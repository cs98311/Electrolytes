// Water Molecules H-bond Analysis

/*

To Run:
gcc Codes/BondCheckerZn.c -o Codes/BondCheckerZn -lm
./Codes/BondCheckerZn
rm Codes/BondCheckerZn

*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

//////////////////////////////////////////////
// Defining Global Variables for System Info//
//////////////////////////////////////////////

// Box length of the cubic system
double box_length;

// Count of no. of clusters
int count = 0;

// Current timestep
int timestep = 0;

// numW = number of Water molecules in the system
// numT = number of TFSI molecules in the system
// numL = number of Lithium ions in the system
// numZ = number of Zinc ions in the system
int numW, numT, numL, numZ;

// Global variables for the counting unique H-bonds/Coordination with H2O,TFSI, Li,Zn
double countW = 0, countA = 0, countL = 0, countZ = 0;

double Wrt[8500][6], Ort[400][7], Frt[300][7], Nrt[1000][6];
double Lrt[2][4], Zrt[4000][4];
int Wcounter = 0, Fcounter = 0, Ocounter = 0, Ncounter = 0;
int Lcounter = 0, Zcounter = 0;

double **exchange;

// int **Adj_mat;

///////////////////////
// DEFINING FUNCTIONS //
///////////////////////

// Function to get magnitude of 2 vectors(of size=3)
double get_magnitude(double m[3])
{
	double mag = 0;
	mag = sqrt(m[0] * m[0] + m[1] * m[1] + m[2] * m[2]);
	return mag;
}

// Function to get dot product of 2 vectors(of size=3)
double dot_product(double m[3], double n[3])
{

	double dt = 0.0;
	int i = 0;
	for (i = 0; i < 3; i++)
	{
		dt += m[i] * n[i];
	}
	return dt;
}

// Function to get angle between 2 vectors(of size=3)
double vector_angle(double v1[3], double v2[3])
{

	double dot = 0, mag1 = 0, mag2 = 0;
	dot = dot_product(v1, v2);
	mag1 = get_magnitude(v1);
	mag2 = get_magnitude(v2);
	return acos((dot) / (mag1 * mag2));
}

// Function to swap two numbers
void swap(int *xp, int *yp)
{
	int temp = *xp;
	*xp = *yp;
	*yp = temp;
}

// Function to apply selection sort to 1D array of size=n
void selectionSort(int arr[], int n)
{
	int i, j, min_idx;
	for (i = 0; i < n - 1; i++)
	{
		min_idx = i;
		for (j = i + 1; j < n; j++)
			if (arr[j] < arr[min_idx])
				min_idx = j;
		swap(&arr[min_idx], &arr[i]);
	}
}

// PBC(Periodic Boundary Condition)
void PBC(double vector[3])
{
	int j = 0;
	for (j = 0; j < 3; j++)
	{
		if (vector[j] > box_length / 2)
		{
			vector[j] -= box_length;
		}
		else if (vector[j] < -box_length / 2)
		{
			vector[j] += box_length;
		}
	}
}

// // DFS(Depth First Search)
// void DFS(int i, int visited[numW], int components[numW])
// {
// 	visited[i] = 1;
// 	components[i] = count;
// 	int j = 0;

// 	for (j = 0; j < numW; j++)
// 	{
// 		if ((Adj_mat[i][j] == 1) && (visited[j] == 0))
// 		{
// 			DFS(j, visited, components);
// 		}
// 	}
// }

// Function to check if H-bond is formed by current water molecule(m) being checked
// with another water oxygen(i) in the system
void Check_HB_Water(double H1[numW][3], double OW[numW][3], int m, int Donor[numW], int Acceptor[numW], int noH, FILE *fEdges, FILE *fOriOw, FILE *fEdgesWT)
{
	// Takes the arrays for Hydrogen and Oxygen coordinates, iteration number(m) of the water
	// molecule currently being czeched, and the Donor and Acceptor HB counter arrays as parameters

	int i = 0, j = 0;

	// v1: is the vector from Water Oxygen(m) to Hydrogen(m) being checked for
	// v2: is the vector from Water Oxygen(m)
	//     to the Acceptor atom(i) (here oxygen of another water molecule)
	// dL: is the potential HB vector from Hydrogen(m) to to the Acceptor atom(i)
	// distance: is the magnitude of the dL ie. the potential HB
	//           Its checked acc. to geometrical criteria (here < 0.22nm)
	// angle: is the acute angle between the vectors v1 and v2
	//        Its checked acc. to geometrical criteria (here < 35 deg)

	double distance = 0, v1[3], v2[3], angle = 0;
	double dL[3];

	// Calculating v1 vector
	for (j = 0; j < 3; j++)
	{
		v1[j] = H1[m][j] - OW[m][j];
	}
	PBC(v1);

	// Starting iterations over the remaining water molecules in the system to check for HB with them
	// by calculating distance and angle with them and applying geometric check

	// Iteration Starts here
	for (i = 0; i < numW; i++)
	{

		// Avoiding the current water molecule for check
		if (i == m)
		{
			continue;
		}

		// Adjusting for PBC
		for (j = 0; j < 3; j++)
		{
			dL[j] = (OW[i][j] - H1[m][j]);

			if (dL[j] > box_length / 2)
			{
				dL[j] -= box_length;
			}
			else if (dL[j] < -box_length / 2)
			{
				dL[j] += box_length;
			}
		}

		// Calculating the v2 vector through vector addition
		for (j = 0; j < 3; j++)
		{
			v2[j] = v1[j] + dL[j];
		}
		PBC(v2);
		// Calculating the magnitude of potential HB vector
		distance = sqrt(dL[0] * dL[0] + dL[1] * dL[1] + dL[2] * dL[2]);

		// Calculating the angle between v1 and v2 for geometric check
		angle = vector_angle(v1, v2);
		// Converting from radians to degrees
		angle *= 57.2958;

		// Applying the geometric check for HB
		if (distance < 0.25 && angle < 35)
		{
			// Adding for the Donor Hydrogen(m)
			Donor[m] += 1;
			// Adding for the Acceptor oxygen(i)
			Acceptor[i] += 1;
			// Counting twice for both donor and acceptor
			countW += 2;
			// Filling the Adjacency Matrix
			// Adj_mat[i][m] = 1;
			// Adj_mat[m][i] = 1;

			Wrt[Wcounter][0] = m;
			Wrt[Wcounter][1] = noH;
			Wrt[Wcounter][2] = i;
			Wrt[Wcounter][3] = distance;
			Wrt[Wcounter][4] = angle;
			Wrt[Wcounter][5] = 1;

			Wcounter += 1;

			if (noH == 1)
			{
				exchange[2 * m][2] = 1;
				exchange[2 * m][3] = i;
				exchange[2 * m][6] = distance;
				exchange[2 * m][7] = angle;
			}

			else if (noH == 2)
			{
				exchange[2 * m + 1][2] = 1;
				exchange[2 * m + 1][3] = i;
				exchange[2 * m + 1][6] = distance;
				exchange[2 * m + 1][7] = angle;
			}

			fprintf(fEdges, "%d,%d\n", m, i);
			fprintf(fEdgesWT, "%d,%d\n", m, i);
		}

		if (distance < 0.5)
		{
			fprintf(fOriOw, "%.5lf %.4lf\n", distance, angle);
		}
	}
	// Iteration over other water molecules ends here
}

// Function to check if H-bond is formed with TFSI oxygen as acceptor
void Check_HB_Otfsi(double H1[numW][3], double OW[numW][3], double Otfsi[4 * numT][3], int m, int HB_with_Otfsi[numW], int Acceptor_Ot[4 * numT], int noH, FILE *fOriOt, FILE *fEdgesWT)
{

	int i = 0, j = 0;
	double distance = 0, v1[3], v2[3], angle = 0;
	double dL[3];

	for (j = 0; j < 3; j++)
	{
		v1[j] = H1[m][j] - OW[m][j];
	}
	PBC(v1);

	for (i = 0; i < 4 * numT; i++)
	{

		for (j = 0; j < 3; j++)
		{
			dL[j] = (Otfsi[i][j] - H1[m][j]);

			if (dL[j] > box_length / 2)
			{
				dL[j] -= box_length;
			}
			else if (dL[j] < -box_length / 2)
			{
				dL[j] += box_length;
			}
		}

		for (j = 0; j < 3; j++)
		{
			v2[j] = v1[j] + dL[j];
		}
		PBC(v2);
		distance = sqrt(dL[0] * dL[0] + dL[1] * dL[1] + dL[2] * dL[2]);

		angle = vector_angle(v1, v2);
		angle *= 57.2958;

		if (distance < 0.25 && angle < 35)
		{
			HB_with_Otfsi[m] += 1;
			Acceptor_Ot[i] += 1;
			countA += 1;

			Ort[Ocounter][0] = m;
			Ort[Ocounter][1] = noH;
			Ort[Ocounter][2] = i % numT;
			Ort[Ocounter][3] = (i / numT) + 1;
			Ort[Ocounter][4] = distance;
			Ort[Ocounter][5] = angle;
			Ort[Ocounter][6] = 1;

			Ocounter += 1;

			if (noH == 1)
			{
				exchange[2 * m][2] = 2;
				exchange[2 * m][3] = i % numT;
				exchange[2 * m][5] = (i / numT) + 1;
				exchange[2 * m][6] = distance;
				exchange[2 * m][7] = angle;
			}

			else if (noH == 2)
			{
				exchange[2 * m + 1][2] = 2;
				exchange[2 * m + 1][3] = i % numT;
				exchange[2 * m + 1][5] = (i / numT) + 1;
				exchange[2 * m + 1][6] = distance;
				exchange[2 * m + 1][7] = angle;
			}

			fprintf(fEdgesWT, "%d,%d\n", m, numW + i % numT);
		}

		if (distance < 0.5)
		{
			fprintf(fOriOt, "%.5lf %.4lf\n", distance, angle);
		}
	}
}

// Function to check if H-bond is formed with TFSI fluorine as acceptor (PBC applied)
void Check_HB_F(double H1[numW][3], double OW[numW][3], double F[6 * numT][3], int m, int HB_with_F[numW], int noH, FILE *fOriF, FILE *fEdgesWT)
{

	int i = 0, j = 0;
	double distance = 0, v1[3], v2[3], angle = 0;
	double dL[3];

	for (j = 0; j < 3; j++)
	{
		v1[j] = H1[m][j] - OW[m][j];
	}
	PBC(v1);
	for (i = 0; i < 6 * numT; i++)
	{

		for (j = 0; j < 3; j++)
		{
			dL[j] = (F[i][j] - H1[m][j]);

			if (dL[j] > box_length / 2)
			{
				dL[j] -= box_length;
			}
			else if (dL[j] < -box_length / 2)
			{
				dL[j] += box_length;
			}
		}

		for (j = 0; j < 3; j++)
		{
			v2[j] = v1[j] + dL[j];
		}
		PBC(v2);

		distance = sqrt(dL[0] * dL[0] + dL[1] * dL[1] + dL[2] * dL[2]);

		angle = vector_angle(v1, v2);
		angle *= 57.2958;

		if (distance < 0.25 && angle < 35)
		{
			HB_with_F[m] += 1;
			countA += 1;

			Frt[Fcounter][0] = m;
			Frt[Fcounter][1] = noH;
			Frt[Fcounter][2] = i % numT;
			Frt[Fcounter][3] = (i / numT) + 1;
			Frt[Fcounter][4] = distance;
			Frt[Fcounter][5] = angle;
			Frt[Fcounter][6] = 1;

			Fcounter += 1;

			if (noH == 1)
			{
				exchange[2 * m][2] = 3;
				exchange[2 * m][3] = i % numT;
				exchange[2 * m][5] = (i / numT) + 1;
				exchange[2 * m][4] = (i / (3 * numT)) + 1;
				exchange[2 * m][6] = distance;
				exchange[2 * m][7] = angle;
			}

			else if (noH == 2)
			{
				exchange[2 * m + 1][2] = 3;
				exchange[2 * m + 1][3] = i % numT;
				exchange[2 * m + 1][5] = (i / numT) + 1;
				exchange[2 * m + 1][4] = (i / (3 * numT)) + 1;
				exchange[2 * m + 1][6] = distance;
				exchange[2 * m + 1][7] = angle;
			}

			fprintf(fEdgesWT, "%d,%d\n", m, numW + i % numT);
		}

		if (distance < 0.5)
		{
			fprintf(fOriF, "%.5lf %.4lf\n", distance, angle);
		}
	}
}

// Function to check if H-bond is formed with TFSI nitrogen as acceptor (PBC applied)
void Check_HB_N(double H1[numW][3], double OW[numW][3], double N1[numT][3], int m, int HB_with_N[numW], int noH, FILE *fOriN, FILE *fEdgesWT)
{

	int i = 0, j = 0;
	double distance = 0, v1[3], v2[3], angle = 0;
	double dL[3];

	for (j = 0; j < 3; j++)
	{
		v1[j] = H1[m][j] - OW[m][j];
	}
	PBC(v1);
	for (i = 0; i < numT; i++)
	{

		for (j = 0; j < 3; j++)
		{
			dL[j] = (N1[i][j] - H1[m][j]);

			if (dL[j] > box_length / 2)
			{
				dL[j] -= box_length;
			}
			else if (dL[j] < -box_length / 2)
			{
				dL[j] += box_length;
			}
		}

		for (j = 0; j < 3; j++)
		{
			v2[j] = v1[j] + dL[j];
		}
		PBC(v2);
		distance = sqrt(dL[0] * dL[0] + dL[1] * dL[1] + dL[2] * dL[2]);

		angle = vector_angle(v1, v2);
		angle *= 57.2958;

		if (distance < 0.25)
		{
			HB_with_N[m] += 1;
			countA += 1;

			Nrt[Ncounter][0] = m;
			Nrt[Ncounter][1] = noH;
			Nrt[Ncounter][2] = i;
			Nrt[Ncounter][3] = distance;
			Nrt[Ncounter][4] = angle;
			Nrt[Ncounter][5] = 1;

			Ncounter += 1;

			if (noH == 1)
			{
				exchange[2 * m][2] = 4;
				exchange[2 * m][3] = i;
				exchange[2 * m][6] = distance;
				exchange[2 * m][7] = angle;
			}

			else if (noH == 2)
			{
				exchange[2 * m + 1][2] = 4;
				exchange[2 * m + 1][3] = i;
				exchange[2 * m + 1][6] = distance;
				exchange[2 * m + 1][7] = angle;
			}

			fprintf(fEdgesWT, "%d,%d\n", m, numW + i % numT);
		}

		if (distance < 0.5)
		{
			fprintf(fOriN, "%.5lf %.4lf\n", distance, angle);
		}
	}
}

// Function to check for coordination of Oxygen(m) of water with Li ion (PBC applied)
void Check_Cood_Li(double OW[numW][3], double LI[numL][3], int m, int Cood_with_Li[numW], int Li_Ow[numL])
{

	int i = 0, j = 0;
	double distance = 0;
	double dL[3];

	for (i = 0; i < numL; i++)
	{

		for (j = 0; j < 3; j++)
		{
			dL[j] = (LI[i][j] - OW[m][j]);

			if (dL[j] > box_length / 2)
			{
				dL[j] -= box_length;
			}
			else if (dL[j] < -box_length / 2)
			{
				dL[j] += box_length;
			}
		}

		distance = sqrt(dL[0] * dL[0] + dL[1] * dL[1] + dL[2] * dL[2]);

		if (distance < 0.28)
		{
			Cood_with_Li[m] += 1;
			Li_Ow[i] += 1;
			countL += 1;

			Lrt[Lcounter][0] = m;
			Lrt[Lcounter][1] = i;
			Lrt[Lcounter][2] = distance;
			Lrt[Lcounter][3] = 1;

			Lcounter += 1;

			// do the exchange thingy
		}
	}
}

// Function to check for coordination of Oxygen(m) of water with Zn ion (PBC applied)
void Check_Cood_Zn(double OW[numW][3], double ZN[numZ][3], int m, int Cood_with_Zn[numW], int Zn_Ow[numZ])
{

	int i = 0, j = 0;
	double distance = 0;
	double dL[3];

	for (i = 0; i < numZ; i++)
	{

		for (j = 0; j < 3; j++)
		{
			dL[j] = (ZN[i][j] - OW[m][j]);

			if (dL[j] > box_length / 2)
			{
				dL[j] -= box_length;
			}
			else if (dL[j] < -box_length / 2)
			{
				dL[j] += box_length;
			}
		}

		distance = sqrt(dL[0] * dL[0] + dL[1] * dL[1] + dL[2] * dL[2]);

		if (distance < 0.28)
		{
			Cood_with_Zn[m] += 1;
			Zn_Ow[i] += 1;
			countZ += 1;

			Zrt[Zcounter][0] = m;
			Zrt[Zcounter][1] = i;
			Zrt[Zcounter][2] = distance;
			Zrt[Zcounter][3] = 1;

			Zcounter += 1;

			// do the exchange thingy
			exchange[2 * m][8] = 1;
			exchange[2 * m][9] = i;
			exchange[2 * m][10] = distance;

			exchange[2 * m + 1][8] = 1;
			exchange[2 * m + 1][9] = i;
			exchange[2 * m + 1][10] = distance;
		}
	}
}

// Function to check for coordination of Oxygen(m) of TFSI with Li ion (PBC applied)
void Check_Cood_Li_Ot(double Ot[4 * numT][3], double LI[numL][3], int m, int Cood_with_Li_Ot[4 * numT], int Li_Ot[numL])
{

	int i = 0, j = 0;
	double distance = 0;
	double dL[3];

	for (i = 0; i < numL; i++)
	{

		for (j = 0; j < 3; j++)
		{
			dL[j] = (LI[i][j] - Ot[m][j]);

			if (dL[j] > box_length / 2)
			{
				dL[j] -= box_length;
			}
			else if (dL[j] < -box_length / 2)
			{
				dL[j] += box_length;
			}
		}

		distance = sqrt(dL[0] * dL[0] + dL[1] * dL[1] + dL[2] * dL[2]);

		if (distance < 0.28)
		{
			Cood_with_Li_Ot[m] += 1;
			Li_Ot[i] += 1;
		}
	}
}

// Function to check for coordination of Oxygen(m) of TFSI with Zn ion (PBC applied)
void Check_Cood_Zn_Ot(double Ot[4 * numT][3], double ZN[numZ][3], int m, int Cood_with_Zn_Ot[4 * numT], int Zn_Ot[numZ])
{

	int i = 0, j = 0;
	double distance = 0;
	double dL[3];

	for (i = 0; i < numZ; i++)
	{

		for (j = 0; j < 3; j++)
		{
			dL[j] = (ZN[i][j] - Ot[m][j]);

			if (dL[j] > box_length / 2)
			{
				dL[j] -= box_length;
			}
			else if (dL[j] < -box_length / 2)
			{
				dL[j] += box_length;
			}
		}

		distance = sqrt(dL[0] * dL[0] + dL[1] * dL[1] + dL[2] * dL[2]);

		if (distance < 0.28)
		{
			Cood_with_Zn_Ot[m] += 1;
			Zn_Ot[i] += 1;
		}
	}
}

//////////////////
// Main Function /
//////////////////

int main()
{

	// Opening the files for Coordinates, System Info, Results,etc.//

	FILE *fI = fopen("SysInfo.txt", "r");

	FILE *fH1 = fopen("Coordinates/H1.txt", "r");
	FILE *fH2 = fopen("Coordinates/H2.txt", "r");
	FILE *fOW = fopen("Coordinates/OW.txt", "r");

	FILE *fN1 = fopen("Coordinates/N1.txt", "r");

	FILE *fLI = fopen("Coordinates/LI.txt", "r");
	FILE *fZN = fopen("Coordinates/ZN.txt", "r");

	// Combined fluorine and oxygen
	FILE *fF = fopen("Coordinates/F.txt", "r");
	FILE *fOt = fopen("Coordinates/Ot.txt", "r");

	// Results storage files
	FILE *fHB_water = fopen("Plots/DataFiles/HB_with_water.txt", "w");
	FILE *fCombined = fopen("Plots/DataFiles/Combined.csv", "w");
	FILE *fCounter = fopen("Water/HBperWater.txt", "a");
	FILE *fHBP = fopen("Water/HBpercent.txt", "a");
	FILE *fHB5 = fopen("Water/BifurcatedHB.txt", "a");
	FILE *fDonor = fopen("Water/Donor.txt", "a");
	FILE *fLi_Ow = fopen("Cation/Li_Ow.txt", "a");
	FILE *fZn_Ow = fopen("Cation/Zn_Ow.txt", "a");
	FILE *fLi_Ot = fopen("Cation/Li_Ot.txt", "a");
	FILE *fZn_Ot = fopen("Cation/Zn_Ot.txt", "a");
	FILE *fH_Ot = fopen("Anion/Ot_HB_accepted.txt", "a");
	FILE *fAccept = fopen("Water/Acceptor.txt", "a");
	FILE *fOt_Li = fopen("Anion/Ot_Li_Coodn.txt", "a");
	FILE *fOt_Zn = fopen("Anion/Ot_Zn_Coodn.txt", "a");
	// FILE *file5 = fopen("Water/MyFile.csv", "a");
	// FILE *file6 = fopen("Water/FreqDistri.csv", "a");
	// FILE *file7 = fopen("Water/FreqDistri.txt", "a");
	// FILE *file8 = fopen("Water/ClusterPercent.txt", "a");
	// FILE *file9 = fopen("Water/ClusterPercent1.csv", "a");
	FILE *fRtime = fopen("RT/RtimesW.txt", "a");
	FILE *fEdges = fopen("Water/Edges.csv", "w");
	FILE *fEdgesWT = fopen("Combined/EdgesT.csv", "w");
	FILE *fOriF = fopen("Orientation/F.txt", "a");
	FILE *fOriN = fopen("Orientation/N.txt", "a");
	FILE *fOriOw = fopen("Orientation/Ow.txt", "a");
	FILE *fOriOt = fopen("Orientation/Ot.txt", "a");

	if (fI == NULL || fH1 == NULL || fH2 == NULL || fOW == NULL || fN1 == NULL ||
		fLI == NULL || fZN == NULL || fF == NULL || fOt == NULL || fHB_water == NULL ||
		fCombined == NULL || fCounter == NULL || fHBP == NULL || fHB5 == NULL ||
		fDonor == NULL || fLi_Ow == NULL || fZn_Ow == NULL || fLi_Ot == NULL || fZn_Ot == NULL ||
		fH_Ot == NULL || fAccept == NULL || fOt_Li == NULL || fOt_Zn == NULL ||
		fRtime == NULL || fEdges == NULL || fEdgesWT == NULL || fOriF == NULL || fOriN == NULL ||
		fOriOw == NULL || fOriOt == NULL)
	{
		perror("Error opening one or more files");
		exit(EXIT_FAILURE);
	}

	// Defining and initializing any used variables
	int i = 0, j = 0, k = 0, m = 0;

	// Scanning the current timestep/iteration
	// int timestep = 0;
	scanf("%d", &timestep);

	// Scanning the System Info into global variables
	fscanf(fI, "%lf %d %d %d %d", &box_length, &numW, &numT, &numL, &numZ);

	// exchange
	// Allocate memory for the rows (array of pointers)
	exchange = (double **)malloc(2 * numW * sizeof(double *));

	if (exchange == NULL)
	{
		printf("Memory allocation failed!\n");
		return 1;
	}

	// Allocate memory for each row (array of double values)
	for (int i = 0; i < 2 * numW; i++)
	{
		exchange[i] = (double *)malloc(11 * sizeof(double));
		if (exchange[i] == NULL)
		{
			printf("Memory allocation to exchange[][] failed for row %d!\n", i);
			return 1;
		}
	}

	// Adj_mat
	//  Allocate memory for the rows (array of pointers)
	// Adj_mat = (int **)malloc(numW * sizeof(int *));

	// if (Adj_mat == NULL)
	// {
	// 	printf("Memory allocation failed!\n");
	// 	return 1;
	// }

	// Allocate memory for each row (array of double values)
	// for (int i = 0; i < numW; i++)
	// {
	// 	Adj_mat[i] = (int *)malloc(numW * sizeof(int));
	// 	if (Adj_mat[i] == NULL)
	// 	{
	// 		printf("Memory allocation to Adj_mat[][] failed for row %d!\n", i);
	// 		return 1;
	// 	}
	// }

	// Defining arrays for all elements of the system
	double H1[numW][3], H2[numW][3], OW[numW][3];
	double F[6 * numT][3], N1[numT][3], Ot[4 * numT][3];
	double F1[numT][3], F2[numT][3], F3[numT][3], F4[numT][3], F5[numT][3], F6[numT][3];
	double O1[numT][3], O2[numT][3], O3[numT][3], O4[numT][3];
	double LI[numL][3], ZN[numZ][3];

	// Defining the arrays used to store H-bond/Coordination info for each water molecule
	int Acceptor[numW], Donor[numW], Total_water_HB[numW], Total_non_water_HB[numW], Total_HB[numW];

	int HB_with_F[numW], HB_with_Otfsi[numW], HB_with_N[numW];

	int Cood_with_Li[numW], Cood_with_Zn[numW], Total_cood_with_cation[numW];

	int Zn_Ow[numZ], Li_Ow[numL], Zn_Ot[numZ], Li_Ot[numL];

	int Cood_with_Li_Ot[4 * numT], Cood_with_Zn_Ot[4 * numT];

	int Acceptor_Ot[4 * numT];

	// This one causes segfault; (5000x5000 array)
	// int Adj_mat[numW][numW];

	int visited[numW], components[numW];

	// Scanning the coordinates from respective files //

	for (i = 0; i < numW; i++)
	{
		for (j = 0; j < 3; j++)
		{
			fscanf(fH1, "%lf", &H1[i][j]);
			fscanf(fH2, "%lf", &H2[i][j]);
			fscanf(fOW, "%lf", &OW[i][j]);
		}
	}

	for (i = 0; i < 6 * numT; i++)
	{
		for (j = 0; j < 3; j++)
		{
			fscanf(fF, "%lf", &F[i][j]);
		}
	}

	for (i = 0; i < 4 * numT; i++)
	{
		for (j = 0; j < 3; j++)
		{
			fscanf(fOt, "%lf", &Ot[i][j]);
		}
	}

	for (i = 0; i < numL; i++)
	{
		for (j = 0; j < 3; j++)
		{
			fscanf(fLI, "%lf", &LI[i][j]);
		}
	}

	for (i = 0; i < numZ; i++)
	{
		for (j = 0; j < 3; j++)
		{
			fscanf(fZN, "%lf", &ZN[i][j]);
		}
	}

	for (i = 0; i < numT; i++)
	{
		for (j = 0; j < 3; j++)
		{
			fscanf(fN1, "%lf", &N1[i][j]);
		}
	}

	// Initializing the array elements to 0 //
	for (i = 0; i < numW; i++)
	{
		Acceptor[i] = 0;
	}

	for (i = 0; i < numW; i++)
	{
		Donor[i] = 0;
	}

	for (i = 0; i < numW; i++)
	{
		Total_water_HB[i] = 0;
	}

	for (i = 0; i < numW; i++)
	{
		Total_non_water_HB[i] = 0;
	}

	for (i = 0; i < numW; i++)
	{
		Total_HB[i] = 0;
	}

	for (i = 0; i < numW; i++)
	{
		HB_with_F[i] = 0;
	}
	for (i = 0; i < numW; i++)
	{
		HB_with_N[i] = 0;
	}
	for (i = 0; i < numW; i++)
	{
		HB_with_Otfsi[i] = 0;
	}

	for (i = 0; i < numW; i++)
	{
		Cood_with_Li[i] = 0;
	}

	for (i = 0; i < numW; i++)
	{
		Cood_with_Zn[i] = 0;
	}

	for (i = 0; i < numW; i++)
	{
		Total_cood_with_cation[i] = 0;
	}

	for (i = 0; i < numZ; i++)
	{
		Zn_Ow[i] = 0;
		Zn_Ot[i] = 0;
	}

	for (i = 0; i < numL; i++)
	{
		Li_Ow[i] = 0;
		Li_Ot[i] = 0;
	}

	for (i = 0; i < 4 * numT; i++)
	{
		Cood_with_Li_Ot[i] = 0;
	}

	for (i = 0; i < 4 * numT; i++)
	{
		Cood_with_Zn_Ot[i] = 0;
	}

	for (i = 0; i < 4 * numT; i++)
	{
		Acceptor_Ot[i] = 0;
	}

	// Adj_mat[][]
	// for (i = 0; i < numW; i++)
	// {
	// 	for (j = 0; j < numW; j++)
	// 	{
	// 		Adj_mat[i][j] = 0;
	// 	}
	// }

	// visited[]
	for (i = 0; i < numW; i++)
	{
		visited[i] = 0;
	}

	// components[]
	for (i = 0; i < numW; i++)
	{
		components[i] = 0;
	}

	// exchange
	for (i = 0; i < 2 * numW; i++)
	{
		for (j = 0; j < 11; j++)
		{
			exchange[i][j] = 0;
		}
	}

	///////////////////////
	// Starting Iterations//
	///////////////////////

	// Iteration checks for each water molecule for all
	// H-bonds with other water, TFSI and Coordination with Zn/Li

	for (m = 0; m < numW; m++)
	{

		exchange[m * 2][0] = m;
		exchange[m * 2][1] = 1;

		exchange[m * 2 + 1][0] = m;
		exchange[m * 2 + 1][1] = 2;

		// Check for H1 hydogen first
		Check_HB_Water(H1, OW, m, Donor, Acceptor, 1, fEdges, fOriOw, fEdgesWT);

		Check_HB_Otfsi(H1, OW, Ot, m, HB_with_Otfsi, Acceptor_Ot, 1, fOriOt, fEdgesWT);

		Check_HB_F(H1, OW, F, m, HB_with_F, 1, fOriF, fEdgesWT);

		Check_HB_N(H1, OW, N1, m, HB_with_N, 2, fOriN, fEdgesWT);

		// Check for H2 hydrogen
		Check_HB_Water(H2, OW, m, Donor, Acceptor, 2, fEdges, fOriOw, fEdgesWT);

		Check_HB_Otfsi(H2, OW, Ot, m, HB_with_Otfsi, Acceptor_Ot, 2, fOriOt, fEdgesWT);

		Check_HB_F(H2, OW, F, m, HB_with_F, 2, fOriF, fEdgesWT);

		Check_HB_N(H2, OW, N1, m, HB_with_N, 2, fOriN, fEdgesWT);

		// Check for OW oxygen coordination with cation(s)
		Check_Cood_Li(OW, LI, m, Cood_with_Li, Li_Ow);
		Check_Cood_Zn(OW, ZN, m, Cood_with_Zn, Zn_Ow);
	}

	for (m = 0; m < 4 * numT; m++)
	{
		Check_Cood_Li_Ot(Ot, LI, m, Cood_with_Li_Ot, Li_Ot);
		Check_Cood_Zn_Ot(Ot, ZN, m, Cood_with_Zn_Ot, Zn_Ot);
	}

	///////////////////////
	// Iterations Finished//
	///////////////////////

	//////////////////////////
	// CONNECTED COMPONENTS //
	//////////////////////////

	// for (i = 0; i < numW; i++)
	// {
	// 	if (visited[i] == 0)
	// 	{
	// 		count += 1;
	// 		DFS(i, visited, components);
	// 	}
	// }

	// int ClusterSizes[count];

	// Initializing ClusterSizes[] to all zero
	// for (i = 0; i < count; i++)
	// {
	// 	ClusterSizes[i] = 0;
	// }

	// for (i = 1; i <= count; i++)
	// {
	// 	for (j = 0; j < numW; j++)
	// 	{
	// 		if (components[j] == i)
	// 		{
	// 			ClusterSizes[i - 1] += 1;
	// 		}
	// 	}
	// }

	// selectionSort(ClusterSizes, count);

	// FILE *fCbig = fopen("Water/Largest.txt", "a");
	// fprintf(fCbig, "%d ", ClusterSizes[count - 1]);
	// fclose(fCbig);

	//////////////////////////
	// Operations on arrays //
	//////////////////////////

	for (i = 0; i < numW; i++)
	{
		Total_water_HB[i] = Donor[i] + Acceptor[i];
	}

	for (i = 0; i < numW; i++)
	{
		Total_non_water_HB[i] = HB_with_N[i] + HB_with_F[i] + HB_with_Otfsi[i];
	}

	for (i = 0; i < numW; i++)
	{
		Total_HB[i] = Total_water_HB[i] + Total_non_water_HB[i];
	}

	for (i = 0; i < numW; i++)
	{
		Total_cood_with_cation[i] = Cood_with_Li[i] + Cood_with_Zn[i];
	}

	////////////////////////
	// Cation Coordinates //
	////////////////////////

	FILE *fCtraj = fopen("Cation/Traj.txt", "a");
	fprintf(fCtraj, "%lf %lf %lf\n", LI[0][0], LI[0][1], LI[0][2]);
	fclose(fCtraj);

	////////////////////////////
	// Exchange Probabilities //
	////////////////////////////

	FILE *fTS1 = fopen("RT/Timestep.txt", "r");
	int TS1 = 0;
	fscanf(fTS1, "%d", &TS1);
	fclose(fTS1);

	double previous[numW * 2][11];

	if (timestep == TS1 + 1)
	{

		FILE *fEx1 = fopen("Exchange/Matrix.txt", "r");
		for (i = 0; i < numW * 2; i++)
		{
			for (j = 0; j < 11; j++)
			{
				fscanf(fEx1, "%lf ", &previous[i][j]);
			}
		}
		fclose(fEx1);

		FILE *fPrevious = fopen("Exchange/Previous.txt", "w");
		for (i = 0; i < numW * 2; i++)
		{
			for (j = 0; j < 11; j++)
			{
				fprintf(fPrevious, "%.2lf ", previous[i][j]);
			}
			fprintf(fPrevious, "\n");
		}
		fclose(fPrevious);

		FILE *fEx3 = fopen("Exchange/Exchanges.txt", "a");
		for (i = 0; i < 2 * numW; i++)
		{

			// 1 : Ow //

			if ((previous[i][2] == 1) && (exchange[i][2] == 1) && (previous[i][3] == exchange[i][3]))
			{
				fprintf(fEx3, "1 ");
			}

			else if ((previous[i][2] == 1) && (exchange[i][2] == 1) && (previous[i][3] != exchange[i][3]))
			{
				fprintf(fEx3, "2 ");
			}

			else if ((previous[i][2] == 1) && (exchange[i][2] == 2))
			{
				fprintf(fEx3, "3 ");
			}

			else if ((previous[i][2] == 1) && (exchange[i][2] == 3))
			{
				fprintf(fEx3, "4 ");
			}

			else if ((previous[i][2] == 1) && (exchange[i][2] == 4))
			{
				fprintf(fEx3, "5 ");
			}

			else if ((previous[i][2] == 1) && (exchange[i][2] == 0))
			{
				fprintf(fEx3, "6 ");
			}

			// 2 : Ot //

			else if ((previous[i][2] == 2) && (exchange[i][2] == 2) && (previous[i][3] == exchange[i][3]) && (previous[i][5] == exchange[i][5]))
			{
				fprintf(fEx3, "7 ");
			}

			else if ((previous[i][2] == 2) && (exchange[i][2] == 2) && (previous[i][3] == exchange[i][3]) && (previous[i][5] != exchange[i][5]))
			{
				fprintf(fEx3, "8 ");
			}

			else if ((previous[i][2] == 2) && (exchange[i][2] == 2) && (previous[i][3] != exchange[i][3]))
			{
				fprintf(fEx3, "9 ");
			}

			else if ((previous[i][2] == 2) && (exchange[i][2] == 1))
			{
				fprintf(fEx3, "10 ");
			}

			else if ((previous[i][2] == 2) && (exchange[i][2] == 3) && (previous[i][3] == exchange[i][3]))
			{
				fprintf(fEx3, "11 ");
			}

			else if ((previous[i][2] == 2) && (exchange[i][2] == 3) && (previous[i][3] != exchange[i][3]))
			{
				fprintf(fEx3, "12 ");
			}

			else if ((previous[i][2] == 2) && (exchange[i][2] == 4) && (previous[i][3] == exchange[i][3]))
			{
				fprintf(fEx3, "13 ");
			}

			else if ((previous[i][2] == 2) && (exchange[i][2] == 4) && (previous[i][3] != exchange[i][3]))
			{
				fprintf(fEx3, "14 ");
			}

			else if ((previous[i][2] == 2) && (exchange[i][2] == 0))
			{
				fprintf(fEx3, "15 ");
			}

			// 3 : F

			else if ((previous[i][2] == 3) && (exchange[i][2] == 3) && (previous[i][3] == exchange[i][3]) && (previous[i][4] == exchange[i][4]) && (previous[i][5] == exchange[i][5]))
			{
				fprintf(fEx3, "16 ");
			}

			else if ((previous[i][2] == 3) && (exchange[i][2] == 3) && (previous[i][3] == exchange[i][3]) && (previous[i][4] == exchange[i][4]) && (previous[i][5] != exchange[i][5]))
			{
				fprintf(fEx3, "17 ");
			}

			else if ((previous[i][2] == 3) && (exchange[i][2] == 3) && (previous[i][3] == exchange[i][3]) && (previous[i][4] != exchange[i][4]))
			{
				fprintf(fEx3, "18 ");
			}

			else if ((previous[i][2] == 3) && (exchange[i][2] == 3) && (previous[i][3] != exchange[i][3]))
			{
				fprintf(fEx3, "19 ");
			}

			else if ((previous[i][2] == 3) && (exchange[i][2] == 1))
			{
				fprintf(fEx3, "20 ");
			}

			else if ((previous[i][2] == 3) && (exchange[i][2] == 2) && (previous[i][3] == exchange[i][3]))
			{
				fprintf(fEx3, "21 ");
			}

			else if ((previous[i][2] == 3) && (exchange[i][2] == 2) && (previous[i][3] != exchange[i][3]))
			{
				fprintf(fEx3, "22 ");
			}

			else if ((previous[i][2] == 3) && (exchange[i][2] == 4) && (previous[i][3] == exchange[i][3]))
			{
				fprintf(fEx3, "23 ");
			}

			else if ((previous[i][2] == 3) && (exchange[i][2] == 4) && (previous[i][3] != exchange[i][3]))
			{
				fprintf(fEx3, "24 ");
			}

			else if ((previous[i][2] == 3) && (exchange[i][2] == 0))
			{
				fprintf(fEx3, "25 ");
			}

			// 4 : N

			else if ((previous[i][2] == 4) && (exchange[i][2] == 4) && (previous[i][3] == exchange[i][3]))
			{
				fprintf(fEx3, "26 ");
			}

			else if ((previous[i][2] == 4) && (exchange[i][2] == 4) && (previous[i][3] != exchange[i][3]))
			{
				fprintf(fEx3, "27 ");
			}

			else if ((previous[i][2] == 4) && (exchange[i][2] == 1))
			{
				fprintf(fEx3, "28 ");
			}

			else if ((previous[i][2] == 4) && (exchange[i][2] == 2) && (previous[i][3] == exchange[i][3]))
			{
				fprintf(fEx3, "29 ");
			}

			else if ((previous[i][2] == 4) && (exchange[i][2] == 2) && (previous[i][3] != exchange[i][3]))
			{
				fprintf(fEx3, "30 ");
			}

			else if ((previous[i][2] == 4) && (exchange[i][2] == 3) && (previous[i][3] == exchange[i][3]))
			{
				fprintf(fEx3, "31 ");
			}

			else if ((previous[i][2] == 4) && (exchange[i][2] == 3) && (previous[i][3] != exchange[i][3]))
			{
				fprintf(fEx3, "32 ");
			}

			else if ((previous[i][2] == 4) && (exchange[i][2] == 0))
			{
				fprintf(fEx3, "33 ");
			}

			// 0 : No Bond

			else if ((previous[i][2] == 0) && (exchange[i][2] == 0))
			{
				fprintf(fEx3, "34 ");
			}

			else if ((previous[i][2] == 0) && (exchange[i][2] == 1))
			{
				fprintf(fEx3, "35 ");
			}

			else if ((previous[i][2] == 0) && (exchange[i][2] == 2))
			{
				fprintf(fEx3, "36 ");
			}

			else if ((previous[i][2] == 0) && (exchange[i][2] == 3))
			{
				fprintf(fEx3, "37 ");
			}

			else if ((previous[i][2] == 0) && (exchange[i][2] == 4))
			{
				fprintf(fEx3, "38 ");
			}
		}
		fclose(fEx3);

		FILE *fEx2 = fopen("Exchange/Matrix.txt", "w");
		for (i = 0; i < numW * 2; i++)
		{
			for (j = 0; j < 11; j++)
			{
				fprintf(fEx2, "%.2lf ", exchange[i][j]);
			}
			fprintf(fEx2, "\n");
		}
		fclose(fEx2);
	}

	else
	{
		FILE *fEx = fopen("Exchange/Matrix.txt", "w");
		for (i = 0; i < numW * 2; i++)
		{
			for (j = 0; j < 11; j++)
			{
				fprintf(fEx, "%.2lf ", exchange[i][j]);
			}
			fprintf(fEx, "\n");
		}
		fclose(fEx);
	}

	////////////////////////
	// H2O Donor-Acceptor //
	////////////////////////

	FILE *fH2Oda = fopen("Water/Donor_Acceptor.txt", "a");

	double d0a0 = 0, d0a1 = 0, d0a2 = 0, d1a0 = 0, d1a1 = 0, d1a2 = 0, d2a0 = 0, d2a1 = 0, d2a2 = 0, totda = 0;

	for (i = 0; i < numW; i++)
	{
		if (Donor[i] == 0 && Acceptor[i] == 0)
			d0a0 += 1;
		else if (Donor[i] == 0 && Acceptor[i] == 1)
			d0a1 += 1;
		else if (Donor[i] == 0 && Acceptor[i] == 2)
			d0a2 += 1;
		else if (Donor[i] == 1 && Acceptor[i] == 0)
			d1a0 += 1;
		else if (Donor[i] == 1 && Acceptor[i] == 1)
			d1a1 += 1;
		else if (Donor[i] == 1 && Acceptor[i] == 2)
			d1a2 += 1;
		else if (Donor[i] == 2 && Acceptor[i] == 0)
			d2a0 += 1;
		else if (Donor[i] == 2 && Acceptor[i] == 1)
			d2a1 += 1;
		else if (Donor[i] == 2 && Acceptor[i] == 2)
			d2a2 += 1;
	}

	totda = d0a0 + d0a1 + d0a2 + d1a0 + d1a1 + d1a2 + d2a0 + d2a1 + d2a2;

	fprintf(fH2Oda, "%.2lf %.2lf %.2lf %.2lf %.2lf %.2lf %.2lf %.2lf %.2lf\n", d0a0 * 100 / totda, d0a1 * 100 / totda, d0a2 * 100 / totda, d1a0 * 100 / totda, d1a1 * 100 / totda, d1a2 * 100 / totda, d2a0 * 100 / totda, d2a1 * 100 / totda, d2a2 * 100 / totda);

	fclose(fH2Oda);

	//////////////////
	// Ot Plot Data //
	//////////////////

	FILE *fOtCompZn = fopen("Anion/OtCompiledZn.txt", "a");

	double t0l0 = 0, t0l1 = 0, t1l0 = 0, t1l1 = 0, tot = 0;

	for (i = 0; i < 4 * numT; i++)
	{
		if (Acceptor_Ot[i] == 0 && Cood_with_Zn_Ot[i] == 0)
			t0l0 += 1;
		else if (Acceptor_Ot[i] == 0 && Cood_with_Zn_Ot[i] == 1)
			t0l1 += 1;
		else if (Acceptor_Ot[i] == 1 && Cood_with_Zn_Ot[i] == 0)
			t1l0 += 1;
		else if (Acceptor_Ot[i] == 1 && Cood_with_Zn_Ot[i] == 1)
			t1l1 += 1;
	}
	tot = t0l1 + t0l0 + t1l1 + t1l0;

	fprintf(fOtCompZn, "%.2lf %.2lf %.2lf %.2lf\n", t0l0 * 100 / tot, t0l1 * 100 / tot, t1l0 * 100 / tot, t1l1 * 100 / tot);

	fclose(fOtCompZn);

	////////////////////
	// RESIDENCE TIME //
	////////////////////

	FILE *fTS = fopen("RT/Timestep.txt", "r");
	int TS = 0, Wsize = 0, Fsize = 0, Osize = 0, Nsize = 0, Lsize = 0, Zsize = 0;
	fscanf(fTS, "%d %d %d %d %d %d %d", &TS, &Wsize, &Fsize, &Osize, &Nsize, &Lsize, &Zsize);
	fclose(fTS);

	/////////
	// H2O //
	/////////

	double Prev[Wsize][6];

	int cc = 0;

	if (timestep == (TS + 1))
	{

		FILE *fH20hb1 = fopen("RT/Wrt.txt", "r");

		for (i = 0; i < Wsize; i++)
		{
			for (j = 0; j < 6; j++)
			{
				fscanf(fH20hb1, "%lf", &Prev[i][j]);
			}
		}

		fclose(fH20hb1);

		for (i = 0; i < Wsize; i++)
		{
			cc = 0;
			for (j = 0; j < Wcounter; j++)
			{
				if (((int)Prev[i][0] == (int)(Wrt[j][0])) && ((int)Prev[i][1] == (int)(Wrt[j][1])) && ((int)Prev[i][2] == (int)(Wrt[j][2])))
				{
					Wrt[j][5] = Prev[i][5] + 1;
					cc = 1;
				}
			}
			if (cc == 0)
			{
				fprintf(fRtime, "%d ", (int)Prev[i][5]);
			}
		}

		FILE *fH20hb2 = fopen("RT/Wrt.txt", "w");
		for (i = 0; i < Wcounter; i++)
		{
			for (j = 0; j < 6; j++)
			{
				fprintf(fH20hb2, "%.2lf ", Wrt[i][j]);
			}
			fprintf(fH20hb2, "\n");
		}
		fclose(fH20hb2);
	}

	else
	{

		FILE *fH20hb3 = fopen("RT/Wrt.txt", "w");
		for (i = 0; i < Wcounter; i++)
		{
			for (j = 0; j < 6; j++)
			{
				fprintf(fH20hb3, "%.2lf ", Wrt[i][j]);
			}
			fprintf(fH20hb3, "\n");
		}
		fclose(fH20hb3);
	}

	//////////////
	// FLUORINE //
	//////////////
	FILE *fRtimeF = fopen("RT/RtimesF.txt", "a");

	double PrevF[Fsize][7];

	cc = 0;

	if (timestep == (TS + 1))
	{

		FILE *fFlu1 = fopen("RT/Frt.txt", "r");
		for (i = 0; i < Fsize; i++)
		{
			for (j = 0; j < 7; j++)
			{
				fscanf(fFlu1, "%lf", &PrevF[i][j]);
			}
		}
		fclose(fFlu1);

		for (i = 0; i < Fsize; i++)
		{
			cc = 0;
			for (j = 0; j < Fcounter; j++)
			{
				if (((int)PrevF[i][0] == (int)(Frt[j][0])) && ((int)PrevF[i][1] == (int)(Frt[j][1])) && ((int)PrevF[i][2] == (int)(Frt[j][2])) && ((int)PrevF[i][3] == (int)(Frt[j][3])))
				{
					Frt[j][6] = PrevF[i][6] + 1;
					cc = 1;
				}
			}
			if (cc == 0)
			{
				fprintf(fRtimeF, "%d ", (int)PrevF[i][6]);
			}
		}

		FILE *fFlu2 = fopen("RT/Frt.txt", "w");
		for (i = 0; i < Fcounter; i++)
		{
			for (j = 0; j < 7; j++)
			{
				fprintf(fFlu2, "%.2lf ", Frt[i][j]);
			}
			fprintf(fFlu2, "\n");
		}
		fclose(fFlu2);
	}

	else
	{

		FILE *fFlu = fopen("RT/Frt.txt", "w");
		for (i = 0; i < Fcounter; i++)
		{
			for (j = 0; j < 7; j++)
			{
				fprintf(fFlu, "%.2lf ", Frt[i][j]);
			}
			fprintf(fFlu, "\n");
		}
		fclose(fFlu);
	}

	////////////
	// OXYGEN //
	////////////
	FILE *fRtimeO = fopen("RT/RtimesO.txt", "a");

	double PrevO[Osize][7];

	cc = 0;

	if (timestep == (TS + 1))
	{

		FILE *fOxy1 = fopen("RT/Ort.txt", "r");
		for (i = 0; i < Osize; i++)
		{
			for (j = 0; j < 7; j++)
			{
				fscanf(fOxy1, "%lf", &PrevO[i][j]);
			}
		}
		fclose(fOxy1);

		for (i = 0; i < Osize; i++)
		{
			cc = 0;
			for (j = 0; j < Ocounter; j++)
			{
				if (((int)PrevO[i][0] == (int)(Ort[j][0])) && ((int)PrevO[i][1] == (int)(Ort[j][1])) && ((int)PrevO[i][2] == (int)(Ort[j][2])) && ((int)PrevO[i][3] == (int)(Ort[j][3])))
				{
					Ort[j][6] = PrevO[i][6] + 1;
					cc = 1;
				}
			}
			if (cc == 0)
			{
				fprintf(fRtimeO, "%d ", (int)PrevO[i][6]);
			}
		}

		FILE *fOxy2 = fopen("RT/Ort.txt", "w");
		for (i = 0; i < Ocounter; i++)
		{
			for (j = 0; j < 7; j++)
			{
				fprintf(fOxy2, "%.2lf ", Ort[i][j]);
			}
			fprintf(fOxy2, "\n");
		}
		fclose(fOxy2);
	}

	else
	{

		FILE *fOxy = fopen("RT/Ort.txt", "w");
		for (i = 0; i < Ocounter; i++)
		{
			for (j = 0; j < 7; j++)
			{
				fprintf(fOxy, "%.2lf ", Ort[i][j]);
			}
			fprintf(fOxy, "\n");
		}
		fclose(fOxy);
	}

	//////////////
	// NITROGEN //
	//////////////
	FILE *fRtimeN = fopen("RT/RtimesN.txt", "a");

	double PrevN[Nsize][7];

	cc = 0;

	if (timestep == (TS + 1))
	{

		FILE *fNxy1 = fopen("RT/Nrt.txt", "r");
		for (i = 0; i < Nsize; i++)
		{
			for (j = 0; j < 6; j++)
			{
				fscanf(fNxy1, "%lf", &PrevN[i][j]);
			}
		}
		fclose(fNxy1);

		for (i = 0; i < Nsize; i++)
		{
			cc = 0;
			for (j = 0; j < Ncounter; j++)
			{
				if (((int)PrevN[i][0] == (int)(Nrt[j][0])) && ((int)PrevN[i][1] == (int)(Nrt[j][1])) && ((int)PrevN[i][2] == (int)(Nrt[j][2])))
				{
					Nrt[j][5] = PrevN[i][5] + 1;
					cc = 1;
				}
			}
			if (cc == 0)
			{
				fprintf(fRtimeN, "%d ", (int)PrevN[i][5]);
			}
		}

		FILE *fNxy2 = fopen("RT/Nrt.txt", "w");
		for (i = 0; i < Ncounter; i++)
		{
			for (j = 0; j < 6; j++)
			{
				fprintf(fNxy2, "%.2lf ", Nrt[i][j]);
			}
			fprintf(fNxy2, "\n");
		}
		fclose(fNxy2);
	}

	else
	{

		FILE *fNxy = fopen("RT/Nrt.txt", "w");
		for (i = 0; i < Ncounter; i++)
		{
			for (j = 0; j < 6; j++)
			{
				fprintf(fNxy, "%.2lf ", Nrt[i][j]);
			}
			fprintf(fNxy, "\n");
		}
		fclose(fNxy);
	}

	/////////////
	// LITHIUM //
	/////////////
	FILE *fRtimeL = fopen("RT/RtimesL.txt", "a");

	double PrevL[Lsize][4];

	cc = 0;

	if (timestep == (TS + 1))
	{

		FILE *fLxy1 = fopen("RT/Lrt.txt", "r");
		for (i = 0; i < Lsize; i++)
		{
			for (j = 0; j < 4; j++)
			{
				fscanf(fLxy1, "%lf", &PrevL[i][j]);
			}
		}
		fclose(fLxy1);

		for (i = 0; i < Lsize; i++)
		{
			cc = 0;
			for (j = 0; j < Lcounter; j++)
			{
				if (((int)PrevL[i][0] == (int)(Lrt[j][0])) && ((int)PrevL[i][1] == (int)(Lrt[j][1])))
				{
					Lrt[j][3] = PrevL[i][3] + 1;
					cc = 1;
				}
			}
			if (cc == 0)
			{
				fprintf(fRtimeL, "%d ", (int)PrevL[i][3]);
			}
		}

		FILE *fLxy2 = fopen("RT/Lrt.txt", "w");
		for (i = 0; i < Lcounter; i++)
		{
			for (j = 0; j < 4; j++)
			{
				fprintf(fLxy2, "%.2lf ", Lrt[i][j]);
			}
			fprintf(fLxy2, "\n");
		}
		fclose(fLxy2);
	}

	else
	{

		FILE *fLxy = fopen("RT/Lrt.txt", "w");
		for (i = 0; i < Lcounter; i++)
		{
			for (j = 0; j < 4; j++)
			{
				fprintf(fLxy, "%.2lf ", Lrt[i][j]);
			}
			fprintf(fLxy, "\n");
		}
		fclose(fLxy);
	}

	//////////
	// ZINC //
	//////////
	FILE *fRtimeZ = fopen("RT/RtimesZ.txt", "a");

	double PrevZ[Zsize][4];

	cc = 0;

	if (timestep == (TS + 1))
	{

		FILE *fZxy1 = fopen("RT/Zrt.txt", "r");
		for (i = 0; i < Zsize; i++)
		{
			for (j = 0; j < 4; j++)
			{
				fscanf(fZxy1, "%lf", &PrevZ[i][j]);
			}
		}
		fclose(fZxy1);

		for (i = 0; i < Zsize; i++)
		{
			cc = 0;
			for (j = 0; j < Zcounter; j++)
			{
				if (((int)PrevZ[i][0] == (int)(Zrt[j][0])) && ((int)PrevZ[i][1] == (int)(Zrt[j][1])))
				{
					Zrt[j][3] = PrevZ[i][3] + 1;
					cc = 1;
				}
			}
			if (cc == 0)
			{
				fprintf(fRtimeZ, "%d ", (int)PrevZ[i][3]);
			}
		}

		FILE *fZxy2 = fopen("RT/Zrt.txt", "w");
		for (i = 0; i < Zcounter; i++)
		{
			for (j = 0; j < 4; j++)
			{
				fprintf(fZxy2, "%.2lf ", Zrt[i][j]);
			}
			fprintf(fZxy2, "\n");
		}
		fclose(fZxy2);
	}

	else
	{

		FILE *fZxy = fopen("RT/Zrt.txt", "w");
		for (i = 0; i < Zcounter; i++)
		{
			for (j = 0; j < 4; j++)
			{
				fprintf(fZxy, "%.2lf ", Zrt[i][j]);
			}
			fprintf(fZxy, "\n");
		}
		fclose(fZxy);
	}

	/////////////
	// LITHIUM //
	/////////////
	FILE *fRtimeD2 = fopen("RT/RtimesD2.txt", "a");
	FILE *fRtimeD1 = fopen("RT/RtimesD1.txt", "a");
	FILE *fRtimeD0 = fopen("RT/RtimesD0.txt", "a");

	int PrevD[numW][2];
	int Drt[numW][2];

	for (i = 0; i < numW; i++)
	{
		Drt[i][0] = Donor[i];
		Drt[i][1] = 1;
	}

	cc = 0;

	if (timestep == (TS + 1))
	{

		FILE *fDxy1 = fopen("RT/Drt.txt", "r");
		for (i = 0; i < numW; i++)
		{
			for (j = 0; j < 2; j++)
			{
				fscanf(fDxy1, "%d", &PrevD[i][j]);
			}
		}
		fclose(fDxy1);

		for (i = 0; i < numW; i++)
		{
			cc = 0;

			if ((PrevD[i][0] == Drt[i][0]))
			{
				Drt[i][1] = PrevD[i][1] + 1;
				cc = 1;
			}

			else
			{
				if (PrevD[i][0] == 2)
				{
					fprintf(fRtimeD2, "%d ", PrevD[i][1]);
				}
				else if (PrevD[i][0] == 1)
				{
					fprintf(fRtimeD1, "%d ", PrevD[i][1]);
				}
				else if (PrevD[i][0] == 0)
				{
					fprintf(fRtimeD0, "%d ", PrevD[i][1]);
				}
			}
		}

		FILE *fDxy2 = fopen("RT/Drt.txt", "w");
		for (i = 0; i < numW; i++)
		{
			for (j = 0; j < 2; j++)
			{
				fprintf(fDxy2, "%d ", Drt[i][j]);
			}
			fprintf(fDxy2, "\n");
		}
		fclose(fDxy2);
	}

	else
	{

		FILE *fDxy = fopen("RT/Drt.txt", "w");
		for (i = 0; i < numW; i++)
		{
			for (j = 0; j < 2; j++)
			{
				fprintf(fDxy, "%d ", Drt[i][j]);
			}
			fprintf(fDxy, "\n");
		}
		fclose(fDxy);
	}

	//////////////////////////////
	// DONOR-ACCEPTOR LIFETIMES //
	//////////////////////////////

	/////////////////////
	FILE *fTS2 = fopen("RT/Timestep.txt", "w");
	fprintf(fTS2, "%d %d %d %d %d %d %d", timestep, Wcounter, Fcounter, Ocounter, Ncounter, Lcounter, Zcounter);
	fclose(fTS2);

	/////////////////////
	//*OUTPUT PRINTING*//
	/////////////////////

	///////////////////
	// Cluster Results /
	///////////////////

	// Frequencies
	// Upper bounds included
	// int one = 0, two = 0, three2five = 0, five2ten = 0, ten2twenty = 0, twenty2fifty = 0, fifty2hundred = 0, hundredplus = 0, threehundplus = 0, fivehundplus = 0;
	// for (i = 0; i < count; i++)
	// {
	// 	if (ClusterSizes[i] == 1)
	// 		one += 1;
	// 	else if (ClusterSizes[i] == 2)
	// 		two += 1;
	// 	else if (ClusterSizes[i] > 2 && ClusterSizes[i] <= 5)
	// 		three2five += 1;
	// 	else if (ClusterSizes[i] > 5 && ClusterSizes[i] <= 10)
	// 		five2ten += 1;
	// 	else if (ClusterSizes[i] > 10 && ClusterSizes[i] <= 20)
	// 		ten2twenty += 1;
	// 	else if (ClusterSizes[i] > 20 && ClusterSizes[i] <= 50)
	// 		twenty2fifty += 1;
	// 	else if (ClusterSizes[i] > 50 && ClusterSizes[i] <= 100)
	// 		fifty2hundred += 1;
	// 	else if (ClusterSizes[i] > 100 && ClusterSizes[i] <= 300)
	// 		hundredplus += 1;
	// 	else if (ClusterSizes[i] > 300 && ClusterSizes[i] <= 500)
	// 		threehundplus += 1;
	// 	if (ClusterSizes[i] > 500)
	// 		fivehundplus += 1;
	// }
	// fprintf(file6, "TS%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n", timestep, one, two, three2five, five2ten, ten2twenty, twenty2fifty, fifty2hundred, hundredplus, threehundplus, fivehundplus);
	// fprintf(file7, "%d %d %d %d %d %d %d %d %d %d\n", one, two, three2five, five2ten, ten2twenty, twenty2fifty, fifty2hundred, hundredplus, threehundplus, fivehundplus);

	// fprintf(file5, "TS%d,", timestep);
	// for (i = count - 1; i > 0; i--)
	// {
	// 	fprintf(file5, "%d,", ClusterSizes[i]);
	// }
	// fprintf(file5, "%d", ClusterSizes[i]);
	// fprintf(file5, "\n");

	// // Percentages
	// int one1 = 0, two1 = 0, three2five1 = 0, five2ten1 = 0, ten2twenty1 = 0, twenty2fifty1 = 0, fifty2hundred1 = 0, hundredplus1 = 0, threehundplus1 = 0, fivehundplus1 = 0;
	// for (i = 0; i < count; i++)
	// {
	// 	if (ClusterSizes[i] == 1)
	// 		one1 += ClusterSizes[i];
	// 	else if (ClusterSizes[i] == 2)
	// 		two1 += ClusterSizes[i];
	// 	else if (ClusterSizes[i] > 2 && ClusterSizes[i] <= 5)
	// 		three2five1 += ClusterSizes[i];
	// 	else if (ClusterSizes[i] > 5 && ClusterSizes[i] <= 10)
	// 		five2ten1 += ClusterSizes[i];
	// 	else if (ClusterSizes[i] > 10 && ClusterSizes[i] <= 20)
	// 		ten2twenty1 += ClusterSizes[i];
	// 	else if (ClusterSizes[i] > 20 && ClusterSizes[i] <= 50)
	// 		twenty2fifty1 += ClusterSizes[i];
	// 	else if (ClusterSizes[i] > 50 && ClusterSizes[i] <= 100)
	// 		fifty2hundred1 += ClusterSizes[i];
	// 	else if (ClusterSizes[i] > 100 && ClusterSizes[i] <= 300)
	// 		hundredplus1 += ClusterSizes[i];
	// 	else if (ClusterSizes[i] > 300 && ClusterSizes[i] <= 500)
	// 		threehundplus1 += ClusterSizes[i];
	// 	if (ClusterSizes[i] > 500)
	// 		fivehundplus1 += ClusterSizes[i];
	// }
	// fprintf(file9, "TS%d,%.2lf,%.2lf,%.2lf,%.2lf,%.2lf,%.2lf,%.2lf,%.2lf,%.2lf,%.2lf\n", timestep, (double)one1 * 100 / numW, (double)two1 * 100 / numW, (double)three2five1 * 100 / numW, (double)five2ten1 * 100 / numW, (double)ten2twenty1 * 100 / numW, (double)twenty2fifty1 * 100 / numW, (double)fifty2hundred1 * 100 / numW, (double)hundredplus1 * 100 / numW, (double)threehundplus1 * 100 / numW, (double)fivehundplus1 * 100 / numW);
	// fprintf(file8, "%.2lf %.2lf %.2lf %.2lf %.2lf %.2lf %.2lf %.2lf %.2lf %.2lf\n", (double)one1 * 100 / numW, (double)two1 * 100 / numW, (double)three2five1 * 100 / numW, (double)five2ten1 * 100 / numW, (double)ten2twenty1 * 100 / numW, (double)twenty2fifty1 * 100 / numW, (double)fifty2hundred1 * 100 / numW, (double)hundredplus1 * 100 / numW, (double)threehundplus1 * 100 / numW, (double)fivehundplus1 * 100 / numW);

	// Printing out results in respective files
	for (i = 0; i < numW; i++)
	{
		fprintf(fHB_water, "%d\n", Total_water_HB[i]);
	}

	////////////////////////////////////
	// Generating a Combined csv file //
	////////////////////////////////////

	// Printing header for combined file
	char *bc0 = "";
	char *bc1 = "Water";
	char *bc2 = "TFSI";
	char *bc3 = "Total";
	char *bc4 = "Li";
	char *bc5 = "Zn";
	char *bc6 = "Total Cation";
	char *bc7 = "F";
	char *bc8 = "N";
	char *bc9 = "O(TFSI)";

	fprintf(fCombined, "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n", bc0, bc1, bc2, bc3, bc4, bc5, bc6, bc7, bc8, bc9);

	// Contains H-bond data for every water with each species
	for (i = 0; i < numW; i++)
	{
		fprintf(fCombined, "%d,", i);
		fprintf(fCombined, "%d,%d,%d,%d,%d,%d,%d,%d,%d\n", Total_water_HB[i], Total_non_water_HB[i], Total_HB[i],
				Cood_with_Li[i], Cood_with_Zn[i], Total_cood_with_cation[i], HB_with_F[i], HB_with_N[i], HB_with_Otfsi[i]);
	}

	//////////////////////////////
	// Printing HBperWater Results //
	//////////////////////////////

	fprintf(fCounter, "%.2lf %.2lf %.2lf %.2lf\n", countW / numW, countA / numW, countL / numW, countZ / numW);

	////////////////////////////////
	// Printing HBpercent results //
	////////////////////////////////
	double Count0 = 0;

	Count0 = 0;
	for (i = 0; i < numW; i++)
	{
		if (Total_water_HB[i] == 0)
		{
			Count0 += 1;
		}
	}
	fprintf(fHBP, "%.2lf   ", Count0 * 100 / numW);

	Count0 = 0;
	for (i = 0; i < numW; i++)
	{
		if (Total_water_HB[i] == 1)
		{
			Count0 += 1;
		}
	}
	fprintf(fHBP, "%.2lf   ", Count0 * 100 / numW);

	Count0 = 0;
	for (i = 0; i < numW; i++)
	{
		if (Total_water_HB[i] == 2)
		{
			Count0 += 1;
		}
	}
	fprintf(fHBP, "%.2lf   ", Count0 * 100 / numW);

	Count0 = 0;
	for (i = 0; i < numW; i++)
	{
		if (Total_water_HB[i] == 3)
		{
			Count0 += 1;
		}
	}
	fprintf(fHBP, "%.2lf   ", Count0 * 100 / numW);

	Count0 = 0;
	for (i = 0; i < numW; i++)
	{
		if (Total_water_HB[i] == 4)
		{
			Count0 += 1;
		}
	}
	fprintf(fHBP, "%.2lf   ", Count0 * 100 / numW);

	fprintf(fHBP, "\n");

	///////////////////////
	// Bifurcated H-Bonds//
	///////////////////////

	for (i = 0; i < numW; i++)
	{
		if ((Donor[i] + Total_non_water_HB[i]) > 2)
		{
			fprintf(fHB5, "%d %d %d %d %d %d\n", timestep, i, Donor[i], HB_with_N[i], HB_with_Otfsi[i], HB_with_F[i]);
		}
	}

	/////////////////////////////
	// Donor-ONLY Probablities //
	////////////////////////////

	double Count1 = 0;
	for (i = 0; i < numW; i++)
	{
		if (Donor[i] == 0)
		{
			Count1 += 1;
		}
	}

	double Count2 = 0;
	for (i = 0; i < numW; i++)
	{
		if ((Donor[i] == 1))
		{
			Count2 += 1;
		}
	}

	double Count3 = 0;
	for (i = 0; i < numW; i++)
	{
		if ((Donor[i] == 2))
		{
			Count3 += 1;
		}
	}

	double Total_Count = Count1 + Count2 + Count3;

	fprintf(fDonor, "%.2lf   %.2lf   %.2lf\n", Count1 / Total_Count, Count2 / Total_Count, Count3 / Total_Count);

	/////////////////////
	// Cation_Ow result /
	/////////////////////

	// Lithium
	double count0 = 0, count1 = 0, count2 = 0, count3 = 0, count4 = 0, count5 = 0, count6 = 0, total_count = 0;
	if (numL != 0)
	{
		for (i = 0; i < numL; i++)
		{
			if (Li_Ow[i] == 0)
				count0 += 1;
			else if (Li_Ow[i] == 1)
				count1 += 1;
			else if (Li_Ow[i] == 2)
				count2 += 1;
			else if (Li_Ow[i] == 3)
				count3 += 1;
			else if (Li_Ow[i] == 4)
				count4 += 1;
			else if (Li_Ow[i] == 5)
				count5 += 1;
			else if (Li_Ow[i] == 6)
				count6 += 1;
		}
		total_count = count0 + count1 + count2 + count3 + count4 + count5 + count6;

		fprintf(fLi_Ow, "%.2lf   %.2lf   %.2lf   %.2lf   %.2lf   %.2lf   %.2lf\n", count0 / total_count, count1 / total_count, count2 / total_count, count3 / total_count, count4 / total_count, count5 / total_count, count6 / total_count);
	}

	// Zinc
	count0 = 0, count1 = 0, count2 = 0, count3 = 0, count4 = 0, count5 = 0, count6 = 0, total_count = 0;
	if (numZ != 0)
	{
		for (i = 0; i < numZ; i++)
		{
			if (Zn_Ow[i] == 0)
				count0 += 1;
			else if (Zn_Ow[i] == 1)
				count1 += 1;
			else if (Zn_Ow[i] == 2)
				count2 += 1;
			else if (Zn_Ow[i] == 3)
				count3 += 1;
			else if (Zn_Ow[i] == 4)
				count4 += 1;
			else if (Zn_Ow[i] == 5)
				count5 += 1;
			else if (Zn_Ow[i] == 6)
				count6 += 1;
		}
		total_count = count0 + count1 + count2 + count3 + count4 + count5 + count6;

		fprintf(fZn_Ow, "%.2lf   %.2lf   %.2lf   %.2lf   %.2lf   %.2lf   %.2lf\n", count0 / total_count, count1 / total_count, count2 / total_count, count3 / total_count, count4 / total_count, count5 / total_count, count6 / total_count);
	}

	//////////////////////
	// Cation_Ot Result //
	//////////////////////

	// Lithium
	count0 = 0, count1 = 0, count2 = 0, count3 = 0, count4 = 0, count5 = 0, count6 = 0, total_count = 0;
	if (numL != 0)
	{
		for (i = 0; i < numL; i++)
		{
			if (Li_Ot[i] == 0)
				count0 += 1;
			else if (Li_Ot[i] == 1)
				count1 += 1;
			else if (Li_Ot[i] == 2)
				count2 += 1;
			else if (Li_Ot[i] == 3)
				count3 += 1;
			else if (Li_Ot[i] == 4)
				count4 += 1;
			else if (Li_Ot[i] == 5)
				count5 += 1;
			else if (Li_Ot[i] == 6)
				count6 += 1;
		}
		total_count = count0 + count1 + count2 + count3 + count4 + count5 + count6;

		fprintf(fLi_Ot, "%.2lf   %.2lf   %.2lf   %.2lf   %.2lf   %.2lf   %.2lf\n", count0 / total_count, count1 / total_count, count2 / total_count, count3 / total_count, count4 / total_count, count5 / total_count, count6 / total_count);
	}

	// Zinc
	count0 = 0, count1 = 0, count2 = 0, count3 = 0, count4 = 0, count5 = 0, count6 = 0, total_count = 0;
	if (numZ != 0)
	{
		for (i = 0; i < numZ; i++)
		{
			if (Zn_Ot[i] == 0)
				count0 += 1;
			else if (Zn_Ot[i] == 1)
				count1 += 1;
			else if (Zn_Ot[i] == 2)
				count2 += 1;
			else if (Zn_Ot[i] == 3)
				count3 += 1;
			else if (Zn_Ot[i] == 4)
				count4 += 1;
			else if (Zn_Ot[i] == 5)
				count5 += 1;
			else if (Zn_Ot[i] == 6)
				count6 += 1;
		}
		total_count = count0 + count1 + count2 + count3 + count4 + count5 + count6;

		fprintf(fZn_Ot, "%.2lf   %.2lf   %.2lf   %.2lf   %.2lf   %.2lf   %.2lf\n", count0 / total_count, count1 / total_count, count2 / total_count, count3 / total_count, count4 / total_count, count5 / total_count, count6 / total_count);
	}

	///////////////////
	// Otfsi H-bonds //
	///////////////////

	count0 = 0, count1 = 0, count2 = 0, count3 = 0, count4 = 0, total_count = 0;

	for (i = 0; i < 4 * numT; i++)
	{
		if (Acceptor_Ot[i] == 0)
			count0 += 1;
		else if (Acceptor_Ot[i] == 1)
			count1 += 1;
		else if (Acceptor_Ot[i] == 2)
			count2 += 1;
	}

	total_count = count0 + count1 + count2;

	fprintf(fH_Ot, "%.2lf   %.2lf   %.2lf\n", count0 / total_count, count1 / total_count, count2 / total_count);

	// H2O-H2O Accepted H-bonds

	count0 = 0, count1 = 0, count2 = 0, count3 = 0, count4 = 0, total_count = 0;

	for (i = 0; i < numW; i++)
	{
		if (Acceptor[i] == 0)
			count0 += 1;
		else if (Acceptor[i] == 1)
			count1 += 1;
		else if (Acceptor[i] == 2)
			count2 += 1;
	}

	total_count = count0 + count1 + count2;

	fprintf(fAccept, "%.2lf   %.2lf   %.2lf\n", count0 / total_count, count1 / total_count, count2 / total_count);

	//////////////////////////////
	// Otfsi-Cation Coordination//
	//////////////////////////////

	// Lithium
	if (numL != 0)
	{
		count0 = 0, count1 = 0, count2 = 0, total_count = 0;

		for (i = 0; i < 4 * numT; i++)
		{
			if (Cood_with_Li_Ot[i] == 0)
				count0 += 1;
			else if (Cood_with_Li_Ot[i] == 1)
				count1 += 1;
			else if (Cood_with_Li_Ot[i] == 2)
				count2 += 1;
		}
		total_count = count0 + count1 + count2;
		fprintf(fOt_Li, "%.2lf   %.2lf   %.2lf\n", count0 / total_count, count1 / total_count, count2 / total_count);
	}

	// Zinc
	if (numZ != 0)
	{
		count0 = 0, count1 = 0, count2 = 0, total_count = 0;

		for (i = 0; i < 4 * numT; i++)
		{
			if (Cood_with_Zn_Ot[i] == 0)
				count0 += 1;
			else if (Cood_with_Zn_Ot[i] == 1)
				count1 += 1;
			else if (Cood_with_Zn_Ot[i] == 2)
				count2 += 1;
		}
		total_count = count0 + count1 + count2;
		fprintf(fOt_Zn, "%.2lf   %.2lf   %.2lf\n", count0 / total_count, count1 / total_count, count2 / total_count);
	}

	// Free the allocated memory
	for (int i = 0; i < 2 * numW; i++)
	{
		free(exchange[i]);
	}
	free(exchange);

	// Free the allocated memory
	// for (int i = 0; i < numW; i++)
	// {
	// 	free(Adj_mat[i]);
	// }
	// free(Adj_mat);

	/////////////////////
	// CLOSING THE FILES//
	/////////////////////

	fclose(fH1);
	fclose(fH2);
	fclose(fOW);

	fclose(fLI);
	fclose(fZN);

	fclose(fI);

	fclose(fN1);
	fclose(fF);
	fclose(fOt);

	fclose(fHB_water);

	fclose(fCombined);

	fclose(fCounter);
	fclose(fHBP);

	fclose(fHB5);

	fclose(fLi_Ow);
	fclose(fZn_Ow);

	fclose(fLi_Ot);
	fclose(fZn_Ot);

	fclose(fH_Ot);
	fclose(fAccept);

	fclose(fOt_Li);
	fclose(fOt_Zn);

	// fclose(file5);
	// fclose(file6);
	// fclose(file7);
	// fclose(file8);
	// fclose(file9);

	fclose(fRtime);
	fclose(fRtimeF);
	fclose(fEdges);
	fclose(fEdgesWT);
	fclose(fOriOw);
	fclose(fOriOt);
	fclose(fOriF);
	fclose(fOriN);

	return 0;
}
