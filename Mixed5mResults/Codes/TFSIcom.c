// Centre of Mass Coordinates of TFSI printed

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

double box_length;
int numW, numT;

int main()
{

	FILE *fC1 = fopen("Coordinates/C1.txt", "r");
	FILE *fC2 = fopen("Coordinates/C2.txt", "r");
	FILE *fN1 = fopen("Coordinates/N1.txt", "r");
	FILE *fF1 = fopen("Coordinates/F1.txt", "r");
	FILE *fF2 = fopen("Coordinates/F2.txt", "r");
	FILE *fF3 = fopen("Coordinates/F3.txt", "r");
	FILE *fF4 = fopen("Coordinates/F4.txt", "r");
	FILE *fF5 = fopen("Coordinates/F5.txt", "r");
	FILE *fF6 = fopen("Coordinates/F6.txt", "r");
	FILE *fS1 = fopen("Coordinates/S1.txt", "r");
	FILE *fS2 = fopen("Coordinates/S2.txt", "r");
	FILE *fO1 = fopen("Coordinates/O1.txt", "r");
	FILE *fO2 = fopen("Coordinates/O2.txt", "r");
	FILE *fO3 = fopen("Coordinates/O3.txt", "r");
	FILE *fO4 = fopen("Coordinates/O4.txt", "r");

	FILE *fI = fopen("SysInfo.txt", "r");

	FILE *fCOM = fopen("Plots/DataFiles/TFSIcom.txt", "w");

	fscanf(fI, "%lf %d %d", &box_length, &numW, &numT);

	int i = 0, j = 0, k = 0, m = 0;

	double C1[numT][3], C2[numT][3], N1[numT][3], F1[numT][3], F2[numT][3], F3[numT][3], F4[numT][3], F5[numT][3], F6[numT][3], S1[numT][3], S2[numT][3], O1[numT][3], O2[numT][3], O3[numT][3], O4[numT][3];

	// Scanning the coordinates
	for (i = 0; i < numT; i++)
	{
		for (j = 0; j < 3; j++)
		{
			fscanf(fC1, "%lf", &C1[i][j]);
		}
	}

	for (i = 0; i < numT; i++)
	{
		for (j = 0; j < 3; j++)
		{
			fscanf(fC2, "%lf", &C2[i][j]);
		}
	}

	for (i = 0; i < numT; i++)
	{
		for (j = 0; j < 3; j++)
		{
			fscanf(fN1, "%lf", &N1[i][j]);
		}
	}

	for (i = 0; i < numT; i++)
	{
		for (j = 0; j < 3; j++)
		{
			fscanf(fF1, "%lf", &F1[i][j]);
		}
	}
	for (i = 0; i < numT; i++)
	{
		for (j = 0; j < 3; j++)
		{
			fscanf(fF2, "%lf", &F2[i][j]);
		}
	}
	for (i = 0; i < numT; i++)
	{
		for (j = 0; j < 3; j++)
		{
			fscanf(fF3, "%lf", &F3[i][j]);
		}
	}
	for (i = 0; i < numT; i++)
	{
		for (j = 0; j < 3; j++)
		{
			fscanf(fF4, "%lf", &F4[i][j]);
		}
	}
	for (i = 0; i < numT; i++)
	{
		for (j = 0; j < 3; j++)
		{
			fscanf(fF5, "%lf", &F5[i][j]);
		}
	}
	for (i = 0; i < numT; i++)
	{
		for (j = 0; j < 3; j++)
		{
			fscanf(fF6, "%lf", &F6[i][j]);
		}
	}

	for (i = 0; i < numT; i++)
	{
		for (j = 0; j < 3; j++)
		{
			fscanf(fS1, "%lf", &S1[i][j]);
		}
	}
	for (i = 0; i < numT; i++)
	{
		for (j = 0; j < 3; j++)
		{
			fscanf(fS2, "%lf", &S2[i][j]);
		}
	}

	for (i = 0; i < numT; i++)
	{
		for (j = 0; j < 3; j++)
		{
			fscanf(fO1, "%lf", &O1[i][j]);
		}
	}
	for (i = 0; i < numT; i++)
	{
		for (j = 0; j < 3; j++)
		{
			fscanf(fO2, "%lf", &O2[i][j]);
		}
	}
	for (i = 0; i < numT; i++)
	{
		for (j = 0; j < 3; j++)
		{
			fscanf(fO3, "%lf", &O3[i][j]);
		}
	}
	for (i = 0; i < numT; i++)
	{
		for (j = 0; j < 3; j++)
		{
			fscanf(fO4, "%lf", &O4[i][j]);
		}
	}

	double Xcom = 0, Ycom = 0, Zcom = 0;

	double Mc = 12, Mf = 19, Mo = 16, Ms = 32, Mn = 14;

	for (i = 0; i < numT; i++)
	{

		Xcom = (Mn * N1[i][0] + Mc * C1[i][0] + Mc * C2[i][0] + Ms * S1[i][0] + Ms * S2[i][0] + Mo * O1[i][0] + Mo * O2[i][0] + Mo * O3[i][0] + Mo * O4[i][0] + Mf * F1[i][0] + Mf * F2[i][0] + Mf * F3[i][0] + Mf * F4[i][0] + Mf * F5[i][0] + Mf * F6[i][0]) / (Mn + 2 * Mc + 6 * Mf + 2 * Ms + 4 * Mo);
		fprintf(fCOM, "%lf ", Xcom);

		Ycom = (Mn * N1[i][1] + Mc * C1[i][1] + Mc * C2[i][1] + Ms * S1[i][1] + Ms * S2[i][1] + Mo * O1[i][1] + Mo * O2[i][1] + Mo * O3[i][1] + Mo * O4[i][1] + Mf * F1[i][1] + Mf * F2[i][1] + Mf * F3[i][1] + Mf * F4[i][1] + Mf * F5[i][1] + Mf * F6[i][1]) / (Mn + 2 * Mc + 6 * Mf + 2 * Ms + 4 * Mo);
		fprintf(fCOM, "%lf ", Ycom);

		Zcom = (Mn * N1[i][2] + Mc * C1[i][2] + Mc * C2[i][2] + Ms * S1[i][2] + Ms * S2[i][2] + Mo * O1[i][2] + Mo * O2[i][2] + Mo * O3[i][2] + Mo * O4[i][2] + Mf * F1[i][2] + Mf * F2[i][2] + Mf * F3[i][2] + Mf * F4[i][2] + Mf * F5[i][2] + Mf * F6[i][2]) / (Mn + 2 * Mc + 6 * Mf + 2 * Ms + 4 * Mo);
		fprintf(fCOM, "%lf\n", Zcom);
	}

	fclose(fI);
	fclose(fCOM);
	fclose(fC1);
	fclose(fC2);
	fclose(fS1);
	fclose(fS2);
	fclose(fO1);
	fclose(fO2);
	fclose(fO3);
	fclose(fO4);
	fclose(fN1);
	fclose(fF1);
	fclose(fF2);
	fclose(fF3);
	fclose(fF4);
	fclose(fF5);
	fclose(fF6);

	return 0;
}