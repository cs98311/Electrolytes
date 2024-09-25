// Produces vectors for all H bonds from Oxy to Oxy

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

double box_length;
int numW;
int count = 0;

double get_magnitude(double m[3])
{
	double mag = 0;
	mag = sqrt(m[0] * m[0] + m[1] * m[1] + m[2] * m[2]);
	return mag;
}

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

double vector_angle(double v1[3], double v2[3])
{

	double dot = 0, mag1 = 0, mag2 = 0;
	dot = dot_product(v1, v2);
	mag1 = get_magnitude(v1);
	mag2 = get_magnitude(v2);
	return acos((dot) / (mag1 * mag2));
}

void swap(double *xp, double *yp)
{
	double temp = *xp;
	*xp = *yp;
	*yp = temp;
}

void selectionSort(double arr[], int n)
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

void Check_HB_H1(double H1[numW][3], double Ox[numW][3], int m, FILE *file4)
{

	int i = 0, j = 0, k = 0, p = 0;
	double distance = 0, v1[3], v2[3], angle = 0;
	double dL[3];

	for (j = 0; j < 3; j++)
	{
		v1[j] = H1[m][j] - Ox[m][j];
	}

	for (i = 0; i < numW; i++)
	{

		if (i == m)
		{
			continue;
		}

		for (j = 0; j < 3; j++)
		{
			dL[j] = (Ox[i][j] - H1[m][j]);

			if (dL[j] > box_length / 2)
			{
				dL[j] -= box_length;
			}
			else if (dL[j] < -box_length / 2)
			{
				dL[j] += box_length;
			}
		}

		double mirrvec[3];
		for (j = 0; j < 3; j++)
		{
			mirrvec[j] = Ox[i][j] - Ox[m][j];
			if (mirrvec[j] > box_length / 2)
				mirrvec[j] -= box_length;
			else if (mirrvec[j] < -box_length / 2)
				mirrvec[j] += box_length;
		}

		for (j = 0; j < 3; j++)
		{
			v2[j] = v1[j] + dL[j];
		}
		// fprintf(file4,"old v2: ");
		double oldv2[3];
		for (j = 0; j < 3; j++)
		{
			oldv2[j] = Ox[i][j] - Ox[m][j];
		}

		distance = sqrt(dL[0] * dL[0] + dL[1] * dL[1] + dL[2] * dL[2]);

		angle = vector_angle(v1, v2);
		angle *= 57.2958;

		if (distance < 0.22 && angle < 35)
		{

			for (j = 0; j < 3; j++)
			{
				fprintf(file4, "%lf,", Ox[m][j]);
			}
			for (j = 0; j < 2; j++)
			{
				fprintf(file4, "%lf,", mirrvec[j]);
			}
			fprintf(file4, "%lf", mirrvec[2]);

			fprintf(file4, "\n");
		}
	}
}

int main()
{

	FILE *file0 = fopen("SysInfo.txt", "r");
	FILE *file1 = fopen("Coordinates/H1.txt", "r");
	FILE *file2 = fopen("Coordinates/H2.txt", "r");
	FILE *file3 = fopen("Coordinates/OW.txt", "r");
	FILE *file4 = fopen("Plots/Datafiles/HBvectorsO2O.txt", "w");

	fscanf(file0, "%lf %d", &box_length, &numW);

	int i = 0, j = 0, k = 0, m = 0, stored = 0;
	double H1[numW][3], H2[numW][3], Ow[numW][3];

	for (i = 0; i < numW; i++)
	{
		for (j = 0; j < 3; j++)
		{
			fscanf(file1, "%lf", &H1[i][j]);
		}
	}

	for (i = 0; i < numW; i++)
	{
		for (j = 0; j < 3; j++)
		{
			fscanf(file2, "%lf", &H2[i][j]);
		}
	}

	for (i = 0; i < numW; i++)
	{
		for (j = 0; j < 3; j++)
		{
			fscanf(file3, "%lf", &Ow[i][j]);
		}
	}

	for (m = 0; m < numW; m++)
	{

		Check_HB_H1(H1, Ow, m, file4);
		Check_HB_H1(H2, Ow, m, file4);
	}

	fclose(file0);
	fclose(file1);
	fclose(file2);
	fclose(file3);
	fclose(file4);

	return 0;
}