// Averaging values over timesteps

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

int main()
{

	int m, n;
	int i = 0, j = 0;

	// n=9;
	// n=5;
	// m=10;
	// m=2001;
	// m=10001;

	scanf("%d %d", &m, &n);

	double Array[m][n];

	char *fname1 = malloc(100), *fname2 = malloc(100);
	scanf("%s %s", fname1, fname2);

	// fname1="ClusterSize/FreqDistriA.txt";
	// fname2="Averages/AvgFreqDistri.csv";

	FILE *file1 = fopen(fname1, "r");
	FILE *file2 = fopen(fname2, "a");

	double avg[n];

	for (i = 0; i < n; i++)
	{
		avg[i] = 0;
	}

	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			fscanf(file1, "%lf", &Array[i][j]);
		}
	}

	for (j = 0; j < n; j++)
	{
		for (i = 0; i < m; i++)
		{
			avg[j] += Array[i][j];
		}
	}

	for (j = 0; j < n - 1; j++)
	{
		fprintf(file2, "%.2lf,", avg[j] / m);
	}

	fprintf(file2, "%.2lf", avg[n - 1] / m);

	fprintf(file2, "\n");

	// free(fname1);

	fclose(file1);
	fclose(file2);

	return 0;
}