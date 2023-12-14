// Header values added to csv file

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

/*

To Run:
gcc Codes/Header.c -o Codes/Header -lm
./Codes/Header
rm Codes/Header

*/

int main()
{

	FILE *file1 = fopen("Water/FreqDistri.csv", "w");
	FILE *file2 = fopen("Averages/AvgFreqDistri1.csv", "w");

	FILE *file3 = fopen("Averages/AvgHBpercent.csv", "w");
	FILE *file4 = fopen("Averages/AvgHBperWater.csv", "w");
	FILE *file5 = fopen("Plots/DataFiles/Combined.csv", "w");

	FILE *file6 = fopen("Water/ClusterPercent1.csv", "w");
	FILE *file7 = fopen("Averages/AvgClusterPercent1.csv", "w");

	FILE *file8 = fopen("Averages/AvgHBabsValues.csv", "w");
	FILE *file9 = fopen("Plots/DataFiles/HBabsValues.txt", "w");

	FILE *file10 = fopen("Water/BifurcatedHB.txt", "w");

	FILE *file11 = fopen("Averages/AvgIsoHBabsValues.csv", "w");
	FILE *file12 = fopen("Plots/DataFiles/IsoHBabsValues.txt", "w");

	FILE *file13 = fopen("Averages/AvgOneHBabsValues.csv", "w");
	FILE *file14 = fopen("Plots/DataFiles/OneHBabsValues.txt", "w");

	FILE *file15 = fopen("Averages/AvgDonorProbs.csv", "w");
	FILE *file16 = fopen("Averages/AvgAcceptorProbsOw.csv", "w");

	FILE *file17 = fopen("Averages/AvgLiCoodOw.csv", "w");
	FILE *file18 = fopen("Averages/AvgLiCoodOt.csv", "w");

	FILE *file19 = fopen("Averages/AvgZnCoodOw.csv", "w");
	FILE *file20 = fopen("Averages/AvgZnCoodOt.csv", "w");

	FILE *file21 = fopen("Averages/AvgAcceptorProbsOt.csv", "w");
	FILE *file22 = fopen("Averages/AvgOtCoodLi.csv", "w");
	FILE *file23 = fopen("Averages/AvgOtCoodZn.csv", "w");

	FILE *file29 = fopen("RT/Timestep.txt", "w");
	FILE *file30 = fopen("Averages/AvgDonorLifetimes.txt", "w");

	fprintf(file29, "1 1 1 1 1 1 1");

	fprintf(file30, "0w,1w,2w\n");

	char *t = "Timestep";
	char *s = "Size";
	char *one = "1";
	char *a = "2";
	char *b = "3-5";
	char *c = "6-10";
	char *d = "11-20";
	char *e = "21-50";
	char *f = "51-100";
	char *g = "101-300";
	char *h = "301-500";
	char *i = ">500";
	char *avg = "Avg";

	fprintf(file1, "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n", t, one, a, b, c, d, e, f, g, h, i);

	fprintf(file2, "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n%s,", s, one, a, b, c, d, e, f, g, h, i, avg);

	fprintf(file7, "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n%s,", s, one, a, b, c, d, e, f, g, h, i, avg);

	fprintf(file6, "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n", t, one, a, b, c, d, e, f, g, h, i);

	char *s1 = "#HB";
	char *a1 = "Zero";
	char *b1 = "One";
	char *c1 = "Two";
	char *d1 = "Three";
	char *e1 = "Four";
	char *avg1 = "Avg";

	fprintf(file3, "%s,%s,%s,%s,%s,%s\n%s,", s1, a1, b1, c1, d1, e1, avg1);

	char *s2 = "#HB/H2O";
	char *a2 = "H2O";
	char *b2 = "TFSI";
	char *c2 = "Li";
	char *d2 = "Zn";
	char *avg2 = "Avg";

	fprintf(file4, "%s,%s,%s,%s,%s\n%s,", s2, a2, b2, c2, d2, avg2);

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

	fprintf(file5, "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n", bc0, bc1, bc2, bc3, bc4, bc5, bc6, bc7, bc8, bc9);

	char *aa1 = "Zero";
	char *aa2 = "One";
	char *aa3 = "Two";

	fprintf(file15, "%s,%s,%s\n", aa1, aa2, aa3);

	fprintf(file16, "%s,%s,%s\n", aa1, aa2, aa3);

	fprintf(file21, "%s,%s,%s\n", aa1, aa2, aa3);

	fprintf(file22, "%s,%s,%s\n", aa1, aa2, aa3);

	fprintf(file23, "%s,%s,%s\n", aa1, aa2, aa3);

	char *a3 = "Zero";
	char *b3 = "One";
	char *c3 = "Two";
	char *d3 = "Three";
	char *e3 = "Four";
	char *f3 = "Five";
	char *g3 = "Six";

	fprintf(file17, "%s,%s,%s,%s,%s,%s,%s\n", a3, b3, c3, d3, e3, f3, g3);
	fprintf(file18, "%s,%s,%s,%s,%s,%s,%s\n", a3, b3, c3, d3, e3, f3, g3);
	fprintf(file19, "%s,%s,%s,%s,%s,%s,%s\n", a3, b3, c3, d3, e3, f3, g3);
	fprintf(file20, "%s,%s,%s,%s,%s,%s,%s\n", a3, b3, c3, d3, e3, f3, g3);

	char *t5 = "Timestep";
	char *s5 = "Size";
	char *one5 = "1";
	char *a5 = "2-5";
	char *b5 = "6-10";
	char *c5 = "11-20";
	char *d5 = "21-40";
	char *e5 = "41-60";
	char *f5 = "61-80";
	char *g5 = "81-100";
	char *h5 = "101-120";
	char *i5 = "121-150";
	char *avg5 = "Avg";

	fclose(file1);
	fclose(file2);
	fclose(file3);
	fclose(file4);
	fclose(file5);
	fclose(file6);
	fclose(file7);
	fclose(file8);
	fclose(file9);
	fclose(file10);
	fclose(file11);
	fclose(file12);
	fclose(file13);
	fclose(file14);
	fclose(file15);
	fclose(file16);
	fclose(file17);
	fclose(file18);
	fclose(file19);
	fclose(file20);
	fclose(file21);
	fclose(file22);
	fclose(file23);
	fclose(file29);
	fclose(file30);

	return 0;
}