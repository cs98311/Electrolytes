#include <iostream>
#include <vector>
#include <array>
#include <fstream>
#include <string>
#include <cmath>

// declare global variables
double BOX_LENGTH{0};
int NUM_H2O{0}, NUM_LI{0}, NUM_TFSI{0}, NUM_ZN{0};
const double H_BOND_DISTANCE_CUTOFF{0.25}, H_BOND_ANGLE_CUTOFF{35};
const double SS1_CUTOFF{0.28}, SS2_CUTOFF{0.54};

// declare functions
void readSystemInfo();
void readCoordinates(const std::string &fileName, double coordinates[][3], int numCoordinates);
double calculateDistance(const double point1[3], const double point2[3]);
void PBC(double vector[3]);
double vector_magnitude(double m[3]);
double dot_product(double m[3], double n[3]);
double angle_between_vectors(double v1[3], double v2[3]);

int main()
{

    // intake the system info
    readSystemInfo();

    // declare arrays to hold coordinates of each element
    double Li[NUM_LI][3], Zn[NUM_ZN][3];
    double H1[NUM_H2O][3], H2[NUM_H2O][3], Ow[NUM_H2O][3];
    double Ot[NUM_TFSI * 4][3], F[NUM_TFSI * 6][3], N[NUM_TFSI][3];

    // input the coordinates
    readCoordinates("Coordinates/LI.txt", Li, NUM_LI);
    readCoordinates("Coordinates/ZN.txt", Zn, NUM_ZN);
    readCoordinates("Coordinates/H1.txt", H1, NUM_H2O);
    readCoordinates("Coordinates/H2.txt", H2, NUM_H2O);
    readCoordinates("Coordinates/OW.txt", Ow, NUM_H2O);
    readCoordinates("Coordinates/F.txt", F, 6 * NUM_TFSI);
    readCoordinates("Coordinates/Ot.txt", Ot, 4 * NUM_TFSI);
    readCoordinates("Coordinates/N1.txt", N, NUM_TFSI);

    // working with m-th Li+. let m=0
    // int m{0};
    for (auto m = 0; m < NUM_LI; m++)
    {

        // check distance from every Ow
        std::vector<int> Ow_SS1;
        std::vector<int> Ow_SS2;
        std::vector<int> Ow_outside;
        for (auto i = 0; i < NUM_H2O; i++)
        {
            double dL[3];
            for (auto k = 0; k < 3; k++)
            {
                dL[k] = Ow[i][k] - Li[m][k];
            }
            PBC(dL);
            double distance = vector_magnitude(dL);

            // check if present in SS1
            if (distance < SS1_CUTOFF)
            {
                Ow_SS1.push_back(i);
            }

            // check if present in SS2
            else if (distance > SS1_CUTOFF && distance < SS2_CUTOFF)
            {
                Ow_SS2.push_back(i);
            }

            else
            {
                Ow_outside.push_back(i);
            }
        }

        // std::vector<int> TFSI_SS1;
        // std::vector<int> TFSI_SS2;
        // std::vector<int> TFSI_outside;

        // check distance from every Ot
        std::vector<int> Ot_SS1;
        std::vector<int> Ot_SS2;
        std::vector<int> Ot_outside;
        for (auto i = 0; i < NUM_TFSI * 4; i++)
        {
            // int tfsi_number = i % 4;
            double dL[3];
            for (auto k = 0; k < 3; k++)
            {
                dL[k] = Ot[i][k] - Li[m][k];
            }
            PBC(dL);
            double distance = vector_magnitude(dL);

            // check if present in SS1
            if (distance < SS1_CUTOFF)
            {
                Ot_SS1.push_back(i);
            }

            // check if present in SS2
            else if (distance > SS1_CUTOFF && distance < SS2_CUTOFF)
            {
                Ot_SS2.push_back(i);
            }
            else
            {
                Ot_outside.push_back(i);
            }
        }

        // check distance from every F
        std::vector<int> F_SS1;
        std::vector<int> F_SS2;
        std::vector<int> F_outside;
        for (auto i = 0; i < NUM_TFSI * 6; i++)
        {
            // int tfsi_number = i % 6;
            double dL[3];
            for (auto k = 0; k < 3; k++)
            {
                dL[k] = F[i][k] - Li[m][k];
            }
            PBC(dL);
            double distance = vector_magnitude(dL);

            // check if present in SS1
            if (distance < SS1_CUTOFF)
            {
                F_SS1.push_back(i);
            }

            // check if present in SS2
            else if (distance > SS1_CUTOFF && distance < SS2_CUTOFF)
            {
                F_SS2.push_back(i);
            }
            else
            {
                F_outside.push_back(i);
            }
        }

        // check distance from every N
        std::vector<int> N_SS1;
        std::vector<int> N_SS2;
        std::vector<int> N_outside;
        for (auto i = 0; i < NUM_TFSI; i++)
        {
            // int tfsi_number = i % 6;
            double dL[3];
            for (auto k = 0; k < 3; k++)
            {
                dL[k] = N[i][k] - Li[m][k];
            }
            PBC(dL);
            double distance = vector_magnitude(dL);

            // check if present in SS1
            if (distance < SS1_CUTOFF)
            {
                N_SS1.push_back(i);
            }

            // check if present in SS2
            else if (distance > SS1_CUTOFF && distance < SS2_CUTOFF)
            {
                N_SS2.push_back(i);
            }
            else
            {
                N_outside.push_back(i);
            }
        }

        // SS1 composition
        std::ofstream comp1Li{"Cation/LiSS1comp.txt", std::ios::app};
        if (!comp1Li)
        {
            std::cerr << "Uh oh, SS1compLi.txt could not be opened for writing!\n";
            return 1;
        }
        comp1Li << Ow_SS1.size() << " " << Ot_SS1.size() << " " << F_SS1.size() << " " << N_SS1.size() << "\n";
        comp1Li.close();

        // SS2 composition
        std::ofstream comp2Li{"Cation/LiSS2comp.txt", std::ios::app};
        if (!comp2Li)
        {
            std::cerr << "Uh oh, SS2comp.txt could not be opened for writing!\n";
            return 1;
        }
        comp2Li << Ow_SS2.size() << " " << Ot_SS2.size() << " " << F_SS2.size() << " " << N_SS2.size() << "\n";
        comp2Li.close();

        // bonding of SS1 water
        int SS1_data[Ow_SS1.size()][9];
        for (auto i = 0; i < Ow_SS1.size(); ++i)
        {
            for (auto j = 0; j < 9; ++j)
            {
                SS1_data[i][j] = 0;
            }
        }
        // fill 1st column of SS1 data with water tag numbers
        for (auto i = 0; i < Ow_SS1.size(); i++)
        {
            SS1_data[i][0] = Ow_SS1[i];
        }

        int SS2_data[Ow_SS2.size()][13];
        for (auto i = 0; i < Ow_SS2.size(); ++i)
        {
            for (auto j = 0; j < 13; ++j)
            {
                SS2_data[i][j] = 0;
            }
        }
        // fill 1st column of SS2 data with water tag numbers
        for (auto i = 0; i < Ow_SS2.size(); i++)
        {
            SS2_data[i][0] = Ow_SS2[i];
        }

        // 1) SS1 H1 --> SS2 Ow
        int i_counter = 0;
        for (const auto &i : Ow_SS1)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H1[i][k] - Ow[i][k];
            }

            int j_counter = 0;
            for (const auto &j : Ow_SS2)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = Ow[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS1_data[i_counter][1]++;
                    SS2_data[j_counter][2]++;
                }
                j_counter++;
            }
            i_counter++;
        }

        // 2) SS1 H2 --> SS2 Ow
        i_counter = 0;
        for (const auto &i : Ow_SS1)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H2[i][k] - Ow[i][k];
            }
            PBC(v1);

            int j_counter = 0;
            for (const auto &j : Ow_SS2)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = Ow[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS1_data[i_counter][1]++;
                    SS2_data[j_counter][2]++;
                }
                j_counter++;
            }
            i_counter++;
        }

        // 3) SS2 H1 --> SS1 Ow
        i_counter = 0;
        for (const auto &i : Ow_SS2)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H1[i][k] - Ow[i][k];
            }
            PBC(v1);

            int j_counter = 0;
            for (const auto &j : Ow_SS1)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = Ow[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS1_data[j_counter][1]++;
                    SS2_data[i_counter][2]++;
                }
                j_counter++;
            }
            i_counter++;
        }

        // 4) SS2 H2 --> SS1 Ow
        i_counter = 0;
        for (const auto &i : Ow_SS2)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H2[i][k] - Ow[i][k];
            }
            PBC(v1);

            int j_counter = 0;
            for (const auto &j : Ow_SS1)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = Ow[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS1_data[j_counter][1]++;
                    SS2_data[i_counter][2]++;
                }
                j_counter++;
            }
            i_counter++;
        }

        // 5) SS1 H --> SS2 TFSI (Ot,F,N)
        // H1 --> Ot
        i_counter = 0;
        for (const auto &i : Ow_SS1)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H1[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : Ot_SS2)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = Ot[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS1_data[i_counter][3]++;
                }
            }
            i_counter++;
        }

        // H2 --> Ot
        i_counter = 0;
        for (const auto &i : Ow_SS1)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H2[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : Ot_SS2)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = Ot[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS1_data[i_counter][3]++;
                }
            }
            i_counter++;
        }

        // H1 --> F
        i_counter = 0;
        for (const auto &i : Ow_SS1)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H1[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : F_SS2)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = F[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS1_data[i_counter][4]++;
                }
            }
            i_counter++;
        }

        // H2 --> F
        i_counter = 0;
        for (const auto &i : Ow_SS1)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H2[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : F_SS2)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = F[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS1_data[i_counter][4]++;
                }
            }
            i_counter++;
        }

        // H1 --> N
        i_counter = 0;
        for (const auto &i : Ow_SS1)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H1[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : N_SS2)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = N[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS1_data[i_counter][5]++;
                }
            }
            i_counter++;
        }

        // H2 --> N
        i_counter = 0;
        for (const auto &i : Ow_SS1)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H2[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : N_SS2)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = N[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS1_data[i_counter][5]++;
                }
            }
            i_counter++;
        }

        // 6) SS1 H --> SS1 TFSI (Ot,F,N)
        // H1 --> Ot
        i_counter = 0;
        for (const auto &i : Ow_SS1)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H1[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : Ot_SS1)
            {
                if (i == j)
                    continue;
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = Ot[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS1_data[i_counter][6]++;
                }
            }
            i_counter++;
        }

        // H2 --> Ot
        i_counter = 0;
        for (const auto &i : Ow_SS1)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H2[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : Ot_SS1)
            {
                if (i == j)
                    continue;
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = Ot[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS1_data[i_counter][6]++;
                }
            }
            i_counter++;
        }

        // H1 --> F
        i_counter = 0;
        for (const auto &i : Ow_SS1)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H1[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : F_SS1)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = F[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS1_data[i_counter][7]++;
                }
            }
            i_counter++;
        }

        // H2 --> F
        i_counter = 0;
        for (const auto &i : Ow_SS1)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H2[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : F_SS1)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = F[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS1_data[i_counter][7]++;
                }
            }
            i_counter++;
        }

        // H1 --> N
        i_counter = 0;
        for (const auto &i : Ow_SS1)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H1[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : N_SS1)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = N[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS1_data[i_counter][8]++;
                }
            }
            i_counter++;
        }

        // H2 --> N
        i_counter = 0;
        for (const auto &i : Ow_SS1)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H2[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : N_SS1)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = N[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS1_data[i_counter][8]++;
                }
            }
            i_counter++;
        }

        // 7) SS1 H1 --> SS1 Ow
        i_counter = 0;
        for (const auto &i : Ow_SS1)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H1[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : Ow_SS1)
            {
                if (i == j)
                    continue;
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = Ow[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS1_data[i_counter][2]++;
                }
            }
            i_counter++;
        }

        // 8) SS1 H2 --> SS1 Ow
        i_counter = 0;
        for (const auto &i : Ow_SS1)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H2[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : Ow_SS1)
            {
                if (i == j)
                    continue;
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = Ow[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS1_data[i_counter][2]++;
                }
            }
            i_counter++;
        }

        // bonding of SS2 water

        // 1) SS2 Ow --> SS2 Ow
        // H1 --> Ow
        i_counter = 0;
        for (const auto &i : Ow_SS2)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H1[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : Ow_SS2)
            {
                if (i == j)
                    continue;

                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = Ow[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS2_data[i_counter][1]++;
                }
            }
            i_counter++;
        }

        // H2 --> Ow
        i_counter = 0;
        for (const auto &i : Ow_SS2)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H2[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : Ow_SS2)
            {
                if (i == j)
                    continue;

                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = Ow[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS2_data[i_counter][1]++;
                }
            }
            i_counter++;
        }

        // 2) SS2 Ow --> SS2 TFSI
        // H1 --> Ot
        i_counter = 0;
        for (const auto &i : Ow_SS2)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H1[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : Ot_SS2)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = Ot[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS2_data[i_counter][4]++;
                }
            }
            i_counter++;
        }

        // H2 --> Ot
        i_counter = 0;
        for (const auto &i : Ow_SS2)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H2[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : Ot_SS2)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = Ot[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS2_data[i_counter][4]++;
                }
            }
            i_counter++;
        }

        // H1 --> F
        i_counter = 0;
        for (const auto &i : Ow_SS2)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H1[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : F_SS2)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = F[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS2_data[i_counter][5]++;
                }
            }
            i_counter++;
        }

        // H2 --> F
        i_counter = 0;
        for (const auto &i : Ow_SS2)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H2[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : F_SS2)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = F[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS2_data[i_counter][5]++;
                }
            }
            i_counter++;
        }

        // H1 --> N
        i_counter = 0;
        for (const auto &i : Ow_SS2)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H1[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : N_SS2)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = N[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS2_data[i_counter][6]++;
                }
            }
            i_counter++;
        }

        // H2 --> N
        i_counter = 0;
        for (const auto &i : Ow_SS2)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H2[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : N_SS2)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = N[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS2_data[i_counter][6]++;
                }
            }
            i_counter++;
        }

        // 3) SS2 Ow --> Outside Ow
        // SS2 H1 --> Outside Ow
        i_counter = 0;
        for (const auto &i : Ow_SS2)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H1[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : Ow_outside)
            {

                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = Ow[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS2_data[i_counter][3]++;
                }
            }
            i_counter++;
        }

        // SS2 H2 --> Outside Ow
        i_counter = 0;
        for (const auto &i : Ow_SS2)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H2[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : Ow_outside)
            {

                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = Ow[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS2_data[i_counter][3]++;
                }
            }
            i_counter++;
        }

        // outside H1 --> SS2 Ow
        i_counter = 0;
        for (const auto &i : Ow_outside)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H1[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : Ow_SS2)
            {

                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = Ow[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS2_data[i_counter][3]++;
                }
            }
            i_counter++;
        }

        // Outside H2 --> SS2 Ow
        i_counter = 0;
        for (const auto &i : Ow_outside)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H2[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : Ow_SS2)
            {

                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = Ow[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS2_data[i_counter][3]++;
                }
            }
            i_counter++;
        }

        // 4) SS2 Ow --> Outside TFSI (Ot,N,F)
        // H1 --> Ot
        i_counter = 0;
        for (const auto &i : Ow_SS2)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H1[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : Ot_outside)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = Ot[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS2_data[i_counter][10]++;
                }
            }
            i_counter++;
        }

        // H2 --> Ot
        i_counter = 0;
        for (const auto &i : Ow_SS2)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H2[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : Ot_outside)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = Ot[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS2_data[i_counter][10]++;
                }
            }
            i_counter++;
        }

        // H1 --> F
        i_counter = 0;
        for (const auto &i : Ow_SS2)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H1[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : F_outside)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = F[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS2_data[i_counter][11]++;
                }
            }
            i_counter++;
        }

        // H2 --> F
        i_counter = 0;
        for (const auto &i : Ow_SS2)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H2[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : F_outside)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = F[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS2_data[i_counter][11]++;
                }
            }
            i_counter++;
        }

        // H1 --> N
        i_counter = 0;
        for (const auto &i : Ow_SS2)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H1[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : N_outside)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = N[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS2_data[i_counter][12]++;
                }
            }
            i_counter++;
        }

        // H2 --> N
        i_counter = 0;
        for (const auto &i : Ow_SS2)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H2[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : N_outside)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = N[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS2_data[i_counter][12]++;
                }
            }
            i_counter++;
        }

        // 5) SS2 Ow --> SS1 TFSI
        // H1 --> Ot
        i_counter = 0;
        for (const auto &i : Ow_SS2)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H1[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : Ot_SS1)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = Ot[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS2_data[i_counter][7]++;
                }
            }
            i_counter++;
        }

        // H2 --> Ot
        i_counter = 0;
        for (const auto &i : Ow_SS2)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H2[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : Ot_SS1)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = Ot[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS2_data[i_counter][7]++;
                }
            }
            i_counter++;
        }

        // H1 --> F
        i_counter = 0;
        for (const auto &i : Ow_SS2)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H1[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : F_SS1)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = F[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS2_data[i_counter][8]++;
                }
            }
            i_counter++;
        }

        // H2 --> F
        i_counter = 0;
        for (const auto &i : Ow_SS2)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H2[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : F_SS1)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = F[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS2_data[i_counter][8]++;
                }
            }
            i_counter++;
        }

        // H1 --> N
        i_counter = 0;
        for (const auto &i : Ow_SS2)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H1[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : N_SS1)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = N[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS2_data[i_counter][9]++;
                }
            }
            i_counter++;
        }

        // H2 --> N
        i_counter = 0;
        for (const auto &i : Ow_SS2)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H2[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : N_SS1)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = N[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS2_data[i_counter][9]++;
                }
            }
            i_counter++;
        }

        // SS1 bonding output
        std::ofstream bond1Li{"Cation/LiSS1bond.txt", std::ios::app};
        if (!bond1Li)
        {
            std::cerr << "Uh oh, SS1bond.txt could not be opened for writing!\n";
            return 1;
        }
        for (auto i = 0; i < Ow_SS1.size(); i++)
        {
            for (auto j = 1; j < 9; j++)
            {
                bond1Li << SS1_data[i][j] << " ";
            }
            bond1Li << "\n";
        }
        bond1Li.close();

        // SS2 bonding output
        std::ofstream bond2Li{"Cation/LiSS2bond.txt", std::ios::app};
        if (!bond2Li)
        {
            std::cerr << "Uh oh, SS2bond.txt could not be opened for writing!\n";
            return 1;
        }
        for (auto i = 0; i < Ow_SS2.size(); i++)
        {
            for (auto j = 1; j < 13; j++)
            {
                bond2Li << SS2_data[i][j] << " ";
            }
            bond2Li << "\n";
        }
        bond2Li.close();
    }

    // working with m-th Zn+2. let m=0
    // int m{0};
    for (auto m = 0; m < NUM_ZN; m++)
    {

        // check distance from every Ow
        std::vector<int> Ow_SS1;
        std::vector<int> Ow_SS2;
        std::vector<int> Ow_outside;
        for (auto i = 0; i < NUM_H2O; i++)
        {
            double dL[3];
            for (auto k = 0; k < 3; k++)
            {
                dL[k] = Ow[i][k] - Zn[m][k];
            }
            PBC(dL);
            double distance = vector_magnitude(dL);

            // check if present in SS1
            if (distance < SS1_CUTOFF)
            {
                Ow_SS1.push_back(i);
            }

            // check if present in SS2
            else if (distance > SS1_CUTOFF && distance < SS2_CUTOFF)
            {
                Ow_SS2.push_back(i);
            }

            else
            {
                Ow_outside.push_back(i);
            }
        }

        // std::vector<int> TFSI_SS1;
        // std::vector<int> TFSI_SS2;
        // std::vector<int> TFSI_outside;

        // check distance from every Ot
        std::vector<int> Ot_SS1;
        std::vector<int> Ot_SS2;
        std::vector<int> Ot_outside;
        for (auto i = 0; i < NUM_TFSI * 4; i++)
        {
            // int tfsi_number = i % 4;
            double dL[3];
            for (auto k = 0; k < 3; k++)
            {
                dL[k] = Ot[i][k] - Zn[m][k];
            }
            PBC(dL);
            double distance = vector_magnitude(dL);

            // check if present in SS1
            if (distance < SS1_CUTOFF)
            {
                Ot_SS1.push_back(i);
            }

            // check if present in SS2
            else if (distance > SS1_CUTOFF && distance < SS2_CUTOFF)
            {
                Ot_SS2.push_back(i);
            }
            else
            {
                Ot_outside.push_back(i);
            }
        }

        // check distance from every F
        std::vector<int> F_SS1;
        std::vector<int> F_SS2;
        std::vector<int> F_outside;
        for (auto i = 0; i < NUM_TFSI * 6; i++)
        {
            // int tfsi_number = i % 6;
            double dL[3];
            for (auto k = 0; k < 3; k++)
            {
                dL[k] = F[i][k] - Zn[m][k];
            }
            PBC(dL);
            double distance = vector_magnitude(dL);

            // check if present in SS1
            if (distance < SS1_CUTOFF)
            {
                F_SS1.push_back(i);
            }

            // check if present in SS2
            else if (distance > SS1_CUTOFF && distance < SS2_CUTOFF)
            {
                F_SS2.push_back(i);
            }
            else
            {
                F_outside.push_back(i);
            }
        }

        // check distance from every N
        std::vector<int> N_SS1;
        std::vector<int> N_SS2;
        std::vector<int> N_outside;
        for (auto i = 0; i < NUM_TFSI; i++)
        {
            // int tfsi_number = i % 6;
            double dL[3];
            for (auto k = 0; k < 3; k++)
            {
                dL[k] = N[i][k] - Zn[m][k];
            }
            PBC(dL);
            double distance = vector_magnitude(dL);

            // check if present in SS1
            if (distance < SS1_CUTOFF)
            {
                N_SS1.push_back(i);
            }

            // check if present in SS2
            else if (distance > SS1_CUTOFF && distance < SS2_CUTOFF)
            {
                N_SS2.push_back(i);
            }
            else
            {
                N_outside.push_back(i);
            }
        }

        // SS1 composition
        std::ofstream comp1Zn{"Cation/ZnSS1comp.txt", std::ios::app};
        if (!comp1Zn)
        {
            std::cerr << "Uh oh, SS1compZn.txt could not be opened for writing!\n";
            return 1;
        }
        comp1Zn << Ow_SS1.size() << " " << Ot_SS1.size() << " " << F_SS1.size() << " " << N_SS1.size() << "\n";
        comp1Zn.close();

        // SS2 composition
        std::ofstream comp2Zn{"Cation/ZnSS2comp.txt", std::ios::app};
        if (!comp2Zn)
        {
            std::cerr << "Uh oh, SS2comp.txt could not be opened for writing!\n";
            return 1;
        }
        comp2Zn << Ow_SS2.size() << " " << Ot_SS2.size() << " " << F_SS2.size() << " " << N_SS2.size() << "\n";
        comp2Zn.close();

        // bonding of SS1 water
        int SS1_data[Ow_SS1.size()][9];
        for (auto i = 0; i < Ow_SS1.size(); ++i)
        {
            for (auto j = 0; j < 9; ++j)
            {
                SS1_data[i][j] = 0;
            }
        }
        // fill 1st column of SS1 data with water tag numbers
        for (auto i = 0; i < Ow_SS1.size(); i++)
        {
            SS1_data[i][0] = Ow_SS1[i];
        }

        int SS2_data[Ow_SS2.size()][13];
        for (auto i = 0; i < Ow_SS2.size(); ++i)
        {
            for (auto j = 0; j < 13; ++j)
            {
                SS2_data[i][j] = 0;
            }
        }
        // fill 1st column of SS2 data with water tag numbers
        for (auto i = 0; i < Ow_SS2.size(); i++)
        {
            SS2_data[i][0] = Ow_SS2[i];
        }

        // 1) SS1 H1 --> SS2 Ow
        int i_counter = 0;
        for (const auto &i : Ow_SS1)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H1[i][k] - Ow[i][k];
            }

            int j_counter = 0;
            for (const auto &j : Ow_SS2)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = Ow[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS1_data[i_counter][1]++;
                    SS2_data[j_counter][2]++;
                }
                j_counter++;
            }
            i_counter++;
        }

        // 2) SS1 H2 --> SS2 Ow
        i_counter = 0;
        for (const auto &i : Ow_SS1)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H2[i][k] - Ow[i][k];
            }
            PBC(v1);

            int j_counter = 0;
            for (const auto &j : Ow_SS2)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = Ow[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS1_data[i_counter][1]++;
                    SS2_data[j_counter][2]++;
                }
                j_counter++;
            }
            i_counter++;
        }

        // 3) SS2 H1 --> SS1 Ow
        i_counter = 0;
        for (const auto &i : Ow_SS2)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H1[i][k] - Ow[i][k];
            }
            PBC(v1);

            int j_counter = 0;
            for (const auto &j : Ow_SS1)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = Ow[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS1_data[j_counter][1]++;
                    SS2_data[i_counter][2]++;
                }
                j_counter++;
            }
            i_counter++;
        }

        // 4) SS2 H2 --> SS1 Ow
        i_counter = 0;
        for (const auto &i : Ow_SS2)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H2[i][k] - Ow[i][k];
            }
            PBC(v1);

            int j_counter = 0;
            for (const auto &j : Ow_SS1)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = Ow[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS1_data[j_counter][1]++;
                    SS2_data[i_counter][2]++;
                }
                j_counter++;
            }
            i_counter++;
        }

        // 5) SS1 H --> SS2 TFSI (Ot,F,N)
        // H1 --> Ot
        i_counter = 0;
        for (const auto &i : Ow_SS1)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H1[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : Ot_SS2)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = Ot[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS1_data[i_counter][3]++;
                }
            }
            i_counter++;
        }

        // H2 --> Ot
        i_counter = 0;
        for (const auto &i : Ow_SS1)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H2[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : Ot_SS2)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = Ot[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS1_data[i_counter][3]++;
                }
            }
            i_counter++;
        }

        // H1 --> F
        i_counter = 0;
        for (const auto &i : Ow_SS1)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H1[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : F_SS2)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = F[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS1_data[i_counter][4]++;
                }
            }
            i_counter++;
        }

        // H2 --> F
        i_counter = 0;
        for (const auto &i : Ow_SS1)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H2[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : F_SS2)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = F[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS1_data[i_counter][4]++;
                }
            }
            i_counter++;
        }

        // H1 --> N
        i_counter = 0;
        for (const auto &i : Ow_SS1)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H1[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : N_SS2)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = N[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS1_data[i_counter][5]++;
                }
            }
            i_counter++;
        }

        // H2 --> N
        i_counter = 0;
        for (const auto &i : Ow_SS1)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H2[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : N_SS2)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = N[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS1_data[i_counter][5]++;
                }
            }
            i_counter++;
        }

        // 6) SS1 H --> SS1 TFSI (Ot,F,N)
        // H1 --> Ot
        i_counter = 0;
        for (const auto &i : Ow_SS1)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H1[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : Ot_SS1)
            {
                if (i == j)
                    continue;
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = Ot[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS1_data[i_counter][6]++;
                }
            }
            i_counter++;
        }

        // H2 --> Ot
        i_counter = 0;
        for (const auto &i : Ow_SS1)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H2[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : Ot_SS1)
            {
                if (i == j)
                    continue;
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = Ot[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS1_data[i_counter][6]++;
                }
            }
            i_counter++;
        }

        // H1 --> F
        i_counter = 0;
        for (const auto &i : Ow_SS1)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H1[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : F_SS1)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = F[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS1_data[i_counter][7]++;
                }
            }
            i_counter++;
        }

        // H2 --> F
        i_counter = 0;
        for (const auto &i : Ow_SS1)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H2[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : F_SS1)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = F[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS1_data[i_counter][7]++;
                }
            }
            i_counter++;
        }

        // H1 --> N
        i_counter = 0;
        for (const auto &i : Ow_SS1)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H1[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : N_SS1)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = N[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS1_data[i_counter][8]++;
                }
            }
            i_counter++;
        }

        // H2 --> N
        i_counter = 0;
        for (const auto &i : Ow_SS1)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H2[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : N_SS1)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = N[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS1_data[i_counter][8]++;
                }
            }
            i_counter++;
        }

        // 7) SS1 H1 --> SS1 Ow
        i_counter = 0;
        for (const auto &i : Ow_SS1)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H1[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : Ow_SS1)
            {
                if (i == j)
                    continue;
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = Ow[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS1_data[i_counter][2]++;
                }
            }
            i_counter++;
        }

        // 8) SS1 H2 --> SS1 Ow
        i_counter = 0;
        for (const auto &i : Ow_SS1)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H2[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : Ow_SS1)
            {
                if (i == j)
                    continue;
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = Ow[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS1_data[i_counter][2]++;
                }
            }
            i_counter++;
        }

        // bonding of SS2 water

        // 1) SS2 Ow --> SS2 Ow
        // H1 --> Ow
        i_counter = 0;
        for (const auto &i : Ow_SS2)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H1[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : Ow_SS2)
            {
                if (i == j)
                    continue;

                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = Ow[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS2_data[i_counter][1]++;
                }
            }
            i_counter++;
        }

        // H2 --> Ow
        i_counter = 0;
        for (const auto &i : Ow_SS2)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H2[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : Ow_SS2)
            {
                if (i == j)
                    continue;

                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = Ow[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS2_data[i_counter][1]++;
                }
            }
            i_counter++;
        }

        // 2) SS2 Ow --> SS2 TFSI
        // H1 --> Ot
        i_counter = 0;
        for (const auto &i : Ow_SS2)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H1[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : Ot_SS2)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = Ot[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS2_data[i_counter][4]++;
                }
            }
            i_counter++;
        }

        // H2 --> Ot
        i_counter = 0;
        for (const auto &i : Ow_SS2)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H2[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : Ot_SS2)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = Ot[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS2_data[i_counter][4]++;
                }
            }
            i_counter++;
        }

        // H1 --> F
        i_counter = 0;
        for (const auto &i : Ow_SS2)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H1[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : F_SS2)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = F[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS2_data[i_counter][5]++;
                }
            }
            i_counter++;
        }

        // H2 --> F
        i_counter = 0;
        for (const auto &i : Ow_SS2)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H2[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : F_SS2)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = F[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS2_data[i_counter][5]++;
                }
            }
            i_counter++;
        }

        // H1 --> N
        i_counter = 0;
        for (const auto &i : Ow_SS2)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H1[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : N_SS2)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = N[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS2_data[i_counter][6]++;
                }
            }
            i_counter++;
        }

        // H2 --> N
        i_counter = 0;
        for (const auto &i : Ow_SS2)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H2[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : N_SS2)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = N[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS2_data[i_counter][6]++;
                }
            }
            i_counter++;
        }

        // 3) SS2 Ow --> Outside Ow
        // SS2 H1 --> Outside Ow
        i_counter = 0;
        for (const auto &i : Ow_SS2)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H1[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : Ow_outside)
            {

                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = Ow[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS2_data[i_counter][3]++;
                }
            }
            i_counter++;
        }

        // SS2 H2 --> Outside Ow
        i_counter = 0;
        for (const auto &i : Ow_SS2)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H2[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : Ow_outside)
            {

                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = Ow[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS2_data[i_counter][3]++;
                }
            }
            i_counter++;
        }

        // outside H1 --> SS2 Ow
        i_counter = 0;
        for (const auto &i : Ow_outside)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H1[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : Ow_SS2)
            {

                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = Ow[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS2_data[i_counter][3]++;
                }
            }
            i_counter++;
        }

        // Outside H2 --> SS2 Ow
        i_counter = 0;
        for (const auto &i : Ow_outside)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H2[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : Ow_SS2)
            {

                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = Ow[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS2_data[i_counter][3]++;
                }
            }
            i_counter++;
        }

        // 4) SS2 Ow --> Outside TFSI (Ot,N,F)
        // H1 --> Ot
        i_counter = 0;
        for (const auto &i : Ow_SS2)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H1[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : Ot_outside)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = Ot[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS2_data[i_counter][10]++;
                }
            }
            i_counter++;
        }

        // H2 --> Ot
        i_counter = 0;
        for (const auto &i : Ow_SS2)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H2[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : Ot_outside)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = Ot[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS2_data[i_counter][10]++;
                }
            }
            i_counter++;
        }

        // H1 --> F
        i_counter = 0;
        for (const auto &i : Ow_SS2)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H1[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : F_outside)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = F[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS2_data[i_counter][11]++;
                }
            }
            i_counter++;
        }

        // H2 --> F
        i_counter = 0;
        for (const auto &i : Ow_SS2)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H2[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : F_outside)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = F[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS2_data[i_counter][11]++;
                }
            }
            i_counter++;
        }

        // H1 --> N
        i_counter = 0;
        for (const auto &i : Ow_SS2)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H1[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : N_outside)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = N[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS2_data[i_counter][12]++;
                }
            }
            i_counter++;
        }

        // H2 --> N
        i_counter = 0;
        for (const auto &i : Ow_SS2)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H2[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : N_outside)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = N[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS2_data[i_counter][12]++;
                }
            }
            i_counter++;
        }

        // 5) SS2 Ow --> SS1 TFSI
        // H1 --> Ot
        i_counter = 0;
        for (const auto &i : Ow_SS2)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H1[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : Ot_SS1)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = Ot[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS2_data[i_counter][7]++;
                }
            }
            i_counter++;
        }

        // H2 --> Ot
        i_counter = 0;
        for (const auto &i : Ow_SS2)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H2[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : Ot_SS1)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = Ot[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS2_data[i_counter][7]++;
                }
            }
            i_counter++;
        }

        // H1 --> F
        i_counter = 0;
        for (const auto &i : Ow_SS2)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H1[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : F_SS1)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = F[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS2_data[i_counter][8]++;
                }
            }
            i_counter++;
        }

        // H2 --> F
        i_counter = 0;
        for (const auto &i : Ow_SS2)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H2[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : F_SS1)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = F[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS2_data[i_counter][8]++;
                }
            }
            i_counter++;
        }

        // H1 --> N
        i_counter = 0;
        for (const auto &i : Ow_SS2)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H1[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : N_SS1)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = N[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS2_data[i_counter][9]++;
                }
            }
            i_counter++;
        }

        // H2 --> N
        i_counter = 0;
        for (const auto &i : Ow_SS2)
        {
            double v1[3];
            for (auto k = 0; k < 3; k++)
            {
                v1[k] = H2[i][k] - Ow[i][k];
            }
            PBC(v1);

            for (const auto &j : N_SS1)
            {
                double v2[3], dL[3];
                for (auto k = 0; k < 3; k++)
                {
                    v2[k] = N[j][k] - Ow[i][k];
                }
                PBC(v2);

                for (auto k = 0; k < 3; k++)
                {
                    dL[k] = v2[k] - v1[k];
                }

                double angle = angle_between_vectors(v1, v2);

                double distance = vector_magnitude(dL);

                if (distance < H_BOND_DISTANCE_CUTOFF && angle < H_BOND_ANGLE_CUTOFF)
                {
                    SS2_data[i_counter][9]++;
                }
            }
            i_counter++;
        }

        // SS1 bonding output
        std::ofstream bond1Zn{"Cation/ZnSS1bond.txt", std::ios::app};
        if (!bond1Zn)
        {
            std::cerr << "Uh oh, ZnSS1bond.txt could not be opened for writing!\n";
            return 1;
        }
        for (auto i = 0; i < Ow_SS1.size(); i++)
        {
            for (auto j = 1; j < 9; j++)
            {
                bond1Zn << SS1_data[i][j] << " ";
            }
            bond1Zn << "\n";
        }
        bond1Zn.close();

        // SS2 bonding output
        std::ofstream bond2Zn{"Cation/ZnSS2bond.txt", std::ios::app};
        if (!bond2Zn)
        {
            std::cerr << "Uh oh, ZnSS2bond.txt could not be opened for writing!\n";
            return 1;
        }
        for (auto i = 0; i < Ow_SS2.size(); i++)
        {
            for (auto j = 1; j < 13; j++)
            {
                bond2Zn << SS2_data[i][j] << " ";
            }
            bond2Zn << "\n";
        }
        bond2Zn.close();
    }

    return 0;
}

//  define functions

void readSystemInfo()
{
    std::ifstream sysInfo("SysInfo.txt");
    if (sysInfo.is_open())
    {
        sysInfo >> BOX_LENGTH >> NUM_H2O >> NUM_TFSI >> NUM_LI >> NUM_ZN;
        sysInfo.close();
    }
    else
    {
        std::cerr << "Unable to open file." << std::endl;
    }
}

void readCoordinates(const std::string &fileName, double coordinates[][3], int numCoordinates)
{
    std::ifstream inputFile(fileName);
    if (inputFile.is_open())
    {
        for (int i = 0; i < numCoordinates; ++i)
        {
            inputFile >> coordinates[i][0] >> coordinates[i][1] >> coordinates[i][2];
        }
        inputFile.close();
    }
    else
    {
        std::cerr << "Unable to open file: " << fileName << std::endl;
    }
}

double calculateDistance(const double point1[3], const double point2[3])
{
    return std::sqrt(std::pow(point1[0] - point2[0], 2) +
                     std::pow(point1[1] - point2[1], 2) +
                     std::pow(point1[2] - point2[2], 2));
}

void PBC(double vector[3])
{
    int j = 0;
    for (j = 0; j < 3; j++)
    {
        if (vector[j] > BOX_LENGTH / 2)
        {
            vector[j] -= BOX_LENGTH;
        }
        else if (vector[j] < -BOX_LENGTH / 2)
        {
            vector[j] += BOX_LENGTH;
        }
    }
}

double vector_magnitude(double m[3])
{
    return std::sqrt(m[0] * m[0] + m[1] * m[1] + m[2] * m[2]);
}

double dot_product(double m[3], double n[3])
{

    double dt = 0.0;
    for (auto i = 0; i < 3; i++)
    {
        dt += m[i] * n[i];
    }
    return dt;
}

double angle_between_vectors(double v1[3], double v2[3])
{
    double dot = dot_product(v1, v2);
    double mag1 = vector_magnitude(v1);
    double mag2 = vector_magnitude(v2);
    return 57.2958 * acos((dot) / (mag1 * mag2));
}
