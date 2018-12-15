#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <cstring>
#include <cstdlib>
#include <math.h>
#include <stdlib.h>
#include <time.h>

using namespace std;

int main(int argc, char *argv[])
{
    int X[200][200];

    // -------------- READ INPUT FILE --------------

    ifstream input_file(argv[1]);

    if (input_file.is_open()){

        int row = 0;

        string line;

        while(getline(input_file, line)){

            int column = 0;

            string element;

            istringstream iss(line);

            while(iss >> element) { X[row][column++] = atoi(element.c_str()); }

            row++;
        }

        input_file.close();

    } else
        throw "Unable to Open Input File!";

    // -------------- -------------- --------------

    // -------------- DENOISE --------------

    int Z[200][200];

    for(int row = 0; row < 200; row++) {
        for(int column = 0; column < 200; column++) {
            Z[row][column] = X[row][column];
        }
    }

    double beta = atof(argv[3]);

    double pi = atof(argv[4]);
    double gamma = 0.5*log((1-pi)/pi);

    long iteration_limit = 1000000;

    srand(time(NULL));

    for(long iteration_number = 1; iteration_number < iteration_limit; iteration_number++) {

        int i = rand() % 200;
        int j = rand() % 200;

        int neighbour_sum = 0;

        for(int k = -1 ; k <= 1; k++) {

            for(int l = -1; l <= 1; l++) {

                int n_i = i + k;
                int n_j = j + l;

                if((n_i > -1 && n_i < 200) && (n_j > -1 && n_j < 200)) {

                    if(!(n_i == i && n_j == j)) {
                        neighbour_sum += Z[n_i][n_j];
                    }
                }
            }

        }

        double delta_E = -2 * gamma * Z[i][j] * X[i][j] - 2 * beta * Z[i][j] * neighbour_sum;

        if ( ((double) rand() / (RAND_MAX)) < exp(delta_E) ) { Z[i][j] = -Z[i][j]; }
    }

    // -------------- WRITE TO OUTPUT FILE --------------

    ofstream output_file(argv[2]);

    if (output_file.is_open()){

        for(int row = 0; row < 200; row++) {

            for(int column = 0; column < 200; column++) {
                output_file << Z[row][column] << " ";
            }

            output_file << "\n";
        }

        output_file.close();

    } else
        throw "Unable to Open Output File!";

    // -------------- -------------- --------------

    return 0;
}
