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
#include <mpi.h>

using namespace std;

int main(int argc, char** argv)
{	

    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	
	int N = 200;
	int S = world_size-1;
	int R = N/S;
	int C = N;

    if(world_rank == 0){
		
	// -------------- PREPARE AN ARRAY TO BE FILLED --------------

		int** X = new int*[N];
		for(int i = 0; i < N; i++) { X[i] = new int[N]; }	


	// -------------- READ INPUT FILE AND FILL THE PREPARED ARRAY --------------

	    ifstream input_file(argv[1]);

	    if(input_file.is_open()){

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

	    } else {
			throw "Unable to Open Input File!";
		}

	// -------------- DISTRIBUTE DATA - TO SLAVES --------------

		for(int s = 0; s < S; s++)
			for(int r = 0; r < R; r++) 
 				MPI_Send(X[s*R + r], N, MPI_INT, s+1, r, MPI_COMM_WORLD);
				
    // -------------- -------------- --------------

    } else {

	// -------------- COLLECT DATA - FROM MASTER --------------
		
		int X[R][N];

		for(int r = 0; r < R; r++)
			MPI_Recv(X[r], N, MPI_INT, 0, r, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	// -------------- DENOISE - INITIALIZE VALUES --------------

		int Z[R][N];

		for(int r = 0; r < R; r++)
			for(int c = 0; c < N; c++)
				Z[r][c] = X[r][c];

		int EXTRA_ROW_ABOVE[N];
		int EXTRA_ROW_BELOW[N];

		double beta = atof(argv[3]);

		double pi = atof(argv[4]);
		double gamma = 0.5*log((1-pi)/pi);

		long iteration_limit = 100000 / S;

		srand(time(NULL));

		for(long iteration_number = 1; iteration_number < iteration_limit; iteration_number++) {

		// -------------- DISTRIBUTE LAST ROWS -------------- 

			if (world_rank != 1)
				MPI_Recv(EXTRA_ROW_ABOVE, N, MPI_INT, world_rank - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			MPI_Send(Z[R-1], N, MPI_INT, world_rank % S + 1, 1, MPI_COMM_WORLD);

			if (world_rank == 1) 
				MPI_Recv(EXTRA_ROW_ABOVE, N, MPI_INT, S, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		// -------------- DISTRIBUTE FIRST ROWS -------------- 

			if (world_rank != S) 
				MPI_Recv(EXTRA_ROW_BELOW, N, MPI_INT, world_rank + 1, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				
			MPI_Send(Z[0], N, MPI_INT, (world_rank - 2 + S) % S + 1, 2, MPI_COMM_WORLD);

			if (world_rank == S) 
				MPI_Recv(EXTRA_ROW_BELOW, N, MPI_INT, 1, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		
		// -------------- DENOISE - PROPOSE BIT FLIP -------------- 

		    int i = rand() % R;
		    int j = rand() % N;

		    int neighbour_sum = 0;

		    for(int k = -1 ; k <= 1; k++) {

		        for(int l = -1; l <= 1; l++) {

		            int n_i = i + k;
		            int n_j = j + l;

		            if((n_i > -1 && n_i < R) && (n_j > -1 && n_j < N)) {

		                if(!(n_i == i && n_j == j)) {
		                    neighbour_sum += Z[n_i][n_j];
		                }
		            }
		        }

		    }

		    double delta_E = -2 * gamma * Z[i][j] * X[i][j] - 2 * beta * Z[i][j] * neighbour_sum;

		    if ( ((double) rand() / (RAND_MAX)) < exp(delta_E) ) { Z[i][j] = -Z[i][j]; }
		}

    // -------------- DISTRIBUTE DATA - TO MASTER --------------	
		
		for(int r = 0; r < R; r++)
			MPI_Send(Z[r], N, MPI_INT, 0, (world_rank-1)*R + r, MPI_COMM_WORLD);
		

    // -------------- -------------- --------------

	}

	MPI_Barrier(MPI_COMM_WORLD);
	
	if(world_rank == 0) {
	
	// -------------- COLLECT DATA - FROM SLAVES --------------
		int Z[N][N];
		
		for(int s = 0; s < S; s++)
			for(int r = 0; r < R; r++)
				MPI_Recv(Z[s*R + r], N, MPI_INT, s+1, s*R + r, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	
	// -------------- WRITE TO OUTPUT FILE --------------

		ofstream output_file(argv[2]);

		if (output_file.is_open()){

		    for(int row = 0; row < N; row++) {

		        for(int column = 0; column < N; column++) {
		            output_file << Z[row][column] << " ";
		        }

		        output_file << "\n";
		    }

		    output_file.close();

		} else {
		    throw "Unable to Open Output File!";
		}

    // -------------- -------------- --------------

	}

    // Finalize the MPI environment. No more MPI calls can be made after this

	MPI_Barrier(MPI_COMM_WORLD);

    MPI_Finalize();

    return 0;
}
