#include <iostream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <string>   
#include <sstream>
#include <fstream>
#include <iomanip>
#include <mpi.h>

#define MC_STEPS     1000000     // NUMBER OF MONTE CARLO CYCLES    ! the 
                                 // effective value includes a N^2 term.
#define THERM_STEPS  50000       // NUMBER OF STEPS FOR THERMALIZATION this is
                                 // around 5% of the total number of steps

using namespace std;
FILE * output;


// RANDOM NUMBER FUNCTION
inline double ran(){
	return ((double) rand()) / RAND_MAX;
}


// FUNCTION TO SET UP THE INITIAL VALUES FOR THE MATRIX AND E
void initialize_matrix(short ** lat, int N, double &E){
	E = 0;
	for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            //(ran() < 0.5) ? lat[i][j] = 1 : lat[i][j] = -1;
            lat[i][j] = -1;
	
	for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
			E += -lat[i][j] * (lat[i][(j-1 == -1 ? N-1 : j-1)]  + 
                               lat[(i-1 == -1 ? N-1 : i-1)][j]);
}


// MAIN FUNCTION OF THE PROGRAM: given the matrix tries to flip a spin and 
// saves the current energy state and total flip counter.
void CountStates(short ** lat, int N, double * deltaE, int * states,
				 int &pivot, int &flips){
	int spin_combination, i, j;
	
	// select random item
	i = (int) (ran()*(N-1));
	j = (int) (ran()*(N-1));
	spin_combination = 2 * lat[i][j] * (lat[i] [(j+1 == N  ? 0   : j+1)] +
										lat[i] [(j-1 == -1 ? N-1 : j-1)] +
										lat[(i+1 == N  ? 0   : i+1)] [j] +
										lat[(i-1 == -1 ? N-1 : i-1)] [j]);
	
	// select if we can accept the move and update parameters
	if(deltaE [spin_combination + 8] >= ran()) {
		lat[i][j] *= -1;
		
		// move the energy flag over the states array and add 1 to flip counter
		pivot += spin_combination;
		flips += 1;
	}
	
	// add 1 to the current state counter
	states[pivot] += 1;
}



int main(int argc, char* argv[]){	
	
	// Matrix and general variables
	short ** lat;
	int * states;	
	int * total_states;
	int N, rank, my_rank, numprocs, pivot, flips, total_flips;
	double E;
	
	// initialization of MPI
	MPI_Init (&argc, &argv);
	MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
	
	// set up the limits of calculation that are process dependent
	int slice_size  = MC_STEPS / numprocs;
	int slice_begin = my_rank * slice_size + 1;
	int slice_end   = (my_rank + 1) * slice_size;
	if (my_rank == numprocs-1)
		slice_end = MC_STEPS;
	
	// initialize seed with a different number for each processor
	srand(time(NULL)-my_rank);
	rank = my_rank;
	
	// create empty output file to append the data later
	if(rank == 0){
		stringstream file1;
		string file_name1;
		const char * file_namec1;
		file1 << "Data/flips.txt";
		file_name1 = file1.str();
		file_namec1 =  file_name1.c_str();
		output = fopen(file_namec1, "w+");
		fclose(output);	
	}
	
	
	N = 20;    // SIZE OF THE MATRIX
	
	// loop over temperatures
	for(double T = 1.0; T < 3.1; T += 0.1){
		
		// create the matrix
		lat = new short * [N];
		for (int i = 0; i < N; i++)
			lat[i] = new short [N];
		
		// create the counter variable for states and flips
		states = new int [4*N*N+1];
		total_states = new int [4*N*N+1];
		flips = 0;
		total_flips = 0;
		
		for (int i = 0; i < 4*N*N+1; i++){
			states[i] = 0;
			total_states[i] = 0;
		}
		
		// compute the energy difference and the corresponding exponential
		// in advance for the current temperature
		double deltaE[17];
		for(int i =-8; i <= 8; i++)
			deltaE[i+8] = 0;
		for(int i =-8; i <= 8; i+=4)
			deltaE[i+8] = exp(-i/(T));
		
		// set up matrix and calculate initial parameters 
		initialize_matrix(lat, N, E);
		pivot = E + 2*N*N;
		
		// small loop for thermalization. no data will be stored from this
		for(int i = slice_begin; i <= slice_end; i++){
			CountStates(lat, N, deltaE, states, pivot, flips);
		}
		
		// erase thermalization data, we don't need it
		for (int i = 0; i < 4*N*N+1; i++){
			states[i] = 0;
			total_states[i] = 0;
		}
		
		// main loop over all the montecarlo cycles. at each loop all the 
		// counters are updated.
        for(int i = slice_begin; i <= slice_end; i++)
			CountStates(lat, N, deltaE, states, pivot, flips);
		
		// merge data from all the processes into one variable
		for (int i = 0; i < 4*N*N+1; i+=4)
			MPI_Reduce(&states[i], &total_states[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&flips, &total_flips, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		
		// compute the final values and write output
		if ( rank == 0){
			stringstream file;
			string file_name;
			const char * file_namec;
			
			// file for energy counters at current temperature
			file << "Data/" << (int) (T*10) << ".txt";
			file_name = file.str();
			file_namec =  file_name.c_str();
			output = fopen(file_namec, "w+");
			
			// write out the counts for each energy level at temperature T
			for (int i = 0; i < 4*N*N+1; i+=4)
				fprintf(output, "%d,   %d\n", i-2*N*N, total_states[i]);
			fclose(output);	
			
			// compute the expectation value and the a posteriori variance 
			double expected = 0;
			double var = 0;
			for (int i = 0; i < 4*N*N+1; i+=4){
				expected += (double)(i-2*N*N)  / MC_STEPS * total_states[i];
				var += (double)(i-2*N*N) * (double)(i-2*N*N) 
						/ MC_STEPS * total_states[i];
			}
			
			// file for total accepted moves, energy expectation value and
			// variance at different temperatures
			stringstream file1;
			string file_name1;
			const char * file_namec1;
			file1 << "Data/flips.txt";
			file_name1 = file1.str();
			file_namec1 =  file_name1.c_str();
			output = fopen(file_namec1, "a");
			fprintf(output, "%.5f,    %d,    %.5f,    %.5f\n", 
					T, total_flips, expected/N/N, (var-expected*expected)/N/N);
			fclose(output);
		}
		
		
		// delete matrix and vectors
		for (int i = 0; i < N; i++)
			delete [] lat[i];
		delete [] lat;
		delete [] states;
		delete [] total_states;
	}
	
	// close MPI before exit
	MPI_Finalize();
	return 0;
}

