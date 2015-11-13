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


// FUNCTION TO SET UP THE INITIAL VALUES FOR THE MATRIX, E AND M
void initialize_matrix(short ** lat, int N, double &E, double &M){
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			//(ran() < 0.5) ? lat[i][j] = 1 : lat[i][j] = -1;
			lat[i][j] = 1;
	
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++){
			M  += lat[i][j];
			E += -lat[i][j] * (lat[i][(j-1 == -1 ? N-1 : j-1)]  + 
							   lat[(i-1 == -1 ? N-1 : i-1)][j]);
	}         
}


// MAIN FUNCTION OF THE PROGRAM: given the matrix tries to flip N^2 and updates
// each time the values for magnetization and energy.
void MonteCarlo(short ** lat, int N, double * deltaE, double &E, double &M){
	int spin_combination, i, j;
	
	for(int n = 0; n < N*N; n++){
		// select random item
		i = (int) (ran()*(N-1));
		j = (int) (ran()*(N-1));
		
		// compute the energy difference of the eventual spin flip
		spin_combination = 2 * lat[i][j] * (lat[i] [(j+1 == N  ? 0   : j+1)] +
											lat[i] [(j-1 == -1 ? N-1 : j-1)] +
											lat[(i+1 == N  ? 0   : i+1)] [j] +
											lat[(i-1 == -1 ? N-1 : i-1)] [j]);
		
		// select if we can accept the move and update parameters
		if(deltaE [spin_combination + 8] >= ran()) {
			lat[i][j] *= -1;
			E  += (double) spin_combination;
			M  += (double) 2*lat[i][j];
		}
	}
}



int main(int argc, char* argv[]){	
	
	// Matrix and general variables
	short ** lat;
	int N_array[7] = {20, 40, 60, 80, 100, 150, 200}; // size of the matrices
	int N, n, rank, my_rank, numprocs;
	
	// variables for the thermodynamical values
	double E, M;
	double E0, dE0, M0, dM0, Ma0;
	double T0, Tf, dT, total_E, total_dE, total_M, total_dM, total_Ma;
	
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
	
	
	// loop over all sizes
	for(n = 0; n < 7; n++){	
		N = N_array[n];
		
		// create the matrixe of the right dimension
		lat = new short * [N];
		for (int i = 0; i < N; i++)
			lat[i] = new short [N];
		
		// variables for the ouput file and filename
		stringstream file;
		string file_name;
		const char * file_namec;
		
		
		// temperature intervals
		T0 = 1.6;   Tf = 3.01; dT = 0.05;
		
		// create empty output file with N as the name 
		if ( rank == 0){
			file << "Data/" << N << ".txt";
			file_name = file.str();
			file_namec =  file_name.c_str();
			output = fopen(file_namec, "w+");
			if(!output) exit(1);
			fclose(output);
		}
		
		// loop over temperatures
		for(double T = T0; T < Tf; T += dT){
			if ( rank == 0) printf("%.3f\n", T);
			
			// set the temperature step to be smaller if we are around Tc
			if(T > 1.98 && T < 2.49) dT = 0.01;
			else dT = 0.05;
			
			// compute the energy difference and the corresponding exponential
			// in advance for the current temperature
			double deltaE[17];
			for(int i =-8; i <= 8; i++)
				deltaE[i+8] = 0;
			for(int i =-8; i <= 8; i+=4)
				deltaE[i+8] = exp(-i/(T));
			
			// set all cumulative variables to 0
			E = M = 0;
			E0 = dE0 = M0 = dM0 = Ma0 = 0;
			total_E = total_dE = total_M = total_dM = total_Ma = 0;
			
			// set up matrix and calculate initial parameters 
			initialize_matrix(lat, N, E, M);
			
			// small loop for thermalization. no data will be stored from this
			for(int i = 0; i < THERM_STEPS; i++)
				MonteCarlo(lat, N, deltaE, E, M);
			
			// main loop over all the montecarlo cycles. at each loop all the 
			// parameters are updated.
			for(int i = slice_begin; i < slice_end; i++){
				MonteCarlo(lat, N, deltaE, E, M);
				E0  += E;
				dE0 += E*E;
				M0  += M;
				dM0 += M*M;
				Ma0 += fabs(M);
			}
			
			// merge data from all the processes into one variable
			MPI_Reduce(&E0,  &total_E,  1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&dE0, &total_dE, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&M0,  &total_M,  1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&dM0, &total_dM, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&Ma0, &total_Ma, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			
			// compute the final values and write output
			if ( rank == 0){
				output = fopen(file_namec, "a");
				
				// normalize to the number of cycles
				total_E  /= (double) MC_STEPS;
				total_M  /= (double) MC_STEPS;
				total_dE /= (double) MC_STEPS;
				total_dM /= (double) MC_STEPS;
				total_Ma /= (double) MC_STEPS;
				
				// compute heat capacity and susceptibility
				total_dE = (total_dE - total_E  * total_E ) /T/T;
				total_Ma = (total_dM - total_Ma * total_Ma) / T;
				
				// output
				fprintf(output, "%.3f    %.5f   %.5f   %.5f   %.5f\n",
						T, total_E/N/N, total_dE/N/N,
						total_M/N/N, total_Ma/N/N);
				fclose(output);
			}
		}
		
		// delete the matrix
		for (int i = 0; i < N; i++)
			delete [] lat[i];
		delete [] lat;
	}
	
	// close MPI before exit
	MPI_Finalize();
	return 0;
}
