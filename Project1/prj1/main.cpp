#include <iostream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <string>   
#include <sstream>

#include "lib.h"

using namespace std;

//The function describing the radial distribution of charge
//     f(x) = 100 e^(-10x)
inline double f(double x) {
	return exp(-10*x) * 100;
}


//The analytical solution of the differential equation
//	f(x) = 1 - (1 - e^(-10))x - e^(-10x)
inline double ua (double x) {
	return 1 - (1 - exp(-10)) * x - exp(-10*x);
}


//Fuction to compute the solution of the differential equation using the gaussian elimination
//on the tridiagonal matrix
//      - u''(x) + f(x) = 0
double * compute_vectors (unsigned int size, double *A, double *B, double *C, double *D, bool speed = 0) {
	double *v = new double[size];
	
	//General tridiagonal matrix solver
	if(!speed){
		//forward substitution
		for (unsigned int i = 1; i < size; i++){
			B[i] = B[i] - A[i] / B[i - 1] * C[i];
			D[i] = D[i] - D[i - 1] * A[i] / B[i - 1];
		}
		//backword substitution
		v[size - 1] = D[size - 1] / B[size - 1];
		for (int i = size - 2; i >= 0; i--)
			v[i] = (D[i] -  C[i] * v[i + 1] ) / B[i];
	}
	
	
	//Optimized 6n algorithm for the particular problem
	else {
		//forward substitution
		for (unsigned int i = 1; i < size; i++){
			B[i] -= 1 / B[i - 1];
			D[i] += D[i - 1] / B[i - 1];
		}	
		//backward substitution
		v[size - 1] = D[size - 1] / B[size - 1];
		for (int i = size - 2; i >= 0; i--)
			v[i] = (D[i] + v[i + 1]) / B[i];
	}	
	return v;
	delete [] v;
}


double * compute_matrix (unsigned int size, double **A, double *B) {
	int *indx = new int[size];
	double d;
    
	ludcmp(A, size, indx, &d);
	lubksb(A, size, indx, B);
    
	delete [] indx;
	return B;
}

//Function to compute the error given a numerical solution and an analytical solution
//	err = | (v - u) / u |
double * compute_error (unsigned int size, double *V, double *A, bool absolute = 0) {
	double *E = new double[size];
	
	if(!absolute)
		for (unsigned int i = 0; i < size; i++)
			E[i] = abs((V[i] - A[i]) / A[i]);
	else
		for (unsigned int i = 0; i < size; i++)
			E[i] = abs(V[i] - A[i]);
	return E;
	delete [] E;
}


//Function to generate one output file containing the x position and the corresponding vector value
//Has flags to add tags if the data series is of an error function or a LU decomposition result
void generate_file (double *v, unsigned int size, double step, bool error = 0, bool LU = 0){
	FILE * output;
	stringstream file;
    
    //If it's a gaussian elimination data series
	if(!LU){
		if (!error)
			file << "Data/" << size << ".txt";
		
		//If it's an error  data series
		else
			file << "Data/" << size << "err.txt";
		string file_name = file.str();
		const char* file_namec =  file_name.c_str();
		output = fopen(file_namec, "w+");
	}
    
    //If it's a LU decomposition data series
	else{
		if (!error)
			file << "Data/" << size << "LU.txt";
		
		//If it's an error  data series
		else
			file << "Data/" << size << "errLU.txt";
		string file_name = file.str();
		const char* file_namec =  file_name.c_str();
		output = fopen(file_namec, "w+");
	}
    
    //Generate the file
	for(unsigned int i = 0; i < size; i++)
		fprintf(output, "%.20lf,  %.20lf\n", step * (i + 1), v[i]);
	fclose(output);
}



int main(int argn, char* argv[]) {
    
    // SECTION 1: VECTOR CALCULATION
	for (unsigned int n = 10; n <= 1000000; n *= 10) {	
		double x_0 = 0.0, x_1 = 1.0;
		double h = (x_1 - x_0) / (n + 1);
		
		double *a = new double[n];
		double *b = new double[n];
		double *c = new double[n];
		double *d = new double[n];
		double *u = new double[n];
		double *ue = new double[n];
		double *e = new double[n];
		
		for(unsigned int i = 0; i < n; i++){
			a[i] = -1;
			b[i] = 2;
			c[i] = -1;
			d[i] = f(x_0 + (i + 1) * h) * h * h;
			ue[i] = ua(x_0 + (i + 1) * h);
		}

		clock_t start, stop;
		start = clock();
		u = compute_vectors(n, a, b, c, d, 1);
		stop = clock();
		float time = ((float)stop-(float)start)/CLOCKS_PER_SEC;
		printf("Vector time %d \t = %.10lf\n", n, time);
		
		generate_file(u, n , h);
		e = compute_error(n, u, ue);
		generate_file(e, n , h);
		delete [] a;
		delete [] b;
		delete [] c;
		delete [] d;
		delete [] u;
		delete [] ue;
		delete [] e;
	}
    
    //SECTION 2: LU DECOMPOSITION
	for (unsigned int n = 10; n <= 1000; n *= 10) {
		double x_0 = 0.0, x_1 = 1.0;
		double h = (x_1 - x_0) / (n + 1);
		
		double ** A;
		double *d = new double[n];
		double *u = new double[n];
		double *ue = new double[n];
		double *e = new double[n];
		
		A = new double*[n];
		for (unsigned int i = 0; i < n; i++) {
			A[i] = new double[n];
			d[i] = f(x_0 + (i + 1) * h) * h * h;
			ue[i] = ua(x_0 + (i + 1) * h);
		}
		
		for (unsigned int i = 0; i < n; i++){
			for (unsigned int j = 0; j < n; j++){
				if(i == j && i != 0 && i != n - 1){
					A[i][j] = 2;
					A[i][j - 1] = -1;
					A[i][j + 1] = -1;
					j++;
				}
				else
					A[i][j] = 0;
			}
		}
		A[0][0] = 2;
		A[0][1] = -1;
		A[n - 1][n - 1] = 2;
		A[n - 1][n - 2] = -1;
		

		clock_t start, stop;
		start = clock();
		u = compute_matrix(n, A, d);
		stop = clock();
		float time = ((float)stop-(float)start)/CLOCKS_PER_SEC;
		printf("Matrix time %d  = %.10lf\n", n, time);
		
		generate_file(u, n, h, 0, 1);
		e = compute_error(n, u, ue);
		generate_file(e, n, h, 1, 1);
		
		for (unsigned int i = 0; i < n; i++)
			delete [] A[i];
		delete [] A;
		delete [] u;
		delete [] ue;
		delete [] e;
	}
	
	
	//SECTION 3: ERROR ANALYSIS
	FILE * output;
	output = fopen("Data/error.txt", "w+");
	
	for (unsigned int n = 10; n <= 1000000; n *= 1.1) {	
		double x_0 = 0.0, x_1 = 1.0;
		double h = (x_1 - x_0) / (n + 1);
		
		double *a = new double[n];
		double *b = new double[n];
		double *c = new double[n];
		double *d = new double[n];
		double *u = new double[n];
		double *ue = new double[n];
		double *e = new double[n];
		
		for(unsigned int i = 0; i < n; i++){
			a[i] = -1;
			b[i] = 2;
			c[i] = -1;
			d[i] = f(x_0 + (i + 1) * h) * h * h;
			ue[i] = ua(x_0 + (i + 1) * h);
		}
		u = compute_vectors(n, a, b, c, d, 1);
		e = compute_error(n, u, ue);
		fprintf(output, "%.20lf,  %.20lf\n", h, e[n/2]);
		
		delete [] a;
		delete [] b;
		delete [] c;
		delete [] d;
		delete [] u;
		delete [] ue;
		delete [] e;
	}
	fclose(output);
    
    return 0;
}

