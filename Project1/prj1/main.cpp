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


//Fuction to compute the solution of the differential equation
//      - u''(x) + f(x) = 0
double * compute_vectors (unsigned int size, double *A, double *B, double *C, double *D) {
    double *v = new double[size];
    double *t = new double[size];
    
    //forward substitution
    C[0] = C[0] / B[0];
    D[0] = D[0] / B[0];
    for (unsigned int i = 1; i < size; i++){
        double temp =  1 / (B[i] - A[i] * C[i - 1]);
	C[i] *= temp;
	D[i] = (D[i] - A[i] * D[i - 1]) * temp;
    }
    
    v[size - 1] = D[size - 1];
    //backword substitution
    for (unsigned int i = size - 2; i > 0; i--){
	v[i] = D[i] -  C[i] * v[i + 1];
    }
    return v;
    delete [] v;
    delete [] t;
}

//Function to compute the error given a numerical solution and an analytical solution
//	err = | (v - u) / u |
double * compute_error (unsigned int size, double *V, double *A) {
    double *E = new double[size];
    
    for (unsigned int i = 0; i < size; i++)
	E[i] = abs((V[i] - A[i]) / A[i]);
    return E;
    delete [] E;
}


//Function to generate one output file containing the x position and the corresponding vector value
void generate_file (double *v, unsigned int size, double step, bool error = 0, bool LU = 0){
    FILE * output;
    stringstream file;
    if(!LU){
	if (!error)
	    file << "Data/" << size - 1<< ".txt";
	else
	    file << "Data/" << size - 1<< "err.txt";
	string file_name = file.str();
	const char* file_namec =  file_name.c_str();
	output = fopen(file_namec, "w+");
    }
    else{
	if (!error)
	    file << "Data/" << size - 1<< "LU.txt";
	else
	    file << "Data/" << size - 1<< "errLU.txt";
	string file_name = file.str();
	const char* file_namec =  file_name.c_str();
	output = fopen(file_namec, "w+");
    }
    for(unsigned int i = 0; i < size; i++)
	fprintf(output, "%.20lf,  %.20lf\n", step * i, v[i]);
    fclose(output);
}

double * compute_matrix (unsigned int size, double **A, double *B) {
    int *indx = new int[size];
    double d;
    
    ludcmp(A, size, indx, &d);
    lubksb(A, size, indx, B);
    
    delete [] indx;
    return B;
}

int main(int argn, char* argv[]) {
    double x_0 = 0.0, x_1 = 1.0;
    
    // SECTION 1: VECTOR CALCULATION
    for (unsigned int n = 10; n <= 10000000; n *= 10) {
        double *a = new double[n];
        double *b = new double[n];
        double *c = new double[n];
	double *d = new double[n];
        double *u = new double[n + 1];
	double *ue = new double[n + 1];
	double *e = new double[n + 1];
	
	double h = (x_1 - x_0) / (n);
	
        for(unsigned int i = 0; i <= n; i++){
            a[i] = -1;
            b[i] = 2;
            c[i] = -1;
	    d[i] = f(x_0 + i * h) * h * h;
	    ue[i] = ua(x_0 + i * h);
        }
        
	clock_t start, stop;
	start = clock();
        u = compute_vectors(n, a, b, c, d);
	stop = clock();
	float time = ((float)stop-(float)start)/CLOCKS_PER_SEC;
	printf("Vector time %d \t = %.10lf\n", n, time);
	
	//u[n] = 0.0;
        generate_file(u, n + 1, h);
	e = compute_error(n, u, ue);
	generate_file(e, n + 1, h, 1);
        delete [] a;
        delete [] b;
        delete [] c;
	delete [] d;
        delete [] u;
	delete [] ue;
	delete [] e;
    }
    
    //SECTION 2: LU DECOMPOSITION
    for (unsigned int n = 10; n <= 1000; n *= 10) 
    {
	double h = (x_1 - x_0) / (n);
	
        double ** A;
	double *d = new double[n];
        double *u = new double[n + 1];
	double *ue = new double[n + 1];
	double *e = new double[n + 1];
	
	A = new double*[n];
	for (unsigned int i = 0; i < n; i++) {
	    A[i] = new double[n];
	    d[i] = f(x_0 + i * h) * h * h;
	    ue[i] = ua(x_0 + i * h);
	}
	
	for (unsigned int i = 0; i < n; i++){
	    for (unsigned int j = 0; j < n; j++){
		if(i == j && i != 0 && i != n - 1){
		    A[i][j] = -2;
		    A[i][j - 1] = 1;
		    A[i][j + 1] = 1;
		    j++;
		}
		else
		    A[i][j] = 0;
	    }
	}
	A[0][0] = -2;
	A[0][1] = 1;
	A[n - 1][n - 1] = -2;
	A[n - 1][n - 2] = 1;
	

	clock_t start, stop;
	start = clock();
        u = compute_matrix(n, A, d);
	stop = clock();
	float time = ((float)stop-(float)start)/CLOCKS_PER_SEC;
	printf("Matrix time %d  = %.10lf\n", n, time);
	
	for (unsigned int i = 0; i < n; i++)
	    u[i] = -u[i];
	//u[n] = 0.0;
	generate_file(u, n + 1, h, 0, 1);
	e = compute_error(n + 1, u, ue);
	generate_file(e, n + 1, h, 1, 1);
	
	for (unsigned int i = 0; i < n; i++)
	    delete [] A[i];
	delete [] A;
        delete [] u;
	delete [] ue;
	delete [] e;
    }
    return 0;
}

