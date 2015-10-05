#include <iostream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <string>   
#include <sstream>

#include "lib.h"

using namespace std;

inline double V_harm(double x, double useless) {
	return x*x;
}


inline double V_2e(double x, double omega) {
	return x*x* omega + 1.0/x;
}


void tqli1(double *d, double *e, int n, double **z)
{
   register int   m,l,iter,i,k;
   double         s,r,p,g,f,dd,c,b;
   for(i = 1; i < n; i++) e[i-1] = e[i];
   e[n] = 0.0;
   for(l = 0; l < n; l++) {
      iter = 0;
      do {
         for(m = l; m < n-1; m++) {
            dd = fabs(d[m]) + fabs(d[m+1]);
            if((double)(fabs(e[m])+dd) == dd) break;
         }
         
         if(m != l) {
            if(iter++ == 30) {
               printf("\n\nToo many iterations in tqli.\n");
               exit(1);
            }
            g = (d[l+1] - d[l])/(2.0 * e[l]);
            r = pythag(g,1.0);
            g = d[m]-d[l]+e[l]/(g+SIGN(r,g));
            s = c = 1.0;
            p = 0.0;
            for(i = m-1; i >= l; i--) {
               f      = s * e[i];
               b      = c*e[i];
               e[i+1] = (r=pythag(f,g));
               if(r == 0.0) {
                  d[i+1] -= p;
                  e[m]    = 0.0;
                  break;
               }
               s      = f/r;
               c      = g/r;
               g      = d[i+1] - p;
               r      = (d[i] - g) * s + 2.0 * c * b;
               d[i+1] = g + (p = s * r);
               g      = c * r - b;
               for(k = 0; k < n; k++) {
                  f         = z[k][i+1];
                  z[k][i+1] = s * z[k][i] + c * f;
                  z[k][i]   = c * z[k][i] - s * f;
               } /* end k-loop */
            } /* end i-loop */
            if(r == 0.0 && i >= l) continue;
            d[l] -= p;
            e[l]  = g;
            e[m]  = 0.0;
         } /* end if-loop for m != 1 */
      } while(m != l);
   } /* end l-loop */
} /* End: function tqli(), (C) Copr. 1986-92 Numerical Recipes Software )%. */



//FUNCTION TO CREATE THE INITIAL MATRIX FOR THE DIMENSIONLES SCHR. EQ.
void generate_matrix(double ** A, unsigned int size, double step, double omega, double (*V)(double, double)){
	
	//fill the matrix with zeros
	for(unsigned int i = 0; i < size; i++)
		for(unsigned int j = 0; j < size; j++)
			A[i][j] = 0;
	
	//set the terms on the three diagonals
	double h1 = 2 / (step * step);     //the main diagonal constant term
	double h2 = - 1 /  (step * step);  //the upper and lower constant terms
	A[0][0] = V(step, omega) + h1;
	A[0][1] = h2;
	for(unsigned int i = 1; i < size - 1; i++){
		A[i][i - 1] = h2;
		A[i][i + 1] = h2;
		A[i][i] = V(step * (i + 1), omega) + h1;
	}
	A[size - 1][size - 1] = V(step * size, omega) + h1;
	A[size - 1][size - 2] = h2;
}



//FUNCTION TO CREATE THE INITIAL VECTORS FOR THE DIMENSIONLES SCHR. EQ.
void generate_vectors(double * A, double * B, unsigned int size, double step, double omega, double (*V)(double, double)){
    
     
    //set the terms on the three diagonals
    double h1 = 2 / (step * step);     //the main diagonal constant term
    double h2 = - 1 /  (step * step);
    for(unsigned int i = 0; i < size; i++)
    {
        A[i] = V(step * (i + 1), omega) + h1;
        B[i] = h2;
    }
}



//FUNCTON TO CREATE AN IDENTITY MATRIX
void generate_identity(double ** A, unsigned int size){
	
	//fill the matrix with zeros
	for(unsigned int i = 0; i < size; i++)
		for(unsigned int j = 0; j < size; j++)
			A[i][j] = 0;
	
	for(unsigned int i = 0; i < size; i++)
		A[i][i] = 1;
}


void find_max_vect(double * A, unsigned int size, unsigned int& a, double& value){
	
	for(unsigned int i = 0; i < size; i++)
		if(A[i] > value){
			value = A[i];
			a = i;
		}
}


void generate_file (double **m, unsigned int size, double step, int index, unsigned int state, double eigenvalue){
	FILE * output;
	stringstream file;
	
	//If it's a gaussian elimination data series
	file << "Data/" << state << ".txt";
	string file_name = file.str();
	const char* file_namec =  file_name.c_str();
	output = fopen(file_namec, "w+");
    fprintf(output, "%d,  %.5lf\n", state, eigenvalue );
	
	//normazlize vectors
	double sum = 0;
	for(unsigned int i = 0; i < size; i++)
		sum += m[i][index] * m[i][index] *  step;
	
	//Generate the file
	for(unsigned int i = 0; i < size; i++)
		fprintf(output, "%.20lf,  %.20lf\n", step * (i + 1), m[i][index] * m[i][index] /sum );
	fclose(output);
}


int main(int argn, char* argv[]) {
	
	//DEFINITION OF MAIN PARAMETERS, MESH SIZE, STEP LENGTH, X-LIMIT
	unsigned int Nl[4] = {100, 150, 200, 250};
    for (int f = 0; f < 4; f++){
        unsigned int N = Nl[f]; 
	double xlim = 7;
	double h = xlim / (N + 1);
	unsigned int n = 5;   //number of eigenvalues to extract

	
	//CONSTRUCTION OF OPERATOR MATRIX
	double * D = new double [N];
    double * U = new double [N];
    generate_vectors(D, U, N, h, 0.1, V_harm);

	
	//CONSTRUCTION OF EIGENVECTORS MATRIX
	double ** R = new double * [N];
	for(unsigned int i = 0; i < N; i++)
		R[i] = new double [N];
	generate_identity(R, N);
	
	/*
	//DIAGONALIZE THE MATRIX WITH THE STANDARD JACOBI METHOD
	double tau, t, c, s;
	double max = 1;
	unsigned int k, l;
	int iter = 0;
	int max_iter = N*N*N;
	while( max >  tolerance && iter < max_iter){
		
		//first we find the largest element in the matrix
		find_max_mat(M, N, &k, &l, &max);
		//we compute the angle needed to stet the largest element to 0
		tau = (M[l][l] - M[k][k]) / (2 * M[k][l]);
		if (tau >= 0) 
			t = 1.0/(tau + sqrt(1.0 + tau * tau));
		else
			t = - 1.0 / (-tau + sqrt(1.0 + tau * tau));
		c = 1 /sqrt(1 + t*t);
		s = c * t;
		
		
		//now we apply the rotation transformation and 
		Jacobi_transform(M, R, N, k, l, s, c);
		iter++;
		printf("%d  %.15lf\n", iter, max);
	}
	
	
	//EIGENVALUES
	unsigned int g;
	double maxv = 0.0;
	
	double * E = new double [n];
	int * idx = new int [n];
	for (unsigned int i = 0; i < n; i++){  //we initialize the vector with the first 5 eigenvalues
		E[i] = M[i][i];
		idx[i] = i;
	}
	for(unsigned int i = 0; i < n; i++)  //we look for the greates value
		if(E[i] > maxv){
			maxv = E[i];
			g = i;
		}
	//extract the 5 smallest eigenvalues by replacing the grates value if one is smaller
	for (unsigned int i = n; i < N; i++)
		if(M[i][i] <= maxv){
			E[g] = M[i][i];
			idx[g] = i;
			maxv = 0;
			for(unsigned int i = 0; i < n; i++)
				if(E[i] > maxv){
					maxv = E[i];
					g = i;
				}
		}
	
	//bubble sort
	int temp_idx;
	double temp_e;
	for(unsigned int j = 0; j < n - 1; j++) 
		for(unsigned int k = 0;  k < n - 1 - j; k++) 
			if(E[k] > E[k + 1]) {
				temp_e = E[k];
				E[k] = E[k + 1];
				E[k + 1] = temp_e;
				temp_idx = idx[k];
				idx[k] = idx[k + 1];
				idx[k + 1] = temp_idx;
			}
	
	for (unsigned int i = 0; i < n; i++){
		printf("%d ->  %.15lf \n", i, E[i]);
		generate_file( R, N, h, idx[i], i, E[i]);
	}
		
	//MEMORY DEALLOCATION AND EXIT STATEMENTS
	for (unsigned int i = 0; i < N; i++){
		delete [] M[i];
		delete [] R[i];
	}
	delete [] M;
	delete [] R;
    */
    
    tqli1(D, U, N, R);
    
    //EIGENVALUES
    unsigned int g;
    double maxv = 0.0;
    
    double * E = new double [n];
    int * idx = new int [n];
    for (unsigned int i = 0; i < n; i++){  //we initialize the vector with the first 5 eigenvalues
        E[i] = D[i];
        idx[i] = i;
    }
    for(unsigned int i = 0; i < n; i++)  //we look for the greates value
        if(E[i] > maxv){
            maxv = E[i];
            g = i;
        }
    //extract the 5 smallest eigenvalues by replacing the grates value if one is smaller
    for (unsigned int i = n; i < N; i++)
        if(D[i] <= maxv){
            E[g] = D[i];
            idx[g] = i;
            maxv = 0;
            for(unsigned int i = 0; i < n; i++)
                if(E[i] > maxv){
                    maxv = E[i];
                    g = i;
                }
        }
    
    //bubble sort
    int temp_idx;
    double temp_e;
    for(unsigned int j = 0; j < n - 1; j++) 
        for(unsigned int k = 0;  k < n - 1 - j; k++) 
            if(E[k] > E[k + 1]) {
                temp_e = E[k];
                E[k] = E[k + 1];
                E[k + 1] = temp_e;
                temp_idx = idx[k];
                idx[k] = idx[k + 1];
                idx[k + 1] = temp_idx;
            }
    
    for (unsigned int i = 0; i < n; i++){
        //printf("%d ->  %.15lf \n", i, E[i]);
        generate_file( R, N, h, idx[i], i, E[i]);
    }
    
    delete [] U;
    delete [] D;
	delete [] E;
	delete [] idx;}
    return 0;
}

