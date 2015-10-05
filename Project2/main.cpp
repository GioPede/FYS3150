#include <iostream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <string>   
#include <sstream>


using namespace std;

inline double V_harm(double x, double useless) {
	return x*x;
}


inline double V_2e(double x, double omega) {
	return x*x* omega + 1.0/x;
}



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


//FUNCTON TO CREATE AN IDENTITY MATRIX
void generate_identity(double ** A, unsigned int size){
	
	//fill the matrix with zeros
	for(unsigned int i = 0; i < size; i++)
		for(unsigned int j = 0; j < size; j++)
			A[i][j] = 0;
	
	for(unsigned int i = 0; i < size; i++)
		A[i][i] = 1;
}


//FUNCTION TO COMPUTE THE SUM OF ALL NON DIAGONAL TERMS
double off (double ** A, unsigned int size){
	
	double sum = 0.0;
	for(unsigned int i = 0; i < size; i++)
		for(unsigned int j = 0; j < size; j++)
			if(i != j)
				sum += abs(A[i][j]);
	
	return sum;
}


//FUNCTION TO FIND THE BIGGEST NON DIAGONAL ELEMENT OF THE MATRIX. 
//RETURNS THE INDEXES AND THE VALUE OF THE ELEMENT.
void find_max_mat(double ** A, unsigned int size, unsigned int *a, unsigned int *b, double *value){

	*value = 0.0;
	for(unsigned int i = 0; i < size; i++)
		for(unsigned int j = 0; j < size; j++)
			if(i != j && abs(A[i][j]) > *value){
				*value = abs(A[i][j]);
				*a = i;
				*b = j;
			}
}


//GIVEN A MATRIX, A SINE AND A COSINE APPLIES A ROTATION TO THE MATRIX
void Jacobi_transform(double ** A, double ** EV, unsigned int size, unsigned int k, unsigned int l, double sin, double cos){
	
	double kk, ll, ik, il;
	kk = A[k][k];
	ll = A[l][l];
	
	//the diagonal values
	A[k][k] = cos * cos * kk - 2.0 * cos * sin * A[k][l] + sin * sin * ll;
	A[l][l] = sin * sin * kk + 2.0 * cos * sin * A[k][l] + cos * cos * ll;
	
	//the pivot values 
	A[k][l] = 0.0;
	A[l][k] = 0.0;  
	for (unsigned int i = 0; i < size; i++ ){
		if ( i != k && i != l ) {
			ik = A[i][k];
			il = A[i][l];
			A[i][k] = cos * ik - sin * il;
			A[k][i] = A[i][k];
			A[i][l] = cos * il + sin * ik;
			A[l][i] = A[i][l];
		}
		
		//we also change the eigenvectors
		double ik = EV[i][k];
		double il = EV[i][l];
		EV[i][k] = cos * ik - sin * il;
		EV[i][l] = cos * il + sin * ik;
	}
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
	
    double * r2 = new double [size]; 
    double * psi2 = new double [size];
    
	//If it's a gaussian elimination data series
	file << "Data/" << state << ".txt";
	string file_name = file.str();
	const char* file_namec =  file_name.c_str();
	output = fopen(file_namec, "w+");
    fprintf(output, "%d,  %.5lf\n", state, eigenvalue );
	
	//normazlize vectors
	double sum = 0;
	for(unsigned int i = 0; i < size; i++){
        r2[i] = step * (i + 1) * step * (i + 1);
        psi2[i] = m[i][index] * m[i][index] / r2[i];
		sum += psi2[i] *  step;
    }
    
	//Generate the file
	for(unsigned int i = 0; i < size; i++)
		fprintf(output, "%.20lf,  %.20lf\n", step * (i + 1), psi2[i] / sum );
	fclose(output);
}


int main(int argn, char* argv[]) {
	
	//DEFINITION OF MAIN PARAMETERS, MESH SIZE, STEP LENGTH, X-LIMIT
	unsigned int N = 250;
	double xlim = 10;
	double h = xlim / (N + 1);
	double tolerance = pow(10, -1);
	unsigned int n = 3;   //number of eigenvalues to extract

	
	//CONSTRUCTION OF OPERATOR MATRIX
	double ** M = new double * [N];
	for(unsigned int i = 0; i < N; i++)
		M[i] = new double [N];
	generate_matrix(M, N, h, 0.1, V_2e);

	
	//CONSTRUCTION OF EIGENVECTORS MATRIX
	double ** R = new double * [N];
	for(unsigned int i = 0; i < N; i++)
		R[i] = new double [N];
	generate_identity(R, N);
	
	
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
		//printf("%d  %.15lf\n", iter, max);
	}
	
	printf("%d \n", iter);
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
	//extract the n smallest eigenvalues by replacing the grates value if one is smaller
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
	delete [] E;
	delete [] idx;
    return 0;
}

