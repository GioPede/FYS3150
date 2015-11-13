#include <iostream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <string>   
#include <sstream>
#include <fstream>
#include <iomanip>
#include <omp.h>
#include <random>
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include <unistd.h>

#define NUM_THREADS 8

short ** M;
short ** E;
int N;
double dz[17];

using namespace std;

inline double ran(){
	//return ((double) generator())/2147483647;
	return ((double) rand()) / RAND_MAX;
}

void boundary(int N){
	for (int i = 1; i <= N; i++){
		M[0][i] = M[N][i];
		M[N + 1][i] = M[1][i];
		M[i][0] = M[i][N];
		M[i][N + 1] = M[i][1];
	}
}

void generate_random(int N){
	srand(time(NULL));
	for (int i = 1; i <= N; i++)
		for (int j = 1; j <= N; j++)
			//(ran() > 0.5 ? M[i][j] = 1 : M[i][j] = -1) ;
			M[i][j] = 1;
}

int t, dt, NS; 
double energy, denergy, st;


void renderScene(void) {
	
	if(t % dt == 0){
		float T = 3.0 - (float) t / (dt*10);
		for( int i =-8; i <= 8; i++) dz[i+8] = 0;
		for( int i =-8; i <= 8; i+=4) dz[i+8] = exp(-i/(T));
		
		energy  /= NS;
		denergy /= NS;
		st = energy*energy / denergy;
		printf("%.10f       %.10f       %.10f\n", T, energy, st);
		NS = 0;
		energy = 0;
		denergy = 0;
	}
	t++;
	NS++;
	int dE, i, j;
	int n;
	#pragma omp parallel for private (dE, i, j)
	for (n = 0; n <= N*N; n++){
			i = (int) (ran()*N) + 1;
			j = (int) (ran()*N) + 1;
			dE = 2 * M[i][j] * (M[i][j+1] + M[i][j-1] + M[i+1][j] + M[i-1][j]);
			if(dz[dE+8] > ran()){
				M[i][j] *= -1;
			}
			boundary(N);
	}
	
	float step = 2.0/N;
	
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
	for (int i = 1; i <= N; i++){
		for (int j = 1; j <= N; j++){
			glBegin(GL_QUADS);
				(M[i][j] == 1 ? glColor3f(1.0f, 0.0f, 0.0f) : glColor3f(0.0f, 0.0f, 1.0f) );
				glVertex2f(-1.0f + step * (i - 1), -1.0f + step * (j - 1));
				glVertex2f(-1.0f + step * (i) , -1.0f + step * (j - 1));
				glVertex2f(-1.0f + step * (i), -1.0f + step * (j));
				glVertex2f(-1.0f + step * (i - 1), -1.0f + step * (j));
			glEnd();
		}
	}
	glutSwapBuffers();
	glutPostRedisplay();
}



int main(int argc, char **argv){
    omp_set_num_threads(NUM_THREADS);  
	
	// init GLUT and create Window
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowPosition(100,100);
	glutInitWindowSize(450,450);
	glutCreateWindow("GLUT");

	// register callbacks
	glutDisplayFunc(renderScene);
	
	N = 150; // SIZE OF THE MATRIX
	
	M = new short * [N + 2];
	for (int i = 0; i < N + 2; i++)
		M[i] = new short [N + 2];
	
	E = new short * [5];
	for (int i = 0; i < 5; i++)
		E[i] = new short [5];
	
	t = 0;
	dt = 200;
	srand(time(NULL));
	generate_random(N);
	glutMainLoop();
	
	
	for (int i = 0; i < N + 2; i++)
		delete [] M[i];
	delete [] M;
}