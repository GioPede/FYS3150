#include <iostream>
#include <math.h>
#include <armadillo>
#include <stdlib.h>
#include "mpi.h"

using namespace std;
using namespace arma;
ofstream ofile;
#define soglia 0.5
#define NUM_THREADS 8

double frand()
{
    double f = (double)rand() / RAND_MAX;
    return f;
}
void upmatrix(int n, double **A){
    //setting all the spins upwards
    for(int j=0;j<n;j++){
        for (int i=0;i<n;i++){
            A[i][j]=1.0;
        }
    }
}
void randommatrix(int n, double ** A){
    //setting randomly the lattice
    for(int i=0;i<n;i++){
       for(int j=0;j<n;j++){
          double temp= frand();
          if(1.0-temp>soglia){
          A[i][j]=1.0;
           }
          else A[i][j]=-1.0;
       }
       }
}
int seekingneighboor(int n,int position, int pm){
    //periodic boundary conditions
    int neighboor=0;
    if (pm==1){
        if (position==n-1){
            neighboor = 0;
        }
        else neighboor= position+1;
    }
    else if(pm == -1){
        if (position==0){
           neighboor = n-1;
        }
        else neighboor= position-1;
    }
   return neighboor;
}

void energy(int n, double& E,double& M, double ** s){
    //calculation of the energy
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            M+=s[i][j];
            if(j==n-1){
                if(i==n-1){
                    E -= s[i][j]*(s[0][j]+s[i][0]);
                }
                 else E -= s[i][j]*(s[i+1][j]+s[i][0]);

            }
       else if(i==n-1){
                E -= s[i][j]*(s[0][j]+s[i][j+1]);
            }
            else{ E -= s[i][j]*(s[i+1][j]+s[i][j+1]);}

        }


     }
}

void TEST_EXIT(double* first, double* second, double mcycle){
    //exit of three vectors for plotting with matplotlib afterwards
    fstream fs;
    ofstream ofs;
    ofs.open("../Ising_model/test.txt", std::ofstream::out | std::ofstream::trunc);
    ofs.close();
   fs.open ("../Ising_model/test.txt", fstream::in | fstream::out | fstream::app);

    for(int i=2;i<mcycle;i*=1.5){
                 fs << i << endl;
    }
    fs<<"change"<<endl;
    for(int i=2;i<mcycle;i*=1.5){
                 fs << first[i] << endl;
    }
    fs<<"change"<<endl;
    for(int i=2;i<mcycle;i*=1.5){
             fs << second[i] << endl;
    }
    fs<<"change"<<endl;
    fs<<"break"<<endl;
    fs.close();

}

void metro_giulio(int mcycle,int n, double E,double M, double**s, double * ene, double *magn){
    //setting temperature
    double T=2.4;
    int i=0;
    //start metropolis algorithm
    for(i=0;i<mcycle;i++){
         int x= (int) ((double)(n-1) * frand());
         int y= (int) ((double)(n-1) * frand());
         double deltaM;
         double deltaE=2.0*s[x][y]*(s[x][seekingneighboor(n,y,-1)]+s[x][seekingneighboor(n,y,1)]+s[seekingneighboor(n,x,-1)][y]+s[seekingneighboor(n,x,1)][y]);

         if(frand()<=exp((-deltaE)/T)){
            s[x][y]= s[x][y]*(-1.0);
            deltaM=2.0*s[x][y];
          }
          else{deltaE=0.0;
              deltaM=0.0;
          }
         M+=deltaM;
         E+=deltaE;
         ene[i]+=E;
         magn[i]+=M/(n*n);
         }
}
void EX_C(){

}
int main(int argc, char** argv)
{
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    int numproc;
    MPI_Comm_size(MPI_COMM_WORLD, &numproc);

    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    //initialise random values
    srand(time(NULL));

    //setting size of the matrix
    int n1=20;

    //create matrix
    double ** A=new double*[n1];
    for(int i=0;i<n1;i++){
        A[i]= new double[n1];
    }
    double E1=0;
    double M1=0;

    //setting number of MC cycles
    int mcycle=pow(10,6);

    double *ene= new double[mcycle];
    double *magn=new double[mcycle];
    double *ene_tot= new double[mcycle];
    double *magn_tot=new double[mcycle];

    for(int i=0; i<mcycle;i++){
        ene[i]=0;
        magn[i]=0;
        ene_tot[i]=0;
        magn_tot[i]=0;

    }

    //setting number of repetitions (multiply it for 4 processors in my computer)
    int N_mean=2500;

    for(int i=0; i<N_mean;i++){

        //decide which kind of matrix
        randommatrix(n1,A);
        //upmatrix(n1,A); //setting the matrix
        E1=0;
        M1=0;
        energy(n1,E1,M1,A);
        metro_giulio(mcycle,n1,E1,M1,A,ene,magn);
    }

    //reduce it
    for(int i=0;i<mcycle;i++){
        MPI_Reduce(&ene[i], &ene_tot[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&magn[i], &magn_tot[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }

    if(my_rank == 0){
        for(int i=0; i<mcycle;i++){
            ene_tot[i]/=((double)N_mean*numproc);
            magn_tot[i]= abs(magn[i])/((double)N_mean*numproc);
        }
        TEST_EXIT(ene_tot,magn_tot,mcycle);
    }

    //kill matrix
    for(int i=0;i<n1;i++){
        delete[] A[i];
    }
    delete [] A;

    //kill vectors
    delete [] ene;
    delete [] magn;
    delete [] magn_tot;
    delete [] ene_tot;

    // Finalize the MPI environment.
    MPI_Finalize();

    return 0;
}
