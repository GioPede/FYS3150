#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <armadillo>

using namespace std;
using namespace arma;


// The function we are considering
inline double f(double x){
    return atan(x);
}


// Functiond sto compute the derivatives
inline double fwd_derive(double h, double x){
    return (f(x+h)-f(x))/h;
}

inline double bwd_derive(double h, double x){
    return (f(x)-f(x-h))/h;
}

inline double tp_derive(double h, double x){
    return (f(x+h)-f(x-h))/(2*h);
}

inline double second_derive(double h, double x){
    return (f(x+h) + f(x-h) -2*f(x))/(h*h);
}


void print_data (unsigned int size, vec step, vec doub, fvec floa, vec seco, fvec fseco)
{
    FILE * output;
    output = fopen("data.txt", "w+");
    for (unsigned int i = 0; i < size; i++)
        fprintf(output, "%.15lf, %.15lf, %.15lf, %.15lf, %.15lf\n", step(i), doub(i), floa(i), seco(i), fseco(i));
    fclose(output);
}


void print_data_err(unsigned int size, vec step, vec doub, fvec floa, vec err, fvec ferr)
{
    FILE * output;
    output = fopen("data_err.txt", "w+");
    for (unsigned int i = 0; i < size; i++)
        fprintf(output, "%.15lf, %.15lf, %.15lf, %.15lf, %.15lf\n", step(i), doub(i), floa(i), err(i), ferr(i));
    fclose(output);
}

int main(int argc, char* argv[])
{
    int n_steps = atoi(argv[1]);
    vec exps = linspace<vec>(0, (-1)*n_steps, (1 + n_steps)*10);
    vec h;
    h.copy_size(exps);

    for(unsigned int i = 0; i < h.n_elem; i++)
        h(i) = pow(10, exps(i));


    vec tp, second, err;
    fvec ftp, fsecond, ferr;

    tp.copy_size(h);
    second.copy_size(h);
    err.copy_size(h);

    ftp.copy_size(h);
    fsecond.copy_size(h);
    ferr.copy_size(h);

    for (unsigned int i = 0; i < h.n_elem; i++)
    {

        tp(i) = tp_derive(h(i), sqrt(2));
        second(i) = second_derive(h(i), sqrt(2));
        ftp(i) = tp_derive(h(i), sqrt(2));
        fsecond(i) = second_derive(h(i), sqrt(2));
    }

    vec sol;
    fvec fsol;

    sol.copy_size(h);   sol.fill(1.0/3);
    fsol.copy_size(h);  fsol.fill(1.0/3);

    print_data(h.n_elem, h, tp, ftp, second, fsecond);
    print_data_err(h.n_elem, h, tp, ftp, abs((tp-sol)/sol), abs((ftp-fsol)/fsol));
}
