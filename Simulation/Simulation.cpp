//
// Created by m4zz31 on 29/10/21.
//
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

//#include <Eigen/Dense>
//using namespace Eigen;


void toy_simulation(int argc, char* argv[]){

    // Ring size initialization
    const int N = 20;

    // Open files
    ofstream cfiles[N];
    for (int  i=0; i<N; i++){
        cfiles[i].open("osc" + to_string(i) + ".txt");
    }

    // Time initialization
    double t_0 = 0;
    double t_max = 1;
    double dt=1e-3;
    double t = t_0;
    const int Nt = (int) (1/dt) + 1;

    // Iterator initialization
    int i,j,k;

    // Container initialization
    double B[N] = {0};
    double W[N] = {1};
    double A[N][Nt] ;

    // Define eqs. parameters
    double K = 3.;
    if (argc>1){
        //K = stod (argv[0], nullptr);
        K = atof(argv[1]);
        cout << "Recieved K: " << K << endl;
    }
    const double K_by_N = K / (double) N;
    double dW = 0.01;
    double TH_0;
    TH_0 = (double) M_PI / (double) 6;
    double sum[N];
    for (i=0;i<N;i++){
        W[i] = W[i] + i * dW;
        A[i][0] = M_PI * (1/ (double) N) * i + TH_0;
    }


    // eq. is:
    //
    // dA_i/dt = w_i + K/N * sum_j(A_j - A_i);
    //
    // therefore the code should be
    // B_i(t) = (w_i + K/N * sum_j(A_j-A_i))
    // A_i(t+dt) = dt * B_i(t) + A_i(t)

    // Start writing and integrating!
    for (j=0; j<Nt-1; j++) {
# pragma omp for
        for (i=0; i<N; i++) {
            for (k=0; k<N; k++){
                sum[i] += sin(A[k][j] - A[i][j]);
            }
            B[i] = W[i] + K_by_N  * sum[i];
        }
        for (i=0; i<N; i++) {
            A[i][j+1] = A[i][j] + dt * B[i];
            cfiles[i] << A[i][j] << "," ;
        }
        t += dt;
    }

    for (int  i=0; i<N; i++){
        cfiles[i].close();
    }



}
