//
// Created by m4zz31 on 2/11/21.
//

#ifndef CPPPROJCT_EULERSOLVER_H
#define CPPPROJCT_EULERSOLVER_H

#include "GeneralSolver.h"

template <typename Equation>
class EulerSolver {
public:
    EulerSolver(){};
    void evolve(double t,
                  double h,
                  double a,
                  std::vector<double> &b,
                  std::vector<double> &c,
                  std::vector<double> &d,
                  Equation &E,//auto F,// ScalarFlow &F,
                  FlowSpecs &Specs,
                  double &answer,
                  double * P);
};



template <typename Equation>
void EulerSolver<Equation>::evolve(double t,
                           double h,
                           double a,
                           std::vector<double> &b,
                           std::vector<double> &c,
                           std::vector<double> &d,
                           Equation &E,//ScalarFlow &F,
                           FlowSpecs &Specs,
                           double &answer,
                           double * P){

    // OVERWRITE: testing number of errors!
//    answer = a+1;
//    return;
    // -------END--OF--DEBUGGIN!

    // Initialize auxiliary vectors
    std::vector<double> T1={0},T2={1},T3={0},T4={1}; // Any val.
    answer = 0;
    if (*(P+0) != 0) {
        // Initialize result
        T1[0] = 0;
        T2[0] = 1;
        T3[0] = 0;
        T4[0] = 1;
        Specs.result = 0;
        E.UpdateFlowSpecs(T1,T2,T3,T4,d.size());
        E.Field(t,a,b,c,d);
        // Populate answer step 1 of 4 :-)
        // it is the reslt weighted by *P
        answer += (T1[0] + T2[0] * a + h * Specs.result) * (*P);
    }
    if (*(P+1) != 0) {
        printf("The current equation does not have method for their derivatives, Euler method fails!\n");
        std::cout<<std::flush;
        exit(1);
    }
    if (*(P+2) != 0) {
        printf("The current equation does not have method for their derivatives, Euler method fails!\n");
        std::cout<<std::flush;
        exit(1);
    }
    if (*(P+3) != 0) {
        printf("The current equation does not have method for their derivatives, Euler method fails!\n");
        std::cout<<std::flush;
        exit(1);
    }

}

#endif //CPPPROJCT_EULERSOLVER_H
