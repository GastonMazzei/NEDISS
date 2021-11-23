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
                  Equation &E,
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
                           Equation &E,
                           FlowSpecs &Specs,
                           double &answer,
                           double * P){


    // Initialize auxiliary vectors to any value
    std::vector<double> T1={0},T2={1},T3={0},T4={1};

    // Values defined for Euler.
    T1[0] = 0;
    T2[0] = 1;
    T3[0] = 0;
    T4[0] = 1;
    E.UpdateFlowSpecs(T1,T2,T3,T4,d.size());

    // Zero Taylor Order
    answer = a;

    // First Taylor Order
    if (*(P+0) != 0) {
        Specs.result = 0;
	E.Reset();
        E.Field(t,a,b,c,d);
        answer +=  h * Specs.result;
    }

    // Second Taylor Order
    if (*(P+1) != 0) {
        Specs.result = 0;
	E.Reset();
        E.d1Field(t,a,b,c,d);
        answer +=  h * h * Specs.result  / 2;
    }

    // Third Taylor Order
    if (*(P+2) != 0) {
        Specs.result = 0;
	E.Reset();
        E.d2Field(t,a,b,c,d);
        answer +=  h * h * h * Specs.result  / 6;    
    }

    // Fourth Taylor Order
    if (*(P+3) != 0) {
        Specs.result = 0;
	E.Reset();
        E.d3Field(t,a,b,c,d);
        answer +=  h * h * h * h * Specs.result / 24;
    }
}

#endif //CPPPROJCT_EULERSOLVER_H
