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
    double evolve(double t,
                  double h,
                  double a,
                  std::vector<double> &b,
                  std::vector<double> &c,
                  std::vector<double> &d,
                  Equation &E,//auto F,// ScalarFlow &F,
                  FlowSpecs &Specs);
};



template <typename Equation>
double EulerSolver<Equation>::evolve(double t,
                           double h,
                           double a,
                           std::vector<double> &b,
                           std::vector<double> &c,
                           std::vector<double> &d,
                           Equation &E,//ScalarFlow &F,
                           FlowSpecs &Specs){

    return 1.0;
}

#endif //CPPPROJCT_EULERSOLVER_H
