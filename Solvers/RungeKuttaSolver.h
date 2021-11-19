//
// Created by m4zz31 on 3/11/21.
//

#ifndef CPPPROJCT_RUNGEKUTTASOLVER_H
#define CPPPROJCT_RUNGEKUTTASOLVER_H

#include "GeneralSolver.h"


template <typename Equation>
class RungeKuttaSolver {
public:
    RungeKuttaSolver(){};
    void evolve(double t,
                  double h,
                  double a,
                  std::vector<double> &b,
                  std::vector<double> &c,
                  std::vector<double> &d,
                  Equation &E,//auto F,// ScalarFlow &F,
                  FlowSpecs &Specs,
                  double &answer);
};



template <typename Equation>
void RungeKuttaSolver<Equation>::evolve(double t,
                                     double h,
                                     double a,
                                     std::vector<double> &b,
                                     std::vector<double> &c,
                                     std::vector<double> &d,
                                     Equation &E,//ScalarFlow &F,
                                     FlowSpecs &Specs,
                                     double &answer){

    answer = 1.1;
}


#endif //CPPPROJCT_RUNGEKUTTASOLVER_H
