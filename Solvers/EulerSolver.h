//
// Created by m4zz31 on 2/11/21.
//

#ifndef CPPPROJCT_EULERSOLVER_H
#define CPPPROJCT_EULERSOLVER_H

#include "GeneralSolvers.h"

class EulerSolver : public GeneralSolver {
public:
    int deg = 1;
    EulerSolver();
    EulerSolver(int deg): deg(deg){};
    double evolve(double a, std::vector<double> &b, std::vector<double> &c, std::vector<double> &d);
};


#endif //CPPPROJCT_EULERSOLVER_H
