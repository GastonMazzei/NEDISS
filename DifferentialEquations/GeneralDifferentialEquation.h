//
// Created by m4zz31 on 2/11/21.
//

#ifndef CPPPROJCT_GENERALDIFFERENTIALEQUATION_H
#define CPPPROJCT_GENERALDIFFERENTIALEQUATION_H

#include <functional>
#include <vector>
#include "../Utils/differential_equations_aux.h"

// Typedef to make things easier in other files
#ifndef SCALAR_FLOW_TYPE
#define SCALAR_FLOW_TYPE
typedef std::function<double(double, double, std::vector<double>, std::vector<double>, std::vector<double>)> ScalarFlow;
#endif

// SEE THE .CPP FOR A DISCUSSION ON IMPLEMENTING PERTURBATION FACTORIES ;_)
// or not.. ^.^ we are currently implementing each Diff Eq rather explicitly :p

// The purpose of this class is to potentially support dumb-typing
// PDEs and SDEs and convert them through a general function "BuildForSolver" :-)
class GeneralDifferentialEquation {
public:
    int type; // 0 is ODE, 1 is PDE, 2 could be SDE :-)
    ScalarFlow Field = scalar_flow_placeholder;
    GeneralDifferentialEquation(int type): type(type) {};
    void BuildForSolver();
};

#endif //CPPPROJCT_GENERALDIFFERENTIALEQUATION_H

