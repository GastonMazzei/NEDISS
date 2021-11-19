//
// Created by m4zz31 on 2/11/21.
//

#ifndef CPPPROJCT_GENERALDIFFERENTIALEQUATION_H
#define CPPPROJCT_GENERALDIFFERENTIALEQUATION_H

#include <functional>
#include <vector>
#include "../Utils/differential_equations_aux.h"

// Typedef to make things easier in other files
//#ifndef SCALAR_FLOW_TYPE
//#define SCALAR_FLOW_TYPE
//typedef std::function<double(double, double, std::vector<double>, std::vector<double>, std::vector<double>)> ScalarFlow;
//#endif


struct FlowSpecs{
    std::vector<double> T1;
    std::vector<double> T2;
    std::vector<double> T3;
    std::vector<double> T4;
    int N, j1, j2, j3, j4;
    double result;
    FlowSpecs(){};
    FlowSpecs(std::vector<double> T1,
              std::vector<double> T2,
              std::vector<double> T3,
              std::vector<double> T4,
              int N): T1(std::move(T1)), T2(std::move(T2)), T3(std::move(T3)), T4(std::move(T4)), N(N){};
};

// SEE THE .CPP FOR A DISCUSSION ON IMPLEMENTING PERTURBATION FACTORIES ;_)
// or not.. ^.^ we are currently implementing each Diff Eq rather explicitly :p


// The purpose of this class is to potentially support dumb-typing
// PDEs and SDEs and convert them through a general function "BuildForSolver" :-)
class GeneralDifferentialEquation {
public:
    int type; // 0 is ODE, 1 is PDE, 2 could be SDE :-)
    bool RequiresBuilding = false; // true  if it requires building :-)
    FlowSpecs Specs;
    GeneralDifferentialEquation(int type): type(type) {};

    void UpdateFlowSpecs(std::vector<double> &T1,
                                 std::vector<double> &T2,
                                 std::vector<double> &T3,
                                 std::vector<double> &T4,
                                 int N);
    void Reset(); // reset result to 0
    void BuildForSolver();
};



#endif //CPPPROJCT_GENERALDIFFERENTIALEQUATION_H

