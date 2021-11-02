//
// Created by m4zz31 on 1/11/21.
//

#ifndef CPPPROJCT_GENERALSOLVERS_H
#define CPPPROJCT_GENERALSOLVERS_H

#include <functional>
#include <vector>

// Typedef to make things easier in other files
//#ifndef SOLVER_TYPE
//#define SOLVER_TYPE
typedef std::function<double(double, std::vector<double>, std::vector<double>, std::vector<double>)> Solver;
//#endif

// Factory that is a thin-wrapper for general solvers
template <typename T>
Solver GeneralSolverFactory();

// Selector that instantiates specific solvers from the factory
Solver SolverSelector(std::string name);

class GeneralSolver{
public:
    double h = 0.01;
    double t = 0;
    GeneralSolver(){};
};

#endif //CPPPROJCT_GENERALSOLVERS_H
