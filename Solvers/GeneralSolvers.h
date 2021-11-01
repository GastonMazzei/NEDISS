//
// Created by m4zz31 on 1/11/21.
//

#ifndef CPPPROJCT_GENERALSOLVERS_H
#define CPPPROJCT_GENERALSOLVERS_H

#include <functional>
#include <vector>

// Typedef to make things easier in other files
typedef std::function<double(double, std::vector<double>, std::vector<double>, std::vector<double>)> Solver;

// Factory that is a thin-wrapper for general solvers
template <typename T>
Solver GeneralSolverFactory();

// Selector that instantiates specific solvers from the factory
Solver SolverSelector(std::string name);

class EulerSolver {
public:
    double evolve(double a, std::vector<double> b, std::vector<double> c, std::vector<double> d);
};

#endif //CPPPROJCT_GENERALSOLVERS_H
