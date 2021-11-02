//
// Created by m4zz31 on 1/11/21.
//

#include  "GeneralSolvers.h"
#include "EulerSolver.h"
#include "../Utils/error.h"

template <typename T>
Solver GeneralSolverFactory(){
    T o;
    return [&](double a, std::vector<double> b, std::vector<double> c, std::vector<double> d) -> double
    {
        return o.evolve(a, b, c, d);
    };
};


Solver SolverSelector(std::string name){
    // Types available should be:
    //  (1) Euler ('eul')
    //  (2) Runge-Kutta-4 ('rk4')

    // for the time being
    if (name != "eul") {
        // raise error here
        error_report("Only support for Euler Integration is available :-(");
    }
    return GeneralSolverFactory<EulerSolver>();
};

