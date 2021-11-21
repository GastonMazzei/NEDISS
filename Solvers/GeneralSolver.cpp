//
// Created by m4zz31 on 1/11/21.
//

#include  "GeneralSolver.h"


#include "EulerSolver.h"
#include "RungeKuttaSolver.h"
#include "../Utils/error.h"

SolverConfig::SolverConfig(){
    // 0 is Euler: ask for the depth
    // 1 is Runge-Kutta: ask for the coefficients
    s=std::stoi(std::getenv("SOLVER"));
    if (s == 0){
        d = std::stoi(std::getenv("SOLVERDEPTH"));
    } else if (s == 1) {
        P[0] = (double) std::stof(std::getenv("K1"));
        P[1] = (double) std::stof(std::getenv("K2"));
        P[2] = (double) std::stof(std::getenv("K3"));
        P[3] = (double) std::stof(std::getenv("K4"));
    }
}
