//
// Created by m4zz31 on 2/11/21.
//

#include <cassert>

#include "GeneralDifferentialEquation.h"

#include "../Utils/error.h"


void GeneralDifferentialEquation::BuildForSolver(){
    if (type != 0){
        error_report("Engine only supports scalar flow types! :-)");
    }
}


void GeneralDifferentialEquation::Reset(){
    Specs.result = 0.;
}


void GeneralDifferentialEquation::UpdateFlowSpecs(std::vector<double> &T1,
                                        std::vector<double> &T2,
                                        std::vector<double> &T3,
                                        std::vector<double> &T4,
                                        int N){
    Specs.N = N;
    Specs.T1 = T1;
    Specs.T2 = T2;
    Specs.T3 = T3;
    Specs.T4 = T4;
    Specs.result = 0.;
    Specs.j1 = 0;
    Specs.j2 = 0;
    Specs.j3 = 0;
    Specs.j4 = 0;
    if (Specs.T1.size() != 1) {
        assert(Specs.T1.size() == Specs.N);
        Specs.j1 = 1;
    }
    if (Specs.T1.size() != 1) {
        assert(Specs.T1.size() == Specs.N);
        Specs.j2 = 1;
    }
    if (Specs.T1.size() != 1) {
        assert(Specs.T1.size() == Specs.N);
        Specs.j3 = 1;
    }
    if (Specs.T1.size() != 1) {
        assert(Specs.T1.size() == Specs.N);
        Specs.j4 = 1;
    }
};


