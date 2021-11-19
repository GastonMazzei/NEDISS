//
// Created by m4zz31 on 2/11/21.
//

#include "GeneralDifferentialEquation.h"
#include "../Utils/error.h"
#include <cassert>

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


// EXAMPLE OF A GENERAL PERTURBATION FACTORY :-)
//
// WE ARE CURRENTLY IMPLEMENTING IT EXPLICITLY FOR EACH DIFFERENTIAL EQUATION :-)
//
//auto GeneralKuramotoFactory(std::vector<double> &T1,
//                            std::vector<double> &T2,
//                            std::vector<double> &T3,
//                            std::vector<double> &T4) {
//    //
//    // General template for applying a perturbation to a field :-)
//    //
//    // This is the case for the Noiseless Kuramoto equation :-)
//    //
//    // General flows cant be accepted as we would need to pass by copy
//    // instead of reference, which would only be acceptable if we change
//    // vectors to arrays, OR actively manage memory here, which is not
//    // necessarily well behaved and its even potentially leaky :-(
//    //
//    // Instead of calling F(t,y)
//    // we can call
//    // F(T1 + T2 * t, T3 + T4 * y)
//    //
//    return [&] (double t, double a, std::vector<double> &b,
//                std::vector<double> &c,
//                std::vector<double> &d) -> double
//    {
//        //create temporally
//        //
//        // c = c + T2, etc...
//        return F(//variables go here);
//    };
//};
//
// EXAMPLE OF A FLOW :-)
//
//double KuramotoKernel (double t, double a, std::vector<double> &b,
//                       std::vector<double> &c,
//                       std::vector<double> &d) {
//    // a: central value
//    // b: own parameters
//    // c: neighbor values
//    // d: interactions
//    double result;
//    result = b[0];
//    for (int i = 0; i < b.size(); i++) {
//        result += d[i] * std::sin(c[i] - a);
//    }
//    return result;
//};
