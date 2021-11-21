//
// Created by m4zz31 on 3/11/21.
//

#ifndef CPPPROJCT_RUNGEKUTTASOLVER_H
#define CPPPROJCT_RUNGEKUTTASOLVER_H

#include "GeneralSolver.h"
#include "../macros/macros.h"


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
                  double &answer,
                  double * P);
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
                                     double &answer,
                                     double * P){
                                    // add here a pointer to REF where we will store the vals ;-)
    if (true){
        printf("[WARNING] Runge Kutta is currently not implemented. This is not fatal but just returns integration placehholder\n");
        std::cout << std::flush;
        answer = INTEGRATION_PLACEHOLDER;
    } else {
        printf("[FATAL] RUNGE KUTTA is not available. Pending communication between processor for vector K1,2,3,4");
        std::cout << std::flush;
        exit(1);
    }
//    //-------------------------------------\
//    // k1 = f(tn,yn)                       |
//    // k2 = f(tn + h / 2, yn + h * k1 / 2) |
//    // k3 = f(tn + h / 2, yn + h * k2 / 2) |
//    // k4 = f(tn + h, yn + h * k3)         |
//    //------------------------------------/
//    // answer = a + h * (P0*k1 + P1*k2 + ...)
//    // an example in wikipedia shows 1/6 (k1 + 2k2 + 2k3 + k4)
//
//    // Initialize auxiliary vectors
//    std::vector<double> T1={0},T2={1},T3={0},T4={1}; // Any val.
//    answer = a; // start building the answer
//
//    double k1,k2,k3,k4;
//
//    // Compute k1 and add it to the answer
//    k1 = 0;
//    T1[0] = 0;
//    T2[0] = 1;
//    T3[0] = 0;
//    T4[0] = 1;
//    Specs.result = 0;
//    E.UpdateFlowSpecs(T1,T2,T3,T4,d.size());
//    E.Field(t,a,b,c,d);
//    k1 = Specs.result;
//    answer += h * (*P) * k1;
//
//    // Try to avoid computing terms we wont use: K2
//    if (((*(P+3) == 0) && (*(P+2) == 0)) && (*(P+1) == 0)){
//        k2 = 0;
//    } else {
//        // Compute k2 and update the answer
//        T1[0] = h / 2 * k1;
//        T2[0] = 1;
//        T3[0] = 0;
//        T4[0] = 1;
//        Specs.result = 0;
//        E.UpdateFlowSpecs(T1,T2,T3,T4,d.size());
//        E.Field(t+h/2,a,b,c,d);
//        k2 = Specs.result;
//        answer += h * (*(P+1)) * k2;
//    }
//
//    // Try to avoid computing terms we wont use: K3
//    if ((*(P+3) == 0) && (*(P+2) == 0)){
//        k3 = 0;
//    } else {
//        // Compute k3 and update the answer
//        T1[0] = h / 2 * k2;
//        T2[0] = 1;
//        T3[0] = 0;
//        T4[0] = 1;
//        Specs.result = 0;
//        E.UpdateFlowSpecs(T1,T2,T3,T4,d.size());
//        E.Field(t+h/2,a,b,c,d);
//        k3 = Specs.result;
//        answer += h * (*(P+2)) * k3;
//    }
//
//    // Try to avoid computing terms we wont use: K4
//    if (*(P+3) == 0){
//        k4 = 0;
//    } else {
//        // Compute k4 and update the answer
//        T1[0] = h * k3;
//        T2[0] = 1;
//        T3[0] = 0;
//        T4[0] = 1;
////        Specs.result = 0;
////        E.UpdateFlowSpecs(T1,T2,T3,T4,d.size());
////        E.Field(t,a,b,c,d);
////        k4 = Specs.result;
//        answer += h * (*(P+3)) * k4;
//    }
//
//    if (*(P+0) != 0) {
//        // Initialize result
//        T1[0] = 0;
//        T2[0] = 1;
//        T3[0] = 0;
//        T4[0] = 1;
//        Specs.result = 0;
//        E.UpdateFlowSpecs(T1,T2,T3,T4,d.size());
//        E.Field(t,a,b,c,d);
//        // Populate answer step 1 of 4 :-)
//        // it is the reslt weighted by *P
//        answer += (T1[0] + T2[0] * a + h * Specs.result) * (*P);
//    }
//    if (*(P+1) != 0) {
//    }
//    if (*(P+2) != 0) {
//    }
//    if (*(P+3) != 0) {
//    }

}

#endif //CPPPROJCT_RUNGEKUTTASOLVER_H
