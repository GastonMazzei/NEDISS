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
                  Equation &E,
                  FlowSpecs &Specs,
                  std::vector<double> &RK1,
                  std::vector<double> &RK2,
                  std::vector<double> &RK3,
                  std::vector<double> &RK4,
                  double &answer,
                  double * P);
    void Term1(double t,
               double h,
               double a,
               std::vector<double> &b,
               std::vector<double> &c,
               std::vector<double> &d,
               Equation &E,
               FlowSpecs &Specs,
               std::vector<double> &RK1,
               double * P);
    void Term2(double t,
               double h,
               double a,
               std::vector<double> &b,
               std::vector<double> &c,
               std::vector<double> &d,
               Equation &E,
               FlowSpecs &Specs,
               std::vector<double> &RK1,
               std::vector<double> &RK2,
               double * P);
    void Term3(double t,
               double h,
               double a,
               std::vector<double> &b,
               std::vector<double> &c,
               std::vector<double> &d,
               Equation &E,
               FlowSpecs &Specs,
               std::vector<double> &RK1,
               std::vector<double> &RK2,
               std::vector<double> &RK3,
               double * P);
    void Term4(double t,
               double h,
               double a,
               std::vector<double> &b,
               std::vector<double> &c,
               std::vector<double> &d,
               Equation &E,
               FlowSpecs &Specs,
               std::vector<double> &RK1,
               std::vector<double> &RK2,
               std::vector<double> &RK3,
               std::vector<double> &RK4,
               double * P);
};



template <typename Equation>
void RungeKuttaSolver<Equation>::evolve(double t,
                                     double h,
                                     double a,
                                     std::vector<double> &b,
                                     std::vector<double> &c,
                                     std::vector<double> &d,
                                     Equation &E,
                                     FlowSpecs &Specs,
                                     std::vector<double> &RK1,
                                     std::vector<double> &RK2,
                                     std::vector<double> &RK3,
                                     std::vector<double> &RK4,
                                     double &answer,
                                     double * P){


    answer = a + (*(P+0)) * h * RK1[0] +
                (*(P+1)) * h * RK2[0] +
                (*(P+2)) * h * RK3[0] +
                (*(P+3)) * h * RK4[0];
}

template <typename Equation>
void RungeKuttaSolver<Equation>::Term1(double t,
                                  double h,
                                  double a,
                                  std::vector<double> &b,
                                  std::vector<double> &c,
                                  std::vector<double> &d,
                                  Equation &E,
                                  FlowSpecs &Specs,
                                  std::vector<double> &RK1,
                                  double * P){
    // Initialize auxiliary vectors to any value
    std::vector<double> T1={0},T2={1},T3={0},T4={1};

    // Values defined for Euler.
    T1[0] = 0;
    T2[0] = 1;
    T3[0] = 0;
    T4[0] = 1;
    E.UpdateFlowSpecs(T1,T2,T3,T4,d.size());
    E.Reset();
    E.Field(t,a,b,c,d);
    RK1[0] =  Specs.result; // (*(P+0)) * h

}


template <typename Equation>
void RungeKuttaSolver<Equation>::Term2(double t,
                                  double h,
                                  double a,
                                  std::vector<double> &b,
                                  std::vector<double> &c,
                                  std::vector<double> &d,
                                  Equation &E,
                                  FlowSpecs &Specs,
                                  std::vector<double> &RK1,
                                  std::vector<double> &RK2,
                                  double * P){
    // Initialize auxiliary vectors to any value
    std::vector<double> T1={0},T2={1},T3={0},T4={1};

    // Values defined for Euler.
    // Here RK1 is the first derivative of the other variables :-)
    T1[0] = h * RK1[0] / 2;
    T2[0] = 1;
    T3.resize(RK1.size()-1);
    for (int i=1; i<RK1.size(); ++i) T3[i-1] = h*RK1[i]/2;
    T4[0] = 1;
    E.UpdateFlowSpecs(T1,T2,T3,T4,d.size());
    // Second Taylor Order
    E.Reset();
    E.Field(t+h/2,a,b,c,d);
    RK2[0] =  Specs.result; // (*(P+1)) * h * h / 2
};

template <typename Equation>
void RungeKuttaSolver<Equation>::Term3(double t,
                                  double h,
                                  double a,
                                  std::vector<double> &b,
                                  std::vector<double> &c,
                                  std::vector<double> &d,
                                  Equation &E,
                                  FlowSpecs &Specs,
                                  std::vector<double> &RK1,
                                  std::vector<double> &RK2,
                                  std::vector<double> &RK3,
                                  double * P){

    // Initialize auxiliary vectors to any value
    std::vector<double> T1={0},T2={1},T3={0},T4={1};
    // Values defined for Euler.
    // Here RK1 and RK2 are the first amd second derivatives of the other variables :-)
    T1[0] = h * RK2[0] / 2;
    T2[0] = 1;
    T3.resize(RK2.size()-1);
    for (int i=1; i<RK2.size(); ++i) T3[i-1] = h*RK2[i]/2;
    T4[0] = 1;
    E.UpdateFlowSpecs(T1,T2,T3,T4,d.size());
    // Second Taylor Order
    E.Reset();
    E.Field(t+h/2,a,b,c,d);
    RK3[0] =  Specs.result; // (*(P+2)) * h * h * h / 6
};


template <typename Equation>
void RungeKuttaSolver<Equation>::Term4(double t,
                                  double h,
                                  double a,
                                  std::vector<double> &b,
                                  std::vector<double> &c,
                                  std::vector<double> &d,
                                  Equation &E,
                                  FlowSpecs &Specs,
                                  std::vector<double> &RK1,
                                  std::vector<double> &RK2,
                                  std::vector<double> &RK3,
                                  std::vector<double> &RK4,
                                  double * P){

    // Initialize auxiliary vectors to any value
    std::vector<double> T1={0},T2={1},T3={0},T4={1};
    // Values defined for Euler.
    // Here RK1, RK2 and RK3 are the first, second and third derivatives of the other variables :-)
    T1[0] = h * RK3[0];
    T2[0] = 1;
    T3.resize(RK1.size()-1);
    for (int i=1; i<RK1.size(); ++i) T3[i-1] = h * RK3[i];
    T4[0] = 1;
    E.UpdateFlowSpecs(T1,T2,T3,T4,d.size());
    // Second Taylor Order
    E.Reset();
    E.Field(t+h,a,b,c,d);
    RK4[0] = Specs.result; // (*(P+3)) * h * h * h * h / 24
};



#endif //CPPPROJCT_RUNGEKUTTASOLVER_H
