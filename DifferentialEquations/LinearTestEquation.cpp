//
// Created by m4zz31 on 20/11/21.
//
#include "LinearTestEquation.h"
#include <cmath>
#include <cassert>
#include <utility>
#include <vector>
#include <iostream>

void LinearTestEquation::Field(double t, double a, std::vector<double> &b,
                              std::vector<double> &c,
                              std::vector<double> &d){
    // a: central value
    // b: own parameters, i.e.
    //  dx/dt = b[0] * x + b[1]
    // T1 & T2 linear offset and factor (respectively) to apply to "a"
    if (b.size()==2) {
        Specs.result += b[0] * (Specs.T1[0] + Specs.T2[0] * a) + b[1];
    } else if (b.size() == 1){
        Specs.result += b[0] * (Specs.T1[0] + Specs.T2[0] * a);
    } else {
        printf("[FATAL] LinearTestEquation accepts only 1 and 2 dimensional parameters!\n");
        std::cout<<std::flush;
        exit(1);
    }
};


void LinearTestEquation::d1Field(double t, double a, std::vector<double> &b,
                                std::vector<double> &c,
                                std::vector<double> &d){
    // a: central value
    // b: own parameters, i.e.
    //  dx/dt = b[0] * x + b[1]
    // T1 & T2 linear offset and factor (respectively) to apply to "a"
    if (b.size()==2) {
        Specs.result += b[0] * (b[0] * (Specs.T1[0] + Specs.T2[0] * a) + b[1]);
    } else if (b.size() == 1){
        Specs.result += b[0] * b[0] * (Specs.T1[0] + Specs.T2[0] * a);
    } else {
        printf("[FATAL] LinearTestEquation accepts only 1 and 2 dimensional parameters!\n");
        std::cout<<std::flush;
        exit(1);
    }
};

void LinearTestEquation::d2Field(double t, double a, std::vector<double> &b,
                                std::vector<double> &c,
                                std::vector<double> &d){
    // a: central value
    // b: own parameters
    // c: neighbor values
    // d: interactions
    // T3 & T4: linear offset and factor (respect.) to apply to "c",
    // as per generalized Runge Kutta's f(y1,...) replaced with f(a+b*y1,...)
    if (b.size()==2) {
        Specs.result += b[0] * b[0] * (b[0] * (Specs.T1[0] + Specs.T2[0] * a) + b[1]);
    } else if (b.size() == 1){
        Specs.result += b[0] * b[0] * b[0] * (Specs.T1[0] + Specs.T2[0] * a);
    } else {
        printf("[FATAL] LinearTestEquation accepts only 1 and 2 dimensional parameters!\n");
        std::cout<<std::flush;
        exit(1);
    }
};

void LinearTestEquation::d3Field(double t, double a, std::vector<double> &b,
                                 std::vector<double> &c,
                                 std::vector<double> &d){
    // a: central value
    // b: own parameters
    // c: neighbor values
    // d: interactions
    // T3 & T4: linear offset and factor (respect.) to apply to "c",
    // as per generalized Runge Kutta's f(y1,...) replaced with f(a+b*y1,...)
    if (b.size()==2) {
        Specs.result += b[0] * b[0] * b[0] * (b[0] * (Specs.T1[0] + Specs.T2[0] * a) + b[1]);
    } else if (b.size() == 1){
        Specs.result += b[0] * b[0] * b[0] * b[0] * (Specs.T1[0] + Specs.T2[0] * a);
    } else {
        printf("[FATAL] LinearTestEquation accepts only 1 and 2 dimensional parameters!\n");
        std::cout<<std::flush;
        exit(1);
    }
};