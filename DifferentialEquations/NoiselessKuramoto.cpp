//
// Created by m4zz31 on 2/11/21.
//

#include "NoiselessKuramoto.h"
#include "GeneralDifferentialEquation.h"
#include <cmath>
#include <cassert>
#include <utility>
#include <vector>

void NoiselessKuramoto::Field(double t, double a, std::vector<double> &b,
                      std::vector<double> &c,
                      std::vector<double> &d){
    // a: central value
    // b: own parameters
    // c: neighbor values
    // d: interactions
    // T3 & T4: linear offset and factor (respect.) to apply to "c",
    // as per generalized Runge Kutta's f(y1,...) replaced with f(a+b*y1,...)
    assert((d.size() == c.size()));
    Specs.result += b[0];
    for (int i = 0; i < b.size(); i++) {
        Specs.result += d[i] * std::sin( (Specs.T3[Specs.j3*i] + Specs.T4[Specs.j4*i] * c[i]) -
                                                 (Specs.T1[0] + Specs.T2[0] * a));
    }
};


void NoiselessKuramoto::d1Field(double t, double a, std::vector<double> &b,
                              std::vector<double> &c,
                              std::vector<double> &d){
    // a: central value
    // b: own parameters
    // c: neighbor values
    // d: interactions
    // T3 & T4: linear offset and factor (respect.) to apply to "c",
    // as per generalized Runge Kutta's f(y1,...) replaced with f(a+b*y1,...)
    printf("First derivative of the NoiselessKuramoto eq. failed to be parametrized");
    exit(1);

};

void NoiselessKuramoto::d2Field(double t, double a, std::vector<double> &b,
                              std::vector<double> &c,
                              std::vector<double> &d){
    // a: central value
    // b: own parameters
    // c: neighbor values
    // d: interactions
    // T3 & T4: linear offset and factor (respect.) to apply to "c",
    // as per generalized Runge Kutta's f(y1,...) replaced with f(a+b*y1,...)
    printf("Second derivative of the NoiselessKuramoto eq. failed to be parametrized");
    exit(1);

};

void NoiselessKuramoto::d3Field(double t, double a, std::vector<double> &b,
                              std::vector<double> &c,
                              std::vector<double> &d){
    // a: central value
    // b: own parameters
    // c: neighbor values
    // d: interactions
    // T3 & T4: linear offset and factor (respect.) to apply to "c",
    // as per generalized Runge Kutta's f(y1,...) replaced with f(a+b*y1,...)
    printf("third derivative of the NoiselessKuramoto eq. failed to be parametrized");
    exit(1);
};