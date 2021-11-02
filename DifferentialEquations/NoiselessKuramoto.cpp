//
// Created by m4zz31 on 2/11/21.
//

#include "NoiselessKuramoto.h"
#include "GeneralDifferentialEquation.h"
#include <cmath>

#include <utility>
#include <vector>

double NoiselessKuramoto::Field(double t, double a, std::vector<double> &b,
                      std::vector<double> &c,
                      std::vector<double> &d){
    // a: central value
    // b: own parameters
    // c: neighbor values
    // d: interactions
    // T3 & T4: linear offset and factor (respect.) to apply to "c"
    // (we are asking they have the same length... is it too much to ask?)
    Specs.result = b[0];
    for (int i = 0; i < b.size(); i++) {
        Specs.result += d[i] * std::sin( (Specs.T3[Specs.j3*i] + Specs.T4[Specs.j4*i] * c[i])- a);
    }
    return Specs.result;
};
