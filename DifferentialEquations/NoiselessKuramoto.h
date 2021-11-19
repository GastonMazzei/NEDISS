//
// Created by m4zz31 on 2/11/21.
//

#ifndef CPPPROJCT_NOISELESSKURAMOTO_H
#define CPPPROJCT_NOISELESSKURAMOTO_H

#include "GeneralDifferentialEquation.h"

//auto GeneralKuramotoFactory(double T1, double T2, double T3, double T4);

class NoiselessKuramoto : public GeneralDifferentialEquation {
public:
    NoiselessKuramoto(): GeneralDifferentialEquation(0){};
    void Field(double t, double a, std::vector<double> &b,
                 std::vector<double> &c,
                 std::vector<double> &d);
};

#endif //CPPPROJCT_NOISELESSKURAMOTO_H
