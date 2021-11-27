//
// Created by m4zz31 on 2/11/21.
//

#ifndef CPPPROJCT_NOISELESSKURAMOTO_H
#define CPPPROJCT_NOISELESSKURAMOTO_H

#include "GeneralDifferentialEquation.h"

class NoiselessKuramoto : public GeneralDifferentialEquation {
public:
    NoiselessKuramoto(): GeneralDifferentialEquation(0){};

    bool requiresCom(int d);

    void Field(double t, double a, std::vector<double> &b,
                 std::vector<double> &c,
                 std::vector<double> &d);

    void d1Field(double t, double a, std::vector<double> &b,
               std::vector<double> &c,
               std::vector<double> &d);

    void d2Field(double t, double a, std::vector<double> &b,
                 std::vector<double> &c,
                 std::vector<double> &d);

    void d3Field(double t, double a, std::vector<double> &b,
                 std::vector<double> &c,
                 std::vector<double> &d);
};

#endif //CPPPROJCT_NOISELESSKURAMOTO_H
