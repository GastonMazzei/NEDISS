//
// Created by m4zz31 on 20/11/21.
//

#ifndef CPPPROJCT_LINEARTESTEQUATION_H
#define CPPPROJCT_LINEARTESTEQUATION_H


#include "GeneralDifferentialEquation.h"

class LinearTestEquation : public GeneralDifferentialEquation {
public:
    LinearTestEquation(): GeneralDifferentialEquation(0){};
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







#endif //CPPPROJCT_LINEARTESTEQUATION_H
