//
// Created by m4zz31 on 20/11/21.
//

#ifndef CPPPROJCT_LINEARTESTEQUATION_H
#define CPPPROJCT_LINEARTESTEQUATION_H


#include "GeneralDifferentialEquation.h"

// df/dx = a * f + b, where the sol is f(x) = b + c * e^(ax)
// it is useful for tests as it does not carry interactions: all nodes
// should produce the same result if given the same initial condition :-)
class LinearTestEquation : public GeneralDifferentialEquation {
public:
    LinearTestEquation(): GeneralDifferentialEquation(0){};

    // for compliance with other equations
    // which do have a non-constant behaviour
    bool requiresCom(int d);

    void Field(double t, double a, 
	       std::vector<double> &b,
               std::vector<double> &c,
               std::vector<double> &d);

    // For this equation, the first (second,third,...) derivatives
    // have a pre-established form that does not require targeting the
    // Field of neighbors, so we write it.
    void d1Field(double t, double a, 
		 std::vector<double> &b,
                 std::vector<double> &c,
                 std::vector<double> &d);

    void d2Field(double t, double a, 
		 std::vector<double> &b,
                 std::vector<double> &c,
                 std::vector<double> &d);

    void d3Field(double t, double a, 
		 std::vector<double> &b,
                 std::vector<double> &c,
                 std::vector<double> &d);
};


#endif //CPPPROJCT_LINEARTESTEQUATION_H
