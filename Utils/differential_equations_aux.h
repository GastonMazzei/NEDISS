//
// Created by m4zz31 on 2/11/21.
//

#ifndef CPPPROJCT_DIFFERENTIAL_EQUATIONS_AUX_H
#define CPPPROJCT_DIFFERENTIAL_EQUATIONS_AUX_H

#include <functional>
#include <vector>

#ifndef SCALAR_FLOW_TYPE
#define SCALAR_FLOW_TYPE
typedef std::function<double(double,double, std::vector<double>, std::vector<double>, std::vector<double>)> ScalarFlow;
#endif

double scalar_flow_placeholder(double,double, std::vector<double>, std::vector<double>, std::vector<double>);

#endif //CPPPROJCT_DIFFERENTIAL_EQUATIONS_AUX_H
