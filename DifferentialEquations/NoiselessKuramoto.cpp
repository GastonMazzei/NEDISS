//
// Created by m4zz31 on 2/11/21.
//

#include "NoiselessKuramoto.h"
#include <cmath>
#include <cassert>
#include <vector>

auto GeneralKuramotoFactory(std::vector<double> &T1,
                            std::vector<double> &T2,
                            std::vector<double> &T3,
                            std::vector<double> &T4,
                            int &N) {
    //
    // General template for applying a perturbation to a field :-)
    //
    // This is the case for the Noiseless Kuramoto equation :-)
    //
    // General flows cant be accepted as we would need to pass by copy
    // instead of reference, which would only be acceptable if we change
    // vectors to arrays, OR actively manage memory here, which is not
    // necessarily well behaved and its even potentially leaky :-(
    //
    // Instead of calling F(t,y)
    // we can call
    // F(T1 + T2 * t, T3 + T4 * y)
    //
    // FOR KURAMOTO: we actually dont care about the T1 and T2 terms ;-)
    //               our system is a closed system ^.^
    //
    int j=0,k=0;
    if (T3.size()!=1){
        assert(T3.size() == N);
        int j=1;
    }
    if (T4.size()!=1){
        assert(T4.size() == N);
        int k=1;
    }
    return [&, j, k] (double t, double a, std::vector<double> &b,
                          std::vector<double> &c,
                          std::vector<double> &d) -> double
    {
        // a: central value
        // b: own parameters
        // c: neighbor values
        // d: interactions
        // T3 & T4: linear offset and factor (respect.) to apply to "c"
        // (we are asking they have the same length... is it too much to ask?)
        double result;
        result = b[0];
        for (int i = 0; i < b.size(); i++) {
            result += d[i] * std::sin((T3[j*i] + T4[k*i] * c[i]) - a);
        }
        return result;
        };
};


NoiselessKuramoto::NoiselessKuramoto(): GeneralDifferentialEquation(0){
    //Field = KuramotoKernel;
};

