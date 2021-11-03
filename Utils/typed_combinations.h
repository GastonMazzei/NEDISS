//
// Created by m4zz31 on 3/11/21.
//

#ifndef CPPPROJCT_TYPED_COMBINATIONS_H
#define CPPPROJCT_TYPED_COMBINATIONS_H

#include "../Solvers/GeneralSolver.h"
#include "../Solvers/EulerSolver.h"
#include "../Solvers/RungeKuttaSolver.h"
#include "../DifferentialEquations/NoiselessKuramoto.h"

template struct GeneralSolver<NoiselessKuramoto, EulerSolver<NoiselessKuramoto>>;
template struct GeneralSolver<NoiselessKuramoto, RungeKuttaSolver<NoiselessKuramoto>>;


#endif //CPPPROJCT_TYPED_COMBINATIONS_H
