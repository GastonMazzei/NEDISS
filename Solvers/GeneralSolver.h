//
// Created by m4zz31 on 1/11/21.
//

#ifndef CPPPROJCT_GENERALSOLVER_H
#define CPPPROJCT_GENERALSOLVER_H

#include <functional>


#include "../GraphClasses/GeneralGraph.h"
#include "../DifferentialEquations/GeneralDifferentialEquation.h"

typedef std::function<double(double, double, std::vector<double>, std::vector<double>, std::vector<double>)> ScalarFlow;
typedef std::function<double(double a, std::vector<double> &b, std::vector<double> &c, std::vector<double> &d, Graph &g, ScalarFlow &F)> SolverOp;

struct TimeStructure{
    double h = 0.01;
    double t = 0;
};


template <typename DIFFEQ, typename SOLVER>
class GeneralSolver{
public:
    // Configuration specific to the solver structure :-)
    double Params[4]; // up to 4 weights for runge-kutta
    //                      for example, for Heun's method it would be {1/2,1/2,0,0}
    //              Note: in this example, k3 and k4 would still be computed if deg=4.
    int deg; // the degree of the integration algorithm

    std::string type;
    GeneralSolver();
    GeneralSolver(std::string valtype);
    GeneralSolver(std::string valtype, int d);
    GeneralSolver(std::string valtype, int d, double params[4]);

    // Instantiate the Differential Equation and prepare stuff if necessary
    DIFFEQ DifferentialEquation;
    SOLVER Solver;
    void PostTemplateProcessing();


    // main function: evolver
    void evolve(double a,
                std::vector<double> &b,
                std::vector<double> &c,
                std::vector<double> &d,
                double &answer);

    // Configuration specific to the time structure :-)
    TimeStructure T;
    void SetStep(double h);
    void SetT0(double t0);
    void EvolveTime();
};

template <typename DIFFEQ, typename SOLVER>
GeneralSolver<DIFFEQ, SOLVER>::GeneralSolver(std::string valtype, int d, double params[4]){
    type = valtype;
    deg = d;
    Params[0] = params[0];
    Params[1] = params[1];
    Params[2] = params[2];
    Params[3] = params[3];
};

template <typename DIFFEQ, typename SOLVER>
GeneralSolver<DIFFEQ, SOLVER>::GeneralSolver(std::string valtype){
    // Types available should be:
    //  (1) Euler ('eu')
    //  (2) Runge-Kutta ('rk')
    std::string methods_str = "(1) Euler ('eu'), (2) Runge-Kutta up to order 4 included ('rk').";
    type = std::move(valtype);
    if (type == "rk") {
        // Heun method's initialization :-)
        deg = 2;
        Params[0] = 0.5;
        Params[1] = 0.5;
        Params[2] = 0;
        Params[3] = 0;
    } else if (type == "eu") {
        deg = 1;
    } else if (type != "eu") {
        error_report("Only support for Integration is available for the following methods:" + methods_str);
    }
};

template <typename DIFFEQ, typename SOLVER>
GeneralSolver<DIFFEQ, SOLVER>::GeneralSolver(){
    type = "eu";
    deg = 1;
}

template <typename DIFFEQ, typename SOLVER>
GeneralSolver<DIFFEQ, SOLVER>::GeneralSolver(std::string valtype, int d) {
    // Types available should be:
    //  (1) Euler ('eu')
    //  (2) Runge-Kutta ('rk')
    std::string methods_str = "(1) Euler ('eu'), (2) Runge-Kutta up to order 4 included ('rk').";
    type = std::move(valtype);
    deg = d;
    if (type == "rk") {
        // Heun method's initialization :-)
        Params[0] = 0.5;
        Params[1] = 0.5;
        Params[2] = 0;
        Params[3] = 0;
    } else if (type != "eu") {
        error_report("Only support for Integration is available for the following methods:" + methods_str);
    }
};

template <typename DIFFEQ, typename SOLVER>
void GeneralSolver<DIFFEQ, SOLVER>::SetStep(double h){
    T.h = h;
}

template <typename DIFFEQ, typename SOLVER>
void GeneralSolver<DIFFEQ, SOLVER>::SetT0(double t0){
    T.t = t0;
}

template <typename DIFFEQ, typename SOLVER>
void GeneralSolver<DIFFEQ, SOLVER>::EvolveTime(){
    T.t += T.h;
}



template <typename DIFFEQ, typename SOLVER>
void GeneralSolver<DIFFEQ, SOLVER>::evolve(double a,
                                             std::vector<double> &b,
                                             std::vector<double> &c,
                                             std::vector<double> &d,
                                             double &answer){
    return Solver.evolve(T.t, T.h, a, b, c, d, DifferentialEquation,
                         DifferentialEquation.Specs, answer);
}

template <typename DIFFEQ, typename SOLVER>
void GeneralSolver<DIFFEQ, SOLVER>::PostTemplateProcessing(){
    if (DifferentialEquation.RequiresBuilding){
        DifferentialEquation.BuildForSolver();
    };
}




#endif //CPPPROJCT_GENERALSOLVER_H