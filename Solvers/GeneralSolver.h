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

struct SolverConfig{
    int s=0;
    int d=0;
    double P[4];
    SolverConfig();
};

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
    bool requires_communication = false;
    GeneralSolver();
    GeneralSolver(std::string valtype);
    GeneralSolver(std::string valtype, int d);
    GeneralSolver(std::string valtype, int d, double * params);

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
GeneralSolver<DIFFEQ, SOLVER>::GeneralSolver(std::string valtype, int d, double * params){
    type = valtype;
    deg = d;
    Params[0] = *(params+0);
    Params[1] = *(params+1);
    Params[2] = *(params+2);
    Params[3] = *(params+3);
    printf("[WARNING] GeneralSoler constructor that  explicitly uses four weights  params assigns 'true' to 'requires_communication', change it if the current method does not require it. Current method: %s with d=%d\n", valtype, d);
    std::cout << std::flush;
    requires_communication = true;
};

template <typename DIFFEQ, typename SOLVER>
GeneralSolver<DIFFEQ, SOLVER>::GeneralSolver(std::string valtype){
    // Types available should be:
    //  (1) Euler ('eu')
    //  (2) Runge-Kutta ('rk')
    std::string methods_str = "(1) Euler ('eu').";
    type = std::move(valtype);
    if (type == "rk") {
        // Heun method's initialization :-)
        requires_communication =  true;
        deg = 2;
        Params[0] = 0.5;
        Params[1] = 0.5;
        Params[2] = 0;
        Params[3] = 0;
    } else if (type == "eu") {
        deg = 1;
        Params[0] = 1;
        Params[1] = 0;
        Params[2] = 0;
        Params[3] = 0;
    } else if (type != "eu") {
        error_report("Only support for Integration is available for the following methods:" + methods_str);
    }
};

template <typename DIFFEQ, typename SOLVER>
GeneralSolver<DIFFEQ, SOLVER>::GeneralSolver(){
    type = "eu";
    deg = 1;
    Params[0] = 1;
    Params[1] = 0;
    Params[2] = 0;
    Params[3] = 0;
    requires_communication = false;
}

template <typename DIFFEQ, typename SOLVER>
GeneralSolver<DIFFEQ, SOLVER>::GeneralSolver(std::string valtype, int d) {
    // Types available should be:
    //  (1) Euler ('eu')
    //  (2) Runge-Kutta ('rk')
    std::string methods_str = "(1) Euler ('eu').";
    type = std::move(valtype);
    deg = d;
    if (type == "rk") {
        if (d == 2) {
            // Heun method's initialization :-)
            Params[0] = 0.5;
            Params[1] = 0.5;
            Params[2] = 0;
            Params[3] = 0;
            requires_communication = true;
        } else {
            printf("[FATAL] GeneralSolver(type='rk', int d) only accepts d=2 for Heun's method, but %d was the input.\n",d);
            std::cout<<std::flush;
            exit(1);
        }
    } else if (type == "eu") {
        // Euler order 1 :-)
        for (int i=0; i<4; ++i){
            if (i+1 <= d){
                Params[i] = 1;
            } else {
                Params[i] = 0;
            }
        }
        if (d!=1) PRINTF_DBG("\n\n[WARNING]\n\nEuler od Order != 1 may not be available for some Equations, i.e. some equation classes may not have defined their field derivatives up to the requested order :O.\n\n");
    } else {
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
                         DifferentialEquation.Specs, answer, &Params[0]);
}

template <typename DIFFEQ, typename SOLVER>
void GeneralSolver<DIFFEQ, SOLVER>::PostTemplateProcessing(){
    if (DifferentialEquation.RequiresBuilding){
        DifferentialEquation.BuildForSolver();
    };
}




#endif //CPPPROJCT_GENERALSOLVER_H