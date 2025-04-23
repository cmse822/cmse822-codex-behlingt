#include "../include/Constants.hpp"
#include "../include/MathFunctions.hpp"

#include <cstdlib>

#include <vector>
#include <array>

#include <iostream>

#include <cmath>

std::array<double, 3> obtainU(double rho, double u, double e){

    std::array<double, 3> output {0,0,0};
    output[0] = rho;
    output[1] = rho * u;
    output[2] = rho * e;

    return(output);
}

std::array<double, 4> obtainProperties(std::array<double, 3> U){

    std::array<double, 4> output {0,0,0,0};
    
    output[1] = U[0];
    output[2] = U[1] / U[0];
    output[3] = U[2] / U[0];

    output[0] = U[0] * CONST::R * obtainT(output[3]);

    return(output);
}

std::array<double, 3> obtainF(std::array<double, 3> &U, std::array<double, 2> &stress){
    //  prop = {P, rho, u, e}
    std::array<double, 4> prop {obtainProperties(U)};

    //Calculating F using the properties we got & the stresses provided;
    std::array<double, 3> output {0,0,0};

    output[0] = prop[1] * prop[2];
    output[1] = (prop[1] * prop[2] * prop[2]) + prop[0] - stress[0];
    output[2] = prop[1] * prop[2] * prop[3] + stress[1];

    return(output);
}

std::array<double, 3> obtainJ(std::array<double, 3> &U, 
                              std::array<double, 3> &Uprev, 
                              std::array<double, 3> &Unext, 
                              std::array<double, 2> &stress){
    //  prop = {P, rho, u, e}
    std::array<double, 4> prop {obtainProperties(U)};
    std::array<double, 4> prop_next {obtainProperties(Unext)};
    std::array<double, 4> prop_prev {obtainProperties(Uprev)};


    //Calculating F using the properties we got & the stresses provided;
    std::array<double, 3> output {0,0,0};

    output[2] = (stress[0] - prop[0]) * centralDifference(prop_prev[2], prop_next[2], PROBLEM::DELTAX);
    
    return(output);
}

std::array<double, 2> obtainStress(std::array<double, 3> &Uprev, 
                                   std::array<double, 3> &Unext) {
    //  prop = {P, rho, u, e}
    std::array<double, 4> prop_next {obtainProperties(Unext)};
    std::array<double, 4> prop_prev {obtainProperties(Uprev)};

    std::array<double, 2> output {0,0};
    output[0] = (4.0/3.0) * CONST::MU * centralDifference(prop_next[2], 
                                                          prop_prev[2], 
                                                          PROBLEM::DELTAX);

    output[1] = (-1.0) * CONST::K * centralDifference(obtainT(prop_next[3]), 
                                                      obtainT(prop_prev[3]), 
                                                      PROBLEM::DELTAX);

    return(output);
}

double obtainT(double e) {
    return(e / CONST::Cv);
}

double obtainE(double T) {
    return(T * CONST::Cv);
}

std::array<double, 3> forwardDifference(std::array<double, 3> &F, std::array<double, 3> &Fnext, double DelX) {
    std::array<double, 3> output {0,0,0};

    output[0] = (Fnext[0] - F[0]) / DelX;
    output[1] = (Fnext[1] - F[1]) / DelX;
    output[2] = (Fnext[2] - F[2]) / DelX;

    return(output);
}

std::array<double, 3> forwardDifference(std::array<double, 3> &F, 
                                        std::array<double, 3> &Fnext, 
                                        std::array<double, 3> &Fnext2,
                                        double DelX) {
    std::array<double, 3> output {0,0,0};

    output[0] = (-7*F[0] + 8*Fnext[0] - Fnext2[0]) / (6*DelX);
    output[1] = (-7*F[1] + 8*Fnext[1] - Fnext2[1]) / (6*DelX);
    output[2] = (-7*F[2] + 8*Fnext[2] - Fnext2[2]) / (6*DelX);
    return(output);
}

std::array<double, 3> backwardDifference(std::array<double, 3> &F, std::array<double, 3> &Fprev, double DelX) {
    std::array<double, 3> output {0,0,0};

    output[0] = (F[0] - Fprev[0]) / DelX;
    output[1] = (F[1] - Fprev[1]) / DelX;
    output[2] = (F[2] - Fprev[2]) / DelX;

    return(output);
}


std::array<double, 3> backwardDifference(std::array<double, 3> &F, 
                                         std::array<double, 3> &Fprev, 
                                         std::array<double, 3> &Fprev2,
                                         double DelX) {
    std::array<double, 3> output {0,0,0};

    output[0] = (7*F[0] - 8*Fprev[0] + Fprev2[0]) / (6*DelX);
    output[1] = (7*F[1] - 8*Fprev[1] + Fprev2[1]) / (6*DelX);
    output[2] = (7*F[2] - 8*Fprev[2] + Fprev2[2]) / (6*DelX);

    return(output);
}


double centralDifference(double Vprev, double Vnext, double DelX) {
    return ((Vnext - Vprev) / (2*DelX));
}

std::array<double, 3> predictorStep(std::array<double, 3> U, 
                                    std::array<double, 3> Fdiff, 
                                    std::array<double, 3> &J,
                                    double timestep){

    std::array<double, 3> output {0,0,0};
    output[0] = U[0] - (timestep * (Fdiff[0] - J[0]));
    output[1] = U[1] - (timestep * (Fdiff[1] - J[1]));
    output[2] = U[2] - (timestep * (Fdiff[2] - J[2]));

    return(output);
}

std::array<double, 3> correctorStep(std::array<double, 3> U, 
                                    std::array<double, 3> &Ubar, 
                                    std::array<double, 3> Fbardiff, 
                                    std::array<double, 3> &Jbar,
                                    double timestep){

    std::array<double, 3> output {0,0,0};

    output[0] = 0.5 * (U[0] + Ubar[0] - (timestep * (Fbardiff[0] - Jbar[0])));
    output[1] = 0.5 * (U[1] + Ubar[1] - (timestep * (Fbardiff[1] - Jbar[1])));
    output[2] = 0.5 * (U[2] + Ubar[2] - (timestep * (Fbardiff[2] - Jbar[2])));

    return(output);
}

std::array<double, 3> artificialViscosity(std::array<double, 3> &U, 
                                          std::array<double, 3> &Uprev, 
                                          std::array<double, 3> &Unext){

    std::array<double, 3> output {0,0,0};

    //First, we want to find the front term;
    std::array<double, 4> U_prop = obtainProperties(U);
    std::array<double, 4> Uprev_prop = obtainProperties(Uprev);
    std::array<double, 4> Unext_prop = obtainProperties(Unext);

    double first_term {0};
    first_term = PROBLEM::CX * std::abs(Unext_prop[0] - 2*U_prop[0] + Uprev_prop[0]) / (Unext_prop[0] + U_prop[0] + Uprev_prop[0]);
    
    output[0] = first_term * (Unext[0] - 2*U[0] + Uprev[0]);
    output[1] = first_term * (Unext[1] - 2*U[1] + Uprev[1]);
    output[2] = first_term * (Unext[2] - 2*U[2] + Uprev[2]);

    return(output);
}

double geta(double T) {return pow( CONST::GAMMA * CONST::R * T, 0.5);}

double obtainTimestep(std::vector<std::array<double, 3>> &fieldU){
    //we have the entire field.
    //first, for each U, we want to calculate that timestep;
    // we then compare to the previous; if new is lower, we insert the new one
    double min_delt {1e10};

    for(int ix {0}; ix < PROBLEM::Nx; ix++){

        //first, get properties for that U;
        std::array<double, 4> prop = obtainProperties(fieldU[ix]);

        //now, determine max nu;
        double nu1 = (4/3) * CONST::MU / prop[1];
        double nu2 = (CONST::MU / prop[1]) * (CONST::GAMMA / CONST::PRANTEL);

        double nu {};
        if(nu1 > nu2) {nu = nu1;}
        else {nu = nu2;}

        double x = PROBLEM::DELTAX; //PROBLEM::ORIGINX + ix * PROBLEM::DELTAX;

        double term1 = std::abs(prop[2]) / x;
        double term2 = geta(obtainT(prop[3])) / x;
        double term3 = (2*nu) / (PROBLEM::DELTAX * PROBLEM::DELTAX);

        double delt = 1 / (term1 + term2 + term3);

        if(delt < min_delt) {min_delt = delt;}
    }

    return(min_delt * PROBLEM::Rt);

}