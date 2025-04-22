#include <cstdlib>

#include <vector>
#include <array>

#include <iostream>

/**
 * @brief Obtaining U vector from initial density, velocity u, and internal energy e.
 */
std::array<double, 3> obtainU(double rho, double u, double e);

/**
 * @brief returns P, rho, u, e
 */
std::array<double, 4> obtainProperties(std::array<double, 3> U);

/**
 * @brief Obtaining F vector from U vector and Stress vector.
 */
std::array<double, 3> obtainF(std::array<double, 3> &U, std::array<double, 2> &stress);

/**
 * @brief Obtaining J vector from U vector and Stress vector.
 */
std::array<double, 3> obtainJ(std::array<double, 3> &U,
                              std::array<double, 3> &Uprev, 
                              std::array<double, 3> &Unext,  
                              std::array<double, 2> &stress);


std::array<double, 2> obtainStress(std::array<double, 3> &Uprev, 
                                   std::array<double, 3> &Unext);

const double obtainT(double e);
const double obtainE(double T);

/**
 * @brief Forward Differencing Function;
 */
std::array<double, 3> forwardDifference(std::array<double, 3> &F, std::array<double, 3> &Fnext, double DelX);
std::array<double, 3> forwardDifference(std::array<double, 3> &F, std::array<double, 3> &Fnext, std::array<double, 3> &Fnext2, double DelX);

/**
 * @brief Backward Differencing Function;
 */
std::array<double, 3> backwardDifference(std::array<double, 3> &F, std::array<double, 3> &Fprev, double DelX);
std::array<double, 3> backwardDifference(std::array<double, 3> &F, std::array<double, 3> &Fprev, std::array<double, 3> &Fprev2, double DelX);


/**
 * @brief Central Differencing Function (single double);
 */
double centralDifference(double Vprev, double Vnext, double DelX);


std::array<double, 3> predictorStep(std::array<double, 3> U, std::array<double, 3> Fdiff, std::array<double, 3> &J, double timestep);

std::array<double, 3> correctorStep(std::array<double, 3> U, std::array<double, 3> &Ubar, std::array<double, 3> Fbardiff, std::array<double, 3> &Jbar, double timestep);

std::array<double, 3> artificialViscosity(std::array<double, 3> &U, std::array<double, 3> &Uprev, std::array<double, 3> &Unext);

double obtainTimestep(std::vector<std::array<double, 3>> &fieldU);