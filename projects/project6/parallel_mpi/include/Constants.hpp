#include <cstdlib>

#ifndef CASE
#define CASE 1
#endif

namespace PROBLEM {

    /**
     * @brief number of cells in x direction
     */
    const int Nx = 500;

    const double Lx = 10;

    const double DELTAX = Lx / Nx;

    const double ORIGINX = -5;
    const double CONTACTX = 0;


    const double NTIMESTEPS = 1000;
    const int OUTINTERVAL = 99999;

    /**
     * @brief timestep size
     */
    // const double TIMESTEP = TIMETARGET/NTIMESTEPS;

#if (CASE==1)
    /**
     * @brief Viscosity parameter Cx.
     */
    const double CX = 0.05;

    const double TIMETARGET = 0.005;

    /**
     * @brief timestep minimum parameter R
     */
    const double Rt = 0.6;
#endif


#if CASE==2

    const double TIMETARGET = 0.002;
    const double CX = 0.1;
    const double Rt = 0.6;

#endif

#if CASE==3
    const double TIMETARGET = 0.005;
    const double CX = 0.3;
    const double Rt = 0.6;
#endif

#if CASE==4
    const double TIMETARGET = 0.005;
    const double CX = 0.3;
    const double Rt = 0.6;
#endif

}

namespace CONST {

#if CASE==1

    const double PRESSURE1 = 1e4;
    const double PRESSURE4 = 1e5;

    const double PR = 3.031;

    const double RHO1 = 0.125;
    const double RHO4 = 1;

    const double U1 = 0;
    const double U4 = 0;

    const double MU = 1.79e-5;

    const double K = 2.63e-2;

    const double R = 287;

    const double GAMMA = 1.4;

    const double PRANTEL = 0.7;

#endif

#if CASE==2
    const double PRESSURE1 = 1e3;
    const double PRESSURE4 = 1e5;

    const double PR = 7.155;

    const double RHO1 = 0.01;
    const double RHO4 = 1;

    const double U1 = 0;
    const double U4 = 0;

    const double MU = 1.79e-5;

    const double K = 2.63e-2;

    const double R = 287;

    const double GAMMA = 1.4;

    const double PRANTEL = 0.7;
#endif

#if CASE==3
    const double PRESSURE1 = 5e4;
    const double PRESSURE4 = 1e5;

    const double PR = 1.452;

    const double RHO1 = 0.75;
    const double RHO4 = 1;

    const double U1 = 0;
    const double U4 = 0;

    const double MU = 1.79e-5;

    const double K = 2.63e-2;

    const double R = 287;

    const double GAMMA = 1.4;

    const double PRANTEL = 0.7;
#endif

#if CASE==4
    const double PRESSURE1 = 1e4;
    const double PRESSURE4 = 1e5;

    const double PR = 3.031;

    const double RHO1 = 0.125;
    const double RHO4 = 1;

    const double U1 = 0;
    const double U4 = 0;

    const double MU = 8.96e-6;

    const double K = 1.83e-1;

    const double R = 4124;

    const double GAMMA = 1.405;

    const double PRANTEL = 0.7;
#endif


const double TEMP1 = CONST::PRESSURE1 / (CONST::RHO1 * CONST::R);
const double TEMP4 = CONST::PRESSURE4 / (CONST::RHO4 * CONST::R);

const double Cv = CONST::R / (CONST::GAMMA - 1);

}