//including header file;
#include "../include/Constants.hpp"
#include "../include/MathFunctions.hpp"

#include <cstdlib>

#include <vector>
#include <array>

#include <iostream>
#include <fstream>
#include <string>

#include <chrono>

/**
 * @brief writes the vector of arrays to a file.
 * Each line represents one value;
 * P, u, rho, T
 */
void writeArrFile(std::vector<std::array<double, 3>> vec_arr, std::string filename) {
    std::ofstream myfile (filename);

    if (myfile.is_open()) {

        for(int i = 0; i < 3; i++) {
            for(int count = 0; count < vec_arr.size(); count ++){
                myfile << vec_arr[count][i];
                if(count < vec_arr.size() - 1){myfile << ",";}
            }
             myfile << '\n';
        }
      myfile.close();
    }
    else {
        std::cout << "failed to open file. " << '\n';
    }
}



int main () {

    //Instantiating arrays for U, F, J

    /**
     * @brief Vector of arrays representing U
     */
    std::vector<std::array<double, 3>> fieldU; fieldU.resize(PROBLEM::Nx, {0,0,0});
    std::vector<std::array<double, 3>> fieldUbar; fieldUbar.resize(PROBLEM::Nx, {0,0,0});

    /**
     * @brief Vector of arrays representing F
     */
    std::vector<std::array<double, 3>> fieldF; fieldF.resize(PROBLEM::Nx, {0,0,0});
    std::vector<std::array<double, 3>> fieldFbar; fieldFbar.resize(PROBLEM::Nx, {0,0,0});


    /**
     * @brief Vector of arrays representing J
     */
    std::vector<std::array<double, 3>> fieldJ; fieldJ.resize(PROBLEM::Nx, {0,0,0});
    std::vector<std::array<double, 3>> fieldJbar; fieldJbar.resize(PROBLEM::Nx, {0,0,0});

    /**
     * @brief Vector of arrays representing Stress; tau: 0, q: 1
     */
    std::vector<std::array<double, 2>> stress; stress.resize(PROBLEM::Nx, {0,0});

    // =====================
    // Initializing values
    // =====================

    int stepcount {0};

    //Iterate through locations in x;
    for(int ix {0}; ix < PROBLEM::Nx; ix++){
        //the current position in x
        double cpos {PROBLEM::ORIGINX + ix * PROBLEM::DELTAX};
        
        //initialize temporary variables
        double rho {0};
        double P {0};
        double T {0};
        double u {0};

        //If we are above the contact;
        if(cpos > PROBLEM::CONTACTX){
            rho = CONST::RHO1;
            P = CONST::PRESSURE1;
            T = CONST::TEMP1;
            u = CONST::U1;
        }

        //If we are below the contact;
        else {
            rho = CONST::RHO4;
            P = CONST::PRESSURE4;
            T = CONST::TEMP4;
            u = CONST::U4;
        }

        //getting the neighboring indices, accounting for BCs
        int px {ix-1};
        int nx {ix+1};
        if(ix == 0) {px = ix;}
        if(ix == PROBLEM::Nx - 1) {nx = ix;}

        //initialize this U, F, J
        fieldU[ix] = obtainU(rho, u, obtainE(T));
        fieldF[ix] = obtainF(fieldU[ix], stress[ix]);
        fieldJ[ix] = obtainJ(fieldU[ix], fieldU[px], fieldU[nx], stress[ix]);
    }

    writeArrFile(fieldU, "output/initial.txt");

    //Now, we want to while loop up to our final time;
    double ctime {0};

    auto start = std::chrono::high_resolution_clock::now();

    std::cout << "=========== Starting ===========" << '\n';

    while(ctime < PROBLEM::TIMETARGET){

        //incrementing time upwards.
        // ctime += PROBLEM::TIMESTEP;
        //using the non-equal timesteps instead;
        double timestep = obtainTimestep(fieldU);

        //if this would put us past the final time, clamp;
        if((ctime + timestep) >= PROBLEM::TIMETARGET){
            timestep = PROBLEM::TIMETARGET - ctime;
            ctime = PROBLEM::TIMETARGET;
        }
        //otherwise, increase ctime.
        else {ctime += timestep;}
        

        //update stress vector to this timestep;
        for(int ix {0}; ix < PROBLEM::Nx; ix++){
            if(ix == 0)                  {stress[ix] = obtainStress(fieldU[ix], fieldU[ix+1]);}
            else if(ix == PROBLEM::Nx-1) {stress[ix] = obtainStress(fieldU[ix-1], fieldU[ix]);}
            else                         {stress[ix] = obtainStress(fieldU[ix-1], fieldU[ix+1]);}
        }

        //updating F and J at this timestep;
        for(int ix {0}; ix < PROBLEM::Nx; ix++){

            int px {ix-1};
            int nx {ix+1};
            if(ix == 0) {px = ix;}
            if(ix == PROBLEM::Nx - 1) {nx = ix;}

            fieldF[ix] = obtainF(fieldU[ix], stress[ix]);
            fieldJ[ix] = obtainJ(fieldU[ix], fieldU[px], fieldU[nx], stress[ix]);
        }

        std::array<double, 3> derivF {0,0,0}; 

        int in {0}; //index of the next U
        int ip {0}; //index of the previous U

        //Performing the Predictor step;
        for(int ix{0}; ix < PROBLEM::Nx; ix++){
            //Obtaining the derivative of F

            //Determine BC behavior;
            if(ix == 0) {derivF = {0,0,0};}
            else {derivF = backwardDifference(fieldF[ix], fieldF[ix-1], PROBLEM::DELTAX);}

            //Performing predictor step;
            fieldUbar[ix] = predictorStep(fieldU[ix], derivF, fieldJ[ix], timestep);

            //determine BC behavior for U;
            if(ix == 0) {ip = ix; in = ix+1;}
            else if(ix == PROBLEM::Nx-1) {ip = ix-1; in = ix;}
            else {ip = ix-1; in = ix+1;}

            //adding artificial viscosity term;
            std::array<double, 3> visc = artificialViscosity(fieldU[ix], fieldU[ip], fieldU[in]);
            fieldUbar[ix][0] += visc[0];
            fieldUbar[ix][1] += visc[1];
            fieldUbar[ix][2] += visc[2];

        }

        //get predictor stress;
        for(int ix {0}; ix < PROBLEM::Nx; ix++){
            if(ix == 0)                  {stress[ix] = obtainStress(fieldUbar[ix], fieldUbar[ix+1]);}
            else if(ix == PROBLEM::Nx-1) {stress[ix] = obtainStress(fieldUbar[ix-1], fieldUbar[ix]);}
            else                         {stress[ix] = obtainStress(fieldUbar[ix-1], fieldUbar[ix+1]);}
        }

        //get predictor F and J;
        for(int ix{0}; ix < PROBLEM::Nx; ix++){

            int px {ix-1};
            int nx {ix+1};
            if(ix == 0) {px = ix;}
            if(ix == PROBLEM::Nx - 1) {nx = ix;}

            fieldFbar[ix] = obtainF(fieldUbar[ix], stress[ix]);
            fieldJbar[ix] = obtainJ(fieldUbar[ix], fieldUbar[px], fieldUbar[nx], stress[ix]);
        }

        std::array<double, 3> derivFbar {0,0,0};
        //Now, we do the corrector step;
        for(int ix{0}; ix < PROBLEM::Nx; ix++){

            if(ix == PROBLEM::Nx-1) {derivFbar = {0,0,0};}
            else {derivFbar = forwardDifference(fieldFbar[ix], fieldFbar[ix+1], PROBLEM::DELTAX);}

            fieldU[ix] = correctorStep(fieldU[ix], fieldUbar[ix], derivFbar, fieldJbar[ix], timestep);

            //determine BC behavior for U;
            if(ix == 0) {ip = ix; in = ix+1;}
            else if(ix == PROBLEM::Nx-1) {ip = ix-1; in = ix;}
            else {ip = ix-1; in = ix+1;}

            //adding artificial viscosity term;
            std::array<double, 3> visc = artificialViscosity(fieldUbar[ix], fieldUbar[ip], fieldUbar[in]);
            fieldU[ix][0] += visc[0];
            fieldU[ix][1] += visc[1];
            fieldU[ix][2] += visc[2];

        }

        //We are donesies!!!
        // We have applied the predictor and corrector steps, so fieldU is up-to-date.

        //putting out a message at a certain cadence;
        // if(stepcount % PROBLEM::OUTINTERVAL == 0) {
        //     std::cout << "time: " << ctime << "/" << PROBLEM::TIMETARGET << " - " << 100 * ctime / PROBLEM::TIMETARGET << "%" << " - it: " << stepcount << '\n'; 
        // }

        stepcount++;

    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    //printing the final timestep.
    //std::cout << "time: " << ctime << "/" << PROBLEM::TIMETARGET << " - " << 100 * ctime / PROBLEM::TIMETARGET << "%" << " - it: " << stepcount << '\n'; 
    std::cout << "time elapsed: " << elapsed.count() << " seconds" << std::endl;
    std::cout << "it/sec: " << stepcount / elapsed.count() << std::endl;
    //writing out the final result to a text file.
    // writeArrFile(fieldU, "output/final.txt");
    std::cout << "=========== Done! ===========" << '\n';

    return 0;
}

