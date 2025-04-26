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

#include <omp.h>

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
        double rho {0}; rho = 0;
        double P {0}; P = 0;
        double T {0}; T = 0;
        double u {0}; u = 0;

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

    //As part of this, I need to swap to pointers (unfortunately).

    std::array<double, 3>* fieldU_data {fieldU.data()}; 
    std::array<double, 3>* fieldUbar_data {fieldUbar.data()};

    std::array<double, 3>* fieldF_data {fieldF.data()};
    std::array<double, 3>* fieldFbar_data {fieldFbar.data()};

    std::array<double, 3>* fieldJ_data {fieldJ.data()};
    std::array<double, 3>* fieldJbar_data {fieldJbar.data()};

    std::array<double, 2>* stress_data {stress.data()};

    std::cout << "=========== Starting ===========" << '\n';

    int Nx = PROBLEM::Nx;


    #pragma omp target data map(tofrom: fieldU_data[0:Nx], fieldUbar_data[0:Nx], \
        fieldF_data[0:Nx], fieldFbar_data[0:Nx], \
        fieldJ_data[0:Nx], fieldJbar_data[0:Nx], \
        stress_data[0:Nx])
    {

    while(ctime < PROBLEM::TIMETARGET){

        //incrementing time upwards.
        double timestep = PROBLEM::TIMESTEP;
        //using the non-equal timesteps instead;
        // double timestep = obtainTimestep(fieldU_data);
        ctime += timestep;

        //if this would put us past the final time, clamp;
        if((ctime) >= PROBLEM::TIMETARGET){
            timestep = PROBLEM::TIMETARGET - ctime;
            ctime = PROBLEM::TIMETARGET;
        }

        //update stress vector to this timestep;
        // handle the two edges on their own, outside the loop;
        stress_data[0]             = obtainStress(fieldU_data[0]            , fieldU_data[1]);
        stress_data[PROBLEM::Nx-1] = obtainStress(fieldU_data[PROBLEM::Nx-2], fieldU_data[PROBLEM::Nx-1]);

        // #pragma omp parallel for
        // #pragma omp target teams distribute parallel for
        #pragma omp target loop
        for(int ix {1}; ix < PROBLEM::Nx-1; ix++){
            stress_data[ix] = obtainStress(fieldU_data[ix-1], fieldU_data[ix+1]);
        }

        //updating F and J at this timestep;
        // #pragma omp parallel for
        // #pragma omp target teams distribute parallel for
        #pragma omp target loop
        for(int ix {0}; ix < PROBLEM::Nx; ix++){

            int px {ix-1};
            int nx {ix+1};
            if(ix == 0) {px = ix;}
            if(ix == PROBLEM::Nx - 1) {nx = ix;}

            fieldF_data[ix] = obtainF(fieldU_data[ix], stress_data[ix]);
            fieldJ_data[ix] = obtainJ(fieldU_data[ix], fieldU_data[px], fieldU_data[nx], stress_data[ix]);
        }


        std::array<double, 3> derivF {0,0,0}; 

        int in {0}; //index of the next U
        int ip {0}; //index of the previous U

        //Performing the Predictor step;
        // #pragma omp parallel for
        // #pragma omp target teams distribute parallel for
        #pragma omp target loop
        for(int ix{0}; ix < PROBLEM::Nx; ix++){
            //Obtaining the derivative of F

            // std::array<double, 3> derivF;
            //Determine BC behavior;
            if(ix == 0) {derivF = {0,0,0};}
            else {derivF = backwardDifference(fieldF_data[ix], fieldF_data[ix-1], PROBLEM::DELTAX);}

            //Performing predictor step;
            fieldUbar_data[ix] = predictorStep(fieldU_data[ix], derivF, fieldJ_data[ix], timestep);

            //determine BC behavior for U;
            if(ix == 0) {ip = ix; in = ix+1;}
            else if(ix == PROBLEM::Nx-1) {ip = ix-1; in = ix;}
            else {ip = ix-1; in = ix+1;}

            //adding artificial viscosity term;
            std::array<double, 3> visc = artificialViscosity(fieldU_data[ix], fieldU_data[ip], fieldU_data[in]);
            fieldUbar_data[ix][0] += visc[0];
            fieldUbar_data[ix][1] += visc[1];
            fieldUbar_data[ix][2] += visc[2];

        }

        //get predictor stress;
        // #pragma omp parallel for
        // #pragma omp target teams distribute parallel for
        #pragma omp target loop
        for(int ix {0}; ix < PROBLEM::Nx; ix++){
            if(ix == 0)                  {stress_data[ix] = obtainStress(fieldUbar_data[ix], fieldUbar_data[ix+1]);}
            else if(ix == PROBLEM::Nx-1) {stress_data[ix] = obtainStress(fieldUbar_data[ix-1], fieldUbar_data[ix]);}
            else                         {stress_data[ix] = obtainStress(fieldUbar_data[ix-1], fieldUbar_data[ix+1]);}
        }

        //get predictor F and J;
        // #pragma omp parallel for
        // #pragma omp target teams distribute parallel for
        #pragma omp target loop
        for(int ix{0}; ix < PROBLEM::Nx; ix++){

            int px {ix-1};
            int nx {ix+1};
            if(ix == 0) {px = ix;}
            if(ix == PROBLEM::Nx - 1) {nx = ix;}

            fieldFbar_data[ix] = obtainF(fieldUbar_data[ix], stress_data[ix]);
            fieldJbar_data[ix] = obtainJ(fieldUbar_data[ix], fieldUbar_data[px], fieldUbar_data[nx], stress_data[ix]);
        }

        std::array<double, 3> derivFbar {0,0,0};
        //Now, we do the corrector step;
        // #pragma omp parallel for
        // #pragma omp target teams distribute parallel for
        #pragma omp target loop
        for(int ix{0}; ix < PROBLEM::Nx; ix++){

            if(ix == PROBLEM::Nx-1) {derivFbar = {0,0,0};}
            else {derivFbar = forwardDifference(fieldFbar_data[ix], fieldFbar_data[ix+1], PROBLEM::DELTAX);}

            fieldU_data[ix] = correctorStep(fieldU_data[ix], fieldUbar_data[ix], derivFbar, fieldJbar_data[ix], timestep);

            //determine BC behavior for U;
            if(ix == 0) {ip = ix; in = ix+1;}
            else if(ix == PROBLEM::Nx-1) {ip = ix-1; in = ix;}
            else {ip = ix-1; in = ix+1;}

            //adding artificial viscosity term;
            std::array<double, 3> visc = artificialViscosity(fieldUbar_data[ix], fieldUbar_data[ip], fieldUbar_data[in]);
            fieldU_data[ix][0] += visc[0];
            fieldU_data[ix][1] += visc[1];
            fieldU_data[ix][2] += visc[2];

        }

        //We are donesies!!!
        // We have applied the predictor and corrector steps, so fieldU is up-to-date.

        //putting out a message at a certain cadence;
        // if(stepcount % PROBLEM::OUTINTERVAL == 0) {
        //     std::cout << "time: " << ctime << "/" << PROBLEM::TIMETARGET << " - " << 100 * ctime / PROBLEM::TIMETARGET << "%" << " - it: " << stepcount << '\n'; 
        // }

        stepcount++;

    }

}

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    //printing the final timestep.
    //std::cout << "time: " << ctime << "/" << PROBLEM::TIMETARGET << " - " << 100 * ctime / PROBLEM::TIMETARGET << "%" << " - it: " << stepcount << '\n'; 
    std::cout << "time elapsed: " << elapsed.count() << " seconds" << std::endl;
    std::cout << "it/sec: " << stepcount / elapsed.count() << std::endl;
    //writing out the final result to a text file.
    writeArrFile(fieldU, "output/final_cpu.txt");
    std::cout << "=========== Done! ===========" << '\n';

    return 0;
}

