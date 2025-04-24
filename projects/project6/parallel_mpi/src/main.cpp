//including header file;
#include "../include/Constants.hpp"
#include "../include/MathFunctions.hpp"
#include "../include/MPIFunctions.hpp"


#include <cstdlib>

#include <vector>
#include <array>

#include <iostream>
#include <fstream>
#include <string>

#include <chrono>

#include <mpi.h>

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



//Instantiating arrays for U, F, J

/**
 * @brief Vector of arrays representing U
 */
std::vector<std::array<double, 3>> fieldU;   
std::vector<std::array<double, 3>> fieldUbar;

/**
 * @brief Vector of arrays representing F
 */
std::vector<std::array<double, 3>> fieldF;    
std::vector<std::array<double, 3>> fieldFbar; 


/**
 * @brief Vector of arrays representing J
 */
std::vector<std::array<double, 3>> fieldJ;   
std::vector<std::array<double, 3>> fieldJbar;

/**
 * @brief Vector of arrays representing Stress; tau: 0, q: 1
 */
std::vector<std::array<double, 2>> stress;


// setting up ghost value arrays;
std::array<std::array<double, 3>, 2> fieldU_ghost;    
std::array<std::array<double, 3>, 2> fieldUbar_ghost; 

std::array<std::array<double, 3>, 2> fieldF_ghost;    
std::array<std::array<double, 3>, 2> fieldFbar_ghost; 

std::array<std::array<double, 3>, 2> fieldJ_ghost;    
std::array<std::array<double, 3>, 2> fieldJbar_ghost; 

std::array<std::array<double, 2>, 2> stress_ghost;    

std::array<int, 2> neighbors;

void exchangeFields(int comm_size, int comm_rank){
    exchangeFieldGhostCells(fieldU, fieldU_ghost, neighbors, comm_size, comm_rank);
    exchangeFieldGhostCells(fieldUbar, fieldUbar_ghost, neighbors, comm_size, comm_rank);

    exchangeFieldGhostCells(fieldF, fieldF_ghost, neighbors, comm_size, comm_rank);
    exchangeFieldGhostCells(fieldFbar, fieldFbar_ghost, neighbors, comm_size, comm_rank);

    exchangeFieldGhostCells(fieldJ, fieldJ_ghost, neighbors, comm_size, comm_rank);
    exchangeFieldGhostCells(fieldJbar, fieldJbar_ghost, neighbors, comm_size, comm_rank);

    exchangeStressGhostCells(stress, stress_ghost, neighbors, comm_size, comm_rank);
}


int main (int argc, char* argv[]) {

    //starts MPI
    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD; //communication object; how we TALK

    //declaring variables
    int comm_rank{};
    int comm_size{};

    //this guy gives us rank (what num) and overall n of processors
    //  wants mem address of ints.
    MPI_Comm_rank(comm, &comm_rank);
    MPI_Comm_size(comm, &comm_size);

    //finding neighbors;
    neighbors[0] = comm_rank-1;
    neighbors[1] = comm_rank+1;

    if(comm_rank == 0) {
        neighbors[0] = -1;
    }

    if(comm_rank == comm_size - 1) {
        neighbors[1] = -1;
    }
    
    //determining the individual lengths;
    int N = 0;
    int localXLength = 0;

    //iterate down to find even split of range;
    for(int i {PROBLEM::Nx}; i >= 1; i--) {
        if(i % comm_size == 0) { //if the range is evenly divisible;
            localXLength = i / comm_size;
            break;
        }
    }

    //determining the starting position;
    int sPos = comm_rank * localXLength;

    //adding the remainder pixels to the last rank;
    if(comm_rank == comm_size-1) {
        localXLength += (PROBLEM::Nx - comm_size * localXLength);
    }

    //printing info for debug;
    if(comm_rank == 0) {
        std::cout << "rank:" << comm_rank << ", length: " << localXLength << ", neighbors: [" << neighbors[0] << ", " << neighbors[1] << "]" << std::endl;
    }
    if(comm_rank == comm_size-1) {
        std::cout << "rank:" << comm_rank << ", length: " << localXLength << ", neighbors: [" << neighbors[0] << ", " << neighbors[1] << "]" << std::endl;
    }

    fieldU   .resize(localXLength, {0,0,0});
    fieldUbar.resize(localXLength, {0,0,0});
    fieldF   .resize(localXLength, {0,0,0});
    fieldFbar.resize(localXLength, {0,0,0});
    fieldJ   .resize(localXLength, {0,0,0});
    fieldJbar.resize(localXLength, {0,0,0});

    stress   .resize(localXLength, {0,0});

    fieldU_ghost[0]    = {0,0,0}; fieldU_ghost[1]    = {0,0,0};
    fieldUbar_ghost[0] = {0,0,0}; fieldUbar_ghost[1] = {0,0,0};
    fieldF_ghost[0]    = {0,0,0}; fieldF_ghost[1]    = {0,0,0};
    fieldFbar_ghost[0] = {0,0,0}; fieldFbar_ghost[1] = {0,0,0};
    fieldJ_ghost[0]    = {0,0,0}; fieldJ_ghost[1]    = {0,0,0};
    fieldJbar_ghost[0] = {0,0,0}; fieldJbar_ghost[1] = {0,0,0};

    stress_ghost[0] = {0,0}; stress_ghost[1] = {0,0};


    // =====================
    // Initializing values
    // =====================

    int stepcount {0};

    //Iterate through locations in x;
    for(int ix {0}; ix < localXLength; ix++){

        //the current position in x
        double cpos {PROBLEM::ORIGINX + (sPos + ix) * PROBLEM::DELTAX};
        
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
        if(ix == localXLength - 1) {nx = ix;}

        //initialize this U, F, J
        fieldU[ix] = obtainU(rho, u, obtainE(T));
        fieldF[ix] = obtainF(fieldU[ix], stress[ix]);
        fieldJ[ix] = obtainJ(fieldU[ix], fieldU[px], fieldU[nx], stress[ix]);
    }

    exchangeFields(comm_size, comm_rank);

    if(comm_rank == 0) {

        std::cout << "=========== Starting ===========" << '\n';
        // writeArrFile(fieldU, "output/initial.txt");
    }

    //Now, we want to while loop up to our final time;
    double ctime {0};

    auto start = std::chrono::high_resolution_clock::now();

    while(ctime < PROBLEM::TIMETARGET){

        //incrementing time upwards.
        //ctime += PROBLEM::TIMESTEP;
        double timestep = PROBLEM::TIMESTEP;
        //using the non-equal timesteps instead;
        // double timestep = obtainTimestep(fieldU);

        //we need to the allreduce to get the minimum timestep across the ranges;

        //if this would put us past the final time, clamp;
        if((ctime + timestep) >= PROBLEM::TIMETARGET){
            timestep = PROBLEM::TIMETARGET - ctime;
            ctime = PROBLEM::TIMETARGET;
        }
        //otherwise, increase ctime.
        else {ctime += timestep;}
        

        //update stress vector to this timestep;
        for(int ix {0}; ix < localXLength; ix++){
            if(ix == 0)                   {stress[ix] = obtainStress(fieldU_ghost[0], fieldU[ix+1]);}
            else if(ix == localXLength-1) {stress[ix] = obtainStress(fieldU[ix-1], fieldU_ghost[1]);}
            else                          {stress[ix] = obtainStress(fieldU[ix-1], fieldU[ix+1]);}
        }

        exchangeStressGhostCells(stress, stress_ghost, neighbors, comm_size, comm_rank);

        std::array<double, 3> fieldU_n;
        std::array<double, 3> fieldU_p;

        //updating F and J at this timestep;
        for(int ix {0}; ix < localXLength; ix++){


            if(ix == 0) {fieldU_n = fieldU_ghost[0];}
            else {fieldU_n = fieldU[ix-1];}

            if(ix == localXLength - 1) {fieldU_p = fieldU_ghost[1];}
            else {fieldU_p = fieldU[ix+1];}

            fieldF[ix] = obtainF(fieldU[ix], stress[ix]);
            fieldJ[ix] = obtainJ(fieldU[ix], fieldU_p, fieldU_n, stress[ix]);
        }

        exchangeFieldGhostCells(fieldF, fieldF_ghost, neighbors, comm_size, comm_rank);
        exchangeFieldGhostCells(fieldJ, fieldJ_ghost, neighbors, comm_size, comm_rank);

        std::array<double, 3> derivF {0,0,0}; 

        //Performing the Predictor step;
        for(int ix{0}; ix < localXLength; ix++){
            //Obtaining the derivative of F

            //Determine BC behavior;
            if(ix == 0) {derivF = backwardDifference(fieldF_ghost[0], fieldF[ix], PROBLEM::DELTAX);}
            else {derivF = backwardDifference(fieldF[ix], fieldF[ix-1], PROBLEM::DELTAX);}

            //Performing predictor step;
            fieldUbar[ix] = predictorStep(fieldU[ix], derivF, fieldJ[ix], timestep);

            //determine BC behavior for U;
            if(ix == 0) {fieldU_n = fieldU_ghost[0];}
            else {fieldU_n = fieldU[ix-1];}

            if(ix == localXLength-1) {fieldU_p = fieldU_ghost[1];}
            else {fieldU_p = fieldU[ix+1];}

            //adding artificial viscosity term;
            std::array<double, 3> visc = artificialViscosity(fieldU[ix], fieldU_p, fieldU_n);
            fieldUbar[ix][0] += visc[0];
            fieldUbar[ix][1] += visc[1];
            fieldUbar[ix][2] += visc[2];

        }

        //updating Ubar.
        exchangeFieldGhostCells(fieldUbar, fieldUbar_ghost, neighbors, comm_size, comm_rank);

        //get predictor stress;
        for(int ix {0}; ix < localXLength; ix++){
            if(ix == 0)                   {stress[ix] = obtainStress(fieldUbar_ghost[0], fieldUbar[ix+1]);}
            else if(ix == localXLength-1) {stress[ix] = obtainStress(fieldUbar[ix-1], fieldUbar_ghost[1]);}
            else                          {stress[ix] = obtainStress(fieldUbar[ix-1], fieldUbar[ix+1]);}
        }

        //updating stress
        exchangeStressGhostCells(stress, stress_ghost, neighbors, comm_size, comm_rank);

        std::array<double, 3> fieldUbar_n;
        std::array<double, 3> fieldUbar_p;

        //get predictor F and J;
        for(int ix{0}; ix < localXLength; ix++){


            if(ix == 0) {fieldUbar_n = fieldUbar_ghost[0];}
            else {fieldUbar_n = fieldUbar[ix-1];}

            if(ix == localXLength - 1) {fieldUbar_p = fieldUbar_ghost[1];}
            else {fieldUbar_p = fieldUbar[ix+1];}

            fieldFbar[ix] = obtainF(fieldUbar[ix], stress[ix]);
            fieldJbar[ix] = obtainJ(fieldUbar[ix], fieldUbar_p, fieldUbar_n, stress[ix]);
        }

        //updating Fbar, Jbar
        exchangeFieldGhostCells(fieldFbar, fieldFbar_ghost, neighbors, comm_size, comm_rank);
        exchangeFieldGhostCells(fieldJbar, fieldJbar_ghost, neighbors, comm_size, comm_rank);


        std::array<double, 3> derivFbar {0,0,0};

        //Now, we do the corrector step;
        for(int ix{0}; ix < localXLength; ix++){

            if(ix == localXLength) {derivFbar = forwardDifference(fieldFbar[ix], fieldFbar_ghost[1], PROBLEM::DELTAX);}
            else {derivFbar = forwardDifference(fieldFbar[ix], fieldFbar[ix+1], PROBLEM::DELTAX);}

            fieldU[ix] = correctorStep(fieldU[ix], fieldUbar[ix], derivFbar, fieldJbar[ix], timestep);

            //determine BC behavior for U;
            if(ix == 0) {fieldUbar_n = fieldUbar_ghost[0];}
            else {fieldUbar_n = fieldUbar[ix-1];}

            if(ix == localXLength-1) {fieldUbar_p = fieldUbar_ghost[1];}
            else {fieldUbar_p = fieldUbar[ix+1];}

            //adding artificial viscosity term;
            std::array<double, 3> visc = artificialViscosity(fieldUbar[ix], fieldUbar_p, fieldUbar_n);
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

    //We now need to combine the data onto rank 0.
    std::vector<std::array<double, 3>> fullFieldU; fullFieldU.resize(PROBLEM::Nx, {0,0,0});
    std::vector<double> fullFieldU_e1, fieldU_e1; fullFieldU_e1.resize(PROBLEM::Nx, 0.0); fieldU_e1.resize(localXLength, 0.0); 
    std::vector<double> fullFieldU_e2, fieldU_e2; fullFieldU_e2.resize(PROBLEM::Nx, 0.0); fieldU_e2.resize(localXLength, 0.0); 
    std::vector<double> fullFieldU_e3, fieldU_e3; fullFieldU_e3.resize(PROBLEM::Nx, 0.0); fieldU_e3.resize(localXLength, 0.0); 

    //need to rearrange this to be conducive to sending via gather...
    //   from vector of arrays to three 1 dim vectors;
    for(int i = 0; i < fieldU.size(); i++) {
        fieldU_e1[i] = fieldU[i][0];
        fieldU_e2[i] = fieldU[i][1];
        fieldU_e3[i] = fieldU[i][2];
    }

    // for(int i = 0; i < fieldU.size(); i++) {
    //     std::cout << fieldU[i][0] << ", ";
    // }
    // std::cout << std::endl;

    MPI_Gather(fieldU_e1.data(), localXLength, MPI_DOUBLE,
               fullFieldU_e1.data(), localXLength, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(fieldU_e2.data(), localXLength, MPI_DOUBLE,
               fullFieldU_e2.data(),  localXLength, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(fieldU_e3.data(), localXLength, MPI_DOUBLE,
               fullFieldU_e3.data(), localXLength, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    //printing the final timestep.
    if(comm_rank == 0) {
        std::cout << "time elapsed: " << elapsed.count() << " seconds" << std::endl;
        std::cout << "it/sec: " << stepcount / elapsed.count() << std::endl;

        //collecting fullFieldU;
        for(int i = 0; i < PROBLEM::Nx; i++) {
            fullFieldU[i][0] = fieldU_e1[i];
            fullFieldU[i][1] = fieldU_e2[i];
            fullFieldU[i][2] = fieldU_e3[i];
        }

        //writing out the final result to a text file.
        writeArrFile(fullFieldU, "output/final.txt");
        std::cout << "=========== Done! ===========" << '\n';
    }


    MPI_Finalize();
    return 0;
}

