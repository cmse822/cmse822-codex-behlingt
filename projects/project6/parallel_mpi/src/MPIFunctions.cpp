#include <cstdlib>

#include <vector>
#include <array>

#include <iostream>
#include <fstream>
#include <string>

#include <chrono>

#include <mpi.h>

std::array<double, 3> negFieldBufferSend {0,0};
std::array<double, 3> posFieldBufferSend {0,0};

std::array<double, 3> negFieldBufferRecv {0,0};
std::array<double, 3> posFieldBufferRecv {0,0};


std::array<double, 2> negStressBufferSend {0,0};
std::array<double, 2> posStressBufferSend {0,0};

std::array<double, 2> negStressBufferRecv {0,0};
std::array<double, 2> posStressBufferRecv {0,0};

std::array<MPI_Request, 2> requests;

/**
 * @brief exchanges Ghost Cells for Field objects.
 */
void exchangeFieldGhostCells(std::vector<std::array<double, 3>> &fieldX, 
                        std::array<std::array<double, 3>, 2> &fieldX_ghost, 
                        std::array<int, 2> neighbors, int comm_size, int comm_rank){

    // std::cout << "[" << comm_rank << "] exchanging..." << std::endl;

    std::fill(requests.begin(), requests.end(), MPI_REQUEST_NULL);

    //for each neighbor, we want to exchange edge point with friends;
    int reqCount = 0;

    //If we have a valid negative neighbor
    //  Exchange negative ghost point
    if(neighbors[0] >= 0) {
        //post receive
        MPI_Irecv(negFieldBufferRecv.data(), 3, MPI_DOUBLE, 
                  neighbors[0], 0, MPI_COMM_WORLD, &requests[reqCount++]);

        negFieldBufferSend = fieldX[0]; //preparing buffer

        MPI_Isend(negFieldBufferSend.data(), 3, MPI_DOUBLE,
                  neighbors[0], 0, MPI_COMM_WORLD, &requests[reqCount++]);
    }

    //If we have a valid positive neighbor
    //  Exchange positive ghost point
    if(neighbors[1] >= 0) {
        //post receive
        MPI_Irecv(posFieldBufferRecv.data(), 3, MPI_DOUBLE, 
                  neighbors[1], 0, MPI_COMM_WORLD, &requests[reqCount++]);

        posFieldBufferSend = fieldX[fieldX.size()-1]; //preparing buffer

        MPI_Isend(posFieldBufferSend.data(), 3, MPI_DOUBLE,
                  neighbors[1], 0, MPI_COMM_WORLD, &requests[reqCount++]);
    }

    // std::cout << "[" << comm_rank << "] made it to waitall" << std::endl;

    //Wait for requests to arrive;
    MPI_Waitall(reqCount, requests.data(), MPI_STATUSES_IGNORE);

    // std::cout << "[" << comm_rank << "] past waitall" << std::endl;

    //unpack buffers;
    if(neighbors[0] >= 0) { //Negative X Buffer (left side);
        fieldX_ghost[0] = negFieldBufferRecv;
        // std::cout << "[" << comm_rank << "] negative successful" << std::endl;
    }
    else { //sneaky trick; if this is the edge, we apply boundary conditions instead;
        fieldX_ghost[0] = fieldX[0];
    }
    if(neighbors[1] >= 0) { //Positive X Buffer (right side);
        fieldX_ghost[1] = posFieldBufferRecv;
        // std::cout << "[" << comm_rank << "] positive successful" << std::endl;
    }
    else { //sneaky trick; if this is the edge, we apply boundary conditions instead;
        fieldX_ghost[1] = fieldX[fieldX.size()-1];
    }

    // std::cout << "[" << comm_rank << "] the end" << std::endl;
}

/**
 * @brief exchanges Ghost Cells for stress objects.
 */
void exchangeStressGhostCells(std::vector<std::array<double, 2>> &stress, 
    std::array<std::array<double, 2>, 2> &stress_ghost, 
    std::array<int, 2> neighbors, int comm_size, int comm_rank){

    std::fill(requests.begin(), requests.end(), MPI_REQUEST_NULL);

    //for each neighbor, we want to exchange edge point with friends;
    int reqCount = 0;

    //If we have a valid negative neighbor
    //  Exchange negative ghost point
    if(neighbors[0] >= 0) {
        //post receive
        MPI_Irecv(negStressBufferRecv.data(), 2, MPI_DOUBLE, 
        neighbors[0], 0, MPI_COMM_WORLD, &requests[reqCount++]);

        negStressBufferSend = stress[0]; //preparing buffer

        MPI_Isend(negStressBufferSend.data(), 2, MPI_DOUBLE,
        neighbors[0], 0, MPI_COMM_WORLD, &requests[reqCount++]);
    }

    //If we have a valid positive neighbor
    //  Exchange positive ghost point
    if(neighbors[1] >= 0) {
    //post receive
    MPI_Irecv(posStressBufferRecv.data(), 2, MPI_DOUBLE, 
        neighbors[1], 0, MPI_COMM_WORLD, &requests[reqCount++]);

        posStressBufferSend = stress[stress.size()-1]; //preparing buffer

        MPI_Isend(posStressBufferSend.data(), 2, MPI_DOUBLE,
        neighbors[1], 0, MPI_COMM_WORLD, &requests[reqCount++]);
    }

    //Wait for requests to arrive;
    MPI_Waitall(reqCount, requests.data(), MPI_STATUSES_IGNORE);

    //unpack buffers;
    if(neighbors[0] >= 0) { //Negative X Buffer (left side);
        stress_ghost[0] = negStressBufferRecv;
    }
    else { //sneaky trick; if this is the edge, we apply boundary conditions instead;
        stress_ghost[0] = stress[0];
    }
    if(neighbors[1] >= 0) { //Positive X Buffer (right side);
        stress_ghost[1] = posStressBufferRecv;
    }
    else { //sneaky trick; if this is the edge, we apply boundary conditions instead;
        stress_ghost[1] = stress[stress.size()-1];
    }
}