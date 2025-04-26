#include <cstdlib>

#include <vector>
#include <array>

#include <iostream>

/**
 * @brief exchanges Ghost Cells for Field objects.
 */
void exchangeFieldGhostCells(std::vector<std::array<double, 3>> &fieldX, 
    std::array<std::array<double, 3>, 2> &fieldX_ghost, 
    std::array<int, 2> neighbors, int comm_size, int comm_rank);

/**
 * @brief exchanges Ghost Cells for stress objects.
 */
void exchangeStressGhostCells(std::vector<std::array<double, 2>> &stress, 
    std::array<std::array<double, 2>, 2> &stress_ghost, 
    std::array<int, 2> neighbors, int comm_size, int comm_rank);