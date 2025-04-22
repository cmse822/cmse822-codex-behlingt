#include <iostream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <omp.h>
#include "mm_utils.hpp"

// Constants
constexpr double TOLERANCE = 0.001;
constexpr int DEF_SIZE = 1000;
constexpr int MAX_ITERS = 100000;
constexpr double LARGE = 1000000.0;
constexpr int OUTINTERVAL = 500;

//very funny...
//  why is the matrix indexed backwards???
int index(int i, int j, int ndim) {
    return(j + i * ndim);
}

int main(int argc, char **argv) {
    int Ndim = (argc == 2) ? std::atoi(argv[1]) : DEF_SIZE;
    std::cout << "Matrix dimension (Ndim) = " << Ndim << std::endl;

    // Use std::vector to allocate storage dynamically.
    std::vector<double> A(Ndim * Ndim);
    std::vector<double> b(Ndim);
    std::vector<double> xnew(Ndim, 0.0);
    std::vector<double> xold(Ndim, 0.0);

    // Generate diagonally dominant matrix A
    initDiagDomNearIdentityMatrix(Ndim, A.data());

    // Initialize b with random values (between 0.0 and 0.5)
    for (int i = 0; i < Ndim; ++i) {
        b[i] = static_cast<double>(std::rand() % 51) / 100.0;
    }


    double start_time = omp_get_wtime();
    double conv = LARGE;
    int iters = 0;

    std::cout << "Starting loop..." << std::endl;

    while ((conv > TOLERANCE) && (iters < MAX_ITERS)) {
        ++iters;

        // Compute new iteration
        // for each element in x;
        for (int i = 0; i < Ndim; ++i) {
            xnew[i] = 0.0;

            //we want to add together all the i != j elements of A*xold[j].
            for (int j = 0; j < Ndim; ++j) {
                if (i != j)
                    xnew[i] += A[index(i,j,Ndim)] * xold[j]; 
            }
            //adding bi, and dividing by aii
            xnew[i] = (b[i] - xnew[i]) / (A[index(i,i,Ndim)]);
        }

        // Compute convergence criterion (Euclidean norm of difference)
        conv = 0.0;

        for (int i = 0; i < Ndim; ++i) {
            double tmp = xnew[i] - xold[i];
            conv += tmp * tmp;
        }
        conv = std::sqrt(conv);

        if((iters % OUTINTERVAL) == 0) {
            std::cout << "it:" << iters << ", conv:" << conv << std::endl;
        }

        // Swap vectors for next iteration
        std::swap(xold, xnew);
    }

    double elapsed_time = omp_get_wtime() - start_time;
    std::cout << "Converged after " << iters << " iterations in "
              << elapsed_time << " seconds with final convergence = "
              << conv << std::endl;

    double err = 0.0, chksum = 0.0;

    for (int i = 0; i < Ndim; ++i) {
        xold[i] = static_cast<double>(0.0);
        for (int j = 0; j < Ndim; ++j)
            xold[i] += A[i * Ndim + j] * xnew[j];

        double diff = xold[i] - b[i];
        chksum += xnew[i];
        err += diff * diff;
    }
    err = static_cast<double>(std::sqrt(err));

    std::cout << "Solution verification: Error = " << err
            << ", Checksum = " << chksum << std::endl;

    if (err > TOLERANCE)
        std::cout << "WARNING: Solution error exceeds tolerance!" << std::endl;

    return 0;
}
