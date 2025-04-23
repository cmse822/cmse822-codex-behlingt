#include <iostream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <omp.h>
#include "mm_utils.hpp"

// Constants
constexpr double TOLERANCE = 0.001; //bumped this from 0.005; I kept getting errors of ~0.001
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

    //need to transpose this A. :)
    std::vector<double> Atemp(Ndim * Ndim);
    for(int i = 0; i < Ndim; i++){
        for(int j = 0; j < Ndim; j++) {
            Atemp[index(i,j, Ndim)] = A[index(j,i, Ndim)];
        }
    }
    std::swap(A, Atemp);

    // Pointers to the data;
    //   needed for the OMP target offload; doesn't support vector class.
    //   OMP <-- strange??? :/
    double* A_ptr {A.data()};
    double* b_ptr {b.data()};
    double* xnew_ptr {xnew.data()};
    double* xold_ptr {xold.data()};


    // Initialize b with random values (between 0.0 and 0.5)
    for (int i = 0; i < Ndim; ++i) {
        b[i] = static_cast<double>(std::rand() % 51) / 100.0;
    }


    double start_time = omp_get_wtime();
    double conv = LARGE;
    int iters = 0;

    std::cout << "Starting loop..." << std::endl;

    #pragma omp target data map(tofrom: xnew_ptr[0:Ndim], xold_ptr[0:Ndim]) \
                            map(to: A_ptr[0:Ndim*Ndim], b_ptr[0:Ndim])
    {
    while ((conv > TOLERANCE) && (iters < MAX_ITERS)) {
        ++iters;

        // Compute new iteration
        // for each element in x;
        // #pragma omp target map(tofrom: xnew_ptr[0:Ndim], xold_ptr[0:Ndim]) \
        //                    map(to: A_ptr[0:Ndim*Ndim], b_ptr[0:Ndim])
        // #pragma omp loop
        #pragma omp target loop
        for (int i = 0; i < Ndim; ++i) {
            xnew_ptr[i] = 0.0;

            //we want to add together all the i != j elements of A*xold[j].
            for (int j = 0; j < Ndim; ++j) {
                // if (i != j)
                //     xnew_ptr[i] += A_ptr[index(i,j,Ndim)] * xold_ptr[j]; 
                xnew_ptr[i] += A_ptr[index(j,i,Ndim)] * xold_ptr[j] * static_cast<double>(i != j);
            }
            //adding bi, and dividing by aii
            xnew_ptr[i] = (b_ptr[i] - xnew_ptr[i]) / (A_ptr[index(i,i,Ndim)]);
        }

        // Compute convergence criterion (Euclidean norm of difference)
        conv = 0.0;
        // #pragma omp target //map(to: xnew_ptr[0:Ndim], xold_ptr[0:Ndim]) \
        //                    map(tofrom: conv)
        #pragma omp target loop map(tofrom: conv) reduction(+:conv)
        // #pragma omp target loop reduction(+:conv)
        for (int i = 0; i < Ndim; ++i) {
            double tmp = xnew_ptr[i] - xold_ptr[i];
            conv += tmp * tmp;
        }
        conv = std::sqrt(conv);

        if((iters % OUTINTERVAL) == 0) {
            std::cout << "it:" << iters << ", conv:" << conv << std::endl;
        }

        // Swap vectors for next iteration
        // std::swap(xold, xnew);
        #pragma omp target loop
        for(int i = 0; i < Ndim; i++) {
            xold_ptr[i] = xnew_ptr[i];
        }
    }
    }

    double elapsed_time = omp_get_wtime() - start_time;
    std::cout << "Converged after " << iters << " iterations in "
              << elapsed_time << " seconds with final convergence = "
              << conv << std::endl;

    double err = 0.0, chksum = 0.0;

    for (int i = 0; i < Ndim; ++i) {
        xold[i] = static_cast<double>(0.0);
        for (int j = 0; j < Ndim; ++j)
            xold[i] += A[j * Ndim + i] * xnew[j];

        double diff = xold[i] - b[i];
        chksum += xnew[i];
        err += diff * diff;
    }
    err = static_cast<double>(std::sqrt(err));

    double totalFlops = 2.0 * static_cast<double>(Ndim)*static_cast<double>(Ndim)*static_cast<double>(iters);
    double flopsPerSecond = (elapsed_time > 0) ? totalFlops / elapsed_time : 0.0;
    std::cout << "FLOP rate: " << flopsPerSecond << " FLOPS/s" << std::endl;
    
    std::cout << "Solution verification: Error = " << err
            << ", Checksum = " << chksum << std::endl;

    if (err > TOLERANCE)
        std::cout << "WARNING: Solution error exceeds tolerance!" << std::endl;

    return 0;
}
