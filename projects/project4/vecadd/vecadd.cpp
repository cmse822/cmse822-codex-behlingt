#include <iostream>
#include <vector>
#include <chrono>
#include <omp.h>

//module load NVHPC/23.7-CUDA-12.1.1
//nvc++ -O2 -mp=multicore -o vecadd_cpu vecadd.cpp 
//nvc++ -O2 -mp=gpu -gpu=ccnative -o vecadd_gpu vecadd.cpp 

//PLOTS TO MAKE; (Jacobian, if it doesn't work use Vector Addition)
// - FLOPS/S vs CPU Threads for different vector sizes (BIIIG factors of 2, 1048576...)
// - FLOPS/S vs Vector Size, GPU vs CPU comparison


// constexpr int N = 10000;
constexpr int N = 10000000;
constexpr double TOL = 1e-4;

int main() {
    std::vector<double> a(N), b(N), c(N), res(N);

    // Fill the arrays using simple for loops.
    for (int i = 0; i < N; ++i) {
        a[i] = static_cast<double>(i);
        b[i] = 2.0f * static_cast<double>(i);
        res[i] = a[i] + b[i];
        c[i] = 0.0;
    }

    double* a_ptr {a.data()};
    double* b_ptr {b.data()};
    double* c_ptr {c.data()};

    // Add two vectors using a simple loop and measure time.
    auto start = std::chrono::high_resolution_clock::now();

    // // Start parallelism
    // // A & B are shared]
    // #pragma omp parallel for shared(a, b, c)
    // for (int i = 0; i < N; ++i) {
    //     c[i] = a[i] + b[i];
    // }
    
    //need to map data from host to device
    //  we can specify tofrom, or only send to, pull from.
    {
    #pragma omp target loop map(to: a_ptr[0:N], b_ptr[0:N]) \
                            map(from: c_ptr[0:N])
    for (int i = 0; i < N; ++i) {
        c_ptr[i] = a_ptr[i] + b_ptr[i];
    }
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    // Test results.
    int err = 0;    
    for (int i = 0; i < N; ++i) {
        double diff = c[i] - res[i];
        if (diff * diff > TOL) {
            ++err;
        }
    }

    std::cout << "Vectors added with " << err << " errors\n";
    std::cout << "Elapsed time: " << elapsed.count() << " seconds\n";
    double flops = 2.0 * N / elapsed.count();
    std::cout << "FLOP rate: " << flops << " FLOP/s\n";
    return 0;
}
