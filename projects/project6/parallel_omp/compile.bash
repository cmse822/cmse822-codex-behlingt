module load NVHPC/23.7-CUDA-12.1.1

nvc++ -O2 -mp=multicore -o p6_omp_cpu src/main.cpp src/MathFunctions.cpp --diag_suppress set_but_not_used
nvc++ -O2 -mp=gpu -gpu=ccnative -o p6_omp_gpu src/main.cpp src/MathFunctions.cpp --diag_suppress set_but_not_used