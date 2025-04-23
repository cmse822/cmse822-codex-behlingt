module load NVHPC/23.7-CUDA-12.1.1

nvc++ -O2 -o p6_serial src/main.cpp src/MathFunctions.cpp --diag_suppress set_but_not_used