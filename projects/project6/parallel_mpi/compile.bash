# module load NVHPC/23.7-CUDA-12.1.1
module load OpenMPI

mpic++ -O2 -o p6_mpi src/main.cpp src/MathFunctions.cpp src/MPIFunctions.cpp