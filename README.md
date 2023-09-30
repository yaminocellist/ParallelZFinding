# zFinding using ROOT::TH1
This .c script can be easily run in ROOT. The main strategy for optimization is to reduce the number of 'if' statements.

# ParallelZFinding
On macOS:

g++-11 -std=c++17 -fopenmp Mocktilt_plan2_para_zfinding.cpp

export OMP_NUM_THREADS=10 (input the number of cores of your mac's CPU)

./a.out 0
