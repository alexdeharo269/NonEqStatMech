/*
Create an Ornstein-Uhlenbeck process simulation in C++17 in 1D, using the Euler-Maruyama method. 
Simulates N particles and computes their Mean Square Displacement (MSD).
*/

#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <numeric>
#include <fstream>
#include <tuple>
#include <chrono>
#include "P2-funcs.h"

#define M_PI 3.14159265358979323846
const double tau=1.0;
const double dt=0.01;
const int Nsteps=1000;
const double D=1.0;
const double L=10.0;  // Box length for PBC
const int Nparticles=500; // Number of particles

int main()
{
    
    // Initialize particles
    std::vector<double> x(Nparticles, 0.0);  // Positions
    std::vector<double> xi(Nparticles, 0.0); // Noise variables
    std::vector<double> x0(Nparticles, 0.0); // Initial positions for MSD
    
    // Open output files
    std::ofstream posfile("P2-1-positions.txt");
    std::ofstream msdfile("P2-1-msd.txt");
    
    
    // 1) Single timed run. 
    // Start timing
    auto start = std::chrono::high_resolution_clock::now();
    run(Nsteps, Nparticles, dt, tau, D, L, true ,posfile,0, msdfile);
    auto end = std::chrono::high_resolution_clock::now();

    

    posfile.close();
    msdfile.close();
    
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Simulation time: " << duration.count() << " ms" << std::endl;
    std::cout << "Simulated " << Nparticles << " particles for " << Nsteps << " steps" << std::endl;
    
    return 0;
}