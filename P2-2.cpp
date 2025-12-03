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
const double dt = 0.1;
const int Nsteps = 1000;
const double D = 1.0;
const double L = 10.0;      // Box length for PBC
const int Nparticles = 500; // Number of particles

int main()
{


    // Open output files
    std::ofstream posfile("P2-2-pos-dummy.txt");
    std::ofstream msdfile("P2-2-msd.txt");

    // 2) Repeat for multiple runs with different tau parameters
    for (double tau : {1.0, 10.0, 100.0})
    {
        // Run simulation
        run(Nsteps, Nparticles, dt, tau, D, L,false, posfile,2, msdfile);
        // Add separator between different tau runs

        msdfile << std::endl;
    }
    // Erase pos file, close msd file
    posfile.close();
    std::remove("P2-2-pos-dummy.txt"); // Delete the temporary position file
    msdfile.close();

    fprintf(stdout, "Simulated %d particles for %d steps with varying tau\n", Nparticles, Nsteps);
    return 0;
}