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
const double D = 1.0;
const double L = 10.0;       // Box length for PBC
const int Nparticles = 2000; // Number of particles

int main()
{
    // Open output files
    std::ofstream posfile("P2-3-pos-dummy.txt");
    std::ofstream corrfile("P2-3-corr.txt");

    std::vector<double> taus = {1.0, 10.0, 100.0};

    for (double tau : taus)
    {
        // CRITICAL: Adjust Nsteps based on tau
        // Need at least: transient (3*tau) + correlation measurement (5*tau) = 8*tau
        // Add extra factor for multiple snapshots: 12*tau
        int Nsteps = static_cast<int>(std::ceil(12.0 * tau / dt));

        std::cout << "Running simulation for tau = " << tau << std::endl;
        std::cout << "  Nsteps = " << Nsteps << " (" << (Nsteps * dt) << " time units)" << std::endl;

        // Run simulation
        run(Nsteps, Nparticles, dt, tau, D, L, false, posfile, 3, corrfile);

        // Add separator between different tau runs
        corrfile << std::endl;
    }

    // Clean up
    posfile.close();
    std::remove("P2-3-pos-dummy.txt"); // Delete the temporary position file
    corrfile.close();

    std::cout << "Simulation complete!" << std::endl;
    std::cout << "Results written to P2-3-corr.txt" << std::endl;

    return 0;
}