// Now we compute the correlation function C(t) = <v(0)v(t)> for different tau values and compare with theory.

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
const double tau=10.0;

int main()
{

    // Open output files
    std::ofstream posfile("P2-4-pos.txt");
    std::ofstream msdfile("P2-4-corr.txt");
    std::ofstream disfile("P2-4-dis.txt");

    // 3) Repeat for multiple runs with different tau parameters

    for (double force : std::vector<double>{0.0,0.001,0.01,0.1,1})
    {

        // Run simulation

        run_with_force(Nsteps, Nparticles, dt, tau, D, L, force, posfile, msdfile, disfile);

        posfile << std::endl;
        msdfile << std::endl;
        disfile << std::endl;
    }
    // Erase pos file, close msd file
    posfile.close();
    msdfile.close();
    disfile.close();

    fprintf(stdout, "Simulated %d particles for %d steps with varying tau\n", Nparticles, Nsteps);
    return 0;
}