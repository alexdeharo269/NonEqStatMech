#include <cmath>
#include <limits>
#include <random>
#include <utility>
#include <vector>
#include "P1-funcs.h"

const double gamma = 1.0;
const double dt = 0.01;
const int L = 100;
const int Nmols = 1000;
const int tau = 100;

int main()
{
    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<> runif(0.0, L);

    // Open output file for D vs Gamma data
    FILE *outfile = fopen("P1_5-data.txt", "w");
    if (!outfile)
    {
        printf("Error opening file!\n");
        return 1;
    }

    fprintf(outfile, "# Gamma D_noPBC D_PBC\n");

    // Test with several Gamma values
    std::vector<double> Gammas = {0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100, 300, 1000};

    for (auto Gamma : Gammas)
    {

        // === Simulation WITHOUT PBC ===
        std::vector<std::pair<double, double>> pos_noPBC;
        std::vector<std::pair<double, double>> initial_pos_noPBC;

        // Initialize positions
        for (int i = 0; i < Nmols; i++)
        {
            double x = runif(rng);
            double y = runif(rng);
            pos_noPBC.push_back(std::make_pair(x, y));
            initial_pos_noPBC.push_back(std::make_pair(x, y));
        }

        // Run simulation without PBC
        int steps = (int)(tau / dt);
        for (int i = 0; i < Nmols; i++)
        {
            auto [x, y] = initial_pos_noPBC[i];
            for (int step = 0; step < steps; step++)
            {
                auto [xupt, yupt] = diffuse(x, y, Gamma, dt);
                x = xupt;
                y = yupt;
            }
            pos_noPBC[i] = std::make_pair(x, y);
        }

        // Compute MSD without PBC
        double msd_noPBC = 0.0;
        for (int i = 0; i < Nmols; i++)
        {
            auto [x0, y0] = initial_pos_noPBC[i];
            auto [xf, yf] = pos_noPBC[i];
            double dx = xf - x0;
            double dy = yf - y0;
            msd_noPBC += dx * dx + dy * dy;
        }
        double D_noPBC = msd_noPBC / (4.0 * tau * Nmols);

        // === Simulation WITH PBC ===
        std::vector<std::tuple<double, double, int, int>> pos_PBC;
        std::vector<std::pair<double, double>> initial_pos_PBC;

        // Initialize positions
        for (int i = 0; i < Nmols; i++)
        {
            double x = runif(rng);
            double y = runif(rng);
            pos_PBC.push_back(std::make_tuple(x, y, 0, 0));
            initial_pos_PBC.push_back(std::make_pair(x, y));
        }

        // Run simulation with PBC
        for (int step = 0; step < steps; step++)
        {
            for (int j = 0; j < Nmols; j++)
            {
                auto [x, y, nx, ny] = pos_PBC[j];

                auto [xupt, yupt] = diffuse(x, y, Gamma, dt);

                // Apply PBC and track crossings
                if (xupt < 0)
                {
                    xupt += L;
                    nx -= 1;
                }
                else if (xupt >= L)
                {
                    xupt -= L;
                    nx += 1;
                }

                if (yupt < 0)
                {
                    yupt += L;
                    ny -= 1;
                }
                else if (yupt >= L)
                {
                    yupt -= L;
                    ny += 1;
                }

                pos_PBC[j] = std::make_tuple(xupt, yupt, nx, ny);
            }
        }

        // Compute MSD with PBC (unwrapped coordinates)
        double msd_PBC = 0.0;
        for (int i = 0; i < Nmols; i++)
        {
            auto [x0, y0] = initial_pos_PBC[i];
            auto [xf, yf, nx, ny] = pos_PBC[i];

            // Unwrap final position
            double xf_unwrapped = xf + nx * L;
            double yf_unwrapped = yf + ny * L;

            // Compute squared displacement
            double dx = xf_unwrapped - x0;
            double dy = yf_unwrapped - y0;
            msd_PBC += dx * dx + dy * dy;
        }
        double D_PBC = msd_PBC / (4.0 * tau * Nmols);

        // Write results to file
        fprintf(outfile, "%.6f %.6f %.6f\n", Gamma, D_noPBC, D_PBC);
        printf("Gamma = %.2f: D_noPBC = %.6f, D_PBC = %.6f\n", Gamma, D_noPBC, D_PBC);
    }

    fclose(outfile);
    printf("\nDone! Results written to P1_5-data.txt\n");

    return 0;
}