#include <cmath>
#include <limits>
#include <random>
#include <utility>

#include <fstream>
#include <iostream>
#define M_PI 3.14159265358979323846


//"mu" is the mean of the distribution, and "sigma" is the standard deviation.
inline std::pair<double, double> generateGaussianNoise(double mu, double sigma)
{
    constexpr double two_pi = 2.0 * M_PI;

    // initialize the random uniform number generator (runif) in a range 0 to 1
    static std::mt19937 rng(std::random_device{}()); // Standard mersenne_twister_engine seeded with rd()
    static std::uniform_real_distribution<> runif(0.0, 1.0);

    // create two random numbers, make sure u1 is greater than zero
    double u1, u2;
    do
    {
        u1 = runif(rng);
    } while (u1 == 0);
    u2 = runif(rng);

    // compute z0 and z1
    auto mag = sigma * sqrt(-2.0 * log(u1));
    auto z0 = mag * cos(two_pi * u2) + mu;
    auto z1 = mag * sin(two_pi * u2) + mu;

    return std::make_pair(z0, z1);
}


inline double wrapPosition(double x, double L)
{
    while (x < 0)
        x += L;
    while (x >= L)
        x -= L;
    return x;
}

// -----------------------------------------------------------------------------
// minimal-change implementation: circular buffer + stationary init + correct theory
// -----------------------------------------------------------------------------

// Compute mean squared displacement
inline double computeMSD(const std::vector<double> &x_unwrapped,
                         const std::vector<double> &x0, int Nparticles)
{
    double msd = 0.0;
    for (int p = 0; p < Nparticles; ++p)
    {
        double displacement = x_unwrapped[p] - x0[p];
        msd += displacement * displacement;
    }
    return msd / static_cast<double>(Nparticles);
}






// Write positions to file
inline void writePositions(std::ofstream &posfile, const std::vector<double> &x, int Nparticles)
{
    for (int p = 0; p < Nparticles; ++p)
    {
        posfile << x[p];
        if (p < Nparticles - 1)
            posfile << " ";
    }
    posfile << std::endl;
}


// Update single particle state
inline void updateParticle(double &x, double &x_unwrapped, double &xi,
                           double dt, double tau, double D, double L, bool pbc)
{
    // Generate Gaussian random number
    auto [z0, z1] = generateGaussianNoise(0.0, 1.0);

    // Update position (unwrapped for MSD)
    x_unwrapped += xi * dt;
    x = x_unwrapped;

    // Apply periodic boundary conditions
    if (pbc)
    {
        x = x - L * floor(x / L);
    }

    // Update Ornstein-Uhlenbeck noise
    xi = xi * (1.0 - dt / tau) + (std::sqrt(2.0 * D * dt) / tau) * z0;
}

// Compute correlation for a specific lag
double computeCorrelationAtLag(
    const std::vector<double> &xi_hist,
    int Nparticles,
    int buf_len,
    int current_step,
    int lag)
{
    if (current_step < lag)
        return 0.0;

    int idx_now = current_step % buf_len;
    int idx_old = (current_step - lag) % buf_len;

    double corr = 0.0;
    for (int p = 0; p < Nparticles; ++p)
    {
        double xi_now = xi_hist[p * buf_len + idx_now];
        double xi_old = xi_hist[p * buf_len + idx_old];
        corr += xi_now * xi_old;
    }
    return corr / static_cast<double>(Nparticles);
}

// Initialize particle states
inline void initializeParticles(std::vector<double> &x, std::vector<double> &x_unwrapped,
                                std::vector<double> &xi, std::vector<double> &xi_hist,
                                int Nparticles, int buf_len, double D, double tau)
{
    double xi_std = std::sqrt(D / tau);
    for (int p = 0; p < Nparticles; ++p)
    {
        auto pr = generateGaussianNoise(0.0, xi_std);
        xi[p] = pr.first;

        // Initialize buffer with independent samples
        for (int j = 0; j < buf_len; ++j)
        {
            auto noise = generateGaussianNoise(0.0, xi_std);
            xi_hist[p * buf_len + j] = noise.first;
        }
    }
}

inline void writeExerciseOutput(std::ofstream &specific_file, int exercise, int i, int lag_max,
                                double dt, double D, double tau,
                                const std::vector<double> &xi, const std::vector<double> &xi_hist,
                                const std::vector<double> &x_unwrapped, const std::vector<double> &x0,
                                int Nparticles, int buf_len, int transient_steps)
{
    if (exercise == 2 || exercise == 0)
    {
        double msd = computeMSD(x_unwrapped, x0, Nparticles);
        specific_file << (i * dt) << " " << msd << std::endl;
    }
    else if (exercise == 3)
    {
        // Only output after transient period
        if (i < transient_steps)
            return;

        // Compute every 500 steps to get good statistics but not too much data
        int output_interval = 500;
        if ((i - transient_steps) % output_interval == 0)
        {
            // Sample lags logarithmically for better visualization
            // This gives more points at early times where decay is fast
            std::vector<int> lags_to_sample;

            // Early times: every dt
            for (int lag = 0; lag <= 50 && lag <= lag_max; lag++)
                lags_to_sample.push_back(lag);

            // Medium times: every 5 dt
            for (int lag = 55; lag <= 200 && lag <= lag_max; lag += 5)
                lags_to_sample.push_back(lag);

            // Late times: every 20 dt
            for (int lag = 220; lag <= lag_max; lag += 20)
                lags_to_sample.push_back(lag);

            for (int lag : lags_to_sample)
            {
                double corr = computeCorrelationAtLag(xi_hist, Nparticles, buf_len, i, lag);
                double theory = (D / tau) * std::exp(-(lag * dt) / tau);
                specific_file << (lag * dt) << " " << corr << " " << theory << " " << tau << std::endl;
            }
            specific_file << std::endl; // Blank line between snapshots
        }
    }
}

// Main simulation loop
inline void run(int Nsteps, int Nparticles, double dt, double tau, double D, double L,
                bool pbc, std::ofstream &posfile, int exercise, std::ofstream &specific_file)
{
    // Need to track correlations up to ~5*tau for full decay
    int max_lag = static_cast<int>(std::ceil(5.0 * tau / dt));
    int buf_len = max_lag + 1;

    // Wait for system to equilibrate (at least 3*tau)
    int transient_steps = static_cast<int>(std::ceil(3.0 * tau / dt));

    std::cout << "Simulation parameters:" << std::endl;
    std::cout << "  tau = " << tau << std::endl;
    std::cout << "  dt = " << dt << std::endl;
    std::cout << "  max_lag = " << max_lag << " steps = " << (max_lag * dt) << " time units" << std::endl;
    std::cout << "  buf_len = " << buf_len << std::endl;
    std::cout << "  transient_steps = " << transient_steps << std::endl;
    std::cout << "  Nsteps = " << Nsteps << std::endl;

    // Initialize state vectors
    std::vector<double> x(Nparticles, 0.0);
    std::vector<double> x_unwrapped(Nparticles, 0.0);
    std::vector<double> xi(Nparticles, 0.0);
    std::vector<double> xi_hist(Nparticles * buf_len, 0.0);

    initializeParticles(x, x_unwrapped, xi, xi_hist, Nparticles, buf_len, D, tau);

    std::vector<double> x0 = x_unwrapped;

    // Time evolution
    for (int i = 0; i < Nsteps; ++i)
    {
        int idx_now = i % buf_len;

        // Update all particles
        for (int p = 0; p < Nparticles; ++p)
        {
            updateParticle(x[p], x_unwrapped[p], xi[p], dt, tau, D, L, pbc);
            xi_hist[p * buf_len + idx_now] = xi[p];
        }

        if (exercise != 3) // Only write positions if not doing correlation
            writePositions(posfile, x, Nparticles);

        writeExerciseOutput(specific_file, exercise, i, max_lag, dt, D, tau,
                            xi, xi_hist, x_unwrapped, x0, Nparticles, buf_len, transient_steps);
    }
}

inline void run_with_force(int Nsteps, int Nparticles, double dt, double tau, double D, double L,
                            double force, std::ofstream &posfile, std::ofstream &msd_file, std::ofstream &disfile)
{
    int runs=10;    

    std::vector<double> msd_t(Nsteps,0.0);
    std::vector<double> dis_t(Nsteps,0.0);

    // Intialize force field
    std::mt19937 rng(std::random_device{}()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> runif(0.0, 1.0);

    //for different run seeds
    for (int run=0;run<runs;run++){

        // Initialize state vectors
        std::vector<double> x(Nparticles, 0.0);
        std::vector<double> xi(Nparticles, 0.0);

        
        std::vector<int> eps(Nparticles,0);
        for (int part=0;part<Nparticles;part++){
            double r=runif(rng);
            if(r<0.5){
                eps[part]=-1;
            }else{
                eps[part]=1;
            }
        }
        // El MSD y el desplazamiento lo computo como media de semillas y lo saco cada tiempo.
        // Time evolution

        
        for (int i = 0; i < Nsteps; ++i)
        {
            
            for (int p = 0; p < Nparticles; ++p)
            {
                // Generate Gaussian random number
                auto [z0, z1] = generateGaussianNoise(0.0, 1.0);

                //
                x[p] += xi[p] * dt+ force*eps[p]*dt;
                // Update Ornstein-Uhlenbeck noise
                xi[p] = xi[p] * (1.0 - dt / tau) + (std::sqrt(2.0 * D * dt) / tau) * z0;
                
                msd_t[i]+= x[p]*x[p]/(Nparticles*runs);
                dis_t[i]+= x[p]/(Nparticles*runs);
            }
        }
    }
    // Write to files
    for (int i = 0; i < Nsteps; ++i)
    {
        msd_file << i * dt << " " << msd_t[i] << " " << tau << std::endl;
        disfile << i*dt << " " << dis_t[i] << " " << tau << std::endl;
    }
}