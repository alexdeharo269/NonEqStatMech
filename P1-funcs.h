#include <cmath>
#include <limits>
#include <random>
#include <utility>
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

inline std::pair<double, double> diffuse(double x, double y, double Gamma, double dt){
    auto [z0, z1] = generateGaussianNoise(0.0, 1.0);
    x += z0 * sqrt(2.0 * Gamma * dt);
    y += z1 * sqrt(2.0 * Gamma * dt);
    return std::make_pair(x,y);}