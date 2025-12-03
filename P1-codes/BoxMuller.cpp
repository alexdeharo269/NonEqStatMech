#include <cmath>
#include <limits>
#include <random>
#include <utility>
#include "P1-funcs.h"
#define M_PI 3.14159265358979323846


int main()
{
    // Example usage
    double mu = 0.0;    // Mean
    double sigma = 1.0; // Standard deviation

    auto [z0, z1] = generateGaussianNoise(mu, sigma);
    // Output the generated Gaussian noise values
    printf("Generated Gaussian Noise: X = %f, Y = %f\n", z0, z1);

    // Now we compute 10e4 numbers and calculate the mean and standard deviation of these quantities in 2D
    const int N=1e4;
    double sumx=0.0, sumy=0.0;
    double sum2x=0.0, sum2y=0.0;
    
    FILE* outfile = fopen("gaussian_stats.txt", "w");
    if (!outfile) {
        printf("Error opening file!\n");
        return 1;
    }
    
    for (int i=0; i<N/2; i++) {
        auto [z0, z1] = generateGaussianNoise(mu, sigma);
        fprintf(outfile, "%f %f\n", z0, z1);
        sumx += z0; sumy += z1;
        sum2x += z0*z0; sum2y += z1*z1;
    }
    
    double meanx = sumx/(N/2);
    double meany = sumy/(N/2);
    double stdx = sqrt(sum2x/(N/2) - meanx*meanx);
    double stdy = sqrt(sum2y/(N/2) - meany*meany);
    
    
    fclose(outfile);
    
    fprintf(stdout, "Computed Mean: X = %f, Y = %f\n", meanx, meany);
    fprintf(stdout, "Computed Std Dev: X = %f, Y = %f\n", stdx, stdy);


    return 0;
}