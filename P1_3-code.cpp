#include <cmath>
#include <limits>
#include <random>
#include <utility>
#include "P1-funcs.h"
#define M_PI 3.14159265358979323846


const double D=0.1;
const double gamma=1.0;
const double dt=0.001;

//Now we compute the Euler-Maruyama algorithm for the brownian motion of a particle in 2D

int main(){

FILE *outfile = fopen("P1_3-data.txt", "w");
if (!outfile)
{
    printf("Error opening file!\n");
    return 1;
}

for (int k=0; k<4;k++){
    double x = 0.0, y = 0.0;
    fprintf(outfile, "%f %f\n", x, y);

    for (int i = 0; i < 1e4; i++)
    {
        fprintf(outfile, "%f %f\n", x, y);
        auto [xupt, yupt] = diffuse(x, y, D, dt);
        x = xupt;
        y = yupt;
    }
fprintf(outfile, "\n");
}
fprintf(stdout, "Done!");
fclose(outfile);
return 0;
}