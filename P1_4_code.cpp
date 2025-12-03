#include <cmath>
#include <limits>
#include <random>
#include <utility>
#include "P1-funcs.h"
#define M_PI 3.14159265358979323846

const double D = 0.1;
const double gamma = 1.0;
const double dt = 0.001;
const int L=100;
const int Nmols=1000;
const int tau=100;

// Now we compute the Euler-Maruyama algorithm for the brownian motion of a particle in 2D
// The twist is that now we use periodic boudary conditions LxL. 
// For each particle, we want to keep track of its position and the number of times it has crossed
// the boundary in each direction.

int main()
{
    //Initilize positions randomly from a uniform distribution
    std::vector<std::tuple<double, double, int, int>> pos;
    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<> runif(0.0, L);
    for (int i=0; i<Nmols; i++){
        double x=runif(rng);
        double y=runif(rng);
        pos.push_back(std::make_tuple(x,y,0,0));
    }

    FILE *outfile = fopen("P1_4-data.txt", "w");
    if (!outfile)
    {
        printf("Error opening file!\n");
        return 1;
    }
    
    //For tau steps, we update the position of each particle and their crossings of the boundary bc of diffusion

    for (int i = 0; i < tau; i++)
    {
        
        //Print periodic positons of all particles x, y 
        printpos(pos, outfile);


        //Update periodic positions of all particles, but keeping track of the number of crossings
        for (int j=0; j<Nmols; j++){
            double x, y;
            int nx, ny;
            auto [x, y , nx, ny]=pos[j];
            
            std::tuple<double, double> upt=diffuse(x,y,D,dt);
            auto [xupt, yupt]=upt;
            //Check if the particle has crossed the boundary 
            if (xupt <0){xupt+=L;nx-=1;}
            else if (xupt >=L){xupt-=L;nx+=1;}
            
            if (yupt <0){yupt+=L;ny-=1;}
            else if (yupt >=L){yupt-=L;ny+=1;}

            pos[j]=std::make_tuple(xupt,yupt,nx,ny);
            
        }
        printpos(pos, outfile);
    }

    fprintf(stdout, "Done!");
    fclose(outfile);
    return 0;
}

void printpos(const std::vector<std::tuple<double, double, int, int>>& pos, FILE* outfile){
    for (const auto& [x, y, nx, ny] : pos) {
        fprintf(outfile, "%f %f ", x, y);
    }
    fprintf(outfile, "\n");
}