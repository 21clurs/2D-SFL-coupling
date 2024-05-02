#ifndef SIM_H
#define SIM_H

#include <math.h>
#include <Eigen/Dense>
#include "mesh.h"

class Sim
{
    public:

        Sim(Mesh& m, int n, float dt); 

        bool outputFrame(std::string filename, std::string filelocation="./out/");

        void step_sim(int frame);

        void step_advect();
        
    private:
        int n;
        float dt;
        float gravity = -1;
        
        Mesh m;


};

#endif