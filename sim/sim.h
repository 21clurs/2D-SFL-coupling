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
        
    private:
        int n;
        float dt;
        float gravity = -1.0;
        
        Mesh m;

        void step_advect();
        void step_HHD();
        void step_BEM();

        Eigen::Vector2d lin_interp(Eigen::Vector2d v_a, Eigen::Vector2d v_b, double q);

        double cross2d(Eigen::Vector2d a, Eigen::Vector2d b);

};

#endif