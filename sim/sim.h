#ifndef SIM_H
#define SIM_H

#include <math.h>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include "mesh.h"
#include "wallobject.h"

class Sim
{
    friend class TestingHelpers;
    public:

        Sim(Mesh& m, int n, float dt); 
        
        void addWall(WallObject* wall);

        bool outputFrame(std::string filename, std::string filelocation="./out/");

        void step_sim(double curr_t);
        void step_advect(double t);
        void step_HHD();
        void step_gravity();
        void step_BEM();
        
        std::vector<Eigen::Vector2d>& get_vels(){ return m.vels; }
        
        void collide();
    private:
        int n;
        float dt;
        Mesh& m;

        double sigma, sigma_SL, sigma_SA;
        double rho;
        Eigen::Vector2d gravity;

        std::vector<WallObject*> walls;

        void remesh();

        Eigen::Vector2d lin_interp(Eigen::Vector2d v_a, Eigen::Vector2d v_b, double q);
        double M_1(double t);
        double M_2(double t);

        double cross2d(Eigen::Vector2d a, Eigen::Vector2d b);

        void step_BEM_BC(Eigen::VectorXd& BC_p, Eigen::VectorXd& BC_dpdn);
        void step_BEM_solve(Eigen::VectorXd& BC_p, Eigen::VectorXd& BC_dpdn, Eigen::VectorXd& p, Eigen::VectorXd& dpdn);
        void step_BEM_gradP(Eigen::VectorXd& BC_p, Eigen::VectorXd& BC_dpdn, Eigen::VectorXd& p, Eigen::VectorXd& dpdn);

};

#endif