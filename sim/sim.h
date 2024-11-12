#ifndef SIM_H
#define SIM_H

#include <math.h>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include "liquidmesh.h"
#include "rigidbody.h"

class Sim
{
    friend class TestingHelpers;
    friend class Scenes;
    public:

        Sim();
        Sim(LiquidMesh& m, int n, float dt); 
        ~Sim(); // need to manually delete stuff in the solids vector...

        void addRigidBody(RigidBody* rigidBody);

        bool outputOscillationX(std::vector<float> oscillation_t, std::vector<float> oscillation_x, std::string filename, std::string filelocation="./out/oscillation_tests/");
        bool outputPostStepV(std::string filename, std::string filelocation="./out/buoyancy_tests/");
        bool outputFrame(std::string filename, std::string filelocation="./out/");

        void step_sim(double curr_t);
        //void step_solidinflux();
        void step_advect(double t);
        //void step_solidinfluxreverse();
        void step_HHD();
        void step_gravity();
        void step_BEM();
        
        std::vector<Eigen::Vector2d>& get_vels(){ return m.vels; }
        
        void collide();

        void genMarkerParticles(double l, double r, double b, double t, double spacing);
        Eigen::Vector2d HHD_interior(Eigen::Vector2d x, double delta); 
        std::vector<Eigen::Vector2d> markerparticles; // TODO: INITIATE THESE

        static bool setAndLoadSimOptions(std::string infileName);
        void run();

        static double cross2d(Eigen::Vector2d a, Eigen::Vector2d b);
    private:
        int n;
        float dt;
        LiquidMesh m;

        int outframe_frequency;

        // simulation parameters
        double sigma, sigma_SL, sigma_SA;
        double rho;
        Eigen::Vector2d gravity;

        // solid objects in the sim
        //std::vector<SolidMesh*> solids;
        // rigid bodies in the sim
        //std::vector<RigidBody*> rigidBodies;

        std::vector<RigidBody*> rigidBodies_scripted;
        std::vector<RigidBody*> rigidBodies_unscripted;

        void remesh(int remesh_itr);

        Eigen::Vector2d lin_interp(Eigen::Vector2d v_a, Eigen::Vector2d v_b, double q);
        double M_1(double t);
        double M_2(double t);

        void step_BEM_BC(Eigen::VectorXd& BC_p, Eigen::VectorXd& BC_dpdn);
        void step_BEM_solve(Eigen::VectorXd& BC_p, Eigen::VectorXd& BC_dpdn, Eigen::VectorXd& p, Eigen::VectorXd& dpdn, std::vector<Eigen::Vector3d>& V_rigidBodies);
        void step_BEM_gradP(Eigen::VectorXd& BC_p, Eigen::VectorXd& BC_dpdn, Eigen::VectorXd& p, Eigen::VectorXd& dpdn);
        void step_BEM_rigidBodyV(std::vector<Eigen::Vector3d> & V_rigidBodies);
};

#endif