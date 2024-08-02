#ifndef RIGID_BODY_H
#define RIGID_BODY_H

#include <math.h>
#include <Eigen/Dense>
#include <iostream>
#include "mesh.h"
#include "liquidmesh.h"

class RigidBody : public Mesh
{
    friend class Sim;
    friend class Scenes;
    friend class LiquidMesh;
    public:
        double area;
        double mass;                    // rho * area
        Eigen::Vector2d com;            // center of mass
        double moi;                     // moment of inertia about the COM

    public:
        // constructors
        RigidBody();
        RigidBody(const std::vector<Eigen::Vector2d>& in_verts, const std::vector<Eigen::Vector2i>& in_faces);
        RigidBody(const std::vector<Eigen::Vector2d>& in_verts, const std::vector<Eigen::Vector2i>& in_faces, const std::vector<Eigen::Vector2d>& in_vels);

        // setters
        void setEps(double e){ epsilon = e; }

        void setRotation(double theta){ rotationTheta = theta; };
        void setTranslation(Eigen::Vector2d t){ translation = t; };

        void setRigidBodyV(Eigen::Vector3d V_in) { V_t = Eigen::Vector2d(V_in.x(), V_in.y()); V_omega = V_in.z(); };
        void setRigidBodyV(Eigen::Vector2d V_t_in, double V_omega_in) { V_t = V_t_in; V_omega = V_omega_in; };
        void setRigidBodyV(double V_t_x, double V_t_y, double V_omega_in) { V_t = Eigen::Vector2d(V_t_x, V_t_y); V_omega = V_omega_in; };

        void setVelFunc(std::function<Eigen::Vector2d(double)> func);

        // getter
        Eigen::Vector3d retrieveRigidBodyV() { return Eigen::Vector3d(V_t.x(), V_t.y(), V_omega); };
        
        void advectFE(double dt);
        void advectFE(double curr_t, double dt);

        void updateVerts();                 // updates the world coordinates (stored in verts), based on the translation & rotation
        void updatePerVertexVels();                  // updates the per-vertex velocities, based on the translational & rotational velocities
        void updateRigidBodyVars();         // something to quickly call to update all of Area, COM, MOI, rotation matrix
    
        void collideAndSnap(LiquidMesh& l);
    public:
        static bool loadMeshFromFile(RigidBody &m, std::string infileName);

    protected:
        u_int rb_index_in_sim;

        double rho;
        // the coordinates of the RigidBody without any of the rotation/translation
        // the vertices of the RigidBody with its transformations applied are stored in verts (field inherited from Mesh)
        std::vector<Eigen::Vector2d> verts_no_transform; 

        Eigen::Vector2d translation;        // RigidBody translation
        double rotationTheta;               // RigidBody angle of rotation around its COM
        
        Eigen::Vector2d V_t;                // RigidBody translational velocity
        double V_omega;                     // RigidBody rotational velocity

        double epsilon = 0.02;
        std::function<Eigen::Vector2d(double)> vel_func;

        void calculateArea();
        void calculateCOM();
        void calculateMOI();
};

#endif