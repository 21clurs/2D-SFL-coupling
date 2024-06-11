#ifndef RIGID_BODY_H
#define RIGID_BODY_H

#include <math.h>
#include <Eigen/Dense>
#include <iostream>
#include "solidmesh.h"

class RigidBody : public SolidMesh
{
    public:
        double theta;                   // RigidBody rotation around its COM
        Eigen::Vector2d translation;    // RigidBody translation
        
        double area;
        Eigen::Vector2d com;            // center of mass
        double moi;                     // moment of inertia about the COM

        // constructors
        RigidBody(const std::vector<Eigen::Vector2d>& in_verts, const std::vector<Eigen::Vector2i>& in_faces);
        RigidBody(const std::vector<Eigen::Vector2d>& in_verts, const std::vector<Eigen::Vector2i>& in_faces, const std::vector<Eigen::Vector2d>& in_vels);

    protected:
        void calculateArea();
        void calculateCOM();
        void calculateMOI();
};

#endif