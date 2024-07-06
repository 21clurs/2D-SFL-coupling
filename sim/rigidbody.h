#ifndef RIGID_BODY_H
#define RIGID_BODY_H

#include <math.h>
#include <Eigen/Dense>
#include <iostream>
#include "solidmesh.h"

class RigidBody : public SolidMesh
{
    public:
        double area;
        Eigen::Vector2d com;            // center of mass
        double moi;                     // moment of inertia about the COM

        // constructors
        RigidBody();
        RigidBody(const std::vector<Eigen::Vector2d>& in_verts, const std::vector<Eigen::Vector2i>& in_faces);
        RigidBody(const std::vector<Eigen::Vector2d>& in_verts, const std::vector<Eigen::Vector2i>& in_faces, const std::vector<Eigen::Vector2d>& in_vels);

        // setters
        void setRotation(double theta);
        void setTranslation(Eigen::Vector2d t){ translation = t; };

        void retrieveCurrentVerts(std::vector<Eigen::Vector2d>& v);

        void updateRigidBodyVars();     // something to quickly call to update all of Area, COM, MOI, rotation matrix

    protected:
        Eigen::Vector2d translation;    // RigidBody translation
        double rotationTheta;           // RigidBody angle of rotation around its COM

        Eigen::Matrix2d rotationMat;    // this is just updated whenever theta is changed...

        void calculateArea();
        void calculateCOM();
        void calculateMOI();

        void updateRotationMat();
};

#endif