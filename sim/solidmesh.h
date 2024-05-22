#ifndef SOLID_MESH_H
#define SOLID_MESH_H

#include <math.h>
#include <Eigen/Dense>
#include <iostream>
#include "mesh.h"
#include "wallobject.h"


class SolidMesh : public Mesh
{
    public:
        // constructors
        SolidMesh(const std::vector<Eigen::Vector2d>& in_verts, const std::vector<Eigen::Vector2i>& in_faces) : Mesh(in_verts, in_faces){}
        void setEps(double e){ epsilon = e; }
        bool checkCollisionAndSnap(Eigen::Vector2d& x);

    private:
        std::vector<Eigen::Vector2d> vert_normals;
        std::vector<Eigen::Vector2d> face_normals;

        double epsilon = 0.01;
};
#endif