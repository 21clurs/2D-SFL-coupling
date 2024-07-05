#ifndef SOLID_MESH_H
#define SOLID_MESH_H

#include <math.h>
#include <Eigen/Dense>
#include <iostream>
#include "mesh.h"
#include "liquidmesh.h"

class SolidMesh : public Mesh
{
    friend class LiquidMesh;
    public:
        static bool loadMeshFromFile(SolidMesh &m, std::string infileName);
    public:
        // constructors
        SolidMesh();
        SolidMesh(const std::vector<Eigen::Vector2d>& in_verts, const std::vector<Eigen::Vector2i>& in_faces);
        SolidMesh(const std::vector<Eigen::Vector2d>& in_verts, const std::vector<Eigen::Vector2i>& in_faces, const std::vector<Eigen::Vector2d>& in_vels);
        void setEps(double e){ epsilon = e; }
        bool checkCollisionAndSnap(Eigen::Vector2d& currpt);
        void setVelFunc(std::function<Eigen::Vector2d(double)> func);
        void advectFE(double curr_t, double dt);
        void collideAndSnap(LiquidMesh& l);
    protected:
        std::vector<Eigen::Vector2d> vert_normals;
        std::vector<Eigen::Vector2d> face_normals;
        Eigen::Vector2d v_effective; // assuming only translation-type movement

        double epsilon = 0.02;

        std::function<Eigen::Vector2d(double)> vel_func;
};
#endif