#ifndef LIQUID_MESH_H
#define LIQUID_MESH_H

#include <Eigen/Dense>
#include <vector>
#include "mesh.h"
#include "wallobject.h"
#include "rectmesh.h"

class LiquidMesh : public Mesh
{
    friend class Sim;
    friend class TestingHelpers;
    protected:
        std::vector<Eigen::Vector2d> vels_solid;
    public:
        // constructors
        LiquidMesh(const std::vector<Eigen::Vector2d>& in_verts, const std::vector<Eigen::Vector2i>& in_faces);
        LiquidMesh(const std::vector<Eigen::Vector2d>& in_verts, const std::vector<Eigen::Vector2i>& in_faces, const std::vector<Eigen::Vector2d>& in_vels);
        
        void set_boundaries(std::vector<bool> air, std::vector<bool> solid, std::vector<bool> triple);
        void update_face_orientations_from_norms(const std::vector<Eigen::Vector2d>& face_normals);

        std::vector<bool> get_solid_faces();

        void laplacian_smoothing();

        void collide_with_wall(WallObject& w);
        void collide_with_rect(RectMesh r);

    private:
        std::vector<bool> is_air;
        std::vector<bool> is_solid;
        std::vector<bool> is_triple;

        void update_triple_points();
};

#endif