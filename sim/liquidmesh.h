#ifndef LIQUID_MESH_H
#define LIQUID_MESH_H

#include <Eigen/Dense>
#include <vector>
#include "mesh.h"
#include "wallobject.h"
#include "solidmesh.h"

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

        std::vector<bool> get_solid_faces();

        void remesh();

        void reset_boundary_types();
        void collide_with_wall(WallObject& w);
        void collide_with_solid(SolidMesh& s);

        void laplacian_smoothing();
        void edge_collapse();
        void edge_split();
    private:
        std::vector<bool> is_air;
        std::vector<bool> is_solid;
        std::vector<bool> is_triple;

        void update_triple_points();

        double minFaceLength;
        double maxFaceLength;

};

#endif