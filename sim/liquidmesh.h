#ifndef LIQUID_MESH_H
#define LIQUID_MESH_H

#include <Eigen/Dense>
#include <vector>
#include "mesh.h"

class LiquidMesh : public Mesh
{
    friend class Sim;
    friend class Scenes;
    friend class SolidMesh;
    friend class RigidBody;
    friend class TestingHelpers;
    protected:
        std::vector<Eigen::Vector2d> vels_solid;
    public:
        // constructors
        LiquidMesh();
        LiquidMesh(const std::vector<Eigen::Vector2d>& in_verts, const std::vector<Eigen::Vector2i>& in_faces);
        LiquidMesh(const std::vector<Eigen::Vector2d>& in_verts, const std::vector<Eigen::Vector2i>& in_faces, const std::vector<Eigen::Vector2d>& in_vels);
        
        void set_boundaries(std::vector<bool> air, std::vector<bool> solid, std::vector<bool> triple);
        void set_boundaries(std::vector<bool> air, std::vector<bool> solid, std::vector<bool> triple, std::vector<bool> corner);
        void set_boundaries_for_vertex(int i, bool air, bool solid, bool triple, bool corner);

        std::vector<bool> get_solid_faces();
        const Eigen::Vector2d calc_vertex_solid_normal(const int vertIndex);

        void remesh();

        void reset_boundary_types();
        void reset_face_length_limits();

        void laplacian_smoothing();
        void edge_collapse();
        void edge_split();

        double signed_min_dist(Eigen::Vector2d x);

        void resize_mesh(size_t n);

        static bool loadMeshFromFile(LiquidMesh &l, std::string infileName);
    protected:
        std::vector<bool> is_air;
        std::vector<bool> is_solid;
        std::vector<bool> is_triple;

        std::vector<bool> is_corner;

        std::vector<uint> per_vertex_rb_contact;

        void update_triple_points();

        double minFaceLength;
        double maxFaceLength;

};

#endif