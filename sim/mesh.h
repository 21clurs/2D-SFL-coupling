#ifndef MESH_H
#define MESH_H

#include <Eigen/Dense>
#include <vector>

class Mesh
{
    friend class Sim;
    friend class TestingHelpers;
    protected:
        std::vector<Eigen::Vector2d> verts;
        std::vector<Eigen::Vector2i> faces;     // assuming these are organized in some orientation [prev vertex, next vertex]
        std::vector<Eigen::Vector2d> vels;

        // these are populated once at construction, and updated when the mesh is modified (i.e. edge splits, collapses)
        std::vector<int> vertsPrevFace;
        std::vector<int> vertsNextFace;

    public:
        // constructor
        Mesh(const std::vector<Eigen::Vector2d>& in_verts, const std::vector<Eigen::Vector2i>& in_faces);
        Mesh(const std::vector<Eigen::Vector2d>& in_verts, const std::vector<Eigen::Vector2i>& in_faces,  const std::vector<Eigen::Vector2d>& in_vels);

        void update_face_orientations_from_norms(const std::vector<Eigen::Vector2d>& face_normals);

        const Eigen::Vector2i verts_from_face(const int faceIndex);
        const Eigen::Vector2i faces_from_vert(const int vertIndex);
        const int prev_neighbor_index(const int vertIndex);
        const int next_neighbor_index(const int vertIndex);
        const Eigen::Vector2d prev_neighbor(const int vertIndex);
        const Eigen::Vector2d next_neighbor(const int vertIndex);

        const int other_face_from_vert(const int vertIndex, const int faceIndex);
        const int other_vert_from_face(const int faceIndex, const int vertIndex);

        double face_length(const int faceIndex);
        double calc_total_face_length();
        double calc_avg_face_length();
        double vert_area(const int vertIndex);

        const Eigen::Vector2d calc_vertex_normal(const int vertIndex);
        const Eigen::Vector2d calc_vertex_tangent(const int vertIndex);
        const Eigen::Vector2d calc_face_normal(const int faceIndex);    // outward normal
        const Eigen::Vector2d calc_face_tangent(const int faceIndex);   // currently, working with orders of faces such that this gives clockwise tangent

        double signed_mean_curvature(const int vertIndex);
        double solid_angle(const int vertIndex);

        double windingNumber(const Eigen::Vector2d& x);
    protected:
        void update_neighbor_face_vecs();
        void swap_face_vertices(const int faceIndex);
        double turning_angle(Eigen::Vector2d a, Eigen::Vector2d b);
};

#endif