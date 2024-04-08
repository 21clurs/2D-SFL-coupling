#ifndef MESH_H
#define MESH_H

#include <Eigen/Dense>
#include <vector>

class Mesh
{
    public:
        std::vector<Eigen::Vector2d> verts;
        std::vector<Eigen::Vector2i> faces;
        std::vector<Eigen::Vector2d> vels;

        Mesh(const std::vector<Eigen::Vector2d>& in_verts, const std::vector<Eigen::Vector2i>& in_faces);
        Mesh(const std::vector<Eigen::Vector2d>& in_verts, const std::vector<Eigen::Vector2i>& in_faces, const std::vector<Eigen::Vector2d>& in_vels);

        const Eigen::Vector2i verts_from_face(const int faceIndex);
        const Eigen::Vector2i faces_from_vert(const int vertIndex);
        const Eigen::Vector2d prev_neighbor(const int vertIndex);
        const Eigen::Vector2d next_neighbor(const int vertIndex);

    private:
        // these are populated once at construction
        std::vector<int> vertsPrevFace;
        std::vector<int> vertsNextFace;

        void update_neighbor_face_vecs();
};

#endif