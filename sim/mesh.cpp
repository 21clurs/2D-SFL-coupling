#include "mesh.h"
#include <iostream>
Mesh::Mesh(const std::vector<Eigen::Vector2d>& in_verts, const std::vector<Eigen::Vector2i>& in_faces):
    verts(in_verts),
    faces(in_faces)
{
    vels = std::vector<Eigen::Vector2d>(verts.size(), Eigen::Vector2d(0.0, 0.0));

    vertsPrevFace = std::vector<int>(verts.size(),0);
    vertsNextFace = std::vector<int>(verts.size(),0);
    update_neighbor_face_vecs();
}

Mesh::Mesh(const std::vector<Eigen::Vector2d>& in_verts, const std::vector<Eigen::Vector2i>& in_faces, const  std::vector<Eigen::Vector2d>& in_vels):
    verts(in_verts),
    faces(in_faces),
    vels(in_vels)
{
    vertsPrevFace = std::vector<int>(verts.size(),0);
    vertsNextFace = std::vector<int>(verts.size(),0);
    update_neighbor_face_vecs();
}   

// returns the indices of the vertices at the endpts of the given face
const Eigen::Vector2i Mesh::verts_from_face(const int faceIndex)
{
    return faces[faceIndex];
}
const Eigen::Vector2i Mesh::faces_from_vert(const int vertIndex)
{
    return Eigen::Vector2i(vertsPrevFace[vertIndex], vertsNextFace[vertIndex]);
}
const Eigen::Vector2d Mesh::prev_neighbor(const int vertIndex)
{
    return verts[faces[vertsPrevFace[vertIndex]][0]];
}
const Eigen::Vector2d Mesh::next_neighbor(const int vertIndex)
{
    return verts[faces[vertsNextFace[vertIndex]][1]];
}

void Mesh::update_neighbor_face_vecs()
{
    int s_index,t_index;
    for (uint i=0; i<faces.size(); i++){
        s_index = faces[i][0];
        t_index = faces[i][1];

        vertsPrevFace[t_index] = i;
        vertsNextFace[s_index] = i;
    }
}