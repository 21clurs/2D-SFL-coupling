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

double Mesh::face_length(const int faceIndex){
    Eigen::Vector2i endpts = verts_from_face(faceIndex);
    return (verts[endpts[0]]-verts[endpts[1]]).norm();
}

double Mesh::vert_area(const int vertIndex){
    return (1.0/2.0)*(face_length(faces_from_vert(vertIndex)[0])+face_length(faces_from_vert(vertIndex)[1]));
}
const Eigen::Vector2d Mesh::calc_vertex_normal(const int vertIndex){
    Eigen::Vector2i adjacent_face_inds = faces_from_vert(vertIndex);

    Eigen::Vector2d n_a = calc_face_normal(adjacent_face_inds[0]);
    Eigen::Vector2d n_b = calc_face_normal(adjacent_face_inds[1]);

    double d_a = face_length(adjacent_face_inds[0]);
    double d_b = face_length(adjacent_face_inds[1]);

    Eigen::Vector2d n = (d_a*n_a + d_b*n_b)/(d_a+d_b);
    n.normalize();
    return n;
}
const Eigen::Vector2d Mesh::calc_face_normal(const int faceIndex){
    Eigen::Vector2i endpts = verts_from_face(faceIndex);
    Eigen::Vector2d x_s = verts[endpts[0]];
    Eigen::Vector2d x_t = verts[endpts[1]];

    // 90-deg clockwise rotation matrix
    Eigen::Matrix2d rotate;
    rotate << 0, 1, -1, 0;
    Eigen::Vector2d n = rotate * (x_t-x_s);
    n.normalize();
    return n;
}
const Eigen::Vector2d Mesh::calc_face_tangent(const int faceIndex){
    // clockwise tangent?
    Eigen::Vector2i endpts = verts_from_face(faceIndex);
    Eigen::Vector2d x_s = verts[endpts[0]];
    Eigen::Vector2d x_t = verts[endpts[1]];

    Eigen::Vector2d t = -(x_t-x_s);
    t.normalize();
    return t;
}
double Mesh::solid_angle(const int vertIndex){
    Eigen::Vector2d t_prev = verts[vertIndex] - prev_neighbor(vertIndex);
    Eigen::Vector2d t_next = next_neighbor(vertIndex) - verts[vertIndex];

    double theta = M_PI - turning_angle(t_prev,t_next);
    return theta/(2*M_PI);
}
double Mesh::signed_mean_curvature(const int vertIndex){
    Eigen::Vector2d t_prev = verts[vertIndex] - prev_neighbor(vertIndex);
    Eigen::Vector2d t_next = next_neighbor(vertIndex) - verts[vertIndex];

    double phi = turning_angle(t_prev, t_next);
    return 2*phi/(t_prev.norm() + t_next.norm());
}
double Mesh::turning_angle(Eigen::Vector2d a, Eigen::Vector2d b){
    a.normalize();
    b.normalize();
    // FIX THIS LATER
    double crossprod = a[0]*b[1] - a[1]*b[0];
    double dotprod = a[0]*b[0] + a[1]*b[1];

    return atan2(crossprod,dotprod);
}