#include "mesh.h"
#include <iostream>
Mesh::Mesh(const std::vector<Eigen::Vector2d>& in_verts, const std::vector<Eigen::Vector2i>& in_faces):
    verts(in_verts),
    faces(in_faces)
{
    vels = std::vector<Eigen::Vector2d>(verts.size(), Eigen::Vector2d(0.0, 0.0));
    vels_solid = std::vector<Eigen::Vector2d>(verts.size(), Eigen::Vector2d(0.0, 0.0));

    vertsPrevFace = std::vector<int>(verts.size(),0);
    vertsNextFace = std::vector<int>(verts.size(),0);
    update_neighbor_face_vecs();

    // defaults to purely liquid-air boundary?
    is_air = std::vector<bool>(verts.size(), true);
    is_solid = std::vector<bool>(verts.size(), false);
    is_triple = std::vector<bool>(verts.size(), false);
}

Mesh::Mesh(const std::vector<Eigen::Vector2d>& in_verts, const std::vector<Eigen::Vector2i>& in_faces, const  std::vector<Eigen::Vector2d>& in_vels):
    verts(in_verts),
    faces(in_faces),
    vels(in_vels)
{
    vels_solid = std::vector<Eigen::Vector2d>(verts.size(), Eigen::Vector2d(0.0, 0.0));

    vertsPrevFace = std::vector<int>(verts.size(),0);
    vertsNextFace = std::vector<int>(verts.size(),0);
    update_neighbor_face_vecs();

    // defaults to purely liquid-air boundary?
    is_air = std::vector<bool>(verts.size(), true);
    is_solid = std::vector<bool>(verts.size(), false);
    is_triple = std::vector<bool>(verts.size(), false);
}   

void Mesh::set_boundaries(std::vector<bool> air, std::vector<bool> solid, std::vector<bool> triple){
    is_air = air;
    is_solid = solid;
    is_triple = triple;
}

void Mesh::update_face_orientations_from_norms(const std::vector<Eigen::Vector2d>& face_normals){
    assert(face_normals.size() == faces.size());
    for (size_t i=0; i<faces.size(); i++){
        if (calc_face_normal(i).dot(face_normals[i]) > 0)
            swap_face_vertices(i);
        assert(calc_face_normal(i).dot(face_normals[i]) <= 0);
    }

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
const int Mesh::prev_neighbor_index(const int vertIndex)
{
    return faces[vertsPrevFace[vertIndex]][0];
}
const int Mesh::next_neighbor_index(const int vertIndex)
{
    return faces[vertsNextFace[vertIndex]][1];
}
const Eigen::Vector2d Mesh::prev_neighbor(const int vertIndex)
{
    return verts[faces[vertsPrevFace[vertIndex]][0]];
}
const Eigen::Vector2d Mesh::next_neighbor(const int vertIndex)
{
    return verts[faces[vertsNextFace[vertIndex]][1]];
}

void Mesh::laplacian_smoothing()
{
    std::vector<Eigen::Vector2d> C = std::vector<Eigen::Vector2d>(verts.size(), Eigen::Vector2d(0.0, 0.0));
    std::vector<Eigen::Vector2d> proj = std::vector<Eigen::Vector2d>(verts.size(), Eigen::Vector2d(0.0, 0.0));
    std::vector<Eigen::Vector2d> verts_smoothed = std::vector<Eigen::Vector2d>(verts.size(), Eigen::Vector2d(0.0, 0.0));

    Eigen::Vector2d n_i;
    for (size_t i=0; i<verts.size(); i++){
        if (is_triple[i]){
            C[i] = verts[i];
            proj[i] = Eigen::Vector2d(0.0, 0.0);
        }
        else{
            n_i = calc_vertex_normal(i);
            C[i] = 0.5*(prev_neighbor(i) + next_neighbor(i));
            proj[i] = (n_i*n_i.transpose())*(verts[i]-C[i]);
        }
        verts_smoothed[i] = C[i] + proj[i];
    }
    verts = verts_smoothed;
}

std::vector<bool> Mesh::get_solid_faces(){
    std::vector<bool> solid_faces = std::vector<bool>(faces.size(), false);
    Eigen::Vector2i endpts;
    for (size_t i=0; i<faces.size(); i++){
        endpts = verts_from_face(i);
        solid_faces[i] = (is_solid[endpts[0]] || is_solid[endpts[1]]);
    }
    return solid_faces;
}

double Mesh::face_length(const int faceIndex){
    Eigen::Vector2i endpts = verts_from_face(faceIndex);
    return (verts[endpts[0]]-verts[endpts[1]]).norm();
}

double Mesh::calc_avg_face_length(){
    double total_face_length = 0;
    for (size_t i=0; i<faces.size(); i++)
        total_face_length += face_length(i);
    return total_face_length/faces.size();
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
const Eigen::Vector2d Mesh::calc_vertex_tangent(const int vertIndex){
    
    Eigen::Vector2i adjacent_face_inds = faces_from_vert(vertIndex);

    Eigen::Vector2d t_a = calc_face_tangent(adjacent_face_inds[0]);
    Eigen::Vector2d t_b = calc_face_tangent(adjacent_face_inds[1]);
    
    double d_a = face_length(adjacent_face_inds[0]);
    double d_b = face_length(adjacent_face_inds[1]);

    Eigen::Vector2d t = (d_a*t_a + d_b*t_b)/(d_a+d_b);
    t.normalize();
    return t;
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
    return theta; ///(2*M_PI);
}
double Mesh::signed_mean_curvature(const int vertIndex){
    Eigen::Vector2d t_prev = verts[vertIndex] - prev_neighbor(vertIndex);
    Eigen::Vector2d t_next = next_neighbor(vertIndex) - verts[vertIndex];

    double phi = turning_angle(t_prev, t_next);
    return 2*phi/(t_prev.norm() + t_next.norm());
}
void Mesh::update_neighbor_face_vecs(){
    int s_index,t_index;
    for (size_t i=0; i<faces.size(); i++){
        s_index = faces[i][0];
        t_index = faces[i][1];

        vertsPrevFace[t_index] = i;
        vertsNextFace[s_index] = i;
    }
}
void Mesh::update_triple_points(){
    for(size_t i=0; i<verts.size(); i++){
        if (is_solid[i] && (is_air[next_neighbor_index(i)] || is_air[prev_neighbor_index(i)])){
            is_triple[i] = true;
            is_solid[i] = false;
        }
        assert((is_solid[i] ^ is_air[i] ^ is_triple[i]) && (is_solid[i] + is_air[i] + is_triple[i]==1));
    }
}
void Mesh::swap_face_vertices(const int faceIndex){
    int tmp = faces[faceIndex][0];
    faces[faceIndex][0] = faces[faceIndex][1];
    faces[faceIndex][1] = tmp;
}
double Mesh::turning_angle(Eigen::Vector2d a, Eigen::Vector2d b){
    a.normalize();
    b.normalize();
    // FIX THIS LATER
    double crossprod = a[0]*b[1] - a[1]*b[0];
    double dotprod = a[0]*b[0] + a[1]*b[1];

    return atan2(crossprod,dotprod);
}
void Mesh::collide_with_wall(WallObject w){
    for(size_t i=0; i<verts.size(); i++){
        //std::cout<<"og: "<<verts[i]<<std::endl;
        if(w.checkCollisionAndSnap(verts[i]) == true){

            //std::cout<<"snapped: "<<verts[i]<<std::endl;
            vels_solid[i] = w.calcEffectiveVel();
            
            is_solid[i] = true; 
            is_air[i] = false;
            is_triple[i] = false;
            //std::cout<<"bLAH "<<i<<std::endl;
        }
    }
    update_triple_points();
}