#include "mesh.h"

Mesh::Mesh(){
    verts = std::vector<Eigen::Vector2d>(0, Eigen::Vector2d(0.0, 0.0));
    faces = std::vector<Eigen::Vector2i>(0, Eigen::Vector2i(0, 0));

    vels = std::vector<Eigen::Vector2d>(verts.size(), Eigen::Vector2d(0.0, 0.0));
    
    vertsPrevFace = std::vector<int>(verts.size(),0);
    vertsNextFace = std::vector<int>(verts.size(),0);
}

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

const int Mesh::other_face_from_vert(const int vertIndex, const int faceIndex){
    Eigen::Vector2i adjacent_faces = faces_from_vert(vertIndex);
    assert(faceIndex == adjacent_faces.x() || faceIndex == adjacent_faces.y());
    int otherFaceIndex = faceIndex == adjacent_faces.x() ? adjacent_faces.y() : adjacent_faces.x();
    return otherFaceIndex;
}

const int Mesh::other_vert_from_face(const int faceIndex, const int vertIndex){
    Eigen::Vector2i adjacent_verts = verts_from_face(faceIndex);
    assert(vertIndex == adjacent_verts.x() || vertIndex == adjacent_verts.y());
    int otherVertIndex = vertIndex == adjacent_verts.x() ? adjacent_verts.y() : adjacent_verts.x();
    return otherVertIndex;
}

double Mesh::face_length(const int faceIndex){
    Eigen::Vector2i endpts = verts_from_face(faceIndex);
    return (verts[endpts[0]]-verts[endpts[1]]).norm();
}
double Mesh::calc_total_face_length(){
    double total_face_length = 0;
    for (size_t i=0; i<faces.size(); i++)
        total_face_length += face_length(i);
    return total_face_length;
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


double Mesh::calc_area(){
    double A = 0;
    Eigen::Vector2d v_i, v_i_next;
    for(size_t i=0; i<verts.size(); i++){
        v_i = verts[i];
        v_i_next = next_neighbor(i);
        A += (v_i.x() - v_i_next.y()) - (v_i_next.x()-v_i.y());
    }
    return 0.5*A;
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
bool Mesh::solid_angle_is_acute(const int vertIndex){
    return (abs(solid_angle(vertIndex)) <= (M_PI/2)+0.1 || abs(solid_angle(vertIndex)) >= 2*M_PI - ((M_PI/2)+0.1));
}
double Mesh::signed_mean_curvature(const int vertIndex){
    Eigen::Vector2d t_prev = verts[vertIndex] - prev_neighbor(vertIndex);
    Eigen::Vector2d t_next = next_neighbor(vertIndex) - verts[vertIndex];

    double phi = turning_angle(t_prev, t_next);
    return 2*phi/(t_prev.norm() + t_next.norm());
}
void Mesh::update_neighbor_face_vecs(){
    int s_index,t_index;
    if (faces.size() != vertsPrevFace.size() || faces.size() != vertsNextFace.size()){
        vertsPrevFace.resize(faces.size());
        vertsNextFace.resize(faces.size());
    }
    for (size_t i=0; i<faces.size(); i++){
        s_index = faces[i][0];
        t_index = faces[i][1];

        vertsPrevFace[t_index] = i;
        vertsNextFace[s_index] = i;
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

double Mesh::windingNumber(const Eigen::Vector2d& p){
    double w = 0;
    for (size_t i=0; i<faces.size(); i++){
        Eigen::Vector2d a = verts[faces[i][0]] - p;
        Eigen::Vector2d b = verts[faces[i][1]] - p;
        w += atan2( (a.x()*b.y()-a.y()*b.x()) , (a.dot(b)) );
    }
    return w/(2*M_PI);
}