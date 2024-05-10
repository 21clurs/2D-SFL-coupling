#include "generatefields.h"

void GenerateField::zero(int n, std::vector<Eigen::Vector2d>& verts, std::vector<Eigen::Vector2d>& vels){
    for (size_t i=0; i<n; i++){
        vels[i].x() = 0;
        vels[i].y() = 0;
    }
}
void GenerateField::constant(int n, std::vector<Eigen::Vector2d>& verts, std::vector<Eigen::Vector2d>& vels){
    for (size_t i=0; i<n; i++){
        vels[i].x() = 1;
        vels[i].y() = 0;
    }
}
void GenerateField::linear_x(int n, std::vector<Eigen::Vector2d>& verts, std::vector<Eigen::Vector2d>& vels){
    for (size_t i=0; i<n; i++){
        vels[i].x() = verts[i].x();
        vels[i].y() = 0;
    }
}
void GenerateField::linear_y(int n, std::vector<Eigen::Vector2d>& verts, std::vector<Eigen::Vector2d>& vels){
    for (size_t i=0; i<n; i++){
        vels[i].x() = 0;
        vels[i].y() = verts[i].y();
    }
}
void GenerateField::harmonic_1(int n, std::vector<Eigen::Vector2d>& verts, std::vector<Eigen::Vector2d>& vels){
    for (size_t i=0; i<n; i++){
        vels[i].x() = verts[i].y();
        vels[i].y() = verts[i].x();
    }
}
void GenerateField::harmonic_2(int n, std::vector<Eigen::Vector2d>& verts, std::vector<Eigen::Vector2d>& vels){
    for (size_t i=0; i<n; i++){
        vels[i].x() = verts[i].x();
        vels[i].y() = -verts[i].y();
    }
}
void GenerateField::harmonic_3(int n, std::vector<Eigen::Vector2d>& verts, std::vector<Eigen::Vector2d>& vels){
    for (size_t i=0; i<n; i++){
        vels[i].x() = -verts[i].x()*verts[i].x() + verts[i].y()*verts[i].y();
        vels[i].y() = 2* verts[i].x() * verts[i].y();
    }
}