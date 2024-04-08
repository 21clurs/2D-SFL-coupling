#ifndef GENERATE_SHAPES
#define GENERATE_SHAPES

#include <math.h>
#include <vector>
#include <Eigen/Dense>

void generateCircle(int n, std::vector<Eigen::Vector2d>& verts, std::vector<Eigen::Vector2i>& faces);
void generateSquare(int n, std::vector<Eigen::Vector2d>& verts, std::vector<Eigen::Vector2i>& faces);
void generateEllipse(int n, std::vector<Eigen::Vector2d>& verts, std::vector<Eigen::Vector2i>& faces);

#endif