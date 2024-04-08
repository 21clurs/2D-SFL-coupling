#ifndef GENERATE_SHAPES
#define GENERATE_SHAPES

#include <math.h>
#include <vector>
#include <Eigen/Dense>

class GenerateShape
{
    public:
        static void circle(int n, std::vector<Eigen::Vector2d>& verts, std::vector<Eigen::Vector2i>& faces);
        static void square(int n, std::vector<Eigen::Vector2d>& verts, std::vector<Eigen::Vector2i>& faces);
        static void ellipse(int n, std::vector<Eigen::Vector2d>& verts, std::vector<Eigen::Vector2i>& faces);
        static void oscillation(int n, std::vector<Eigen::Vector2d>& verts, std::vector<Eigen::Vector2i>& faces);
        static void semicircle_v(int n, std::vector<Eigen::Vector2d>& verts, std::vector<Eigen::Vector2i>& faces);
        static void semicircle_h(int n, std::vector<Eigen::Vector2d>& verts, std::vector<Eigen::Vector2i>& faces);
};
#endif