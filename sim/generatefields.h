#ifndef GENERATE_FIELDS
#define GENERATE_FIELDS

#include <math.h>
#include <vector>
#include <Eigen/Dense>

class GenerateField
{
    public:
        static void zero(int n, std::vector<Eigen::Vector2d>& verts, std::vector<Eigen::Vector2d>& vels);           // (0,0)
        static void constant(int n, std::vector<Eigen::Vector2d>& verts, std::vector<Eigen::Vector2d>& vels);       // (1,0)
        static void linear_x(int n, std::vector<Eigen::Vector2d>& verts, std::vector<Eigen::Vector2d>& vels);       // (x,0)
        static void linear_y(int n, std::vector<Eigen::Vector2d>& verts, std::vector<Eigen::Vector2d>& vels);       // (y,0)
        static void harmonic_1(int n, std::vector<Eigen::Vector2d>& verts, std::vector<Eigen::Vector2d>& vels);     // (y,x)
        static void harmonic_2(int n, std::vector<Eigen::Vector2d>& verts, std::vector<Eigen::Vector2d>& vels);     // (x,-y)
        static void harmonic_3(int n, std::vector<Eigen::Vector2d>& verts, std::vector<Eigen::Vector2d>& vels);     // (-x^2+y^2,2xy)
};
#endif