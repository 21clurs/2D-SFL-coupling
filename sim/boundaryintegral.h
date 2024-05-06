#ifndef BOUNDARY_INTEGRAL_H
#define BOUNDARY_INTEGRAL_H

#include <math.h>
#include <Eigen/Dense>

class BoundaryIntegral
{
    public:
        static const std::vector<Eigen::Vector2d> gaussian_quadrature(int N=4);     // [-1, 1]
        static const std::vector<Eigen::Vector2d> tanh_sinh_quadrature(int N=9);    // [-1, 1]

        static double G (const Eigen::Vector2d & x, const Eigen::Vector2d  & y) { return -(1 / (2 * M_PI)) * log((x - y).norm()); }
};
#endif