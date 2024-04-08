#ifndef BOUNDARY_INTEGRAL_H
#define BOUNDARY_INTEGRAL_H

#include <math.h>
#include <Eigen/Dense>

class BoundaryIntegral
{
    public:
        static void gaussian_quadrature(int N=4);
        static void tanh_sinh_quadrature(int N=9);

        static double G (const Eigen::Vector2d & x, const Eigen::Vector2d  & y) { return -(1 / (2 * M_PI)) * log((x - y).norm()); }
};
#endif