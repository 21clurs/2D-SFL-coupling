#ifndef RECT_MESH_H
#define RECT_MESH_H

#include <math.h>
#include <Eigen/Dense>
#include <iostream>
#include <algorithm>

class RectMesh
{
    public:

        double left_x;
        double right_x;
        double bottom_y;
        double top_y;

        // constructors
        RectMesh(double l, double r, double b, double t);

        void setEps(double e);
        bool checkCollisionAndSnap(Eigen::Vector2d& x);

    private:
        double left_x_effective;
        double right_x_effective;
        double bottom_y_effective;
        double top_y_effective;

        std::vector<Eigen::Vector2d> vert_normals;
        std::vector<Eigen::Vector2d> face_normals;

        double epsilon = 0.01;

        void update_effective_coords();
};
#endif