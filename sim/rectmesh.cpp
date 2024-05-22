#include "rectmesh.h"

RectMesh::RectMesh(double l, double r, double b, double t):
    left_x(l), right_x(r), bottom_y(b), top_y(t)
{
    update_effective_coords();
}

void RectMesh::setEps(double e){ epsilon = e; }

bool RectMesh::checkCollisionAndSnap(Eigen::Vector2d& x){
    if (x.x()>left_x_effective && x.x()<right_x_effective && x.y()>bottom_y_effective && x.y()<top_y_effective){
        double dx_l = abs(x.x()-left_x);
        double dx_r = abs(x.x()-right_x);
        double dy_b = abs(x.y()-bottom_y);
        double dy_t = abs(x.y()-top_y);

        std::vector<double> d = {dx_l, dx_r, dy_b, dy_t};
        std::vector<double>::iterator min_d = std::min_element(d.begin(), d.end());
        if(dx_l == *min_d){
            x.x() = left_x;
        }else if(dx_r == *min_d){
            x.x() = right_x;
        }else if(dy_b == *min_d){
            x.y() = bottom_y;
        }else if(dy_t == *min_d){
            x.y() = top_y;
        }
        return true;
    }
    return false;
}

void RectMesh::update_effective_coords(){
    left_x_effective = left_x-epsilon;
    right_x_effective = right_x+epsilon;
    bottom_y_effective = bottom_y-epsilon;
    top_y_effective = top_y+epsilon;
}