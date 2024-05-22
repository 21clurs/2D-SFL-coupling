#ifndef WALL_OBJECT_H
#define WALL_OBJECT_H

#include <math.h>
#include <Eigen/Dense>
#include <iostream>

class WallObject
{
    public:

        Eigen::Vector2d endptA;
        Eigen::Vector2d endptB;
        
        WallObject(double dt, Eigen::Vector2d& a, Eigen::Vector2d& b);
        
        void setVelFunc(std::function<Eigen::Vector2d(double)> func);
        void setEps(double e);
        Eigen::Vector2d calcEffectiveVel();
        void advectFE(double curr_t);

        bool checkCollisionAndSnap(Eigen::Vector2d& x);

    private:
        // "wall" line is defined by two points a and b
        // the side of the line that is "wall" vs not is defined by the order of a and b
        // "inside the wall" is determined as opposite the direction of the clockwise normal to vector (b-a)
        double dt;
        Eigen::Vector2d endptA_prev;
        Eigen::Vector2d endptB_prev;
        Eigen::Vector2d n;
        Eigen::Vector2d endptA_effective;
        Eigen::Vector2d endptB_effective;
        double epsilon = 0.01;

        void update_effective_endpts();

        std::function<Eigen::Vector2d(double)> vel_func;
};
#endif