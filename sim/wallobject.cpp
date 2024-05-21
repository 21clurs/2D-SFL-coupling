#include "wallobject.h"

WallObject::WallObject(double dt, Eigen::Vector2d& a, Eigen::Vector2d& b): endptA(a), endptB(b), dt(dt), endptA_prev(a), endptB_prev(b){
    // 90-deg clockwise rotation matrix
    Eigen::Matrix2d rotate;
    rotate << 0, 1, -1, 0;
    n = rotate * (endptB-endptA);
    n.normalize();

    update_effective_endpts();
}

void WallObject::setVelFunc(std::function<Eigen::Vector2d(double)> func){ vel_func = func; }

void WallObject::setEps(double e){ epsilon = e; }

Eigen::Vector2d WallObject::calcEffectiveVel(){
    return (endptA-endptA_prev)/dt;
}

void WallObject::advectFE(double curr_t){
    endptA_prev = endptA;
    endptB_prev = endptB;

    endptA = endptA + vel_func(curr_t)*dt;
    endptB = endptB + vel_func(curr_t)*dt;

    update_effective_endpts();
    //std::cout<<"Advect wall! "<<endptA.x()<<","<<endptA.y()<<std::endl;
}

bool WallObject::checkCollisionAndSnap(Eigen::Vector2d& x){
    double d_to_effective = (x.x()-endptA_effective.x())*(endptB_effective.y()-endptA_effective.y()) - (x.y()-endptA_effective.y())*(endptB_effective.x()-endptA_effective.x());
    //std::cout<<d<<std::endl;
    //std::cout<<"hello: "<<x[0]<<", "<<x[1]<<","<<d_to_effective<<" "<<endptB_effective.x()<<","<<endptB_effective.y()<<std::endl;
    if (d_to_effective<=0){
        double d = ((x.x()-endptA.x())*(endptB.y()-endptA.y()) - (x.y()-endptA.y())*(endptB.x()-endptA.x())) / (endptB-endptA).norm();
        //std::cout<<"hello again: "<<d<<" "<<n<<std::endl;
        x = x - d*n;
        //x.y() = endptA.y();
        return true;
    }
    return false;
}

void WallObject::update_effective_endpts(){
    endptA_effective = endptA+epsilon*n;
    endptB_effective = endptB+epsilon*n;
}