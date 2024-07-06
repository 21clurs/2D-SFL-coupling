#include "rigidbody.h"

RigidBody::RigidBody() : SolidMesh(){
    rotationTheta = 0;
    translation = Eigen::Vector2d(0,0);
    updateRigidBodyVars();
}

RigidBody::RigidBody(const std::vector<Eigen::Vector2d>& in_verts, const std::vector<Eigen::Vector2i>& in_faces) : SolidMesh(in_verts, in_faces){
    rotationTheta = 0;
    translation = Eigen::Vector2d(0,0);
    updateRigidBodyVars();
}

RigidBody::RigidBody(const std::vector<Eigen::Vector2d>& in_verts, const std::vector<Eigen::Vector2i>& in_faces, const std::vector<Eigen::Vector2d>& in_vels) : SolidMesh(in_verts, in_faces, in_vels){
    rotationTheta = 0;
    translation = Eigen::Vector2d(0,0);
    updateRigidBodyVars();
}

void RigidBody::setRotation(double theta){
    rotationTheta = theta;
    updateRotationMat();
}

void RigidBody::updateRigidBodyVars(){
    calculateArea();
    calculateCOM();
    calculateMOI();
    updateRotationMat();
}

void RigidBody::calculateArea(){
    double sum = 0;
    for (size_t i = 0; i<faces.size(); i++){
        Eigen::Vector2d v_a = verts[verts_from_face(i)[0]];
        Eigen::Vector2d v_b = verts[verts_from_face(i)[1]];
        sum += v_a.x()*v_b.y() - v_b.x()*v_a.y();
    }
    area = 0.5*sum;
}

void RigidBody::calculateCOM(){
    double sum_x = 0;
    double sum_y = 0;
    for (size_t i = 0; i<faces.size(); i++){
        Eigen::Vector2d v_a = verts[verts_from_face(i)[0]];
        Eigen::Vector2d v_b = verts[verts_from_face(i)[1]];
        sum_x += (v_a.x()+v_b.x()) * (v_a.x()*v_b.y() - v_b.x()*v_a.y());
        sum_y += (v_a.y()+v_b.y()) * (v_a.x()*v_b.y() - v_b.x()*v_a.y());
    }
    com = Eigen::Vector2d( (1/(6*area))*sum_x, (1/(6*area))*sum_y );
}

void RigidBody::calculateMOI(){
    double I_x = 0;
    double I_y = 0;
    for (size_t i = 0; i<faces.size(); i++){
        Eigen::Vector2d v_a = verts[verts_from_face(i)[0]];
        Eigen::Vector2d v_b = verts[verts_from_face(i)[1]];

        // translate such that COM is at origin
        v_a -= com;
        v_b -= com;

        I_x += (pow(v_a.y(),2) + v_a.y()*v_b.y() + pow(v_b.y(),2)) * (v_a.x()*v_b.y() - v_b.x()*v_a.y());
        I_y += (pow(v_a.x(),2) + v_a.x()*v_b.x() + pow(v_b.x(),2)) * (v_a.x()*v_b.y() - v_b.x()*v_a.y());
    }
    I_x /= 12;
    I_y /= 12;
    moi = I_x + I_y; // placeholder
}

void RigidBody::updateRotationMat(){
    // rotates clockwise by theta
    rotationMat << cos(rotationTheta), -sin(rotationTheta),
                sin(rotationTheta), cos(rotationTheta);
}

void RigidBody::retrieveCurrentVerts(std::vector<Eigen::Vector2d>& v){
    assert(v.size() == verts.size());
    for (size_t i = 0; i<verts.size(); i++){
        v[i] = rotationMat*(verts[i]-com) + com + translation;
    }
}