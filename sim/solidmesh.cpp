#include "solidmesh.h"

SolidMesh::SolidMesh(const std::vector<Eigen::Vector2d>& in_verts, const std::vector<Eigen::Vector2i>& in_faces) : Mesh(in_verts, in_faces){ 
    v_effective = Eigen::Vector2d(0,0);
    vel_func = [](double t)->Eigen::Vector2d{ return Eigen::Vector2d(0, 0); };
}
SolidMesh::SolidMesh(const std::vector<Eigen::Vector2d>& in_verts, const std::vector<Eigen::Vector2i>& in_faces, const std::vector<Eigen::Vector2d>& in_vels) : Mesh(in_verts, in_faces, in_vels){
    v_effective = Eigen::Vector2d(0,0); 
    vel_func = [](double t)->Eigen::Vector2d{ return Eigen::Vector2d(0, 0); };
}

bool SolidMesh::checkCollisionAndSnap(Eigen::Vector2d& curr_pt){
    // iterate through faces
    double min_d;
    Eigen::Vector2d nearest_pt;
    for (size_t i =0; i<faces.size(); i++){
        // check if  curr_pt is on consistent side of all faces (i.e. "inside")
        Eigen::Vector2d ptA = verts[faces[i].x()];
        Eigen::Vector2d ptB = verts[faces[i].y()];

        // finding (signed) distance from curr_pt to line, regardless of endpoints
        double d_orth = ((curr_pt.x()-ptA.x())*(ptB.y()-ptA.y()) - (curr_pt.y()-ptA.y())*(ptB.x()-ptA.x())) / ((ptB-ptA).norm());

        //std::cout<<"("<<ptA[0]<<", "<<ptA[1]<<")"<<"("<<ptB[0]<<", "<<ptB[1]<<")"<<" "<<d_orth<<std::endl;

        if( d_orth>0 && abs(d_orth)>epsilon){
            // this point is not inside the solid mesh, can end here
            return false;
        } else{    
            // finding closest point on face to curr_pt
            Eigen::Vector2d u = ptB-ptA;
            //std::cout<<"u: "<<u.x()<<","<<u.y()<<std::endl;
            Eigen::Vector2d v = curr_pt-ptA;
            //std::cout<<"v: "<<v.x()<<","<<v.y()<<std::endl;
            double t = (u.dot(v)/u.dot(u));
            //std::cout<<"t: "<<t<<std::endl;
            Eigen::Vector2d proj_pt;
            
            if (t>=0 && t<=1){
                proj_pt = (1-t)*ptA + t*ptB;
            } else{ // project to nearest endpoint
                //std::cout<<"?"<<std::endl;
                double g0 = (ptA - curr_pt).squaredNorm();
                double g1 = (ptB - curr_pt).squaredNorm();
                proj_pt = g0 < g1 ? ptA : ptB;
            }

            //std::cout<<"proj pt: "<<proj_pt[0]<<", "<<proj_pt[1]<<" : "<<t<<std::endl;
        
            // if this is reached, we are (at least currently) inside the solid mesh, or within the epsilon 'snapping' distance
            if (i==0 || (proj_pt-curr_pt).norm()<=abs(min_d)){
                min_d = (proj_pt-curr_pt).norm();
                nearest_pt = proj_pt;// need to figure out how to do this
            }
        }
    }
    // if this is reached, then currpt is colliding with the solid mesh
    // snap the point to the nearest face
    curr_pt = nearest_pt;
    return true;
}

void SolidMesh::setVelFunc(std::function<Eigen::Vector2d(double)> func){ vel_func = func; }

void SolidMesh::advectFE(double curr_t, double dt){
    for (size_t i=0; i< verts.size(); i++){
        verts[i] += vel_func(curr_t)*dt;
    }
    v_effective = vel_func(curr_t)*dt;
}