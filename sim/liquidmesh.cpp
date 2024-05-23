#include "liquidmesh.h"
#include <iostream>

LiquidMesh::LiquidMesh(const std::vector<Eigen::Vector2d>& in_verts, const std::vector<Eigen::Vector2i>& in_faces) :
    Mesh(in_verts, in_faces) 
{
    vels_solid = std::vector<Eigen::Vector2d>(verts.size(), Eigen::Vector2d(0.0, 0.0));

    // defaults to purely liquid-air boundary?
    is_air = std::vector<bool>(verts.size(), true);
    is_solid = std::vector<bool>(verts.size(), false);
    is_triple = std::vector<bool>(verts.size(), false);
}

LiquidMesh::LiquidMesh(const std::vector<Eigen::Vector2d>& in_verts, const std::vector<Eigen::Vector2i>& in_faces, const std::vector<Eigen::Vector2d>& in_vels) :
    Mesh(in_verts, in_faces, in_vels) 
{
    vels_solid = std::vector<Eigen::Vector2d>(verts.size(), Eigen::Vector2d(0.0, 0.0));
    
    // defaults to purely liquid-air boundary?
    is_air = std::vector<bool>(verts.size(), true);
    is_solid = std::vector<bool>(verts.size(), false);
    is_triple = std::vector<bool>(verts.size(), false);
}

void LiquidMesh::set_boundaries(std::vector<bool> air, std::vector<bool> solid, std::vector<bool> triple){
    is_air = air;
    is_solid = solid;
    is_triple = triple;
}

void LiquidMesh::laplacian_smoothing()
{
    std::vector<Eigen::Vector2d> C = std::vector<Eigen::Vector2d>(verts.size(), Eigen::Vector2d(0.0, 0.0));
    std::vector<Eigen::Vector2d> proj = std::vector<Eigen::Vector2d>(verts.size(), Eigen::Vector2d(0.0, 0.0));
    std::vector<Eigen::Vector2d> verts_smoothed = std::vector<Eigen::Vector2d>(verts.size(), Eigen::Vector2d(0.0, 0.0));

    Eigen::Vector2d n_i;
    for (size_t i=0; i<verts.size(); i++){
        if (is_triple[i]){
            C[i] = verts[i];
            proj[i] = Eigen::Vector2d(0.0, 0.0);
        }
        else{
            n_i = calc_vertex_normal(i);
            C[i] = 0.5*(prev_neighbor(i) + next_neighbor(i));
            proj[i] = (n_i*n_i.transpose())*(verts[i]-C[i]);
        }
        verts_smoothed[i] = C[i] + proj[i];
    }
    verts = verts_smoothed;
}

std::vector<bool> LiquidMesh::get_solid_faces(){
    std::vector<bool> solid_faces = std::vector<bool>(faces.size(), false);
    Eigen::Vector2i endpts;
    for (size_t i=0; i<faces.size(); i++){
        endpts = verts_from_face(i);
        solid_faces[i] = (is_solid[endpts[0]] || is_solid[endpts[1]]);
    }
    return solid_faces;
}

void LiquidMesh::update_triple_points(){
    for(size_t i=0; i<verts.size(); i++){
        if (is_solid[i] && (is_air[next_neighbor_index(i)] || is_air[prev_neighbor_index(i)])){
            is_triple[i] = true;
            is_solid[i] = false;
        }
        assert((is_solid[i] ^ is_air[i] ^ is_triple[i]) && (is_solid[i] + is_air[i] + is_triple[i]==1));
    }
}

void LiquidMesh::collide_with_wall(WallObject& w){
    for(size_t i=0; i<verts.size(); i++){
        //std::cout<<"og: "<<verts[i]<<std::endl;
        if(w.checkCollisionAndSnap(verts[i]) == true){

            //std::cout<<"snapped: "<<verts[i]<<std::endl;
            vels_solid[i] = w.calcEffectiveVel();
            
            is_solid[i] = true; 
            is_air[i] = false;
            is_triple[i] = false;
            //std::cout<<"bLAH "<<i<<std::endl;
        }
    }
    update_triple_points();
}

void LiquidMesh::collide_with_solid(SolidMesh& s){
    for(size_t i=0; i<verts.size(); i++){
        //std::cout<<"og: "<<verts[i]<<std::endl;
        if(s.checkCollisionAndSnap(verts[i]) == true){

            //std::cout<<"snapped: "<<verts[i]<<std::endl;
            vels_solid[i] = Eigen::Vector2d(0,0);//placeholder
            
            is_solid[i] = true; 
            is_air[i] = false;
            is_triple[i] = false;
            //std::cout<<"bLAH "<<i<<std::endl;
        }
    }
    update_triple_points();
}