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

    minFaceLength = 0.7*calc_avg_face_length();
    maxFaceLength = 1.3*calc_avg_face_length();
}

LiquidMesh::LiquidMesh(const std::vector<Eigen::Vector2d>& in_verts, const std::vector<Eigen::Vector2i>& in_faces, const std::vector<Eigen::Vector2d>& in_vels) :
    Mesh(in_verts, in_faces, in_vels) 
{
    vels_solid = std::vector<Eigen::Vector2d>(verts.size(), Eigen::Vector2d(0.0, 0.0));
    
    // defaults to purely liquid-air boundary?
    is_air = std::vector<bool>(verts.size(), true);
    is_solid = std::vector<bool>(verts.size(), false);
    is_triple = std::vector<bool>(verts.size(), false);

    minFaceLength = 0.7*calc_avg_face_length();
    maxFaceLength = 1.4*calc_avg_face_length();
}

void LiquidMesh::set_boundaries(std::vector<bool> air, std::vector<bool> solid, std::vector<bool> triple){
    is_air = air;
    is_solid = solid;
    is_triple = triple;
}


void LiquidMesh::remesh(){
    edge_split();
    edge_collapse();
    laplacian_smoothing();
}

void LiquidMesh::laplacian_smoothing()
{
    std::vector<Eigen::Vector2d> C = std::vector<Eigen::Vector2d>(verts.size(), Eigen::Vector2d(0.0, 0.0));
    std::vector<Eigen::Vector2d> proj = std::vector<Eigen::Vector2d>(verts.size(), Eigen::Vector2d(0.0, 0.0));
    std::vector<Eigen::Vector2d> verts_smoothed = std::vector<Eigen::Vector2d>(verts.size(), Eigen::Vector2d(0.0, 0.0));

    Eigen::Vector2d n_i;
    for (size_t i=0; i<verts.size(); i++){
        // we do not smooth/move triple points
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

void LiquidMesh::edge_collapse(){
    //std::cout<<"Min edge length: "<<minFaceLength<<std::endl;
    size_t n_collapse = 0;
    size_t n_old = faces.size();
    std::vector<bool> verts_to_delete(verts.size(), false);
    std::vector<bool> faces_to_delete(faces.size(), false);
    for (size_t i=0; i<n_old; i++){
    //std::cout<<"Edge length: "<<face_length(i)<<std::endl;
        if (face_length(i) < minFaceLength && !faces_to_delete[i]){

            int endpt_a_ind = faces[i].x();
            int endpt_b_ind = faces[i].y();

            int endpt_a_face = other_face_from_vert(endpt_a_ind, i);
            int endpt_b_face = other_face_from_vert(endpt_b_ind, i);

            if (faces_to_delete[endpt_a_face] || faces_to_delete[endpt_b_face] ){
                continue;
                // maybe not the best logic, but maybe reasonable?
            }

            //std::cout<<"WAHH"<<std::endl;
            n_collapse++;

            if (is_triple[endpt_b_ind] && !is_triple[endpt_a_ind]){
                // merge current face with endpt_a_face
                faces[i].x() = other_vert_from_face(endpt_a_face, endpt_a_ind); //check this
                verts_to_delete[endpt_a_ind] = true;
                faces_to_delete[endpt_a_face] = true;
            } else if (is_triple[endpt_a_ind] && !is_triple[endpt_b_ind]){
                // merge current face with endpt_b_face
                faces[i].y() = other_vert_from_face(endpt_b_face, endpt_b_ind); //check this
                verts_to_delete[endpt_b_ind] = true;
                faces_to_delete[endpt_b_face] = true;
            } else { //either both endpoints are triple points, or both are not
                // -- a -- b --
                if (face_length(endpt_a_face)>face_length(endpt_b_face)) {
                    // merge current face with endpt_b_face
                    faces[i].y() = other_vert_from_face(endpt_b_face, endpt_b_ind); //check this
                    verts_to_delete[endpt_b_ind] = true;
                    faces_to_delete[endpt_b_face] = true;
                } else if (face_length(endpt_a_face)<face_length(endpt_b_face)){
                    // merge current face with endpt_a_face
                    faces[i].x() = other_vert_from_face(endpt_a_face, endpt_a_ind); //check this
                    verts_to_delete[endpt_a_ind] = true;
                    faces_to_delete[endpt_a_face] = true;
                } else{
                    // collapse into the midpoint of the face
                    // implemented by moving endpoint a of my current face, and then merging current face with endpt_b_face
                    verts[endpt_a_ind] = 0.5*(verts[endpt_a_ind]+verts[endpt_b_ind]);
                    faces[i].y() = other_vert_from_face(endpt_b_face, endpt_b_ind); //check this
                    verts_to_delete[endpt_b_ind] = true;
                    faces_to_delete[endpt_b_face] = true;
                }
            }
            
            // if we are also not already marked for deletion
            // we are face A
            // so we determine which face to glom it into 
            // by finding the longer one-- face B
            // modify face A to now have the end point of face B
            // 'mark' the vertex that needs to be deleted for deletion
            // 'mark' face B for deletion
        }
    }
    /*std::cout<<"Verts: "<<std::endl;
    for(size_t i=0; i<verts.size(); i++){
        std::cout<<"("<<verts[i].x()<<","<<verts[i].y()<<")"<<std::endl;
    }
    std::cout<<"Faces post initial finangling: "<<std::endl;
    for(size_t i=0; i<faces.size(); i++){
        std::cout<<"("<<faces[i].x()<<","<<faces[i].y()<<")"<<std::endl;
    }
    std::cout<<"Faces to delete: "<<std::endl;
    for(size_t i=0; i<faces_to_delete.size(); i++){
        std::cout<<"("<<faces_to_delete[i]<<")"<<std::endl;
    }
    std::cout<<"Verts to delete: "<<std::endl;
    for(size_t i=0; i<verts_to_delete.size(); i++){
        std::cout<<"("<<verts_to_delete[i]<<")"<<std::endl;
    }*/

    if (n_collapse > 0){
        std::vector<Eigen::Vector2d> verts_remeshed;
        verts_remeshed.reserve(n_old-n_collapse);
        std::vector<Eigen::Vector2d> vels_remeshed;
        vels_remeshed.reserve(n_old-n_collapse);
        std::vector<Eigen::Vector2d> vels_solid_remeshed;
        vels_solid_remeshed.reserve(n_old-n_collapse);

        std::vector<bool> is_air_remeshed;
        is_air_remeshed.reserve(n_old-n_collapse);
        std::vector<bool> is_solid_remeshed;
        is_solid_remeshed.reserve(n_old-n_collapse);
        std::vector<bool> is_triple_remeshed;
        is_triple_remeshed.reserve(n_old-n_collapse);

        std::unordered_map<size_t, int> verts_map; // maps old

        for(size_t i=0; i<n_old; i++){
            if (!verts_to_delete[i]){
                verts_remeshed.emplace_back(verts[i]);
                verts_map[i] = verts_remeshed.size()-1;
                vels_remeshed.emplace_back(vels[i]);
                vels_solid_remeshed.emplace_back(vels_solid[i]);

                is_air_remeshed.push_back(is_air[i]);
                is_solid_remeshed.push_back(is_solid[i]);
                is_triple_remeshed.push_back(is_triple[i]);
            } else{
                verts_map[i] = -1;
            }
        }

        /*std::cout<<"Verts map: "<<std::endl;
        for(auto tmp: verts_map){
            std::cout<<"("<<tmp.first<<"->"<<tmp.second<<")"<<std::endl;
        }*/

        std::vector<Eigen::Vector2i> faces_remeshed;
        faces_remeshed.reserve(n_old-n_collapse);
        for(size_t i=0; i<n_old; i++){
            if(!faces_to_delete[i]){
                faces_remeshed.emplace_back(verts_map[faces[i][0]],verts_map[faces[i][1]]);
                //faces_remeshed[i][0] = verts_map[faces_remeshed[i][0]];
                //faces_remeshed[i][1] = verts_map[faces_remeshed[i][1]];
            }
        }
        verts = verts_remeshed;
        faces = faces_remeshed;
        vels = vels_remeshed;
        vels_solid = vels_solid_remeshed;

        is_air = is_air_remeshed;
        is_solid = is_solid_remeshed;
        is_triple = is_triple_remeshed;
    }
    assert(faces.size() == n_old - n_collapse);
    assert(verts.size() == n_old - n_collapse);
    assert(vels.size() == verts.size());
    assert(vels_solid.size() == verts.size());

    assert(is_air.size() == verts.size());
    assert(is_solid.size() == verts.size());
    assert(is_triple.size() == verts.size());

    update_neighbor_face_vecs();
    
}
void LiquidMesh::edge_split(){
    //std::cout<<"Max edge length: "<<maxFaceLength<<std::endl;
    size_t n_split = 0;
    size_t n_old = faces.size();
    for (size_t i=0; i<n_old; i++){
        if (face_length(i) > maxFaceLength){
            n_split++;

            verts.emplace_back( 0.5 * (verts[verts_from_face(i).x()] + verts[verts_from_face(i).y()]) );
            vels.emplace_back( 0.5 * (vels[verts_from_face(i).x()] + vels[verts_from_face(i).y()]) );
            vels_solid.emplace_back(0.5* (vels_solid[verts_from_face(i).x()] + vels_solid[verts_from_face(i).y()]) );

            /*std::cout<<"Vert placed at: "<<verts[verts.size()-1][0]<<","<<verts[verts.size()-1][1]<<std::endl;
            std::cout<<"Endpt: "<<verts[verts_from_face(i).x()][0]<<","<<verts[verts_from_face(i).x()][1]<<std::endl;
            std::cout<<"Endpt: "<<verts[verts_from_face(i).y()][0]<<","<<verts[verts_from_face(i).y()][1]<<std::endl;*/

            int newVertIndex = verts.size()-1;
            faces.emplace_back(Eigen::Vector2i(newVertIndex, faces[i].y()));
            faces[i].y() = newVertIndex;

            if(is_solid[verts_from_face(i).x()] || is_solid[verts_from_face(i).y()]){
                is_air.push_back(false);
                is_solid.push_back(true);
                is_triple.push_back(false);
            } else{
                is_air.push_back(true);
                is_solid.push_back(false);
                is_triple.push_back(false);
            }
        }
    }
    assert(faces.size() == n_old + n_split);
    assert(verts.size() == n_old + n_split);
    assert(vels.size() == verts.size());
    assert(vels_solid.size() == verts.size());

    assert(is_air.size() == verts.size());
    assert(is_solid.size() == verts.size());
    assert(is_triple.size() == verts.size());
    
    update_neighbor_face_vecs();
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

void LiquidMesh::reset_boundary_types(){
    is_air = std::vector<bool>(is_air.size(), true);
    is_solid = std::vector<bool>(is_solid.size(), false);
    is_triple = std::vector<bool>(is_triple.size(), false);
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

            // vels[i] = w.calcEffectiveVel(); //?
        }
    }
}

void LiquidMesh::collide_with_solid(SolidMesh& s){
    for(size_t i=0; i<verts.size(); i++){
        //std::cout<<"og: "<<verts[i]<<std::endl;
        if(s.checkCollisionAndSnap(verts[i]) == true){

            //std::cout<<"snapped: "<<verts[i]<<std::endl;
            vels_solid[i] = s.v_effective;
            
            is_solid[i] = true; 
            is_air[i] = false;
            is_triple[i] = false;

            //vels[i] =  s.v_effective; //?
        }
    }
}