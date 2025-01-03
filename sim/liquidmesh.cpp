#include "liquidmesh.h"
#include <fstream>
#include <iostream>
#include "simoptions.h"

LiquidMesh::LiquidMesh() : Mesh(){
    vels_solid = std::vector<Eigen::Vector2d>(verts.size(), Eigen::Vector2d(0.0, 0.0));// defaults to purely liquid-air boundary?
    
    is_air = std::vector<bool>(verts.size(), true);
    is_solid = std::vector<bool>(verts.size(), false);
    is_triple = std::vector<bool>(verts.size(), false);

    is_corner = std::vector<bool>(verts.size(), false);

    per_vertex_rb_contact = std::vector<uint>(verts.size(), 0);
}

LiquidMesh::LiquidMesh(const std::vector<Eigen::Vector2d>& in_verts, const std::vector<Eigen::Vector2i>& in_faces) :
    Mesh(in_verts, in_faces) 
{
    vels_solid = std::vector<Eigen::Vector2d>(verts.size(), Eigen::Vector2d(0.0, 0.0));

    // defaults to purely liquid-air boundary?
    is_air = std::vector<bool>(verts.size(), true);
    is_solid = std::vector<bool>(verts.size(), false);
    is_triple = std::vector<bool>(verts.size(), false);

    is_corner = std::vector<bool>(verts.size(), false);

    per_vertex_rb_contact = std::vector<uint>(verts.size(), 0);

    reset_face_length_limits();
}

LiquidMesh::LiquidMesh(const std::vector<Eigen::Vector2d>& in_verts, const std::vector<Eigen::Vector2i>& in_faces, const std::vector<Eigen::Vector2d>& in_vels) :
    Mesh(in_verts, in_faces, in_vels) 
{
    vels_solid = std::vector<Eigen::Vector2d>(verts.size(), Eigen::Vector2d(0.0, 0.0));
    
    // defaults to purely liquid-air boundary?
    is_air = std::vector<bool>(verts.size(), true);
    is_solid = std::vector<bool>(verts.size(), false);
    is_triple = std::vector<bool>(verts.size(), false);

    is_corner = std::vector<bool>(verts.size(), false);

    per_vertex_rb_contact = std::vector<uint>(verts.size(), 0);

    reset_face_length_limits();
}

void LiquidMesh::resize_mesh(size_t n){
    verts.resize(n);
    faces.resize(n);
    vels.resize(n);

    vertsPrevFace.resize(n);
    vertsNextFace.resize(n);

    vels_solid.resize(n);
    
    is_air.resize(n);
    is_solid.resize(n);
    is_triple.resize(n);

    is_corner.resize(n);

    per_vertex_rb_contact.resize(n);
}

void LiquidMesh::set_boundaries(std::vector<bool> air, std::vector<bool> solid, std::vector<bool> triple){
    is_air = air;
    is_solid = solid;
    is_triple = triple;
}

void LiquidMesh::set_boundaries(std::vector<bool> air, std::vector<bool> solid, std::vector<bool> triple, std::vector<bool> corner){
    is_air = air;
    is_solid = solid;
    is_triple = triple;
    is_corner = corner;
}

void LiquidMesh::set_boundaries_for_vertex(int i, bool air, bool solid, bool triple, bool corner){
    is_air[i] = air;
    is_solid[i] = solid;
    is_triple[i] = triple;

    is_corner[i] = corner;
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
        // we do not smooth/move triple points or points at corners
        if (is_triple[i]|| is_corner[i]){
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
        int endpt_a_ind = faces[i].x();
        int endpt_a_face = other_face_from_vert(endpt_a_ind, i);
        if (face_length(endpt_a_face) == 0 && !faces_to_delete[endpt_a_face]){ 
            // NEED to delete this, and should not be disruptive to do so
            n_collapse++;

            //int endpt_a_ind = faces[i].x();
            //int endpt_a_face = other_face_from_vert(endpt_a_ind, i);
            
            // merge current face with endpt_a_face
            faces[i].x() = other_vert_from_face(endpt_a_face, endpt_a_ind);
            verts_to_delete[endpt_a_ind] = true;
            faces_to_delete[endpt_a_face] = true;

            //std::cout<<"A"<<endpt_a_face<<std::endl;

        } else if (face_length(i) < minFaceLength && !faces_to_delete[i]){
            // slightly more normal case

            endpt_a_ind = faces[i].x();
            int endpt_b_ind = faces[i].y();

            endpt_a_face = other_face_from_vert(endpt_a_ind, i);
            int endpt_b_face = other_face_from_vert(endpt_b_ind, i);

            // don't bother deleting on this iteration if either of its neighbor faces are already marked
            // this is just because I don't feel like going through all the cases there
            if (faces_to_delete[endpt_a_face] || faces_to_delete[endpt_b_face] ){
                continue;
            } else if (verts_to_delete[endpt_a_ind] || verts_to_delete[endpt_b_ind] ){
                continue;
            } 
            // we also do not delete faces around the corners currently
            if (is_corner[endpt_a_ind] || is_corner[endpt_b_ind] ){
                continue;
            }

            n_collapse++;

            if (is_triple[endpt_b_ind] && !is_triple[endpt_a_ind]){
                // merge current face with endpt_a_face
                faces[i].x() = other_vert_from_face(endpt_a_face, endpt_a_ind); //check this
                verts_to_delete[endpt_a_ind] = true;
                faces_to_delete[endpt_a_face] = true;
                //std::cout<<"B"<<endpt_a_face<<std::endl;
            } else if (is_triple[endpt_a_ind] && !is_triple[endpt_b_ind]){
                // merge current face with endpt_b_face
                faces[i].y() = other_vert_from_face(endpt_b_face, endpt_b_ind); //check this
                verts_to_delete[endpt_b_ind] = true;
                faces_to_delete[endpt_b_face] = true;
                //std::cout<<"C"<<endpt_b_face<<std::endl;
            } else { //either both endpoints are triple points, or both are not
                // -- a -- b --
                if (face_length(endpt_a_face)>face_length(endpt_b_face)) {
                    // merge current face with endpt_b_face
                    faces[i].y() = other_vert_from_face(endpt_b_face, endpt_b_ind); //check this
                    verts_to_delete[endpt_b_ind] = true;
                    faces_to_delete[endpt_b_face] = true;
                    //std::cout<<"D"<<endpt_b_face<<std::endl;
                } else if (face_length(endpt_a_face)<face_length(endpt_b_face)){
                    // merge current face with endpt_a_face
                    faces[i].x() = other_vert_from_face(endpt_a_face, endpt_a_ind); //check this
                    verts_to_delete[endpt_a_ind] = true;
                    faces_to_delete[endpt_a_face] = true;
                    //std::cout<<"E"<<endpt_a_face<<std::endl;
                } else{
                    // collapse into the midpoint of the face
                    // implemented by moving endpoint a of my current face, and then merging current face with endpt_b_face
                    verts[endpt_a_ind] = 0.5*(verts[endpt_a_ind]+verts[endpt_b_ind]);
                    faces[i].y() = other_vert_from_face(endpt_b_face, endpt_b_ind); //check this
                    verts_to_delete[endpt_b_ind] = true;
                    faces_to_delete[endpt_b_face] = true;
                    //std::cout<<"F"<<endpt_b_face<<std::endl;
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

        std::vector<bool> is_corner_remeshed;
        is_corner_remeshed.reserve(n_old-n_collapse);

        std::unordered_map<size_t, int> verts_map; // maps old to new verts

        for(size_t i=0; i<n_old; i++){
            if (!verts_to_delete[i]){
                verts_remeshed.emplace_back(verts[i]);
                verts_map[i] = verts_remeshed.size()-1;
                vels_remeshed.emplace_back(vels[i]);
                vels_solid_remeshed.emplace_back(vels_solid[i]);

                is_air_remeshed.push_back(is_air[i]);
                is_solid_remeshed.push_back(is_solid[i]);
                is_triple_remeshed.push_back(is_triple[i]);

                is_corner_remeshed.push_back(is_corner[i]);
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

        is_corner = is_corner_remeshed;

        per_vertex_rb_contact.resize(verts.size());
    }
    assert(faces.size() == n_old - n_collapse);
    assert(verts.size() == n_old - n_collapse);
    assert(vels.size() == verts.size());
    assert(vels_solid.size() == verts.size());

    assert(is_air.size() == verts.size());
    assert(is_solid.size() == verts.size());
    assert(is_triple.size() == verts.size());

    assert(is_corner.size() == verts.size());

    assert(per_vertex_rb_contact.size() == verts.size());

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

                is_corner.push_back(false);
            } else{
                is_air.push_back(true);
                is_solid.push_back(false);
                is_triple.push_back(false);

                is_corner.push_back(false);
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

    assert(is_corner.size() == verts.size());

    per_vertex_rb_contact.resize(verts.size());
    assert(per_vertex_rb_contact.size() == verts.size());
    
    update_neighbor_face_vecs();
}

std::vector<bool> LiquidMesh::get_solid_faces(){
    std::vector<bool> solid_faces = std::vector<bool>(faces.size(), false);
    Eigen::Vector2i endpts;
    for (size_t i=0; i<faces.size(); i++){
        endpts = verts_from_face(i);
        // for now this interprets  anything with two triple points at the ends as a solid face...
        //solid_faces[i] = ( (is_solid[endpts[0]]||is_triple[endpts[0]]) && (is_solid[endpts[1]||is_triple[endpts[1]]]) );
        solid_faces[i] = (is_solid[endpts[0]] || is_solid[endpts[1]]) || (is_triple[endpts[0]] && is_triple[endpts[1]]);
    }
    return solid_faces;
}

const Eigen::Vector2d LiquidMesh::calc_vertex_solid_normal(const int vertIndex){
    Eigen::Vector2i adjacent_face_inds = faces_from_vert(vertIndex);
    // figure out which face is the solid face
    // return normal of that face
    // this is NOT super robust. I think I miss a LOT of edge cases here, especially to do with the triple points
    Eigen::Vector2d n_solid;
    if (is_air[other_vert_from_face(adjacent_face_inds[0], vertIndex)]){
        n_solid = calc_face_normal(adjacent_face_inds[1]);
    } else {
        n_solid = calc_face_normal(adjacent_face_inds[0]);
    }
    n_solid.normalize();
    return n_solid;
}

double LiquidMesh::vert_solid_area(const int vertIndex){
    Eigen::Vector2i adjacent_face_inds = faces_from_vert(vertIndex);
    // figure out which face is the solid face
    // return half of area of that face
    // this is NOT super robust. I think I miss a LOT of edge cases here, especially to do with the triple points
    double area_solid;
    if (is_air[other_vert_from_face(adjacent_face_inds[0], vertIndex)]){
        area_solid = 0.5 * face_length(adjacent_face_inds[1]);
    } else {
        area_solid = 0.5 * face_length(adjacent_face_inds[0]);
    }
    return area_solid;
}

void LiquidMesh::update_triple_points(){
    for(size_t i=0; i<verts.size(); i++){
        if (is_solid[i] && (is_air[next_neighbor_index(i)] || is_air[prev_neighbor_index(i)])){
            is_triple[i] = true;
            is_solid[i] = false;
            is_air[i] = false;
        }
        assert((is_solid[i] ^ is_air[i] ^ is_triple[i]) && (is_solid[i] + is_air[i] + is_triple[i]==1));
    }
}

void LiquidMesh::reset_boundary_types(){
    is_air = std::vector<bool>(is_air.size(), true);
    is_solid = std::vector<bool>(is_solid.size(), false);
    is_triple = std::vector<bool>(is_triple.size(), false);

    is_corner = std::vector<bool>(is_corner.size(), false);
}
void LiquidMesh::reset_face_length_limits(){
    minFaceLength = SimOptions::doubleValue("mesh-edge-min-ratio")*calc_avg_face_length();
    maxFaceLength = SimOptions::doubleValue("mesh-edge-max-ratio")*calc_avg_face_length();
}
// returns the closest distance to the liquid mesh
// returns a positive number if inside the mesh
// returns a negative number if outside the mesh
double LiquidMesh::signed_min_dist(Eigen::Vector2d x){
    double min_d;
    Eigen::Vector2d nearest_pt;

    // iterate through faces
    for (size_t i = 0; i<faces.size(); i++){

        // retrieve endpoints of current face, these should be oriented
        Eigen::Vector2d ptA = verts[faces[i][0]];
        Eigen::Vector2d ptB = verts[faces[i][1]];

        // finding closest point on the line defined by the face to curr_pt
        Eigen::Vector2d u = ptB-ptA;        //std::cout<<"u: "<<u.x()<<","<<u.y()<<std::endl;
        Eigen::Vector2d v = x-ptA;    //std::cout<<"v: "<<v.x()<<","<<v.y()<<std::endl;
        double t = (u.dot(v)/u.dot(u));     //std::cout<<"t: "<<t<<std::endl;
        // finding closest point
        Eigen::Vector2d proj_pt;
        if (t>=0 && t<=1){ 
            // closest point is within the segment length
            proj_pt = (1-t)*ptA + t*ptB;
        } else{ 
            // if outside the segment length, project to nearest endpoint
            double g0 = (ptA - x).squaredNorm();
            double g1 = (ptB - x).squaredNorm();
            proj_pt = g0 < g1 ? ptA : ptB;
        }
    
        // logic to get minimal distance to a segment in the SolidMesh
        if (i==0 || (proj_pt-x).norm()<=abs(min_d)){
            min_d = (proj_pt-x).norm();
            nearest_pt = proj_pt;
        }
    }

    // if inside
    if (windingNumber(x)>1.0e-8){
        return abs(min_d);
    } else{
        return -1*abs(min_d);
    }
}

bool LiquidMesh::loadMeshFromFile(LiquidMesh &m, std::string infileName){
    // load liquid mesh file
    std::ifstream infile(infileName);
    if (!infile.is_open()) { 
        std::cerr << "Unable to open liquid file "<<infileName<<"!" << std::endl; 
        assert(!"Unable to open liquid mesh file!");
    }

    std::vector<Eigen::Vector2d> v;
    std::vector<Eigen::Vector2i> f;
    std::vector<Eigen::Vector2d> u;

    std::string line;
    while(!infile.eof()){
        std::getline(infile, line);
        std::stringstream ss(line);

        std::string linetype;
        ss >> linetype;
        if (linetype == "#" || linetype == "//" || linetype == "" || ss.eof())    // skip comment lines and empty lines
            continue;
        
        if (linetype.compare("v") == 0){
            double a,b;
            ss >> a;
            ss >> b;
            v.emplace_back(Eigen::Vector2d(a,b));
        } else if (linetype.compare("f") == 0){
            int a,b;
            ss >> a;
            ss >> b;
            f.emplace_back(Eigen::Vector2i(a,b));
        } else if (linetype.compare("u") == 0){
            double a,b;
            ss >> a;
            ss >> b;
            u.emplace_back(Eigen::Vector2d(a,b));
        }
        else {
            std::cerr << "Invalid line in "<<infileName<<"! Skipping line..." << std::endl;
        }
    }
    infile.close();
    
    assert(v.size() == f.size());

    m.resize_mesh(v.size());

    m.verts = v;
    m.faces = f;
    m.update_neighbor_face_vecs();

    // inputting per-vertex velocity is optional
    if (u.size() > 0){
        assert(u.size() == v.size());
        m.vels = u;
    }

    return true;
}
