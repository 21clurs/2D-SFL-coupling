#include "solidmesh.h"
#include <fstream>

SolidMesh::SolidMesh() : Mesh(){ 
    v_effective = Eigen::Vector2d(0,0);
    vel_func = [](double t)->Eigen::Vector2d{ return Eigen::Vector2d(0, 0); };
}
SolidMesh::SolidMesh(const std::vector<Eigen::Vector2d>& in_verts, const std::vector<Eigen::Vector2i>& in_faces) : Mesh(in_verts, in_faces){ 
    v_effective = Eigen::Vector2d(0,0);
    vel_func = [](double t)->Eigen::Vector2d{ return Eigen::Vector2d(0, 0); };
}
SolidMesh::SolidMesh(const std::vector<Eigen::Vector2d>& in_verts, const std::vector<Eigen::Vector2i>& in_faces, const std::vector<Eigen::Vector2d>& in_vels) : Mesh(in_verts, in_faces, in_vels){
    v_effective = Eigen::Vector2d(0,0); 
    vel_func = [](double t)->Eigen::Vector2d{ return Eigen::Vector2d(0, 0); };
}

bool SolidMesh::loadMeshFromFile(SolidMesh &m, std::string infileName){
    // load sim options file
    std::ifstream infile(infileName);
    if (!infile.is_open()) { 
        std::cerr << "Unable to open solid file "<<infileName<<"!" << std::endl; 
        assert(!"Unable to open options file!");
    }

    std::vector<Eigen::Vector2d> v;
    std::vector<Eigen::Vector2i> f;

    std::string line;
    while(!infile.eof()){
        std::getline(infile, line);
        std::stringstream ss(line);

        std::string linetype;
        ss >> linetype;
        if (linetype == "#" || linetype == "" || ss.eof())    // skip comment lines and empty lines
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
        }
        else {
            std::cerr << "Invalid line in "<<infileName<<"! Skipping line..." << std::endl;
        }
    }
    infile.close();
    
    assert(v.size() == f.size());

    m.verts = v;
    m.faces = f;
    m.update_neighbor_face_vecs();
    
    m.vels.resize(v.size());

    return true;
}

bool SolidMesh::checkCollisionAndSnap(Eigen::Vector2d& curr_pt){
    double min_d;
    Eigen::Vector2d nearest_pt;

    // iterate through faces
    for (size_t i = 0; i<faces.size(); i++){

        // retrieve endpoints of current face, these should be oriented
        Eigen::Vector2d ptA = verts[faces[i][0]];
        Eigen::Vector2d ptB = verts[faces[i][1]];

        // first check to snap to 'sharp' features in geometry
        // if we snap here, skip everything else
        if ((curr_pt-ptA).norm()<epsilon){
            curr_pt = ptA;
            return true;
        } else if ((curr_pt-ptB).norm()<epsilon){
            curr_pt = ptB;
            return true;
        }

        // finding closest point on the line defined by the face to curr_pt
        Eigen::Vector2d u = ptB-ptA;        //std::cout<<"u: "<<u.x()<<","<<u.y()<<std::endl;
        Eigen::Vector2d v = curr_pt-ptA;    //std::cout<<"v: "<<v.x()<<","<<v.y()<<std::endl;
        double t = (u.dot(v)/u.dot(u));     //std::cout<<"t: "<<t<<std::endl;
        // finding closest point
        Eigen::Vector2d proj_pt;
        if (t>=0 && t<=1){ 
            // closest point is within the segment length
            proj_pt = (1-t)*ptA + t*ptB;
        } else{ 
            // if outside the segment length, project to nearest endpoint
            double g0 = (ptA - curr_pt).squaredNorm();
            double g1 = (ptB - curr_pt).squaredNorm();
            proj_pt = g0 < g1 ? ptA : ptB;
        }
    
        // logic to get minimal distance to a segment in the SolidMesh
        if (i==0 || (proj_pt-curr_pt).norm()<=abs(min_d)){
            min_d = (proj_pt-curr_pt).norm();
            nearest_pt = proj_pt;
        }
    }

    // since winding number when exactly on the segment is a little funky
    // we define 'collision' as when min_d is less than our epsilon OR when the winding number is nonzero
    // as nonzero winding number indicates that we are intersecting the mesh
    if(min_d <= epsilon || windingNumber(curr_pt)>1.0e-8){
        // snap the point to the nearest face
        curr_pt = nearest_pt;
        return true;
    } else{
        return false;
    }
}

void SolidMesh::collideAndSnap(LiquidMesh& l){
    // for each vertex in the given LiquidMesh
    for(size_t j=0; j<l.verts.size(); j++){
        Eigen::Vector2d curr_pt = l.verts[j];
        bool pt_snapped_to_solid_corner = false;

        double min_d;
        Eigen::Vector2d nearest_pt;

        // iterate through faces
        for (size_t i = 0; i<faces.size(); i++){

            // retrieve endpoints of current face, these should be oriented
            Eigen::Vector2d ptA = verts[faces[i][0]];
            Eigen::Vector2d ptB = verts[faces[i][1]];

            // first check to snap to 'sharp' features in geometry, if so, we can skip the rest of the loop
            if ((curr_pt-ptA).norm()<epsilon){
                min_d = (curr_pt-ptA).norm();
                nearest_pt = ptA;
                pt_snapped_to_solid_corner = true;
                break;
            } else if ((curr_pt-ptB).norm()<epsilon){
                min_d = (curr_pt-ptB).norm();
                nearest_pt = ptB;
                pt_snapped_to_solid_corner = true;
                break;
            } else {
                // finding closest point on the line defined by the face to curr_pt
                Eigen::Vector2d u = ptB-ptA;        //std::cout<<"u: "<<u.x()<<","<<u.y()<<std::endl;
                Eigen::Vector2d v = curr_pt-ptA;    //std::cout<<"v: "<<v.x()<<","<<v.y()<<std::endl;
                double t = (u.dot(v)/u.dot(u));     //std::cout<<"t: "<<t<<std::endl;
                // finding closest point
                Eigen::Vector2d proj_pt;
                if (t>=0 && t<=1){ 
                    // closest point is within the segment length
                    proj_pt = (1-t)*ptA + t*ptB;
                } else{ 
                    // if outside the segment length, project to nearest endpoint
                    double g0 = (ptA - curr_pt).squaredNorm();
                    double g1 = (ptB - curr_pt).squaredNorm();
                    proj_pt = g0 < g1 ? ptA : ptB;
                }
            
                // logic to get minimal distance to a segment in the SolidMesh
                if (i==0 || (proj_pt-curr_pt).norm()<=abs(min_d)){
                    min_d = (proj_pt-curr_pt).norm();
                    nearest_pt = proj_pt;
                }
            }
        }

        // since winding number when precisely on the segment is a bit misbehaved
        // we define 'collision' as when min_d is less than our epsilon OR when the winding number is nonzero
        // as a nonzero winding number indicates that we are intersecting the mesh
        if(min_d <= epsilon || windingNumber(curr_pt)>1.0e-8){
            // snap the point to the nearest face
            l.verts[j] = nearest_pt;
            // set effective velocity
            l.vels_solid[j] = v_effective;
            l.set_boundaries_for_vertex(j, false, true, false, pt_snapped_to_solid_corner);
        }
    }
}

void SolidMesh::setVelFunc(std::function<Eigen::Vector2d(double)> func){ vel_func = func; }

void SolidMesh::advectFE(double curr_t, double dt){
    for (size_t i=0; i< verts.size(); i++){
        verts[i] += vel_func(curr_t)*dt;
    }
    v_effective = vel_func(curr_t);//*dt;
}

void SolidMesh::advectFE(double dt){
    assert(verts.size() == vels.size());
    for (size_t i=0; i< verts.size(); i++){
        v_effective = vels[i]*dt;
        verts[i] += vels[i]*dt;
    }
}