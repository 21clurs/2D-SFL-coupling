#include "rigidbody.h"

#include <fstream>
#include "sim.h"

RigidBody::RigidBody() : Mesh(){
    rho = INFINITY;
    verts_no_transform = std::vector<Eigen::Vector2d>(0, Eigen::Vector2d(0.0, 0.0));
    rotationTheta = 0;
    translation = Eigen::Vector2d(0,0);
    V_t = Eigen::Vector2d(0,0);
    V_omega = 0;
    updateRigidBodyVars();
    updateVerts();
    mass = rho*area;
    vel_func = [](double t)->Eigen::Vector2d{ return Eigen::Vector2d(0, 0); };
}

RigidBody::RigidBody(const std::vector<Eigen::Vector2d>& in_verts, const std::vector<Eigen::Vector2i>& in_faces) : Mesh(in_verts, in_faces){
    rho = INFINITY;
    verts_no_transform = in_verts;
    rotationTheta = 0;
    translation = Eigen::Vector2d(0,0);
    V_t = Eigen::Vector2d(0,0);
    V_omega = 0;
    updateRigidBodyVars();
    updateVerts();
    mass = rho*area;
    vel_func = [](double t)->Eigen::Vector2d{ return Eigen::Vector2d(0, 0); };
}

RigidBody::RigidBody(const std::vector<Eigen::Vector2d>& in_verts, const std::vector<Eigen::Vector2i>& in_faces, const std::vector<Eigen::Vector2d>& in_vels) : Mesh(in_verts, in_faces, in_vels){
    rho = INFINITY;
    verts_no_transform = in_verts;
    rotationTheta = 0;
    translation = Eigen::Vector2d(0,0);
    V_t = Eigen::Vector2d(0,0);
    V_omega = 0;
    updateRigidBodyVars();
    updateVerts();
    mass = rho*area;
    vel_func = [](double t)->Eigen::Vector2d{ return Eigen::Vector2d(0, 0); };
}

void RigidBody::setVelFunc(std::function<Eigen::Vector2d(double)> func){ vel_func = func; }

void RigidBody::advectFE(double dt){
    setTranslation(translation + V_t*dt);
    setRotation(rotationTheta +  V_omega*dt);
    updateVerts();
    //calculateCOM();
}

void RigidBody::advectFE(double curr_t, double dt){
    setTranslation(translation + vel_func(curr_t)*dt);
    updateVerts();
    //calculateCOM();
}

void RigidBody::updateRigidBodyVars(){
    calculateArea();
    calculateCOM();
    calculateMOI();
}

void RigidBody::calculateArea(){
    double sum = 0;
    for (size_t i = 0; i<faces.size(); i++){
        Eigen::Vector2d v_a = verts_no_transform[verts_from_face(i)[0]];
        Eigen::Vector2d v_b = verts_no_transform[verts_from_face(i)[1]];
        sum += v_a.x()*v_b.y() - v_b.x()*v_a.y();
    }
    area = abs(0.5*sum);
}

void RigidBody::calculateCOM(){
    double sum_x = 0;
    double sum_y = 0;
    for (size_t i = 0; i<faces.size(); i++){
        Eigen::Vector2d v_a = verts_no_transform[verts_from_face(i)[0]];
        Eigen::Vector2d v_b = verts_no_transform[verts_from_face(i)[1]];
        sum_x += (v_a.x()+v_b.x()) * (v_a.x()*v_b.y() - v_b.x()*v_a.y());
        sum_y += (v_a.y()+v_b.y()) * (v_a.x()*v_b.y() - v_b.x()*v_a.y());
    }
    com = Eigen::Vector2d( (1/(6*area))*sum_x, (1/(6*area))*sum_y );
}

void RigidBody::calculateMOI(){
    double I_x = 0;
    double I_y = 0;
    for (size_t i = 0; i<faces.size(); i++){
        Eigen::Vector2d v_a = verts_no_transform[verts_from_face(i)[0]];
        Eigen::Vector2d v_b = verts_no_transform[verts_from_face(i)[1]];

        // translate such that COM is at origin
        v_a -= com;
        v_b -= com;

        I_x += (pow(v_a.y(),2) + v_a.y()*v_b.y() + pow(v_b.y(),2)) * (v_a.x()*v_b.y() - v_b.x()*v_a.y());
        I_y += (pow(v_a.x(),2) + v_a.x()*v_b.x() + pow(v_b.x(),2)) * (v_a.x()*v_b.y() - v_b.x()*v_a.y());
    }
    I_x /= 12;
    I_y /= 12;
    moi = rho * abs(I_x + I_y);
}

void RigidBody::updateVerts(){
    if (verts.size()!=verts_no_transform.size()){
        verts.resize(verts_no_transform.size());
    }
    for (size_t i = 0; i<verts.size(); i++){
        verts[i] = Eigen::Rotation2Dd(rotationTheta)*(verts_no_transform[i]-com) + com + translation;
    }
}

void RigidBody::updatePerVertexVels(){
    // sets the per-vertex velocities (in vels) based on the translational and rotational velocities V
    
    assert(verts.size() == vels.size());
    for(size_t i=0; i<verts.size(); i++){
        Eigen::Vector2d r_i = verts[i] - com;
        vels[i] = V_t + V_omega*Eigen::Vector2d(-r_i.y(),r_i.x());
    }
    
    // another way to do this is to construct J^T matrix
    /*
    assert(verts.size() == vels.size());
    Eigen::MatrixX3d J_transpose_n(verts.size(),3);
    Eigen::MatrixX3d J_transpose_t(verts.size(),3);
    for(size_t i=0; i<verts.size(); i++){
        double a_i = vert_area(i);
        Eigen::Vector2d r_i = verts[i] - com;

        Eigen::Vector2d n_i = calc_vertex_normal(i);
        J_transpose_n.row(i) = a_i * Eigen::Vector3d(n_i.x(), n_i.y(), Sim::cross2d(verts[i] - com, n_i));

        Eigen::Vector2d t_i = calc_vertex_tangent(i);
        J_transpose_t.row(i) = a_i * Eigen::Vector3d(t_i.x(), t_i.y(), Sim::cross2d(verts[i] - com, t_i));
    }
    Eigen::VectorXd u_n = J_transpose_n * retrieveRigidBodyV();
    Eigen::VectorXd u_t = J_transpose_t * retrieveRigidBodyV();
    for(size_t i=0; i<verts.size(); i++){
        Eigen::Vector2d n_i = calc_vertex_normal(i);
        Eigen::Vector2d t_i = calc_vertex_tangent(i);
        vels[i] = u_n(i)*n_i + u_t(i)*t_i;
    }
    */
       
} 

// essentially the same as the one in SolidMesh, but with one little tag...
// TODO: make this better later? either restructure or consolidate stuff
void RigidBody::collideAndSnap(LiquidMesh& l){
    // for each vertex in the given LiquidMesh
    for(size_t j=0; j<l.verts.size(); j++){
        Eigen::Vector2d curr_pt = l.verts[j];
        bool pt_snapped_to_solid_corner = false;

        double min_d;
        Eigen::Vector2d nearest_pt;
        int nearest_face;

        // iterate through faces of rigid body
        for (size_t i = 0; i<faces.size(); i++){

            // retrieve endpoints of current face, these should be oriented
            int ptA_ind = faces[i][0];
            int ptB_ind = faces[i][1];

            Eigen::Vector2d ptA = verts[ptA_ind];
            Eigen::Vector2d ptB = verts[ptB_ind];

            Eigen::Vector2d ptA_vel = vels[ptA_ind];
            Eigen::Vector2d ptB_vel = vels[ptB_ind];

            // first check to snap to 'sharp' features in geometry, if so, we can skip the rest of the loop
            if ((curr_pt-ptA).norm()<epsilon && solid_angle_is_acute(ptA_ind)){
                min_d = (curr_pt-ptA).norm();
                nearest_pt = ptA;
                pt_snapped_to_solid_corner = true;
                break;
            } else if ((curr_pt-ptB).norm()<epsilon && solid_angle_is_acute(ptB_ind)){
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
                    // snap to solid corner here too probably
                }
            
                // logic to get minimal distance to a segment in the SolidMesh
                if (i==0 || (proj_pt-curr_pt).norm()<=abs(min_d)){
                    min_d = (proj_pt-curr_pt).norm();
                    nearest_pt = proj_pt;
                }
            }
        }

        // Collision detected!!
        // since winding number when precisely on the segment is a bit misbehaved
        // we define 'collision' as when min_d is less than our epsilon OR when the winding number is nonzero
        // as a nonzero winding number indicates that we are intersecting the mesh
        if(min_d <= epsilon || windingNumber(curr_pt)>1.0e-8){
            // snap the point to the nearest face
            l.verts[j] = nearest_pt;
            // set effective velocity
            l.per_vertex_rb_contact[j] = rb_sim_id;
            l.vels_solid[j] = V_t + V_omega*Eigen::Vector2d(-(l.verts[j]-com).y(),(l.verts[j]-com).x());
            l.set_boundaries_for_vertex(j, false, true, false, pt_snapped_to_solid_corner);
        }
    }
}

bool RigidBody::loadMeshFromFile(RigidBody &m, std::string infileName){
    // load rigid body mesh file
    std::ifstream infile(infileName);
    if (!infile.is_open()) { 
        std::cerr << "Unable to open solid file "<<infileName<<"!" << std::endl; 
        assert(!"Unable to open solid mesh file!");
    }

    std::vector<Eigen::Vector2d> v;
    std::vector<Eigen::Vector2i> f;

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
        } else if (linetype.compare("theta") == 0){
            ss >> m.rotationTheta;
        } else if (linetype.compare("translation") == 0){
            double a,b;
            ss >> a;
            ss >> b;
            m.translation = Eigen::Vector2d(a,b);
        } else if (linetype.compare("u_rotational") == 0){
            ss >> m.V_omega;
        } else if (linetype.compare("u_translational") == 0){
            double a,b;
            ss >> a;
            ss >> b;
            m.V_t = Eigen::Vector2d(a,b);
        } else if (linetype.compare("rho") == 0){
            ss >> m.rho;
        } else if (linetype.compare("clamp_translation_x") == 0){
            ss >> m.clamp_translation_x;
        }else if (linetype.compare("clamp_translation_y") == 0){
            ss >> m.clamp_translation_y;
        }else if (linetype.compare("clamp_rotation") == 0){
            ss >> m.clamp_rotation;
        }
        else {
            std::cerr << "Invalid line in "<<infileName<<"! Skipping line..." << std::endl;
        }
    }
    infile.close();
    
    assert(v.size() == f.size());

    m.verts_no_transform = v;
    m.verts = v;
    m.faces = f;
    m.update_neighbor_face_vecs();

    m.vels.resize(v.size());
    
    m.updateRigidBodyVars();
    m.updateVerts();
    m.updatePerVertexVels();

    m.mass = m.rho*m.area;

    return true;
}