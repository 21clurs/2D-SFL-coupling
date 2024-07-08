#include "rigidbody.h"

#include <fstream>
#include "sim.h"

RigidBody::RigidBody() : SolidMesh(){
    verts_no_transform = std::vector<Eigen::Vector2d>(0, Eigen::Vector2d(0.0, 0.0));
    rotationTheta = 0;
    translation = Eigen::Vector2d(0,0);
    V_t = Eigen::Vector2d(0,0);
    V_omega = 0;
    updateRigidBodyVars();
    updateVerts();
}

RigidBody::RigidBody(const std::vector<Eigen::Vector2d>& in_verts, const std::vector<Eigen::Vector2i>& in_faces) : SolidMesh(in_verts, in_faces){
    verts_no_transform = in_verts;
    rotationTheta = 0;
    translation = Eigen::Vector2d(0,0);
    V_t = Eigen::Vector2d(0,0);
    V_omega = 0;
    updateRigidBodyVars();
    updateVerts();
}

RigidBody::RigidBody(const std::vector<Eigen::Vector2d>& in_verts, const std::vector<Eigen::Vector2i>& in_faces, const std::vector<Eigen::Vector2d>& in_vels) : SolidMesh(in_verts, in_faces, in_vels){
    verts_no_transform = in_verts;
    rotationTheta = 0;
    translation = Eigen::Vector2d(0,0);
    V_t = Eigen::Vector2d(0,0);
    V_omega = 0;
    updateRigidBodyVars();
    updateVerts();
}

void RigidBody::setRotation(double theta){
    rotationTheta = theta;
    updateRotationMat();
}

void RigidBody::advectFE(double dt){
    setTranslation(translation + V_t*dt);
    setRotation(rotationTheta +  V_omega*dt);
    updateVerts();
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
        Eigen::Vector2d v_a = verts_no_transform[verts_from_face(i)[0]];
        Eigen::Vector2d v_b = verts_no_transform[verts_from_face(i)[1]];
        sum += v_a.x()*v_b.y() - v_b.x()*v_a.y();
    }
    area = 0.5*sum;
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
    moi = I_x + I_y; // placeholder
}

void RigidBody::updateRotationMat(){
    // rotates clockwise by theta
    rotationMat << cos(rotationTheta), -sin(rotationTheta),
                sin(rotationTheta), cos(rotationTheta);
}

void RigidBody::updateVerts(){
    if (verts.size()!=verts_no_transform.size()){
        verts.resize(verts_no_transform.size());
    }
    for (size_t i = 0; i<verts.size(); i++){
        verts[i] = rotationMat*(verts_no_transform[i]-com) + com + translation;
    }
}

void RigidBody::updatePerVertexVels(){
    // sets the per-vertex velocities (in vels) based on the translational and rotational velocities V
    assert(verts.size() == vels.size());
    for(size_t i=0; i<verts.size(); i++){
        Eigen::Vector2d r_i = verts[i] - com;
        vels[i] = V_t + V_omega*Eigen::Vector2d(-r_i.y(),r_i.x());
    }
    /*
    // another way to do this is to construct J^T matrix
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

bool RigidBody::loadMeshFromFile(RigidBody &m, std::string infileName){
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
        } else if (linetype.compare("theta") == 0){
            ss >> m.rotationTheta;
        } else if (linetype.compare("translation") == 0){
            int a,b;
            ss >> a;
            ss >> b;
            m.translation = Eigen::Vector2d(a,b);
        } else if (linetype.compare("u_rotational") == 0){
            ss >> m.V_omega;
        } else if (linetype.compare("u_translational") == 0){
            int a,b;
            ss >> a;
            ss >> b;
            m.V_t = Eigen::Vector2d(a,b);
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

    return true;
}