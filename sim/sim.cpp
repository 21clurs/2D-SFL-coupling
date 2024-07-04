#include "sim.h"

#include <iostream>
#include <fstream> 
#include <stdexcept>
#include "boundaryintegral.h"
#include "simoptions.h"
#include "scenes.h"

using Eigen::MatrixXd;
using Eigen::Vector2d;
using Eigen::Vector2i;

// adapted from Options::ParseOptionFile() from Da 2016 code
bool Sim::setAndLoadSimOptions(std::string infileName){
    // set a bunch of default sim options
    SimOptions::addStringOption ("scene", "circle");
    SimOptions::addDoubleOption ("time-step", 0.01);
    SimOptions::addDoubleOption ("simulation-time", 10.0);
    SimOptions::addDoubleOption ("gravity", 0);

    SimOptions::addDoubleOption ("sigma-sl", 1);
    SimOptions::addDoubleOption ("sigma-la", 1);
    SimOptions::addDoubleOption ("sigma-sa", 1);
    SimOptions::addDoubleOption ("rho", 1);

    SimOptions::addIntegerOption ("mesh-size-n", 32);
    SimOptions::addIntegerOption ("mesh-remesh-on", 0);
    SimOptions::addIntegerOption ("mesh-remesh-iters", 0);
    SimOptions::addDoubleOption ("mesh-edge-max-ratio", 1.1);
    SimOptions::addDoubleOption ("mesh-edge-min-ratio", 0.9);

    SimOptions::addDoubleOption ("mesh-collision-epsilon", 0.9);

    SimOptions::addStringOption ("mesh-initial-velocity", "zero");

    SimOptions::addBooleanOption ("marker-particles-on", false);

    // scene/original liquid mesh shape specific parameters
    SimOptions::addDoubleOption ("radius", 1); // circle
    SimOptions::addDoubleOption ("width", 1); // rectangle
    SimOptions::addDoubleOption ("height", 1); // rectangle
    
    // load sim options file
    SimOptions::loadSimOptions(infileName);
}

void Sim::run(){
    // set up scene
    Scenes::scene(this, SimOptions::strValue("scene"), SimOptions::strValue("mesh-initial-velocity"));
    // collide liquid mesh with all the solids and such
    collide();
    // generate marker particles after collision
    if (SimOptions::boolValue("marker-particles-on")){
        genMarkerParticles(-1.5, 1.5, -0.5, 0, 0.1);
    }

    double dt = SimOptions::doubleValue("time-step");
    int frames = (int) SimOptions::doubleValue("simulation-time")/dt;
    for (int i=0; i<frames; i++){
        // sim stuff
        outputFrame(std::to_string(i)+".txt");
        step_sim(i*dt);
        
        // progress messages
        std::cout<<"Simulation steps "<<i+1<<"/"<<frames<<" complete."<<"\r";
        std::cout.flush();

        if (i==30){
        //testSolid.setVelFunc([](double t)->Vector2d{ return Vector2d(0, -(1.0/4.0)*sin(t*4.0)); });
        } 
        if (i==344){
        //testSolid.setVelFunc([](double t)->Vector2d{ return Vector2d(0, 0); });
        }
    }
    std::cout << std::endl; 
}

Sim::Sim(){
    n = SimOptions::intValue("mesh-size-n");
    dt = SimOptions::doubleValue("time-step");
    m = LiquidMesh();
    sigma = SimOptions::doubleValue("sigma-la");
    sigma_SL = SimOptions::doubleValue("sigma-sl");
    sigma_SA = SimOptions::doubleValue("sigma-sa");
    rho = SimOptions::doubleValue("rho");
    gravity = Eigen::Vector2d({0.0, SimOptions::doubleValue("gravity")});
    
    p = Eigen::VectorXd::Zero(n);
    dpdn = Eigen::VectorXd::Zero(n);
    markerparticles = {};
}

Sim::Sim(LiquidMesh& m, int n, float dt):
    n(n),
    dt(dt),
    m(m),
    sigma(0.5),
    sigma_SL(1.0),
    sigma_SA(1.0),
    rho(1.0),
    gravity(Eigen::Vector2d({0.0, -5.0}))
{
    p = Eigen::VectorXd::Zero(n);
    dpdn = Eigen::VectorXd::Zero(n);
    markerparticles = {};
}

void Sim::addSolid(SolidMesh* solid){
    solids.emplace_back(solid);
}

bool Sim::outputFrame(std::string filename, std::string filelocation){
    std::ofstream file(filelocation+filename);
    for (size_t i=0; i<m.verts.size(); i++){
        char vtype;
        if (m.is_corner[i])
            vtype = 'c';
        else if (m.is_air[i])
            vtype = 'a';
        else if (m.is_solid[i])
            vtype = 's';
        else
            vtype = 't';
        file<<"v "<<m.verts[i][0]<<" "<<m.verts[i][1]<<" "<<vtype<<std::endl;
    }
    file<<std::endl;

    // outward vertex norms
    for (size_t i=0; i<m.verts.size(); i++)
        file<<"vn "<<(m.calc_vertex_normal(i))[0]<<" "<<(m.calc_vertex_normal(i))[1]<<std::endl;
    file<<std::endl;
    
    // vertex tangents -- should be oriented as 90deg anticlockwise rotation from the outward norms
    for (size_t i=0; i<m.faces.size(); i++)
        file<<"vt "<<(m.calc_vertex_tangent(i))[0]<<" "<<(m.calc_vertex_tangent(i))[1]<<std::endl;
    file<<std::endl;

    // vels
    for (size_t i=0; i<m.faces.size(); i++)
        file<<"u "<<m.vels[i][0]<<" "<<m.vels[i][1]<<std::endl;
    file<<std::endl;

    // faces
    for (size_t i=0; i<m.faces.size(); i++)
        file<<"f "<<m.faces[i][0]<<" "<<m.faces[i][1]<<std::endl;
    file<<std::endl;

    // solids
    size_t count = 0;
    for (size_t i=0; i<solids.size(); i++){
        // solid verts
        for (size_t j=0; j<solids[i]->verts.size(); j++){
            file<<"v "<<solids[i]->verts[j][0]<<" "<<solids[i]->verts[j][1]<<std::endl;
        }
        // solid faces
        for (size_t j=0; j<solids[i]->faces.size(); j++){
            file<<"f "<< m.faces.size()+ solids[i]->faces[j][0] + count<<" "<< m.faces.size()+ solids[i]->faces[j][1] + count<<std::endl;
        }
        count += solids[i]->verts.size();
    }
    file<<std::endl;

    // marker particles
    for (size_t i=0; i<markerparticles.size(); i++){
        file<<"p "<<markerparticles[i][0]<<" "<<markerparticles[i][1]<<std::endl;
    }
    file<<std::endl;

    // marker particle velocities
    for (size_t i=0; i<markerparticles.size(); i++){
        Eigen::Vector2d tmp = HHD_FD(markerparticles[i], 0.01);
        file<<"pv "<<tmp[0]<<" "<<tmp[1]<<std::endl;
    }
    file.close();
    

    return true;
}

void Sim::step_sim(double curr_t){
    /*std::cout<<"pre remesh vec"<<std::endl;
    for (size_t i=0; i<m.verts.size(); i++){
        std::cout<<m.verts[i][0]<<", "<<m.verts[i][1]<<std::endl;
    }*/
    /*
    std::cout<<"is_air: ";
    for (size_t i=0; i<m.verts.size(); i++){
        std::cout<<m.is_air[i]<<", ";
    }
    std::cout<<std::endl;

    std::cout<<"is_solid: ";
    for (size_t i=0; i<m.verts.size(); i++){
        std::cout<<m.is_solid[i]<<", ";
    }
    std::cout<<std::endl;

    std::cout<<"is_triple: ";
    for (size_t i=0; i<m.verts.size(); i++){
        std::cout<<m.is_triple[i]<<", ";
    }
    std::cout<<std::endl;*/
    step_advect(curr_t);

    //step_solidinfluxreverse();

    collide();
    remesh();
    //std::cout<<"post remesh vec "<<m.verts.size()<<std::endl;
    /*for (size_t i=0; i<m.verts.size(); i++){
        std::cout<<"("<<m.verts[i][0]<<", "<<m.verts[i][1]<<")";
    }
    std::cout<<std::endl;*/
    //collide(); // remesh will collide things

    //outputFrame(std::to_string(1)+".txt");

    //std::cout<<"post collide vec "<<m.verts.size()<<std::endl;
    /*for (size_t i=0; i<m.verts.size(); i++){
        std::cout<<"("<<m.verts[i][0]<<", "<<m.verts[i][1]<<")";
    }
    std::cout<<std::endl;*/
    step_HHD();
    //std::cout<<"post HHD vels "<<m.vels.size()<<std::endl;
    /*for (size_t i=0; i<m.vels.size(); i++){
        std::cout<<m.vels[i][0]<<", "<<m.vels[i][1]<<std::endl;
    }*/
    step_gravity(); // body forces only applied on div-free velocity field.

    step_BEM();
    //std::cout<<"post BEM vels "<<m.vels.size()<<std::endl;
    /*for (size_t i=0; i<m.vels.size(); i++){
        std::cout<<m.vels[i][0]<<", "<<m.vels[i][1]<<std::endl;
    }*/

    //step_solidinflux();
}

/*
void Sim::step_solidinflux(){
    std::vector<Eigen::Vector2d> nSolid(m.verts.size(), Eigen::Vector2d::Zero());
    std::vector<bool> solid_faces = m.get_solid_faces();
    for (size_t i=0; i<m.faces.size(); i++){
        if (solid_faces[i] == true){
            Eigen::Vector2d n = m.calc_face_normal(i);
            Eigen::Vector2i vCurrFace = m.verts_from_face(i);
            nSolid[vCurrFace[0]] = (nSolid[vCurrFace[0]]+n).normalized();
            nSolid[vCurrFace[1]] = (nSolid[vCurrFace[1]]+n).normalized();
        }
    }
    for (size_t i=0; i<m.verts.size(); i++){
        if (m.is_solid[i]){ //what about triple point?
            Eigen::Vector2d n = nSolid[i];
            assert(n.squaredNorm() != 0);
            m.vels[i] -= m.vels[i].dot(n)*n;
            m.vels[i] += m.vels_solid[i];
        }
    }
}
*/

void Sim::step_advect(double t){
    // advect the liquid mesh
    for (size_t i=0; i<n; i++){
        m.verts[i] = m.verts[i] + m.vels[i]*dt;
    }
    // advect the scripted solids
    for (size_t i=0; i<solids.size(); i++){
        solids[i]->advectFE(t-dt*30, dt);
    }
    // "advect" the marker particles
    for (size_t i=0; i<markerparticles.size(); i++){
        markerparticles[i] = markerparticles[i] + HHD_FD(markerparticles[i],0.0001)*dt;
    }
}

/*
void Sim::step_solidinfluxreverse(){
    std::vector<Eigen::Vector2d> nSolid(m.verts.size(), Eigen::Vector2d::Zero());
    std::vector<bool> solid_faces = m.get_solid_faces();
    for (size_t i=0; i<m.faces.size(); i++){
        if (solid_faces[i] == true){
            Eigen::Vector2d n = m.calc_face_normal(i);
            Eigen::Vector2i vCurrFace = m.verts_from_face(i);
            nSolid[vCurrFace[0]] = (nSolid[vCurrFace[0]]+n).normalized();
            nSolid[vCurrFace[1]] = (nSolid[vCurrFace[1]]+n).normalized();
        }
    }
    for (size_t i=0; i<m.verts.size(); i++){
        if (m.is_solid[i]){ //what about triple point?
            Eigen::Vector2d n = nSolid[i];
            assert(n.squaredNorm() != 0);
            m.vels[i] = m.vels[i] - m.vels[i].dot(n)*n;// -  gravity*n ;
        }
    }
}
*/

void Sim::step_HHD(){
    const std::vector<Eigen::Vector2d> quadrature_GQ = BoundaryIntegral::gaussian_quadrature();
    const std::vector<Eigen::Vector2d> quadrature_DE = BoundaryIntegral::tanh_sinh_quadrature();

    const double negOneOver2pi = -1.0/(2.0*M_PI);

    std::vector<double> dPhidn(n, 0.0);
    std::vector<double> curlA(n, 0.0);

    for (size_t i=0; i<m.verts.size(); i++){
        int xi_a, xi_b;
        // iterate over each adjacent edge (only 2 in this case)
        for (size_t p=0; p<2; p++){
            int f_i = m.faces_from_vert(i)[p];

            Eigen::Vector2d n_x = m.calc_face_normal(f_i);
            Eigen::Vector2d t_x = m.calc_face_tangent(f_i);
            
            // shuffling around the endpoint indices
            // such that xi_a is the 'current' vertex
            if (m.verts_from_face(f_i)[0] == i){
                xi_a = m.verts_from_face(f_i)[0];
                xi_b = m.verts_from_face(f_i)[1];
            } else {
                xi_a = m.verts_from_face(f_i)[1];
                xi_b = m.verts_from_face(f_i)[0];
            }

            double II_gradPhi = 0;
            double charge = 0;
            double II_curlA = 0;

            double jacobian = 1.0/2.0; // a little sketchy here but ok lets see...

            // iterate over quadrature points on the adjacent faces
            for (size_t qi=0; qi<quadrature_DE.size(); qi++){
                double qik = quadrature_DE[qi].x();
                double qiw = quadrature_DE[qi].y();

                double theta_i = (1-qik)/2;

                Eigen::Vector2d x = lin_interp(m.verts[xi_a], m.verts[xi_b], qik);
                Eigen::Vector2d v_x = lin_interp(m.vels[xi_a], m.vels[xi_b], qik);

                double I_gradPhi = 0;
                double I_curlA = 0;

                // iterating over all faces (surface integral over the entire boundary)
                for (size_t f_j=0; f_j<m.faces.size(); f_j++){
                    Eigen::Vector2d n_y = m.calc_face_normal(f_j);

                    double f_len = m.face_length(f_j);

                    if (f_i == f_j){ //singularity case 
                        // gradPhi
                        // contribution is 0 for this case

                        // curlA
                        // split around the point
                        for (size_t k=0; k<2; k++){
                            double val_curlA = 0;
                            double len = (x-m.verts[m.verts_from_face(f_j)[k]]).norm();
                            for (size_t qj=0; qj<quadrature_DE.size(); qj++){
                                double qjk = quadrature_DE[qj].x();
                                double qjw = quadrature_DE[qj].y();
                                
                                Eigen::Vector2d y = lin_interp(x, m.verts[m.verts_from_face(f_j)[k]], qjk);
                                Eigen::Vector2d v_y = lin_interp(v_x, m.vels[m.verts_from_face(f_j)[k]], qjk);

                                if (p==0)
                                    val_curlA += jacobian * qjw * cross2d(n_y,v_y) * BoundaryIntegral::G(x,y) ;
                                else if (p==1)
                                    val_curlA += -1 * jacobian * qjw * cross2d(n_y,v_y) * BoundaryIntegral::G(x,y) ; // weirdsign flip
                            }
                            I_curlA += val_curlA * len;
                        }
                    } else{
                        // gradPhi
                        double val_gradPhi = 0;
                        for (size_t qj=0; qj<quadrature_GQ.size(); qj++){
                            double qjk = quadrature_GQ[qj].x();
                            double qjw = quadrature_GQ[qj].y();
                            
                            Eigen::Vector2d y = lin_interp(m.verts[m.verts_from_face(f_j)[0]], m.verts[m.verts_from_face(f_j)[1]], qjk);
                            Eigen::Vector2d v_y = lin_interp(m.vels[m.verts_from_face(f_j)[0]], m.vels[m.verts_from_face(f_j)[1]], qjk);

                            Eigen::Vector2d dx = x-y;
                            double dxn = dx.norm();
                            double dGdnx = negOneOver2pi * dx.dot(n_x) / (dxn*dxn);

                            val_gradPhi += jacobian * qjw * v_y.dot(n_y) * dGdnx ;
                        }
                        I_gradPhi += val_gradPhi * f_len;

                        // curlA
                        double val_curlA = 0;
                        for (size_t qj=0; qj<quadrature_DE.size(); qj++){
                            double qjk = quadrature_DE[qj].x();
                            double qjw = quadrature_DE[qj].y();

                            Eigen::Vector2d y = lin_interp(m.verts[m.verts_from_face(f_j)[0]], m.verts[m.verts_from_face(f_j)[1]], qjk);
                            Eigen::Vector2d v_y = lin_interp(m.vels[m.verts_from_face(f_j)[0]], m.vels[m.verts_from_face(f_j)[1]], qjk);

                            Eigen::Vector2d dx = x-y;
                            double dxn = dx.norm();
                            double dGdtx = negOneOver2pi * dx.dot(t_x) / (dxn*dxn);

                            val_curlA += jacobian * qjw * cross2d(n_y,v_y) * dGdtx ;
                        }
                        I_curlA += theta_i * val_curlA * f_len;
                    }
                }
                II_gradPhi += theta_i * jacobian * qiw * I_gradPhi;
                charge += theta_i * jacobian * qiw * (-negOneOver2pi)*m.solid_angle(i)*v_x.dot(n_x);
                II_curlA += jacobian * qiw * I_curlA;

                //std::cout<<"I_gradPhi_"<<f_i<<"_"<<qik<<": "<<I_gradPhi<<", "<<(-negOneOver2pi)*m.solid_angle(i)*v_x.dot(n_x)<<std::endl;
                //std::cout<<v_x.dot(n_x)<<std::endl;
            }
            dPhidn[i] += (II_gradPhi + charge)*m.face_length(f_i);
            curlA[i] += II_curlA*m.face_length(f_i);
        }
    }
    /*
    std::cout<<"dPhidn: "<<std::endl;
    for (size_t i=0;i<dPhidn.size();i++){
        std::cout<<dPhidn[i]<<",";
    }
    std::cout<<std::endl;
    std::cout<<"curlA: "<<std::endl;
    for (size_t i=0;i<curlA.size();i++){
        std::cout<<curlA[i]<<",";
    }
    std::cout<<std::endl;
    */

    // constructing tangential projection
    std::vector<Eigen::Vector2d> V_t = std::vector<Eigen::Vector2d>(n,Eigen::Vector2d(0.0,0.0));
    std::vector<Eigen::Vector2d> V_n = std::vector<Eigen::Vector2d>(n,Eigen::Vector2d(0.0,0.0));
    for (size_t i=0; i<m.verts.size(); i++){

        Eigen::Vector2d n_i = m.calc_vertex_normal(i);

        Eigen::Matrix2d P = Eigen::Matrix2d::Identity()-n_i*n_i.transpose();
        V_t[i] = P*m.vels[i];

        dPhidn[i] /= m.vert_area(i);
        curlA[i] /= m.vert_area(i);

        V_n[i] = (dPhidn[i] + curlA[i])*n_i;
    }
    for (size_t i=0; i<m.verts.size(); i++){
        m.vels[i] = V_t[i] + V_n[i];
    }
}

void Sim::step_gravity(){
    for (size_t i=0; i<n; i++){
        m.vels[i] += gravity*dt;
    }
}

void Sim::step_BEM(){
    int N = m.verts.size();
    Eigen::VectorXd BC_p = Eigen::VectorXd::Zero(N);
    Eigen::VectorXd BC_dpdn = Eigen::VectorXd::Zero(N);
    step_BEM_BC(BC_p, BC_dpdn);
    
    // reset p and dpdn, and make correct size
    p = Eigen::VectorXd::Zero(N);
    dpdn = Eigen::VectorXd::Zero(N);

    step_BEM_solve(BC_p, BC_dpdn, p, dpdn);
    
    step_BEM_gradP(BC_p, BC_dpdn, p,dpdn);
}

void Sim::step_BEM_BC(Eigen::VectorXd& BC_p, Eigen::VectorXd& BC_dpdn){
    // Setting BCs
    int N = m.verts.size();
    std::vector<bool> face_is_solid = m.get_solid_faces();

    double triple_junction_virtual_width = 0.25 * m.calc_avg_face_length();

    //Eigen::VectorXd BC_p = Eigen::VectorXd::Zero(N);
    //Eigen::VectorXd BC_dpdn = Eigen::VectorXd::Zero(N);

    for (size_t i=0; i<N; i++){
        Eigen::Vector2i incident_faces = m.faces_from_vert(i);
        // 'air' vertex: if it is explicitly air vertex or neither incident face is a solid face
        if (m.is_air[i] || (!face_is_solid[incident_faces[0]] && !face_is_solid[incident_faces[1]])){
            double H_i = m.signed_mean_curvature(i);
            
            BC_p[i] = sigma * H_i * dt;     // known
            BC_dpdn[i] = 0;                 // unknown
        }
        // 'solid' vertex: if it is explicity solid vertex or both incident faces are solid faces
        else if (m.is_solid[i] || (face_is_solid[incident_faces[0]] && face_is_solid[incident_faces[1]])){
            Eigen::Vector2d n_i = m.calc_vertex_normal(i);

            BC_p[i] = 0;
            //BC_dpdn[i] = n_i.dot(m.vels[i] + gravity*dt - m.vels_solid[i]) * rho;
            BC_dpdn[i] = n_i.dot(m.vels[i] - m.vels_solid[i]) * rho;
        }
        // triple junction
        else{
            Eigen::Vector2d n_solid_outward;
            int solid_vert, air_vert;
            if (face_is_solid[incident_faces[0]]){
                n_solid_outward = m.calc_face_normal(incident_faces[0]);
                
                // TODO: make this nicer somehow? messy...
                Eigen::Vector2i solid_face_verts = m.verts_from_face(incident_faces[0]);
                solid_vert = (i == solid_face_verts[0])? solid_face_verts[1]: solid_face_verts[0];
                Eigen::Vector2i air_face_verts = m.verts_from_face(incident_faces[1]);
                air_vert = (i == air_face_verts[0])? air_face_verts[1]: air_face_verts[0];
            } else if (face_is_solid[incident_faces[1]]){
                n_solid_outward = m.calc_face_normal(incident_faces[1]);

                // TODO: make this nicer somehow? messy...
                Eigen::Vector2i solid_face_verts = m.verts_from_face(incident_faces[1]);
                solid_vert = (i == solid_face_verts[0])? solid_face_verts[1]: solid_face_verts[0];
                Eigen::Vector2i air_face_verts = m.verts_from_face(incident_faces[0]);
                air_vert = (i == air_face_verts[0])? air_face_verts[1]: air_face_verts[0];
            } else
                throw std::logic_error("This should not have been reached-- triple junction should fall into one of the above statements");

            Eigen::Vector2d t_solid_outward = m.verts[i] - m.verts[solid_vert];
            t_solid_outward.normalize();
            Eigen::Vector2d t_air_outward = m.verts[i] - m.verts[air_vert];
            t_air_outward.normalize();
            /*
            std::cout<<"t_solid_outward at "<<i<<": "<<t_solid_outward<<std::endl;
            std::cout<<"t_air_outward at "<<i<<": "<<t_air_outward<<std::endl;
            */
            Eigen::Vector2d st_force_SL = -t_solid_outward * sigma * sigma_SL;
            Eigen::Vector2d st_force_SA = t_solid_outward * sigma * sigma_SA;
            Eigen::Vector2d st_force_LA = -t_air_outward * sigma;

            Eigen::Vector2d st_force_combined = st_force_SL + st_force_SA + st_force_LA;
            /*
            std::cout<<"st_force_combined: "<<st_force_combined<<std::endl;
            std::cout<<"triple_junction_virtual_width: "<<triple_junction_virtual_width<<std::endl;
            */
            double pressure_jump_normal_to_triple_junction = st_force_combined.dot(-t_solid_outward) / triple_junction_virtual_width;

            BC_p[i] = pressure_jump_normal_to_triple_junction * dt;
            //BC_dpdn[i] = n_solid_outward.dot(m.vels[i] + gravity*dt - m.vels_solid[i]) * rho; // Note: technically, this is dpdn_s only, not dpdn
            BC_dpdn[i] = n_solid_outward.dot(m.vels[i] - m.vels_solid[i]) * rho; 
        }
    }

    assert(BC_p == BC_p);
    assert(BC_dpdn == BC_dpdn);
        
    //std::cout<<"BC_p: "<<BC_p<<std::endl;
    //std::cout<<"BC_dpdn: "<<BC_dpdn<<std::endl;
}

void Sim::step_BEM_solve(Eigen::VectorXd& BC_p, Eigen::VectorXd& BC_dpdn, Eigen::VectorXd& p, Eigen::VectorXd& dpdn){
    // BEM
    int N = m.verts.size();
    std::vector<bool> face_is_solid = m.get_solid_faces();

    /*
    std::cout<<"x:";
    for (size_t i=0; i<N; i++){
        std::cout<<m.verts[i]<<std::endl;
        std::cout<<"Solid: "<<m.is_solid[i]<<", Air: "<<m.is_air[i]<<", Triple: "<<m.is_triple[i]<<std::endl;
    }
    std::cout<<std::endl;

    std::cout<<"p_input"<<BC_p<<std::endl;
    std::cout<<"dpdn_input"<<BC_dpdn<<std::endl;
    */

    Eigen::MatrixXd collocA_solid = Eigen::MatrixXd::Zero(N, N);
    Eigen::MatrixXd collocA_air = Eigen::MatrixXd::Zero(N, N);
    Eigen::MatrixXd collocB_solid = Eigen::MatrixXd::Zero(N, N);
    Eigen::MatrixXd collocB_air = Eigen::MatrixXd::Zero(N, N);

    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(N, N);
    Eigen::MatrixXd rhs_per_vertex = Eigen::MatrixXd::Zero(N, N);
    
    // BEM coeffs/assembly
    double negOneOver2pi = -1.0/(2.0*M_PI);

    const int quad_N = 4;
    const std::vector<Eigen::Vector2d> quad_GQ = BoundaryIntegral::gaussian_quadrature(quad_N);
    
    // iterate through vertices-- reference point
    for (size_t i=0; i<N; i++){
        int prev_i = m.prev_neighbor_index(i);
        int next_i = m.next_neighbor_index(i);
        
        double d_prev = (m.verts[i] - m.verts[prev_i]).norm();
        double d_next = (m.verts[i] - m.verts[next_i]).norm();

        Eigen::Vector2i neighbor_faces = m.faces_from_vert(i);

        // iterate through vertices-- field point
        for (size_t j=0; j<N; j++){
            Eigen::Vector2d y_j = m.verts[j];

            double G_quad_prev = 0;
            double G_quad_next = 0;
            double dGdn_quad_prev = 0;
            double dGdn_quad_next = 0;

            // 2 singular integrals
            if (i==j){
                G_quad_prev = (d_prev/2.0) * (log(d_prev)-3.0/2.0) * negOneOver2pi;
                dGdn_quad_prev = 0; // placeholder, we deal with this later w/row-sum
                G_quad_next = (d_next/2.0) * (log(d_next)-3.0/2.0) * negOneOver2pi;
                dGdn_quad_next = 0; // placeholder, we deal with this later w/row-sum
            } else{
                if (j==prev_i){
                    G_quad_prev = (d_prev/2.0) * (log(d_prev)-1.0/2.0) * negOneOver2pi;
                    dGdn_quad_prev = 0;
                } else{
                    // generic case -- no singularities
                    Eigen::Vector2d n_x = m.calc_face_normal(neighbor_faces[0]);
                    /*if(i==0 && j==2){
                        std::cout<<"PREV VERTEX:"<<std::endl;
                        std::cout<<"from: "<<m.verts[i]<<std::endl;
                        std::cout<<"to: "<<m.verts[prev_i]<<std::endl;
                    }*/
                    for(size_t k=0; k<quad_N; k++){
                        double qk = quad_GQ[k].x();
                        double qw = quad_GQ[k].y();
                        Eigen::Vector2d x_i = lin_interp(m.verts[prev_i], m.verts[i], qk);
                        double Mi = M_2(qk);
                        G_quad_prev += qw * BoundaryIntegral::G(x_i, y_j) * Mi;
                        dGdn_quad_prev += qw * BoundaryIntegral::dGdnx(x_i, y_j, n_x) * Mi;
                        /*if(i==0 && j==2){
                            std::cout<<"x_"<<k<<": "<<x_i<<std::endl;
                            std::cout<<"y_"<<k<<": "<<y_j<<std::endl;
                            std::cout<<"M_"<<k<<": "<<Mi<<std::endl;
                            std::cout<<"qk: "<<qk<<std::endl;
                            std::cout<<"qw: "<<qw<<std::endl;
                            std::cout<<"G(x,y): "<<BoundaryIntegral::G(x_i, y_j)<<std::endl;
                        }*/
                    }
                    G_quad_prev *= 0.5*d_prev;
                    dGdn_quad_prev *= 0.5*d_prev;
                    /*if(i==0 && j==2){
                        std::cout<<"i=0,j=2 G_prev: "<<G_quad_prev<<std::endl;
                        std::cout<<"i=0,j=2 dGdn_prev: "<<dGdn_quad_prev<<std::endl;
                    }*/
                }
                if (j==next_i){
                    G_quad_next = (d_next/2.0) * (log(d_next)-1.0/2.0) * negOneOver2pi;
                } else{
                    // generic case -- no singularities
                    Eigen::Vector2d n_x = m.calc_face_normal(neighbor_faces[1]);
                    /*if(i==0 && j==2){
                        std::cout<<"NEXT VERTEX:"<<std::endl;
                        std::cout<<"from: "<<m.verts[i]<<std::endl;
                        std::cout<<"to: "<<m.verts[next_i]<<std::endl;
                    }*/
                    for(size_t k=0; k<quad_N; k++){
                        double qk = quad_GQ[k].x();
                        double qw = quad_GQ[k].y();
                        Eigen::Vector2d x_i = lin_interp(m.verts[i], m.verts[next_i], qk);
                        double Mi = M_1(qk);
                        G_quad_next += qw * BoundaryIntegral::G(x_i, y_j) * Mi;
                        dGdn_quad_next += qw * BoundaryIntegral::dGdnx(x_i, y_j, n_x) * Mi;

                        /*if(i==0 && j==2){
                            std::cout<<"x_"<<k<<": "<<x_i<<std::endl;
                            std::cout<<"y_"<<k<<": "<<y_j<<std::endl;
                            std::cout<<"M_"<<k<<": "<<Mi<<std::endl;
                            std::cout<<"qk: "<<qk<<std::endl;
                            std::cout<<"qw: "<<qw<<std::endl;
                            std::cout<<"G(x,y): "<<BoundaryIntegral::G(x_i, y_j)<<std::endl;
                        }*/
                    }
                    G_quad_next *= 0.5*d_next;
                    dGdn_quad_next *= 0.5*d_next;
                    /*if(i==0 && j==2){
                        std::cout<<"i=0,j=2 G_next: "<<G_quad_next<<std::endl;
                        std::cout<<"i=0,j=2 dGdn_next: "<<G_quad_next<<std::endl;
                    }*/
                }
            }

            // organizing into collocA and collocB
            // prev face
            if (face_is_solid[neighbor_faces[0]]){
                collocA_solid(j,i) += dGdn_quad_prev;
                collocB_solid(j,i) += G_quad_prev;
            } else{
                collocA_air(j,i) += dGdn_quad_prev;
                collocB_air(j,i) += G_quad_prev;
            }
            // next face
            if (face_is_solid[neighbor_faces[1]]){
                collocA_solid(j,i) += dGdn_quad_next;
                collocB_solid(j,i) += G_quad_next;
            } else{
                collocA_air(j,i) += dGdn_quad_next;
                collocB_air(j,i) += G_quad_next;
            }
        }
    }
    Eigen::MatrixXd collocA = collocA_solid + collocA_air;

    // row-sum to find diagonal entries
    for (size_t i=0; i<N; i++){
        double omega_i = m.solid_angle(i) * negOneOver2pi;
        collocA(i,i) = omega_i - collocA.row(i).sum(); 
    }

    /*
    std::cout<<"collocA: "<<collocA<<std::endl;
    std::cout<<"collocB_air: "<<collocB_air<<std::endl;
    std::cout<<"collocB_solid: "<<collocB_solid<<std::endl;
    */

    // re-arranging
    for (size_t i=0; i<N; i++){
        double omega_i = m.solid_angle(i) * negOneOver2pi;
        if(m.is_air[i]){
            A.col(i) = -(collocB_solid.col(i) + collocB_air.col(i));
            rhs_per_vertex.col(i) = -collocA.col(i) * BC_p[i];
            rhs_per_vertex(i,i) += omega_i*BC_p[i];
        } else if (m.is_solid[i]){
            A.col(i) = collocA.col(i);
            rhs_per_vertex.col(i) = (collocB_solid.col(i) + collocB_air.col(i)) * BC_dpdn[i];
            A(i,i) += -omega_i;
        } else if (m.is_triple[i]) {
            A.col(i) = -(collocB_air.col(i));
            rhs_per_vertex.col(i) = -collocA.col(i) * BC_p[i] + collocB_solid.col(i) * BC_dpdn[i];
            rhs_per_vertex(i,i) += omega_i*BC_p[i];
        } else {
            throw std::logic_error("This should not have been reached-- all vertices should be categorized as air, solid, or triple");
        }
    } 
    Eigen::VectorXd rhs = rhs_per_vertex.rowwise().sum();

    // linear solve
    // DIRECT - LU
    //Eigen::VectorXd soln = A.lu().solve(rhs);
    // ITERATIVE - BICGSTAB
    Eigen::BiCGSTAB<Eigen::MatrixXd> solver(A);
    Eigen::VectorXd soln = solver.solve(rhs);
    
    /*
    std::cout<<"A: "<<A<<std::endl;
    std::cout<<"rhs: "<<rhs<<std::endl;
    std::cout<<"rhs_per_vertex: "<<rhs_per_vertex<<std::endl;
    std::cout<<"x: "<<soln<<std::endl;
    */

    // assembly
    //Eigen::VectorXd p = BC_p;
    //Eigen::VectorXd dpdn = BC_dpdn;
    for (size_t i=0; i<N; i++){
        if(m.is_air[i]){
            p[i] = BC_p[i];
            dpdn[i] = soln[i];
        }
        else if(m.is_solid[i]){
            p[i] = soln[i];
            dpdn[i] = BC_dpdn[i];
        }
        else if(m.is_triple[i]){
            p[i] = BC_p[i];
            dpdn[i] = soln[i]; // technically only dpdn_a component of dpdn, but we store it all together
        }
        else
            throw std::logic_error("This should not have been reached-- all vertices should be categorized as air, solid, or triple");
    }

    assert(p == p);
    assert(dpdn == dpdn);
    
}

void Sim::step_BEM_gradP(Eigen::VectorXd& BC_p, Eigen::VectorXd& BC_dpdn, Eigen::VectorXd& p, Eigen::VectorXd& dpdn){
    // gradP and update velocities
    int N = m.verts.size();
    std::vector<bool> face_is_solid = m.get_solid_faces();

    std::vector<Eigen::Vector2d> dv = std::vector<Eigen::Vector2d>(N,Eigen::Vector2d(0.0,0.0));
    // linearly interpolate p from vertices to get dpdt per face
    Eigen::VectorXd dpdt_face(N);
    for (size_t i=0; i<N; i++){
        dpdt_face[i] = (p[m.verts_from_face(i)[1]] - p[m.verts_from_face(i)[0]])/m.face_length(i);
    }
    //std::cout<<"dpdt_face: "<<dpdt_face<<std::endl;
    
    // we average back onto the vertices to get dpdt per vertex
    Eigen::VectorXd dpdt(N);
    for (size_t i=0; i<N; i++){
        Eigen::Vector2d t_i = -m.calc_vertex_tangent(i); // TODO: figure out why this needs to be sign flipped, what am I doing here that is wrong hmmm
        //std::cout<<"t_"<<i<<": "<<t_i<<std::endl;
        Eigen::Vector2d n_i = m.calc_vertex_normal(i);

        Eigen::Vector2i neighbor_faces = m.faces_from_vert(i);

        if (!m.is_triple[i]){
            dpdt[i] = (dpdt_face[neighbor_faces[0]]*m.face_length(neighbor_faces[0]) + dpdt_face[neighbor_faces[1]]*m.face_length(neighbor_faces[1])) / (m.face_length(neighbor_faces[0]) + m.face_length(neighbor_faces[1]));
            dv[i] = -(1.0/rho) * (dpdt[i]*t_i + dpdn[i]*n_i);
        } else{
            //double dpdn_S, dpdn_A;
            double dpdt_SL, dpdt_LA;
            Eigen::Vector2d n_S, n_A;
            Eigen::Vector2d t_SL, t_LA;
            if (face_is_solid[neighbor_faces[0]]){
                t_SL = m.calc_face_tangent(neighbor_faces[0]); //want it to face outward from vertex i
                n_S = m.calc_face_normal(neighbor_faces[0]);
                dpdt_SL = dpdt_face[neighbor_faces[0]];
                //dpdn_S = BC_dpdn[i];
                
                t_LA = -m.calc_face_tangent(neighbor_faces[1]);
                n_A = m.calc_face_normal(neighbor_faces[1]);
                dpdt_LA = dpdt_face[neighbor_faces[1]];
                //dpdn_A = dpdn[i];
            } else {
                t_SL = -m.calc_face_tangent(neighbor_faces[1]);
                n_S = m.calc_face_normal(neighbor_faces[1]);
                dpdt_SL = dpdt_face[neighbor_faces[1]];
                //dpdn_S = BC_dpdn[i];

                t_LA = m.calc_face_tangent(neighbor_faces[0]); //want it to face outward from vertex i
                n_A = m.calc_face_normal(neighbor_faces[0]);
                dpdt_LA = dpdt_face[neighbor_faces[0]];
                //dpdn_A = dpdn[i];
            }
            /*
            std::cout<<"for vert "<<i<<":"<<std::endl;
            std::cout<<"t_SL: "<<t_SL<<std::endl;
            std::cout<<"n_S: "<<n_S<<std::endl;
            std::cout<<"dpdt_SL: "<<dpdt_SL<<std::endl;
            std::cout<<"t_LA: "<<t_LA<<std::endl;
            std::cout<<"n_A: "<<n_A<<std::endl;
            std::cout<<"dpdt_LA: "<<dpdt_LA<<std::endl;
            */
            Eigen::MatrixXd gradp_proj = Eigen::MatrixXd::Zero(4, 2);
            Eigen::VectorXd gradp_proj_rhs = Eigen::VectorXd::Zero(4);
            gradp_proj.row(0) = n_S.transpose();        gradp_proj_rhs[0] = BC_dpdn[i];                 // dpdn on the solid side
            gradp_proj.row(1) = n_A.transpose();        gradp_proj_rhs[1] = dpdn[i];                    // dpdn on the air side
            gradp_proj.row(2) = t_SL.transpose();        gradp_proj_rhs[2] = dpdt_SL;  // dpdt on the solid side
            gradp_proj.row(3) = t_LA.transpose();        gradp_proj_rhs[3] = dpdt_LA;    // dpdn on the air side
            
            // TODO: look over this later
            Eigen::Vector2d gradp = (gradp_proj.transpose()*gradp_proj).partialPivLu().solve(gradp_proj.transpose()*gradp_proj_rhs); // stole this from 3D SFL
            // project the pressure gradient to remove the solid normal component, because that component should be the prescribed BC_dpdn[i]  
            //gradp = (Eigen::Matrix2d::Identity() - n_S * n_S.transpose()) * gradp + n_S * BC_dpdn[i];
            /*
            std::cout<<"gradp_proj: "<<gradp_proj<<std::endl;
            std::cout<<"gradp_proj_rhs: "<<gradp_proj_rhs<<std::endl;
            std::cout<<"gradp: "<<gradp<<std::endl;
            */
            dv[i] = -(1.0/rho) * gradp;
        }
    }
    
    //std::cout<<"dpdt: "<<dpdt<<std::endl;
    
    /*for (size_t i=0; i<m.verts.size(); i++){
        std::cout<<dv[i][0]*-rho<<", "<<dv[i][1]*-rho<<std::endl;
    }*/
    // accelerate to end of step velocity
    for (size_t i=0; i<m.verts.size(); i++){
        m.vels[i] += dv[i];
    }
}

Eigen::Vector2d Sim::HHD_FD(Eigen::Vector2d x, double delta){
    // gonna do HHD here again yay to evaluate/approximate velocity at x
    // except no need to worry about singularities???
    std::vector<Eigen::Vector2d> x_deltas = {
        Eigen::Vector2d(x.x() + delta, x.y()),
        Eigen::Vector2d(x.x() - delta, x.y()),
        Eigen::Vector2d(x.x(), x.y() + delta),
        Eigen::Vector2d(x.x(), x.y() - delta)
    };
    std::vector<double> phi_gammas(x_deltas.size(),0.0);
    std::vector<double> A_gammas(x_deltas.size(),0.0);
    for(size_t i=0; i<x_deltas.size(); i++){
        phi_gammas[i] = BIE_Phi(x_deltas[i]);
        A_gammas[i] = BIE_A(x_deltas[i]);
    }
    Eigen::Vector2d gradPhi_FD(
        (phi_gammas[0]-phi_gammas[1])/(2*delta),
        (phi_gammas[2]-phi_gammas[3])/(2*delta)
    );
    Eigen::Vector2d curlA_FD(
        (A_gammas[2]-A_gammas[3])/(2*delta),
        -(A_gammas[0]-A_gammas[1])/(2*delta)
    );
    
    return -(curlA_FD - gradPhi_FD); //idk why this sign needs to be flipped here???
}

double Sim::BIE_Phi(Eigen::Vector2d x){
    // evaluating Phi_Gamma at point x
    const std::vector<Eigen::Vector2d> quadrature_GQ = BoundaryIntegral::gaussian_quadrature();
    double jacobian = 0.5;

    double phi = 0;
    for (size_t i=0; i<m.faces.size(); i++){
        double f_len = m.face_length(i);
        Eigen::Vector2d n_y = m.calc_face_normal(i);

        double phi_face = 0;
        for (size_t qi = 0; qi<quadrature_GQ.size(); qi++){
            double qik = quadrature_GQ[qi].x();
            double qiw = quadrature_GQ[qi].y();
            
            Eigen::Vector2d y = lin_interp(m.verts[m.verts_from_face(i)[0]], m.verts[m.verts_from_face(i)[1]], qik);
            Eigen::Vector2d v_y = lin_interp(m.vels[m.verts_from_face(i)[0]], m.vels[m.verts_from_face(i)[1]], qik);

            double G = BoundaryIntegral::G(x,y);

            phi_face += qiw * (n_y.dot(v_y)) * G;
        }
        phi += phi_face * (jacobian * f_len);
    }
    return phi;
}

double Sim::BIE_A(Eigen::Vector2d x){
    // evaluating A_Gamma at point x
    const std::vector<Eigen::Vector2d> quadrature_GQ = BoundaryIntegral::gaussian_quadrature();
    double jacobian = 0.5;

    double A = 0;
    for (size_t i=0; i<m.faces.size(); i++){
        double f_len = m.face_length(i);
        Eigen::Vector2d n_y = m.calc_face_normal(i);

        double A_face = 0;
        for (size_t qi = 0; qi<quadrature_GQ.size(); qi++){
            double qik = quadrature_GQ[qi].x();
            double qiw = quadrature_GQ[qi].y();
            
            Eigen::Vector2d y = lin_interp(m.verts[m.verts_from_face(i)[0]], m.verts[m.verts_from_face(i)[1]], qik);
            Eigen::Vector2d v_y = lin_interp(m.vels[m.verts_from_face(i)[0]], m.vels[m.verts_from_face(i)[1]], qik);

            double G = BoundaryIntegral::G(x,y);

            A_face += qiw * (cross2d(n_y,v_y)) * G;
        }
        A += A_face * (jacobian * f_len);
    }
    return A;
}

Eigen::Vector2d Sim::lin_interp(Eigen::Vector2d v_a, Eigen::Vector2d v_b, double q){
    // maps q \in [-1, 1] to x \in [a, b]
    return v_a*((1-q)/2) + v_b*((1+q)/2);
}

double Sim::M_1(double t){
    return (1.0-t)/2.0;
}
double Sim::M_2(double t){
    return (1.0+t)/2.0;
}

double Sim::cross2d(Eigen::Vector2d a, Eigen::Vector2d b){
   return a.x()*b.y() - a.y()*b.x();
}

void Sim::remesh(){
    for(size_t i=0; i<6; i++){
        m.remesh();
        collide();
    }
    n = m.verts.size();
}

void Sim::collide(){

    // collide liquid mesh with each solid
    // then recalibrate triple points

    m.reset_boundary_types();
    for (size_t i=0; i<solids.size(); i++){
        solids[i]->collideAndSnap(m);
    }
    m.update_triple_points();
}
void Sim::genMarkerParticles(double l, double r, double b, double t, double spacing){
    double min_dist_to_liquid_allowed = 0.05;

    std::vector<Eigen::Vector2d> tmpMarkers;
    tmpMarkers.reserve(int(((r-l)/spacing)*((t-b)/spacing)));

    for (double x=l; x<=r; x+=spacing){
        for (double y=b; y<=t; y+=spacing){
            Eigen::Vector2d marker(x,y);
            double dist_to_liquid = m.signed_min_dist(marker);
            if ( dist_to_liquid>min_dist_to_liquid_allowed){
                tmpMarkers.push_back(marker);
            }
        }
    }

    markerparticles.resize(tmpMarkers.size());
    for(size_t i = 0; i<tmpMarkers.size(); i++){
        markerparticles[i] = tmpMarkers[i];
    }
}