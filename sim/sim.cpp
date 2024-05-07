#include "sim.h"

#include <iostream>
#include <fstream> 
#include <stdexcept>
#include "generateshapes.h"
#include "boundaryintegral.h"

using Eigen::MatrixXd;
using Eigen::Vector2d;
using Eigen::Vector2i;

Sim::Sim(Mesh& m, int n, float dt):
    n(n),
    dt(dt),
    m(m),
    sigma(1.0),
    sigma_SL(1.0),
    sigma_SA(1.0),
    rho(1.0),
    gravity(Eigen::Vector2d({0.0, -1.0}))
{
}


bool Sim::outputFrame(std::string filename, std::string filelocation){
    std::ofstream file(filelocation+filename);
    for (uint i=0; i<m.verts.size(); i++)
        file<<"v "<<m.verts[i][0]<<" "<<m.verts[i][1]<<std::endl;
    file<<std::endl;

    // outward vertex norms
    for (uint i=0; i<m.verts.size(); i++)
        file<<"vn "<<(m.calc_vertex_normal(i))[0]<<" "<<(m.calc_vertex_normal(i))[1]<<std::endl;
    file<<std::endl;
    
    // faces
    for (uint i=0; i<m.faces.size(); i++)
        file<<"f "<<m.faces[i][0]<<" "<<m.faces[i][1]<<std::endl;
    file.close();

    return true;
}

void Sim::step_sim(int frame){
    step_advect();
    remesh();
    step_HHD();
    step_BEM();
}

void Sim::step_advect(){
    for (size_t i=0; i<n; i++){
        m.verts[i] = m.verts[i] + m.vels[i]*dt;
    }
}

void Sim::step_HHD(){
    const std::vector<Eigen::Vector2d> quadrature_GQ = BoundaryIntegral::gaussian_quadrature();
    const std::vector<Eigen::Vector2d> quadrature_DE = BoundaryIntegral::tanh_sinh_quadrature();

    const double negOneOver2pi = -1.0/(2*M_PI);

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

            // iterate over quadrature points
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
                                    val_curlA += jacobian * qjw * cross2d(v_y,n_y) * BoundaryIntegral::G(x,y) ;
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
                charge += theta_i * jacobian * qiw * m.solid_angle(i)*v_x.dot(n_x);
                II_curlA += jacobian * qiw*I_curlA;
            }
            dPhidn[i] += (II_gradPhi + charge)*m.face_length(f_i);
            curlA[i] += II_curlA*m.face_length(f_i);
        }
    }
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

void Sim::step_BEM(){
    
    // Setting BCs
    int N = m.verts.size();
    double triple_junction_virtual_width = 0.25 * m.calc_avg_face_length();

    Eigen::VectorXd BC_p = Eigen::VectorXd::Zero(N);
    Eigen::VectorXd BC_dpdn = Eigen::VectorXd::Zero(N);

    std::vector<bool> face_is_solid = m.get_solid_faces();

    Eigen::Vector2d v_solid =  Eigen::Vector2d({0.0, 0.0}); //placeholder

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
            BC_dpdn[i] = n_i.dot(m.vels[i] + gravity*dt - v_solid) * rho;
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
            }
            else if (face_is_solid[incident_faces[1]]){
                n_solid_outward = m.calc_face_normal(incident_faces[1]);

                // TODO: make this nicer somehow? messy...
                Eigen::Vector2i solid_face_verts = m.verts_from_face(incident_faces[1]);
                solid_vert = (i == solid_face_verts[0])? solid_face_verts[1]: solid_face_verts[0];
                Eigen::Vector2i air_face_verts = m.verts_from_face(incident_faces[0]);
                air_vert = (i == air_face_verts[0])? air_face_verts[1]: air_face_verts[0];
            }
            else
                throw std::logic_error("This should not have been reached-- triple junction should fall into one of the above statements");

            Eigen::Vector2d t_solid_outward = m.verts[i] - m.verts[solid_vert];
            t_solid_outward.normalize();
            Eigen::Vector2d t_air_outward = m.verts[i] - m.verts[air_vert];
            t_air_outward.normalize();

            Eigen::Vector2d st_force_SL = -t_solid_outward * sigma * sigma_SL;
            Eigen::Vector2d st_force_SA = t_solid_outward * sigma * sigma_SA;
            Eigen::Vector2d st_force_LA = -t_air_outward * sigma;

            Eigen::Vector2d st_force_combined = st_force_SL + st_force_SA + st_force_LA;
            double pressure_jump_normal_to_triple_junction = st_force_combined.dot(-t_solid_outward) / triple_junction_virtual_width;

            BC_p[i] = pressure_jump_normal_to_triple_junction * dt;
            BC_dpdn[i] = n_solid_outward.dot(m.vels[i] + gravity*dt - v_solid) * rho;
        }
    }

    // BEM

}

Eigen::Vector2d Sim::lin_interp(Eigen::Vector2d v_a, Eigen::Vector2d v_b, double q){
    // maps q \in [-1, 1] to x \in [a, b]
    return v_a*((1-q)/2) + v_b*((1+q)/2);
}

double Sim::cross2d(Eigen::Vector2d a, Eigen::Vector2d b){
   return a.x()*b.y() - a.y()*b.x();
}

void Sim::remesh(){
    m.laplacian_smoothing();
}
