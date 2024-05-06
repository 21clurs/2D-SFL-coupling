#include "sim.h"

#include <iostream>
#include <fstream> 
#include "generateshapes.h"
#include "boundaryintegral.h"

using Eigen::MatrixXd;
using Eigen::Vector2d;
using Eigen::Vector2i;

Sim::Sim(Mesh& m, int n, float dt):
    n(n),
    dt(dt),
    m(m)
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

    //for(size_t i =0; i<quadrature_DE.size(); i++){
    //    std::cout<<quadrature_DE[i].x()<<","<<quadrature_DE[i].y()<<std::endl;
    //}

    // i for outer -- when looping through verts, and then the faces adjacent to those verts
    // j for inner

    // iterate over vertices

    const double negOneOver2pi = -1.0/(2*M_PI);

    std::vector<double> dPhidn(n, 0.0);
    std::vector<double> curlA(n, 0.0);

    for (size_t i=0; i<m.verts.size(); i++){
        //std::cout<<"i"<<": "<<i<<std::endl;
        int xi_a, xi_b;
        // iterate over each adjacent edge (only 2 in this case)
        for (size_t p=0; p<2; p++){
            int f_i = m.faces_from_vert(i)[p];

            Eigen::Vector2d n_x = m.calc_face_normal(f_i);
            Eigen::Vector2d t_x = m.calc_face_tangent(f_i);
            /*if (i==0 ){
                std::cout<<"f_i"<<": "<<f_i<<std::endl;
                //std::cout<<"t_x"<<": "<<t_x<<std::endl;
            }*/
            
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

            double jacobian = 1.0/2.0; //super sketch but ok lets see...

            // iterate over quadrature points
            for (size_t qi=0; qi<quadrature_DE.size(); qi++){
                double qik = quadrature_DE[qi].x();
                double qiw = quadrature_DE[qi].y();

                double theta_i = (1-qik)/2;

                Eigen::Vector2d x = lin_interp(m.verts[xi_a], m.verts[xi_b], qik);
                Eigen::Vector2d v_x = lin_interp(m.vels[xi_a], m.vels[xi_b], qik);
                /*if (i==0 ){
                    std::cout<<"quad_x"<<": "<<x<<std::endl;
                }*/
                //std::cout<<"quad_pt_v_"<<i+1<<": "<<v_x<<std::endl;

                double I_gradPhi = 0;
                double I_curlA = 0;

                // iterating over all faces (surface integral over the entire boundary)
                for (size_t f_j=0; f_j<m.faces.size(); f_j++){
                    Eigen::Vector2d n_y = m.calc_face_normal(f_j);
                    //Eigen::Vector2d t_y = m.calc_face_tangent(f_j);
                    /*if (i==0 ){
                        std::cout<<"f_j:"<<f_j<<std::endl;
                        std::cout<<"n_y"<<": "<<n_y<<std::endl;
                    }*/

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
                                /*if (i==0 &&qi==0){
                                    std::cout<<"y"<<": "<<y<<std::endl;
                                    std::cout<<"val"<<": "<<cross2d(n_y,v_y) * BoundaryIntegral::G(x,y)<<std::endl;
                                    std::cout<<"val_curlA * len"<<": "<<val_curlA * len<<std::endl;
                                }*/
                            }
                            I_curlA += val_curlA * len; // or f_len, like it is in my matlab?
                            /*if (i==0 ){
                                std::cout<<"SING_val_curlA"<<"_k_"<<k<<": "<<val_curlA<<std::endl;
                                std::cout<<"SING_val_curlA * len"<<"_k_"<<k<<": "<< val_curlA * len<<std::endl;
                                std::cout<<"SING_I_curlA"<<"_k_"<<k<<": "<< I_curlA <<std::endl;
                            }*/
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
                        /*if (i==0 ){
                            std::cout<<"val_gradPhi"<<": "<<val_gradPhi<<std::endl;
                            std::cout<<"val_gradPhi * len"<<": "<< val_gradPhi * f_len<<std::endl;
                            std::cout<<"I_gradPhi"<<": "<< I_gradPhi <<std::endl;
                        }*/

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
                        //if (i==0 && qi==0)
                            //std::cout<<"val_curlA_f_j"<<f_j<<": "<<val_curlA<<std::endl;
                    }
                }
                II_gradPhi += theta_i * jacobian * qiw * I_gradPhi;
                II_curlA += jacobian * qiw*I_curlA;
                /*if (i==0 ){
                    std::cout<<"II_gradPhi"<<": "<< II_gradPhi <<std::endl;
                }*/
                charge += theta_i * jacobian * qiw * m.solid_angle(i)*v_x.dot(n_x);
            }
            dPhidn[i] += (II_gradPhi + charge)*m.face_length(f_i);
            curlA[i] += II_curlA*m.face_length(f_i);
            /*if (i==0 ){
                std::cout<<"dPhidn_i=0"<<": "<< dPhidn[i] <<std::endl;
                std::cout<<"--------------------------" <<std::endl;
            }*/
        }
    }

    for(size_t i =0; i<n; i++){
       //std::cout<<"dPhidn_"<<i+1<<": "<<dPhidn[i]<<std::endl;
       //std::cout<<"curlA_"<<i+1<<": "<<curlA[i]<<std::endl;
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
/*
    for(size_t i =0; i<n; i++){
       std::cout<<"dPhidn_"<<i+1<<": "<<dPhidn[i]<<std::endl;
    }
    for(size_t i =0; i<n; i++){
       std::cout<<"curlA_"<<i+1<<": "<<curlA[i]<<std::endl;
    }*/
}

void Sim::step_BEM(){
    // BC
}

Eigen::Vector2d Sim::lin_interp(Eigen::Vector2d v_a, Eigen::Vector2d v_b, double q){
    // maps q \in [-1, 1] to x \in [a, b]
    return v_a*((1-q)/2) + v_b*((1+q)/2);
}

double Sim::cross2d(Eigen::Vector2d a, Eigen::Vector2d b){
   return a.x()*b.y() - a.y()*b.x();
}