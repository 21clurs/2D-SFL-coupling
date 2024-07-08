#include "testinghelpers.h"
#include "scenes.h"

void TestingHelpers::testHHDErrorTables(std::string shape, std::string fieldType){
    std::vector<size_t> N = {4, 8, 16, 32, 64, 128, 256, 512};
    std::vector<Eigen::Vector2d> errs(N.size());
    for (size_t i=0; i<N.size(); i++){
        testHHD(shape, fieldType, N[i], errs[i]);
        // progress messages
        std::cout<<"Test "<<i+1<<"/"<<N.size()<<" complete."<<"\r";
        std::cout.flush();
    }
    std::cout<<std::endl;

    std::cout<<"N"<<"\t"<<"err_mean"<<"\t"<<"err_max"<<"\t"<<"err_mean_ratios"<<"\t"<<"err_max_ratios"<<std::endl;
    std::cout<<"----------------------------------------------------------------"<<std::endl;
    for (size_t i=0; i<N.size(); i++){
        if (i>0){
            std::cout<<N[i]<<"\t"<<errs[i].x()<<"\t"<<errs[i].y()<<"\t"<<errs[i].x()/errs[i-1].x()<<"\t"<<errs[i].y()/errs[i-1].y()<<std::endl;
        } else{
            std::cout<<N[i]<<"\t"<<errs[i].x()<<"\t"<<errs[i].y()<<"\t"<<"n/a"<<"\t"<<"n/a"<<std::endl;
        }
    }
}

void TestingHelpers::testHHD(std::string shape, std::string fieldType, int n, Eigen::Vector2d& errs){
    // set up test mesh
    std::vector<Eigen::Vector2d> verts(n);
    std::vector<Eigen::Vector2i> faces(n);
    LiquidMesh m(verts,faces);
    Scenes::setupSceneShape(m, shape);

    // set up test velocities -- this test only makes sense when this is a harmonic velocity field
    std::vector<Eigen::Vector2d> vels(n);
    Scenes::setupSceneVelocities(verts, vels, fieldType);

    // set up expected/theoretical result
    std::vector<Eigen::Vector2d> expected_post_HHD_vels(n);
    for (size_t i=0; i<n; i++){
        // a very basic expected velocity :)
        expected_post_HHD_vels[i] = vels[i];
    }

    // set up test simulation
    LiquidMesh mesh(verts,faces,vels);
    float dt = 1.0/60.0;
    Sim s(mesh, n, dt);
    
    // step HHD
    s.step_HHD();

    // calculate errors
    Eigen::VectorXd err_norms(n);
    for (size_t i=0; i<n; i++){
        err_norms[i] = (s.get_vels()[i] - expected_post_HHD_vels[i]).norm();
    }

    // return errors...
    errs[0] = err_norms.mean();// mean error
    errs[1] = err_norms.cwiseAbs().maxCoeff();// maximum error
}

void TestingHelpers::testBEMErrorTables(std::string shape, std::string fieldType){
    std::vector<size_t> N = {8, 16, 32, 64, 128, 256, 512, 1024};
    std::vector<std::vector<Eigen::Vector2d>> errs(N.size(), std::vector<Eigen::Vector2d>(3,Eigen::Vector2d::Zero()));
    for (size_t i=0; i<N.size(); i++){
        testBEM(shape, fieldType, N[i], errs[i]);
        // progress messages
        std::cout<<"Test "<<i+1<<"/"<<N.size()<<" complete."<<"\r";
        std::cout.flush();
    }
    std::cout<<std::endl;

    for (size_t i=0; i<3; i++){
        switch(i) {
            case 0:
                std::cout<<"Error table for p:"<<std::endl;
                break;
            case 1:
                std::cout<<"Error table for dpdn:"<<std::endl;
                break;
            case 2:
                std::cout<<"Error table for dpdn_a:"<<std::endl;
                break;
        }
        std::cout<<"N"<<"\t"<<"err_mean"<<"\t"<<"err_max"<<"\t"<<"err_mean_ratios"<<"\t"<<"err_max_ratios"<<std::endl;
        std::cout<<"----------------------------------------------------------------"<<std::endl;
        for (size_t j=0; j<N.size(); j++){
            if (j>0){
                std::cout<<N[j]<<"\t"<<errs[j][i].x()<<"\t"<<errs[j][i].y()<<"\t"<<errs[j][i].x()/errs[j-1][i].x()<<"\t"<<errs[j][i].y()/errs[j-1][i].y()<<std::endl;
            } else{
                std::cout<<N[j]<<"\t"<<errs[j][i].x()<<"\t"<<errs[j][i].y()<<"\t"<<"n/a"<<"\t"<<"n/a"<<std::endl;
            }
        }
    }
}

void TestingHelpers::testBEM(std::string shape, std::string fieldType, int n, std::vector<Eigen::Vector2d>& errs){    
    // set up test shape
    std::vector<Eigen::Vector2d> verts(n);
    std::vector<Eigen::Vector2i> faces(n);
    LiquidMesh m(verts,faces);
    Scenes::setupSceneShape(m, shape);

    // set BC type
    std::vector<bool> is_air(n);
    std::vector<bool> is_solid(n);
    std::vector<bool> is_triple(n);
    std::string BC_type = "Fun!";
    if (BC_type.compare("Dirichlet") == 0){
        std::fill(is_air.begin(), is_air.end(), true);
        std::fill(is_solid.begin(), is_solid.end(), false);
        std::fill(is_triple.begin(), is_triple.end(), false);
    } else if (BC_type.compare("Mixed") == 0){
        for (size_t i=0; i<n; i++){
            if (i<n/2){
                is_air[i] = false;
                is_solid[i] = true;
                is_triple[i] = false;
            } else{
                is_air[i] = true;
                is_solid[i] = false;
                is_triple[i] = false;
            }
        }
    } else if (BC_type.compare("Fun!") == 0){
        setFunBC(shape, n, is_air, is_solid, is_triple);
    }

    // set up test simulation
    LiquidMesh mesh(verts,faces);
    mesh.set_boundaries(is_air, is_solid, is_triple);
    float dt = 1.0/60.0;
    Sim s(mesh, n, dt);

    // set input and expected BC
    Eigen::VectorXd p_input = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd dpdn_input = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd p_expected(n);
    Eigen::VectorXd dpdn_expected(n);
    std::vector<bool> face_is_solid = mesh.get_solid_faces();
    for( size_t i=0; i<n; i++){
        if (fieldType.compare("1") == 0){
            // p = y
            p_expected[i] = mesh.verts[i].y();
            if (is_triple[i]){
                for (int face: mesh.faces_from_vert(i)){
                    if (!face_is_solid[face]){
                        dpdn_expected[i] = mesh.calc_face_normal(face).dot(Eigen::Vector2d(0.0,1.0));
                    } else{
                        dpdn_input[i] = mesh.calc_face_normal(face).dot(Eigen::Vector2d(0.0,1.0));
                    }
                }
            } else {
                dpdn_expected[i] = mesh.calc_vertex_normal(i).dot(Eigen::Vector2d(0.0,1.0));
            }
        } else if (fieldType.compare("2") == 0){
            // p = x+y
            p_expected[i] = mesh.verts[i].x() + mesh.verts[i].y();
            if (is_triple[i]){
                for (int face: mesh.faces_from_vert(i)){
                    if (!face_is_solid[face]){
                        dpdn_expected[i] = mesh.calc_face_normal(face).dot(Eigen::Vector2d(1.0,1.0));
                    } else{
                        dpdn_input[i] = mesh.calc_face_normal(face).dot(Eigen::Vector2d(1.0,1.0));
                    }
                }
            } else {
                dpdn_expected[i] = mesh.calc_vertex_normal(i).dot(Eigen::Vector2d(1.0,1.0));
            }
        } else if (fieldType.compare("3") == 0){
            // p = y^3-3yx^2
            p_expected[i] = pow(mesh.verts[i].y(),3) -  3 * mesh.verts[i].y() * pow(mesh.verts[i].x(),2);
            if (is_triple[i]){
                for (int face: mesh.faces_from_vert(i)){
                    if (!face_is_solid[face]){
                        dpdn_expected[i] = mesh.calc_face_normal(face).dot(Eigen::Vector2d(-6.0*mesh.verts[i].x()*mesh.verts[i].y(), 3.0*pow(mesh.verts[i].y(),2)-3.0*pow(mesh.verts[i].x(),2)));
                    } else{
                        dpdn_input[i] = mesh.calc_face_normal(face).dot(Eigen::Vector2d(-6.0*mesh.verts[i].x()*mesh.verts[i].y(), 3.0*pow(mesh.verts[i].y(),2)-3.0*pow(mesh.verts[i].x(),2)));
                    }
                }
            } else {
                dpdn_expected[i] = mesh.calc_vertex_normal(i).dot(Eigen::Vector2d(-6.0*mesh.verts[i].x()*mesh.verts[i].y(), 3.0*pow(mesh.verts[i].y(),2)-3.0*pow(mesh.verts[i].x(),2)));
            }
        } else if (fieldType.compare("4") == 0){
            // p = x^3-3xy^2
            p_expected[i] = pow(mesh.verts[i].x(),3) -  3 * mesh.verts[i].x() * pow(mesh.verts[i].y(),2);
            if (is_triple[i]){
                for (int face: mesh.faces_from_vert(i)){
                    if (!face_is_solid[face]){
                        dpdn_expected[i] = mesh.calc_face_normal(face).dot(Eigen::Vector2d(3.0*pow(mesh.verts[i].x(),2)-3.0*pow(mesh.verts[i].y(),2), -6.0*mesh.verts[i].x()*mesh.verts[i].y()));
                    } else{
                        dpdn_input[i] = mesh.calc_face_normal(face).dot(Eigen::Vector2d(3.0*pow(mesh.verts[i].x(),2)-3.0*pow(mesh.verts[i].y(),2), -6.0*mesh.verts[i].x()*mesh.verts[i].y()));
                    }
                }
            } else {
                dpdn_expected[i] = mesh.calc_vertex_normal(i).dot(Eigen::Vector2d(3.0*pow(mesh.verts[i].x(),2)-3.0*pow(mesh.verts[i].y(),2), -6.0*mesh.verts[i].x()*mesh.verts[i].y()));
            }
        }

        if (is_solid[i]){
            dpdn_input[i] = dpdn_expected[i];
        } else if (is_air[i]){
            p_input[i] = p_expected[i];
        } else if (is_triple[i]){
            p_input[i] = p_expected[i];
        }
    }
    
    // step BEM
    Eigen::VectorXd p_test = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd dpdn_test = Eigen::VectorXd::Zero(n);
    s.step_BEM_solve(p_input, dpdn_input, p_test, dpdn_test);

    // calculate errors
    Eigen::VectorXd p_err_norms = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd dpdn_err_norms = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd dpdna_err_norms = Eigen::VectorXd::Zero(n);

    for (size_t i=0; i<n; i++){
        if(is_solid[i])
            p_err_norms[i] = abs(p_expected[i] - p_test[i]);
        else if(is_air[i])
            dpdn_err_norms[i] = abs(dpdn_expected[i] - dpdn_test[i]);
        else if(is_triple[i])
            dpdna_err_norms[i] = abs(dpdn_expected[i] - dpdn_test[i]);
    }

    errs[0][0] = p_err_norms.cwiseAbs().sum() / std::count(is_solid.begin(), is_solid.end(), true);
    errs[0][1] = p_err_norms.cwiseAbs().maxCoeff();

    errs[1][0] = dpdn_err_norms.cwiseAbs().sum() / std::count(is_air.begin(), is_air.end(), true);
    errs[1][1] = dpdn_err_norms.cwiseAbs().maxCoeff();

    errs[2][0] = dpdna_err_norms.cwiseAbs().sum() / std::count(is_triple.begin(), is_triple.end(), true);
    errs[2][1] = dpdna_err_norms.cwiseAbs().maxCoeff();
}

void TestingHelpers::setFunBC(std::string shape, int n, std::vector<bool>& is_air, std::vector<bool>& is_solid, std::vector<bool>& is_triple){
    for (size_t i=0; i<n; i++){
        if (shape.compare("semicircle_v")==0){
            if (i==0 || i==n/2){
                is_air[i] = false;
                is_solid[i] = false;
                is_triple[i] = true;
            }
            else if (i<n/2){
                is_air[i] = true;
                is_solid[i] = false;
                is_triple[i] = false;
            }
            else if(i>n/2){
                is_air[i] = false;
                is_solid[i] = true;
                is_triple[i] = false;
            }
        } else if (shape.compare("semicircle_h")==0){
            if (i==0 || i==n/2){
                is_air[i] = false;
                is_solid[i] = false;
                is_triple[i] = true;
            }
            else if (i<n/2){
                is_air[i] = true;
                is_solid[i] = false;
                is_triple[i] = false;
            }
            else if(i>n/2){
                is_air[i] = false;
                is_solid[i] = true;
                is_triple[i] = false;
            }
        } else if (shape.compare("donut")==0){
            double ratio = 2.0;
            int n_inner = (n/(ratio+1))+1;

            if (i<n_inner){
                is_air[i] = false;
                is_solid[i] = true;
                is_triple[i] = false;
            }
            else {
                is_air[i] = true;
                is_solid[i] = false;
                is_triple[i] = false;
            }
        }
    }
}
