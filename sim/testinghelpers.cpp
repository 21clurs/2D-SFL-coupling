#include "testinghelpers.h"

void TestingHelpers::generateVField(std::string fieldType, int n, std::vector<Eigen::Vector2d>& verts, std::vector<Eigen::Vector2d>& vels){
    for (size_t i=0; i<n; i++){
        if (fieldType.compare("zero") == 0){
            // (0,0)
            vels[i].x() = 0;
            vels[i].y() = 0;
        } else if (fieldType.compare("constant") == 0){
            // (1,0)
            vels[i].x() = 1;
            vels[i].y() = 0;
        } else if (fieldType.compare("linear_x") == 0){
            // (x,0)
            vels[i].x() = verts[i].x();
            vels[i].y() = 0;
        } else if (fieldType.compare("linear_y") == 0){
            // (0,y)
            vels[i].x() = 0;
            vels[i].y() = verts[i].y();
        } else if (fieldType.compare("harmonic_1") == 0){
            // (y,x)
            vels[i].x() = verts[i].y();
            vels[i].y() = verts[i].x();
        } else if (fieldType.compare("harmonic_2") == 0){
            // (x,-y)
            vels[i].x() = verts[i].x();
            vels[i].y() = -verts[i].y();
        } else if (fieldType.compare("harmonic_3") == 0){
            // (-x^2+y^2,2xy)
            vels[i].x() = -verts[i].x()*verts[i].x() + verts[i].y()*verts[i].y();
            vels[i].y() = 2 * verts[i].x() * verts[i].y();
        } else {
            std::cout<<"please enter a valid test velocity field!"<<std::endl;
        }
    }
}

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
    // set up test shape
    std::vector<Eigen::Vector2d> verts(n);
    std::vector<Eigen::Vector2i> faces(n);
    genShape(shape, n, verts, faces);

    // set up test velocities -- this test only makes sense when this is a harmonic velocity field
    std::vector<Eigen::Vector2d> vels(n);
    TestingHelpers::generateVField(fieldType,n,verts,vels);

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
    genShape(shape, n, verts, faces);

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


void TestingHelpers::genShape(std::string shape, int n, std::vector<Eigen::Vector2d>& verts, std::vector<Eigen::Vector2i>& faces){
    if (shape.compare("circle")==0){
        float theta = 2*M_PI/n;
        double r=1;
        Eigen::Vector2d center(0.0,0.0);

        for (size_t i=0; i<n; i++){
            verts[i](0) = r*cos(theta*i) + center(0);
            verts[i](1) = r*sin(theta*i) + center(1);

            faces[i] = Eigen::Vector2i(i,(i+1)%n);
        }
    } else if (shape.compare("square")==0){
        assert(((void)"n is a multiple of 4 when generating a square", n%4==0));
        
        int nPerSide = n/4;
        float sideLength = 1;
        float delta = sideLength/nPerSide;

        for(size_t i=0; i<nPerSide; i++){
            // bottom
            verts[i] = Eigen::Vector2d((-sideLength/2) + i*delta, -sideLength/2);
            // right
            verts[nPerSide+i] = Eigen::Vector2d(sideLength/2, (-sideLength/2) + i*delta);
            // top
            verts[2*nPerSide+i] = Eigen::Vector2d((sideLength/2)-i*delta, sideLength/2);
            // left
            verts[3*nPerSide+i] = Eigen::Vector2d(-sideLength/2, (sideLength/2)-i*delta);
        }
        for(size_t i=0; i<n; i++)
            faces[i] = Eigen::Vector2i(i,(i+1)%n);
    } else if (shape.compare("rectangle")==0){ 
        assert(((void)"n is a multiple of 2 when generating a rectangle", n%2==0));
        float w = 3;
        float h = 0.5;
        float ratio = w/h;

        float bottom_y = -0.5;
        
        int nVertical =  (n/2)/(ratio+1);
        int nHorizontal = (n/2)-nVertical;

        // bottom
        float delta = w/nHorizontal;
        for(size_t i=0; i<nHorizontal; i++){
            verts[i] = Eigen::Vector2d((-w/2) + i*delta, bottom_y);
        }
        // right
        delta = h/nVertical;
        for(size_t i=0; i<nVertical; i++){
            verts[i+nHorizontal] = Eigen::Vector2d((w/2), bottom_y + i*delta);
        }
        // top
        delta = w/nHorizontal;
        for(size_t i=0; i<nHorizontal; i++){
            verts[i+nHorizontal+nVertical] = Eigen::Vector2d((w/2)-i*delta, h+bottom_y);
        }
        // left
        delta = h/nVertical;
        for(size_t i=0; i<nVertical; i++){
            verts[i+2*nHorizontal+nVertical] = Eigen::Vector2d(-(w/2), (h+bottom_y)-i*delta);
        }
        for(size_t i=0; i<n; i++)
            faces[i] = Eigen::Vector2i(i,(i+1)%n);
    } else if (shape.compare("ellipse")==0){
        float theta = 2*M_PI/n;
    
        double r=1;
        Eigen::Vector2d center(0.0,0.0);

        for (size_t i=0; i<n; i++){
            verts[i](0) = 2*r*cos(theta*i) + center(0);
            verts[i](1) = r*sin(theta*i) + center(1);

            faces[i] = Eigen::Vector2i(i,(i+1)%n);
        }
    } else if (shape.compare("oscillation")==0){
        float theta = 2*M_PI/n;

        double a = 1.0/3.0;
        double eps = 0.05*a;
        int m = 2; // the mode of oscillation we are interested in

        double r;
        for (size_t i=0; i<n; i++){
            r = a + eps*cos(m*theta*i);

            verts[i](0) = r*cos(theta*i);
            verts[i](1) = r*sin(theta*i);

            faces[i] = Eigen::Vector2i(i,(i+1)%n);
        }
    } else if (shape.compare("semicircle_v")==0){
        float theta = 2*M_PI/n;

        double r=1;
        Eigen::Vector2d center(0.0,0.0);

        for (size_t i=0; i<n; i++){
            verts[i](0) = r*cos(theta*i - M_PI/2) + center(0);
            verts[i](1) = r*sin(theta*i - M_PI/2) + center(1);

            faces[i] = Eigen::Vector2i(i,(i+1)%n);
        }

        double delta = 2*r/(n/2);
        for (size_t i = 0; i<n/2; i++){
            verts[i+n/2](0) = 0;
            verts[i+n/2](1) = r - delta*i;
        }
    } else if (shape.compare("semicircle_h")==0){
        assert(((void)"n is a multiple of 2 when generating a semicircle", n%4==0));
    
        float theta = 2*M_PI/n;

        double r=1;
        Eigen::Vector2d center(0.0,0.0);

        for (size_t i=0; i<n; i++){
            verts[i](0) = r*cos(theta*i) + center(0);
            verts[i](1) = r*sin(theta*i) + center(1);

            faces[i] = Eigen::Vector2i(i,(i+1)%n);
        }

        double delta = 2*r/(n/2);
        for (size_t i = 0; i<n/2; i++){
            verts[i+n/2](0) = -r + delta*i;
            verts[i+n/2](1) = 0;
        }
    } else if (shape.compare("donut")==0){
        double ratio = 2.0;
        int n_inner = (n/(ratio+1))+1;
        int n_outer = n-n_inner;

        Eigen::Vector2d center(0.0,0.0);

        float theta_inner = 2*M_PI/n_inner;
        double r_inner = 0.5;
        for (size_t i=0; i<n_inner; i++){
            verts[i](0) = r_inner*cos(theta_inner*i) + center(0);
            verts[i](1) = r_inner*sin(theta_inner*i) + center(1);

            faces[i] = Eigen::Vector2i((i+1)%n_inner,i);
        }

        float theta_outer = 2*M_PI/n_outer;
        double r_outer=r_inner*ratio;
        for (size_t i=0; i<n_outer; i++){
            verts[n_inner+i](0) = r_outer*cos(theta_outer*i) + center(0);
            verts[n_inner+i](1) = r_outer*sin(theta_outer*i) + center(1);

            faces[n_inner+i] = Eigen::Vector2i(n_inner+i,n_inner+(i+1)%n_outer);
        }
    } else if (shape.compare("square_donut")==0){
        double ratio = 4.0;
        int n_inner = (n/(ratio+1));
        int n_outer = n-n_inner;

        assert(((void)"n_inner is a multiple of 4 when generating a square donut", n_inner%4==0));
        
        int nPerSide_inner = n_inner/4;
        float sideLength_inner = 0.5;
        float delta_inner = sideLength_inner/nPerSide_inner;

        int nPerSide_outer = n_outer/4;
        float sideLength_outer = sideLength_inner*ratio;
        float delta_outer = sideLength_outer/nPerSide_outer;

        // populate outer square
        for(size_t i=0; i<nPerSide_outer; i++){
            // bottom
            verts[i] = Eigen::Vector2d((-sideLength_outer/2) + i*delta_outer, -sideLength_outer/2);
            // right
            verts[nPerSide_outer+i] = Eigen::Vector2d(sideLength_outer/2, (-sideLength_outer/2) + i*delta_outer);
            // top
            verts[2*nPerSide_outer+i] = Eigen::Vector2d((sideLength_outer/2)-i*delta_outer, sideLength_outer/2);
            // left
            verts[3*nPerSide_outer+i] = Eigen::Vector2d(-sideLength_outer/2, (sideLength_outer/2)-i*delta_outer);
        }
        for(size_t i=0; i<n_outer; i++)
            faces[i] = Eigen::Vector2i(i,(i+1)%n_outer);

        // populate inner square
        for(size_t i=0; i<nPerSide_inner; i++){
            // bottom
            verts[n_outer + i] = Eigen::Vector2d((-sideLength_inner/2) + i*delta_inner, -sideLength_inner/2);
            // right
            verts[n_outer + nPerSide_inner+i] = Eigen::Vector2d(sideLength_inner/2, (-sideLength_inner/2) + i*delta_inner);
            // top
            verts[n_outer + 2*nPerSide_inner+i] = Eigen::Vector2d((sideLength_inner/2)-i*delta_inner, sideLength_inner/2);
            // left
            verts[n_outer + 3*nPerSide_inner+i] = Eigen::Vector2d(-sideLength_inner/2, (sideLength_inner/2)-i*delta_inner);
        }
        for(size_t i=0; i<n_inner; i++)
            faces[n_outer + i] = Eigen::Vector2i(n_outer + (i+1)%n_inner, n_outer + i);
    } else {
        std::cout<<"\""<< shape<<"\" is not a valid shape type. Please enter a valid shape type!"<<std::endl;
        return;
    }
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
