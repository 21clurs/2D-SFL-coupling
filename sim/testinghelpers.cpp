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

void TestingHelpers::testErrorTables(std::string testType, std::string shape, std::string fieldType){
    std::vector<size_t> N = {8, 16, 32, 64, 128, 256, 512, 1024};
    std::vector<Eigen::Vector2d> errs(N.size());
    for (size_t i=0; i<N.size(); i++){
        if (testType.compare("HHD") == 0){
            testHHD(shape, fieldType, N[i], errs[i]);
        } else if (testType.compare("BEM") == 0){
            testBEM(shape, fieldType, N[i], errs[i]);
        } else {
            std::cout<<"please enter a valid test type!"<<std::endl;
        }
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
    TestingHelpers::generateVField("harmonic_1",n,verts,vels);

    // set up expected/theoretical result
    std::vector<Eigen::Vector2d> expected_post_HHD_vels(n);
    for (size_t i=0; i<n; i++){
        // a very basic expected velocity :)
        expected_post_HHD_vels[i] = vels[i];
    }

    // set up test simulation
    Mesh mesh(verts,faces,vels);
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


void TestingHelpers::testBEM(std::string shape, std::string fieldType, int n, Eigen::Vector2d& errs){
    // set up test shape
    std::vector<Eigen::Vector2d> verts(n);
    std::vector<Eigen::Vector2i> faces(n);
    genShape(shape, n, verts, faces);

    // set BC type
    std::vector<bool> is_air(n);
    std::vector<bool> is_solid(n);
    std::vector<bool> is_triple(n);
    std::string BC_type = "Dirichlet";
    if (BC_type.compare("Dirichlet") == 0){
        std::fill(is_air.begin(), is_air.end(), true);
        std::fill(is_solid.begin(), is_solid.end(), false);
        std::fill(is_triple.begin(), is_triple.end(), false);
    }

    // set up test simulation
    Mesh mesh(verts,faces);
    mesh.set_boundaries(is_air, is_solid, is_triple);
    float dt = 1.0/60.0;
    Sim s(mesh, n, dt);

    // set input and expected BC
    Eigen::VectorXd p_input = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd dpdn_input = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd p_expected(n);
    Eigen::VectorXd dpdn_expected(n);
    for( size_t i=0; i<n; i++){
        if (fieldType.compare("1") == 0){
            // p = y
            p_expected[i] = mesh.verts[i].y();
            dpdn_expected[i] = mesh.calc_vertex_normal(i).dot(Eigen::Vector2d(0.0,1.0));
        } else if (fieldType.compare("2") == 0){
            // p = x+y
            p_expected[i] = mesh.verts[i].x() + mesh.verts[i].y();
            dpdn_expected[i] = mesh.calc_vertex_normal(i).dot(Eigen::Vector2d(1.0,1.0));
        } else if (fieldType.compare("3") == 0){
            // p = y^3-3yx^2
            p_expected[i] = pow(mesh.verts[i].y(),3) -  mesh.verts[i].y() * pow(mesh.verts[i].x(),2);
            dpdn_expected[i] = mesh.calc_vertex_normal(i).dot(Eigen::Vector2d(-6.0*mesh.verts[i].x()*mesh.verts[i].y(), 3.0*pow(mesh.verts[i].y(),2)-3.0*pow(mesh.verts[i].x(),2)));
        } else if (fieldType.compare("4") == 0){
            // p = x^3-3xy^2
            p_expected[i] = pow(mesh.verts[i].x(),3) -  mesh.verts[i].x() * pow(mesh.verts[i].y(),2);
            dpdn_expected[i] = mesh.calc_vertex_normal(i).dot(Eigen::Vector2d(3.0*pow(mesh.verts[i].x(),2)-3.0*pow(mesh.verts[i].y(),2), -6.0*mesh.verts[i].x()*mesh.verts[i].y()));
        }

        if (is_solid[i]){
            dpdn_input[i] = dpdn_expected[i];
        } else if (is_air[i]){
            p_input[i] = p_expected[i];
        }
    }
    // step BEM
    Eigen::VectorXd p_test = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd dpdn_test = Eigen::VectorXd::Zero(n);
    s.step_BEM_solve(p_input, dpdn_input, p_test, dpdn_test);

    // calculate errors
    Eigen::VectorXd err_norms(n);
    for (size_t i=0; i<n; i++){
        err_norms[i] = abs(dpdn_expected[i] - dpdn_test[i]);
    }

    // return errors...
    errs[0] = err_norms.mean();// mean error
    errs[1] = err_norms.cwiseAbs().maxCoeff();// maximum error
}


void TestingHelpers::genShape(std::string shape, int n, std::vector<Eigen::Vector2d>& verts, std::vector<Eigen::Vector2i>& faces){
    if (shape.compare("circle")==0){
        float theta = 2*M_PI/n;
        double r=1;
        Eigen::Vector2d center(0.0,0.0);

        for (uint i=0; i<n; i++){
            verts[i](0) = r*cos(theta*i) + center(0);
            verts[i](1) = r*sin(theta*i) + center(1);

            faces[i] = Eigen::Vector2i(i,(i+1)%n);
        }
    } else if (shape.compare("square")==0){
        assert(((void)"n is a multiple of 4 when generating a square", n%4==0));
        
        int nPerSide = n/4;
        float sideLength = 1;
        float delta = sideLength/nPerSide;

        for(uint i=0; i<nPerSide; i++){
            // bottom
            verts[i] = Eigen::Vector2d((-sideLength/2) + i*delta, -sideLength/2);
            // right
            verts[nPerSide+i] = Eigen::Vector2d(sideLength/2, (-sideLength/2) + i*delta);
            // top
            verts[2*nPerSide+i] = Eigen::Vector2d((sideLength/2)-i*delta, sideLength/2);
            // left
            verts[3*nPerSide+i] = Eigen::Vector2d(-sideLength/2, (sideLength/2)-i*delta);
        }
        for(uint i=0; i<n; i++)
            faces[i] = Eigen::Vector2i(i,(i+1)%n);
    } else if (shape.compare("ellipse")==0){
        float theta = 2*M_PI/n;
    
        double r=1;
        Eigen::Vector2d center(0.0,0.0);

        for (uint i=0; i<n; i++){
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
        for (uint i=0; i<n; i++){
            r = a + eps*cos(m*theta*i);

            verts[i](0) = r*cos(theta*i);
            verts[i](1) = r*sin(theta*i);

            faces[i] = Eigen::Vector2i(i,(i+1)%n);
        }
    } else if (shape.compare("semicircle_v")==0){
        float theta = 2*M_PI/n;

        double r=1;
        Eigen::Vector2d center(0.0,0.0);

        for (uint i=0; i<n; i++){
            verts[i](0) = r*cos(theta*i - M_PI/2) + center(0);
            verts[i](1) = r*sin(theta*i - M_PI/2) + center(1);

            faces[i] = Eigen::Vector2i(i,(i+1)%n);
        }

        double delta = 2*r/(n/2);
        for (uint i = 0; i<n/2; i++){
            verts[i+n/2](0) = 0;
            verts[i+n/2](1) = r - delta*i;
        }
    } else if (shape.compare("semicircle_h")==0){
        assert(((void)"n is a multiple of 2 when generating a semicircle", n%4==0));
    
        float theta = 2*M_PI/n;

        double r=1;
        Eigen::Vector2d center(0.0,0.0);

        for (uint i=0; i<n; i++){
            verts[i](0) = r*cos(theta*i) + center(0);
            verts[i](1) = r*sin(theta*i) + center(1);

            faces[i] = Eigen::Vector2i(i,(i+1)%n);
        }

        double delta = 2*r/(n/2);
        for (uint i = 0; i<n/2; i++){
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
        for (uint i=0; i<n_inner; i++){
            verts[i](0) = r_inner*cos(theta_inner*i) + center(0);
            verts[i](1) = r_inner*sin(theta_inner*i) + center(1);

            faces[i] = Eigen::Vector2i((i+1)%n_inner,i);
        }

        float theta_outer = 2*M_PI/n_outer;
        double r_outer=r_inner*ratio;
        for (uint i=0; i<n_outer; i++){
            verts[n_inner+i](0) = r_outer*cos(theta_outer*i) + center(0);
            verts[n_inner+i](1) = r_outer*sin(theta_outer*i) + center(1);

            faces[n_inner+i] = Eigen::Vector2i(n_inner+i,n_inner+(i+1)%n_outer);
        }
    } else {
        std::cout<<"please enter a valid shape type!"<<std::endl;
        return;
    }
}