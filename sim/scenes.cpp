#include "scenes.h"
#include "simoptions.h"
#include "rigidbody.h"

namespace {
    void gen_ellipse(int n, const Eigen::Vector2d& center, double a, double b, std::vector<Eigen::Vector2d> &v, std::vector<Eigen::Vector2i> &f ){
        float theta = 2*M_PI/n;

        v.resize(n);
        f.resize(n);

        for (size_t i=0; i<n; i++){
            v[i].x() = a*cos(theta*i) + center.x();
            v[i].y() = b*sin(theta*i) + center.y();

            f[i][0] = i;
            f[i][1] = (i+1)%n;
        }
    }
    void gen_rectangle(int n, const Eigen::Vector2d& center, double w, double h, std::vector<Eigen::Vector2d> &v, std::vector<Eigen::Vector2i> &f ){
        assert(((void)"n is a multiple of 2 when generating a rectangle", n%2==0));
        
        v.resize(n);
        f.resize(n);

        float ratio = w/h;

        float bottom_y = center.y() - h/2;
        float left_x = center.x() - w/2;

        int nVertical = (n/2)/(ratio+1);
        int nHorizontal = (n/2)-nVertical;

        // bottom
        float delta = w/nHorizontal;
        for(size_t i=0; i<nHorizontal; i++){
            v[i] = Eigen::Vector2d(left_x + i*delta, bottom_y);
        }
        // right
        delta = h/nVertical;
        for(size_t i=0; i<nVertical; i++){
            v[i+nHorizontal] = Eigen::Vector2d(left_x + w, bottom_y + i*delta);
        }
        // top
        delta = w/nHorizontal;
        for(size_t i=0; i<nHorizontal; i++){
            v[i+nHorizontal+nVertical] = Eigen::Vector2d((left_x+w)-i*delta, h+bottom_y);
        }
        // left
        delta = h/nVertical;
        for(size_t i=0; i<nVertical; i++){
            v[i+2*nHorizontal+nVertical] = Eigen::Vector2d(left_x, (h+bottom_y)-i*delta);
        }

        // face indices
        for(size_t i=0; i<n; i++){
            f[i] = Eigen::Vector2i(i,(i+1)%n);
        }
    }
    void gen_semicircle_v(int n, const Eigen::Vector2d& center, double a, double b, std::vector<Eigen::Vector2d> &v, std::vector<Eigen::Vector2i> &f ){
        float theta = 2*M_PI/n;

        v.resize(n);
        f.resize(n);

        for (size_t i=0; i<n; i++){
            v[i].x() = a*cos(theta*i - M_PI/2) + center.x();
            v[i].y() = b*sin(theta*i - M_PI/2) + center.y();

            f[i][0] = i;
            f[i][1] = (i+1)%n;
        }

        double delta = 2*b/(n/2);
        for (size_t i=0; i<n/2; i++){
            v[i+n/2].x() = 0;
            v[i+n/2].y() = b - delta*i;
        }
    }
    void gen_semicircle_h(int n, const Eigen::Vector2d& center, double a, double b, std::vector<Eigen::Vector2d> &v, std::vector<Eigen::Vector2i> &f ){
        float theta = 2*M_PI/n;

        v.resize(n);
        f.resize(n);

        for (size_t i=0; i<n; i++){
            v[i].x() = a*cos(theta*i) + center.x();
            v[i].y() = b*sin(theta*i) + center.y();

            f[i][0] = i;
            f[i][1] = (i+1)%n;
        }

        double delta = 2*a/(n/2);
        for (size_t i=0; i<n/2; i++){
            v[i+n/2].x() = -a + delta*i;
            v[i+n/2].y() = 0;
        }
    }
    void gen_donut(int n, const Eigen::Vector2d& center, double r_inner, double r_outer, std::vector<Eigen::Vector2d> &v, std::vector<Eigen::Vector2i> &f ){
        v.resize(n);
        f.resize(n);
        
        double ratio = r_outer/r_inner;

        int n_inner = (n/(ratio+1))+1;
        int n_outer = n-n_inner;

        float theta_inner = 2*M_PI/n_inner;
        for (size_t i=0; i<n_inner; i++){
            v[i].x() = r_inner*cos(theta_inner*i) + center.x();
            v[i].y() = r_inner*sin(theta_inner*i) + center.y();

            f[i].x() = (i+1)%n_inner;
            f[i].y() = i;
        }

        float theta_outer = 2*M_PI/n_outer;
        for (size_t i=0; i<n_outer; i++){
            v[n_inner+i].x() = r_outer*cos(theta_outer*i) + center.x();
            v[n_inner+i].y() = r_outer*sin(theta_outer*i) + center.y();

            f[n_inner+i].x() = n_inner + i;
            f[n_inner+i].y() = n_inner + (i+1)%n_outer;
        }
    }
    void gen_square_donut(int n, const Eigen::Vector2d& center, double sideLength_inner, double sideLength_outer, std::vector<Eigen::Vector2d> &v, std::vector<Eigen::Vector2i> &f ){
        
        double ratio = sideLength_outer/sideLength_inner;

        int n_inner = n/(ratio+1);
        int n_outer = n-n_inner;

        n_inner = n_inner + (n_inner%4);
        n_outer = n_outer + (n_outer%4);
        v.resize(n_inner + n_outer);
        f.resize(n_inner + n_outer);

        assert(((void)"n_inner is a multiple of 4 when generating a square donut", n_inner%4==0));
        assert(((void)"n_outer is a multiple of 4 when generating a square donut", n_outer%4==0));
        
        int nPerSide_inner = n_inner/4;
        float delta_inner = sideLength_inner/nPerSide_inner;

        int nPerSide_outer = n_outer/4;
        float delta_outer = sideLength_outer/nPerSide_outer;

        // populate outer square
        for(size_t i=0; i<nPerSide_outer; i++){
            // bottom
            v[i] = Eigen::Vector2d((-sideLength_outer/2) + i*delta_outer, -sideLength_outer/2);
            // right
            v[nPerSide_outer+i] = Eigen::Vector2d(sideLength_outer/2, (-sideLength_outer/2) + i*delta_outer);
            // top
            v[2*nPerSide_outer+i] = Eigen::Vector2d((sideLength_outer/2)-i*delta_outer, sideLength_outer/2);
            // left
            v[3*nPerSide_outer+i] = Eigen::Vector2d(-sideLength_outer/2, (sideLength_outer/2)-i*delta_outer);
        }
        for(size_t i=0; i<n_outer; i++)
            f[i] = Eigen::Vector2i(i,(i+1)%n_outer);

        // populate inner square
        for(size_t i=0; i<nPerSide_inner; i++){
            // bottom
            v[n_outer + i] = Eigen::Vector2d((-sideLength_inner/2) + i*delta_inner, -sideLength_inner/2);
            // right
            v[n_outer + nPerSide_inner+i] = Eigen::Vector2d(sideLength_inner/2, (-sideLength_inner/2) + i*delta_inner);
            // top
            v[n_outer + 2*nPerSide_inner+i] = Eigen::Vector2d((sideLength_inner/2)-i*delta_inner, sideLength_inner/2);
            // left
            v[n_outer + 3*nPerSide_inner+i] = Eigen::Vector2d(-sideLength_inner/2, (sideLength_inner/2)-i*delta_inner);
        }
        for(size_t i=0; i<n_inner; i++)
            f[n_outer + i] = Eigen::Vector2i(n_outer + (i+1)%n_inner, n_outer + i);
    }
    void gen_generic_donut(int n, const Eigen::Vector2d& center, double r_outer, std::vector<Eigen::Vector2d> &v, std::vector<Eigen::Vector2i> &f){
        // assumes at least one rigid body that this donut spawns around
        assert(SimOptions::intValue("num-rb") >= 1);
        RigidBody *m = new RigidBody();
        RigidBody::loadMeshFromFile(*m, SimOptions::strValue("rigid-body-file-1"));

        int n_inner = m->verts.size();
        int n_outer = n;
        n_outer = n_outer + n_outer%4;

        n = n_inner + n_outer;

        // creates a rough square donut around an arbitrary rigid body shape in the middle
        v.resize(n);
        f.resize(n);
        
        // creates a rough circular donut around an arbitrary rigid body shape in the middle
        v.resize(n);
        f.resize(n);

        for (size_t i=0; i<n_inner; i++){
            v[i] = m->verts[i];
            f[i] = Eigen::Vector2i(m->faces[i].y(),m->faces[i].x());
        }

        // then the circle
        float theta = 2*M_PI/n_outer;
        for (size_t i=0; i<n_outer; i++){
            v[i+n_inner].x() = r_outer*cos(theta*i) + center.x();
            v[i+n_inner].y() = r_outer*sin(theta*i) + center.y();

            f[i+n_inner][0] = n_inner + i;
            f[i+n_inner][1] = n_inner + (i+1)%n_outer;
        }
    }
    void gen_generic_square_donut(int n, const Eigen::Vector2d& center, double sideLength_outer, std::vector<Eigen::Vector2d> &v, std::vector<Eigen::Vector2i> &f){
        // assumes at least one rigid body that this donut spawns around
        assert(SimOptions::intValue("num-rb") >= 1);
        RigidBody *m = new RigidBody();
        RigidBody::loadMeshFromFile(*m, SimOptions::strValue("rigid-body-file-1"));

        int n_inner = m->verts.size();
        int n_outer = n;
        n_outer = n_outer + n_outer%4;

        n = n_inner + n_outer;

        // creates a rough square donut around an arbitrary rigid body shape in the middle
        v.resize(n);
        f.resize(n);

        for (size_t i=0; i<n_inner; i++){
            v[i] = m->verts[i];
            f[i] = Eigen::Vector2i(m->faces[i].y(),m->faces[i].x());
        }

        // then the square
        int nPerSide = n_outer/4;
        double delta = sideLength_outer/nPerSide;
        // populate outer square
        for(size_t i=0; i<nPerSide; i++){
            // bottom
            v[n_inner + i] = Eigen::Vector2d((-sideLength_outer/2) + i*delta, -sideLength_outer/2);
            // right
            v[n_inner + nPerSide+i] = Eigen::Vector2d(sideLength_outer/2, (-sideLength_outer/2) + i*delta);
            // top
            v[n_inner + 2*nPerSide+i] = Eigen::Vector2d((sideLength_outer/2)-i*delta, sideLength_outer/2);
            // left
            v[n_inner + 3*nPerSide+i] = Eigen::Vector2d(-sideLength_outer/2, (sideLength_outer/2)-i*delta);
        }
        for(size_t i=0; i<n_outer; i++)
            f[n_inner + i] = Eigen::Vector2i(n_inner + i, n_inner + (i+1)%n_outer);
    }
}

void Scenes::scene(Sim * const &sim, const std::string & scenename, const std::string & initialvelocity){
    int N = sim->n;
    sim->m.resize_mesh(N);

    setupSceneShape(sim->m, scenename);
    setupSceneVelocities(sim->m.verts, sim->m.vels, initialvelocity);

    setupSceneSolids(sim);
}

void Scenes::sceneFromFile(Sim * const &sim, const std::string & filename, const std::string & initialvelocity){
    int N = sim->n;
    sim->m.resize_mesh(N);

    setupSceneShapeFromFile(sim->m, filename);
    setupSceneVelocities(sim->m.verts, sim->m.vels, initialvelocity);
    setupSceneSolids(sim);
}

void Scenes::setupSceneSolids(Sim * const &sim){
    for (int i=1; i<=SimOptions::intValue("num-rb"); i++){
        RigidBody *r = new RigidBody();
        RigidBody::loadMeshFromFile(*r, SimOptions::strValue("rigid-body-file-"+std::to_string(i)));
        sim->addRigidBody(r);
    }
}

void Scenes::setupSceneShape(LiquidMesh& m, const std::string & scenename){
    int N = m.verts.size();

    std::vector<Eigen::Vector2d> v(N, Eigen::Vector2d(0,0));
    std::vector<Eigen::Vector2i> f(N, Eigen::Vector2i(0,0));

    if (scenename == "circle"){

        double r = SimOptions::doubleValue("radius");

        gen_ellipse(N, Eigen::Vector2d(0.0,0.0), r, r, v, f);
    } else if (scenename == "rectangle"){
        double w = SimOptions::doubleValue("width");
        double h = SimOptions::doubleValue("height");

        gen_rectangle(N, Eigen::Vector2d(0.0,0.0), w, h, v, f);

    } else if (scenename == "ellipse"){
        double a = SimOptions::doubleValue("axis-horizontal");
        double b = SimOptions::doubleValue("axis-vertical");

        gen_ellipse(N, Eigen::Vector2d(0.0,0.0), a, b, v, f);

    } else if (scenename == "oscillation_test"){
        float theta = 2*M_PI/N;
        
        double a = 1.0/3.0;
        double eps = 0.05*a;
        int m = 2; // the mode of oscillation we are interested in

        double r;

        v.resize(N);
        f.resize(N);

        for (size_t i=0; i<N; i++){
            r = a + eps*cos(m*theta*i);

            v[i].x() = r*cos(theta*i);
            v[i].y() = r*sin(theta*i);

            f[i] = Eigen::Vector2i(i,(i+1)%N);
        }
    } else if (scenename == "semicircle_horizontal") {
        double r = SimOptions::doubleValue("radius");

        gen_semicircle_h(N, Eigen::Vector2d(0.0,0.0), r, r, v, f);

    } else if (scenename == "semicircle_vertical"){
        double r = SimOptions::doubleValue("radius");

        gen_semicircle_v(N, Eigen::Vector2d(0.0,0.0), r, r, v, f);

    } else if (scenename == "donut"){
        double r_outer = SimOptions::doubleValue("radius-outer");
        double r_inner = SimOptions::doubleValue("radius-inner");

        gen_donut(N, Eigen::Vector2d(0.0,0.0), r_inner, r_outer, v, f);

    } else if (scenename == "square_donut"){
        double d_outer = SimOptions::doubleValue("size-outer");
        double d_inner = SimOptions::doubleValue("size-inner");

        gen_square_donut(N, Eigen::Vector2d(0.0,0.0), d_inner, d_outer, v, f);

    } else if (scenename == "generic_donut"){
        double r_outer = SimOptions::doubleValue("radius-outer");

        gen_generic_donut(N, Eigen::Vector2d(0.0,0.0), r_outer, v, f);

    } else if (scenename == "generic_square_donut"){
        double d_outer = SimOptions::doubleValue("size-outer");

        gen_generic_square_donut(N, Eigen::Vector2d(0.0,0.0), d_outer, v, f);

    } else {
        std::cout<<"\""<< scenename<<"\" is not a valid scene. Please enter a valid scene!"<<std::endl;
        return;
    }

    if (N != v.size()){
        std::cout<<"Mesh was resized to have "<<v.size()<<" verts in order to better setup scene."<<std::endl;
        N = v.size();
        m.resize_mesh(N);
    }

    std::vector<Eigen::Vector2d> u(N, Eigen::Vector2d(0,0));

    std::vector<Eigen::Vector2d> v_solid(N, Eigen::Vector2d(0,0));

    std::vector<bool> air(N, true);
    std::vector<bool> solid(N, false);
    std::vector<bool> triple(N, false);
    std::vector<bool> corner(N, false);

    for (size_t i = 0; i < v.size(); i++) 
        m.verts[i] = Eigen::Vector2d (v[i].x(), v[i].y());
    for (size_t i = 0; i < f.size(); i++) 
        m.faces[i] = Eigen::Vector2i (f[i][0], f[i][1]);

    m.update_neighbor_face_vecs();
    m.reset_face_length_limits();
    
    for (size_t i = 0; i < v_solid.size(); i++) 
        m.vels_solid[i] = Eigen::Vector2d (v_solid[i][0], v_solid[i][1]);    

    m.set_boundaries(air, solid, triple, corner);
    
}

void Scenes::setupSceneShapeFromFile(LiquidMesh& m, const std::string & filename){
    LiquidMesh::loadMeshFromFile(m, filename);
    
    int N = m.verts.size();

    m.reset_boundary_types();
    m.vels = std::vector<Eigen::Vector2d>(N, Eigen::Vector2d(0,0));
    m.vels_solid = std::vector<Eigen::Vector2d>(N, Eigen::Vector2d(0,0));
    m.is_corner = std::vector<bool>(N, false);
    m.reset_face_length_limits(); // this is fairly important to set!

    // double checking
    assert(N == m.verts.size());
    assert(N == m.vels.size());
    assert(N == m.faces.size());
    assert(N == m.vertsPrevFace.size());
    assert(N == m.vertsNextFace.size());
    assert(N == m.vels_solid.size());
    assert(N == m.is_air.size());
    assert(N == m.is_solid.size());
    assert(N == m.is_triple.size());
    assert(N == m.is_corner.size());
}

void Scenes::setupSceneVelocities(std::vector<Eigen::Vector2d> & verts, std::vector<Eigen::Vector2d> & vels, const std::string & initialvelocity){
    int N = vels.size();

    std::vector<Eigen::Vector2d> u(N, Eigen::Vector2d(0,0));

    if (initialvelocity == "zero"){
        // (0,0)
        for (size_t i=0; i<N; i++){
            u[i].x() = 0;
            u[i].y() = 0;
        }
    } else if (initialvelocity == "constant"){
        // (1,0)
        for (size_t i=0; i<N; i++){
            u[i].x() = 1;
            u[i].y() = 0;
        }
    } else if (initialvelocity == "linear_x"){
        // (x,0)
        for (size_t i=0; i<N; i++){
            u[i].x() = verts[i].x();
            u[i].y() = 0;
        }
    } else if (initialvelocity == "linear_y"){
        // (0,y)
        for (size_t i=0; i<N; i++){
            u[i].x() = 0;
            u[i].y() = verts[i].y();
        }
    } else if (initialvelocity == "rotational_cw"){
        // (y,-x)
        for (size_t i=0; i<N; i++){
            u[i].x() = verts[i].y();
            u[i].y() = -verts[i].x();
        }
    } else if (initialvelocity == "rotational_ccw"){
        // (-y,x)
        for (size_t i=0; i<N; i++){
            u[i].x() = -verts[i].y();
            u[i].y() = verts[i].x();
        }
    } else if (initialvelocity == "harmonic_1"){
        // (y,x)
        for (size_t i=0; i<N; i++){
            u[i].x() = verts[i].y();
            u[i].y() = verts[i].x();
        }
    } else if (initialvelocity == "harmonic_2"){
        // (x,-y)
        for (size_t i=0; i<N; i++){
            u[i].x() = verts[i].x();
            u[i].y() = -verts[i].y();
        }
    } else if (initialvelocity == "harmonic_3"){
        // (-x^2+y^2,2xy)
        for (size_t i=0; i<N; i++){
            u[i].x() = -verts[i].x()*verts[i].x() + verts[i].y()*verts[i].y();
            u[i].y() = 2 * verts[i].x() * verts[i].y();
        }
    } else {
        std::cout<<"please enter a valid test velocity field!"<<std::endl;
        std::cout<<"defaulting to 0 velocity field!"<<std::endl;
    }

    assert(vels.size() == N);
    for (size_t i = 0; i < u.size(); i++) 
        vels[i] = Eigen::Vector2d (u[i][0], u[i][1]); 
}