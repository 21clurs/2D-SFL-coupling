#include "scenes.h"
#include "simoptions.h"

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
        // topdelta = w/nHorizontal;
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
}

void Scenes::scene(Sim &sim, const std::string & scenename){
    int N = sim.n;

    std::vector<Eigen::Vector2d> v;
    std::vector<Eigen::Vector2i> f;
    std::vector<Eigen::Vector2d> u(N, Eigen::Vector2d(0,0));

    if (scenename == "circle"){

        double r = SimOptions::doubleValue("radius");

        gen_ellipse(N, Eigen::Vector2d(0.0,0.0), r, r, v, f);
        

    } else if (scenename == "rectangle"){
        double w = SimOptions::doubleValue("width");
        double h = SimOptions::doubleValue("height");

        gen_rectangle(N, Eigen::Vector2d(0.0,0.0), w, h, v, f);

    } else if (scenename == "ellipse"){
        double a = SimOptions::doubleValue("axis_horizontal");
        double b = SimOptions::doubleValue("axis_vertical");

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
    } else if (scenename == "droplet_on_surface_horizontal") {
        double r = SimOptions::doubleValue("radius");

    } else if (scenename == "droplet_on_surface_vertical"){
        double r = SimOptions::doubleValue("radius");

    } else if (scenename == "donut"){
        double r_outer = SimOptions::doubleValue("radius_outer");
        double r_inner = SimOptions::doubleValue("radius_inner");

    } else if (scenename == "donut_square"){
        double d_outer = SimOptions::doubleValue("size_outer");
        double d_inner = SimOptions::doubleValue("size_inner");

    } else if (scenename == "cup"){
        double w = SimOptions::doubleValue("width");
        double h = SimOptions::doubleValue("height");

    } else if (scenename == "cup_with_block"){
        double w_cup = SimOptions::doubleValue("cup_width");
        double h_cup = SimOptions::doubleValue("cup_height");
        double w_block = SimOptions::doubleValue("block_width");
        double h_block = SimOptions::doubleValue("block_height");

    } else {

    }

    sim.m.verts.resize(v.size());
    for (size_t i = 0; i < v.size(); i++) 
        sim.m.verts[i] = Eigen::Vector2d (v[i].x(), v[i].y());
    sim.m.faces.resize(f.size()); 
    for (size_t i = 0; i < f.size(); i++) 
        sim.m.faces[i] = Eigen::Vector2i (f[i][0], f[i][1]);
    sim.m.vels.resize(u.size());
    for (size_t i = 0; i < u.size(); i++) 
        sim.m.vels[i] = Eigen::Vector2d (u[i][0], u[i][1]);
    sim.m.reset_face_length_limits();

}

