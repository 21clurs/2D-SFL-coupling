#include "generateshapes.h"
void GenerateShape::circle(int n, std::vector<Eigen::Vector2d>& verts, std::vector<Eigen::Vector2i>& faces){
  float theta = 2*M_PI/n;
  
  double r=1;
  Eigen::Vector2d center(0.0,0.0);

  for (uint i=0; i<n; i++){
      verts[i](0) = r*cos(theta*i) + center(0);
      verts[i](1) = r*sin(theta*i) + center(1);

      faces[i] = Eigen::Vector2i(i,(i+1)%n);
  }
}
void GenerateShape::square(int n, std::vector<Eigen::Vector2d>& verts, std::vector<Eigen::Vector2i>& faces){
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
}
void GenerateShape::ellipse(int n, std::vector<Eigen::Vector2d>& verts, std::vector<Eigen::Vector2i>& faces){
    float theta = 2*M_PI/n;
    
    double r=1;
    Eigen::Vector2d center(0.0,0.0);

    for (uint i=0; i<n; i++){
        verts[i](0) = 2*r*cos(theta*i) + center(0);
        verts[i](1) = r*sin(theta*i) + center(1);

        faces[i] = Eigen::Vector2i(i,(i+1)%n);
    }
}
void GenerateShape::oscillation(int n, std::vector<Eigen::Vector2d>& verts, std::vector<Eigen::Vector2i>& faces){
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
}
void GenerateShape::semicircle_h(int n, std::vector<Eigen::Vector2d>& verts, std::vector<Eigen::Vector2i>& faces){
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
}
void GenerateShape::semicircle_v(int n, std::vector<Eigen::Vector2d>& verts, std::vector<Eigen::Vector2i>& faces){
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

}
