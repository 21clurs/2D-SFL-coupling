#include "generateshapes.h"
void generateCircle(int n, std::vector<Eigen::Vector2d>& verts, std::vector<Eigen::Vector2i>& faces){
  float theta = 2*M_PI/n;
  
  int r=1;
  Eigen::Vector2d center(0.0,0.0);

  for (uint i=0; i<n; i++){
      verts[i](0) = r*cos(theta*i) + center(0);
      verts[i](1) = r*sin(theta*i) + center(1);

      faces[i] = Eigen::Vector2i(i,(i+1)%n);
  }
}
void generateSquare(int n, std::vector<Eigen::Vector2d>& verts, std::vector<Eigen::Vector2i>& faces){
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
void generateEllipse(int n, std::vector<Eigen::Vector2d>& verts, std::vector<Eigen::Vector2i>& faces){
    float theta = 2*M_PI/n;
    
    int r=1;
    Eigen::Vector2d center(0.0,0.0);

    for (uint i=0; i<n; i++){
        verts[i](0) = 2*r*cos(theta*i) + center(0);
        verts[i](1) = r*sin(theta*i) + center(1);

        faces[i] = Eigen::Vector2i(i,(i+1)%n);
    }
}
