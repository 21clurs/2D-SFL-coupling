#include <iostream>
#include <Eigen/Dense>
#include <math.h>

#include "sim.h"
#include "mesh.h"
#include "generateshapes.h"

using Eigen::MatrixXd;
using Eigen::Vector2d;
using Eigen::Vector2i;

int main(){

  //Mesh a;
  MatrixXd m(2,2);
  m(0,0) = 3;
  m(1,0) = 2.5;
  m(0,1) = -1;
  m(1,1) = m(1,0) + m(0,1);
  //std::cout << m << std::endl;

  int n = 32;
  std::vector<Vector2d> verts = std::vector<Vector2d>(n,Vector2d(0.0,0.0));
  std::vector<Vector2i> faces = std::vector<Vector2i>(n,Vector2i(0,0));
  std::vector<Vector2d> vels = std::vector<Vector2d>(n,Vector2d(1.0,0.0));
  GenerateShape::circle(n,verts,faces);
  Mesh circMesh(verts,faces,vels);
  for (int i=0; i<n; i++){
    //std::cout << circMesh.next_neighbor(i)<< std::endl;
    //std::cout << circMesh.vels[i]<< std::endl;
  }

  float dt = 1.0/10.0;

  Sim s(circMesh, n, dt);

  //s.outputFrame(circMesh,"tester.txt");

  for (int i=0; i<50; i++){
    s.step_sim(i);
    s.outputFrame(std::to_string(i)+".txt");
  }

  return 0;
}
