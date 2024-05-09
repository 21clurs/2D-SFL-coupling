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
  int n = 32;
  std::vector<Vector2d> verts = std::vector<Vector2d>(n,Vector2d(0.0,0.0));
  std::vector<Vector2i> faces = std::vector<Vector2i>(n,Vector2i(0,0));
  GenerateShape::circle(n,verts,faces);
  std::vector<Vector2d> vels = std::vector<Vector2d>(n,Vector2d(0.0,0.0));
  for (int i=0; i<n; i++){
    vels[i].x() = verts[i].x();
  }
  Mesh circMesh(verts,faces,vels);
  for (int i=0; i<n; i++){
    //std::cout << circMesh.next_neighbor(i)<< std::endl;
    //std::cout << circMesh.vels[i]<< std::endl;
  }

  float dt = 1.0/60.0;

  Sim s(circMesh, n, dt);

  //s.outputFrame(circMesh,"tester.txt");

  int num_frames = 30;
  for (int i=0; i<num_frames; i++){
    // sim stuff
    s.outputFrame(std::to_string(i)+".txt");
    s.step_sim(i);
    
    // progress messages
    std::cout<<"Simulation steps "<<i+1<<"/"<<num_frames<<" complete."<<"\r";
    std::cout.flush();
  }
  std::cout << std::endl;

  return 0;
}
