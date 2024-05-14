#include <iostream>
#include <Eigen/Dense>
#include <math.h>

#include "testinghelpers.h"
#include "sim.h"
#include "mesh.h"

using Eigen::MatrixXd;
using Eigen::Vector2d;
using Eigen::Vector2i;

int main(){
  /*
  int n = 32;

  std::vector<Vector2d> verts(n);
  std::vector<Vector2i> faces(n);
  GenerateShape::donut(n,verts,faces);
  std::vector<Vector2d> vels(n);
  TestingHelpers::generateVField("zero",n,verts,vels);
  Mesh mesh(verts,faces,vels);
  for (int i=0; i<n; i++){
    //std::cout << mesh.next_neighbor(i)<< std::endl;
    //std::cout << mesh.vels[i]<< std::endl;
  }

  float dt = 1.0/60.0;

  Sim s(mesh, n, dt);

  //s.outputFrame(circMesh,"tester.txt");
  
  
  int num_frames = 50;
  for (int i=0; i<num_frames; i++){
    // sim stuff
    s.outputFrame(std::to_string(i)+".txt");
    s.step_sim(i);
    
    // progress messages
    std::cout<<"Simulation steps "<<i+1<<"/"<<num_frames<<" complete."<<"\r";
    std::cout.flush();
  }
  std::cout << std::endl;
  */

  TestingHelpers::testBEMErrorTables("donut","1");

  return 0;
}