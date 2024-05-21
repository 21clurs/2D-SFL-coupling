#include <iostream>
#include <Eigen/Dense>
#include <math.h>

#include "testinghelpers.h"
#include "sim.h"
#include "mesh.h"
#include "wallobject.h"

using Eigen::MatrixXd;
using Eigen::Vector2d;
using Eigen::Vector2i;

int main(){
  int n = 128;

  std::vector<Vector2d> verts(n);
  std::vector<Vector2i> faces(n);
  TestingHelpers::genShape("square",n,verts,faces);
  std::vector<Vector2d> vels(n);
  TestingHelpers::generateVField("zero",n,verts,vels);
  Mesh mesh(verts,faces,vels);
  for (int i=0; i<n; i++){
    //std::cout << mesh.next_neighbor(i)<< std::endl;
    //std::cout << mesh.vels[i]<< std::endl;
  }

  float dt = 1.0/240.0;

  Sim s(mesh, n, dt);


  //TestingHelpers::testHHDErrorTables("donut","1");
  
  // moving surface
  /*
  Vector2d a(-2.0,0);
  Vector2d b(2.0,0);
  WallObject w(dt, b, a);
  w.setVelFunc([](double t)->Vector2d{ return Vector2d(0, (1.0/4.0)*cos(t*2.0-M_PI_2)); });
  //s.addWall(&w);
  //s.collide(); // snap everything/set boundary flags properly, etc.
  */


  // cup for 2x2 square
  Vector2d a(-1.0,1.5);
  Vector2d b(-1.0,-1.0);
  Vector2d c(1.0,-1.0);
  Vector2d d(1.0,1.5);
  WallObject w1(dt, b, a);
  WallObject w2(dt, c, b);
  WallObject w3(dt, d, c);
  w1.setVelFunc([](double t)->Vector2d{ return Vector2d(0.0,0.0); });
  w2.setVelFunc([](double t)->Vector2d{ return Vector2d(0.0,0.0); });
  w3.setVelFunc([](double t)->Vector2d{ return Vector2d(0.0,0.0); });
  s.addWall(&w1);
  s.addWall(&w2);
  s.addWall(&w3);
  s.collide(); // snap everything/set boundary flags properly, etc.

  int num_frames = 800;
  for (int i=0; i<num_frames; i++){
    // sim stuff
    s.outputFrame(std::to_string(i)+".txt");
    s.step_sim(i*dt);
    
    // progress messages
    std::cout<<"Simulation steps "<<i+1<<"/"<<num_frames<<" complete."<<"\r";
    std::cout.flush();

    //std::cout << "FRAME: "<<i<< std::endl;
  }
  std::cout << std::endl;

  return 0;
}