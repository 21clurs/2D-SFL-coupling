#include <iostream>
#include <Eigen/Dense>
#include <math.h>

#include "liquidmesh.h"
#include "sim.h"
#include "solidmesh.h"
#include "testinghelpers.h"
#include "wallobject.h"


using Eigen::MatrixXd;
using Eigen::Vector2d;
using Eigen::Vector2i;

int main(){
  int n = 100;

  std::vector<Vector2d> verts(n);
  std::vector<Vector2i> faces(n);
  TestingHelpers::genShape("rectangle",n,verts,faces);
  std::vector<Vector2d> vels(n);
  TestingHelpers::generateVField("zero",n,verts,vels);
  LiquidMesh mesh(verts,faces,vels);
  for (int i=0; i<n; i++){
    //std::cout << mesh.next_neighbor(i)<< std::endl;
    //std::cout << mesh.vels[i]<< std::endl;
  }

  float dt = 1.0/200.0;

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


  // cup
  Vector2d a(-1.5,1.5);
  Vector2d b(-1.5,-0.5);
  Vector2d c(1.5,-0.5);
  Vector2d d(1.5,1.5);
  WallObject w1(dt, b, a);
  WallObject w2(dt, c, b);
  WallObject w3(dt, d, c);
  w1.setVelFunc([](double t)->Vector2d{ return Vector2d(0.0,0.0); });
  w2.setVelFunc([](double t)->Vector2d{ return Vector2d(0.0,0.0); });
  w3.setVelFunc([](double t)->Vector2d{ return Vector2d(0.0,0.0); });
  s.addWall(&w1);
  s.addWall(&w2);
  s.addWall(&w3);


  // a little block in my cup
  //int n_inner = 8;
  //std::vector<Vector2d> testverts(n_inner);
  //std::vector<Vector2i> testfaces(n_inner);
  //TestingHelpers::genShape("circle",n_inner,testverts,testfaces);
  std::vector<Vector2d> testverts = {
    Vector2d(-0.5,-0.05),
    Vector2d(0.5,-0.05),
    Vector2d(0.5,1.0),
    Vector2d(-0.5,1.0)
  };
  std::vector<Vector2i> testfaces = {
    Vector2i(0,1),
    Vector2i(1,2),
    Vector2i(2,3),
    Vector2i(3,0)
  };
  SolidMesh testSolid(testverts,testfaces);
  s.addSolid(&testSolid);


  s.collide(); // snap everything/set boundary flags properly, etc.

  int num_frames = 800;
  for (int i=0; i<num_frames; i++){
    // sim stuff
    s.outputFrame(std::to_string(i)+".txt");
    s.step_sim(i*dt);
    
    // progress messages
    std::cout<<"Simulation steps "<<i+1<<"/"<<num_frames<<" complete."<<"\r";
    std::cout.flush();

    if (i==30){
      testSolid.setVelFunc([](double t)->Vector2d{ return Vector2d(0, -(1.0/2.0)*sin(t*8.0)); });
    }
    //std::cout << "FRAME: "<<i<< std::endl;
  }
  std::cout << std::endl;

  // quick and dirty solid mesh test
  /*
  std::vector<Vector2d> testverts = {
    Vector2d(-0.5,-0.5),
    Vector2d(0.5,-0.5),
    Vector2d(0.5,0.5),
    Vector2d(-0.5,0.5)
  };
  std::vector<Vector2i> testfaces = {
    Vector2i(0,1),
    Vector2i(1,2),
    Vector2i(2,3),
    Vector2i(3,0)
  };
  */
  
  //SolidMesh testSolid(testverts,testfaces);
  /*
  Vector2d testPoint = Vector2d(0.2501,-.25);
  std::cout<<"Is point "<<testPoint[0]<<","<<testPoint[1]<<" in our solid mesh?"<<std::endl;
  std::cout<<testSolid.checkCollisionAndSnap(testPoint)<<std::endl;
  std::cout<<testPoint[0]<<" "<<testPoint[1]<<std::endl;*/
}