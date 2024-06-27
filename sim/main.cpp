#include <iostream>
#include <Eigen/Dense>
#include <math.h>

#include "liquidmesh.h"
#include "sim.h"
#include "solidmesh.h"
#include "testinghelpers.h"
#include "rigidbody.h"
#include "scenes.h"


using Eigen::MatrixXd;
using Eigen::Vector2d;
using Eigen::Vector2i;

int main(){
  float dt = 1.0/240.0;
  int n = 128; // check back later why 10 on a rect is funky

  std::vector<Vector2d> verts(n);
  std::vector<Vector2i> faces(n);
  std::vector<Vector2d> vels(n);
  LiquidMesh mesh(verts,faces,vels);
  Sim s(mesh, n, dt);

  Scenes::scene(s, "oscillation_test");
/*
  // dysfunctional because I don't deal with concave shapes (yet?)...
  std::vector<Vector2d> cupverts = {
    Vector2d(1.5,1.5),
    Vector2d(1.5,-0.5),
    Vector2d(-1.5,-0.5),
    Vector2d(-1.5,1.5),
    Vector2d(-1.7,1.5),
    Vector2d(-1.7,-0.7),
    Vector2d(1.7,-0.7),
    Vector2d(1.7,1.5),
  };
  std::vector<Vector2i> cupfaces = {
    Vector2i(0,1),
    Vector2i(1,2),
    Vector2i(2,3),
    Vector2i(3,4),
    Vector2i(4,5),
    Vector2i(5,6),
    Vector2i(6,7),
    Vector2i(7,0),
  };
  SolidMesh cup(cupverts,cupfaces);
  s.addSolid(&cup);
*/

/*
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
*/
  s.collide(); // snap everything/set boundary flags properly, etc.
  int num_frames = 1000;
  for (int i=0; i<num_frames; i++){
    // sim stuff
    s.outputFrame(std::to_string(i)+".txt");
    s.step_sim(i*dt);
    
    // progress messages
    std::cout<<"Simulation steps "<<i+1<<"/"<<num_frames<<" complete."<<"\r";
    std::cout.flush();

    if (i==30){
      //testSolid.setVelFunc([](double t)->Vector2d{ return Vector2d(0, -(1.0/4.0)*sin(t*4.0)); });
    } 
    if (i==344){
      //testSolid.setVelFunc([](double t)->Vector2d{ return Vector2d(0, 0); });
    }
  }
  std::cout << std::endl;

  return 0;
}