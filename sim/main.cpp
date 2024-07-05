#include <fstream>
#include <iostream>
#include <Eigen/Dense>
#include <math.h>

#include "liquidmesh.h"
#include "sim.h"
#include "solidmesh.h"
#include "testinghelpers.h"
#include "rigidbody.h"


using Eigen::MatrixXd;
using Eigen::Vector2d;
using Eigen::Vector2i;

using namespace std;

int main(int argc, char **argv){
  if (argc <= 1){
    cout<<"Please input an valid filename or path."<<endl;
    return 1;
  }

  // trying to open input file
  string inputFileName = argv[1];
  Sim::setAndLoadSimOptions(inputFileName);

  Sim s_noparam;
  s_noparam.run();
  /*
  for (size_t i=0; i<s.markerparticles.size(); i++){
    Eigen::Vector2d v = s.HHD_FD(s.markerparticles[i],0.01);
    cout<<"("<<s.markerparticles[i].x()<<","<<s.markerparticles[i].y()<<"): ["<<v.x()<<","<<v.y()<<"]"<<endl;
  }
  */
  
  /*
  cout<<"CIRCLE HHD"<<endl;
  TestingHelpers::testHHDErrorTables("circle","constant");
  TestingHelpers::testHHDErrorTables("circle","harmonic_1");
  TestingHelpers::testHHDErrorTables("circle","harmonic_3");

  cout<<"SEMICIRCLE HHD"<<endl;
  TestingHelpers::testHHDErrorTables("semicircle_h","constant");
  TestingHelpers::testHHDErrorTables("semicircle_h","harmonic_1");
  TestingHelpers::testHHDErrorTables("semicircle_h","harmonic_3");
  
  cout<<"CIRCLE BEM"<<endl;
  TestingHelpers::testBEMErrorTables("circle","1");
  TestingHelpers::testBEMErrorTables("circle","2");
  TestingHelpers::testBEMErrorTables("circle","3");

  cout<<"SEMICIRCLE BEM"<<endl;
  TestingHelpers::testBEMErrorTables("semicircle_h","1");
  TestingHelpers::testBEMErrorTables("semicircle_h","2");
  TestingHelpers::testBEMErrorTables("semicircle_h","3");
  */

  /*vector<Vector2d> testverts = {
    Vector2d(-3.0,-1.0),
    Vector2d(3.0,-1.0),
    Vector2d(3.0,0.0),
    Vector2d(-3.0,0.0)
  };
  vector<Vector2i> testfaces = {
    Vector2i(0,1),
    Vector2i(1,2),
    Vector2i(2,3),
    Vector2i(3,0)
  };
  SolidMesh testSolid(testverts,testfaces);
  testSolid.setVelFunc([](double t)->Vector2d{ return Vector2d(0, (1.0/4.0)*cos(t*2.0-M_PI_2)); });
  s.addSolid(&testSolid);*/
  

  // a little block in my cup
  //int n_inner = 8;
  //vector<Vector2d> testverts(n_inner);
  //vector<Vector2i> testfaces(n_inner);
  //TestingHelpers::genShape("circle",n_inner,testverts,testfaces);
  /*vector<Vector2d> testverts = {
    Vector2d(-0.5,-0.05),
    Vector2d(0.5,-0.05),
    Vector2d(0.5,1.0),
    Vector2d(-0.5,1.0)
  };
  vector<Vector2i> testfaces = {
    Vector2i(0,1),
    Vector2i(1,2),
    Vector2i(2,3),
    Vector2i(3,0)
  };
  SolidMesh testSolid(testverts,testfaces);*/
  //s.addSolid(&testSolid);

  // moving surface
  /*Vector2d a(-2.0,0);
  Vector2d b(2.0,0);
  Object w(dt, b, a);
  w.setVelFunc([](double t)->Vector2d{ return Vector2d(0, (1.0/4.0)*cos(t*2.0-M_PI_2)); });
  s.addWall(&w);*/
  
  

  // quick and dirty rigid body tests
  vector<Vector2d> squareverts = {
    Vector2d(0,0),
    Vector2d(1,0),
    Vector2d(1,1),
    Vector2d(0,1),
    /*Vector2d(0.5,0.5),
    Vector2d(0.8,0.5),
    Vector2d(0.8,0.2),
    Vector2d(0.5,0.2)*/
  };
  vector<Vector2i> squarefaces = {
    Vector2i(0,1),
    Vector2i(1,2),
    Vector2i(2,3),
    Vector2i(3,0),
    /*Vector2i(4,5),
    Vector2i(5,6),
    Vector2i(6,7),
    Vector2i(7,4)*/
  };
  /*
  RigidBody testRB(squareverts,squarefaces);
  cout<<"area: "<<testRB.area<<endl;
  cout<<"COM: "<<testRB.com[0]<<","<<testRB.com[1]<<endl;
  cout<<"MOI: "<<testRB.moi<<endl;
  */
}