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




  float dt = 1.0/240.0;
  int n = 128; // check back later why 10 on a rect is funky

  vector<Vector2d> verts(n);
  vector<Vector2i> faces(n);
  TestingHelpers::genShape("rectangle",n,verts,faces);
  vector<Vector2d> vels(n);
  TestingHelpers::generateVField("harmonic_3",n,verts,vels);
  LiquidMesh mesh(verts,faces,vels);

  Sim s(mesh, n, dt);
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
  
  // cup using solids
  /*vector<Vector2d> cupLeft = {
    Vector2d(-1.5,1.5),
    Vector2d(-1.7,1.5),
    Vector2d(-1.7,-0.7),
    Vector2d(-1.5,-0.7)
  };
  vector<Vector2d> cupBottom = {
    Vector2d(-1.7,-0.5),
    Vector2d(-1.7,-0.7),
    Vector2d(1.7,-0.7),
    Vector2d(1.7,-0.5)
  };
  vector<Vector2d> cupRight = {
    Vector2d(1.5,1.5),
    Vector2d(1.5,-0.7),
    Vector2d(1.7,-0.7),
    Vector2d(1.7,1.5)
  };
  vector<Vector2i> cupFaces = {
    Vector2i(0,1),
    Vector2i(1,2),
    Vector2i(2,3),
    Vector2i(3,0)
  };
  SolidMesh cupL(cupLeft,cupFaces);
  SolidMesh cupB(cupBottom,cupFaces);
  SolidMesh cupR(cupRight,cupFaces);
  s.addSolid(&cupL);
  s.addSolid(&cupB);
  s.addSolid(&cupR);*/
  
  // dysfunctional because I don't deal with concave shapes (yet?)...
  /*vector<Vector2d> cupverts = {
    Vector2d(1.5,1.5),
    Vector2d(1.5,-0.5),
    Vector2d(-1.5,-0.5),
    Vector2d(-1.5,1.5),
    Vector2d(-1.7,1.5),
    Vector2d(-1.7,-0.7),
    Vector2d(1.7,-0.7),
    Vector2d(1.7,1.5),
  };
  vector<Vector2i> cupfaces = {
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
  s.addSolid(&cup);*/
/*
  Vector2d testVert1 = Vector2d(-1.5,1.5);
  cout<<"Testing... winding number of "<<testVert1.x()<<","<<testVert1.y()<<": "<<cup.windingNumber(testVert1)<<endl;
  Vector2d testVert2 = Vector2d(-1.505,1.495);
  cout<<"Testing... winding number of "<<testVert2.x()<<","<<testVert2.y()<<": "<<cup.windingNumber(testVert2)<<endl;
  Vector2d testVert3 = Vector2d(-1.4,1.0);
  cout<<"Testing... winding number of "<<testVert3.x()<<","<<testVert3.y()<<": "<<cup.windingNumber(testVert3)<<endl;
*/

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
  
  /*vector<Vector2d> splitVerts = {
    Vector2d(-2,0),
    Vector2d(2,0),
    Vector2d(0.1,0.5),
    Vector2d(0.05,0.5),
    Vector2d(0,0.5),
    Vector2d(-0.05,0.5),
    Vector2d(-0.1,0.5)
  };
  vector<Vector2i> splitFaces = {
    Vector2i(0,1),
    Vector2i(1,2),
    Vector2i(2,3),
    Vector2i(3,4),
    Vector2i(4,5),
    Vector2i(5,6),
    Vector2i(6,0)
  };
  LiquidMesh splitMesh(splitVerts,splitFaces);
  cout<<"Mesh avg length: "<<splitMesh.calc_avg_face_length()<<endl;
  Sim s(splitMesh, splitVerts.size(), dt);
  */

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

  // quick and dirty solid mesh test
  /*
  vector<Vector2d> testverts = {
    Vector2d(-0.5,-0.5),
    Vector2d(0.5,-0.5),
    Vector2d(0.5,0.5),
    Vector2d(-0.5,0.5)
  };
  vector<Vector2i> testfaces = {
    Vector2i(0,1),
    Vector2i(1,2),
    Vector2i(2,3),
    Vector2i(3,0)
  };
  */
  
  //SolidMesh testSolid(testverts,testfaces);
  /*
  Vector2d testPoint = Vector2d(0.2501,-.25);
  cout<<"Is point "<<testPoint[0]<<","<<testPoint[1]<<" in our solid mesh?"<<endl;
  cout<<testSolid.checkCollisionAndSnap(testPoint)<<endl;
  cout<<testPoint[0]<<" "<<testPoint[1]<<endl;*/

  // checking out edge splitting logic...
  /*vector<Vector2d> splitVerts = {
    Vector2d(-1,0),
    Vector2d(1,0),
    Vector2d(0.5,0.2),
    Vector2d(0,0.2),
    Vector2d(-0.5,0.2)
  };
  vector<Vector2i> splitFaces = {
    Vector2i(0,1),
    Vector2i(1,2),
    Vector2i(2,3),
    Vector2i(3,4),
    Vector2i(4,0)
  };
  LiquidMesh splitMesh(splitVerts,splitFaces);
  cout<<"Mesh avg length: "<<splitMesh.calc_avg_face_length()<<endl;
  splitMesh.edge_split();
  cout<<"Verts: "<<endl;
  for(size_t i=0; i<splitMesh.verts.size(); i++){
    cout<<"("<<splitMesh.verts[i].x()<<","<<splitMesh.verts[i].y()<<")"<<endl;
  }
  cout<<"Faces: "<<endl;
  for(size_t i=0; i<splitMesh.faces.size(); i++){
    cout<<"("<<splitMesh.faces[i].x()<<","<<splitMesh.faces[i].y()<<")"<<endl;
  }
  cout<<"Verts prev face: "<<endl;
  for(size_t i=0; i<splitMesh.vertsPrevFace.size(); i++){
    cout<<i<<"->("<<splitMesh.vertsPrevFace[i]<<")"<<endl;
  }
  cout<<"Verts next face: "<<endl;
  for(size_t i=0; i<splitMesh.vertsNextFace.size(); i++){
    cout<<i<<"->("<<splitMesh.vertsNextFace[i]<<")"<<endl;
  }*/
  /*vector<Vector2d> splitVerts = {
    Vector2d(-1,0),
    Vector2d(1,0),
    Vector2d(0.5,0.2),
    Vector2d(0,0.2),
    Vector2d(-0.5,0.2)
  };
  vector<Vector2i> splitFaces = {
    Vector2i(0,1),
    Vector2i(1,2),
    Vector2i(2,3),
    Vector2i(3,4),
    Vector2i(4,0)
  };
  LiquidMesh splitMesh(splitVerts,splitFaces);
  cout<<"Mesh avg length: "<<splitMesh.calc_avg_face_length()<<endl;
  splitMesh.edge_split();
  cout<<"Verts: "<<endl;
  for(size_t i=0; i<splitMesh.verts.size(); i++){
    cout<<"("<<splitMesh.verts[i].x()<<","<<splitMesh.verts[i].y()<<")"<<endl;
  }
  cout<<"Faces: "<<endl;
  for(size_t i=0; i<splitMesh.faces.size(); i++){
    cout<<"("<<splitMesh.faces[i].x()<<","<<splitMesh.faces[i].y()<<")"<<endl;
  }
  cout<<"Verts prev face: "<<endl;
  for(size_t i=0; i<splitMesh.vertsPrevFace.size(); i++){
    cout<<i<<"->("<<splitMesh.vertsPrevFace[i]<<")"<<endl;
  }
  cout<<"Verts next face: "<<endl;
  for(size_t i=0; i<splitMesh.vertsNextFace.size(); i++){
    cout<<i<<"->("<<splitMesh.vertsNextFace[i]<<")"<<endl;
  }*/

  /*vector<Vector2d> collapseVerts = {
    Vector2d(-1,0),
    Vector2d(1,0),
    Vector2d(0,1),
    Vector2d(-0.1,1)
  };
  vector<Vector2i> collapseFaces = {
    Vector2i(0,1),
    Vector2i(1,2),
    Vector2i(2,3),
    Vector2i(3,0)
  };
  LiquidMesh collapseMesh(collapseVerts,collapseFaces);
  cout<<"Mesh avg length: "<<collapseMesh.calc_avg_face_length()<<endl;
  collapseMesh.edge_collapse();
  cout<<"Verts: "<<endl;
  for(size_t i=0; i<collapseMesh.verts.size(); i++){
    cout<<"("<<collapseMesh.verts[i].x()<<","<<collapseMesh.verts[i].y()<<")"<<endl;
  }
  cout<<"Faces: "<<endl;
  for(size_t i=0; i<collapseMesh.faces.size(); i++){
    cout<<"("<<collapseMesh.faces[i].x()<<","<<collapseMesh.faces[i].y()<<")"<<endl;
  }
  cout<<"Verts prev face: "<<endl;
  for(size_t i=0; i<collapseMesh.vertsPrevFace.size(); i++){
    cout<<i<<"->("<<collapseMesh.vertsPrevFace[i]<<")"<<endl;
  }
  cout<<"Verts next face: "<<endl;
  for(size_t i=0; i<collapseMesh.vertsNextFace.size(); i++){
    cout<<i<<"->("<<collapseMesh.vertsNextFace[i]<<")"<<endl;
  }*/

}