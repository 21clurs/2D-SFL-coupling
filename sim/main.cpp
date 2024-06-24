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

int main(){
  float dt = 1.0/240.0;
  int n = 128; // check back later why 10 on a rect is funky

  std::vector<Vector2d> verts(n);
  std::vector<Vector2i> faces(n);
  TestingHelpers::genShape("rectangle",n,verts,faces);
  std::vector<Vector2d> vels(n);
  TestingHelpers::generateVField("zero",n,verts,vels);
  LiquidMesh mesh(verts,faces,vels);

  Sim s(mesh, n, dt);
  
  /*
  std::cout<<"CIRCLE HHD"<<std::endl;
  TestingHelpers::testHHDErrorTables("circle","constant");
  TestingHelpers::testHHDErrorTables("circle","harmonic_1");
  TestingHelpers::testHHDErrorTables("circle","harmonic_3");

  std::cout<<"SEMICIRCLE HHD"<<std::endl;
  TestingHelpers::testHHDErrorTables("semicircle_h","constant");
  TestingHelpers::testHHDErrorTables("semicircle_h","harmonic_1");
  TestingHelpers::testHHDErrorTables("semicircle_h","harmonic_3");
  
  std::cout<<"CIRCLE BEM"<<std::endl;
  TestingHelpers::testBEMErrorTables("circle","1");
  TestingHelpers::testBEMErrorTables("circle","2");
  TestingHelpers::testBEMErrorTables("circle","3");

  std::cout<<"SEMICIRCLE BEM"<<std::endl;
  TestingHelpers::testBEMErrorTables("semicircle_h","1");
  TestingHelpers::testBEMErrorTables("semicircle_h","2");
  TestingHelpers::testBEMErrorTables("semicircle_h","3");
  */

  /*std::vector<Vector2d> testverts = {
    Vector2d(-3.0,-1.0),
    Vector2d(3.0,-1.0),
    Vector2d(3.0,0.0),
    Vector2d(-3.0,0.0)
  };
  std::vector<Vector2i> testfaces = {
    Vector2i(0,1),
    Vector2i(1,2),
    Vector2i(2,3),
    Vector2i(3,0)
  };
  SolidMesh testSolid(testverts,testfaces);
  testSolid.setVelFunc([](double t)->Vector2d{ return Vector2d(0, (1.0/4.0)*cos(t*2.0-M_PI_2)); });
  s.addSolid(&testSolid);*/
  
  // cup using solids
  /*std::vector<Vector2d> cupLeft = {
    Vector2d(-1.5,1.5),
    Vector2d(-1.7,1.5),
    Vector2d(-1.7,-0.7),
    Vector2d(-1.5,-0.7)
  };
  std::vector<Vector2d> cupBottom = {
    Vector2d(-1.7,-0.5),
    Vector2d(-1.7,-0.7),
    Vector2d(1.7,-0.7),
    Vector2d(1.7,-0.5)
  };
  std::vector<Vector2d> cupRight = {
    Vector2d(1.5,1.5),
    Vector2d(1.5,-0.7),
    Vector2d(1.7,-0.7),
    Vector2d(1.7,1.5)
  };
  std::vector<Vector2i> cupFaces = {
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
/*
  Vector2d testVert1 = Vector2d(-1.5,1.5);
  std::cout<<"Testing... winding number of "<<testVert1.x()<<","<<testVert1.y()<<": "<<cup.windingNumber(testVert1)<<std::endl;
  Vector2d testVert2 = Vector2d(-1.505,1.495);
  std::cout<<"Testing... winding number of "<<testVert2.x()<<","<<testVert2.y()<<": "<<cup.windingNumber(testVert2)<<std::endl;
  Vector2d testVert3 = Vector2d(-1.4,1.0);
  std::cout<<"Testing... winding number of "<<testVert3.x()<<","<<testVert3.y()<<": "<<cup.windingNumber(testVert3)<<std::endl;
*/

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
  
  /*std::vector<Vector2d> splitVerts = {
    Vector2d(-2,0),
    Vector2d(2,0),
    Vector2d(0.1,0.5),
    Vector2d(0.05,0.5),
    Vector2d(0,0.5),
    Vector2d(-0.05,0.5),
    Vector2d(-0.1,0.5)
  };
  std::vector<Vector2i> splitFaces = {
    Vector2i(0,1),
    Vector2i(1,2),
    Vector2i(2,3),
    Vector2i(3,4),
    Vector2i(4,5),
    Vector2i(5,6),
    Vector2i(6,0)
  };
  LiquidMesh splitMesh(splitVerts,splitFaces);
  std::cout<<"Mesh avg length: "<<splitMesh.calc_avg_face_length()<<std::endl;
  Sim s(splitMesh, splitVerts.size(), dt);
  */

  // moving surface
  /*Vector2d a(-2.0,0);
  Vector2d b(2.0,0);
  Object w(dt, b, a);
  w.setVelFunc([](double t)->Vector2d{ return Vector2d(0, (1.0/4.0)*cos(t*2.0-M_PI_2)); });
  s.addWall(&w);*/
  
  
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
      testSolid.setVelFunc([](double t)->Vector2d{ return Vector2d(0, -(1.0/4.0)*sin(t*4.0)); });
    } 
    if (i==344){
      //testSolid.setVelFunc([](double t)->Vector2d{ return Vector2d(0, 0); });
    }
  }
  std::cout << std::endl;

  // quick and dirty rigid body tests
  std::vector<Vector2d> squareverts = {
    Vector2d(0,0),
    Vector2d(1,0),
    Vector2d(1,1),
    Vector2d(0,1),
    /*Vector2d(0.5,0.5),
    Vector2d(0.8,0.5),
    Vector2d(0.8,0.2),
    Vector2d(0.5,0.2)*/
  };
  std::vector<Vector2i> squarefaces = {
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
  std::cout<<"area: "<<testRB.area<<std::endl;
  std::cout<<"COM: "<<testRB.com[0]<<","<<testRB.com[1]<<std::endl;
  std::cout<<"MOI: "<<testRB.moi<<std::endl;
  */

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

  // checking out edge splitting logic...
  /*std::vector<Vector2d> splitVerts = {
    Vector2d(-1,0),
    Vector2d(1,0),
    Vector2d(0.5,0.2),
    Vector2d(0,0.2),
    Vector2d(-0.5,0.2)
  };
  std::vector<Vector2i> splitFaces = {
    Vector2i(0,1),
    Vector2i(1,2),
    Vector2i(2,3),
    Vector2i(3,4),
    Vector2i(4,0)
  };
  LiquidMesh splitMesh(splitVerts,splitFaces);
  std::cout<<"Mesh avg length: "<<splitMesh.calc_avg_face_length()<<std::endl;
  splitMesh.edge_split();
  std::cout<<"Verts: "<<std::endl;
  for(size_t i=0; i<splitMesh.verts.size(); i++){
    std::cout<<"("<<splitMesh.verts[i].x()<<","<<splitMesh.verts[i].y()<<")"<<std::endl;
  }
  std::cout<<"Faces: "<<std::endl;
  for(size_t i=0; i<splitMesh.faces.size(); i++){
    std::cout<<"("<<splitMesh.faces[i].x()<<","<<splitMesh.faces[i].y()<<")"<<std::endl;
  }
  std::cout<<"Verts prev face: "<<std::endl;
  for(size_t i=0; i<splitMesh.vertsPrevFace.size(); i++){
    std::cout<<i<<"->("<<splitMesh.vertsPrevFace[i]<<")"<<std::endl;
  }
  std::cout<<"Verts next face: "<<std::endl;
  for(size_t i=0; i<splitMesh.vertsNextFace.size(); i++){
    std::cout<<i<<"->("<<splitMesh.vertsNextFace[i]<<")"<<std::endl;
  }*/
  /*std::vector<Vector2d> splitVerts = {
    Vector2d(-1,0),
    Vector2d(1,0),
    Vector2d(0.5,0.2),
    Vector2d(0,0.2),
    Vector2d(-0.5,0.2)
  };
  std::vector<Vector2i> splitFaces = {
    Vector2i(0,1),
    Vector2i(1,2),
    Vector2i(2,3),
    Vector2i(3,4),
    Vector2i(4,0)
  };
  LiquidMesh splitMesh(splitVerts,splitFaces);
  std::cout<<"Mesh avg length: "<<splitMesh.calc_avg_face_length()<<std::endl;
  splitMesh.edge_split();
  std::cout<<"Verts: "<<std::endl;
  for(size_t i=0; i<splitMesh.verts.size(); i++){
    std::cout<<"("<<splitMesh.verts[i].x()<<","<<splitMesh.verts[i].y()<<")"<<std::endl;
  }
  std::cout<<"Faces: "<<std::endl;
  for(size_t i=0; i<splitMesh.faces.size(); i++){
    std::cout<<"("<<splitMesh.faces[i].x()<<","<<splitMesh.faces[i].y()<<")"<<std::endl;
  }
  std::cout<<"Verts prev face: "<<std::endl;
  for(size_t i=0; i<splitMesh.vertsPrevFace.size(); i++){
    std::cout<<i<<"->("<<splitMesh.vertsPrevFace[i]<<")"<<std::endl;
  }
  std::cout<<"Verts next face: "<<std::endl;
  for(size_t i=0; i<splitMesh.vertsNextFace.size(); i++){
    std::cout<<i<<"->("<<splitMesh.vertsNextFace[i]<<")"<<std::endl;
  }*/

  /*std::vector<Vector2d> collapseVerts = {
    Vector2d(-1,0),
    Vector2d(1,0),
    Vector2d(0,1),
    Vector2d(-0.1,1)
  };
  std::vector<Vector2i> collapseFaces = {
    Vector2i(0,1),
    Vector2i(1,2),
    Vector2i(2,3),
    Vector2i(3,0)
  };
  LiquidMesh collapseMesh(collapseVerts,collapseFaces);
  std::cout<<"Mesh avg length: "<<collapseMesh.calc_avg_face_length()<<std::endl;
  collapseMesh.edge_collapse();
  std::cout<<"Verts: "<<std::endl;
  for(size_t i=0; i<collapseMesh.verts.size(); i++){
    std::cout<<"("<<collapseMesh.verts[i].x()<<","<<collapseMesh.verts[i].y()<<")"<<std::endl;
  }
  std::cout<<"Faces: "<<std::endl;
  for(size_t i=0; i<collapseMesh.faces.size(); i++){
    std::cout<<"("<<collapseMesh.faces[i].x()<<","<<collapseMesh.faces[i].y()<<")"<<std::endl;
  }
  std::cout<<"Verts prev face: "<<std::endl;
  for(size_t i=0; i<collapseMesh.vertsPrevFace.size(); i++){
    std::cout<<i<<"->("<<collapseMesh.vertsPrevFace[i]<<")"<<std::endl;
  }
  std::cout<<"Verts next face: "<<std::endl;
  for(size_t i=0; i<collapseMesh.vertsNextFace.size(); i++){
    std::cout<<i<<"->("<<collapseMesh.vertsNextFace[i]<<")"<<std::endl;
  }*/

}