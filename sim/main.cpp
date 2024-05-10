#include <iostream>
#include <Eigen/Dense>
#include <math.h>

#include "sim.h"
#include "mesh.h"
#include "generateshapes.h"
#include "generatefields.h"

using Eigen::MatrixXd;
using Eigen::Vector2d;
using Eigen::Vector2i;

void test_error_tables();
void test_HHD_error(int n, Eigen::Vector2d& errs);

int main(){
  int n = 32;
  std::vector<Vector2d> verts(n);
  std::vector<Vector2i> faces(n);
  GenerateShape::circle(n,verts,faces);
  std::vector<Vector2d> vels(n);
  GenerateField::harmonic_1(n,verts,vels);
  Mesh mesh(verts,faces,vels);
  for (int i=0; i<n; i++){
    //std::cout << mesh.next_neighbor(i)<< std::endl;
    //std::cout << mesh.vels[i]<< std::endl;
  }

  float dt = 1.0/60.0;

  Sim s(mesh, n, dt);

  //s.outputFrame(circMesh,"tester.txt");
  /*
  int num_frames = 240;
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

 test_error_tables();

  return 0;
}

void test_error_tables(){
  std::vector<size_t> N = {4, 8, 16, 32, 64, 128, 256};
  std::vector<Eigen::Vector2d> errs(N.size());
  for (size_t i=0; i<N.size(); i++){
    test_HHD_error(N[i], errs[i]);

    // progress messages
    std::cout<<"Test "<<i+1<<"/"<<N.size()<<" complete."<<"\r";
    std::cout.flush();
  }
  std::cout<<std::endl;

  std::cout<<"N"<<"\t"<<"err_mean"<<"\t"<<"err_max"<<"\t"<<"err_mean_ratios"<<"\t"<<"err_max_ratios"<<std::endl;
  std::cout<<"----------------------------------------------------------------"<<std::endl;
  for (size_t i=0; i<N.size(); i++){
      if (i>0){
        std::cout<<N[i]<<"\t"<<errs[i].x()<<"\t"<<errs[i].y()<<"\t"<<errs[i].x()/errs[i-1].x()<<"\t"<<errs[i].y()/errs[i-1].y()<<std::endl;
      } else{
        std::cout<<N[i]<<"\t"<<errs[i].x()<<"\t"<<errs[i].y()<<"\t"<<"n/a"<<"\t"<<"n/a"<<std::endl;
      }
  }
}

void test_HHD_error(int n, Eigen::Vector2d& errs){
  // set up test shape
  std::vector<Vector2d> verts(n);
  std::vector<Vector2i> faces(n);
  GenerateShape::circle(n,verts,faces);

  // set up test velocities
  // this test only makes sense when this is a harmonic velocity field
  std::vector<Vector2d> vels(n);
  GenerateField::harmonic_1(n,verts,vels);

  // set up expected/theoretical result
  std::vector<Vector2d> expected_post_HHD_vels(n);
  for (size_t i=0; i<n; i++){
    // a very basic expected velocity
    expected_post_HHD_vels[i] = vels[i];
  }

  // set up test simulation
  Mesh mesh(verts,faces,vels);
  float dt = 1.0/60.0;
  Sim s(mesh, n, dt);
  
  // step HHD
  s.step_HHD();

  // calculate errors
  Eigen::VectorXd err_norms(n);
  for (size_t i=0; i<n; i++){
    err_norms[i] = (s.get_vels()[i] - expected_post_HHD_vels[i]).norm();
  }

  // return errors...
  errs[0] = err_norms.mean();// mean error
  errs[1] = err_norms.cwiseAbs().maxCoeff();// maximum error
}