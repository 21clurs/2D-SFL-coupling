#include <fstream>
#include <iostream>
#include <Eigen/Dense>
#include <math.h>

#include "liquidmesh.h"
#include "sim.h"
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

  // try to open and load input file
  string inputFileName = argv[1];
  Sim::setAndLoadSimOptions(inputFileName);

  Sim s_noparam;
  s_noparam.run();

}