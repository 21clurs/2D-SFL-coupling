#ifndef TESTING_HELPERS_H
#define TESTING_HELPERS_H

#include <math.h>
#include <vector>
#include <algorithm> 
#include <Eigen/Dense>
#include <iostream>
#include "liquidmesh.h"
#include "sim.h"

class TestingHelpers
{
    public:
        static void generateVField(std::string fieldType, int n, std::vector<Eigen::Vector2d>& verts, std::vector<Eigen::Vector2d>& vels);
        static void testHHDErrorTables(std::string shape, std::string fieldType);
        static void testBEMErrorTables(std::string shape, std::string fieldType);
        static void testHHD(std::string shape, std::string fieldType, int n, Eigen::Vector2d& errs);
        static void testBEM(std::string shape, std::string fieldType, int n, std::vector<Eigen::Vector2d>& errs);
        static void genShape(std::string shape, int n, std::vector<Eigen::Vector2d>& verts, std::vector<Eigen::Vector2i>& faces);
        static void setFunBC(std::string shape, int n, std::vector<bool>& is_air, std::vector<bool>& is_solid, std::vector<bool>& is_triple);
};
#endif