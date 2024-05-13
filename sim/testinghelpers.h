#ifndef GENERATE_FIELDS
#define GENERATE_FIELDS

#include <math.h>
#include <vector>
#include <Eigen/Dense>
#include <iostream>
#include "mesh.h"
#include "sim.h"

class TestingHelpers
{
    public:
        static void generateVField(std::string fieldType, int n, std::vector<Eigen::Vector2d>& verts, std::vector<Eigen::Vector2d>& vels);
        static void testErrorTables(std::string testType, std::string shape, std::string fieldType);
        static void testHHD(std::string shape, std::string fieldType, int n, Eigen::Vector2d& errs);
        static void testBEM(std::string shape, std::string fieldType, int n, Eigen::Vector2d& errs);
        static void genShape(std::string shape, int n, std::vector<Eigen::Vector2d>& verts, std::vector<Eigen::Vector2i>& faces);
};
#endif