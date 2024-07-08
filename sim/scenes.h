// scenes.h
// 06/08/2024

#ifndef SCENES_H
#define SCENES_H

#include <iostream>
#include "sim.h"
#include "testinghelpers.h"

class Scenes{
    friend class TestingHelpers;
    public:
        static void scene(Sim * const &sim, const std::string & scenename, const std::string & initialvelocity);
    protected:
        static void setupSceneShape(LiquidMesh& m, const std::string & scenename);
        static void setupSceneVelocities(std::vector<Eigen::Vector2d> & verts, std::vector<Eigen::Vector2d> & vels, const std::string & initialvelocity);
};

#endif