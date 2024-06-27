// scenes.h
// 06/08/2024

#ifndef SCENES_H
#define SCENES_H

#include <iostream>
#include "sim.h"

class Scenes{
    public:
        static void scene(Sim &sim, const std::string & scenename);
};

#endif