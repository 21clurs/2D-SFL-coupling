#include "output.h"

#include <iostream>
#include <fstream> 

void outputFrame(Mesh &m, std::string filename, std::string filelocation){
    std::ofstream file(filelocation+filename);
    for (uint i=0; i<m.verts.size(); i++)
        file<<"v "<<m.verts[i][0]<<" "<<m.verts[i][1]<<std::endl;
    file<<std::endl;

    // outward vertex norms
    for (uint i=0; i<m.verts.size(); i++)
        file<<"vn "<<(m.calc_vertex_normal(i))[0]<<" "<<(m.calc_vertex_normal(i))[1]<<std::endl;
    file<<std::endl;
    
    // faces
    for (uint i=0; i<m.faces.size(); i++)
        file<<"f "<<m.faces[i][0]<<" "<<m.faces[i][1]<<std::endl;
    file.close();
}
