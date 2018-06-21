//*****************************************************
//*** AsFem: A simple finite element method program ***
//*** Author: Yang Bai                              ***
//*** Bug report: walkandthinker@gmail.com          ***
//*** QQ group: 797998860                           ***
//*****************************************************
//*** Define the base class for mesh                ***
//*** All the mesh class should inherit from this   ***
//*****************************************************

#ifndef ASFEM_MESHTOOLS_H
#define ASFEM_MESHTOOLS_H

#include <iostream>
#include <string>
#include <iomanip>
#include "petsc.h"

//
#include "StringUtils/StringUtils.h"

#include "Mesh/Mesh.h"

using namespace std;

class MeshTools
{
public:
    MeshTools();

    void CreateMesh(double xmin,double xmax,int ne,string elmttype,Mesh &mesh);

    void CreateMesh(double xmin,double xmax,
                    double ymin,double ymax,
                    int nx,int ny,string elmttype,
                    Mesh &mesh);

};


#endif //ASFEM_MESHBASE_H