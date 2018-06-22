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

    void CreateMesh(Mesh &mesh,double xmin=0.0,double xmax=1.0,int nx=10,string elmttype="edge2");

    void CreateMesh(Mesh &mesh,
                    double xmin=0.0,double xmax=1.0,
                    double ymin=0.0,double ymax=1.0,
                    int nx=2,int ny=2,string elmttype="quad4");

    void CreateMesh(Mesh &mesh,
                    double xmin=0.0,double xmax=1.0,
                    double ymin=0.0,double ymax=1.0,
                    double zmin=0.0,double zmax=1.0,
                    int nx=2,int ny=2,int nz=2,string elmttype="hex8");

};


#endif //ASFEM_MESHBASE_H