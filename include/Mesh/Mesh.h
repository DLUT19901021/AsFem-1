//*****************************************************
//*** AsFem: A simple finite element method program ***
//*** Author: Yang Bai                              ***
//*** Bug report: walkandthinker@gmail.com          ***
//*** QQ group: 797998860                           ***
//*****************************************************
//*** Define final mesh                             ***
//***  inherit from mesh base                       ***
//*****************************************************

#ifndef ASFEM_MESH_H
#define ASFEM_MESH_H

#include<iostream>
#include<vector>

#include "petsc.h"

#include "Mesh/Mesh1D.h"
#include "Mesh/Mesh2D.h"

using namespace std;

class Mesh
{
public:
    Mesh(int dim=3,int nnodesperelmt=27,int ndofspernode=1);
    void Release();

    int Get1DMeshNum() const { return Mesh1DList.size();}
    int Get2DMeshNum() const { return Mesh2DList.size();}

    double IthNodeJthCoords(int i,int j) const;
    int IthConnJthIndex(int e,int j) const;

    int GetIthVTKCellType() const;

    void Init();
    void Add1DMesh(Mesh1D &mesh1d);
    void Add2DMesh(Mesh2D &mesh2d);

    void PrintMeshInfo(string str="") const;
    void PrintMeshDetailInfo(string str="") const;


private:
    bool Has1DMesh=false,Has2DMesh=false;
    vector<Mesh1D> Mesh1DList;
    vector<Mesh2D> Mesh2DList;

    bool IsInit=false;
    double Xmin,Xmax,Ymin,Ymax,Zmin,Zmax;
    int Nx,Ny,Nz;
    int nDim,nNodes,nElmts,nMaxNodesPerElmt;
    int nDofsPerNode;
    vector<int> Conn;
    vector<double> NodeCoords;

};


#endif //ASFEM_MESH_H
