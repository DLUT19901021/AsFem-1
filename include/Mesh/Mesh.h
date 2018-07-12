//*************************************************************
//* This file is part of the AsFem framework                  *
//* https://github.com/walkandthinker/AsFem                   *
//* All rights reserved, see COPYRIGHT for full restrictions  *
//* Licensed under LGPL 3.0, please see LICENSE for details   *
//* https://www.gnu.org/licenses/gpl-3.0.html                 *
//*************************************************************
//*** AsFem: A simple finite element method program         ***
//*** Author: Yang Bai @CopyRight                           ***
//*** Bug report: walkandthinker@gmail.com                  ***
//*** QQ group: 797998860                                   ***
//*************************************************************
//*** Created by Y. Bai on 12.06.18.                        ***
//*** Define the mesh class of mesh                         ***
//*************************************************************

#ifndef ASFEM_MESH_H
#define ASFEM_MESH_H

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#include "petsc.h"

using namespace std;

class Mesh
{
public:
    Mesh();

    void SetDim(int dim);
    void SetXmin(double xmin);
    void SetXmax(double xmax);
    void SetYmin(double ymin);
    void SetYmax(double ymax);
    void SetZmin(double zmin);
    void SetZmax(double zmax);

    void SetNx(int nx);
    void SetNy(int ny);
    void SetNz(int nz);

    void SetMeshType(string meshtype);

    bool IsRequiredInfoComplete();

    //* For mesh generation---> built-in mesh generation

    void CreateMesh();

    // TODO: implement gmsh and abaqus mesh support!!!
    void ImportGmsh(string filename);
    void ImportAbaqusInp(string filename);;
    void ImportNetgen(string filename);

    int GetIthConnJthIndex(int i, int j) const;
    double GetIthNodeJthCoord(int i, int j) const;
    void PrintMeshShortInfo() const;
    void PrintMeshDetailedInfo() const;

    //* For boundary mesh information
    int GetSideBCIthConnJthIndex(string sidename,int i,int j) const;
    int GetBCIthConnJthIndex(int i,int j);
    double GetSideBCIthNodeJthCoord(string sidename,int i,int j) const;
    double GetBCIthNodeJthCoord(int i,int j) const;

    //* For VTK information
    int GetVTKCellType() const { return VTKCellType;}

private:
    void Create1DMesh();
    void Create2DMesh();
    void Create3DMesh();

private:
    int nDims,Nx,Ny,Nz;
    string MeshType;
    bool IsDimSet,IsNxSet,IsNySet,IsNzSet;
    bool IsMeshTypeSet;

private:
    // For mesh information
    bool IsMeshCreated;
    double Xmin,Xmax,Ymin,Ymax,Zmin,Zmax;
    int IsXminSet,IsXmaxSet,IsYminSet,IsYmaxSet,IsZminSet,IsZmaxSet;
    vector<int> Conn;
    vector<double> NodeCoords;// should be Nx4
    int nNodes,nElmts,nNodesPerElmt;

    // For boundary mesh information
    //* left---> for x=xmin boundary
    //* right--> for x=xmax boundary
    //* bottom-> for y=ymin boundary
    //* top----> for y=ymax boundary
    //* back---> for z=zmin boundary
    //* front--> for z=zmax boundary
    bool IsBCMeshSplit;
    vector<int> BCConn;
    vector<int> LeftBCConn,RightBCConn;
    vector<int> BottomBCConn,TopBCConn;
    vector<int> BackBCConn,FrontBCConn;
    int nNodesPerBCElmt;

private:
    // For VTK information
    int VTKCellType;


};

#endif //ASFEM_MESH_H
