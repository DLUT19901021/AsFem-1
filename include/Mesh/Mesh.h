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

    int GetDims() const { return nDims;}
    int Get1DMeshNum() const { return Mesh1DList.size();}
    int Get2DMeshNum() const { return Mesh2DList.size();}

    int GetDofsNum() const { return nDofs;}
    int GetNodesNum() const { return nNodes;}
    int GetElmtsNum() const { return nElmts;}

    int GetNodesNumPerBCElmt() const { return nNodesPerBCElmt;}
    int GetNodesNumPerElmt() const { return nMaxNodesPerElmt;}
    int GetDofsNumPerNode() const { return nDofsPerNode;}

    int GetSideNodesNum(string sidename) const;
    int GetSideElmtsNum(string sidename) const;
    vector<int> GetIthBCElmtConn(string sidename,int e) const;

    pair<string,vector<int>> GetSideSet(string sidename) const;
    vector<pair<string,vector<int>>> GetBoundaryElmtSet() const { return BoundaryElmtSet;}

    double IthNodeJthCoords(int i,int j) const;
    int IthConnJthIndex(int e,int j) const;

    int GetIthVTKCellType() const {return VTKCellType;}

    void Init();
    void Add1DMesh(Mesh1D &mesh1d);
    void Add2DMesh(Mesh2D &mesh2d);

    void PrintMeshInfo(string str="") const;
    void PrintMeshDetailInfo(string str="") const;


private:
    bool Has1DMesh=false,Has2DMesh=false;
    vector<Mesh1D> Mesh1DList;
    vector<Mesh2D> Mesh2DList;
    int VTKCellType;

    bool IsInit=false;
    double Xmin,Xmax,Ymin,Ymax,Zmin,Zmax;
    int Nx,Ny,Nz;
    int nDims,nNodes,nElmts,nMaxNodesPerElmt;
    int nDofsPerNode,nDofs;
    vector<int> Conn;
    vector<double> NodeCoords;

    // For boundary mesh information
    pair<string,vector<int>> LeftBCElmtSet,RightBCElmtSet;
    pair<string,vector<int>> BottomBCElmtSet,TopBCElmtSet;
    pair<string,vector<int>> BackBCElmtSet,FrontBCElmtSet;
    vector<pair<string,vector<int>> > BoundaryElmtSet;
    int nNodesPerBCElmt,nBCNodes,nBCElmts;

};


#endif //ASFEM_MESH_H
