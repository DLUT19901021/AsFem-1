//*****************************************************
//*** AsFem: A simple finite element method program ***
//*** Author: Yang Bai                              ***
//*** Bug report: walkandthinker@gmail.com          ***
//*** QQ group: 797998860                           ***
//*****************************************************
//*** Define 2D mesh                                ***
//***  inherit from mesh base                       ***
//*****************************************************

#ifndef ASFEM_MESH2D_H
#define ASFEM_MESH2D_H

#include "StringUtils/StringUtils.h"
#include "MeshBase.h"
#include "petsc.h"

class Mesh2D:public MeshBase
{
public:
    Mesh2D(double xmin=0.0,double xmax=1.0,double ymin=0.0,double ymax=1.0,int nx=2,int ny=2,string elmttype="quad4");

    virtual void PrintMeshInfo(string str="") const override ;
    virtual void PrintMeshDetailInfo(string str="") const override ;

    virtual bool CreateMesh() override ;
    virtual int GetNodesNum() const override { return nNodes;}
    virtual int GetElmtsNum() const override { return nElmts;}
    virtual int GetNodesNumPerElmt() const override { return nNodesPerElmt;}
    virtual int GetVTKCellType() const override { return VTKCellType;}
    virtual string GetElmtType() const override { return ElmtType;}

    virtual double IthNodeJthCoords(int i,int j) const override ;
    virtual int IthConnJthIndex(int e,int j) const override ;

    void Release();

private:
    double Xmin,Xmax,Ymin,Ymax;
    int Nx,Ny;
    int P;

    // So, BoundaryNode=Left+Right+Bottom+Top
    int *BoundaryNodeIndex,nBoundaryNodeIndex;

    int *LeftEdgeNodeIndex,nLeftEdgeNodeIndex;
    int *RightEdgeNodeIndex,nRightEdgeNodeIndex;

    int *BottomEdgeNodeIndex,nBottomEdgeNodeIndex;
    int *TopEdgeNodeIndex,nTopEdgeNodeIndex;
};


#endif //ASFEM_MESH2D_H
