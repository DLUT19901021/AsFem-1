//*****************************************************
//*** AsFem: A simple finite element method program ***
//*** Author: Yang Bai                              ***
//*** Bug report: walkandthinker@gmail.com          ***
//*** QQ group: 797998860                           ***
//*****************************************************
//*** Define 1D mesh                                ***
//***                                               ***
//*****************************************************

#ifndef ASFEM_MESH1D_H
#define ASFEM_MESH1D_H

#include "StringUtils/StringUtils.h"
#include "MeshBase.h"
#include "petsc.h"


class Mesh1D:public MeshBase
{
public:
    Mesh1D(double xmin=0.0,double xmax=1.0,int nx=10,string elmttype="edge2");

    virtual bool CreateMesh() override ;

    virtual double IthNodeJthCoords(int i,int j) const override ;
    virtual int IthConnJthIndex(int e,int j) const override ;

    virtual void PrintMeshInfo(string str="") const override ;
    virtual void PrintMeshDetailInfo(string str="") const override ;


    virtual int GetNodesNum() const override { return nNodes;}
    virtual int GetElmtsNum() const override { return nElmts;}
    virtual int GetNodesNumPerElmt() const override { return nNodesPerElmt;}
    virtual int GetVTKCellType() const override { return VTKCellType;}
    virtual string GetElmtType() const override { return ElmtType;}

    int GetXmin() const {return Xmin;}
    int GetXmax() const {return Xmax;}
    bool IsMeshGenerated() const {return MeshGenerated;}
    
    void Release();
private:
    double Xmax,Xmin;
    int Nx,P;
    int BoundaryNodeIndex[2],nBoundaryNodeIndex;


};


#endif //ASFEM_MESH1D_H
