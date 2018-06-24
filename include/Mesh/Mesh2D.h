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
#include "Mesh1D.h"
#include "petsc.h"

class Mesh1D;

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

    double GetXmin() const {return Xmin;}
    double GetXmax() const {return Xmax;}
    double GetYmin() const {return Ymin;}
    double GetYmax() const {return Ymax;}

    int GetBCNodesNumPerElmt() const { return nNodesPerBCElmt;}
    int GetBCNodesNum() const { return nBCNodes;}
    int GetBCElmtsNum() const { return nBCElmts;}

    int GetSideNodesNum(string sidename) const;
    int GetSideElmtsNum(string sidename) const;
    vector<int> GetIthBCElmtConn(string sidename,int e) const;


    pair<string,vector<int>> GetSideSet(string sidename) const;
    vector<pair<string,vector<int>>> GetBoundaryElmtSet() const { return BoundaryElmtSet;}

    bool IsMeshGenerated() const {return MeshGenerated;}
    
    void Release();

private:
    double Xmin,Xmax,Ymin,Ymax;
    int Nx,Ny;
    int P;

    void SplitBoundaryMesh();
    pair<string,vector<int>> LeftBCElmtSet,RightBCElmtSet;
    pair<string,vector<int>> BottomBCElmtSet,TopBCElmtSet;
    vector<pair<string,vector<int>> > BoundaryElmtSet;
    int nNodesPerBCElmt,nBCNodes,nBCElmts;

};


#endif //ASFEM_MESH2D_H
