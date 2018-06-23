//*****************************************************
//*** AsFem: A simple finite element method program ***
//*** Author: Yang Bai                              ***
//*** Bug report: walkandthinker@gmail.com          ***
//*** QQ group: 797998860                           ***
//*****************************************************
//*** Define the base class for mesh                ***
//*** All the mesh class should inherit from this   ***
//*****************************************************

#ifndef ASFEM_MESHBASE_H
#define ASFEM_MESHBASE_H

#include <iostream>
#include <string>
#include <vector>
#include <iomanip>


using namespace std;

class MeshBase
{
public:
    virtual void PrintMeshInfo(string str="") const =0;
    virtual void PrintMeshDetailInfo(string str="") const=0;

    virtual bool CreateMesh() =0;
    virtual int GetNodesNum() const =0;
    virtual int GetElmtsNum() const =0;
    virtual int GetNodesNumPerElmt() const =0;
    virtual int GetVTKCellType() const =0;
    virtual string GetElmtType() const =0;

    virtual double IthNodeJthCoords(int i,int j) const =0;
    virtual int IthConnJthIndex(int e,int j) const =0;

    int nNodes,nElmts,nNodesPerElmt,nDims;
    double Xmin,Xmax,Ymin,Ymax,Zmin,Zmax;
    int VTKCellType;
    string ElmtType;
    vector<double> NodeCoords;
    vector<int> Conn;
    bool MeshGenerated;

};


#endif //ASFEM_MESHBASE_H
