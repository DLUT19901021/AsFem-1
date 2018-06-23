//*****************************************************
//*** AsFem: A simple finite element method program ***
//*** Author: Yang Bai                              ***
//*** Bug report: walkandthinker@gmail.com          ***
//*** QQ group: 797998860                           ***
//*****************************************************
//*** This class handle the dofs of asfem           ***
//***  i.e. local to global dof map and so on       ***
//*****************************************************

#ifndef ASFEM_DOFHANDLER_H
#define ASFEM_DOFHANDLER_H

#include<iostream>
#include<vector>

#include "petsc.h"

// AsFem's own header file
#include "Mesh/Mesh.h"

class DofHandler
{
public:
    DofHandler();

    bool CreateLocalToGlobalDofMap(Mesh &mesh);
    void GetLocalDofMap(const int e,int &ndofs,int *rInd,int *cInd);

    void Release();

    void PrintDofMap() const;
private:
    bool HasDofMap=false;
    int nDofs,nNodes,nElmts;
    int nDofsPerNode,nNodesPerElmts;
    int nDofsPerElmt;

    vector<int> GlobalDofMap;

};


#endif