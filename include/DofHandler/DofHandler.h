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
#include <set>

#include "petsc.h"

// AsFem's own header file
#include "Mesh/Mesh.h"

class DofHandler
{
public:
    DofHandler();

    bool CreateLocalToGlobalDofMap(Mesh &mesh,int ndofspernode);
    void GetLocalDofMap(const int e,int &ndofs,int (&rInd)[500],int (&cInd)[500]) const;
    void GetLocalBCDofMap(string sidename,const int e,int &ndofs,int (&Ind)[500]) const;

    void Release();

    void PrintDofMap() const;
private:
    bool HasDofMap=false;
    bool HasBCDofMap=false;
    int nDofs,nNodes,nElmts;
    int nDofsPerNode,nNodesPerElmts;
    int nDofsPerElmt;

    vector<int> GlobalDofMap;

    pair<string,vector<int>> LeftSideDofMap,RightSideDofMap;
    pair<string,vector<int>> BottomSideDofMap,TopSideDofMap;
    vector< pair<string,vector<int>> > GlobalBCDofMap;

    int nDofsPerBCElmt,nNodesPerBCElmt;

};


#endif