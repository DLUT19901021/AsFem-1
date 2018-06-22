//*****************************************************
//*** AsFem: A simple finite element method program ***
//*** Author: Yang Bai                              ***
//*** Bug report: walkandthinker@gmail.com          ***
//*** QQ group: 797998860                           ***
//*****************************************************
//*** Define final mesh                             ***
//***  inherit from mesh base                       ***
//*****************************************************

#include "Mesh/Mesh.h"

Mesh::Mesh(int dim,int nnodesperelmt,int ndofspernode)
{
    if(dim<1||dim>3)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: dim=%2d is invalid for a FEM mesh!!!\n",dim);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }

    if(nnodesperelmt<2||nnodesperelmt>27)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: nNodesPerElmt=%d is invalid for a FEM mesh!!!\n",nnodesperelmt);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }

    if(ndofspernode<1||ndofspernode>10)
    {
        // Currently, only 10 dofs is supported for each node
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: nDofsPerNode=%d is invalid in current version!!!\n",ndofspernode);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        max dofs per node<=10!!!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }

    Has1DMesh=false;
    Has2DMesh=false;
    IsInit=false;
    nDim=dim;nNodes=0;nElmts=0;
    nMaxNodesPerElmt=nnodesperelmt;
    nDofsPerNode=ndofspernode;

    Mesh1DList.clear();
    Mesh2DList.clear();

    Conn.clear();
    NodeCoords.clear();

}