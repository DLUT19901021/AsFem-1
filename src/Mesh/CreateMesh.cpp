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
//*** Implement mesh generation in AsFem                    ***
//*************************************************************

#include "Mesh/Mesh.h"

void Mesh::CreateMesh()
{
    if(IsRequiredInfoComplete())
    {
        if(nDims==1)
        {
            Create1DMesh();
        }
        else if(nDims==2)
        {
            Create2DMesh();
        }
        else
        {
            Create3DMesh();
        }
    }
    else
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: mesh info is not complete!!!    ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        mesh generation failed!!!       ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }
}
//*****************************************************
void Mesh::Create1DMesh()
{
    IsMeshCreated=false;

    int i,j,e,P;

    VTKCellType=4;
    if(MeshType=="edge2")
    {
        P=1;
        VTKCellType=3;
    }
    if(MeshType=="edge3") P=2;
    if(MeshType=="edge4") P=3;

    nElmts=Nx;
    nNodesPerElmt=P+1;
    nNodes=nElmts*P+1;

    double dx=(Xmax-Xmin)/(nNodes-1.0);


    NodeCoords.resize(nNodes*4,0.0);
    Conn.resize(nElmts*nNodesPerElmt,0);





    //TODO: it seems I need to consider about the mesh split by
    //      different cores, and the ghost effects!!!
    for(i=0;i<nNodes;++i)
    {
        NodeCoords[i*4  ]=1.0;
        NodeCoords[i*4+1]=Xmin+i*dx;
        NodeCoords[i*4+2]=0.0;
        NodeCoords[i*4+3]=0.0;
    }




    //TODO: make mesh generation also parallel!
    /*
    rankn=nElmts/size;
    iStart=rank*rankn;
    iEnd=(rank+1)*rankn;
    if(rank==size-1) iEnd=nElmts;
    */


    for(e=0;e<nElmts;++e)
    {
        for(j=1;j<=nNodesPerElmt;++j)
        {
            Conn[e*nNodesPerElmt+j-1]=e*P+j;
        }
    }

    IsMeshCreated=true;
    // Now we start to creat boundary element sets
    pair<string,int> LeftSet=pair<string,int>("left",1);
    pair<string,int> RightSet=pair<string,int>("right",nNodes);
    //BoundaryElmtSet.clear();
    //BoundaryElmtSet.push_back(LeftSet);
    //BoundaryElmtSet.push_back(RightSet);
}

//*****************************************************
// For 2D mesh generation
//*****************************************************
void Mesh::Create2DMesh()
{

}

//*****************************************************
// For 3D mesh generation
//*****************************************************
void Mesh::Create3DMesh()
{

}

