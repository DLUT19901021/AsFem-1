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
    nDims=dim;nNodes=0;nElmts=0;
    nMaxNodesPerElmt=nnodesperelmt;
    nDofsPerNode=ndofspernode;

    Xmin=Ymin=Zmin=1.0e12;
    Xmax=Ymax=Zmax=1.0e-12;
    Nx=0;Ny=0;Nz=0;

    Mesh1DList.clear();
    Mesh2DList.clear();

    Conn.clear();
    NodeCoords.clear();

}

void Mesh::Release()
{
    if(IsInit)
    {
        Mesh1DList.clear();
        Mesh2DList.clear();
        Conn.clear();
        NodeCoords.clear();
    }
}

//*************************************
void Mesh::Add1DMesh(Mesh1D &mesh1d)
{
    if(Has2DMesh)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: already has 2D mesh, can't insert 1d mesh!!!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        currently only one kind of mesh is supported!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }

    if(nDims!=1)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: you can't insert 2/3D mesh into 1d case!!!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        currently only one kind of mesh is supported!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }

    if(!mesh1d.IsMeshGenerated())
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: you can't insert this 1d mesh to the mesh list!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        since no mesh has been generated for it     !!!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }

    Mesh1DList.push_back(mesh1d);
    if(mesh1d.GetXmin()<Xmin) Xmin=mesh1d.GetXmin();
    if(mesh1d.GetXmax()>Xmax) Xmax=mesh1d.GetXmax();
    Has1DMesh=true;
    nNodes+=mesh1d.GetNodesNum();
    nElmts+=mesh1d.GetElmtsNum();
}
void Mesh::Add2DMesh(Mesh2D &mesh2d)
{
    if(Has1DMesh)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: already has 1D mesh, can't insert 2d mesh!!!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        currently only one kind of mesh is supported!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }

    if(nDims!=2)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: you can't insert 1/3D mesh into 2d case!!!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        currently only one kind of mesh is supported!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }

    if(!mesh2d.IsMeshGenerated())
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: you can't insert this 2d mesh to the mesh list!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        since no mesh has been generated for it     !!!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }

    Mesh2DList.push_back(mesh2d);
    if(mesh2d.GetXmin()<Xmin) Xmin=mesh2d.GetXmin();
    if(mesh2d.GetXmax()>Xmax) Xmax=mesh2d.GetXmax();
    if(mesh2d.GetYmin()<Ymin) Ymin=mesh2d.GetYmin();
    if(mesh2d.GetYmax()>Ymax) Ymax=mesh2d.GetYmax();
    Has2DMesh=true;
    nNodes+=mesh2d.GetNodesNum();
    nElmts+=mesh2d.GetElmtsNum();
}

//********************************************
void Mesh::Init()
{
    // TODO: now only one mesh block is 
    //       supported, in the future
    //       multiple blocks should be supported!
    IsInit=false;

    NodeCoords=vector<double>(4*nNodes,0.0);
    NodeCoords.resize(4*nNodes);
    Conn=vector<int>(nElmts*nMaxNodesPerElmt,0);
    Conn.resize(nElmts*nMaxNodesPerElmt);

    if(Has1DMesh)
    {
        for(int i=1;i<=Mesh1DList[0].GetNodesNum();i++)
        {
            NodeCoords[(i-1)*4+0]=Mesh1DList[0].IthNodeJthCoords(i,0);
            NodeCoords[(i-1)*4+1]=Mesh1DList[0].IthNodeJthCoords(i,1);
            NodeCoords[(i-1)*4+2]=Mesh1DList[0].IthNodeJthCoords(i,2);
            NodeCoords[(i-1)*4+3]=Mesh1DList[0].IthNodeJthCoords(i,3);
        }
        for(int e=1;e<=Mesh1DList[0].GetElmtsNum();e++)
        {
            for(int j=1;j<=Mesh1DList[0].GetNodesNumPerElmt();j++)
            {
                Conn[(e-1)*nMaxNodesPerElmt+j-1]=Mesh1DList[0].IthConnJthIndex(e,j);
            }
        }
        IsInit=true;
    }
    else if(Has2DMesh)
    {
        for(int i=1;i<=Mesh2DList[0].GetNodesNum();i++)
        {
            NodeCoords[(i-1)*4+0]=Mesh2DList[0].IthNodeJthCoords(i,0);
            NodeCoords[(i-1)*4+1]=Mesh2DList[0].IthNodeJthCoords(i,1);
            NodeCoords[(i-1)*4+2]=Mesh2DList[0].IthNodeJthCoords(i,2);
            NodeCoords[(i-1)*4+3]=Mesh2DList[0].IthNodeJthCoords(i,3);
        }
        for(int e=1;e<=Mesh2DList[0].GetElmtsNum();e++)
        {
            for(int j=1;j<=Mesh2DList[0].GetNodesNumPerElmt();j++)
            {
                Conn[(e-1)*nMaxNodesPerElmt+j-1]=Mesh2DList[0].IthConnJthIndex(e,j);
            }
        }
        VTKCellType=Mesh2DList[0].GetVTKCellType();
        IsInit=true;
    }
}

//***********************************************
double Mesh::IthNodeJthCoords(int i,int j) const
{
    if(!IsInit)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't get node's coordinates, mesh hasn't been generated!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        CreateMesh should be called before this!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }

    if(i<1||i>nNodes)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: i=%6d is out of nNodes=%6d\n",i,nNodes);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }
    if(j<0||j>3)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: j=%2d is out 0~3\n",i,nNodes);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        j should be 0->w,1->x,2->y,3->z!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }

    return NodeCoords[4*(i-1)+j];
}

int Mesh::IthConnJthIndex(int e,int j) const
{
    if(!IsInit)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't get connectivity info, mesh hasn't been generated!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        CreateMesh should be called before this!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }

    if(e<1||e>nElmts)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: e=%8d is out of nElmts=%8d\n",e,nNodes);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }
    if(j<1||j>nMaxNodesPerElmt)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: j=%2d is out of nNodesPerElmt=%2d\n",j,nMaxNodesPerElmt);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        j should be 1~%2d!\n",nMaxNodesPerElmt);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }

    return Conn[(e-1)*nMaxNodesPerElmt+j-1];

}

//*****************************************
void Mesh::PrintMeshInfo(string str) const
{
    if(!IsInit)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't print mesh info, mesh hasn't been generated!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        CreateMesh should be called before this!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }


    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
    if(str.size()>1)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** str= %s\n",str.c_str());
    }
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Mesh information:                      ***\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   nDims=%2d                             ***\n",nDims);
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   nNodes=%10d                    ***\n",nNodes);
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   nElmts=%10D                    ***\n",nElmts);
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   nNodesPerElmt=%3d                    ***\n",nMaxNodesPerElmt);
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");

}
void Mesh::PrintMeshDetailInfo(string str) const
{
    int i,e;
    double x,y,z,w;
    if(!IsInit)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't print mesh info, mesh hasn't been generated!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        CreateMesh should be called before this!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }


    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
    if(str.size()>1)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** str= %s\n",str.c_str());
    }
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Mesh detailed information:              ***\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   nDims=%2d                             ***\n",nDims);
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   nNodes=%10d                    ***\n",nNodes);
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   nElmts=%10D                    ***\n",nElmts);
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   nNodesPerElmt=%3d                    ***\n",nMaxNodesPerElmt);
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Nodes' coordinates:                    ***\n");
    for(i=1;i<=nNodes;++i)
    {
        w=IthNodeJthCoords(i,0);
        x=IthNodeJthCoords(i,1);
        y=IthNodeJthCoords(i,2);
        z=IthNodeJthCoords(i,3);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** %6d-th node:",i);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"x=%12.5e,y=%12.5e,z=%12.5e,w=%6.2f\n",x,y,z,w);
    }

    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Element's connectivity:                ***\n");
    for(e=1;e<=nElmts;++e)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** e=%6d:",e);
        for(i=1;i<=nMaxNodesPerElmt;++i)
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"%5d ",IthConnJthIndex(e,i));
        }
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\n");
    }
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");

}

