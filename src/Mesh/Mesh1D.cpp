//
// Created by Y. Bai on 17.06.18.
//

#include "Mesh/Mesh1D.h"

Mesh1D::Mesh1D(double xmin, double xmax, int nx,string elmttype)
{
    if(xmin>=xmax)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: xmin=%12.6e is larger than xmax=12.6e!!!\n",xmin,xmax);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }

    if(nx<2)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: nx=%6d is too less for a FEM problem!!!\n",nx);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }

    string str=RemoveSpace(elmttype);
    elmttype=StrToLower(str);
    if(elmttype.find("edge2")==string::npos&&
       elmttype.find("edge3")==string::npos&&
       elmttype.find("edge4")==string::npos)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: unsupported mesh type !!!\n",nx);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        only type=edge2, edge3 or edge4 is supported for 1D mesh !!!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }

    Xmin=xmin;Xmax=xmax;
    Nx=nx;
    nNodes=0;nNodesPerElmt=0;nElmts=0;
    MeshGenerated=false;
    ElmtType=elmttype;

    BoundaryElmtSet.clear();
}

void Mesh1D::Release()
{
    if(MeshGenerated)
    {
        NodeCoords.clear();
        Conn.clear();
        BoundaryElmtSet.clear();
    }
}
//**********************************************
bool Mesh1D::CreateMesh()
{
    int iStart,iEnd,rankn;


    int i,j,e;
    MeshGenerated=false;

    VTKCellType=4;
    if(ElmtType=="edge2")
    {
        P=1;
        VTKCellType=3;
    }
    if(ElmtType=="edge3") P=2;
    if(ElmtType=="edge4") P=3;

    nElmts=Nx;
    nNodesPerElmt=P+1;
    nNodes=nElmts*P+1;

    double dx=(Xmax-Xmin)/(nNodes-1.0);

    
    NodeCoords.resize(nNodes*4,0.0);
    Conn.resize(nElmts*nNodesPerElmt,0);

    

    

    //TODO: it seems I need to consider about the mesh split by
    //      different cores, and the ghost effects!!!
    iStart=0;iEnd=nNodes;
    for(i=iStart;i<iEnd;++i)
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

    iStart=0;iEnd=nElmts;
    for(e=iStart;e<iEnd;++e)
    {
        for(j=1;j<=nNodesPerElmt;++j)
        {
            Conn[e*nNodesPerElmt+j-1]=e*P+j;
        }
    }

    MeshGenerated=true;
    // Now we start to creat boundary element sets
    pair<string,int> LeftSet=pair<string,int>("left",1);
    pair<string,int> RightSet=pair<string,int>("right",nNodes);
    BoundaryElmtSet.clear();
    BoundaryElmtSet.push_back(LeftSet);
    BoundaryElmtSet.push_back(RightSet);

    return MeshGenerated;
}

// Get element information

double Mesh1D::IthNodeJthCoords(int i,int j) const
{
    if(!MeshGenerated)
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

int Mesh1D::IthConnJthIndex(int e,int j) const
{
    if(!MeshGenerated)
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
    if(j<1||j>nNodesPerElmt)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: j=%2d is out of nNodesPerElmt=%2d\n",j,nNodesPerElmt);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        j should be 1~%2d!\n",nNodesPerElmt);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }

    return Conn[(e-1)*nNodesPerElmt+j-1];
}

//****************************************
void Mesh1D::PrintMeshInfo(string str) const
{
    if(!MeshGenerated)
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
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** 1D mesh info:                          ***\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   element order=%2d                     ***\n",P);
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   nNodes=%10d                    ***\n",nNodes);
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   nElmts=%10D                    ***\n",nElmts);
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   nNodesPerElmt=%3d                    ***\n",nNodesPerElmt);
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");

}

void Mesh1D::PrintMeshDetailInfo(string str) const
{
    int i,e;
    double x,y,z,w;
    if(!MeshGenerated)
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
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** 1D mesh detail info:                   ***\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   element order=%2d                     ***\n",P);
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   nNodes=%10d                    ***\n",nNodes);
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   nElmts=%10D                    ***\n",nElmts);
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   nNodesPerElmt=%3d                    ***\n",nNodesPerElmt);
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
        for(i=1;i<=nNodesPerElmt;++i)
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"%5d ",IthConnJthIndex(e,i));
        }
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\n");
    }
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
}