//*****************************************************
//*** AsFem: A simple finite element method program ***
//*** Author: Yang Bai                              ***
//*** Bug report: walkandthinker@gmail.com          ***
//*** QQ group: 797998860                           ***
//*****************************************************
//*** This class handle the dofs of asfem           ***
//***  i.e. local to global dof map and so on       ***
//*****************************************************

#include "DofHandler/DofHandler.h"

DofHandler::DofHandler()
{
    nDofs=0;nNodes=0;nElmts=0;
    HasDofMap=false;
    GlobalDofMap.clear();
}

void DofHandler::Release()
{
    if(HasDofMap)
    {
        GlobalDofMap.clear();
    }
}

//******************************
bool DofHandler::CreateLocalToGlobalDofMap(Mesh &mesh)
{
    HasDofMap=false;
    GlobalDofMap.clear();

    nNodes=mesh.GetNodesNum();
    nElmts=mesh.GetElmtsNum();
    nNodesPerElmts=mesh.GetNodesNumPerElmt();
    nDofsPerNode=mesh.GetDofsNumPerNode();
    nDofs=nNodes*nDofsPerNode;
    nDofsPerElmt=nNodesPerElmts*nDofsPerNode;

    GlobalDofMap=vector<int>(nElmts*nDofsPerNode*nNodesPerElmts,0);

    int e,i,j,k,iInd;
    for(e=1;e<=nElmts;e++)
    {
        for(i=1;i<=nNodesPerElmts;i++)
        {
            k=mesh.IthConnJthIndex(e,i);
            for(j=1;j<=nDofsPerNode;j++)
            {
                iInd=(k-1)*nDofsPerNode+j;
                GlobalDofMap[(e-1)*nDofsPerElmt+(i-1)*nDofsPerNode+j-1]=iInd;
            }
        }
    }
    HasDofMap=true;

    return HasDofMap;
}

//*******************************************
void DofHandler::GetLocalDofMap(const int e,int &ndofs,int *rInd,int *cInd)
{
    
    if(!HasDofMap)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't get dof map info, DofHandler hasn't been generated!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        CreateMesh should be called before this!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }

    if(e<1||e>nElmts)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: e=%8d is out of nElmts=%8d\n",e,nElmts);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }
    

    ndofs=nDofsPerElmt;
    for(int i=1;i<=nDofsPerElmt;i++)
    {
        rInd[i-1]=GlobalDofMap[(e-1)*nDofsPerElmt+i-1];
        cInd[i-1]=GlobalDofMap[(e-1)*nDofsPerElmt+i-1];
    }

}

//********************************
void DofHandler::PrintDofMap() const
{
    if(!HasDofMap)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't print dof map info, DofHandler hasn't been generated!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should create mesh first, then generate the dof map\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        then AsFem can print dof map for you!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"** Dof map information:                 ******\n");
    for(int e=1;e<=nElmts;e++)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** e=%6d: ",e);
        for(int i=1;i<=nDofsPerElmt;i++)
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"%5d ",GlobalDofMap[(e-1)*nDofsPerElmt+i-1]);
        }
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\n");
    }
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
    
}
