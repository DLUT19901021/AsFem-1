//*****************************************************
//*** AsFem: A simple finite element method program ***
//*** Author: Yang Bai @CopyRight                   ***
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
    HasDofMap=false;HasBCDofMap=false;
    GlobalDofMap.clear();
    GlobalBCDofMap.clear();
}

void DofHandler::Release()
{
    if(HasDofMap)
    {
        GlobalDofMap.clear();
        GlobalBCDofMap.clear();
        LeftSideDofMap.second.clear();
        RightSideDofMap.second.clear();
        BottomSideDofMap.second.clear();
        TopSideDofMap.second.clear();
    }
}

//******************************
bool DofHandler::CreateLocalToGlobalDofMap(Mesh &mesh,int ndofspernode)
{
    HasDofMap=false;
    GlobalDofMap.clear();

    nNodes=mesh.GetNodesNum();
    nElmts=mesh.GetElmtsNum();
    nNodesPerElmts=mesh.GetNodesNumPerElmt();
    nDofsPerNode=ndofspernode;
    nDofs=nNodes*nDofsPerNode;
    nDofsPerElmt=nNodesPerElmts*nDofsPerNode;

    GlobalDofMap.resize(nElmts*nDofsPerElmt,0);

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

    // Now to create boundary's dof map

    nNodesPerBCElmt=mesh.GetNodesNumPerBCElmt();
    nDofsPerBCElmt=nNodesPerBCElmt*nDofsPerNode;
    vector<int> SideDofMap;
    if(mesh.GetDims()==2)
    {
        // for left side
        pair<string,vector<int>> LeftNodeSet=mesh.GetSideSet("left");
        SideDofMap.resize(nDofsPerNode*LeftNodeSet.second.size(),0);
        for(i=1;i<=LeftNodeSet.second.size();i++)
        {
            for(j=1;j<=nDofsPerNode;j++)
            {
                SideDofMap[(i-1)*nDofsPerNode+j-1]=(LeftNodeSet.second[i-1]-1)*nDofsPerNode+j;
            }
        }
        LeftSideDofMap=make_pair("left",SideDofMap);
        SideDofMap.clear();

        // for right side
        pair<string,vector<int>> RightNodeSet=mesh.GetSideSet("right");
        SideDofMap.resize(nDofsPerNode*RightNodeSet.second.size(),0);
        for(i=1;i<=RightNodeSet.second.size();i++)
        {
            for(j=1;j<=nDofsPerNode;j++)
            {
                SideDofMap[(i-1)*nDofsPerNode+j-1]=(RightNodeSet.second[i-1]-1)*nDofsPerNode+j;
            }
        }
        RightSideDofMap=make_pair("right",SideDofMap);
        SideDofMap.clear();

        // for bottom side
        pair<string,vector<int>> BottomNodeSet=mesh.GetSideSet("bottom");
        SideDofMap.resize(nDofsPerNode*BottomNodeSet.second.size(),0);
        for(i=1;i<=BottomNodeSet.second.size();i++)
        {
            for(j=1;j<=nDofsPerNode;j++)
            {
                SideDofMap[(i-1)*nDofsPerNode+j-1]=(BottomNodeSet.second[i-1]-1)*nDofsPerNode+j;
            }
        }
        BottomSideDofMap=make_pair("bottom",SideDofMap);
        SideDofMap.clear();

        // for top side
        pair<string,vector<int>> TopNodeSet=mesh.GetSideSet("top");
        SideDofMap.resize(nDofsPerNode*TopNodeSet.second.size(),0);
        for(i=1;i<=TopNodeSet.second.size();i++)
        {
            for(j=1;j<=nDofsPerNode;j++)
            {
                SideDofMap[(i-1)*nDofsPerNode+j-1]=(TopNodeSet.second[i-1]-1)*nDofsPerNode+j;
            }
        }
        TopSideDofMap=make_pair("top",SideDofMap);
        SideDofMap.clear();

    }

    GlobalBCDofMap.push_back(LeftSideDofMap);
    GlobalBCDofMap.push_back(RightSideDofMap);
    GlobalBCDofMap.push_back(BottomSideDofMap);
    GlobalBCDofMap.push_back(TopSideDofMap);

    HasBCDofMap=true;

    return HasDofMap;
}

//*******************************************
void DofHandler::GetLocalDofMap(const int e,int &ndofs,int (&rInd)[500],int (&cInd)[500]) const
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

void DofHandler::GetLocalBCDofMap(string sidename,const int e,int &ndofs,int (&Ind)[500]) const
{
    if(!HasBCDofMap)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't get bc's dof map info, DofHandler hasn't been generated!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        CreateMesh should be called before this!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }

    ndofs=nDofsPerBCElmt;// the length of Ind array

    if(sidename=="left")
    {
        if(e<1||e>LeftSideDofMap.second.size()/nDofsPerBCElmt)
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: e=%d is out of boundary element's range!\n",e);
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
            PetscFinalize();
            abort();
        }
        for(int i=1;i<=nDofsPerBCElmt;i++)
        {
            Ind[i-1]=LeftSideDofMap.second[(e-1)*nDofsPerBCElmt+i-1];
        }
    }
    else if(sidename=="right")
    {
        if(e<1||e>RightSideDofMap.second.size()/nDofsPerBCElmt)
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: e=%d is out of boundary element's range!\n",e);
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
            PetscFinalize();
            abort();
        }
        for(int i=1;i<=nDofsPerBCElmt;i++)
        {
            Ind[i-1]=RightSideDofMap.second[(e-1)*nDofsPerBCElmt+i-1];
        }
    }
    else if(sidename=="bottom")
    {
        if(e<1||e>BottomSideDofMap.second.size()/nDofsPerBCElmt)
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: e=%d is out of boundary element's range!\n",e);
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
            PetscFinalize();
            abort();
        }
        for(int i=1;i<=nDofsPerBCElmt;i++)
        {
            Ind[i-1]=BottomSideDofMap.second[(e-1)*nDofsPerBCElmt+i-1];
        }
    }
    else if(sidename=="top")
    {
        if(e<1||e>TopSideDofMap.second.size()/nDofsPerBCElmt)
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: e=%d is out of boundary element's range!\n",e);
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
            PetscFinalize();
            abort();
        }
        for(int i=1;i<=nDofsPerBCElmt;i++)
        {
            Ind[i-1]=TopSideDofMap.second[(e-1)*nDofsPerBCElmt+i-1];
        }
    }
    else
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: unsupported sidename(=%s)!\n",sidename.c_str());
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        please use left,right,top,bottom,back,front!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
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
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Dof map information:                   ***\n");
    for(int e=1;e<=nElmts;e++)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** e=%6d: ",e);
        for(int i=1;i<=nDofsPerElmt;i++)
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"%5d ",GlobalDofMap[(e-1)*nDofsPerElmt+i-1]);
        }
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\n");
    }

    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Boundary dof map information:          ***\n");

    int dofind[500],len;
    for(int i=0;i<GlobalBCDofMap.size();i++)
    {
        string sidename=GlobalBCDofMap[i].first;
        for(int e=1;e<=GlobalBCDofMap[i].second.size()/nDofsPerBCElmt;e++)
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** %8s side->",sidename.c_str());
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"e=%3d:",e);

            GetLocalBCDofMap(sidename,e,len,dofind);
            for(int j=0;j<len;j++)
            {
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"%5d ",dofind[j]);
            }
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\n");
        }
    }
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
    
}
