//*****************************************************
//*** AsFem: A simple finite element method program ***
//*** Author: Yang Bai @CopyRight                   ***
//*** Bug report: walkandthinker@gmail.com          ***
//*** QQ group: 797998860                           ***
//*****************************************************
//*** This class defines the core parts of AsFem    ***
//***  basical and elemental calculation of FEM     ***
//***  element loop and matrix assembly code        ***
//***  should be put here                           ***
//*****************************************************

#include "FE/FE.h"

FE::FE(Mesh &mesh,DofHandler &dofHandler)
{
    DIM=mesh.GetDims();
    int p=1;
    if(DIM==1)
    {
        p=mesh.GetNodesNumPerElmt()-1;
    }
    else if(DIM==2)
    {
        if(mesh.GetNodesNumPerElmt()>=8)
        {
            p=2;
        }
    }
    else if(DIM==3)
    {
        if(mesh.GetNodesNumPerElmt()>=20)
        {
            p=2;
        }
    }

    qrule=new Qrule(DIM,p);
}

void FE::FormKR(Mesh &mesh,DofHandler &dofHandler,EquationSystem &equationSystem)
{
    int e,i,j;
    int ii,jj;
    double k[270][270],rhs[270];
    double elCoords[27][4]={0.0};
    double elU[27][2]={0.0};//0->u,1->v
    PetscInt elConn[27]={0};


    VecScatterCreateToAll(equationSystem.U,&Uscatter,&Useq);
    VecScatterBegin(Uscatter,equationSystem.U,Useq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(Uscatter,equationSystem.U,Useq,INSERT_VALUES,SCATTER_FORWARD);

    VecScatterCreateToAll(equationSystem.V,&Vscatter,&Vseq);
    VecScatterBegin(Uscatter,equationSystem.V,Useq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(Uscatter,equationSystem.V,Vseq,INSERT_VALUES,SCATTER_FORWARD);

    // split the element into rank-rated index, then it can be carried parallel
    PetscMPIInt rank,size;
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    MPI_Comm_size(PETSC_COMM_WORLD,&size);

    int rankne=mesh.GetElmtsNum()/size;
    int eStart=rank*rankne;
    int eEnd=(rank+1)*rankne;
    if(rank==size-1) eEnd=mesh.GetElmtsNum();

    for(e=eStart;e<eEnd;e++)
    {
        nNodesPerElmt=mesh.GetNodesNumPerElmt();

        for(i=1;i<=nNodesPerElmt;++i)
        {
            j=mesh.IthConnJthIndex(e+1,i);
            elConn[i-1]=j-1;
            elCoords[i-1][0]=0.0;
            elCoords[i-1][1]=mesh.IthNodeJthCoords(j,1);
            elCoords[i-1][2]=mesh.IthNodeJthCoords(j,2);
            elCoords[i-1][3]=mesh.IthNodeJthCoords(j,3);

            for(ii=1;ii<=mesh.GetDofsNumPerNode();ii++)
            {
                jj=(j-1)*mesh.GetDofsNumPerNode()+ii-1;
                VecGetValues(Useq,1,&jj,&elU[(i-1)*mesh.GetDofsNumPerNode()+ii-1][0]);
                VecGetValues(Vseq,1,&jj,&elU[(i-1)*mesh.GetDofsNumPerNode()+ii-1][1]);
            }
        }

    }
}

void FE::QpLoop(const double (&elCoords)[27][4], const double (&elU)[27][2],
                double (&k)[270][270],double (&rhs)[270])
{
    double xi,eta,zeta;
    double shp[27][4],xsj;
    int I,J;
    for(qp=1;qp<=qrule->GetQPointsNum();qp++)
    {
        if(DIM==1)
        {
            xi=qrule->GetComponent(qp,1);
            Shp1D(DIM,nNodesPerElmt,xi,elCoords,shp,xsj,true);
        }
        else if(DIM==2)
        {
            xi=qrule->GetComponent(qp,1);
            eta=qrule->GetComponent(qp,2);
            Shp2D(DIM,nNodesPerElmt,xi,eta,elCoords,shp,xsj,true);
        }
        else if(DIM==3)
        {
            xi=qrule->GetComponent(qp,1);
            eta=qrule->GetComponent(qp,2);
            zeta=qrule->GetComponent(qp,3);
            Shp3D(DIM,nNodesPerElmt,xi,eta,zeta,elCoords,shp,xsj,true);
        }
        JxW=xsj*qrule->GetComponent(qp,0);

        for(I=0;1<nNodesPerElmt;I++)
        {
            ComputeResidual();
            for(J=0;J<nNodesPerElmt;J++)
            {
                ComputeJacobian();
            }
        }
    }
}


