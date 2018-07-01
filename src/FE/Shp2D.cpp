//*****************************************************
//*** AsFem: A simple finite element method program ***
//*** Author: Yang Bai @CopyRight                   ***
//*** Bug report: walkandthinker@gmail.com          ***
//*** QQ group: 797998860                           ***
//*****************************************************
//*** Define 2D shape functions for N and dNdx      ***
//*****************************************************

#include "FE/FE.h"

void FE::Shp2D(const int &ndim,
               const int &nnodes,
               const double &xi,
               const double &eta,
               const double (&Coords)[27][4],
               double (&shp)[27][4],
               double &DetJac,
               bool IsNature)
{
    int i;
    double dxdxi, dydxi;
    double dxdeta, dydeta;
    double temp;
    double XJac[2][2],Jac[2][2];


    if (ndim!=2)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: Wrong dimension case !!!        ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: Wrong dimension(dim=%d)!!!",ndim);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error happens in Shp2D!!!              ***\n");
        abort();
    }

    if (nnodes<3 || nnodes>9)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error:Wrong nodes on current element(nNodes=%d)!!!\n",nnodes);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem only support 2D-4,8,9 nodes mesh!***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error happens in Shp2D!!!              ***\n");
        abort();
    }



    if(nnodes==4)
    {
        // 2D-4Nodes rectangle element
        shp[0][0]=(1.0-xi)*(1.0-eta)/4.0;
        shp[1][0]=(1.0+xi)*(1.0-eta)/4.0;
        shp[2][0]=(1.0+xi)*(1.0+eta)/4.0;
        shp[3][0]=(1.0-xi)*(1.0+eta)/4.0;

        shp[0][1]= (eta-1.0)/4.0;
        shp[0][2]= (xi -1.0)/4.0;

        shp[1][1]= (1.0-eta)/4.0;
        shp[1][2]=-(1.0+xi )/4.0;

        shp[2][1]= (1.0+eta)/4.0;
        shp[2][2]= (1.0+xi )/4.0;

        shp[3][1]=-(1.0+eta)/4.0;
        shp[3][2]= (1.0-xi )/4.0;
    }
    else if(nnodes==8)
    {
        // 2D-8Nodes rectangle element
        shp[0][0]=(1.0-xi)*(1.0-eta)*(-xi-eta-1.0)/4.0;
        shp[1][0]=(1.0+xi)*(1.0-eta)*( xi-eta-1.0)/4.0;
        shp[2][0]=(1.0+xi)*(1.0+eta)*( xi+eta-1.0)/4.0;
        shp[3][0]=(1.0-xi)*(1.0+eta)*(-xi+eta-1.0)/4.0;
        shp[4][0]=(1.0-xi*xi)*(1.0-eta    )/2.0;
        shp[5][0]=(1.0+xi   )*(1.0-eta*eta)/2.0;
        shp[6][0]=(1.0-xi*xi)*(1.0+eta    )/2.0;
        shp[7][0]=(1.0-xi   )*(1.0-eta*eta)/2.0;

        // derivatives over xi and eta
        shp[0][1]=(1.0-eta)*(2.0*xi+eta)/4.0;
        shp[0][2]=(1.0-xi )*(xi+2.0*eta)/4.0;

        shp[1][1]=(1.0-eta)*(2.0*xi-eta)/4.0;
        shp[1][2]=(1.0+xi )*(2.0*eta-xi)/4.0;

        shp[2][1]=(1.0+eta)*(2.0*xi+eta)/4.0;
        shp[2][2]=(1.0+xi )*(xi+2.0*eta)/4.0;

        shp[3][1]=(1.0+eta)*(2.0*xi-eta)/4.0;
        shp[3][2]=(1.0-xi )*(2.0*eta-xi)/4.0;

        shp[4][1]=xi*(eta-1.0);
        shp[4][2]=(xi*xi-1.0)/2.0;

        shp[5][1]=(1.0-eta*eta)/2.0;
        shp[5][2]=-(1.0+xi)*eta;

        shp[6][1]=-xi*(1.0+eta);
        shp[6][2]=(1.0-xi*xi)/2.0;

        shp[7][1]=(eta*eta-1.0)/2.0;
        shp[7][2]=(xi-1.0)*eta;
    }
    else if(nnodes==9)
    {

        // 2D-9Nodes rectangle element
        shp[0][0]=(xi*xi-xi )*(eta*eta-eta)/4.0;
        shp[1][0]=(xi*xi+xi )*(eta*eta-eta)/4.0;
        shp[2][0]=(xi*xi+xi )*(eta*eta+eta)/4.0;
        shp[3][0]=(xi*xi-xi )*(eta*eta+eta)/4.0;
        shp[4][0]=(1.0-xi*xi)*(eta*eta-eta)/2.0;
        shp[5][0]=(xi*xi+xi )*(1.0-eta*eta)/2.0;
        shp[6][0]=(1.0-xi*xi)*(eta*eta+eta)/2.0;
        shp[7][0]=(xi*xi-xi )*(1.0-eta*eta)/2.0;
        shp[8][0]=(1.0-xi*xi)*(1.0-eta*eta);

        shp[0][1]=(2.0*xi-1.0)*(eta*eta-eta)/4.0;
        shp[0][2]=(xi*xi-xi  )*(2.0*eta-1.0)/4.0;

        shp[1][1]=(2.0*xi+1.0)*(eta*eta-eta)/4.0;
        shp[1][2]=(xi*xi+xi  )*(2.0*eta-1.0)/4.0;

        shp[2][1]=(2.0*xi+1.0)*(eta*eta+eta)/4.0;
        shp[2][2]=(xi*xi+xi  )*(2.0*eta+1.0)/4.0;

        shp[3][1]=(2.0*xi-1.0)*(eta*eta+eta)/4.0;
        shp[3][2]=(xi*xi-xi  )*(2.0*eta+1.0)/4.0;

        shp[4][1]=-xi*(eta*eta-eta);
        shp[4][2]=(1.0-xi*xi )*(2.0*eta-1.0)/2.0;

        shp[5][1]=(2.0*xi+1.0)*(1.0-eta*eta)/2.0;
        shp[5][2]=-(xi*xi+xi )*eta;

        shp[6][1]=-xi*(eta*eta+eta);
        shp[6][2]=(1.0-xi*xi )*(2.0*eta+1.0)/2.0;

        shp[7][1]=(2.0*xi-1.0)*(1.0-eta*eta)/2.0;
        shp[7][2]=-(xi*xi-xi )*eta;

        shp[8][1]=-2.0*xi*(1.0-eta*eta);
        shp[8][2]=-2.0*eta*(1.0-xi*xi);
    }


    // compute jacob transform matrix
    dxdxi=0.0;dxdeta=0.0;
    dydxi=0.0;dydeta=0.0;
    for(i=0;i<nnodes;i++)
    {
        dxdxi +=Coords[i][1]*shp[i][1];
        dxdeta+=Coords[i][1]*shp[i][2];

        dydxi +=Coords[i][2]*shp[i][1];
        dydeta+=Coords[i][2]*shp[i][2];
    }


    Jac[0][0]= dxdxi;Jac[0][1]= dydxi;
    Jac[1][0]=dxdeta;Jac[1][1]=dydeta;


    DetJac=Jac[0][0]*Jac[1][1]-Jac[0][1]*Jac[1][0];

    if(fabs(DetJac)<1.e-13)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: Singular in element(Shp2D)!!!   ***\n");
        abort();
    }
    XJac[0][0]= Jac[1][1]/DetJac;
    XJac[0][1]=-Jac[0][1]/DetJac;
    XJac[1][0]=-Jac[1][0]/DetJac;
    XJac[1][1]= Jac[0][0]/DetJac;



    if(!IsNature)
    {
        for(i=0;i<nnodes;++i)
        {
            temp=XJac[0][0]*shp[i][1]+XJac[0][1]*shp[i][2];
            shp[i][2]=XJac[1][0]*shp[i][1]+XJac[1][1]*shp[i][2];
            shp[i][1]=temp;
        }
    }
}

