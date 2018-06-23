//*****************************************************
//*** AsFem: A simple finite element method program ***
//*** Author: Yang Bai                              ***
//*** Bug report: walkandthinker@gmail.com          ***
//*** QQ group: 797998860                           ***
//*****************************************************
//*** Define class for gauss integration            ***
//***   generate gauss point and weight for inte    ***
//***                                               ***
//*****************************************************

#include "FE/Qrule.h"

Qrule::Qrule(int dim,int n)
{
    if(dim<1||dim>3)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: dim=%2d is invalid for gauss integration!!!\n",dim);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }

    if(n<2||n>4)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: number of gauss point=%d is invalid for a FEM mesh!!!\n",n);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }

    ngp=n;
    nGaussPoints=ngp;
    if(dim==2) nGaussPoints=ngp*ngp;
    if(dim==3) nGaussPoints=ngp*ngp*ngp;

    nDims=dim;

    GenerateGaussPoints();
}

void Qrule::GenerateGaussPoints()
{
    int i,j,k,l;
    weight.resize(nGaussPoints);
    coords.resize(nGaussPoints*nDims);


    if(ngp==2)
    {
        if(nDims==1)
        {
            weight[0]=p1w[0];weight[1]=p1w[1];
            coords[0]=p1xi[0];coords[1]=p1xi[1];
        }
        else if(nDims==2)
        {
            k=0;
            for(j=0;j<ngp;j++)
            {
                for(i=0;i<ngp;i++)
                {
                    weight[k]=p1w[i]*p1w[j];
                    coords[k*nDims+0]=p1xi[i];
                    coords[k*nDims+1]=p1xi[j];
                    k+=1;
                }
            }
        }
        else if(nDims==3)
        {
            l=0;
            for(k=0;k<ngp;k++)
            {
                for(j=0;j<ngp;j++)
                {
                    for(i=0;i<ngp;i++)
                    {
                        weight[l]=p1w[i]*p1w[j]*p1w[k];
                        coords[l*nDims+0]=p1xi[i];
                        coords[l*nDims+1]=p1xi[j];
                        coords[l*nDims+2]=p1xi[k];
                        l+=1;
                    }
                }
            }
        }
    }
    else if(ngp==3)
    {
        if(nDims==1)
        {
            weight[0]=p2w[0];weight[1]=p2w[1];weight[2]=p2w[2];
            coords[0]=p2xi[0];coords[1]=p2xi[1];coords[2]=p2xi[2];
        }
        else if(nDims==2)
        {
            k=0;
            for(j=0;j<ngp;j++)
            {
                for(i=0;i<ngp;i++)
                {
                    weight[k]=p2w[i]*p2w[j];
                    coords[k*nDims+0]=p2xi[i];
                    coords[k*nDims+1]=p2xi[j];
                    k+=1;
                }
            }
        }
        else if(nDims==3)
        {
            l=0;
            for(k=0;k<ngp;k++)
            {
                for(j=0;j<ngp;j++)
                {
                    for(i=0;i<ngp;i++)
                    {
                        weight[l]=p2w[i]*p2w[j]*p2w[k];
                        coords[l*nDims+0]=p2xi[i];
                        coords[l*nDims+1]=p2xi[j];
                        coords[l*nDims+2]=p2xi[k];
                        l+=1;
                    }
                }
            }
        }
    }
    else if(ngp==4)
    {
        if(nDims==1)
        {
            weight[0]=p3w[0];weight[1]=p3w[1];weight[2]=p3w[2];weight[3]=p3w[3];
            coords[0]=p3xi[0];coords[1]=p3xi[1];coords[2]=p3xi[2];coords[3]=p3xi[3];
        }
        else if(nDims==2)
        {
            k=0;
            for(j=0;j<ngp;j++)
            {
                for(i=0;i<ngp;i++)
                {
                    weight[k]=p3w[i]*p3w[j];
                    coords[k*nDims+0]=p3xi[i];
                    coords[k*nDims+1]=p3xi[j];
                    k+=1;
                }
            }
        }
        else if(nDims==3)
        {
            l=0;
            for(k=0;k<ngp;k++)
            {
                for(j=0;j<ngp;j++)
                {
                    for(i=0;i<ngp;i++)
                    {
                        weight[l]=p3w[i]*p3w[j]*p3w[k];
                        coords[l*nDims+0]=p3xi[i];
                        coords[l*nDims+1]=p3xi[j];
                        coords[l*nDims+2]=p3xi[k];
                        l+=1;
                    }
                }
            }
        }
    }
}

//****************************************
double Qrule::GetComponent(int i,int j) const
{
    if(i<1||i>nGaussPoints)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: i=%d is invalid for a gauss point!!!\n",i);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }
    if(j<0||j>nDims)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: j=%d is an invalid component of gauss point!!!\n",j);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }

    if(j==0)
    {
        return weight[i-1];
    }
    else
    {
        return coords[(i-1)*nDims+j-1];
    }
}

//***************************************
void Qrule::PrintGaussPointInfo() const
{
    int i,j;
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Gauss point information:               ***\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   ngp=%2d, nGauss=%2d                    ***\n",ngp,nGaussPoints);
    for(i=1;i<=nGaussPoints;++i)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** i=%2d: ",i);
        for(j=1;j<=nDims;j++)
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"%10.3e ",GetComponent(i,j));
        }
        PetscSynchronizedPrintf(PETSC_COMM_WORLD," %8.2e  ***\n",GetComponent(i,0));
    }
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
}