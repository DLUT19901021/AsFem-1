//*****************************************************
//*** AsFem: A simple finite element method program ***
//*** Author: Yang Bai @CopyRight                   ***
//*** Bug report: walkandthinker@gmail.com          ***
//*** QQ group: 797998860                           ***
//*****************************************************
//***  the kernel for laplace operation             ***
//***  k*div*(grad[u])=F                            ***
//***                                               ***
//*****************************************************

#include "Kernels/Laplace.h"

Laplace::Laplace(int dim,string kernelname,string dofname,int dofindex)
{
    KernelName=kernelname;
    DofName=dofname;
    DofIndex=dofindex;
    DIM=dim;
}

double Laplace::computeQpResidual()
{
    residual=0.0;
    for(int i=1;i<=DIM;i++)
    {
        residual-=_K*_grad_u[i]*_grad_testI[i];
    }
    residual-=_F*_testI;
    return residual;
}

double Laplace::computeQpJacobian()
{
    jacobian=0.0;
    for(int i=1;i<=DIM;i++)
    {
        jacobian+=-_K*_grad_phiJ[i]*_grad_testI[i];
    }
    return jacobian;
}

