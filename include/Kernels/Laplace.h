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

#ifndef ASFEM_LAPLACE_H
#define ASFEM_LAPLACE_H

#include "Kernels/Kernel.h"

class Laplace:public Kernel
{
public:
    Laplace(int dim,string kernelname,string dofname,int dofindex);

protected:
    virtual double computeQpResidual() override ;
    virtual double computeQpJacobian() override ;

private:
    double _F,_K;

    double residual,jacobian;

};


#endif //ASFEM_LAPLACE_H
