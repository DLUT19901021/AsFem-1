//*****************************************************
//*** AsFem: A simple finite element method program ***
//*** Author: Yang Bai @CopyRight                   ***
//*** Bug report: walkandthinker@gmail.com          ***
//*** QQ group: 797998860                           ***
//*****************************************************
//*** This class defines the core parts of AsFem    ***
//***  the kernels for FEM calculation              ***
//***  if you want write your own kernel            ***
//***  just copy the template and modifiy it        ***
//*****************************************************

#ifndef ASFEM_KERNELS_H
#define ASFEM_KERNELS_H

#include <iostream>
#include <vector>
#include <string>



using namespace std;

class Kernel
{
public:
    int GetDofIndex() const { return DofIndex;}
    string GetDofName() const { return DofName;}
    string GetKernelName() const { return KernelName;}

    virtual double computeQpResidual()=0;
    virtual double computeQpJacobian()=0;
    //virtual double computeQpOffDiagJacobian();
    //virtual void getQpMaterialProperty();

    double _grad_u[4],_u;
    double _testI,_grad_testI[4];
    double _phiJ,_grad_phiJ[4];

    string KernelName;
    string DofName;
    string MaterialName;
    int DofIndex;
    int DIM;


};


#endif //ASFEM_KERNELS_H
