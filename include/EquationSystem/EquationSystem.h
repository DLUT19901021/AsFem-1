//*****************************************************
//*** AsFem: A simple finite element method program ***
//*** Author: Yang Bai                              ***
//*** Bug report: walkandthinker@gmail.com          ***
//*** QQ group: 797998860                           ***
//*****************************************************
//*** Define class for system equations             ***
//***   where Ax=F required info is defined         ***
//***                                               ***
//*****************************************************

#ifndef ASFEM_EQUATIONSYSTEM_H
#define ASFEM_EQUATIONSYSTEM_H

#include <iostream>
#include "petsc.h"

using namespace std;

class EquationSystem
{
public:
    EquationSystem(int ndofs);
    Mat AMATRX;
    Vec RHS;

    void Initilaze();

private:
    bool IsInit;
    int nDofs;
    Vec U,U0;
    KSP ksp;
};

#endif