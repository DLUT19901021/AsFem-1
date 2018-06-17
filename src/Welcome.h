//*****************************************************
//*** AsFem: A simple finite element method program ***
//*** Author: Yang Bai                              ***
//*** Bug report: walkandthinker@gmail.com          ***
//*** QQ group: 797998860                           ***
//*****************************************************
//*** Print the welcome message in console          ***
//***  show a short introduction for AsFem          ***
//*****************************************************

#ifndef ASFEM_WELCOME_H
#define ASFEM_WELCOME_H

#include "petsc.h"

void Welcome(double version=0.1)
{
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Welcome to use AsFem ^_^               ***\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** A Simple Finite Element Method program ***\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Version: %6.1f                        ***\n",version);
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Author: Yang Bai @ CopyRight           ***\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Email: walkandthinker@gmail.com        ***\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** License: GPL-3.0                       ***\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** QQ group: 797998860                    ***\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Please feel free to use and discuss    ***\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
}

#endif //ASFEM_WELCOME_H
