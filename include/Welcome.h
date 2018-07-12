//*************************************************************
//* This file is part of the AsFem framework                  *
//* https://github.com/walkandthinker/AsFem                   *
//* All rights reserved, see COPYRIGHT for full restrictions  *
//* Licensed under LGPL 3.0, please see LICENSE for details   *
//* https://www.gnu.org/licenses/gpl-3.0.html                 *
//*************************************************************
//*** AsFem: A simple finite element method program         ***
//*** Author: Yang Bai @CopyRight                           ***
//*** Bug report: walkandthinker@gmail.com                  ***
//*** QQ group: 797998860                                   ***
//*************************************************************
//*** Created by Y. Bai on 12.01.18.                        ***
//*************************************************************

#ifndef ASFEM_WELCOME_H
#define ASFEM_WELCOME_H

#include "petsc.h"

void Welcome(double version=0.1)
{
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Welcome to use ASFEM                   ***\n");
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
