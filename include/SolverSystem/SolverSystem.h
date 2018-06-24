//*****************************************************
//*** AsFem: A simple finite element method program ***
//*** Author: Yang Bai @CopyRight                   ***
//*** Bug report: walkandthinker@gmail.com          ***
//*** QQ group: 797998860                           ***
//*****************************************************
//*** This class define the solver system of AsFem  ***
//***  solver based on PETSc's ksp solver           ***
//*****************************************************

#ifndef ASFEM_SOLVERSYSTEM_H
#define ASFEM_SOLVERSYSTEM_H

#include <iostream>

#include "petsc.h"

#include "EquationSystem/EquationSystem.h"

class SolverSystem
{
public:
    SolverSystem(MPI_Comm comm=PETSC_COMM_WORLD);

    bool Solve(Mat &A,Vec &dU,Vec &RHS);

    int GetKSPIterations() const;



    void Release();

private:
    KSP ksp;
    PetscErrorCode ierr;
};

#endif //ASFEM_SOLVERSYSTEM_H
