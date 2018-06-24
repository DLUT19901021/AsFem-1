//*****************************************************
//*** AsFem: A simple finite element method program ***
//*** Author: Yang Bai @CopyRight                   ***
//*** Bug report: walkandthinker@gmail.com          ***
//*** QQ group: 797998860                           ***
//*****************************************************
//*** This class define the solver system of AsFem  ***
//***  solver based on PETSc's ksp solver           ***
//*****************************************************

#include "SolverSystem/SolverSystem.h"

SolverSystem::SolverSystem(MPI_Comm comm)
{
    KSPCreate(comm,&ksp);

    KSPSetTolerances(ksp,
                     PETSC_DEFAULT,
                     PETSC_DEFAULT,
                     PETSC_DEFAULT,
                     1000000);

    KSPSetFromOptions(ksp);
    /*
    KSPSetUp(ksp);// TODO: it seems call this can lead to some PETSs errors?
     */
}

bool SolverSystem::Solve(Mat &A, Vec &dU, Vec &RHS)
{
    ierr=KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
    ierr=KSPSolve(ksp,RHS,dU);CHKERRQ(ierr);
}

int SolverSystem::GetKSPIterations() const
{
    PetscInt iters;
    KSPGetIterationNumber(ksp,&iters);
    return iters;
}

void SolverSystem::Release()
{
    KSPDestroy(&ksp);
}
