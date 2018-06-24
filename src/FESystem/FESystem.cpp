//*****************************************************
//*** AsFem: A simple finite element method program ***
//*** Author: Yang Bai @CopyRight                   ***
//*** Bug report: walkandthinker@gmail.com          ***
//*** QQ group: 797998860                           ***
//*****************************************************
//*** This class defines the core parts of AsFem    ***
//***  all the kernels, materials, boundarys        ***
//***  all the fem calculation related code         ***
//***  should be put here                           ***
//*****************************************************

#include "FESystem/FESystem.h"

FESystem::FESystem(Mesh &mesh,
                   DofHandler &dofHandler,
                   EquationSystem &equationSystem,
                   SolverSystem &solverSystem)
{
    KernelList.clear();
    MaterialKernelList.clear();
    BoundaryKernelList.clear();
    AuxKernelList.clear();
}

void FESystem::Run()
{
    AssembleFESystem();
    //ApplyBoundaryCondtion();
}

void FESystem::AssembleFESystem()
{
    int e,i,j,nNodesPerElmt;
    double k[27][27],rhs[27];
    double elCoords[27][4],elU={0.0};
    PetscInt elConn[27]={0};
}