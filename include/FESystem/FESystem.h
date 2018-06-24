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

#ifndef ASFEM_FESYSTEM_H
#define ASFEM_FESYSTEM_H

#include "Mesh/Mesh.h"
#include "DofHandler/DofHandler.h"
#include "EquationSystem/EquationSystem.h"
#include "SolverSystem/SolverSystem.h"

class FESystem
{
public:
    FESystem(Mesh &mesh,
             DofHandler &dofHandler,
             EquationSystem &equationSystem,
             SolverSystem &solverSystem);

    enum {Steady,Transient};

    void RegisterKernel();// If you have your own kernel, register it here

    void RegisterMaterial();// If you have your own material, register it here
    void RegisterBoundary();// If you have your own boundary condition, register it here

    void RegisterAuxKernel();

    void ComputeResidual();// calculate local RHS
    void ComputeJacobian();// calculate local K
    void ComputeMaterial();// calcuate local material properties
    void ComputeAuxVariable();//

    void AssembleFESystem();// assemble local K and RHS to FESystem

    void ApplyBoundaryCondtion();

    void Run();

private:
    vector<int> KernelList;
    vector<int> MaterialKernelList;
    vector<int> BoundaryKernelList;
    vector<int> AuxKernelList;
};


#endif //ASFEM_FESYSTEM_H
