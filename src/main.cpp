//*****************************************************
//*** AsFem: A simple finite element method program ***
//*** Author: Yang Bai                              ***
//*** Bug report: walkandthinker@gmail.com          ***
//*** QQ group: 797998860                           ***
//*****************************************************

// Include standard header file from C++
#include <iostream>
// Include required header file from PETSc
#include "petsc.h"

// Include AsFem's own header file
#include "Welcome.h"

#include "Mesh/Mesh2D.h"
#include "Mesh/Mesh.h"

#include "DofHandler/DofHandler.h"

#include "EquationSystem/EquationSystem.h"

#include "SolverSystem/SolverSystem.h"

#include "FESystem/FESystem.h"

using namespace std;

int main(int argc,char *argv[])
{
    double version=0.1;
    PetscErrorCode ierr;
    ierr=PetscInitialize(&argc,&argv,NULL,NULL);CHKERRQ(ierr);

    Welcome(version);

    Mesh2D mesh2D(0.,1.,0.,1.,2,2,"quad8");
    mesh2D.CreateMesh();

    Mesh mesh(2,8,2);

    mesh.Add2DMesh(mesh2D);
    mesh.Init();

    mesh.PrintMeshInfo();

    mesh.PrintMeshDetailInfo();

    DofHandler dofHandler;
    dofHandler.CreateLocalToGlobalDofMap(mesh);
    dofHandler.PrintDofMap();

    EquationSystem equationSystem(mesh.GetDofsNum(),mesh.GetDofsNumPerNode());
    equationSystem.Init();
    equationSystem.AddSolutionNameAndIndex("disp_x",1);
    equationSystem.AddSolutionNameAndIndex("disp_y",2);
    equationSystem.SetSolutionName();
    equationSystem.PrintSolutionNameMap();


    SolverSystem solver(PETSC_COMM_WORLD);


    FESystem feSystem(mesh,dofHandler,equationSystem,solver);

    feSystem.Run();

	ierr=PetscFinalize();CHKERRQ(ierr);
	return ierr;
} 
