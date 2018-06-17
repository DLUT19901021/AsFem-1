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
#include "Mesh/Mesh1D.h"

using namespace std;

int main(int argc,char *argv[])
{
    double version=0.1;
    PetscErrorCode ierr;
    ierr=PetscInitialize(&argc,&argv,NULL,NULL);CHKERRQ(ierr);

    Welcome(version);


	Mesh1D mesh1D(0.0,1.0,10,"edge3");

	mesh1D.CreateMesh();
	//mesh1D.PrintMeshInfo();
    mesh1D.PrintMeshDetailInfo();

	ierr=PetscFinalize();CHKERRQ(ierr);
	return ierr;
} 
