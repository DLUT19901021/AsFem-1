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



using namespace std;

int main(int argc,char *argv[])
{
    double version=0.1;
    PetscErrorCode ierr;
    ierr=PetscInitialize(&argc,&argv,NULL,NULL);CHKERRQ(ierr);

    Welcome(version);

	ierr=PetscFinalize();CHKERRQ(ierr);
	return ierr;
} 
