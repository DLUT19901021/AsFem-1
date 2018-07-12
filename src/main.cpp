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

#include <iostream>

#include "petsc.h"

//* For AsFem's own header file
#include "Welcome.h"
#include "Mesh/Mesh.h"

using namespace std;


int main(int args,char *argv[])
{
    double version=0.1;

    PetscErrorCode ierr;

    //* Initializing
    ierr=PetscInitialize(&args,&argv,NULL,NULL);CHKERRQ(ierr);

    Welcome(version);

    Mesh mesh;

    mesh.SetDim(3);
    mesh.SetNx(1);

    mesh.SetXmax(1.0);
    mesh.SetXmin(0.0);

    mesh.SetYmin(1.0);
    mesh.SetYmax(2.0);

    mesh.SetZmax(2.0);
    mesh.SetZmin(.5);


    ierr=PetscFinalize();CHKERRQ(ierr);
    return ierr;

}
