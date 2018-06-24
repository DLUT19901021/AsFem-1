//*****************************************************
//*** AsFem: A simple finite element method program ***
//*** Author: Yang Bai @CopyRight                   ***
//*** Bug report: walkandthinker@gmail.com          ***
//*** QQ group: 797998860                           ***
//*****************************************************
//*** This class defines the input file read system ***
//***  all the input file informaiton should be     ***
//***  read by this class                           ***
//*****************************************************

#ifndef ASFEM_INPUTSYSTEM_H
#define ASFEM_INPUTSYSTEM_H

#include <iostream>
#include <iomanip>

#include "petsc.h"

// For AsFem's own header file
#include "Mesh/Mesh.h"

#include "EquationSystem/EquationSystem.h"

using namespace std;

class InputSystem
{
public:
    InputSystem(int argc,char *argv[]);

    bool ReadMesh(Mesh &mesh);
    bool ReadDofsName(EquationSystem &equationSystem);
    bool ReadBoundaryCondition();
    bool
private:
    string FileName;

};


#endif //ASFEM_INPUTSYSTEM_H
