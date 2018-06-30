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

#include "InputSystem/InputSystem.h"



InputSystem::InputSystem(int argc, char **argv)
{
    string str;
    if(argc>1)
    {
        str=argv[1];
        if(str.at(0)!='-'&&(str.at(1)!='i'||str.at(1)!='I'))
        {
            FileName=str;
        }
        else
        {
            FileName=str.substr(2,str.length());
        }
        in.open(FileName.c_str(),ios::in);
        if(!in.is_open())
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Wrong inp file name!!!                 ***\n");
            PetscFinalize();
            abort();
        }
    }
    else
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***----------------------------------------***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***--- Start input step -------------------***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***---   Input the inp file name:");
        cin >> FileName;
        in.open(FileName.c_str(), ios::in);

        while (!in.is_open())
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***----------------------------------------***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***--- Wrong input file name!!!            ***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***---   Input the inp file name:");
            cin >> FileName;
            in.open(FileName.c_str(), ios::in);
        }
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***----------------------------------------***\n");
    }
}

void InputSystem::CloseFile()
{
    if(in.is_open())
    {
        in.close();
    }
}