//*****************************************************
//*** AsFem: A simple finite element method program ***
//*** Author: Yang Bai @CopyRight                   ***
//*** Bug report: walkandthinker@gmail.com          ***
//*** QQ group: 797998860                           ***
//*****************************************************
//*** We read the name of each dofs  from inp file  ***
//*****************************************************

#include "InputSystem/InputSystem.h"

bool InputSystem::ReadDofsName(EquationSystem &equationSystem)
{
    //*****************************
    //*** Format:
    //***     [variables]
    //***       name=ux uy
    //***     []
    //*****************************
    string str,line,line0="variables";
    vector<string> namelist;
    vector<int> indexlist;
    int linenum,blockstartlinenum,i;


    linenum=0;
    namelist.clear();
    indexlist.clear();
    if(IsBracketMatch(in,line0,blockstartlinenum))
    {
        GotoLine(in,blockstartlinenum);linenum=blockstartlinenum-1;
        getline(in,line);linenum+=1;// read [variables]

        // read name of each dofs
        getline(in,line);linenum+=1;
        if(line.find("name=")!=string::npos)
        {
            i=line.find("=");
            str=line.substr(i+1,line.size()-i+1);
            namelist=SplitStr(str,' ');
        }
        else
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: 'name=' can't found in line=%-3d ***\n",linenum);
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given name=ux uy..   ***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            return false;
        }
    }
    if(namelist.size()<1)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find dofs name in line=%-3d***\n",linenum);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given name=ux uy..   ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        return false;
    }
    else
    {
        for(i=0;i<namelist.size();i++)
        {
            indexlist.push_back(i+1);
        }

        equationSystem.SetSolutionNameFromVector(namelist,indexlist);
        return true;
    }

}

