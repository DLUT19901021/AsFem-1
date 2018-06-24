//*****************************************************
//*** AsFem: A simple finite element method program ***
//*** Author: Yang Bai @CopyRight                   ***
//*** Bug report: walkandthinker@gmail.com          ***
//*** QQ group: 797998860                           ***
//*****************************************************
//*** This class handle Ax=F system equations       ***
//***  the final solution also be stored here       ***
//***  you must install PETSc library               ***
//***  this class depend on PETSc heavily!!!        ***
//*****************************************************

#ifndef ASFEM_EQUATIONSYSTEM_H
#define ASFEM_EQUATIONSYSTEM_H

#include <iostream>
#include <limits>
#include <vector>
#include <string>
#include <algorithm>
#include "petsc.h"

using namespace std;

class EquationSystem
{
public:
    EquationSystem(const int dofs,int dofspernode=1);

    PetscErrorCode Init();
    Mat AMATRIX;
    Vec RHS,dU,U;

    void AddSolutionNameAndIndex(string name,int order);
    void SetSolutionNameFromVector(vector<string> names,vector<int> orders);
    void SetSolutionName();

    string GetIthDofsName(const int i) const;

    void ReInitEquationSystem();
    PetscErrorCode Release();

    void PrintSolutionNameMap(string str="") const;

private:
    bool IsInit=false;
    PetscErrorCode ierr;
    int nDofs,nDofsPerNode;
    vector<string> solution_name_list;
    vector<int> solution_index_list;
    vector<pair<string,int>> solution_name_map;
    bool SolutionHasName=false;
    Vec U0;

};


#endif //ASFEM_EQUATIONSYSTEM_H
