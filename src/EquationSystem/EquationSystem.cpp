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

#include "EquationSystem/EquationSystem.h"

EquationSystem::EquationSystem(const int dofs,int dofspernode)
{

    if(dofs<2||dofs>MaxDofs)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: dofs number=%-8d is not supported in current version!\n",dofs);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        the maximum dofs of current version is:%-8d!\n",MaxDofs);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }
    if(dofspernode<1||dofspernode>10)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: nDofsPerNode=%d is invalid in current version!!!\n",dofspernode);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        max dofs per node<=10!!!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }

    nDofs=dofs;nDofsPerNode=dofspernode;
    solution_name_list.clear();
    solution_index_list.clear();
    IsInit=false;
    SolutionHasName=false;
}

EquationSystem::EquationSystem()
{

    nDofs=0;nDofsPerNode=0;
    solution_name_list.clear();
    solution_index_list.clear();
    IsInit=false;
    SolutionHasName=false;
}

PetscErrorCode EquationSystem::Init()
{
    if(IsInit)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Warning: equation system is already initialized!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***          AsFem will do nothing for you!\n");
        return 0;
    }

    ierr=MatCreate(PETSC_COMM_WORLD,&AMATRIX);CHKERRQ(ierr);
    ierr=MatSetSizes(AMATRIX,PETSC_DECIDE,PETSC_DECIDE,nDofs,nDofs);CHKERRQ(ierr);
    ierr=MatSetFromOptions(AMATRIX);CHKERRQ(ierr);
    ierr=MatSetUp(AMATRIX);CHKERRQ(ierr);

    // for vector
    ierr=VecCreate(PETSC_COMM_WORLD,&U0);CHKERRQ(ierr);
    ierr=VecSetSizes(U0,PETSC_DECIDE,nDofs);CHKERRQ(ierr);
    ierr=VecSetFromOptions(U0);CHKERRQ(ierr);
    ierr=VecSetUp(U0);CHKERRQ(ierr);

    ierr=VecDuplicate(U0,&RHS);CHKERRQ(ierr);
    ierr=VecDuplicate(U0,&dU);CHKERRQ(ierr);
    ierr=VecDuplicate(U0,&U);CHKERRQ(ierr);

    // initialize
    ierr=VecSet(U0,0.0);CHKERRQ(ierr);
    ierr=VecSet(RHS,0.0);CHKERRQ(ierr);

    nDofsPerNode=solution_name_list.size();
}

PetscErrorCode EquationSystem::Release()
{
    ierr=MatDestroy(&AMATRIX);CHKERRQ(ierr);

    ierr=VecDestroy(&U0);CHKERRQ(ierr);
    ierr=VecDestroy(&U);CHKERRQ(ierr);
    ierr=VecDestroy(&dU);CHKERRQ(ierr);
    ierr=VecDestroy(&RHS);CHKERRQ(ierr);

    solution_name_list.clear();
    solution_index_list.clear();
}

void EquationSystem::ReInitEquationSystem()
{
    // Warning: here this function should and must only
    //          initializing the AMATRIX and RHS,
    //          don't do operating on other vec!!!

    // so this should be called before the N-R iteration!

    MatZeroEntries(AMATRIX);// must be sure AMATRIX is already initalized
    VecSet(RHS,0.0);
}

void EquationSystem::AddSolutionNameAndIndex(string name, int order)
{
    if(SolutionHasName)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: solution already has name!!!\n",name,order);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        can't add name and order to the list!!!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }
    bool IsNameIn=false,IsOrderIn=false;
    for(unsigned int i=0;i<solution_name_list.size();i++)
    {
        if(name==solution_name_list[i])
        {
            IsNameIn=true;
            break;
        }
    }

    for(unsigned int i=0;i<solution_index_list.size();i++)
    {
        if(order==solution_index_list[i])
        {
            IsOrderIn=true;
            break;
        }
    }

    if(IsNameIn || IsOrderIn)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: name=%8s or order=%2d is already in the list!!\n",name.c_str(),order);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        unique name and oder should be given!!!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }
    else
    {
        solution_name_list.push_back(name);
        solution_index_list.push_back(order);
    }
}

void EquationSystem::SetSolutionNameFromVector(vector<string> names, vector<int> orders)
{
    if(SolutionHasName)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: solution already has name!!!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        can't add name and order to the list!!!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }

    bool IsNameUnique,IsOrderUnique;

    auto its=unique(names.begin(),names.end());
    IsNameUnique=(its==names.end());

    auto it=unique(orders.begin(),orders.end());
    IsOrderUnique=(it==orders.end());

    if((!IsNameUnique) || (!IsOrderUnique))
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: either name or oder list isn't unique\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        unique name and oder should be given!!!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }
    else if(names.size()!=orders.size())
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: len(names)=%2d isn't equal len(orders)=%2d\n",names.size(),orders.size());
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        equal length names_list and orders_list should be given!!!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }
    if(names.size()<nDofsPerNode || orders.size()<nDofsPerNode)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: neighter name or oder list is long enough to set solution name!!!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        please add more name and order values to the list!!!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }
    else
    {
        pair<string,int> temp;
        for(unsigned int i=0;i<orders.size();i++)
        {
            temp=make_pair(names[i],orders[i]);
            solution_name_map.push_back(temp);
        }
        SolutionHasName=true;
        nDofsPerNode=orders.size();
    }
}

void EquationSystem::SetSolutionName()
{
    if(SolutionHasName)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: solution already has name!!!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        can't add name and order to the list!!!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }
    if(solution_index_list.size()<nDofsPerNode || solution_name_list.size()<nDofsPerNode)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: neighter name or oder list is long enough to set solution name!!!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        please add more name and order values to the list!!!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }

    bool IsNameUnique,IsOrderUnique;

    auto its=unique(solution_name_list.begin(),solution_name_list.end());
    IsNameUnique=(its==solution_name_list.end());

    auto it=unique(solution_index_list.begin(),solution_index_list.end());
    IsOrderUnique=(it==solution_index_list.end());

    if((!IsNameUnique) || (!IsOrderUnique))
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: either name or oder list isn't unique\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        unique name and oder should be given!!!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }

    pair<string,int> temp;
    for(unsigned int i=0;i<solution_index_list.size();i++)
    {
        temp=make_pair(solution_name_list[i],solution_index_list[i]);
        solution_name_map.push_back(temp);
    }
    SolutionHasName=true;
    nDofsPerNode=solution_index_list.size();
}

//*****************************************
string EquationSystem::GetIthDofsName(const int i) const
{
    if(!SolutionHasName)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: solution has no name!!!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should set name for it first!!!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }

    if(i<1||i>nDofsPerNode)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: i=%d is invalid in for a dof index!!!\n",i);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        1 to %2d is expected!!!\n",nDofsPerNode);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }


    for(unsigned int j=0;j<solution_name_map.size();j++)
    {
        if(solution_name_map[j].second==i)
        {
            return solution_name_map[j].first;
        }
    }
}

void EquationSystem::PrintSolutionNameMap(string str) const
{
    if(!SolutionHasName)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: solution has no name!!!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should set name for it first!!!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }

    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
    if(str.size()>1)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** str= %s\n",str.c_str());
    }
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Solution system information:           ***\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   number of solution=%2d                ***\n",nDofsPerNode);

    string name;
    for(int i=1;i<=nDofsPerNode;i++)
    {
        name=GetIthDofsName(i);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   solution[%2d]<------>%8s         ***\n",i,name.c_str());
    }

    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");

}

