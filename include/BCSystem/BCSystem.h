//*****************************************************
//*** AsFem: A simple finite element method program ***
//*** Author: Yang Bai @CopyRight                   ***
//*** Bug report: walkandthinker@gmail.com          ***
//*** QQ group: 797998860                           ***
//*****************************************************
//*** This class defines the boundary conditon      ***
//***  all the boundary information read from       ***
//***  input file should be stored here             ***
//*****************************************************

#ifndef ASFEM_BCSYSTEM_H
#define ASFEM_BCSYSTEM_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

using namespace std;

class BCSystem
{
public:
    BCSystem();

private:
    vector<pair<string,double>> DirichletBCList;
    vector<pair<string,double>> NeumannBCList;
    vector<pair<string,int>> RobinBCList;
};

#endif //ASFEM_BCSYSTEM_H
