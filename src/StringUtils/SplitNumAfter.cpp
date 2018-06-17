//*****************************************************
//*** AsFem: A simple finite element method program ***
//*** Author: Yang Bai                              ***
//*** Bug report: walkandthinker@gmail.com          ***
//*** QQ group: 797998860                           ***
//*****************************************************
//*** Function to split number from string after pos***
//*****************************************************

#include "StringUtils.h"

vector<double> SplitNumAfter(string instr,int pos)
{
    string str;
    str=instr.substr(pos,instr.size()-pos+1);
    return SplitNum(str);
}