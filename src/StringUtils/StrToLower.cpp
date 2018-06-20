//*****************************************************
//*** AsFem: A simple finite element method program ***
//*** Author: Yang Bai                              ***
//*** Bug report: walkandthinker@gmail.com          ***
//*** QQ group: 797998860                           ***
//*****************************************************
//*** Function to remove space of input string      ***
//***  i.e: "Yang"-->"yang"                         ***
//*****************************************************

#include <algorithm>
#include "StringUtils/StringUtils.h"

string StrToLower(string instr)
{
    string outstr=instr;
    transform(outstr.begin(),outstr.end(),outstr.begin(),::tolower);
    return outstr;
}
