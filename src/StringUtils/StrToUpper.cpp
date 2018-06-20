//*****************************************************
//*** AsFem: A simple finite element method program ***
//*** Author: Yang Bai                              ***
//*** Bug report: walkandthinker@gmail.com          ***
//*** QQ group: 797998860                           ***
//*****************************************************
//*** Function to convert string to upper case      ***
//***  i.e: "yang"-->"YANG"                         ***
//*****************************************************

#include <algorithm>
#include "StringUtils/StringUtils.h"

string StrToUpper(string instr)
{
    string outstr=instr;
    transform(outstr.begin(),outstr.end(),outstr.begin(),::toupper);
    return outstr;
}
