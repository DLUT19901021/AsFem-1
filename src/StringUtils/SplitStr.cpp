//*****************************************************
//*** AsFem: A simple finite element method program ***
//*** Author: Yang Bai                              ***
//*** Bug report: walkandthinker@gmail.com          ***
//*** QQ group: 797998860                           ***
//*****************************************************
//*** Function to split str via input char symbol   ***
//*****************************************************

#include <sstream>
#include "StringUtils.h"

vector<string> SplitStr(string instr,char symbol)
{
    vector<string> outstr;
    stringstream ss(instr);
    string tok;

    while(getline(ss,tok,symbol))
    {
        outstr.push_back(tok);
    }
    return outstr;
}
