//*****************************************************
//*** AsFem: A simple finite element method program ***
//*** Author: Yang Bai                              ***
//*** Bug report: walkandthinker@gmail.com          ***
//*** QQ group: 797998860                           ***
//*****************************************************
//*** Function to remove space of input string      ***
//***  i.e: "yang bai"-->"yangbai"                  ***
//*****************************************************

#include "StringUtils/StringUtils.h"

string RemoveSpace(string instr)
{
    if(instr.size()<=1)
    {
        return instr;
    }
    unsigned int i,length;
    char ch[1000];
    length=0;
    for(i=0;i<instr.size();i++)
    {
        if(instr.at(i)!=' ')
        {
            ch[length]=instr.at(i);
            length+=1;
        }
    }

    string outstr;
    outstr.clear();
    for(i=0;i<length;i++)
    {
        if(ch[i]=='\n')
        {
            break;
        }
        outstr+=ch[i];
    }
    return outstr;
}