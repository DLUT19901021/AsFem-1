//*****************************************************
//*** AsFem: A simple finite element method program ***
//*** Author: Yang Bai @CopyRight                   ***
//*** Bug report: walkandthinker@gmail.com          ***
//*** QQ group: 797998860                           ***
//*****************************************************
//*** Reset the ifstream pointer to expected line   ***
//*****************************************************

#include "StringUtils/StringUtils.h"

void GotoLine(ifstream &in,int linenum)
{
    string line;
    in.clear();
    in.seekg(ios::beg);
    for(int i=0;i<linenum-1;i++)
    {
        getline(in,line);
    }

}