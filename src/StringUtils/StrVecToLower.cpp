//*****************************************************
//*** AsFem: A simple finite element method program ***
//*** Author: Yang Bai                              ***
//*** Bug report: walkandthinker@gmail.com          ***
//*** QQ group: 797998860                           ***
//*****************************************************
//*** Function convert string vector to lower case  ***
//*****************************************************

#include <algorithm>
#include "StringUtils/StringUtils.h"

vector<string> StrVecToLower(vector<string> instrvec)
{
    vector<string> outstrvec=instrvec;
    for(unsigned int i=0;i<outstrvec.size();i++)
    {
        string str=outstrvec[i];
        transform(str.begin(),str.end(),str.begin(),::tolower);
        outstrvec[i]=str;
    }
    return outstrvec;
}
