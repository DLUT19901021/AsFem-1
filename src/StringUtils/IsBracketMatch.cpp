//*****************************************************
//*** AsFem: A simple finite element method program ***
//*** Author: Yang Bai @CopyRight                   ***
//*** Bug report: walkandthinker@gmail.com          ***
//*** QQ group: 797998860                           ***
//*****************************************************
//*** We read mesh from input file                  ***
//*****************************************************

#include "StringUtils/StringUtils.h"

bool IsBracketMatch(ifstream &in,string &bracketstr,int &startline)
{
    string line,str;
    bool FoundHead=false,FoundEnd=false;
    string line0='['+bracketstr+']';
    int endline=0;
    startline=0;
    int linenum=0;
    in.clear();
    in.seekg(0, ios::beg);
    while(!in.eof())
    {
        getline(in,line);linenum+=1;
        str=RemoveSpace(line);
        if(str.compare(line0)==0)
        {
            FoundHead=true;
            startline=linenum;
            while(!in.eof())
            {
                getline(in,line);linenum+=1;
                str=RemoveSpace(line);
                if(str.compare("[]")==0)
                {
                    FoundEnd= true;
                    endline=linenum;
                    break;
                }
            }
        }
    }

    if((FoundHead && FoundEnd)&&(endline>startline))
    {
        //cout<<"start="<<startline<<":end="<<endline<<endl;
        return true;
    }
    else
    {
        return false;
    }
}
