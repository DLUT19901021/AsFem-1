//*****************************************************
//*** AsFem: A simple finite element method program ***
//*** Author: Yang Bai                              ***
//*** Bug report: walkandthinker@gmail.com          ***
//*** QQ group: 797998860                           ***
//*****************************************************
//*** Several function to operator string           ***
//***  i.e: upper to lower                          ***
//***       low to upper                            ***
//***       split number from string                ***
//*****************************************************

#ifndef ASFEM_STRINGUTILS_H
#define ASFEM_STRINGUTILS_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

string RemoveSpace(string instr);
string StrToLower(string instr);
vector<string> StrVecToLower(vector<string> instrvec);

string StrToUpper(string instr);
vector<string> StrVecToUpper(vector<string> instrvec);

vector<string> SplitStr(string instr,char symbol);

vector<double> SplitNum(string instr); //
vector<double> SplitNumAfter(string instr,int pos);

bool IsBracketMatch(ifstream &in,string &bracketstr,int &startline);

void GotoLine(ifstream &in,int linenum);

#endif //ASFEM_STRINGUTILS_H
