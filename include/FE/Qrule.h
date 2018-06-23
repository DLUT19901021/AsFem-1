//*****************************************************
//*** AsFem: A simple finite element method program ***
//*** Author: Yang Bai                              ***
//*** Bug report: walkandthinker@gmail.com          ***
//*** QQ group: 797998860                           ***
//*****************************************************
//*** Define class for gauss integration            ***
//***   generate gauss point and weight for inte    ***
//***                                               ***
//*****************************************************

#ifndef ASFEM_QRULE_H
#define ASFEM_QRULE_H

#include <iostream>
#include <vector>
#include <cmath>

#include "petsc.h"

using namespace std;

class Qrule
{
public:
    Qrule(int dim,int n=2);

    int GetQPointsNum() const {return nGaussPoints;}
    double GetComponent(int i,int j) const;

    void PrintGaussPointInfo() const;

private:
    vector<double> weight,coords;
    int nGaussPoints,ngp,nDims;

    void GenerateGaussPoints();

    const double p1xi[2]={-0.577350269189625764509148780502,
                           0.577350269189625764509148780502};
    const double p1w[2]={1.0,1.0};

    const double p2xi[3]={-0.774596669241483377035853079956,
                           0.0,
                           0.774596669241483377035853079956};
    const double p2w[3]={0.555555555555555555555555555556,
                         0.888888888888888888888888888889,
                         0.555555555555555555555555555556};

    const double p3xi[4]={-0.861136311594052575223946488893,
                          -0.339981043584856264802665759103,
                           0.339981043584856264802665759103,
                           0.861136311594052575223946488893};
    const double p3w[4]={0.347854845137453857373063949222,
                         0.652145154862546142626936050778,
                         0.652145154862546142626936050778,
                         0.347854845137453857373063949222};                                                            
};

#endif