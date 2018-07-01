//*****************************************************
//*** AsFem: A simple finite element method program ***
//*** Author: Yang Bai @CopyRight                   ***
//*** Bug report: walkandthinker@gmail.com          ***
//*** QQ group: 797998860                           ***
//*****************************************************
//*** This class defines the core parts of AsFem    ***
//***  basical and elemental calculation of FEM     ***
//***  element loop and matrix assembly code        ***
//***  should be put here                           ***
//*****************************************************

#ifndef ASFEM_FE_H
#define ASFEM_FE_H

#include <iostream>
#include <string>
#include <vector>
#include "petsc.h"

// AsFem's own header file
#include "FE/Qrule.h"
#include "Mesh/Mesh.h"
#include "DofHandler/DofHandler.h"
#include "EquationSystem/EquationSystem.h"

#include "Kernels/Kernel.h"

using namespace std;

class FE
{
public:
    FE(Mesh &mesh,DofHandler &dofHandler);

    void InitLocalKR();
    void FormKR(Mesh &mesh,DofHandler &dofHandler,EquationSystem &equationSystem);
    void QpLoop(const double (&elCoords)[27][4],const double (&elU)[27][2],
                double (&k)[270][270],double (&rhs)[270]);
    void AssembleLocalToGlobal();

    void KernelListLoop();


private:
    int DIM,nNodesPerElmt;
    int Lint,qp;
    string MeshType;
    double JxW;
    Qrule *qrule;

    VecScatter Uscatter,Vscatter;
    Vec Useq,Vseq;


private:
    void ComputeResidual(vector<Kernel> &KerenlsList,);
    void ComputeJacobian(vector<Kernel> &KerenlsList);


private:
    void Shp1D(const int &ndim,const int &nnodes,const double &xi,const double (&Coords)[27][4],double (&shp)[27][4],double &DetJac,bool IsNature);
    void Shp2D(const int &ndim,const int &nnodes,const double &xi,const double &eta,const double (&Coords)[27][4],double (&shp)[27][4],double &DetJac,bool IsNature);
    void Shp3D(const int &ndim,const int &nnodes,const double &xi,const double &eta,const double &zeta,const double (&Coords)[27][4],double (&shp)[27][4],double &DetJac,bool IsNature);


};


#endif //ASFEM_FE_H
