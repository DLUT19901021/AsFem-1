//*****************************************************
//*** AsFem: A simple finite element method program ***
//*** Author: Yang Bai @CopyRight                   ***
//*** Bug report: walkandthinker@gmail.com          ***
//*** QQ group: 797998860                           ***
//*****************************************************
//*** Define 2D shape functions for N and dNdx      ***
//*****************************************************

#include "FE/FE.h"

void FE::Shp3D(const int &ndim,
               const int &nnodes,
               const double &xi,
               const double &eta,
               const double &zeta,
               const double (&Coords)[27][4],
               double (&shp)[27][4],
               double &DetJac,
               bool IsNature)
{
    if(ndim!=3)
    {
        cout<<"dim="<<ndim<<" is invalid for 3D shape function!!!"<<endl;
        abort();
    }

    int i;
    double dxdxi, dydxi, dzdxi;
    double dxdeta, dydeta, dzdeta;
    double dxdzeta,dydzeta,dzdzeta;
    double xi0,eta0,zeta0;
    double xi1,eta1,zeta1;



    if (nnodes!=8&&nnodes!=20&&nnodes!=27)
    {
        cout<<"*** Error:Wrong nodes on current element(nNodes="<<nnodes<<")!!!" << endl;
        cout<<"*** AsFem only support 3D-8,20,27 nodes hex-element!!!"<<endl;
        cout<<"*** Error happen in Shp3D!!!"<<endl;
        abort();
    }

    //TODO: change Matrix3d to double type array
    /*
    Matrix3d Jac,XJac;
    //shp.setZero();

    if(nnodes==8)
    {
        shp[0][0]= (1-xi)*(1-eta)*(1-zeta)/8.0;
        shp[0][1]=-(1-eta)*(1-zeta)/8.0;
        shp[0][2]=-(1-xi)*(1-zeta)/8.0;
        shp[0][3]=-(1-xi)*(1-eta)/8.0;

        shp[1][0]= (1+xi)*(1-eta)*(1-zeta)/8.0;
        shp[1][1]= (1-eta)*(1-zeta)/8.0;
        shp[1][2]=-(1+xi)*(1-zeta)/8.0;
        shp[1][3]=-(1+xi)*(1-eta)/8.0;

        shp[2][0]= (1+xi)*(1+eta)*(1-zeta)/8.0;
        shp[2][1]= (1+eta)*(1-zeta)/8.0;
        shp[2][2]= (1+xi)*(1-zeta)/8.0;
        shp[2][3]=-(1+xi)*(1+eta)/8.0;

        shp[3][0]= (1-xi)*(1+eta)*(1-zeta)/8.0;
        shp[3][1]=-(1+eta)*(1-zeta)/8.0;
        shp[3][2]= (1-xi)*(1-zeta)/8.0;
        shp[3][3]=-(1-xi)*(1+eta)/8.0;

        shp[4][0]= (1-xi)*(1-eta)*(1+zeta)/8.0;
        shp[4][1]=-(1-eta)*(1+zeta)/8.0;
        shp[4][2]=-(1-xi)*(1+zeta)/8.0;
        shp[4][3]= (1-xi)*(1-eta)/8.0;

        shp[5][0]= (1+xi)*(1-eta)*(1+zeta)/8.0;
        shp[5][1]= (1-eta)*(1+zeta)/8.0;
        shp[5][2]=-(1+xi)*(1+zeta)/8.0;
        shp[5][3]= (1+xi)*(1-eta)/8.0;

        shp[6][0]= (1+xi)*(1+eta)*(1+zeta)/8.0;
        shp[6][1]= (1+eta)*(1+zeta)/8.0;
        shp[6][2]= (1+xi)*(1+zeta)/8.0;
        shp[6][3]= (1+xi)*(1+eta)/8.0;

        shp[7][0]= (1-xi)*(1+eta)*(1+zeta)/8.0;
        shp[7][1]=-(1+eta)*(1+zeta)/8.0;
        shp[7][2]= (1-xi)*(1+zeta)/8.0;
        shp[7][3]= (1-xi)*(1+eta)/8.0;
    }
    else if(nnodes==20)
    {
        double XI[]={0.0,-1.0,1.0,1.0,-1.0,-1.0,1.0,1.0,-1.0};
        double ETA[]={0.0,-1.0,-1.0,1.0,1.0,-1.0,-1.0,1.0,1.0};
        double ZETA[]={0.0,-1.0,-1.0,-1.0,-1.0,1.0,1.0,1.0,1.0};
        double M1[]={0.0,-1.0,1.0,-1.0,1.0};
        double M2[]={0.0,-1.0,-1.0,1.0,1.0};
        // for 8-corner node
        for(i=1;i<=8;++i)
        {
            xi0=XI[i]*xi;
            eta0=ETA[i]*eta;
            zeta0=ZETA[i]*zeta;

            xi1=0.5+0.5*xi0;
            eta1=0.5+0.5*eta0;
            zeta1=0.5+0.5*zeta0;

            shp[i-1][0]=xi1*eta1*zeta1*(xi0+eta0+zeta0-2.0);
            shp[i-1][1]=eta1*zeta1*(xi+0.5*XI[i]*(eta0+zeta0-1.0));
            shp[i-1][2]=xi1*zeta1*(eta+0.5*ETA[i]*(zeta0+xi0-1.0));
            shp[i-1][3]=xi1*eta1*(zeta+0.5*ZETA[i]*(xi0+eta0-1.0));
        }
        // for middle point
        int ix[]={0,9,11,13,15};
        int iy[]={0,12,10,16,14};
        int iz[]={0,17,18,20,19};
        for(i=1;i<=4;++i)
        {
            xi1=(1.0-xi*xi)*0.25;
            eta1=1.0+M1[i]*eta;
            zeta1=1.0+M2[i]*zeta;
            shp[ix[i]-1][0]=xi1*eta1*zeta1;
            shp[ix[i]-1][1]=-0.5*xi*eta1*zeta1;
            shp[ix[i]-1][2]=xi1*M1[i]*zeta1;
            shp[ix[i]-1][3]=xi1*eta1*M2[i];

            xi1=1.0+M1[i]*xi;
            eta1=(1.0-eta*eta)*0.25;
            shp[iy[i]-1][0]= xi1*eta1*zeta1;
            shp[iy[i]-1][1]= M1[i]*eta1*zeta1;
            shp[iy[i]-1][2]=-0.5*xi1*eta*zeta1;
            shp[iy[i]-1][3]= xi1*eta1*M2[i];

            eta1=1.0+M2[i]*eta;
            zeta1=(1.0-zeta*zeta)*0.25;
            shp[iz[i]-1][0]= xi1*eta1*zeta1;
            shp[iz[i]-1][1]= M1[i]*eta1*zeta1;
            shp[iz[i]-1][2]= xi1*M2[i]*zeta1;
            shp[iz[i]-1][3]=-0.5*xi1*eta1*zeta;
        }
    }
    else if(nnodes==27)
    {
        /*
        int ix[]={0,1,2,2,1,1,2,2,1,3,2,3,1,3,2,3,1,1,2,2,1,3,3,1,2,3,3,3};
        int iy[]={0,1,1,2,2,1,1,2,2,1,3,2,3,1,3,2,3,1,1,2,2,3,3,3,3,1,2,3};
        int iz[]={0,1,1,1,1,2,2,2,2,1,1,1,1,2,2,2,2,3,3,3,3,1,2,3,3,3,3,3};
        double xi0[]  ={0.0,    0.5*xi*(xi-1.0),    0.5*xi*(xi+1.0),    1.0-xi*xi};
        double eta0[] ={0.0,  0.5*eta*(eta-1.0),  0.5*eta*(eta+1.0),  1.0-eta*eta};
        double zeta0[]={0.0,0.5*zeta*(zeta-1.0),0.5*zeta*(zeta+1.0),1.0-zeta*zeta};
        double xi1[]  ={0.0,  xi-0.5,  xi+0.5,  -2.0*xi};
        double eta1[] ={0.0, eta-0.5, eta+0.5, -2.0*eta};
        double zeta1[]={0.0,zeta-0.5,zeta+0.5,-2.0*zeta};
        for(i=1;i<=27;++i)
        {
            shp(i-1,0)=xi0[ix[i]]*eta0[iy[i]]*zeta0[iz[i]];
            shp(i-1,1)=xi1[ix[i]]*eta0[iy[i]]*zeta0[iz[i]];
            shp(i-1,2)=xi0[ix[i]]*eta1[iy[i]]*zeta0[iz[i]];
            shp(i-1,3)=xi0[ix[i]]*eta0[iy[i]]*zeta1[iz[i]];
        }
         */

        /*
        shp[0][0]= xi*(xi-1)*eta*(eta-1)*zeta*(zeta-1)/8.0;
        shp[0][1]= (2*xi-1)*eta*(eta-1)*zeta*(zeta-1)/8.0;
        shp[0][2]= xi*(xi-1)*(2*eta-1)*zeta*(zeta-1)/8.0;
        shp[0][3]= xi*(xi-1)*eta*(eta-1)*(2*zeta-1)/8.0;

        shp[1][0]= xi*(xi+1)*eta*(eta-1)*zeta*(zeta-1)/8.0;
        shp[1][1]= (2*xi+1)*eta*(eta-1)*zeta*(zeta-1)/8.0;
        shp[1][2]= xi*(xi+1)*(2*eta-1)*zeta*(zeta-1)/8.0;
        shp[1][3]= xi*(xi+1)*eta*(eta-1)*(2*zeta-1)/8.0;

        shp[2][0] = xi*(xi+1)*eta*(eta+1)*zeta*(zeta-1)/8.0;
        shp[2][1] = (2*xi+1)*eta*(eta+1)*zeta*(zeta-1)/8.0;
        shp[2][2] = xi*(xi+1)*(2*eta+1)*zeta*(zeta-1)/8.0;
        shp[2][3] = xi*(xi+1)*eta*(eta+1)*(2*zeta-1)/8.0;

        shp[3][0]= xi*(xi-1)*eta*(eta+1)*zeta*(zeta-1)/8.0;
        shp[3][1]= (2*xi-1)*eta*(eta+1)*zeta*(zeta-1)/8.0;
        shp[3][2]= xi*(xi-1)*(2*eta+1)*zeta*(zeta-1)/8.0;
        shp[3][3]= xi*(xi-1)*eta*(eta+1)*(2*zeta-1)/8.0;

        shp[4][0]= xi*(xi-1)*eta*(eta-1)*zeta*(zeta+1)/8.0;
        shp[4][1]= (2*xi-1)*eta*(eta-1)*zeta*(zeta+1)/8.0;
        shp[4][2]= xi*(xi-1)*(2*eta-1)*zeta*(zeta+1)/8.0;
        shp[4][3]= xi*(xi-1)*eta*(eta-1)*(2*zeta+1)/8.0;

        shp[5][0]= xi*(xi+1)*eta*(eta-1)*zeta*(zeta+1)/8.0;
        shp[5][1]= (2*xi+1)*eta*(eta-1)*zeta*(zeta+1)/8.0;
        shp[5][2]= xi*(xi+1)*(2*eta-1)*zeta*(zeta+1)/8.0;
        shp[5][3]= xi*(xi+1)*eta*(eta-1)*(2*zeta+1)/8.0;

        shp[6][0]= xi*(xi+1)*eta*(eta+1)*zeta*(zeta+1)/8.0;
        shp[6][1]= (2*xi+1)*eta*(eta+1)*zeta*(zeta+1)/8.0;
        shp[6][2]= xi*(xi+1)*(2*eta+1)*zeta*(zeta+1)/8.0;
        shp[6][3]= xi*(xi+1)*eta*(eta+1)*(2*zeta+1)/8.0;

        shp[7][0]= xi*(xi-1)*eta*(eta+1)*zeta*(zeta+1)/8.0;
        shp[7][1]= (2*xi-1)*eta*(eta+1)*zeta*(zeta+1)/8.0;
        shp[7][2]= xi*(xi-1)*(2*eta+1)*zeta*(zeta+1)/8.0;
        shp[7][3]= xi*(xi-1)*eta*(eta+1)*(2*zeta+1)/8.0;

        shp[8][0]= (1-xi*xi)*eta*(eta-1)*zeta*(zeta-1)/4.0;
        shp[8][1]=-xi*eta*(eta-1)*zeta*(zeta-1)/2.0;
        shp[8][2]= (1-xi*xi)*(2*eta-1)*zeta*(zeta-1)/4.0;
        shp[8][3]= (1-xi*xi)*eta*(eta-1)*(2*zeta-1)/4.0;

        shp[9][0]= xi*(xi+1)*(1-eta*eta)*zeta*(zeta-1)/4.0;
        shp[9][1]= (2*xi+1)*(1-eta*eta)*zeta*(zeta-1)/4.0;
        shp[9][2]=-xi*(xi+1)*eta*zeta*(zeta-1)/2.0;
        shp[9][3]= xi*(xi+1)*(1-eta*eta)*(2*zeta-1)/4.0;

        shp[10][0]= (1-xi*xi)*eta*(eta+1)*zeta*(zeta-1)/4.0;
        shp[10][1]=-xi*eta*(eta+1)*zeta*(zeta-1)/2.0;
        shp[10][2]= (1-xi*xi)*(2*eta+1)*zeta*(zeta-1)/4.0;
        shp[10][3]= (1-xi*xi)*eta*(eta+1)*(2*zeta-1)/4.0;

        shp[11][0]= xi*(xi-1)*(1-eta*eta)*zeta*(zeta-1)/4.0;
        shp[11][1]= (2*xi-1)*(1-eta*eta)*zeta*(zeta-1)/4.0;
        shp[11][2]=-xi*(xi-1)*eta*zeta*(zeta-1)/2.0;
        shp[11][3]= xi*(xi-1)*(1-eta*eta)*(2*zeta-1)/4.0;

        shp[12][0]= (1-xi*xi)*eta*(eta-1)*zeta*(zeta+1)/4.0;
        shp[12][1]=-xi*eta*(eta-1)*zeta*(zeta+1)/2.0;
        shp[12][2]= (1-xi*xi)*(2*eta-1)*zeta*(zeta+1)/4.0;
        shp[12][3]= (1-xi*xi)*eta*(eta-1)*(2*zeta+1)/4.0;

        shp[13][0]= xi*(xi+1)*(1-eta*eta)*zeta*(zeta+1)/4.0;
        shp[13][1]= (2*xi+1)*(1-eta*eta)*zeta*(zeta+1)/4.0;
        shp[13][2]=-xi*(xi+1)*eta*zeta*(zeta+1)/2.0;
        shp[13][3]= xi*(xi+1)*(1-eta*eta)*(2*zeta+1)/4.0;

        shp[14][0]= (1-xi*xi)*eta*(eta+1)*zeta*(zeta+1)/4.0;
        shp[14][1]=-xi*eta*(eta+1)*zeta*(zeta+1)/2.0;
        shp[14][2]= (1-xi*xi)*(2*eta+1)*zeta*(zeta+1)/4.0;
        shp[14][3]= (1-xi*xi)*eta*(eta+1)*(2*zeta+1)/4.0;

        shp[15][0]= xi*(xi-1)*(1-eta*eta)*zeta*(zeta+1)/4.0;
        shp[15][1]= (2*xi-1)*(1-eta*eta)*zeta*(zeta+1)/4.0;
        shp[15][2]=-xi*(xi-1)*eta*zeta*(zeta+1)/2.0;
        shp[15][3]= xi*(xi-1)*(1-eta*eta)*(2*zeta+1)/4.0;

        shp[16][0]= xi*(xi-1)*eta*(eta-1)*(1-zeta*zeta)/4.0;
        shp[16][1]= (2*xi-1)*eta*(eta-1)*(1-zeta*zeta)/4.0;
        shp[16][2]= xi*(xi-1)*(2*eta-1)*(1-zeta*zeta)/4.0;
        shp[16][3]=-xi*(xi-1)*eta*(eta-1)*zeta/2.0;

        shp[17][0]= xi*(xi+1)*eta*(eta-1)*(1-zeta*zeta)/4.0;
        shp[17][1]= (2*xi+1)*eta*(eta-1)*(1-zeta*zeta)/4.0;
        shp[17][2]= xi*(xi+1)*(2*eta-1)*(1-zeta*zeta)/4.0;
        shp[17][3]=-xi*(xi+1)*eta*(eta-1)*zeta/2.0;

        shp[18][0]= xi*(xi+1)*eta*(eta+1)*(1-zeta*zeta)/4.0;
        shp[18][1]= (2*xi+1)*eta*(eta+1)*(1-zeta*zeta)/4.0;
        shp[18][2]= xi*(xi+1)*(2*eta+1)*(1-zeta*zeta)/4.0;
        shp[18][3]=-xi*(xi+1)*eta*(eta+1)*zeta/2.0;

        shp[19][0]= xi*(xi-1)*eta*(eta+1)*(1-zeta*zeta)/4.0;
        shp[19][1]= (2*xi-1)*eta*(eta+1)*(1-zeta*zeta)/4.0;
        shp[19][2]= xi*(xi-1)*(2*eta+1)*(1-zeta*zeta)/4.0;
        shp[19][3]=-xi*(xi-1)*eta*(eta+1)*zeta/2.0;

        shp[20][0]= xi*(xi-1)*(1-eta*eta)*(1-zeta*zeta)/2.0;
        shp[20][1]= (2*xi-1)*(1-eta*eta)*(1-zeta*zeta)/2.0;
        shp[20][2]=-xi*(xi-1)*eta*(1-zeta*zeta);
        shp[20][3]=-xi*(xi-1)*(1-eta*eta)*zeta;

        shp[21][0]= xi*(xi+1)*(1-eta*eta)*(1-zeta*zeta)/2.0;
        shp[21][1]= (2*xi+1)*(1-eta*eta)*(1-zeta*zeta)/2.0;
        shp[21][2]=-xi*(xi+1)*eta*(1-zeta*zeta);
        shp[21][3]=-xi*(xi+1)*(1-eta*eta)*zeta;

        shp[22][0]= (1-xi*xi)*eta*(eta-1)*(1-zeta*zeta)/2.0;
        shp[22][1]=-xi*eta*(eta-1)*(1-zeta*zeta);
        shp[22][2]= (1-xi*xi)*(2*eta-1)*(1-zeta*zeta)/2.0;
        shp[22][3]=-(1-xi*xi)*eta*(eta-1)*zeta;

        shp[23][0]= (1-xi*xi)*eta*(eta+1)*(1-zeta*zeta)/2.0;
        shp[23][1]=-xi*eta*(eta+1)*(1-zeta*zeta);
        shp[23][2]= (1-xi*xi)*(2*eta+1)*(1-zeta*zeta)/2.0;
        shp[23][3]=-(1-xi*xi)*eta*(eta+1)*zeta;

        shp[24][0]= (1-xi*xi)*(1-eta*eta)*zeta*(zeta-1)/2.0;
        shp[24][1]=-xi*(1-eta*eta)*zeta*(zeta-1);
        shp[24][2]=-(1-xi*xi)*eta*zeta*(zeta-1);
        shp[24][3]= (1-xi*xi)*(1-eta*eta)*(2*zeta-1)/2.0;

        shp[25][0]= (1-xi*xi)*(1-eta*eta)*zeta*(zeta+1)/2.0;
        shp[25][1]=-xi*(1-eta*eta)*zeta*(zeta+1);
        shp[25][2]=-(1-xi*xi)*eta*zeta*(zeta+1);
        shp[25][3]= (1-xi*xi)*(1-eta*eta)*(2*zeta+1)/2.0;

        shp[26][0]= (1-xi*xi)*(1-eta*eta)*(1-zeta*zeta);
        shp[26][1]=-2*xi*(1-eta*eta)*(1-zeta*zeta);
        shp[26][2]=-2*(1-xi*xi)*eta*(1-zeta*zeta);
        shp[26][3]=-2*(1-xi*xi)*(1-eta*eta)*zeta;
    }

    dxdxi=0.0;dxdeta=0.0;dxdzeta=0.0;
    dydxi=0.0;dydeta=0.0;dydzeta=0.0;
    dzdxi=0.0;dzdeta=0.0;dzdzeta=0.0;
    for(i=0;i<nnodes;++i)
    {
        dxdxi  +=shp[i][1]*Coords[i][1];
        dydxi  +=shp[i][1]*Coords[i][2];
        dzdxi  +=shp[i][1]*Coords[i][3];

        dxdeta +=shp[i][2]*Coords[i][1];
        dydeta +=shp[i][2]*Coords[i][2];
        dzdeta +=shp[i][2]*Coords[i][3];

        dxdzeta+=shp[i][3]*Coords[i][1];
        dydzeta+=shp[i][3]*Coords[i][2];
        dzdzeta+=shp[i][3]*Coords[i][3];
    }

    Jac<<dxdxi,dydxi,dzdxi,
            dxdeta,dydeta,dzdeta,
            dxdzeta,dydzeta,dzdzeta;


    DetJac=Jac.determinant();

    if(fabs(DetJac)<1.e-13)
    {
        cout<<"*** Error: Singular element(happen in Shp3D.cpp)!!!"<<endl;
        abort();
    }

    if(!IsNature)
    {
        XJac=Jac.inverse();
        double temp1,temp2;
        for(i=0;i<nnodes;++i)
        {
            temp1   = XJac.coeff(0,0)*shp[i][1]
                      +XJac.coeff(0,1)*shp[i][2]
                      +XJac.coeff(0,2)*shp[i][3];

            temp2   = XJac.coeff(1,0)*shp[i][1]
                      +XJac.coeff(1,1)*shp[i][2]
                      +XJac.coeff(1,2)*shp[i][3];

            shp[i][3]= XJac.coeff(2,0)*shp[i][1]
                       +XJac.coeff(2,1)*shp[i][2]
                       +XJac.coeff(2,2)*shp[i][3];

            shp[i][1]=temp1;
            shp[i][2]=temp2;
        }
    }
    */
}
