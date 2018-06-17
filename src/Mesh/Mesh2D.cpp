//*****************************************************
//*** AsFem: A simple finite element method program ***
//*** Author: Yang Bai                              ***
//*** Bug report: walkandthinker@gmail.com          ***
//*** QQ group: 797998860                           ***
//*****************************************************
//*** Define 2D mesh                                ***
//***  inherit from mesh base                       ***
//*****************************************************

#include "Mesh2D.h"

Mesh2D::Mesh2D(double xmin, double xmax, double ymin, double ymax, int nx, int ny, string elmttype)
{
    if(xmin>=xmax)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: xmin=%12.6e is larger than xmax=12.6e!!!\n",xmin,xmax);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }

    if(ymin>=ymax)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: ymin=%12.6e is larger than ymax=12.6e!!!\n",ymin,ymax);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }

    if(nx<1)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: nx=%6d is too less for a FEM problem!!!\n",nx);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }
    if(ny<1)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: ny=%6d is too less for a FEM problem!!!\n",ny);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }

    string str=RemoveSpace(elmttype);
    elmttype=StrToLower(str);
    if(elmttype.find("quad4")==string::npos&&
       elmttype.find("quad8")==string::npos&&
       elmttype.find("quad9")==string::npos)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: unsupported mesh type !!!\n",nx);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        only type=edge2, edge3 or edge4 is supported for 1D mesh !!!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }

    Xmin=xmin;Xmax=xmax;
    Ymin=ymin;Ymax=ymax;

    Nx=nx;Ny=ny;
    nNodes=0;nNodesPerElmt=0;nElmts=0;
    MeshGenerated=false;
    ElmtType=elmttype;

    nBoundaryNodeIndex=0;
    nLeftEdgeNodeIndex=0;nRightEdgeNodeIndex=0;
    nBoundaryNodeIndex=0;nTopEdgeNodeIndex=0;

    VTKCellType=9;

    MeshGenerated=false;
}

void Mesh2D::Release()
{
    if(MeshGenerated)
    {
        delete[] NodeCoords;
        delete[] Conn;
        delete[] BoundaryNodeIndex;

        /*
         * TODO: assign node index or split mesh into boundary mesh
         *       no idea currently, but I should find solution
        delete[] LeftEdgeNodeIndex;
        delete[] RightEdgeNodeIndex;

        delete[] BottomEdgeNodeIndex;
        delete[] TopEdgeNodeIndex;
         */
    }
}

//************************************************
bool Mesh2D::CreateMesh()
{
    double dx,dy;
    int e,i,j,k,kk;
    int i1,i2,i3,i4,i5,i6,i7,i8,i9;

    if(ElmtType=="quad4")
    {
        dx=(Xmax-Xmin)/Nx;
        dy=(Ymax-Ymin)/Ny;


        nElmts=Nx*Ny;
        nNodes=(Nx+1)*(Ny+1);
        nNodesPerElmt=4;
        P=1;

        NodeCoords=new double[nNodes*4];
        Conn=new int[nElmts*nNodesPerElmt];

        if(Nx==1||Ny==1)
        {
            BoundaryNodeIndex=new int[nNodes];
            nBoundaryNodeIndex=nNodes;
        }
        else
        {
            i=2*Nx+2*Ny;
            BoundaryNodeIndex=new int[i];
            nBoundaryNodeIndex=i;
        }

        kk=0;
        for(j=1;j<=Ny+1;++j)
        {
            for(i=1;i<=Nx+1;++i)
            {
                k=(j-1)*(Nx+1)+i;


                NodeCoords[(k-1)*4  ]=1.0;
                NodeCoords[(k-1)*4+1]=Xmin+(i-1)*dx;
                NodeCoords[(k-1)*4+2]=Ymin+(j-1)*dy;
                NodeCoords[(k-1)*4+3]=0.0;

                if(j==1||j==Ny+1)
                {
                    kk+=1;
                    BoundaryNodeIndex[kk-1]=k;
                }
                else
                {
                    if(i==1||i==Nx+1)
                    {
                        kk+=1;
                        BoundaryNodeIndex[kk-1]=k;
                    }
                }
            }
        }

        // Create Connectivity matrix
        kk=0;
        for(j=1;j<=Ny;j++)
        {
            for(i=1;i<=Nx;i++)
            {
                e=(j-1)*Nx+i;
                i1=(j-1)*(Nx+1)+i;
                i2=i1+1;
                i3=i2+Nx+1;
                i4=i3-1;


                Conn[(e-1)*nNodesPerElmt+1-1]=i1;
                Conn[(e-1)*nNodesPerElmt+2-1]=i2;
                Conn[(e-1)*nNodesPerElmt+3-1]=i3;
                Conn[(e-1)*nNodesPerElmt+4-1]=i4;

            }
        }

        VTKCellType=9;
        MeshGenerated=true;
    }
    else if(ElmtType=="quad8")
    {
        // for 2D-8 nodes mesh
        dx=(Xmax-Xmin)/(2.0*Nx);
        dy=(Ymax-Ymin)/(2.0*Ny);

        nElmts=Nx*Ny;
        nNodes=(2*Nx+1)*(2*Ny+1)-nElmts;
        nNodesPerElmt=8;
        P=2;


        NodeCoords=new double[nNodes*4];
        Conn=new int[nElmts*nNodesPerElmt];

        if(Nx==1||Ny==1)
        {
            BoundaryNodeIndex=new int[nNodes-(nElmts-1)];
            nBoundaryNodeIndex=nNodes-(nElmts-1);

        }
        else
        {
            i=2*(2*Nx+1)+2*(2*Ny+1-2);

            BoundaryNodeIndex=new int[i];
            nBoundaryNodeIndex=i;
        }


        kk=0;
        for(j=1;j<=Ny;j++)
        {
            // for bottom line of each element
            for(i=1;i<=2*Nx+1;i++)
            {
                k=(j-1)*(2*Nx+1+Nx+1)+i;


                NodeCoords[(k-1)*4  ]=1.0;
                NodeCoords[(k-1)*4+1]=Xmin+(i-1)*dx;
                NodeCoords[(k-1)*4+2]=Ymin+(j-1)*2*dy;
                NodeCoords[(k-1)*4+3]=0.0;
                if(j==1)
                {
                    kk+=1;
                    BoundaryNodeIndex[kk-1]=k;
                }
                else
                {
                    if(i==1||i==2*Nx+1)
                    {
                        kk+=1;
                        BoundaryNodeIndex[kk-1]=k;
                    }
                }
            }
            // for middle line of each element
            for(i=1;i<=Nx+1;i++)
            {
                k=(j-1)*(2*Nx+1+Nx+1)+2*Nx+1+i;


                NodeCoords[(k-1)*4  ]=1.0;
                NodeCoords[(k-1)*4+1]=Xmin+(i-1)*2*dx;
                NodeCoords[(k-1)*4+2]=Ymin+(j-1)*2*dy+dy;
                NodeCoords[(k-1)*4+3]=0.0;

                if(i==1||i==Nx+1)
                {
                    kk+=1;
                    BoundaryNodeIndex[kk-1]=k;
                }
            }
        }
        // for the last top line
        j=Ny+1;
        for(i=1;i<=2*Nx+1;i++)
        {
            k=(j-1)*(2*Nx+1+Nx+1)+i;


            NodeCoords[(k-1)*4  ]=1.0;
            NodeCoords[(k-1)*4+1]=Xmin+(i-1)*dx;
            NodeCoords[(k-1)*4+2]=Ymin+(j-1)*2*dy;
            NodeCoords[(k-1)*4+3]=0.0;

            kk+=1;
            BoundaryNodeIndex[kk-1]=k;
        }


        // Create Connectivity matrix
        kk=0;
        for(j=1;j<=Ny;j++)
        {
            for(i=1;i<=Nx;i++)
            {
                e=(j-1)*Nx+i;
                i1=(j-1)*(2*Nx+1+Nx+1)+2*i-1;
                i2=i1+2;
                i3=i2+(2*Nx+1+Nx+1);
                i4=i3-2;

                i5=i1+1;
                i6=i2+(2*Nx+1)-i;
                i7=i3-1;
                i8=i1+(2*Nx+1)-(i-1);


                Conn[(e-1)*nNodesPerElmt+1-1]=i1;
                Conn[(e-1)*nNodesPerElmt+2-1]=i2;
                Conn[(e-1)*nNodesPerElmt+3-1]=i3;
                Conn[(e-1)*nNodesPerElmt+4-1]=i4;
                Conn[(e-1)*nNodesPerElmt+5-1]=i5;
                Conn[(e-1)*nNodesPerElmt+6-1]=i6;
                Conn[(e-1)*nNodesPerElmt+7-1]=i7;
                Conn[(e-1)*nNodesPerElmt+8-1]=i8;
            }
        }

        VTKCellType=23;
        MeshGenerated=true;
    }
    else if(ElmtType=="quad9")
    {
        dx=(Xmax-Xmin)/(2.0*Nx);
        dy=(Ymax-Ymin)/(2.0*Ny);

        nElmts=Nx*Ny;
        nNodes=(2*Nx+1)*(2*Ny+1);
        nNodesPerElmt=9;
        P=2;


        NodeCoords=new double[nNodes*4];
        Conn=new int[nElmts*nNodesPerElmt];


        if(Nx==1||Ny==1)
        {
            BoundaryNodeIndex=new int[nNodes-2*nElmts+1];
            nBoundaryNodeIndex=nNodes-2*nElmts+1;
        }
        else
        {
            i=2*(2*Nx+1)+2*(2*Ny+1-2);
            BoundaryNodeIndex=new int[i];
            nBoundaryNodeIndex=i;
        }


        kk=0;
        for(j=1;j<=2*Ny+1;j++)
        {
            for(i=1;i<=2*Nx+1;i++)
            {
                k=(j-1)*(2*Nx+1)+i;


                NodeCoords[(k-1)*4  ]=1.0;
                NodeCoords[(k-1)*4+1]=Xmin+(i-1)*dx;
                NodeCoords[(k-1)*4+2]=Ymin+(j-1)*dy;
                NodeCoords[(k-1)*4+3]=0.0;

                if(j==1||j==2*Ny+1)
                {
                    kk+=1;
                    BoundaryNodeIndex[kk-1]=k;
                }
                else
                {
                    if(i==1||i==2*Nx+1)
                    {
                        kk+=1;
                        BoundaryNodeIndex[kk-1]=k;
                    }
                }
            }
        }


        // Create Connectivity matrix
        kk=0;
        for(j=1;j<=Ny;j++)
        {
            for(i=1;i<=Nx;i++)
            {
                e=(j-1)*Nx+i;
                i1=(j-1)*2*(2*Nx+1)+2*i-1;
                i2=i1+2;
                i3=i2+2*(2*Nx+1);
                i4=i3-2;

                i5=i1+1;
                i6=i2+(2*Nx+1);
                i7=i3-1;
                i8=i1+(2*Nx+1);
                i9=i8+1;


                Conn[(e-1)*nNodesPerElmt+1-1]=i1;
                Conn[(e-1)*nNodesPerElmt+2-1]=i2;
                Conn[(e-1)*nNodesPerElmt+3-1]=i3;
                Conn[(e-1)*nNodesPerElmt+4-1]=i4;
                Conn[(e-1)*nNodesPerElmt+5-1]=i5;
                Conn[(e-1)*nNodesPerElmt+6-1]=i6;
                Conn[(e-1)*nNodesPerElmt+7-1]=i7;
                Conn[(e-1)*nNodesPerElmt+8-1]=i8;
                Conn[(e-1)*nNodesPerElmt+9-1]=i9;

            }
        }

        VTKCellType=28;
        MeshGenerated=true;
    }

    return MeshGenerated;
}

//***************************************
double Mesh2D::IthNodeJthCoords(int i, int j) const
{
    if(!MeshGenerated)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't get node's coordinates, mesh hasn't been generated!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        CreateMesh should be called before this!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }

    if(i<1||i>nNodes)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: i=%6d is out of nNodes=%6d\n",i,nNodes);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }
    if(j<0||j>3)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: j=%2d is out 0~3\n",i,nNodes);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        j should be 0->w,1->x,2->y,3->z!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }

    return NodeCoords[4*(i-1)+j];
}

int Mesh2D::IthConnJthIndex(int e, int j) const
{
    if(!MeshGenerated)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't get connectivity info, mesh hasn't been generated!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        CreateMesh should be called before this!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }

    if(e<1||e>nElmts)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: e=%8d is out of nElmts=%8d\n",e,nNodes);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }
    if(j<1||j>nNodesPerElmt)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: j=%2d is out of nNodesPerElmt=%2d\n",j,nNodesPerElmt);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        j should be 1~%2d!\n",nNodesPerElmt);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }

    return Conn[(e-1)*nNodesPerElmt+j-1];
}

//**************************************************
void Mesh2D::PrintMeshInfo(string str) const
{
    if(!MeshGenerated)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't print mesh info, mesh hasn't been generated!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        CreateMesh should be called before this!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }


    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
    if(str.size()>1)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** str= %s\n",str.c_str());
    }
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** 1D mesh info:                          ***\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   element order=%2d                     ***\n",P);
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   nNodes=%10d                    ***\n",nNodes);
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   nElmts=%10D                    ***\n",nElmts);
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   nNodesPerElmt=%3d                    ***\n",nNodesPerElmt);
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");

}

void Mesh2D::PrintMeshDetailInfo(string str) const
{
    int i,e;
    double x,y,z,w;
    if(!MeshGenerated)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't print mesh info, mesh hasn't been generated!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        CreateMesh should be called before this!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }


    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
    if(str.size()>1)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** str= %s\n",str.c_str());
    }
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** 1D mesh detail info:                   ***\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   element order=%2d                     ***\n",P);
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   nNodes=%10d                    ***\n",nNodes);
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   nElmts=%10D                    ***\n",nElmts);
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   nNodesPerElmt=%3d                    ***\n",nNodesPerElmt);
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Nodes' coordinates:                    ***\n");
    for(i=1;i<=nNodes;++i)
    {
        w=IthNodeJthCoords(i,0);
        x=IthNodeJthCoords(i,1);
        y=IthNodeJthCoords(i,2);
        z=IthNodeJthCoords(i,3);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** %6d-th node:",i);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"x=%12.5e,y=%12.5e,z=%12.5e,w=%6.2f\n",x,y,z,w);
    }

    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Element's connectivity:                ***\n");
    for(e=1;e<=nElmts;++e)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** e=%6d:",e);
        for(i=1;i<=nNodesPerElmt;++i)
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"%5d ",IthConnJthIndex(e,i));
        }
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\n");
    }
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");

}
