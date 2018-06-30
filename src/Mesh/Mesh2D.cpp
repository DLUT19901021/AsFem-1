//*****************************************************
//*** AsFem: A simple finite element method program ***
//*** Author: Yang Bai                              ***
//*** Bug report: walkandthinker@gmail.com          ***
//*** QQ group: 797998860                           ***
//*****************************************************
//*** Define 2D mesh                                ***
//***  inherit from mesh base                       ***
//*****************************************************

#include "Mesh/Mesh2D.h"

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

    Nx=nx;Ny=ny;nDims=2;
    nNodes=0;nNodesPerElmt=0;nElmts=0;
    MeshGenerated=false;
    ElmtType=elmttype;

    BoundaryElmtSet.clear();

    VTKCellType=9;

    MeshGenerated=false;
}

void Mesh2D::Release()
{
    if(MeshGenerated)
    {
        NodeCoords.clear();
        Conn.clear();

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

        
        NodeCoords.resize(nNodes*4,0.0);
        Conn.resize(nElmts*nNodesPerElmt);



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


        NodeCoords.resize(nNodes*4,0.0);
        Conn.resize(nElmts*nNodesPerElmt,0);



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
            }
            // for middle line of each element
            for(i=1;i<=Nx+1;i++)
            {
                k=(j-1)*(2*Nx+1+Nx+1)+2*Nx+1+i;


                NodeCoords[(k-1)*4  ]=1.0;
                NodeCoords[(k-1)*4+1]=Xmin+(i-1)*2*dx;
                NodeCoords[(k-1)*4+2]=Ymin+(j-1)*2*dy+dy;
                NodeCoords[(k-1)*4+3]=0.0;

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


        NodeCoords.resize(nNodes*4,0.0);
        Conn.resize(nElmts*nNodesPerElmt,0);



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

    SplitBoundaryMesh();

    return MeshGenerated;
}

// split the boundary mesh into the BoundaryElmtSet
void Mesh2D::SplitBoundaryMesh()
{
    int i,j,k,e;
    if(!MeshGenerated)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't split boundary mesh, mesh hasn't been generated!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        CreateMesh should be called before this!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }

    vector<int> left,right,bottom,top;


    if(ElmtType=="quad4")
    {
        nNodesPerBCElmt=2;
        nBCElmts=Ny*2+Nx*2;
        nBCNodes=0;
        // for left side element
        i=1;left.clear();
        for(j=1;j<=Ny;j++)
        {
            e=(j-1)*Nx+i;
            left.push_back(IthConnJthIndex(e,4));
            left.push_back(IthConnJthIndex(e,1));
            nBCNodes+=2;
        }
        // for right side element
        i=Nx;right.clear();
        for(j=1;j<=Ny;j++)
        {
            e=(j-1)*Nx+i;
            right.push_back(IthConnJthIndex(e,2));
            right.push_back(IthConnJthIndex(e,3));
            nBCNodes+=2;
        }
        // for bottom element
        bottom.clear();
        j=1;
        for(i=1;i<=Nx;i++)
        {
            e=(j-1)*Nx+i;
            bottom.push_back(IthConnJthIndex(e,1));
            bottom.push_back(IthConnJthIndex(e,2));
            nBCNodes+=2;
        }
        // for bottom element
        top.clear();
        j=Ny;
        for(i=1;i<=Nx;i++)
        {
            e=(j-1)*Nx+i;
            top.push_back(IthConnJthIndex(e,3));
            top.push_back(IthConnJthIndex(e,4));
            nBCNodes+=2;
        }
    }
    else if(ElmtType=="quad8"||ElmtType=="quad9")
    {
        nNodesPerBCElmt=3;
        nBCElmts=Nx*2+Ny*2;
        nBCNodes=0;
        // 4***7***3
        // |       *
        // 8       6
        // |       |
        // 1***5***2
        // for left side set
        left.clear();
        i=1;
        for(j=1;j<=Ny;j++)
        {
            e=(j-1)*Nx+i;
            left.push_back(IthConnJthIndex(e,4));
            left.push_back(IthConnJthIndex(e,8));
            left.push_back(IthConnJthIndex(e,1));
            nBCNodes+=3;
        }
        // for right side set
        right.clear();
        i=Nx;
        for(j=1;j<=Ny;j++)
        {
            e=(j-1)*Nx+i;
            right.push_back(IthConnJthIndex(e,2));
            right.push_back(IthConnJthIndex(e,6));
            right.push_back(IthConnJthIndex(e,3));
            nBCNodes+=3;
        }
        // for bottom side set
        bottom.clear();
        j=1;
        for(i=1;i<=Nx;i++)
        {
            e=(j-1)*Nx+i;
            bottom.push_back(IthConnJthIndex(e,1));
            bottom.push_back(IthConnJthIndex(e,5));
            bottom.push_back(IthConnJthIndex(e,2));
            nBCNodes+=3;
        }
        // for top side set
        top.clear();
        j=Ny;
        for(i=1;i<=Nx;i++)
        {
            e=(j-1)*Nx+i;
            top.push_back(IthConnJthIndex(e,3));
            top.push_back(IthConnJthIndex(e,7));
            top.push_back(IthConnJthIndex(e,4));
            nBCNodes+=3;
        }
    }

    LeftBCElmtSet=make_pair("left",left);
    RightBCElmtSet=make_pair("right",right);

    BottomBCElmtSet=make_pair("bottom",bottom);
    TopBCElmtSet=make_pair("top",top);

    BoundaryElmtSet.push_back(LeftBCElmtSet);
    BoundaryElmtSet.push_back(RightBCElmtSet);
    BoundaryElmtSet.push_back(BottomBCElmtSet);
    BoundaryElmtSet.push_back(TopBCElmtSet);
}

pair<string,vector<int>> Mesh2D::GetSideSet(string sidename) const
{
    if(sidename=="left")
    {
        return LeftBCElmtSet;
    }
    else if(sidename=="right")
    {
        return RightBCElmtSet;
    }
    else if(sidename=="bottom")
    {
        return BottomBCElmtSet;
    }
    else if(sidename=="top")
    {
        return TopBCElmtSet;
    }
    else
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: unsupported sidename(=%s)!\n",sidename.c_str());
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        please use left,right,top,bottom,back,front!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }
};
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
int Mesh2D::GetSideNodesNum(string sidename) const
{
    if(sidename=="left")
    {
        return LeftBCElmtSet.second.size();
    }
    else if(sidename=="right")
    {
        return RightBCElmtSet.second.size();
    }
    else if(sidename=="bottom")
    {
        return BottomBCElmtSet.second.size();
    }
    else if(sidename=="top")
    {
        return TopBCElmtSet.second.size();
    }
    else
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: unsupported sidename(=%s)!\n",sidename.c_str());
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        please use left,right,top,bottom,back,front!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }
}

int Mesh2D::GetSideElmtsNum(string sidename) const
{
    if(sidename=="left")
    {
        return LeftBCElmtSet.second.size()/nNodesPerBCElmt;
    }
    else if(sidename=="right")
    {
        return RightBCElmtSet.second.size()/nNodesPerBCElmt;
    }
    else if(sidename=="bottom")
    {
        return BottomBCElmtSet.second.size()/nNodesPerBCElmt;
    }
    else if(sidename=="top")
    {
        return TopBCElmtSet.second.size()/nNodesPerBCElmt;
    }
    else
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: unsupported sidename(=%s)!\n",sidename.c_str());
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        please use left,right,top,bottom,back,front!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }
}

vector<int> Mesh2D::GetIthBCElmtConn(string sidename,int e) const
{
    vector<int> conn;
    conn.resize(nNodesPerBCElmt);
    if(sidename=="left")
    {
        if(e<1||e>LeftBCElmtSet.second.size()/nNodesPerBCElmt)
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: e=%d is out of boundary element's range!\n",e);
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
            PetscFinalize();
            abort();
        }
        for(int i=1;i<=nNodesPerBCElmt;i++)
        {
            conn[i-1]=LeftBCElmtSet.second[(e-1)*nNodesPerBCElmt+i-1];
        }
    }
    else if(sidename=="right")
    {
        if(e<1||e>RightBCElmtSet.second.size()/nNodesPerBCElmt)
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: e=%d is out of boundary element's range!\n",e);
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
            PetscFinalize();
            abort();
        }
        for(int i=1;i<=nNodesPerBCElmt;i++)
        {
            conn[i-1]=RightBCElmtSet.second[(e-1)*nNodesPerBCElmt+i-1];
        }
    }
    else if(sidename=="bottom")
    {
        if(e<1||e>BottomBCElmtSet.second.size()/nNodesPerBCElmt)
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: e=%d is out of boundary element's range!\n",e);
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
            PetscFinalize();
            abort();
        }
        for(int i=1;i<=nNodesPerBCElmt;i++)
        {
            conn[i-1]=BottomBCElmtSet.second[(e-1)*nNodesPerBCElmt+i-1];
        }
    }
    else if(sidename=="top")
    {
        if(e<1||e>TopBCElmtSet.second.size()/nNodesPerBCElmt)
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: e=%d is out of boundary element's range!\n",e);
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
            PetscFinalize();
            abort();
        }
        for(int i=1;i<=nNodesPerBCElmt;i++)
        {
            conn[i-1]=TopBCElmtSet.second[(e-1)*nNodesPerBCElmt+i-1];
        }
    }
    else
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: unsupported sidename(=%s)!\n",sidename.c_str());
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        please use left,right,top,bottom,back,front!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }

    return conn;
}

//**************************************************

//**************************************************
void Mesh2D::PrintMeshInfo(string str) const
{
    int i,j;
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

    i=LeftBCElmtSet.second.size();j=LeftBCElmtSet.second.size()/nNodesPerBCElmt;
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   Left side nodes=%5d, elmts=%5d   ***\n",i,j);

    i=RightBCElmtSet.second.size();j=RightBCElmtSet.second.size()/nNodesPerBCElmt;
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   Right side nodes=%5d, elmts=%5d  ***\n",i,j);

    i=BottomBCElmtSet.second.size();j=BottomBCElmtSet.second.size()/nNodesPerBCElmt;
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   Bottom side nodes=%5d, elmts=%5d ***\n",i,j);

    i=TopBCElmtSet.second.size();j=TopBCElmtSet.second.size()/nNodesPerBCElmt;
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   Top side nodes=%5d, elmts=%5d    ***\n",i,j);

    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");

}

void Mesh2D::PrintMeshDetailInfo(string str) const
{
    int i,j,e;
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


    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Boundary element's connectivity:       ***\n");
    string sidename;
    for(i=1;i<=BoundaryElmtSet.size();++i)
    {
        sidename=BoundaryElmtSet[i-1].first;

        for(e=1;e<=BoundaryElmtSet[i-1].second.size()/nNodesPerBCElmt;++e)
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** %s side->",sidename.c_str());
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"e=%3d:",e);
            for(j=1;j<=nNodesPerBCElmt;j++)
            {
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"%5d ",BoundaryElmtSet[i-1].second[(e-1)*nNodesPerBCElmt+j-1]);
            }
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\n");
        }
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\n");
    }
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");

}
