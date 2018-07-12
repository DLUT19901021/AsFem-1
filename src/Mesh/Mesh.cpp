//*************************************************************
//* This file is part of the AsFem framework                  *
//* https://github.com/walkandthinker/AsFem                   *
//* All rights reserved, see COPYRIGHT for full restrictions  *
//* Licensed under LGPL 3.0, please see LICENSE for details   *
//* https://www.gnu.org/licenses/gpl-3.0.html                 *
//*************************************************************
//*** AsFem: A simple finite element method program         ***
//*** Author: Yang Bai @CopyRight                           ***
//*** Bug report: walkandthinker@gmail.com                  ***
//*** QQ group: 797998860                                   ***
//*************************************************************
//*** Created by Y. Bai on 12.06.18.                        ***
//*** Define the mesh class of AsFem                        ***
//*************************************************************

#include "Mesh/Mesh.h"

Mesh::Mesh()
{
    Xmin=0.0;Xmax=0.0;
    Ymin=0.0;Ymax=0.0;
    Zmin=0.0;Zmax=0.0;
    Conn.clear();
    NodeCoords.clear();
    nNodes=0;nElmts=0;nNodesPerElmt=0;

    IsDimSet=false;IsNxSet=false;IsNySet=false;IsNzSet=false;
    IsMeshTypeSet=false;

    IsXminSet=false;IsXmaxSet=false;
    IsYminSet=false;IsYmaxSet=false;
    IsZminSet=false;IsZmaxSet=false;

    IsMeshCreated=false;
    IsBCMeshSplit=false;

    VTKCellType=-1;
}

//*************************************************
// Set mesh information, before mesh generation ***
//*************************************************
void Mesh::SetDim(int dim)
{
    if(IsDimSet)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Warning: dim=%2d is already given!      ***\n",nDims);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem will do nothing for you !        ***\n");
    }
    else
    {
        if(dim<1||dim>3)
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: dim=%2d is an invalid value!!!   ***\n",dim);
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscFinalize();
            abort();
        }
        nDims=dim;
        IsDimSet=true;
    }
}

void Mesh::SetNx(int nx)
{
    if(IsNxSet)
    {

        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Warning: Nx=%2d is already given!       ***\n",Nx);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem will do nothing for you !        ***\n");
    }
    else
    {
        if(nx<0)
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: Nx=%2d is an invalid value!!!    ***\n",nx);
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscFinalize();
            abort();
        }
        Nx=nx;
        IsNxSet=true;
    }
}
void Mesh::SetNy(int ny)
{
    if(IsNySet)
    {

        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Warning: Ny=%2d is already given!       ***\n",Ny);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem will do nothing for you !        ***\n");
    }
    else
    {
        if(ny<0)
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: Ny=%2d is an invalid value!!!    ***\n",ny);
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscFinalize();
            abort();
        }
        Ny=ny;
        IsNySet=true;
    }
}
void Mesh::SetNz(int nz)
{
    if(IsNzSet)
    {

        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Warning: Nz=%2d is already given!       ***\n",Nz);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem will do nothing for you !        ***\n");
    }
    else
    {
        if(nz<0)
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: Nz=%2d is an invalid value!!!    ***\n",nz);
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscFinalize();
            abort();
        }
        Nz=nz;
        IsNzSet=true;
    }
}
//************************************************
void Mesh::SetXmin(double xmin)
{
    if(IsXmaxSet)
    {
        if(xmin>=Xmax)
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: xmin=%.3f is larger than xmax  ***\n",xmin);
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscFinalize();
            abort();
        }
        else
        {
            Xmin=xmin;
            IsXminSet=true;
        }
    }
    else
    {
        Xmin=xmin;
        IsXminSet=true;
    }
}
void Mesh::SetXmax(double xmax)
{
    if(IsXmaxSet)
    {
        if(xmax<=Xmin)
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: xmax=%.3f is smaller than xmin ***\n",xmax);
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscFinalize();
            abort();
        }
        else
        {
            Xmax=xmax;
            IsXmaxSet=true;
        }
    }
    else
    {
        Xmax=xmax;
        IsXmaxSet=true;
    }
}
//************************************************
void Mesh::SetYmin(double ymin)
{
    if(IsDimSet)
    {
        if(nDims<2)
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: dim=%2d ,Y-settings not allowed  ***\n",nDims);
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscFinalize();
            abort();
        }
    }
    if(IsYmaxSet)
    {
        if(ymin>=Ymax)
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: ymin=%.3f is larger than ymax  ***\n",ymin);
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscFinalize();
            abort();
        }
        else
        {
            Ymin=ymin;
            IsYminSet=true;
        }
    }
    else
    {
        Ymin=ymin;
        IsYminSet=true;
    }
}
void Mesh::SetYmax(double ymax)
{
    if(IsDimSet)
    {
        if(nDims<2)
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: dim=%2d ,Y-settings not allowed  ***\n",nDims);
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscFinalize();
            abort();
        }
    }

    if(IsYminSet)
    {
        if(ymax<=Ymin)
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: ymax=%.3f is smaller than ymin ***\n",ymax);
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscFinalize();
            abort();
        }
        else
        {
            Ymax=ymax;
            IsYmaxSet=true;
        }
    }
    else
    {
        Ymax=ymax;
        IsYmaxSet=true;
    }
}
//************************************************
void Mesh::SetZmin(double zmin)
{
    if(IsDimSet)
    {
        if(nDims<3)
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: dim=%2d ,Z-settings not allowed  ***\n",nDims);
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscFinalize();
            abort();
        }
    }

    if(IsZmaxSet)
    {
        if(zmin>=Zmax)
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: zmin=%.3f is larger than zmax  ***\n",zmin);
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscFinalize();
            abort();
        }
        else
        {
            Zmin=zmin;
            IsZminSet=true;
        }
    }
    else
    {
        Zmin=zmin;
        IsZminSet=true;
    }
}
void Mesh::SetZmax(double zmax)
{
    if(IsDimSet)
    {
        if(nDims<3)
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: dim=%2d ,Z-settings not allowed  ***\n",nDims);
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscFinalize();
            abort();
        }
    }

    if(IsZminSet)
    {
        if(zmax<=Zmin)
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: zmax=%.3f is smaller than zmin ***\n",zmax);
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscFinalize();
            abort();
        }
        else
        {
            Zmax=zmax;
            IsZmaxSet=true;
        }
    }
    else
    {
        Zmax=zmax;
        IsZmaxSet=true;
    }
}

//*************************************
void Mesh::SetMeshType(string meshtype)
{
    if(IsDimSet)
    {
        if(nDims==1)
        {
            if(meshtype.compare("edge2")==0 ||
               meshtype.compare("edge3")==0 ||
               meshtype.compare("edge4")==0)
            {
                MeshType=meshtype;
                IsMeshTypeSet=true;
            }
            else
            {
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: unsupported mesh type !!!       ***\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        meshtype=edge2,3,4 is expected! ***\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                PetscFinalize();
                abort();
            }
        }
        else if(nDims==2)
        {
            if(meshtype.compare("quad4")==0 ||
               meshtype.compare("quad8")==0 ||
               meshtype.compare("quad9")==0)
            {
                MeshType=meshtype;
                IsMeshTypeSet=true;
            }
            else
            {
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: unsupported mesh type!!!        ***\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        meshtype=quad4,8,9 is expected! ***\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                PetscFinalize();
                abort();
            }
        }
        else if(nDims==3)
        {
            if(meshtype.compare("hex8")==0 ||
               meshtype.compare("hex20")==0 ||
               meshtype.compare("hex27")==0)
            {
                MeshType=meshtype;
                IsMeshTypeSet=true;
            }
            else
            {
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: unsupported mesh type!!!        ***\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        meshtype=hex8,20,27 is expected!***\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                PetscFinalize();
                abort();
            }
        }
    }
    else
    {
        if(meshtype.compare("edge2")==0 ||
           meshtype.compare("edge3")==0 ||
           meshtype.compare("edge4")==0)
        {
            MeshType=meshtype;
            IsMeshTypeSet=true;
        }
        else if(meshtype.compare("quad4")==0 ||
                meshtype.compare("quad8")==0 ||
                meshtype.compare("quad9")==0)
        {
            MeshType=meshtype;
            IsMeshTypeSet=true;
        }
        else if(meshtype.compare("hex8")==0 ||
                meshtype.compare("hex20")==0 ||
                meshtype.compare("hex27")==0)
        {
            MeshType=meshtype;
            IsMeshTypeSet=true;
        }
        else
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: unsupported mesh type !!!       ***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** For 1D:meshtype=edge2,3,4 is expected! ***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** For 2D:meshtype=quad4,8,9 is expected! ***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** For 3D:meshtype=hex8,20,27 is expected!***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscFinalize();
            abort();

        }
    }
}

//**********************************
bool Mesh::IsRequiredInfoComplete()
{
    if(IsDimSet)
    {
        if(nDims==1)
        {
            if(IsMeshTypeSet && IsXminSet && IsXmaxSet && IsNySet)
            {
                return true;
            }
            else
            {
                return false;
            }
        }
        else if(nDims==2)
        {
            if(IsMeshTypeSet&&
               IsXminSet && IsXmaxSet &&
               IsYminSet && IsYmaxSet &&
               IsNxSet && IsNySet)
            {
                return true;
            }
            else
            {
                return false;
            }
        }
        else if(nDims==3)
        {
            if(IsMeshTypeSet &&
               IsXminSet && IsXmaxSet&&
               IsYminSet && IsYmaxSet&&
               IsZminSet && IsZmaxSet&&
               IsNxSet && IsNySet && IsNzSet)
            {
                return true;
            }
            else
            {
                return false;
            }
        }
        return false;
    }
    else
    {
        return false;
    }
}