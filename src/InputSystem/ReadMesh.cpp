//*****************************************************
//*** AsFem: A simple finite element method program ***
//*** Author: Yang Bai @CopyRight                   ***
//*** Bug report: walkandthinker@gmail.com          ***
//*** QQ group: 797998860                           ***
//*****************************************************
//*** We read mesh from input file                  ***
//*****************************************************

#include "InputSystem/InputSystem.h"
#include "StringUtils/StringUtils.h"

#include "petsc.h"

bool InputSystem::ReadMesh(Mesh &mesh)
{
    int linenum=0;
    string str,line,line0="mesh";
    int blockstartlinenum;
    vector<double> numbers;
    int dim;
    double xmin,xmax,ymin,ymax,zmin,zmax;
    string meshtype;
    bool ReadMeshSuccess=false;
    int nx,ny,nz;

    // Read the first comment line
    getline(in,line);linenum+=1;
    if(IsBracketMatch(in,line0,blockstartlinenum))
    {

        GotoLine(in,blockstartlinenum);linenum=blockstartlinenum-1;
        getline(in,line);linenum+=1;// read [mesh]


        // read type
        getline(in,line);linenum+=1;
        str=RemoveSpace(line);
        if(str.compare(0,10,"type=asfem")==0)
        {
            // read built-in mesh
            getline(in,line);linenum+=1;// read dim
            str=RemoveSpace(line);
            if(str.compare(0,4,"dim=")!=0)
            {
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find 'dim=' line-%-3d      ***\n",linenum);
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given dim=1,2,3 !!!  ***\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                return false;
            }
            numbers=SplitNum(line);
            if(numbers.size()<1)
            {
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find number in line-%-3d  ***\n",linenum);
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given dim=1,2,3 !!!  ***\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                return false;
            }
            else
            {
                dim=int(numbers[0]);
                if(dim<1||dim>3)
                {
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: dim=%2d is invalid in line-%-3d!!!***\n",dim,linenum);
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given dim=1,2,3 !!!  ***\n");
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                    return false;
                }
                else if(dim==1)
                {
                    getline(in,line);linenum+=1;// read xmin
                    str=RemoveSpace(line);
                    if(str.compare(0,5,"xmin=")==0)
                    {
                        numbers=SplitNum(line);
                        if(numbers.size()<1)
                        {
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find number in line-%-3d  ***\n",linenum);
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given dim=1,2,3 !!!  ***\n");
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                            return false;
                        }
                        xmin=numbers[0];
                    }
                    else
                    {
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find xmin=  in line-%3d***\n",linenum);
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given xmin=value !!! ***\n");
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                        return false;
                    }

                    //***for xmax
                    getline(in,line);linenum+=1;// read xmax
                    str=RemoveSpace(line);
                    if(str.compare(0,5,"xmax=")==0)
                    {
                        numbers=SplitNum(line);
                        if(numbers.size()<1)
                        {
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find number in line-%-3d  ***\n",linenum);
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given dim=1,2,3 !!!  ***\n");
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                            return false;
                        }
                        xmax=numbers[0];
                    }
                    else
                    {
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find xmax=  in line-%3d***\n",linenum);
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given xmax=value !!! ***\n");
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                        return false;
                    }

                    if(xmin>xmax)
                    {
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: invalid mesh information        ***\n",linenum);
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        xmin=%8.3f>xmax=%8.3f !!! ***\n",xmin,xmax);
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        please check your input file    ***\n");
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                        return false;
                    }

                    //***for nx
                    getline(in,line);linenum+=1;// read nx
                    str=RemoveSpace(line);
                    if(str.compare(0,3,"nx=")==0)
                    {
                        numbers=SplitNum(line);
                        if(numbers.size()<1)
                        {
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find number in line-%-3d  ***\n",linenum);
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given nx=int value   ***\n");
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                            return false;
                        }
                        nx=numbers[0];
                        if(nx<1)
                        {
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: nx=%5d(line-%3d) is invalid !!***\n",nx,linenum);
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given nx=value(>0)   ***\n");
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                            return false;
                        }
                    }
                    else
                    {
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: 'nx=' can't found in line=%-3d   ***\n",linenum);
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given nx=value(>0)   ***\n");
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                        return false;
                    }
                    //**** Read mesh type
                    getline(in,line);linenum+=1;//
                    str=RemoveSpace(line);
                    if(str.compare(0,9,"meshtype=")==0)
                    {
                        if(str.find("edge2")!=string::npos)
                        {
                            meshtype="edge2";
                        }
                        else if(str.find("edge3")!=string::npos)
                        {
                            meshtype="edge3";
                        }
                        else if(str.find("edge4")!=string::npos)
                        {
                            meshtype="edge4";
                        }
                        else
                        {
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: unsupported mesh type(line-%-3d) ***\n",linenum);
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        meshtype=edge2,3,4 is expected!!***\n");
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                            return false;
                        }
                    }

                    Mesh1D mesh1D(xmin,xmax,nx,meshtype);
                    mesh1D.CreateMesh();
                    mesh.Add1DMesh(mesh1D);
                    mesh.Init();
                    mesh.SetDims(dim);
                    return true;
                }
                else if(dim==2)
                {
                    getline(in,line);linenum+=1;// read xmin
                    str=RemoveSpace(line);
                    if(str.compare(0,5,"xmin=")==0)
                    {
                        numbers=SplitNum(line);
                        if(numbers.size()<1)
                        {
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find number in line-%-3d  ***\n",linenum);
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given dim=1,2,3 !!!  ***\n");
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                            return false;
                        }
                        xmin=numbers[0];
                    }
                    else
                    {
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find 'xmin='  in line-%-3d ***\n",linenum);
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given xmin=value !!! ***\n");
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                        return false;
                    }

                    //***for xmax
                    getline(in,line);linenum+=1;// read xmax
                    str=RemoveSpace(line);
                    if(str.compare(0,5,"xmax=")==0)
                    {
                        numbers=SplitNum(line);
                        if(numbers.size()<1)
                        {
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find number in line-%-3d  ***\n",linenum);
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given dim=1,2,3 !!!  ***\n");
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                            return false;
                        }
                        xmax=numbers[0];
                    }
                    else
                    {
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find 'xmax='  in line-%-3d ***\n",linenum);
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given xmax=value !!! ***\n");
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                        return false;
                    }

                    if(xmin>xmax)
                    {
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: invalid mesh information        ***\n",linenum);
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        xmin=%8.3f>xmax=%8.3f !!! ***\n",xmin,xmax);
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        please check your input file    ***\n");
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                        return false;
                    }
                    //****************************
                    getline(in,line);linenum+=1;// read ymin
                    str=RemoveSpace(line);
                    if(str.compare(0,5,"ymin=")==0)
                    {
                        numbers=SplitNum(line);
                        if(numbers.size()<1)
                        {
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find number in line-%-3d  ***\n",linenum);
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given dim=1,2,3 !!!  ***\n");
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                            return false;
                        }
                        ymin=numbers[0];
                    }
                    else
                    {
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find 'ymin=' in line-%3d  ***\n",linenum);
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given ymin=value !!! ***\n");
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                        return false;
                    }

                    //***for ymax
                    getline(in,line);linenum+=1;// read ymax
                    str=RemoveSpace(line);
                    if(str.compare(0,5,"ymax=")==0)
                    {
                        numbers=SplitNum(line);
                        if(numbers.size()<1)
                        {
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find number in line-%-3d  ***\n",linenum);
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given dim=1,2,3 !!!  ***\n");
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                            return false;
                        }
                        ymax=numbers[0];
                    }
                    else
                    {
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find 'ymax=' in line-%3d  ***\n",linenum);
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given ymax=value !!! ***\n");
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                        return false;
                    }

                    if(ymin>ymax)
                    {
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: invalid mesh information        ***\n",linenum);
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        ymin=%8.3f>ymax=%8.3f !!! ***\n",ymin,ymax);
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        please check your input file    ***\n");
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                        return false;
                    }

                    //***for nx
                    getline(in,line);linenum+=1;// read nx
                    str=RemoveSpace(line);
                    if(str.compare(0,3,"nx=")==0)
                    {
                        numbers=SplitNum(line);
                        if(numbers.size()<1)
                        {
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find number in line-%-3d  ***\n",linenum);
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given nx=int value   ***\n");
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                            return false;
                        }
                        nx=numbers[0];
                        if(nx<1)
                        {
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: nx=%5d(line-%3d) is invalid !!***\n",nx,linenum);
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given nx=value(>0)   ***\n");
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                            return false;
                        }
                    }
                    else
                    {
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: 'nx=' can't found in line=%-3d   ***\n",linenum);
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given nx=value(>0)   ***\n");
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                        return false;
                    }
                    //***for ny
                    getline(in,line);linenum+=1;// read ny
                    str=RemoveSpace(line);
                    if(str.compare(0,3,"ny=")==0)
                    {
                        numbers=SplitNum(line);
                        if(numbers.size()<1)
                        {
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find number in line-%-3d  ***\n",linenum);
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given nx=int value   ***\n");
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                            return false;
                        }
                        ny=numbers[0];
                        if(ny<1)
                        {
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: ny=%5d(line-%3d) is invalid !!***\n",nx,linenum);
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given ny=value(>0)   ***\n");
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                            return false;
                        }
                    }
                    else
                    {
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: 'ny=' can't found in line=%-3d   ***\n",linenum);
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given ny=value(>0)   ***\n");
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                        return false;
                    }

                    //**** Read mesh type
                    getline(in,line);linenum+=1;//
                    str=RemoveSpace(line);
                    if(str.compare(0,9,"meshtype=")==0)
                    {
                        if(str.find("quad4")!=string::npos)
                        {
                            meshtype="quad4";
                        }
                        else if(str.find("quad8")!=string::npos)
                        {
                            meshtype="quad8";
                        }
                        else if(str.find("quad9")!=string::npos)
                        {
                            meshtype="quad9";
                        }
                        else
                        {
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: unsupported mesh type(line-%-3d) ***\n",linenum);
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***       'meshtype=quad4,8,9'is expected!!***\n");
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                            return false;
                        }
                    }
                    else
                    {
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error:can't find 'meshtype=' in line%-3d***\n",linenum);
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***       'meshtype=quad4,8,9'is expected!!***\n");
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                        return false;
                    }
                    Mesh2D mesh2D(xmin,xmax,ymin,ymax,nx,ny,meshtype);
                    mesh2D.CreateMesh();
                    mesh.Add2DMesh(mesh2D);
                    mesh.Init();
                }
            }

        }
    }
    else
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***----------------------------------------***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***--- can't find matched [mesh]/[] block--***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***--- please check your input file!!!   --***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!!!                          ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

}

