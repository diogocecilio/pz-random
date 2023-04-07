#pragma once
#include <cmath>
#include <set>
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include "pzgmesh.h"
#include "pzcmesh.h"
#include <pzgmesh.h>
#include <pzcmesh.h>
#include "pzstack.h"
#include "TPZVTKGeoMesh.h"
#include "pzanalysis.h"
#include "pzbndcond.h"
#include "pzfunction.h"
#include <pzskylstrmatrix.h> //symmetric skyline matrix storage
#include <pzstepsolver.h> //for TPZStepSolver
#include "tpzdohrsubstruct.h"
#include "tpzdohrmatrix.h"
#include "tpzdohrprecond.h"
#include "pzdohrstructmatrix.h"
#include <sstream>
#include "pzmaterial.h"
#include <sys/stat.h>
#include <thread>
#include "KLMaterial.h"
#include "tpzgeoelrefpattern.h"
#include <TPZGeoLinear.h>
#include "KLStrMatrix.h"
#include "KLAnalysis.h"
#include "pzelasmatwithmem.h"
#include "pzelasmat.h"
#include "pzelasticmem.h"
#include <random>
#include "TPZYCMohrCoulombPV.h"
#include "pzelastoplastic2D.h"
#include "pzelastoplasticmem.h"
#include "TPZPlasticStepPV.h"
#include "readgidmesh.h"
#include "pzgeopoint.h"
#include <thread>
#include "pzdarcymatwithmem.h"
#include "darcytools.h"
#include "plasticitytools.h"
#include "fieldtools.h"


TPZGeoMesh * CreateGMeshRef ( int ref,int ref2,string file );
TPZGeoMesh * CreateGMesh ( int ref,string file );

void SolveP();

void CreateField();

int main()
{


    SolveP();



    return 0;
}

void CreateField()
{
    int ref=2;
    int ref2=0;
    int porder=1;
    string file ="/home/diogo/projects/pz-random/data/tri-struc.msh";
    TPZGeoMesh *gmesh1= CreateGMeshRef ( ref,ref2,file );
    TPZGeoMesh *gmesh2= CreateGMeshRef ( ref,ref2,file );
    REAL lx=20.;
    REAL ly=2.;
    int M=50;
    int type=3;
    FieldTools fields(gmesh1,gmesh2,lx, ly, M, type, porder);

    fields.ComputeField (  );
}

void SolveP()
{
    int ref=2;
    int ref2=1;
    int porder=2;
    string file ="/home/diogo/projects/pz-random/data/tri-struc.msh";

    //TPZGeoMesh *gmeshpalstic = CreateGMeshRef ( ref,ref2,file );
    TPZGeoMesh *gmeshpalstic = CreateGMesh ( ref,file );
    TPZGeoMesh *gmeshdarcy = CreateGMesh ( ref,file );
    TPZGeoMesh *gmeshfields1 = CreateGMesh ( ref,file );
    TPZGeoMesh *gmeshfields2 = CreateGMesh ( ref,file );

    REAL gammaagua=10.;
    REAL gammasolo=20.;
    int numthreads=10;

    REAL H=10;
    REAL Hw=10.;
    REAL Ht=40.;
    REAL gammaW=10.;
    DarcyTools darcytools(gmeshdarcy,H,Hw,Ht,gammaW,porder);

    darcytools.SolveDarcyProlem();

    PlasticityTools plastictools(gmeshpalstic,gammaagua,gammasolo,numthreads,porder);

    darcytools.SetFlux(plastictools.fcmesh);

    int porderfield=1;
    REAL lx=20.;
    REAL ly=2.;
    int M=50;
    int type=3;
    FieldTools fields(gmeshfields1,gmeshfields2,lx, ly, M, type, porderfield);
    TPZManVector<TPZCompMesh*,2> vecmesh;
    if(false)
    {
        fields.ComputeField (  );
        return;
    }
    else
    {
        vecmesh = fields.SettingCreateFilds();
    }


    plastictools.TransferSolutionFrom(vecmesh,2);

    REAL tolfs =0.01;
    int numiterfs =20;
    REAL tolres = 1.e-6;
    int numiterres =20;
    REAL l =0.2;
    REAL lambda0=0.61;

    REAL fs = plastictools.Solve(tolfs,numiterfs,tolres,numiterres,l,lambda0);

    plastictools.PostPlasticity();

}
TPZGeoMesh * CreateGMeshRef ( int ref,int ref2,string file )
{

    TPZGeoMesh *gmesh  =  new TPZGeoMesh();

    gmesh->SetDimension ( 2 );

    readgidmesh read;

    std::vector<std::vector<int>> meshtopol;
    std::vector<std::vector<double>> meshcoords;

    read.ReadMesh2 ( meshtopol,meshcoords,file );


    cout << "a" << endl;
    int ncoords = meshcoords.size();
    gmesh->NodeVec().Resize ( ncoords );

    TPZVec<REAL> coord ( 2 );
    for ( int inode=0; inode<ncoords; inode++ )
    {
        coord[0] = meshcoords[inode][1];
        coord[1] = meshcoords[inode][2];
        gmesh->NodeVec() [inode] = TPZGeoNode ( inode, coord, *gmesh );
    }

    int sz = meshtopol.size();
    TPZVec <long> TopoQuad ( 4 );
    TPZVec <long> TopoTri ( 3 );
    TPZVec <long> TopoLine ( 2 );
    for ( int iel=0; iel<meshtopol.size(); iel++ )
    {
        if ( meshtopol[iel].size() ==4 )
        {
            TopoTri[0] =meshtopol[iel][1]-1;
            TopoTri[1] =meshtopol[iel][2]-1;
            TopoTri[2] =meshtopol[iel][3]-1;
            new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> ( iel, TopoTri, 1,*gmesh );

        }
        else if ( meshtopol[iel].size() ==5 )
        {

            TopoQuad[0] =meshtopol[iel][1]-1;
            TopoQuad[1] =meshtopol[iel][2]-1;
            TopoQuad[2] =meshtopol[iel][3]-1;
            TopoQuad[3] =meshtopol[iel][4]-1;
            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> ( iel, TopoQuad, 1,*gmesh );

        }
        else if ( meshtopol[iel].size() ==3 )
        {
            TopoLine[0] = meshtopol[iel][1]-1;
            TopoLine[1] = meshtopol[iel][2]-1;
            REAL x0 = meshcoords[TopoLine[0]][1];
            REAL y0 = meshcoords[TopoLine[0]][2];
            REAL xf = meshcoords[TopoLine[1]][1];
            REAL yf = meshcoords[TopoLine[1]][2];
            REAL tol=1.e-3;
            if ( ( fabs ( ( y0-0 ) ) <tol && fabs ( ( yf-0 ) ) <tol ) )
            {
                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -1, *gmesh );//bottom
            }
            else if ( ( fabs ( ( x0-75 ) ) <tol && fabs ( ( xf-75 ) ) <tol ) )
            {
                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -2, *gmesh );//rigth
            }
            else if ( ( fabs ( ( y0-30 ) ) <tol && fabs ( ( yf-30 ) ) <tol ) )
            {
                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -3, *gmesh );//toprigth
            }
            else if ( ( fabs ( ( y0-40 ) ) <tol && fabs ( ( yf-40 ) ) <tol ) )
            {
                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -4, *gmesh );//topleft
            }
            else if ( ( fabs ( ( x0-0 ) ) <tol && fabs ( ( xf-0 ) ) <tol ) )
            {
                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -5, *gmesh );//left
            }
            else if ( ( fabs ( ( xf-x0 ) ) >tol && fabs ( ( yf-y0 ) ) >tol ) )
            {
                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -6, *gmesh );//ramp
            }
            else
            {
                cout<< "bc element not found."<<endl;
                cout<< "x0 = " << x0 << " y0 = "<< y0 << endl;
                cout<< "xf = " << xf << " yf = "<< yf << endl;
                DebugStop();
            }

        }

    }

    TPZVec<long> TopolPoint ( 1 );
    TopolPoint[0]=61;
    new TPZGeoElRefPattern< pzgeom::TPZGeoPoint> ( sz+1, TopolPoint, 2, *gmesh );

    TopoTri[0] =61;
    TopoTri[1] =50;
    TopoTri[2] =53;
    new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> ( sz+1, TopoTri, 3,*gmesh );

    TopoTri[0] =53;
    TopoTri[1] =50;
    TopoTri[2] =38;
    new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> ( sz+1, TopoTri, 4,*gmesh );

    TopoTri[0] =43;
    TopoTri[1] =53;
    TopoTri[2] =38;
    new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> ( sz+1, TopoTri, 5,*gmesh );

    TopoTri[0] =38;
    TopoTri[1] =50;
    TopoTri[2] =36;
    new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> ( sz+1, TopoTri, 6,*gmesh );


    gmesh->BuildConnectivity();
    cout << "c" << endl;
    for ( int d = 0; d<ref; d++ )
    {
        int nel = gmesh->NElements();
        TPZManVector<TPZGeoEl *> subels;
        for ( int iel = 0; iel<nel; iel++ )
        {
            TPZGeoEl *gel = gmesh->ElementVec() [iel];
            gel->Divide ( subels );
        }
    }

    set<int> SETmatRefDir11;
    SETmatRefDir11.insert ( 2 );
    SETmatRefDir11.insert ( 3 );
    SETmatRefDir11.insert ( 4 );
    SETmatRefDir11.insert ( 5 );
    SETmatRefDir11.insert ( 6 );
    for ( int j = 0; j < ref2; j++ )
    {
        long nel = gmesh->NElements();
        for ( long iref = 0; iref < nel; iref++ )
        {
            TPZVec<TPZGeoEl*> filhos;
            TPZGeoEl * gelP11 = gmesh->ElementVec() [iref];
            TPZRefPatternTools::RefineUniformIfNeighMat ( gelP11, SETmatRefDir11 );
        }
    }

    std::ofstream files ( "teste-mesh.vtk" );
    TPZVTKGeoMesh::PrintGMeshVTK ( gmesh,files,false );
    cout << "d" << endl;
    return gmesh;


}
TPZGeoMesh * CreateGMesh ( int ref,string file )
{

    TPZGeoMesh *gmesh  =  new TPZGeoMesh();

    gmesh->SetDimension ( 2 );

    readgidmesh read;

    std::vector<std::vector<int>> meshtopol;
    std::vector<std::vector<double>> meshcoords;

    read.ReadMesh2 ( meshtopol,meshcoords,file );


    cout << "a" << endl;
    int ncoords = meshcoords.size();
    gmesh->NodeVec().Resize ( ncoords );

    TPZVec<REAL> coord ( 2 );
    for ( int inode=0; inode<ncoords; inode++ )
    {
        coord[0] = meshcoords[inode][1];
        coord[1] = meshcoords[inode][2];
        gmesh->NodeVec() [inode] = TPZGeoNode ( inode, coord, *gmesh );
    }

    int sz = meshtopol.size();
    TPZVec <long> TopoQuad ( 4 );
    TPZVec <long> TopoTri ( 3 );
    TPZVec <long> TopoLine ( 2 );
    for ( int iel=0; iel<meshtopol.size(); iel++ )
    {
        if ( meshtopol[iel].size() ==4 )
        {
            TopoTri[0] =meshtopol[iel][1]-1;
            TopoTri[1] =meshtopol[iel][2]-1;
            TopoTri[2] =meshtopol[iel][3]-1;
            new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> ( iel, TopoTri, 1,*gmesh );

        }
        else if ( meshtopol[iel].size() ==5 )
        {

            TopoQuad[0] =meshtopol[iel][1]-1;
            TopoQuad[1] =meshtopol[iel][2]-1;
            TopoQuad[2] =meshtopol[iel][3]-1;
            TopoQuad[3] =meshtopol[iel][4]-1;
            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> ( iel, TopoQuad, 1,*gmesh );

        }
        else if ( meshtopol[iel].size() ==3 )
        {
            TopoLine[0] = meshtopol[iel][1]-1;
            TopoLine[1] = meshtopol[iel][2]-1;
            REAL x0 = meshcoords[TopoLine[0]][1];
            REAL y0 = meshcoords[TopoLine[0]][2];
            REAL xf = meshcoords[TopoLine[1]][1];
            REAL yf = meshcoords[TopoLine[1]][2];
            REAL tol=1.e-3;
            if ( ( fabs ( ( y0-0 ) ) <tol && fabs ( ( yf-0 ) ) <tol ) )
            {
                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -1, *gmesh );//bottom
            }
            else if ( ( fabs ( ( x0-75 ) ) <tol && fabs ( ( xf-75 ) ) <tol ) )
            {
                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -2, *gmesh );//rigth
            }
            else if ( ( fabs ( ( y0-30 ) ) <tol && fabs ( ( yf-30 ) ) <tol ) )
            {
                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -3, *gmesh );//toprigth
            }
            else if ( ( fabs ( ( y0-40 ) ) <tol && fabs ( ( yf-40 ) ) <tol ) )
            {
                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -4, *gmesh );//topleft
            }
            else if ( ( fabs ( ( x0-0 ) ) <tol && fabs ( ( xf-0 ) ) <tol ) )
            {
                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -5, *gmesh );//left
            }
            else if ( ( fabs ( ( xf-x0 ) ) >tol && fabs ( ( yf-y0 ) ) >tol ) )
            {
                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -6, *gmesh );//ramp
            }
            else
            {
                cout<< "bc element not found."<<endl;
                cout<< "x0 = " << x0 << " y0 = "<< y0 << endl;
                cout<< "xf = " << xf << " yf = "<< yf << endl;
                DebugStop();
            }

        }

    }


    gmesh->BuildConnectivity();
    cout << "c" << endl;
    for ( int d = 0; d<ref; d++ )
    {
        int nel = gmesh->NElements();
        TPZManVector<TPZGeoEl *> subels;
        for ( int iel = 0; iel<nel; iel++ )
        {
            TPZGeoEl *gel = gmesh->ElementVec() [iel];
            gel->Divide ( subels );
        }
    }

    std::ofstream files ( "teste-mesh.vtk" );
    TPZVTKGeoMesh::PrintGMeshVTK ( gmesh,files,false );
    cout << "d" << endl;
    return gmesh;


}



