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

void SolveDeterministic();
void MonteCarlo ( int a,int b );
void MonteCarloShearRed ( int a,int b );
void SolveMultiThread ( int a0,int b0,int nthreads );

TPZGeoMesh * CreateRefMesh();

int main()
{


    //SolveP();

   // MonteCarloShearRed ( 0,1000 );
    SolveMultiThread(0,1000,10);

    // SolveDeterministic();
    return 0;
}

void SolveMultiThread ( int a0,int b0,int nthreads )
{
    int samples=b0-a0;
    std::vector <std::thread> threadsmat1,threadsmat2;

    int delta = int ( samples/nthreads );
    int a=a0;
    int b=a+delta;

    for ( int i=0; i<nthreads; i++ )
    {
        std::cout << "a = "<< a <<std::endl;
        std::cout << "b = "<< b <<std::endl;

        string file ="/home/diogo/projects/pz-random/data/tri-struc.msh";
        int porder=1;
        int ref=2;
        bool createfield=false;

        TPZGeoMesh *gmesh1 = CreateGMesh ( ref,file );
        TPZGeoMesh *gmesh2 = CreateGMesh ( ref,file );

        std::thread threadx ( MonteCarloShearRed,a,b );

        threadsmat1.push_back ( std::move ( threadx ) );

        a=b+1;
        b+=delta;
    }

    for ( auto &threadx: threadsmat1 ) threadx.join();
}
void SolveDeterministic()
{
    int ref=1;
    int porder=1;
    //string file ="/home/diogo/projects/pz-random/data/tri-struc2x1.msh";//2:1
    //string file ="/home/diogo/projects/pz-random/data/tri.msh";//2:1
    string file ="/home/diogo/projects/pz-random/data/tri-struc-v2.msh";//1:1

    TPZGeoMesh *gmeshpalstic = CreateGMeshRef ( ref,2,file );
    TPZGeoMesh *gmeshdarcy = CreateGMesh ( ref,file );
    TPZGeoMesh *gmeshfields1 = CreateGMesh ( ref,file );
    TPZGeoMesh *gmeshfields2 = CreateGMesh ( ref,file );
    TPZGeoMesh *gmeshfields3 = CreateGMesh ( ref,file );

    TPZManVector<TPZGeoMesh*,3> gmeshvec ( 3 );
    gmeshvec[0]=gmeshfields1;
    gmeshvec[1]=gmeshfields2;
    gmeshvec[2]=gmeshfields3;

    REAL gammaagua=10.;
    REAL gammasolo=20.;
    REAL coes=10.;
    REAL phi=30.;
    REAL E=20000.;
    REAL nu =0.49;
    int numthreads=10;

    REAL H=10;
    REAL Hw=10.;
    REAL Ht=40.;


    PlasticityTools plastictools ( gmeshpalstic,E, nu, coes,phi,gammaagua*0,gammasolo,numthreads,porder );

    string vtkd="darcydeter.vtk";

    DarcyTools darcytools ( gmeshdarcy,H,Hw,Ht,gammaagua,porder+1 );
    darcytools.SolveDarcyProlem();
    //darcytools.SetFlux ( plastictools.fcmesh );
    darcytools.PostDarcy ( vtkd );

    REAL tolfs =0.01;
    int numiterfs =20;
    REAL tolres = 1.e-6;
    int numiterres =20;
    REAL l =0.2;
    REAL lambda0=0.61;

    // REAL fs = plastictools.Solve ( tolfs,numiterfs,tolres,numiterres,l,lambda0 );

    for ( int iref=0; iref<2; iref++ )
    {
        std::set<long> elindices;
        REAL fs = plastictools.ShearRed();
        plastictools.ComputeElementDeformation();
        plastictools.DivideElementsAbove ( 0.005,elindices );
        plastictools.ComputeElementDeformation();
    }

    std::ofstream files ( "printrefmesh.vtk" );
    TPZVTKGeoMesh::PrintGMeshVTK ( plastictools.fcmesh->Reference(),files,false );

    //TPZGeoMesh
    PlasticityTools plastictools2 ( plastictools.fcmesh->Reference(),E, nu, coes,phi,gammaagua*0,gammasolo,numthreads,porder+1 );
    darcytools.SetFlux ( plastictools2.fcmesh );
    plastictools2.ShearRed();
    //REAL fs = plastictools.GravityIncrease();
    string vtkp="plasticdeter.vtk";
    plastictools2.PostPlasticity ( vtkp );

}

void MonteCarloShearRed ( int a,int b )
{
    int ref=2;
    int ref2=2;
    int porder=1;
    int porderfield=1;

    string file ="/home/diogo/projects/pz-random/data/tri-struc-v2.msh";


    TPZGeoMesh *gmeshfields1 = CreateGMesh ( ref,file );
    TPZGeoMesh *gmeshfields2 = CreateGMesh ( ref,file );
    TPZGeoMesh *gmeshfields3 = CreateGMesh ( ref,file );


    REAL gammaagua=0.;
    REAL gammasolo=20.;
    REAL coes=10.;
    REAL phi=30.*M_PI/180.;
    REAL E=20000.;
    REAL nu =0.49;
    int numthreads=15;

    REAL H=10;
    REAL Hw=10.;
    REAL Ht=40.;


    int samples =1000;
    REAL lx=20.;
    REAL ly=2.;
    int M=50;
    int type=3;
    string coesfile="coesao.txt";
    FieldTools fieldscoes ( gmeshfields1,lx, ly, M, type, porderfield,coesfile );

    lx=20.;
    ly=2.;
    M=50;
    type=3;
    string phifile="atrito.txt";
    FieldTools fieldphi ( gmeshfields2,lx, ly, M, type, porderfield,phifile );

    lx=20.;
    ly=20.;
    M=50;
    type=4;
    string permfile="permeability.txt";
    FieldTools fieldpermeability ( gmeshfields3,lx, ly, M, type, porderfield,permfile );

    TPZManVector<TPZCompMesh*,3> vecmesh ( 3 );
    if ( false )
    {
        REAL mean=10.;
        REAL cov=0.3;
        fieldscoes.ComputeField ( mean,cov,samples );


        mean=30.*M_PI/180.;
        cov=0.2;
        fieldphi.ComputeField ( mean,cov,samples );


        mean=1.;
        cov=0.2;
        fieldpermeability.ComputeField ( mean,cov,samples );


        vecmesh[0]= fieldscoes.GetCMesh();
        vecmesh[1]= fieldphi.GetCMesh();
        vecmesh[2]=fieldpermeability.GetCMesh();

        return;
    }
    else
    {
        vecmesh[0] = fieldscoes.SettingCreateFild();
        vecmesh[1] = fieldphi.SettingCreateFild();
        vecmesh[2] = fieldpermeability.SettingCreateFild();
    }

    for ( int imc=a; imc<b; imc++ )
    {


        TPZGeoMesh *gmeshdarcy = CreateGMesh ( ref,file );
        TPZGeoMesh *gmeshpalstic = CreateGMeshRef ( ref,ref2,file );
        string saidafs = "post/fs";
        string vtk1 = "postvtk/saidamontecarloplasticity";
        string vtk2 = "postvtk/saidamontecarlodarcy";
        auto var=to_string ( imc );
        saidafs+=var;
        vtk1+=var;
        vtk2+=var;
        saidafs+=".dat";
        vtk1+=".vtk";
        vtk2+=".vtk";
        ofstream out ( saidafs );
        std::cout << "imc = " <<  imc << std::endl;


        PlasticityTools plastictoolstemp ( gmeshpalstic,E, nu, coes,phi,gammaagua,gammasolo,numthreads,porder );

        int iref;
        for ( iref=0; iref<2; iref++ )
        {

            std::set<long> elindices;
            REAL fsdummy = plastictoolstemp.ShearRed ( vecmesh,imc );
            plastictoolstemp.ComputeElementDeformation();
            plastictoolstemp.PRefineElementsAbove ( 0.0015,porder+iref,elindices );
            plastictoolstemp.DivideElementsAbove ( 0.0015,elindices );
            gmeshpalstic=plastictoolstemp.fcmesh->Reference();

        }

        DarcyTools darcytools ( gmeshdarcy,H,Hw,Ht,gammaagua,porder+iref);
        darcytools.TransferSolutionFrom ( vecmesh,imc );
        darcytools.SolveDarcyProlem();
        darcytools.PostDarcy ( vtk2 );
        cout << "end refining, solving..."<<endl;

        string meshref = "printrefmesh";
        meshref+=var;
        meshref+=".vtk";
        std::ofstream files ( meshref );
        TPZVTKGeoMesh::PrintGMeshVTK ( gmeshpalstic,files,false );

   //     darcytools.SetFlux ( plastictoolstemp.fcmesh );
        plastictoolstemp.TransferSolutionFrom ( vecmesh,imc );

        REAL fs = plastictoolstemp.ShearRed ( vecmesh,imc );

        plastictoolstemp.PostPlasticity ( vtk1 );
        out << fs << std::endl;
    }
}

// void MonteCarlo ( int a,int b )
// {
//     int ref=2;
//     int ref2=2;
//     int porder=2;
//     string file ="/home/diogo/projects/pz-random/data/tri-struc.msh";
//
//     TPZGeoMesh *gmeshpalstic = CreateGMeshRef ( ref,ref2,file );
//     //TPZGeoMesh *gmeshpalstic = CreateGMesh ( ref,file );
//     TPZGeoMesh *gmeshdarcy = CreateGMesh ( ref,file );
//     TPZGeoMesh *gmeshfields1 = CreateGMesh ( ref,file );
//     TPZGeoMesh *gmeshfields2 = CreateGMesh ( ref,file );
//     TPZGeoMesh *gmeshfields3 = CreateGMesh ( ref,file );
//
//     TPZManVector<TPZGeoMesh*,3> gmeshvec ( 3 );
//     gmeshvec[0]=gmeshfields1;
//     gmeshvec[1]=gmeshfields2;
//     gmeshvec[2]=gmeshfields3;
//
//     REAL gammaagua=10.;
//     REAL gammasolo=20.;
//     REAL coes=10.;
//     REAL phi=30.;
//     REAL E=20000.;
//     REAL nu =0.49;
//     int numthreads=10;
//
//     REAL H=10;
//     REAL Hw=10.;
//     REAL Ht=40.;
//
//     DarcyTools darcytools ( gmeshdarcy,H,Hw,Ht,gammaagua,porder );
//     PlasticityTools plastictools ( gmeshpalstic,E, nu, coes,phi,gammaagua,gammasolo,numthreads,porder );
//
//
//     REAL lx=20.;
//     REAL ly=2.;
//     int M=50;
//     int type=3;
//     int porderfield=1;
//
//     FieldTools fields ( gmeshvec,lx, ly, M, type, porderfield );
//     TPZManVector<TPZCompMesh*,3> vecmesh;
//     if ( false )
//     {
//         fields.ComputeField ( );
//         return;
//     }
//     else
//     {
//         vecmesh = fields.SettingCreateFilds();
//     }
//
//
//     for ( int imc=a; imc<b; imc++ )
//     {
//         string saidafs = "post/fs";
//         string vtk1 = "postvtk/saidamontecarloplasticity";
//         string vtk2 = "postvtk/saidamontecarlodarcy";
//         auto var=to_string ( imc );
//         saidafs+=var;
//         vtk1+=var;
//         vtk2+=var;
//         saidafs+=".dat";
//         vtk1+=".vtk";
//         vtk2+=".vtk";
//         ofstream out ( saidafs );
//         std::cout << "imc = " <<  imc << std::endl;
//         darcytools.TransferSolutionFrom ( vecmesh,imc );
//         darcytools.SolveDarcyProlem();
//         darcytools.SetFlux ( plastictools.fcmesh );
//         darcytools.PostDarcy ( vtk2 );
//
//         plastictools.TransferSolutionFrom ( vecmesh,imc );
//
//         REAL tolfs =0.01;
//         int numiterfs =20;
//         REAL tolres = 1.e-6;
//         int numiterres =20;
//         REAL l =0.2;
//         REAL lambda0=0.61;
//
//         REAL fs = plastictools.Solve ( tolfs,numiterfs,tolres,numiterres,l,lambda0 );
//
//         plastictools.PostPlasticity ( vtk1 );
//         out << fs << std::endl;
//     }
// }

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
    TopolPoint[0]=17;
    new TPZGeoElRefPattern< pzgeom::TPZGeoPoint> ( sz+1, TopolPoint, 2, *gmesh );

//     TopoTri[0] =61;
//     TopoTri[1] =50;
//     TopoTri[2] =53;
//     new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> ( sz+1, TopoTri, 3,*gmesh );
//
//     TopoTri[0] =53;
//     TopoTri[1] =50;
//     TopoTri[2] =38;
//     new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> ( sz+1, TopoTri, 4,*gmesh );
//
//     TopoTri[0] =43;
//     TopoTri[1] =53;
//     TopoTri[2] =38;
//     new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> ( sz+1, TopoTri, 5,*gmesh );
//
//     TopoTri[0] =38;
//     TopoTri[1] =50;
//     TopoTri[2] =36;
//     new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> ( sz+1, TopoTri, 6,*gmesh );

//     TPZVec<long> TopolPoint ( 1 );
//     TopolPoint[0]=37;
//     new TPZGeoElRefPattern< pzgeom::TPZGeoPoint> ( sz+1, TopolPoint, 2, *gmesh );
//
//     TopoTri[0] =37;
//     TopoTri[1] =46;
//     TopoTri[2] =42;
//     new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> ( sz+2, TopoTri, 3,*gmesh );



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
//   SETmatRefDir11.insert ( 3 );
//    SETmatRefDir11.insert ( 4 );
//    SETmatRefDir11.insert ( 5 );
    //   SETmatRefDir11.insert ( 6 );
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

    std::ofstream files ( "teste-meshref.vtk" );
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



