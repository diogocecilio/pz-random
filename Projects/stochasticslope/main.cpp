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

string filelocation = "/home/diogo/projects/pz-random/data";


TPZGeoMesh * CreateGMeshRef ( int ref,int ref2,string file );
TPZGeoMesh * CreateGMesh ( int ref,string file );

void SolveDeterministic();
void MonteCarlo ( int a,int b );
void MonteCarlo2 ( int a,int b );
void MonteCarloShearRed ( int a,int b );
void SolveMultiThread ( int a0,int b0,int nthreads );
void SolveMultiThread ( void ( *function ) ( int,int ), int a0,int b0,int nthreads );
void SolveDeterministic ( bool gimsrm=1 );
TPZGeoMesh * CreateGMesh ( int ref,string file,std::vector<double> coordbc  );

TPZGeoMesh * CreateRefMesh();

int main()
{
    chrono::steady_clock sc;
    auto start = sc.now();

    //SolveMultiThread(MonteCarlo2,400,500,16);
    //MonteCarlo2 ( 0,100 );
    SolveDeterministic ( 0 );

    auto end = sc.now();
    auto time_span = static_cast<chrono::duration<double>> ( end - start );
    cout << " | total time in iterative process =  " << time_span.count() << std::endl;


    // MonteCarloShearRed ( 500,1000 );
    //SolveMultiThread(MonteCarloShearRed,100,200,2);
    //  SolveMultiThread(MonteCarlo,0,1000,12);
    // SolveDeterministic();
    //int srm=1;
    //int gim=0;
    //  SolveDeterministic ( 0 );
    return 0;
}







void MonteCarlo2 ( int a,int b )
{
    int ref=3;
    int porder=2;
    int porderfield=1;
    int nref=3;

    string file =filelocation;
    file+="/tri-struc-v2.msh";

    TPZGeoMesh *gmeshfield = CreateGMesh ( ref,file );
    TPZGeoMesh *gmeshfield2 = CreateGMesh ( ref,file );
    TPZGeoMesh *gmeshfield3 = CreateGMesh ( ref,file );

    REAL gammaagua=0.;
    REAL gammasolo=20.;
    REAL coes=10.;
    REAL phi=30.*M_PI/180.;
    REAL E=20000.;
    REAL nu =0.49;
    int numthreads=10;

    int samples =1000;
    REAL lx=20.;
    REAL ly=2.;
    int M=150;
    int type=3;
    string file1  =  "cohes.txt";
    FieldTools field ( gmeshfield,lx, ly, M, type, porderfield );//aqui leak

    lx=20.;
    ly=2.;
    M=150;
    type=3;
    string file2  =  "fric.txt";
    FieldTools fieldphi ( gmeshfield2,lx, ly, M, type, porderfield );//aqui leak

    lx=20.;
    ly=20.;
    M=150;
    type=4;
    string permfile="permeability.txt";
    FieldTools fieldpermeability ( gmeshfield3,lx, ly, M, type, porderfield,permfile );//aqui leak


    TPZManVector<TPZCompMesh*,3> vecmesh ( 3 );
    if ( false )
    {
        REAL mean1=10.;
        REAL cov1=0.3;
        REAL mean2=30.*M_PI/180.;
        REAL cov2=0.2;
        REAL crossfac=-0.5;
        //ComputeFields ( REAL mean1,REAL mean2,REAL cov1, REAL cov2,int samples,REAL crossfac,string file1,string file2)
        field.ComputeFields ( mean1, mean2, cov1,  cov2, samples, crossfac, file1, file2 );

        mean1=1.;
        cov1=0.2;
        fieldpermeability.ComputeField ( mean1,cov1,samples );

        return;
    }
    else
    {
        vecmesh[0] = field.SettingCreateFild ( file1 );
        vecmesh[1] = fieldphi.SettingCreateFild ( file2 );
        vecmesh[2] = fieldpermeability.SettingCreateFild ( permfile ); //corrigir file
    }

    for ( int imc=a; imc<b; imc++ )
    {


        TPZGeoMesh *gmeshdarcy = CreateGMesh ( ref,file );
        TPZGeoMesh *gmeshpalstic = CreateGMesh ( ref,file );
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


        REAL tolfs =0.01;
        int numiterfs =30;
        REAL tolres = 1.e-6;
        int numiterres =30;
        REAL l =0.1;
        REAL lambda0=0.1;


        REAL fsr=0.;
        int iref;
        for ( iref=1; iref<=nref; iref++ )
        {

            // TPZGeoMesh *gmeshpalstic2 = CreateGMesh ( ref,file );

            PlasticityTools plastictoolstemp ( gmeshpalstic,E, nu, coes,phi,gammaagua,gammasolo,numthreads,porder );
            std::set<long> elindices;
            // REAL fsdummy = plastictoolstemp.ShearRed ( vecmesh,imc );
            plastictoolstemp.TransferSolutionFrom ( vecmesh,imc );
            fsr=plastictoolstemp.Solve ( tolfs,numiterfs,tolres,numiterres,l,lambda0 );
            plastictoolstemp.ComputeElementDeformation();
            //plastictoolstemp.PRefineElementsAbove ( 0.0015,porder+1,elindices );
            plastictoolstemp.DivideElementsAbove ( 0.0015,elindices );
            gmeshpalstic=plastictoolstemp.fcmesh->Reference();

            // delete gmeshpalstic2;


        }
        tolfs =0.01;
        numiterfs =20;
        tolres = 1.e-6;
        numiterres =20;
        l =0.1;
        lambda0=0.2;
        cout << "Fim do refinamento da malha"<< endl;
        // DarcyTools darcytools ( gmeshdarcy,H,Hw,Ht,gammaagua,porder+iref );
        // darcytools.TransferSolutionFrom ( vecmesh,imc );
        // darcytools.SolveDarcyProlem();
        // darcytools.PostDarcy ( vtk2 );
        cout << "end refining, solving..."<<endl;

        string meshref = "post/printrefmesh";
        meshref+=var;
        meshref+=".vtk";
        std::ofstream files ( meshref );
        TPZVTKGeoMesh::PrintGMeshVTK ( gmeshpalstic,files,true );

        PlasticityTools plastictools ( gmeshpalstic,E, nu, coes,phi,gammaagua,gammasolo,numthreads,porder );
        //darcytools.SetFlux ( plastictools.fcmesh );
        plastictools.TransferSolutionFrom ( vecmesh,imc );

        //REAL fs = plastictools.ShearRed ( vecmesh,imc );
        bool converge;
        REAL fs = plastictools.Solve ( tolfs,numiterfs,tolres,numiterres,l,lambda0,converge );

        plastictools.PostPlasticity ( vtk1 );
        out << fs << std::endl;
        out << converge << std::endl;

        delete gmeshpalstic;
        delete gmeshdarcy;
    }
}

void MonteCarlo ( int a,int b )
{
    int ref=3;
    int porder=2;
    int porderfield=1;
    int nref=3;


    string file =filelocation;
    //file+="/tri-struc.msh";
    file+="/tri-struc-v2.msh";
    //string file ="/home/diogo/projects/pz-random/data/tri-struc-v2.msh";


    TPZGeoMesh *gmeshfields1 = CreateGMesh ( ref,file );
    TPZGeoMesh *gmeshfields2 = CreateGMesh ( ref,file );
    TPZGeoMesh *gmeshfields3 = CreateGMesh ( ref,file );


    REAL gammaagua=0.;
    REAL gammasolo=20.;
    REAL coes=10.;
    REAL phi=30.*M_PI/180.;
    REAL E=20000.;
    REAL nu =0.49;
    int numthreads=10;

    REAL H=10;
    REAL Hw=10.;
    REAL Ht=40.;


    int samples =1000;
    REAL lx=20.;
    REAL ly=2.;
    int M=150;
    int type=3;
    string coesfile="coesao-v2.txt";
    FieldTools fieldscoes ( gmeshfields1,lx, ly, M, type, porderfield,coesfile );//aqui leak

    lx=20.;
    ly=2.;
    M=150;
    type=3;
    string phifile="atrito-v2.txt";
    FieldTools fieldphi ( gmeshfields2,lx, ly, M, type, porderfield,phifile );//aqui leak

    lx=20.;
    ly=20.;
    M=150;
    type=4;
    string permfile="permeability-v2.txt";
    FieldTools fieldpermeability ( gmeshfields3,lx, ly, M, type, porderfield,permfile );//aqui leak

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
        TPZGeoMesh *gmeshpalstic = CreateGMesh ( ref,file );
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


        REAL tolfs =0.01;
        int numiterfs =20;
        REAL tolres = 1.e-5;
        int numiterres =20;
        REAL l =0.1;
        REAL lambda0=0.1;


        int iref;
        for ( iref=1; iref<=nref; iref++ )
        {

            // TPZGeoMesh *gmeshpalstic2 = CreateGMesh ( ref,file );

            PlasticityTools plastictoolstemp ( gmeshpalstic,E, nu, coes,phi,gammaagua,gammasolo,numthreads,porder );
            std::set<long> elindices;
            // REAL fsdummy = plastictoolstemp.ShearRed ( vecmesh,imc );
            plastictoolstemp.TransferSolutionFrom ( vecmesh,imc );
            plastictoolstemp.Solve ( tolfs,numiterfs,tolres,numiterres,l,lambda0 );
            plastictoolstemp.ComputeElementDeformation();
            //plastictoolstemp.PRefineElementsAbove ( 0.0015,porder+1,elindices );
            plastictoolstemp.DivideElementsAbove ( 0.0015,elindices );
            gmeshpalstic=plastictoolstemp.fcmesh->Reference();

            // delete gmeshpalstic2;


        }
        tolfs =0.01;
        numiterfs =20;
        tolres = 1.e-6;
        numiterres =20;
        l =0.2;
        lambda0=0.1;
        cout << "Fim do refinamento da malha"<< endl;
        // DarcyTools darcytools ( gmeshdarcy,H,Hw,Ht,gammaagua,porder+iref );
        // darcytools.TransferSolutionFrom ( vecmesh,imc );
        // darcytools.SolveDarcyProlem();
        // darcytools.PostDarcy ( vtk2 );
        cout << "end refining, solving..."<<endl;

        string meshref = "post/printrefmesh";
        meshref+=var;
        meshref+=".vtk";
        std::ofstream files ( meshref );
        TPZVTKGeoMesh::PrintGMeshVTK ( gmeshpalstic,files,true );

        PlasticityTools plastictools ( gmeshpalstic,E, nu, coes,phi,gammaagua,gammasolo,numthreads,porder );
        //darcytools.SetFlux ( plastictools.fcmesh );
        plastictools.TransferSolutionFrom ( vecmesh,imc );

        //REAL fs = plastictools.ShearRed ( vecmesh,imc );
        bool converge;
        REAL fs = plastictools.Solve ( tolfs,numiterfs,tolres,numiterres,l,lambda0,converge );

        plastictools.PostPlasticity ( vtk1 );
        out << fs << std::endl;
        out << converge << std::endl;

        delete gmeshpalstic;
        delete gmeshdarcy;
    }
}


void SolveDeterministic ( bool gimsrm )
{

    int ref=4;
    int porder=2;
    int nref=1;

    string file =filelocation;
    //file+="/tri-struc.msh";
    //file+="/mesh2x1cho.msh";
    file+="/mesh2x1choeq.msh";

    std::vector<double> coordbc(3);
    coordbc[0]=30.;coordbc[1]=5.;coordbc[2]=5.;
    //coordbc[0]=75.;coordbc[1]=30.;coordbc[2]=10.;
    TPZGeoMesh *gmeshpalstic = CreateGMesh ( ref,file,coordbc );
    TPZGeoMesh *gmeshdarcy = CreateGMesh ( ref,file,coordbc );

    REAL gammaagua=0.;
    REAL gammasolo=20.;
    REAL coes=23.;
    REAL phi=1.*M_PI/180.;
    REAL E=20000.;
    REAL nu =0.49;
    int numthreads=10;

    REAL tolfs =0.01;
    int numiterfs =100;
    REAL tolres = 1.e-6;
    int numiterres =30;
    REAL l =0.02;
    REAL lambda0=0.1;


    int iref;

    cout << "Iniciando o refinamento da malha... " << endl;
    std::set<long> elindices2;
    for ( iref=0; iref<nref; iref++ )
    {
        PlasticityTools plastictoolstemp ( gmeshpalstic,E, nu, coes,phi,gammaagua,gammasolo,numthreads,porder );
        std::set<long> elindices;
        // REAL fsdummy = plastictoolstemp.ShearRed ( );
        plastictoolstemp.Solve ( tolfs,numiterfs,tolres,numiterres,l,lambda0 );
        plastictoolstemp.ComputeElementDeformation();
        //plastictoolstemp.PRefineElementsAbove ( 0.0015,porder+iref,elindices );
        plastictoolstemp.DivideElementsAbove ( 0.0015,elindices );
        gmeshpalstic=plastictoolstemp.fcmesh->Reference();
        gmeshdarcy=plastictoolstemp.fcmesh->Reference();
        // elindices2=elindices;

    }

    tolfs =0.01;
    numiterfs =100;
    tolres = 1.e-6;
    numiterres =30;
    lambda0=0.1;
    cout << "Fim do refinamento da malha"<< endl;

    PlasticityTools plastictoolstemp2 ( gmeshpalstic,E, nu, coes,phi,gammaagua,gammasolo,numthreads,porder );
    //PlasticityTools plastictools ( gmeshpalstic,E, nu, coes,phi,gammaagua,gammasolo,numthreads,porder );
    int setflux=0;
    if ( setflux )
    {

        REAL H=10;
        REAL Hw=0;
        REAL Ht=40.;

        DarcyTools darcytools ( gmeshdarcy,H,Hw,Ht,gammaagua,porder );
        darcytools.SolveDarcyProlem();
        darcytools.SetFlux ( plastictoolstemp2.fcmesh );
        string vtk2="daArcyx.vtk";
        darcytools.PostDarcy ( vtk2 );

    }

    if ( gimsrm )
    {
        cout << "Saindo do refinando e resolvendo na malha fina com o metodo da reducao da resistencia.. " << endl;
        plastictoolstemp2.ShearRed ( );
    }
    else
    {
        cout << "Saindo do refinando e resolvendo na malha fina com o metodo da plasticidade incremental com o Arc-Length.. " << endl;
        plastictoolstemp2.Solve ( tolfs,numiterfs,tolres,numiterres,l,lambda0 );
    }


    cout << "pos-processando.. " << iref << endl;
    string meshref = "malharefinada.vtk";
    std::ofstream files ( meshref );
    TPZVTKGeoMesh::PrintGMeshVTK ( gmeshpalstic,files,true );


    string vtk1 = "printdeterr.vtk";
    plastictoolstemp2.PostPlasticity ( vtk1 );




}

void SolveMultiThread ( void ( *function ) ( int,int ), int a0,int b0,int nthreads )
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

        std::thread threadx ( function,a,b );

        threadsmat1.push_back ( std::move ( threadx ) );

        a=b+1;
        b+=delta;
    }

    for ( auto &threadx: threadsmat1 ) threadx.join();
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
    //TopolPoint[0]=61;//tri-struc
    TopolPoint[0]=17;//tri-struc2
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
    TPZVTKGeoMesh::PrintGMeshVTK ( gmesh,files,true );
    cout << "d" << endl;
    return gmesh;


}



TPZGeoMesh * CreateGMesh ( int ref,string file,std::vector<double> coordbc )
{


    REAL L=coordbc[0];
    REAL h1=coordbc[1];
    REAL h2=coordbc[2];

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
            else if ( ( fabs ( ( x0-L ) ) <tol && fabs ( ( xf-L ) ) <tol ) )
            {
                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -2, *gmesh );//rigth
            }
            else if ( ( fabs ( ( y0-h1 ) ) <tol && fabs ( ( yf-h1 ) ) <tol ) )
            {
                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -3, *gmesh );//toprigth
            }
            else if ( ( fabs ( ( y0-(h1+h2) ) ) <tol && fabs ( ( yf-(h1+h2) ) ) <tol ) )
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

void MonteCarloShearRed ( int a,int b )
{
    int ref=2;
    int ref2=2;
    int porder=2;
    int porderfield=1;
    int nref=2;


    string file =filelocation;
    file+="/tri-struc-v2.msh";


    TPZGeoMesh *gmeshfields1 = CreateGMesh ( ref,file );
    TPZGeoMesh *gmeshfields2 = CreateGMesh ( ref,file );
    TPZGeoMesh *gmeshfields3 = CreateGMesh ( ref,file );


    REAL gammaagua=0.;
    REAL gammasolo=20.;
    REAL coes=10.;
    REAL phi=30.*M_PI/180.;
    REAL E=20000.;
    REAL nu =0.49;
    int numthreads=10;

    REAL H=10;
    REAL Hw=10.;
    REAL Ht=40.;


    int samples =1000;
    REAL lx=20.;
    REAL ly=2.;
    int M=100;
    int type=3;
    string coesfile="coesao.txt";
    FieldTools fieldscoes ( gmeshfields1,lx, ly, M, type, porderfield,coesfile );

    lx=20.;
    ly=2.;
    M=100;
    type=3;
    string phifile="atrito.txt";
    FieldTools fieldphi ( gmeshfields2,lx, ly, M, type, porderfield,phifile );

    lx=20.;
    ly=20.;
    M=100;
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




        int iref;
        for ( iref=1; iref<=nref; iref++ )
        {

            PlasticityTools plastictoolstemp ( gmeshpalstic,E, nu, coes,phi,gammaagua,gammasolo,numthreads,porderfield );
            std::set<long> elindices;
            REAL fsdummy = plastictoolstemp.ShearRed ( vecmesh,imc );
            plastictoolstemp.ComputeElementDeformation();
            //plastictoolstemp.PRefineElementsAbove ( 0.0015,porder+iref,elindices );
            plastictoolstemp.DivideElementsAbove ( 0.0015,elindices );
            gmeshpalstic=plastictoolstemp.fcmesh->Reference();

        }

        // DarcyTools darcytools ( gmeshdarcy,H,Hw,Ht,gammaagua,porder+iref );
        // darcytools.TransferSolutionFrom ( vecmesh,imc );
        // darcytools.SolveDarcyProlem();
        // darcytools.PostDarcy ( vtk2 );
        cout << "end refining, solving..."<<endl;

        string meshref = "post/printrefmesh";
        meshref+=var;
        meshref+=".vtk";
        std::ofstream files ( meshref );
        TPZVTKGeoMesh::PrintGMeshVTK ( gmeshpalstic,files,true );

        PlasticityTools plastictools ( gmeshpalstic,E, nu, coes,phi,gammaagua,gammasolo,numthreads,porder );
        //darcytools.SetFlux ( plastictools.fcmesh );
        plastictools.TransferSolutionFrom ( vecmesh,imc );

        REAL fs = plastictools.ShearRed ( vecmesh,imc );

        plastictools.PostPlasticity ( vtk1 );
        out << fs << std::endl;
    }
}

