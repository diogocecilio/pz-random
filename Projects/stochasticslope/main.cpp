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
#define EIGEN_USE_MKL_ALL
string filelocation = "/home/diogo/projects/pz-random/data";

TPZGeoMesh *  CreateGMeshSlope ( int ref );
TPZGeoMesh * CreateGMeshRef ( int ref,int ref2,string file );
TPZGeoMesh * CreateGMesh ( int ref,string file );
TPZFMatrix<REAL> ExportField(TPZCompMesh * vecmesh,int isol);
void SolveDeterministic();
void MonteCarlo ( int a,int b );
void MonteCarlo2 ( int a,int b );
void MonteCarloShearRed ( int a,int b );
void SolveMultiThread ( int a0,int b0,int nthreads );
void SolveMultiThread ( void ( *function ) ( int,int ), int a0,int b0,int nthreads );
void SolveDeterministic ( bool gimsrm=1 );
TPZGeoMesh * CreateGMesh ( int ref,string file,std::vector<double> coordbc  );
void  ReadFile ( std::string file,TPZFMatrix<REAL> &out );
template <class T>
std::vector<T> str_vec ( std::vector<std::string> &vs );
void SettingCreateFild (string file,TPZCompMesh * cmesh );

TPZGeoMesh * CreateRefMesh();
//string filelocation2 = "/home/diogo/projects/pz-random-build-release/Projects/stochasticslope/";
void Testing();

int main()
{
//     chrono::steady_clock sc;
//     auto start = sc.now();
//
//    // SolveMultiThread(MonteCarlo2,500,1000,12);
//       MonteCarlo2 ( 0,1 );
//    //SolveDeterministic ( 0);
//
//     auto end = sc.now();
//     auto time_span = static_cast<chrono::duration<double>> ( end - start );
//     cout << " | total time in iterative process =  " << time_span.count() << std::endl;

//Testing();
    // MonteCarloShearRed ( 500,1000 );
    //SolveMultiThread(MonteCarloShearRed,100,200,2);
    //  SolveMultiThread(MonteCarlo,0,1000,12);
    // SolveDeterministic();
    //int srm=1;
    //int gim=0;
    //  SolveDeterministic ( 0 );
MonteCarlo2 ( 0,2 );
    for(int i =0;i<2;i++)
    {
        //ExportField();
    }


    return 0;
}

TPZGeoMesh *  CreateGMeshSlope ( int ref )
{
        const std::string name ( "Darcy Flow Slope" );

        TPZGeoMesh *gmesh  =  new TPZGeoMesh();

        gmesh->SetName ( name );

        gmesh->SetDimension ( 2 );

        TPZVec<REAL> coord ( 2 );

        vector<vector<double>> co= {
                {0., 0.}, {75., 0.}, {75., 30.},{45., 30.},
                {35., 40.}, {0.,40.},{35./3., 40.},{2 * 35/3., 40.},
                {30., 40.},{30., 30.}, {60.,30.},{2* 35./3.,2* 35/3.},
                {45., 2* 35/3.},{35./3., 35/3.}, {60., 35./3.}
        };
        vector<vector<int>> topol = {
                {0,  1,  14, 13},{1,  2,  10, 14}, {14, 10, 3,  12},
                {13, 14, 12, 11},{11, 12, 3,  9}, {9,  3,  4,  8},
                {11, 9,  8,  7},{13, 11, 7, 6},{0, 13,  6, 5}
        };

        gmesh->NodeVec().Resize ( co.size() );

        for ( int inode=0; inode<co.size(); inode++ ) {
                coord[0] = co[inode][0];
                coord[1] = co[inode][1];
                gmesh->NodeVec() [inode] = TPZGeoNode ( inode, coord, *gmesh );
        }
        TPZVec <long> TopoQuad ( 4 );
        for ( int iel=0; iel<topol.size(); iel++ ) {
                TopoQuad[0] = topol[iel][0];
                TopoQuad[1] = topol[iel][1];
                TopoQuad[2] =	topol[iel][2];
                TopoQuad[3] = topol[iel][3];
                new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> ( iel, TopoQuad, 1,*gmesh );
        }




        int id = topol.size();
        TPZVec <long> TopoLine ( 2 );
        TopoLine[0] = 0;
        TopoLine[1] = 1;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -1, *gmesh );//bottom

        id++;
        TopoLine[0] = 1;
        TopoLine[1] = 2;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -2, *gmesh );//rigth

        id++;
        TopoLine[0] = 2;
        TopoLine[1] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -3, *gmesh );//top-rigth

        id++;
        TopoLine[0] = 4;
        TopoLine[1] = 5;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -4, *gmesh ); //top-left

        id++;
        TopoLine[0] = 3;
        TopoLine[1] = 4;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -6, *gmesh ); //ramp

        id++;
        TopoLine[0] = 5;
        TopoLine[1] = 0;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -5, *gmesh ); //left





        gmesh->BuildConnectivity();
        for ( int d = 0; d<ref; d++ ) {
                int nel = gmesh->NElements();
                TPZManVector<TPZGeoEl *> subels;
                for ( int iel = 0; iel<nel; iel++ ) {
                        TPZGeoEl *gel = gmesh->ElementVec() [iel];
                        gel->Divide ( subels );
                }
        }

        string meshref = "gmesh.vtk";
        std::ofstream files ( meshref );
        TPZVTKGeoMesh::PrintGMeshVTK ( gmesh,files,true );
        return gmesh;
}

void  ReadFile ( std::string file,TPZFMatrix<REAL> &out )
{
    string line,line2, temp;

    ifstream myfile ( file );

    std::vector<std::vector<double>> coords;
    int counter=0;
    while ( getline ( myfile, line ) )
    {
        std::vector<string> tokens;
        istringstream iss ( line );
        while ( iss >> temp ) tokens.push_back ( temp );
        std::vector<double> input_doub_temp= str_vec<double> ( tokens );
        coords.push_back ( input_doub_temp );
        counter++;
    }

    out.Resize ( coords.size(),coords[1].size() );
    for ( int iel=0; iel<coords.size(); iel++ )
    {
        for ( int elnode=0; elnode<coords[iel].size(); elnode++ )
        {
            out ( iel,elnode ) = coords[iel][elnode];
        }
    }
}

template <class T>
std::vector<T> str_vec ( std::vector<std::string> &vs )
{
    std::vector<T> ret;
    for ( std::vector<string>::iterator it = vs.begin(); it != vs.end(); ++it )
    {
        istringstream iss ( *it );
        T temp;
        iss >> temp;
        ret.push_back ( temp );
    }
    return ret;
}




void SettingCreateFild (string file,TPZCompMesh * cmesh )
{

    string outco=filelocation2;
    outco+=file;
    cout << "outco = " << outco;
    TPZFMatrix<REAL> readco;
    ReadFile ( outco,readco );
    cmesh->LoadSolution ( readco );

}

void Testing()
{

    int ref=2;
    int porderfield=1;

    string file =filelocation;
    file+="/tri-struc-v2.msh";

    TPZGeoMesh *gmeshfield = CreateGMesh ( ref,file );
    TPZGeoMesh *gmeshfield1 = CreateGMesh ( ref,file );
    TPZGeoMesh *gmeshfield2 = CreateGMesh ( ref,file );

    REAL lx=20.;
    REAL ly=2.;
    int M=150;
    int type=3;
    FieldTools field ( gmeshfield,lx, ly, M, type, porderfield );//aqui leak
    FieldTools field1 ( gmeshfield1,lx, ly, M, type, porderfield );//aqui leak

    lx=8.;
    ly=2.;//indiferente para tipo 4
    M=150;
    type=4;
    FieldTools field2 ( gmeshfield2,lx, ly, M, type, porderfield );//aqui leak


    TPZCompMesh* meshcoes = field.SettingCreateFild ( "cohes.txt" );
    TPZCompMesh* meshphi = field1.SettingCreateFild ( "fric.txt" );
    TPZCompMesh* meshperm = field2.SettingCreateFild ("perm.txt" );


    TPZGeoMesh *gmeshpalstic = CreateGMesh ( ref,file );

    REAL E=20000.;
    REAL nu=0.49;
    REAL coes = 10.;
    REAL phi=30*M_PI/180.;
    REAL gammaagua=0.;
    REAL gammasolo=20.;
    int numthreads=10;
    int porder=2;


    PlasticityTools plastictool( gmeshpalstic,E, nu, coes,phi,gammaagua,gammasolo,numthreads,porder );

   TPZCompMesh *newcmesh(plastictool.fcmesh);

    SettingCreateFild ( "cohes.txt",newcmesh );

    plastictool.TransferingSol(newcmesh);

    //newcmesh->Print(cout);

}




void MonteCarlo2 ( int a,int b )
{
        bool flux=false;
        bool variableparam=true;
        if(flux==false)
        {
            cout<< "\n\n NO uncetainty in permability is taking into account.\n\n";
        }else{
            cout<< "\n\n uncetainty in permability IS taking into account.\n\n";
        }

        if(variableparam==false)
        {
            cout<< "\n\n NO uncetainty in cohesion or friction angle are taking into account.\n\n";
        }else{
            cout<< "\n\n uncetainty in cohesion and friction angle ARE taking into account.\n\n";
        }


    int ref=4;
    int porder=2;
    int porderfield=1;
    int nref=0;

    string file =filelocation;
    file+="/tri-struc-v2.msh";
    std::vector<double> coordbc(3);

    coordbc[0]=75.;coordbc[1]=30.;coordbc[2]=10.;

//     TPZGeoMesh *gmeshfield = CreateGMesh ( ref,file );
//     TPZGeoMesh *gmeshfield1 = CreateGMesh ( ref,file );
//     TPZGeoMesh *gmeshfield2 = CreateGMesh ( ref,file );
    TPZGeoMesh *gmeshfield = CreateGMeshSlope ( ref );
    TPZGeoMesh *gmeshfield1 = CreateGMeshSlope ( ref );
    TPZGeoMesh *gmeshfield2 = CreateGMeshSlope ( ref );
    REAL gammaagua=0.;
    if(flux)
    {
        gammaagua=10.;
    }

    REAL gammasolo=20.;
    REAL coes=10.;
    REAL phi=30.*M_PI/180.;
    REAL E=20000.;
    REAL nu =0.49;
    int numthreads=12;

    int samples =1000;
    REAL lx=20.;
    REAL ly=2.;
    int M=150;
    int type=3;
    FieldTools field ( gmeshfield,lx, ly, M, type, porderfield );//aqui leak
    FieldTools field1 ( gmeshfield1,lx, ly, M, type, porderfield );//aqui leak


    lx=8.;
    ly=2.;//indiferente para tipo 4
    M=150;
    type=4;
    FieldTools field2 ( gmeshfield2,lx, ly, M, type, porderfield );//aqui leak



    TPZVec<REAL> mean(2),cov(2);
    TPZVec<string> filefields(2);
    TPZVec<REAL> meanperm(1),covperm(1);
    TPZVec<string> filefieldsperm(1);
    mean[0]=10.;mean[1]=30*M_PI/180.;
    cov[0]=0.3;cov[1]=0.2;
    filefields[0]="cohes.txt";filefields[1]="fric.txt";

    meanperm[0]=1;
    filefieldsperm[0]="perm.txt";
    covperm[0]=0.5;
    TPZManVector<TPZCompMesh*,3> vecmesh ( 3 );
    if ( true )
    {
        field.ComputeFieldExport ( mean ,   cov , filefields, samples );
        //field2.ComputeFields ( meanperm ,   covperm , filefieldsperm, samples );
        return;
    }
    else
    {
        vecmesh[0] = field.SettingCreateFild ( filefields[0] );
        vecmesh[1] = field1.SettingCreateFild ( filefields[1] );
        vecmesh[2] = field2.SettingCreateFild ( filefieldsperm[0] );
    }

    std::ofstream print ( "teste2" );

    TPZElastoPlasticAnalysis * anal = new TPZElastoPlasticAnalysis(vecmesh[0]);

    anal->AssembleResidual();
    anal->LoadSolution(anal->Rhs());

//       TPZFMatrix<REAL> pzdata = ExportField(anal->Mesh(),0 );
//
//
//           int row=pzdata.Rows();
//     int cols = pzdata.Cols();
//
//     for ( int i=0; i<row; i++ )
//     {
//         for ( int j=0; j<cols; j++ )
//         {
//             print << pzdata ( i,j ) << " ";
//         }
//         print<< std::endl;
//     }
//
//     print<<"end"<< endl;

      return;
    for ( int imc=0; imc<=2; imc++ )
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
        REAL tolres = 1.e-3;
        int numiterres =30;
        REAL l =0.5;
        REAL lambda0=0.5;

        int iref;


        REAL H=10.;
        REAL Hw=10.;
        REAL Ht=40.;

        DarcyTools darcytools ( gmeshdarcy,H,Hw,Ht,gammaagua,porder);
       // PlasticityTools plastictoolstempx ( gmeshpalstic,E, nu, coes,phi,gammaagua,gammasolo,numthreads,porder );






        continue;
        for ( iref=1; iref<=nref; iref++ )
        {

            // TPZGeoMesh *gmeshpalstic2 = CreateGMesh ( ref,file );



            PlasticityTools plastictoolstemp ( gmeshpalstic,E, nu, coes,phi,gammaagua,gammasolo,numthreads,porder );



            std::set<long> elindices;
            DarcyTools darcytools ( gmeshdarcy,H,Hw,Ht,gammaagua,porder);
            if(flux)
            {
                cout << "refing..setting flux.."<< endl;
                darcytools.TransferSolutionFrom ( vecmesh,imc );
                darcytools.SolveDarcyProlem();
                darcytools.SetFlux ( plastictoolstemp.fcmesh );
            }

            if(variableparam)
            {
                cout << "refing..setting cohesion and friction angle.."<< endl;
                plastictoolstemp.TransferSolutionFrom ( vecmesh,imc );
            }
            cout << "solving.."<< endl;
            plastictoolstemp.Solve ( tolfs,numiterfs,tolres,numiterres,l,lambda0 );
            plastictoolstemp.ComputeElementDeformation();
            //plastictoolstemp.PRefineElementsAbove ( 0.0015,porder+1,elindices );
            plastictoolstemp.DivideElementsAbove ( 0.0015,elindices );
            gmeshpalstic=plastictoolstemp.fcmesh->Reference();

            // delete gmeshpalstic2;


        }
        PlasticityTools plastictools ( gmeshpalstic,E, nu, coes,phi,gammaagua,gammasolo,numthreads,porder );
        tolfs =0.01;
        numiterfs =20;
        tolres = 1.e-6;
        numiterres =20;
        l =0.5;
        lambda0=0.5;
        cout << "END REFING.."<< endl;
        //DarcyTools darcytools ( gmeshdarcy,H,Hw,Ht,gammaagua,porder );
        if(flux)
        {
            cout << " setting flux"<< endl;
            darcytools.TransferSolutionFrom ( vecmesh,imc );
            darcytools.SolveDarcyProlem();
            darcytools.SetFlux ( plastictools.fcmesh );
            darcytools.PostDarcy ( vtk2 );
        }


        string meshref = "post/printrefmesh";
        meshref+=var;
        meshref+=".vtk";
        std::ofstream files ( meshref );
        TPZVTKGeoMesh::PrintGMeshVTK ( gmeshpalstic,files,true );


        if(variableparam)
        {
            cout << "..setting cohesion and friction..."<<endl;
            plastictools.TransferSolutionFrom ( vecmesh,imc );
        }

        cout << "..solving with arclength..."<<endl;
        bool converge;
        REAL fs = plastictools.Solve ( tolfs,numiterfs,tolres,numiterres,l,lambda0,converge );

        cout << "..postprocessing..."<<endl;
        plastictools.PostPlasticity ( vtk1 );
        out << fs << std::endl;
        out << converge << std::endl;

        delete gmeshpalstic;
        delete gmeshdarcy;
    }
}


TPZFMatrix<REAL> ExportField(TPZCompMesh * vecmesh,int isol)
{

    std::vector<std::vector<double>> datastd;
    std::ofstream print ( "teste" );

    int nels =  vecmesh->NElements();

    TPZMaterial *mat1 = vecmesh->FindMaterial ( 1 );
    if ( !mat1 )
    {
        DebugStop();
    }

    for ( int iel=0; iel<nels; iel++ )
    {
        TPZCompEl *cel = vecmesh->ElementVec() [iel];
        if ( !cel )
        {
            continue;
        }
        TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *> ( cel );
        if ( !intel )
        {
            continue;
        }

        TPZGeoEl *gel = cel->Reference();

        TPZMaterialData data;
        data.fNeedsSol = true;
        intel->InitMaterialData ( data );
        if ( true )
        {

            REAL co[8][2] = {{-1,-1},{1,-1},{1,1},{-1,1},{0,-1},{1,0},{0,1},{-1,0}};
            for ( int p = 0; p<1; p++ )
            {


                TPZManVector<REAL,3> qsi ( 2,0. ),xco ( 3,0. );
                qsi[0] = 0;
                qsi[1] =0;
                //gel->X ( qsi, xco );
                intel->ComputeRequiredData ( data, qsi );
                intel->ComputeSolution(qsi,data);
                 REAL x,y,z,solu;
                 x=data.x[0];
                 y=data.x[1];
                 z=data.x[2];
                 solu=data.sol[0][0];
//
                 if(solu>1.e-3)
                 {
                 std::vector<double> localdata(4);
                 print << x << " "<< y<< " " << z << " "<< solu <<std::endl;
                 localdata[0]=x;
                 localdata[1]=y;
                 localdata[2]=z;
                 localdata[3]=solu;
                 datastd.push_back(localdata);
                 }
            }
        }
    }

    int sz=datastd.size();
    TPZFMatrix<REAL> pzdata(sz,4);
    for(int i=0;i<sz;i++)
    {
        pzdata(i,0)=datastd[i][0];
        pzdata(i,1)=datastd[i][1];
        pzdata(i,2)=datastd[i][2];
        pzdata(i,3)=datastd[i][3];
    }
    return pzdata;
}



void MonteCarlo ( int a,int b )
{
    int ref=3;
    int porder=2;
    int porderfield=1;
    int nref=2;


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

    int ref=1;
    int porder=3;
    int nref=3;

    string file =filelocation;
    file+="/tri-struc-v2.msh";
   // file+="/mesh2x1cho.msh";
   // file+="/mesh2x1choeq.msh";
    //file+="/gidmesh.msh";
    std::vector<double> coordbc(3);
    coordbc[0]=30.;coordbc[1]=5.;coordbc[2]=5.;
    coordbc[0]=75.;coordbc[1]=30.;coordbc[2]=10.;
    //coordbc[0]=75.;coordbc[1]=30.;coordbc[2]=5.;
    TPZGeoMesh *gmeshpalstic = CreateGMesh ( ref,file,coordbc );
    //CreateGMesh ( int ref,string file,std::vector<double> coordbc )
    TPZGeoMesh *gmeshdarcy = CreateGMesh ( ref,file,coordbc  );

    REAL gammaagua=0*10.;
    REAL gammasolo=20.;
   // REAL coes=10.;
   // REAL phi=30.*M_PI/180.;
    REAL coes=10.;
    REAL phi=30*M_PI/180.;
    REAL E=20000.;
    REAL nu =0.49;
    int numthreads=20;

    REAL tolfs =0.1;
    int numiterfs =100;
    REAL tolres = 1.e-6;
    int numiterres =10;
    REAL l =1;
    REAL lambda0=0.1;


    int iref;
    lambda0=0.1;
    l=0.5;
    bool converge;
    cout << "Iniciando o refinamento da malha... " << endl;
    std::set<long> elindices2;
    for ( iref=0; iref<nref; iref++ )
    {


        PlasticityTools plastictoolstemp ( gmeshpalstic,E, nu, coes,phi,gammaagua,gammasolo,numthreads,porder );
        std::set<long> elindices;
         REAL fsdummy = plastictoolstemp.ShearRed ( );
        // plastictoolstemp.Solve ( tolfs,numiterfs,tolres,numiterres,l,lambda0 );
        //plastictoolstemp.IterativeProcessHybridArcLength( tolfs, numiterfs, tolres, numiterres, l, lambda0, converge);
        plastictoolstemp.ComputeElementDeformation();
        //plastictoolstemp.PRefineElementsAbove ( 0.0015,porder+iref,elindices );
        plastictoolstemp.DivideElementsAbove ( 0.001,elindices );
        gmeshpalstic=plastictoolstemp.fcmesh->Reference();
        gmeshdarcy=plastictoolstemp.fcmesh->Reference();
        // elindices2=elindices;

    }

    cout << "pos-processando.. " << iref << endl;
    string meshref = "malharefinada.vtk";
    std::ofstream files ( meshref );
    TPZVTKGeoMesh::PrintGMeshVTK ( gmeshpalstic,files,true );

    l =0.1;
    lambda0=0.1;
    tolfs =0.01;
    numiterfs =100;
    tolres = 1.e-6;
    numiterres =30;
    //lambda0=0.1;
    cout << "Fim do refinamento da malha"<< endl;

    PlasticityTools plastictoolstemp2 ( gmeshpalstic,E, nu, coes,phi,gammaagua,gammasolo,numthreads,porder );
    //PlasticityTools plastictools ( gmeshpalstic,E, nu, coes,phi,gammaagua,gammasolo,numthreads,porder );
    int setflux=1;
    if ( false )
    {

        cout << "Resolvendo problema de fluxo"<< endl;
        REAL H=10;
        REAL Hw=10;
        REAL Ht=40.;

        DarcyTools darcytools ( gmeshdarcy,H,Hw,Ht,gammaagua,porder );
        darcytools.SolveDarcyProlem();
        darcytools.SetFlux ( plastictoolstemp2.fcmesh );
        string vtk2="daArcyx.vtk";
        darcytools.PostDarcy ( vtk2 );

    }

    if ( true )
    {
        cout << "Saindo do refinando e resolvendo na malha fina com o metodo da reducao da resistencia.. " << endl;
        plastictoolstemp2.ShearRed ( );
    }
    else
    {
        cout << "Saindo do refinando e resolvendo na malha fina com o metodo da plasticidade incremental com o Arc-Length.. " << endl;
        plastictoolstemp2.Solve ( tolfs,numiterfs,tolres,numiterres,l,lambda0 );
    }





    string vtk1 = "saida-talude.vtk";
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

