
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

typedef TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> LEMC;

typedef   TPZMatElastoPlastic2D <LEMC, TPZElastoPlasticMem > plasticmat;

TPZGeoMesh * CreateGMesh(int ref,string file);

TPZGeoMesh * CreateGMeshRef(int ref,string file);

TPZCompMesh * CreateCompMesh(TPZGeoMesh * gmesh,int porder);

TPZCompMesh * CreateCompMeshKL(TPZGeoMesh * gmesh,int porder);

TPZElastoPlasticAnalysis * CreateAnalysis(TPZCompMesh *cmesh);

REAL GravityIncrease ( TPZCompMesh * cmesh );

void  CreatePostProcessingMesh ( TPZPostProcAnalysis * PostProcess,TPZCompMesh * cmesh );

void PostProcessVariables ( TPZStack<std::string> &scalNames, TPZStack<std::string> &vecNames );

TPZCompMesh * ComputeField(TPZGeoMesh *gmesh,int porder);

TPZFMatrix<REAL> CreateLogNormalRandomField(TPZFMatrix<REAL> PHI, REAL mean, REAL cov,int samples,string outdata);

void PrintMat ( std::string out,TPZFMatrix<REAL> mat );

REAL SolveElastoplastic(TPZManVector<TPZCompMesh*,2> vecmesh,int imc, string vtk,int porder,int ref);

void MonteCarlo();

void LoadingRamp ( TPZCompMesh * cmesh,  REAL factor,REAL gammasolo, REAL gammaagua);

void  ReadFile ( std::string file,TPZFMatrix<REAL> &out );

template <class T>
std::vector<T> str_vec ( std::vector<std::string> &vs );

void SolveDeterministic();

void TransferSolution(TPZManVector<TPZCompMesh*,2> vecmesh,TPZCompMesh * cmesh,int isol);
void ComputeSolution ( TPZCompEl *cel, TPZFMatrix<REAL> &phi,TPZSolVec &sol );
int main()
{


//      SolveDeterministic();
// //
//      return 0;

    MonteCarlo();

    return 0;
}

void SolveDeterministic()
{
    int ref=1;
    int porder=1;
    string file ="/home/diogo/projects/pz-random/data/mesh-tri.msh";
    TPZGeoMesh *gmesh = CreateGMeshRef(ref,file);



    TPZCompMesh * cmesh = CreateCompMesh(gmesh,porder);

    REAL factor =1.;
    REAL gammasolo=20.;
    REAL gammaagua=0.;
    LoadingRamp ( cmesh,factor  , gammasolo,  gammaagua);

    TPZElastoPlasticAnalysis *analysis =  CreateAnalysis(cmesh);

   GravityIncrease ( cmesh );

//         REAL tolfs =0.01;
//     int numiterfs =40;
//     REAL tolres = 1.e-6;
//     int numiterres =20;
//     REAL l =0.1;
//     REAL lambda0=1.;
//    REAL fs = analysis->IterativeProcessArcLength(tolfs,numiterfs,tolres,numiterres,l,lambda0);

     TPZPostProcAnalysis * postprocdeter = new TPZPostProcAnalysis();
     CreatePostProcessingMesh ( postprocdeter, cmesh );

    TPZVec<int> PostProcMatIds ( 1,1 );

    TPZStack<std::string> PostProcVars, scalNames, vecNames;

    PostProcessVariables ( scalNames, vecNames );

    string vtk = "postprocessdeter.vtk";
    postprocdeter->DefineGraphMesh ( 2,scalNames,vecNames,vtk );

    postprocdeter->PostProcess ( 0 );
}

void MonteCarlo()
{
    int ref=2;
    int porder=1;
    //string file ="/home/diogo/projects/pz-random/data/mesh.msh";
    string file ="/home/diogo/projects/pz-random/data/mesh-tri3.msh";
    TPZGeoMesh *gmesh1 = CreateGMesh(ref,file);
      TPZGeoMesh *gmesh2 = CreateGMesh(ref,file);

    TPZCompMesh * cmeshfield;
    TPZCompMesh*  cmeshfield2;
    if(1)
    {
         cmeshfield =  ComputeField(gmesh1,porder);
         return;
    }else{
         cmeshfield =  CreateCompMeshKL(gmesh1,porder);
         cmeshfield2 =  CreateCompMeshKL(gmesh2,porder);
         string outco="/home/diogo/projects/pz-random-build-release/Projects/randomfemplasticityv2/coesao.txt";
         TPZFMatrix<REAL> readco;
         ReadFile ( outco,readco );
         cmeshfield->LoadSolution(readco);

         TPZElastoPlasticAnalysis * analysis1x  = new TPZElastoPlasticAnalysis(cmeshfield);
         analysis1x->LoadSolution(readco);


         string outphi="/home/diogo/projects/pz-random-build-release/Projects/randomfemplasticityv2/atrito.txt";
         TPZFMatrix<REAL> readphi;
         ReadFile ( outphi,readphi );
         cmeshfield2->LoadSolution(readphi);

        TPZElastoPlasticAnalysis * analysis2x  = new TPZElastoPlasticAnalysis(cmeshfield2);
        analysis2x->LoadSolution(readphi);
    }

    TPZManVector<TPZCompMesh*,2> vecmesh(2);
    vecmesh[0]=cmeshfield;
    vecmesh[1]=cmeshfield2;

    int samples = cmeshfield->Solution().Cols();

    std::cout << "samples = " <<  samples << std::endl;

    ofstream out("fs2.dat");

    string vtk = "saidamontecarlo";
    for(int imc=0;imc<samples;imc++)
    {
        auto var=to_string(imc);
        var+=".vtk";
        vtk+=var;
        std::cout << "imc = " <<  imc << std::endl;
        REAL fs= SolveElastoplastic(vecmesh,imc,vtk,porder,ref);
        out << fs << std::endl;
    }


}

REAL SolveElastoplastic(TPZManVector<TPZCompMesh*,2>  vecmesh, int imc,string vtk,int porder,int ref)
{

    int porderint=1;//porder nao importa?
    int refint =ref;//ref deve ser maior ou igual ao ref do kl?
    //TPZGeoMesh *gmesh = rcmesh1->Reference();
    //string file1 ="/home/diogo/projects/pz-random/data/mesh.msh";
    string file1 ="/home/diogo/projects/pz-random/data/mesh-tri3.msh";
    TPZGeoMesh *gmesh = CreateGMeshRef(refint,file1);

    TPZCompMesh * cmesh = CreateCompMesh(gmesh,porderint);

    REAL factor =1.;
    REAL gammasolo=20.;
    REAL gammaagua=0.;
    LoadingRamp ( cmesh,factor  , gammasolo,  gammaagua);

    TransferSolution(vecmesh,cmesh,imc);
    //SetRandomField(rcmesh1,rcmesh2,cmesh,imc);

    TPZElastoPlasticAnalysis *analysis =  CreateAnalysis(cmesh);

     REAL fs = GravityIncrease ( cmesh );


//     REAL tolfs =0.1;
//     int numiterfs =40;
//     REAL tolres = 1.e-3;
//     int numiterres =20;
//     REAL l =0.1;
//     REAL lambda0=1.;
//    REAL fs = analysis->IterativeProcessArcLength(tolfs,numiterfs,tolres,numiterres,l,lambda0);

//     TPZVec<REAL> xd(3);
//     xd[0]=34.99;
//     xd[1]= 39.99;
//     REAL sol = fabs(FindSol(cmesh,xd));

//    std::cout << "top displacement = " << sol << std::endl;

     TPZPostProcAnalysis * postprocdeter = new TPZPostProcAnalysis();
     CreatePostProcessingMesh ( postprocdeter, cmesh );

    TPZVec<int> PostProcMatIds ( 1,1 );

    TPZStack<std::string> PostProcVars, scalNames, vecNames;

    PostProcessVariables ( scalNames, vecNames );

    postprocdeter->DefineGraphMesh ( 2,scalNames,vecNames,vtk );

    postprocdeter->PostProcess ( 0 );

    return fs;
}


TPZCompMesh * ComputeField(TPZGeoMesh *gmesh,int porder)
{

    TPZCompMesh * cmesh = CreateCompMeshKL(gmesh,porder);

    KLAnalysis * klanal = new KLAnalysis ( cmesh );

    KLMaterial *mat = dynamic_cast<KLMaterial*> (cmesh->MaterialVec()[1]);

    klanal->SetExpansionOrder ( mat->GetExpansioOrder() );

    klanal->Solve();



    TPZFMatrix<REAL> eigenfunctions = klanal->Solution();
    REAL mean=10.;
    REAL cov=0.3;
    int samples=1000;
    string outdata = "coesao.txt";

    TPZFMatrix<REAL> field = CreateLogNormalRandomField(eigenfunctions, mean, cov, samples, outdata);


    mean=30.*M_PI/180.;
    cov=0.2;
    string outdataatrito = "atrito.txt";
    TPZFMatrix<REAL> fieldphi = CreateLogNormalRandomField(eigenfunctions, mean, cov, samples, outdataatrito);

    cmesh->LoadSolution(field);

    string file1 ="coesao2.vtk";

    klanal->Post ( file1,2,0 );

    return cmesh;
}

TPZFMatrix<REAL> CreateLogNormalRandomField(TPZFMatrix<REAL> PHI, REAL mean, REAL cov,int samples,string outdata)
{
    TPZFMatrix<REAL>  PHIt;

     PHI.Transpose ( &PHIt );

     int M = PHI.Cols();

     cout << "M = " << M <<endl;

     std::normal_distribution<double> distribution ( 0., 1. );

    TPZFMatrix<REAL>  THETA ( M, samples, 0. );
    for ( int isample = 0; isample < samples; isample++ ) {
        for ( int irdvar = 0; irdvar < M; irdvar++ ) {
            std::random_device rd{};
            std::mt19937 generator{ rd() };
            REAL xic = distribution ( generator );
            REAL xiphi = distribution ( generator );
            THETA(irdvar,isample) = xic;
        }
    }


      TPZFMatrix<REAL>  hhat;

      PHI.Multiply( THETA, hhat );

    REAL sdev = cov * mean;
    REAL xi = sqrt ( log ( 1 + pow ( ( sdev / mean ),2 ) ) );
    REAL lambda = log ( mean ) - xi * xi / 2.;
    for ( int i = 0; i < hhat.Rows(); i++ ) {
        for ( int j = 0; j < hhat.Cols(); j++ ) {
			REAL temp =  hhat(i,j);
            hhat(i,j) = exp ( lambda + xi * temp );
        }
    }

    //string out="testefield.txt";

    PrintMat ( outdata,hhat);

     return hhat;
}

TPZGeoMesh * CreateGMesh(int ref,string file)
{

    TPZGeoMesh *gmesh  =  new TPZGeoMesh();

    gmesh->SetDimension ( 2 );

    readgidmesh read;

    std::vector<std::vector<int>> meshtopol;
    std::vector<std::vector<double>> meshcoords;

    read.ReadMesh2(meshtopol,meshcoords,file);


    cout << "a" << endl;
    int ncoords = meshcoords.size();
    gmesh->NodeVec().Resize ( ncoords );

    TPZVec<REAL> coord(2);
    for ( int inode=0; inode<ncoords; inode++ ) {
        coord[0] = meshcoords[inode][1];
        coord[1] = meshcoords[inode][2];
        gmesh->NodeVec() [inode] = TPZGeoNode ( inode, coord, *gmesh );
    }

        int sz = meshtopol.size();
        TPZVec <long> TopoQuad ( 4 );
        TPZVec <long> TopoTri ( 3 );
        TPZVec <long> TopoLine ( 2 );
        for ( int iel=0; iel<meshtopol.size(); iel++ ) {
            if ( meshtopol[iel].size() ==4 ) {
                TopoTri[0] =meshtopol[iel][1]-1;
                TopoTri[1] =meshtopol[iel][2]-1;
                TopoTri[2] =meshtopol[iel][3]-1;
                new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> ( iel, TopoTri, 1,*gmesh );

            }else if(meshtopol[iel].size() ==5){

                TopoQuad[0] =meshtopol[iel][1]-1;
                TopoQuad[1] =meshtopol[iel][2]-1;
                TopoQuad[2] =meshtopol[iel][3]-1;
                TopoQuad[3] =meshtopol[iel][4]-1;
                new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> ( iel, TopoQuad, 1,*gmesh );

            }else if(meshtopol[iel].size() ==3){
                TopoLine[0] = meshtopol[iel][1]-1;
                TopoLine[1] = meshtopol[iel][2]-1;
                REAL x0 = meshcoords[TopoLine[0]][1];
                REAL y0 = meshcoords[TopoLine[0]][2];
                REAL xf = meshcoords[TopoLine[1]][1];
                REAL yf = meshcoords[TopoLine[1]][2];
                REAL tol=1.e-3;
                if((fabs((y0-0))<tol && fabs((yf-0))<tol))
                {
                    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -1, *gmesh );//bottom
                }
                else if((fabs((x0-75))<tol && fabs((xf-75))<tol))
                {
                    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -2, *gmesh );//rigth
                }
                else if((fabs((y0-30))<tol && fabs((yf-30))<tol))
                {
                    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -3, *gmesh );//toprigth
                }
                else if((fabs((y0-40))<tol && fabs((yf-40))<tol))
                {
                    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -4, *gmesh );//topleft
                }
                else if((fabs((x0-0))<tol && fabs((xf-0))<tol))
                {
                    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -5, *gmesh );//left
                }
                else if((fabs((xf-x0))>tol && fabs((yf-y0))>tol))
                {
                    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -6, *gmesh );//ramp
                }
                else{
                    cout<< "bc element not found."<<endl;
                    cout<< "x0 = " << x0 << " y0 = "<< y0 << endl;
                    cout<< "xf = " << xf << " yf = "<< yf << endl;
                    DebugStop();
                }

            }

        }

        TPZVec<long> TopolPoint(1);
        TopolPoint[0]=28;
        new TPZGeoElRefPattern< pzgeom::TPZGeoPoint> ( sz+1, TopolPoint, 2, *gmesh );

    gmesh->BuildConnectivity();
    cout << "c" << endl;
    for ( int d = 0; d<ref; d++ ) {
        int nel = gmesh->NElements();
        TPZManVector<TPZGeoEl *> subels;
        for ( int iel = 0; iel<nel; iel++ ) {
            TPZGeoEl *gel = gmesh->ElementVec() [iel];
            gel->Divide ( subels );
        }
    }

    std::ofstream files ( "teste-mesh.vtk" );
    TPZVTKGeoMesh::PrintGMeshVTK ( gmesh,files,false );
cout << "d" << endl;
    return gmesh;


}

TPZGeoMesh * CreateGMeshRef(int ref,string file)
{

    TPZGeoMesh *gmesh  =  new TPZGeoMesh();

    gmesh->SetDimension ( 2 );

    readgidmesh read;

    std::vector<std::vector<int>> meshtopol;
    std::vector<std::vector<double>> meshcoords;

    read.ReadMesh2(meshtopol,meshcoords,file);


    cout << "a" << endl;
    int ncoords = meshcoords.size();
    gmesh->NodeVec().Resize ( ncoords );

    TPZVec<REAL> coord(2);
    for ( int inode=0; inode<ncoords; inode++ ) {
        coord[0] = meshcoords[inode][1];
        coord[1] = meshcoords[inode][2];
        gmesh->NodeVec() [inode] = TPZGeoNode ( inode, coord, *gmesh );
    }

        int sz = meshtopol.size();
        TPZVec <long> TopoQuad ( 4 );
        TPZVec <long> TopoTri ( 3 );
        TPZVec <long> TopoLine ( 2 );
        for ( int iel=0; iel<meshtopol.size(); iel++ ) {
            if ( meshtopol[iel].size() ==4 ) {
                TopoTri[0] =meshtopol[iel][1]-1;
                TopoTri[1] =meshtopol[iel][2]-1;
                TopoTri[2] =meshtopol[iel][3]-1;
                new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> ( iel, TopoTri, 1,*gmesh );

            }else if(meshtopol[iel].size() ==5){

                TopoQuad[0] =meshtopol[iel][1]-1;
                TopoQuad[1] =meshtopol[iel][2]-1;
                TopoQuad[2] =meshtopol[iel][3]-1;
                TopoQuad[3] =meshtopol[iel][4]-1;
                new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> ( iel, TopoQuad, 1,*gmesh );

            }else if(meshtopol[iel].size() ==3){
                TopoLine[0] = meshtopol[iel][1]-1;
                TopoLine[1] = meshtopol[iel][2]-1;
                REAL x0 = meshcoords[TopoLine[0]][1];
                REAL y0 = meshcoords[TopoLine[0]][2];
                REAL xf = meshcoords[TopoLine[1]][1];
                REAL yf = meshcoords[TopoLine[1]][2];
                REAL tol=1.e-3;
                if((fabs((y0-0))<tol && fabs((yf-0))<tol))
                {
                    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -1, *gmesh );//bottom
                }
                else if((fabs((x0-75))<tol && fabs((xf-75))<tol))
                {
                    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -2, *gmesh );//rigth
                }
                else if((fabs((y0-30))<tol && fabs((yf-30))<tol))
                {
                    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -3, *gmesh );//toprigth
                }
                else if((fabs((y0-40))<tol && fabs((yf-40))<tol))
                {
                    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -4, *gmesh );//topleft
                }
                else if((fabs((x0-0))<tol && fabs((xf-0))<tol))
                {
                    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -5, *gmesh );//left
                }
                else if((fabs((xf-x0))>tol && fabs((yf-y0))>tol))
                {
                    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, TopoLine, -6, *gmesh );//ramp
                }
                else{
                    cout<< "bc element not found."<<endl;
                    cout<< "x0 = " << x0 << " y0 = "<< y0 << endl;
                    cout<< "xf = " << xf << " yf = "<< yf << endl;
                    DebugStop();
                }

            }

        }

        TPZVec<long> TopolPoint(1);
        TopolPoint[0]=90;
        new TPZGeoElRefPattern< pzgeom::TPZGeoPoint> ( sz+1, TopolPoint, 2, *gmesh );

    gmesh->BuildConnectivity();
    cout << "c" << endl;
    for ( int d = 0; d<ref; d++ ) {
        int nel = gmesh->NElements();
        TPZManVector<TPZGeoEl *> subels;
        for ( int iel = 0; iel<nel; iel++ ) {
            TPZGeoEl *gel = gmesh->ElementVec() [iel];
            gel->Divide ( subels );
        }
    }

        set<int> SETmatRefDir11;
        SETmatRefDir11.insert(2);
        for(int j = 0; j < 6; j++)
		{
			long nel = gmesh->NElements();
			for (long iref = 0; iref < nel; iref++)
			{
				TPZVec<TPZGeoEl*> filhos;
				TPZGeoEl * gelP11 = gmesh->ElementVec()[iref];
                TPZRefPatternTools::RefineUniformIfNeighMat(gelP11, SETmatRefDir11);
			}
		}

    std::ofstream files ( "teste-mesh.vtk" );
    TPZVTKGeoMesh::PrintGMeshVTK ( gmesh,files,false );
cout << "d" << endl;
    return gmesh;


}

TPZCompMesh * CreateCompMesh(TPZGeoMesh * gmesh,int porder)
{
  unsigned int dim  = 2;
    const std::string name ( "ElastoPlastic COMP MESH Footing Problem " );

    TPZCompMesh * cmesh =  new TPZCompMesh ( gmesh );
    cmesh->SetName ( name );
    cmesh->SetDefaultOrder ( porder );
    cmesh->SetDimModel ( dim );

    // Mohr Coulomb data
	REAL mc_cohesion    = 10.;//kpa
    REAL mc_phi         = 30.*M_PI/180.;
    REAL mc_psi         = mc_phi;

    /// ElastoPlastic Material using Mohr Coulomb
    // Elastic predictor
    TPZElasticResponse ER;
    REAL nu = 0.49;
    REAL E = 20000.;

    LEMC plasticstep;
    //TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> LEMC;
    //TPZPlasticStepVoigt<TPZMohrCoulombVoigt, TPZElasticResponse> LEMC;
    ER.SetUp ( E, nu );
    plasticstep.fER =ER;
    // LEMC.SetElasticResponse( ER );
    plasticstep.fYC.SetUp ( mc_phi, mc_psi, mc_cohesion, ER );

    int PlaneStrain = 1;

    plasticmat * material = new plasticmat ( 1,PlaneStrain );
    //plasticmatcrisfield * material = new plasticmatcrisfield ( 1,PlaneStrain );
    material->SetPlasticity ( plasticstep );

    material->SetId ( 1 );

    material->SetLoadFactor(1.);
    material->SetWhichLoadVector(0);//option to compute the total internal force vecor fi=(Bt sigma+ N (b+gradu))

    cmesh->InsertMaterialObject ( material );

    TPZFMatrix<STATE> val1 ( 2,2,0. );
    TPZFMatrix<STATE>  val2 ( 2,1,0. );
    int directionadirichlet =3;
    int newman =1;
    val2 ( 0,0 ) = 1;
    val2 ( 1,0 ) = 1;
    auto * bc_bottom = material->CreateBC ( material, -1,directionadirichlet, val1, val2 );//bottom
    val2 ( 0,0 ) = 1;
    val2 ( 1,0 ) = 0;
    auto * bc_rigth = material->CreateBC ( material, -2, directionadirichlet, val1, val2 );//rigth
    val2 ( 0,0 ) = 1;
    val2 ( 1,0 ) = 0;
    auto * bc_left = material->CreateBC ( material, -5, directionadirichlet, val1, val2 );//left

    val2 ( 0,0 ) = 0;
    val2 ( 1,0 ) = 0;
    auto * bc_topr = material->CreateBC ( material, -3,newman, val1, val2 );//top rigth
    val2 ( 0,0 ) = 0;
    val2 ( 1,0 ) = 0;
    auto * bc_topl = material->CreateBC ( material, -4, newman, val1, val2 );//top left
    val2 ( 0,0 ) = 0;
    val2 ( 1,0 ) = 0;
    auto * bc_ramp = material->CreateBC ( material, -6, newman, val1, val2 );//ramp


    cmesh->InsertMaterialObject ( bc_bottom );
    cmesh->InsertMaterialObject ( bc_rigth );
    cmesh->InsertMaterialObject ( bc_left );
    cmesh->InsertMaterialObject ( bc_topr );
    cmesh->InsertMaterialObject ( bc_topl );
    cmesh->InsertMaterialObject ( bc_ramp );


    //cmesh->InsertMaterialObject ( top );
    cmesh->SetAllCreateFunctionsContinuousWithMem();

    cmesh->AutoBuild();

    return cmesh;



}

TPZElastoPlasticAnalysis * CreateAnalysis(TPZCompMesh *cmesh)
{

    int numthreads=10;

    TPZElastoPlasticAnalysis * analysis =  new TPZElastoPlasticAnalysis ( cmesh ); // Create analysis

    TPZSkylineStructMatrix matskl ( cmesh );
    //TPZFStructMatrix matskl ( cmesh );

    matskl.SetNumThreads ( numthreads );

    analysis->SetStructuralMatrix ( matskl );

    ///Setting a direct solver
    TPZStepSolver<STATE> step;
    step.SetDirect ( ELDLt );
    //step.SetDirect ( ECholesky );
    analysis->SetSolver ( step );



    return analysis;

}

#include "pzfstrmatrix.h"
#include <pzgeopoint.h>
#include <pzgeopoint.h>

void  CreatePostProcessingMesh ( TPZPostProcAnalysis * PostProcess,TPZCompMesh * cmesh )
{
    if ( PostProcess->ReferenceCompMesh() != cmesh ) {

        PostProcess->SetCompMesh ( cmesh );

        TPZVec<int> PostProcMatIds ( 1,1 );
        TPZStack<std::string> PostProcVars, scalNames, vecNames;
        PostProcessVariables ( scalNames, vecNames );

        for ( int i=0; i<scalNames.size(); i++ ) {
            PostProcVars.Push ( scalNames[i] );
        }
        for ( int i=0; i<vecNames.size(); i++ ) {
            PostProcVars.Push ( vecNames[i] );
        }
        //
        PostProcess->SetPostProcessVariables ( PostProcMatIds, PostProcVars );
        TPZFStructMatrix structmatrix ( PostProcess->Mesh() );
        structmatrix.SetNumThreads ( 0 );
        PostProcess->SetStructuralMatrix ( structmatrix );
    }
    //
    //Chamar com o analysis e nao com o postanalysis pois tem o acumulo de sols
    PostProcess->TransferSolution();

}

/// Get the post processing variables
void PostProcessVariables ( TPZStack<std::string> &scalNames, TPZStack<std::string> &vecNames )
{
scalNames.Push ( "StrainVol" );
    scalNames.Push ( "StrainXX" );
    scalNames.Push ( "StrainYY" );
    scalNames.Push ( "StrainZZ" );
    scalNames.Push ( "StrainXY" );
    scalNames.Push ( "StrainXZ" );
    scalNames.Push ( "StrainYZ" );

    scalNames.Push ( "ElStrainVol" );
    scalNames.Push ( "ElStrainXX" );
    scalNames.Push ( "ElStrainYY" );
    scalNames.Push ( "ElStrainZZ" );
    scalNames.Push ( "ElStrainXY" );
    scalNames.Push ( "ElStrainXZ" );
    scalNames.Push ( "ElStrainYZ" );

    scalNames.Push ( "PlStrainVol" );
    scalNames.Push ( "PlStrainXX" );
    scalNames.Push ( "PlStrainYY" );
    scalNames.Push ( "PlStrainZZ" );
    scalNames.Push ( "PlStrainXY" );
    scalNames.Push ( "PlStrainXZ" );
    scalNames.Push ( "PlStrainYZ" );

    scalNames.Push ( "PlStrainSqJ2" );
    scalNames.Push ( "PlStrainSqJ2El" );
    scalNames.Push ( "PlAlpha" );

    scalNames.Push ( "DisplacementX" );
    scalNames.Push ( "DisplacementY" );
    scalNames.Push ( "DisplacementZ" );
    vecNames.Push ( "DisplacementTotal" );


    scalNames.Push ( "YieldSurface1" );
    scalNames.Push ( "YieldSurface2" );
    scalNames.Push ( "YieldSurface3" );

    scalNames.Push ( "POrder" );
    scalNames.Push ( "NSteps" );
    scalNames.Push ( "Cohesion" );
    scalNames.Push ( "FrictionAngle" );
    scalNames.Push ( "FluxX" );
    scalNames.Push ( "FluxY" );
    vecNames.Push ( "Flux" );
    vecNames.Push ( "PrincipalStress" );
    scalNames.Push ( "Pressure" );


}

TPZCompMesh * CreateCompMeshKL(TPZGeoMesh * gmesh,int porder)
{

    int id=1;
    int dim = gmesh->Dimension();
    TPZCompMesh * cmesh = new TPZCompMesh ( gmesh );

    REAL lx=20.;
    REAL ly =2.;
    REAL lz=1.;

    int type=3;
    int expansionorder=150;

    KLMaterial * klmat = new KLMaterial ( id,lx,ly,lz,dim,type,expansionorder );

    cmesh->SetDefaultOrder ( porder );
    cmesh->SetDimModel ( dim );
    cmesh->InsertMaterialObject ( klmat );
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->AutoBuild();
    return cmesh;


}

void PrintMat ( std::string out,TPZFMatrix<REAL> mat )
{
    std::ofstream print ( out );
    int row=mat.Rows();
    int cols = mat.Cols();

    for ( int i=0; i<row; i++ ) {
        for ( int j=0; j<cols; j++ ) {
            print << mat ( i,j ) << " ";
        }
        print<< std::endl;
    }


}
REAL FindSol(TPZCompMesh *cmesh,TPZVec<REAL> &xd) {

    TPZGeoMesh * gmesh =  cmesh->Reference();
    int dim = gmesh->Dimension();
//    TPZVec<REAL> xd(dim,0.);
    TPZVec<REAL> mpt(dim,0.);

//     xd[0] =-0.5;
//     xd[1] = 0.5;
    int id;
    int nels = gmesh->NElements();
    for(int iel=0; iel<nels; iel++) {


        TPZVec<REAL> qsi(dim,0.);
        //TPZGeoEl * gel = gmesh->ElementVec()[iel];

        //TPZCompEl * cel = gel->Reference();

        TPZCompEl *cel = cmesh->Element ( iel );
        TPZGeoEl * gel = cel->Reference();
        TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *> ( cel);

        if(gel->MaterialId()<0)continue;
        bool check = gel->ComputeXInverse(xd,qsi,1.e-5);

        if(check==true)
        {


// 	   		cout << "elemento encontrado"<<endl;
// 			cout << "index do elemento = " <<gel->Index() <<endl;
            int nnodes=gel->NNodes();
            for(int inode=0; inode<nnodes; inode++)
            {
                TPZVec<REAL> co(3);
                TPZGeoNode* node = gel->NodePtr(inode);
                node->GetCoordinates(co);
                id = node->Id();

                if(fabs(co[0]-xd[0])<1.e-3 && fabs(co[1]-xd[1])<1.e-3)
                {
 					//cout << "node id = "<<id <<endl;
 					cout << " Coordinates "<< endl;
 					cout << " x = "<< co[0] << ",  y = " << co[1] << endl;
 					cout << " qsi = "<< qsi << endl;
                   //cel = gel->Reference();
                    //TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *> ( cel );

                    TPZMaterialData data;

                    data.fNeedsSol = true;

                    intel->InitMaterialData ( data );


                    intel->ComputeRequiredData ( data, qsi );

                    //~data.sol.Print(std::cout);
                    //uy
                    REAL uy=data.sol[0][1];// = cmesh->Solution()(2*id+1,0) ;

                    //cout << " sol = "<<data.sol[0] <<endl;
                    return uy;
                }
            }


        }
    }

}

void LoadingRamp ( TPZCompMesh * cmesh,  REAL factor,REAL gammasolo, REAL gammaagua)
{
    plasticmat * body= dynamic_cast<plasticmat *> ( cmesh->FindMaterial ( 1 ) );
    //plasticmatcrisfield * body= dynamic_cast<plasticmatcrisfield *> ( cmesh->FindMaterial ( 1 ) );
    TPZManVector<REAL, 3> force ( 3,0. );
    force[1]=(gammaagua-gammasolo);
    //force[1]=(-gammasolo);
    body->SetLoadFactor(factor);
    body->SetBodyForce ( force );

}
void  ReadFile ( std::string file,TPZFMatrix<REAL> &out )
{
    string line,line2, temp;

    ifstream myfile ( file );

    std::vector<std::vector<double>> coords;
    int counter=0;
    while ( getline ( myfile, line ) ) {
        std::vector<string> tokens;
        istringstream iss ( line );
        while ( iss >> temp ) tokens.push_back ( temp );
        std::vector<double> input_doub_temp= str_vec<double> ( tokens );
        coords.push_back ( input_doub_temp );
        counter++;
    }

    out.Resize ( coords.size(),coords[1].size() );
    for ( int iel=0; iel<coords.size(); iel++ ) {
        for ( int elnode=0; elnode<coords[iel].size(); elnode++ ) {
            out ( iel,elnode ) = coords[iel][elnode];
        }
    }
}

template <class T>
std::vector<T> str_vec ( std::vector<std::string> &vs )
{
    std::vector<T> ret;
    for ( std::vector<string>::iterator it = vs.begin(); it != vs.end(); ++it ) {
        istringstream iss ( *it );
        T temp;
        iss >> temp;
        ret.push_back ( temp );
    }
    return ret;
}


REAL GravityIncrease ( TPZCompMesh * cmesh )
{

    REAL FS=0.1,FSmax=1000.,FSmin=0.,tol=0.01;
    int neq = cmesh->NEquations();
    int maxcount=100;
    TPZFMatrix<REAL> displace ( neq,1 ),displace0 ( neq,1 );

    int counterout = 0;
    REAL factor =1.;
    REAL gammasolo=20.;
    REAL gammaagua=0.;
    LoadingRamp ( cmesh,factor  , gammasolo,  gammaagua);
    REAL norm = 1000.;
    REAL tol2 = 1.e-2;
    int NumIter = 100;
    bool linesearch = true;
    bool checkconv = false;

    do {

        std::cout << "FS = " << FS  <<" | Load step = " << counterout << " | Rhs norm = " << norm  << std::endl;
        LoadingRamp ( cmesh,FS  , gammasolo,  gammaagua);
		//SetSuportPressure(cmesh,FS);

        TPZElastoPlasticAnalysis  * anal = CreateAnalysis ( cmesh );
        chrono::steady_clock sc;
        auto start = sc.now();
        int iters;
        bool conv = anal->IterativeProcess ( cout, tol2, NumIter,linesearch,checkconv,iters);
		//bool conv =anal->IterativeProcess(cout, tol2, NumIter,linesearch,checkconv);
        auto end = sc.now();
        auto time_span = static_cast<chrono::duration<double>> ( end - start );
        cout << "| total time in iterative process =  " << time_span.count()<< std::endl;
        //anal->IterativeProcess ( outnewton, tol2, NumIter);

        norm = Norm ( anal->Rhs() );

        if ( conv==false) {
            cmesh->LoadSolution(displace0);
            //cmesh->Solution().Zero();
            FSmax = FS;
            FS = ( FSmin + FSmax ) / 2.;

        } else {
            // uy+=findnodalsol(cmesh);
            displace0 = anal->Solution();
            FSmin = FS;
            anal->AcceptSolution();
            FS = 1. / ( ( 1. / FSmin + 1. / FSmax ) / 2. );
        }
        // cout << "|asdadadasd =  " << std::endl;
        counterout++;

    }  while ( (( FSmax - FSmin ) / FS > tol && counterout<maxcount) );
    TPZElastoPlasticAnalysis  * anal = CreateAnalysis ( cmesh );
    anal->AcceptSolution();

    return FS;
}
void TransferSolution(TPZManVector<TPZCompMesh*,2> vecmesh,TPZCompMesh * cmesh,int isol)
{
    bool debug=false;
    TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> (cmesh->MaterialVec()[1]);
    int nels =  cmesh->NElements();

    TPZMaterial *mat1 = vecmesh[0]->FindMaterial(1);
    if (!mat1) {
        DebugStop();
    }
    TPZMaterial *mat2 = vecmesh[1]->FindMaterial(1);
    if (!mat2) {
        DebugStop();
    }

    TPZAdmChunkVector<TPZElastoPlasticMem>  &mem = pMatWithMem2->GetMemory();
    if (pMatWithMem2) {
        pMatWithMem2->SetUpdateMem(true);
    }
    else {
        DebugStop();
    }

    for(int iel=0;iel<nels;iel++)
    {

        TPZCompEl *cel = cmesh->ElementVec()[iel];
        if (!cel) {
            continue;
        }
        TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
        if (!intel) {
            continue;
        }
        if (intel->Material() != pMatWithMem2) {
            continue;
        }

    if (debug) {
        std::stringstream sout;
        TPZGeoEl *gel = cel->Reference();
        gel->Print(sout);
        REAL co[8][2] = {{-1,-1},{1,-1},{1,1},{-1,1},{0,-1},{1,0},{0,1},{-1,0}};
        for (int p = 0; p<8; p++) {
            TPZManVector<REAL,3> par(2,0.),x(3,0.);
            par[0] = co[p][0];
            par[1] = co[p][1];
            gel->X(par, x);
            cout << "point " << p << "co " << x << std::endl;
        }
    }

    TPZGeoMesh *gmesh1 = vecmesh[0]->Reference();
    TPZGeoMesh *gmesh2 = vecmesh[1]->Reference();


    const TPZIntPoints &intpoints = intel->GetIntegrationRule();
    int nint = intpoints.NPoints();
    TPZManVector<REAL,3> point(3,0.);
    TPZMaterialData data1,data2,data;
    intel->InitMaterialData(data);
    data.fNeedsSol = false;
    data1.fNeedsSol = true;
    data2.fNeedsSol = true;

    long elementid1 = 0;
    long elementid2 = 0;
    TPZManVector<REAL,3> qsi(2,0.);

    for (long ip =0; ip<nint; ip++) {
        REAL weight;
        intpoints.Point(ip, point, weight);
        data.intLocPtIndex = ip;
        intel->ComputeRequiredData(data, point);

        if(debug)
        {
            int memoryindex = data.intGlobPtIndex;
            std::stringstream sout;
            cout << "Local integration point index " << data.intLocPtIndex << std::endl;
            pMatWithMem2->PrintMem(sout,memoryindex);
        }

        TPZGeoEl *gel1 = gmesh1->FindElement(data.x, qsi, elementid1,2);
        if (!gel1) {
            DebugStop();
        }
        TPZGeoEl *gel2 = gmesh2->FindElement(data.x, qsi, elementid2,2);
        if (!gel2) {
            DebugStop();
        }
        TPZCompEl *cel1 = gel1->Reference();
        if (!cel1) {
            DebugStop();
        }
        TPZCompEl *cel2 = gel2->Reference();
        if (!cel2) {
            DebugStop();
        }

        TPZInterpolationSpace *intel1 = dynamic_cast<TPZInterpolationSpace *>(cel1);
        if (!intel1) {
            DebugStop();
        }

        TPZInterpolationSpace *intel2 = dynamic_cast<TPZInterpolationSpace *>(cel2);
        if (!intel2) {
            DebugStop();
        }

        data1.fNeedsSol = true;
        intel1->InitMaterialData(data1);
        data2.fNeedsSol = true;
        intel1->InitMaterialData(data2);

        intel1->ComputeRequiredData(data1, qsi);
        intel2->ComputeRequiredData(data2, qsi);

        if(debug)
        {
            REAL diff = dist(data1.x,data.x);
            if(diff > 1.e-6)
            {
                std::cout << "Point not found " << data2.x << std::endl;
                //DebugStop();
            }
        }

        int indexplastic =data.intGlobPtIndex;


        REAL coesao=data1.sol[isol][0];
        REAL atrito=data2.sol[isol][0];

//         TPZSolVec sol1,sol2;
//         ComputeSolution ( cel1, data1.phi,sol1);
//         ComputeSolution ( cel2, data2.phi,sol2);
//
//         REAL coesao=sol1[isol][0];
//         REAL atrito=sol2[isol][0];

        mem[indexplastic].fPlasticState.fmatprop.Resize ( 2 );
        mem[indexplastic].fPlasticState.fmatprop[0] = coesao;//coesao
        mem[indexplastic].fPlasticState.fmatprop[1] = atrito;//atritointerno

    }
    pMatWithMem2->SetUpdateMem(false);

    }

}
void ComputeSolution ( TPZCompEl *cel, TPZFMatrix<REAL> &phi,TPZSolVec &sol )
{

    const int numdof = cel->Material()->NStateVariables();


    int counter=0;

    std::set<long> cornerconnectlist;
    cel->BuildCornerConnectList(cornerconnectlist);
    int ncon = cornerconnectlist.size();

    TPZFMatrix<STATE> &MeshSol = cel->Mesh()->Solution();

    long numbersol = MeshSol.Cols();
    sol.Resize ( numbersol );
    for ( long is=0 ; is<numbersol; is++ ) {
        sol[is].Resize ( numdof );
        sol[is].Fill ( 0. );

    }
    TPZBlock<STATE> &block = cel->Mesh()->Block();
    long iv = 0, d;
    for ( int in=0; in<ncon; in++ ) {
        TPZConnect *df = &cel->Connect ( in );
        /** @brief Returns the Sequence number of the connect object */
        /** If the \f$ sequence number == -1 \f$ this means that the node is unused */
        long dfseq = df->SequenceNumber();

        /**
        * @brief Returns block dimension
        * @param block_diagonal Inquired block_diagonal
        */
        int dfvar = block.Size ( dfseq );
        /**
        * @brief Returns the position of first element block dependent on matrix diagonal
        * @param block_diagonal Inquired block_diagonal
        */
        long pos = block.Position ( dfseq );
        for ( int jn=0; jn<dfvar; jn++ ) {
            //cout << "\n pos+jn = " << pos+jn << endl;
            //cout << " ids[in] = " << ids[in] << endl;

            for ( int is=0; is<numbersol; is++ ) {

                sol[is][iv%numdof] += ( STATE ) phi ( iv/numdof,0 ) *MeshSol ( pos+jn,is );
            }
            iv++;
        }
    }

}//metho
