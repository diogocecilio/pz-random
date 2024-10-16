
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

REAL GravityIncrease ( TPZCompMesh * cmesh );
typedef   TPZMatElastoPlastic2D <LEMC, TPZElastoPlasticMem > plasticmat;

TPZGeoMesh * CreateGMeshQuadrado(int ref);

TPZGeoMesh *  CreateGMesh ( int ref );

TPZCompMesh * CreateCompMeshMatWithMem(TPZGeoMesh * gmesh,int porder);

TPZCompMesh * CreateCompMeshKL(TPZGeoMesh * gmesh,int porder);

TPZElastoPlasticAnalysis * CreateAnalysis(TPZCompMesh *cmesh);

void  CreatePostProcessingMesh ( TPZPostProcAnalysis * PostProcess,TPZCompMesh * cmesh );

void PostProcessVariables ( TPZStack<std::string> &scalNames, TPZStack<std::string> &vecNames );

void SetRandomField(TPZCompMesh * rcmesh1,TPZCompMesh * rcmesh2,TPZCompMesh * cmesh,int isol);

TPZCompMesh * ComputeField(TPZGeoMesh *gmesh,int porder);

TPZFMatrix<REAL> CreateLogNormalRandomField(TPZFMatrix<REAL> PHI, REAL mean, REAL cov,int samples,string outdata);

void PrintMat ( std::string out,TPZFMatrix<REAL> mat );

REAL SolveElastoplastic(TPZCompMesh * rcmesh1,TPZCompMesh * rcmesh2,int imc, string vtk);

REAL FindSol(TPZCompMesh *cmesh,TPZVec<REAL> &xd);

void MonteCarlo();

void LoadingRamp ( TPZCompMesh * cmesh,  REAL factor,REAL gammasolo, REAL gammaagua);

TPZCompMesh * ComputeDummyCMesh(TPZGeoMesh *gmesh,int porder);

void  ReadFile ( std::string file,TPZFMatrix<REAL> &out );

template <class T>
std::vector<T> str_vec ( std::vector<std::string> &vs );

void SolveDeterministic();

void ComputeSolution ( TPZCompEl *cel, TPZFMatrix<REAL> &phi,TPZSolVec &sol );

int main()
{
//     TPZGeoMesh *gmesh = CreateGMesh(0);
//
//     return 0;

//     set<int> a;
//     a.insert(1);
//     a.insert(2);
//     a.insert(1);
//     a.insert(3);
//
//     for(auto it = a.begin(); it!=a.end();it++)
//     {
//         cout<<endl<<*it;
//     }
//
//     return 0;

 //    SolveDeterministic();
//
 //    return 0;
//     int ref=3;
//     int porder=1;
//     TPZGeoMesh *gmesh = CreateGMeshQuadrado(ref);
//
//     TPZCompMesh * cmeshfield =  ComputeField(gmesh,porder);
//
//     SolveElastic(cmeshfield,0);

    MonteCarlo();

    return 0;
}

void SolveDeterministic()
{
    int ref=1;
    int porder=2;
    TPZGeoMesh *gmesh = CreateGMeshQuadrado(ref);

   // TPZGeoMesh *gmesh = CreateGMesh(ref);

    TPZCompMesh * cmesh = CreateCompMeshMatWithMem(gmesh,porder);

    REAL factor =1.;
    REAL gammasolo=20.;
    REAL gammaagua=0.;
    LoadingRamp ( cmesh,factor  , gammasolo,  gammaagua);

    TPZElastoPlasticAnalysis *analysis =  CreateAnalysis(cmesh);

//   GravityIncrease ( cmesh );

        REAL tolfs =0.01;
    int numiterfs =40;
    REAL tolres = 1.e-6;
    int numiterres =20;
    REAL l =0.1;
    REAL lambda0=1.;
    bool converge;
   REAL fs = analysis->IterativeProcessArcLength(tolfs,numiterfs,tolres,numiterres,l,lambda0,converge);

    TPZVec<REAL> xd(3);
    xd[0]=35.;
    xd[1]= 40.;
    REAL sol = fabs(FindSol(cmesh,xd));

    std::cout << "top displacement = " << sol << std::endl;

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
    int ref=1;
    int porder=2;
    TPZGeoMesh *gmesh = CreateGMeshQuadrado(ref);
//    TPZGeoMesh *gmesh = CreateGMesh(ref);

    TPZCompMesh * cmeshfield;
    TPZCompMesh*  cmeshfield2;
    if(false)
    {
         cmeshfield =  ComputeField(gmesh,porder);
         return;
    }else{
         cmeshfield =  ComputeDummyCMesh(gmesh,porder);
         cmeshfield2 =  ComputeDummyCMesh(gmesh,porder);
         string outco="/home/diogo/projects/pz-random-build-release/Projects/randomfemplasticity/coesao.txt";
         TPZFMatrix<REAL> readco;
         ReadFile ( outco,readco );
         cmeshfield->LoadSolution(readco);

         string outphi="/home/diogo/projects/pz-random-build-release/Projects/randomfemplasticity/atrito.txt";
         TPZFMatrix<REAL> readphi;
         ReadFile ( outphi,readphi );
         cmeshfield2->LoadSolution(readphi);
    }

   // cmeshfield->Solution().Print(std::cout);
    //return;
    int samples = cmeshfield->Solution().Cols();

    std::cout << "samples = " <<  samples << std::endl;

    ofstream out("fs2.dat");


    string vtk = "saidamontecarlo.vtk";
    for(int imc=0;imc<samples;imc++)
    {
        std::cout << "imc = " <<  imc << std::endl;
        REAL fs= SolveElastoplastic(cmeshfield,cmeshfield2,imc,vtk);
        out << fs << std::endl;
    }


}


REAL SolveElastoplastic(TPZCompMesh * rcmesh1,TPZCompMesh * rcmesh2,int imc,string vtk)
{

    TPZGeoMesh *gmesh = rcmesh1->Reference();

    int porder = rcmesh1->GetDefaultOrder();

    TPZCompMesh * cmesh = CreateCompMeshMatWithMem(gmesh,porder);

    REAL factor =1.;
    REAL gammasolo=20.;
    REAL gammaagua=0.;
    LoadingRamp ( cmesh,factor  , gammasolo,  gammaagua);

    SetRandomField(rcmesh1,rcmesh2,cmesh,imc);

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

//      TPZPostProcAnalysis * postprocdeter = new TPZPostProcAnalysis();
//      CreatePostProcessingMesh ( postprocdeter, cmesh );
//
//     TPZVec<int> PostProcMatIds ( 1,1 );
//
//     TPZStack<std::string> PostProcVars, scalNames, vecNames;
//
//     PostProcessVariables ( scalNames, vecNames );
//
//     postprocdeter->DefineGraphMesh ( 2,scalNames,vecNames,vtk );
//
//     postprocdeter->PostProcess ( 0 );

    return fs;
}

// REAL SolveElastoplastic(TPZCompMesh * rcmesh1,TPZCompMesh * rcmesh2,int imc)
// {
//
//     TPZGeoMesh *gmesh = rcmesh1->Reference();
//
//     int porder = rcmesh1->GetDefaultOrder();
//
//     TPZCompMesh * cmesh = CreateCompMeshMatWithMem(gmesh,porder);
//
//     REAL factor =1.;
//     REAL gammasolo=20.;
//     REAL gammaagua=0.;
//     LoadingRamp ( cmesh,factor  , gammasolo,  gammaagua);
//
//     SetRandomField(rcmesh1,rcmesh2,cmesh,imc);
//
//     TPZElastoPlasticAnalysis *analysis =  CreateAnalysis(cmesh);
//     REAL tolfs =0.1;
//     int numiterfs =40;
//     REAL tolres = 1.e-3;
//     int numiterres =20;
//     REAL l =0.1;
//     REAL lambda0=1.;
//     analysis->IterativeProcessArcLength(tolfs,numiterfs,tolres,numiterres,l,lambda0);
//
//     TPZVec<REAL> xd(3);
//     xd[0]=35.;
//     xd[1]= 40.;
//     REAL sol = fabs(FindSol(cmesh,xd));
//
//     std::cout << "top displacement = " << sol << std::endl;
//
//      TPZPostProcAnalysis * postprocdeter = new TPZPostProcAnalysis();
//      CreatePostProcessingMesh ( postprocdeter, cmesh );
//
//     TPZVec<int> PostProcMatIds ( 1,1 );
//
//     TPZStack<std::string> PostProcVars, scalNames, vecNames;
//
//     PostProcessVariables ( scalNames, vecNames );
//
//     string vtk = "postprocess.vtk";
//     postprocdeter->DefineGraphMesh ( 2,scalNames,vecNames,vtk );
//
//     postprocdeter->PostProcess ( 0 );
//
//     return sol;
// }


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

    string file1 ="coesao.vtk";

    klanal->Post ( file1,2,0 );

    return cmesh;
}

TPZCompMesh * ComputeDummyCMesh(TPZGeoMesh *gmesh,int porder)
{

    TPZCompMesh * cmesh = CreateCompMeshKL(gmesh,porder);

    KLAnalysis * klanal = new KLAnalysis ( cmesh );

    KLMaterial *mat = dynamic_cast<KLMaterial*> (cmesh->MaterialVec()[1]);


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

void SetRandomField(TPZCompMesh * rcmesh1,TPZCompMesh * rcmesh2,TPZCompMesh * cmesh,int isol)
{
    int nels = cmesh->NElements();
    int nels0 = rcmesh1->NElements();

    if ( nels!=nels0 ) {
        //cout << "nels "<< nels << " | nels0 = "<< nels0 << " nelsgmesh = " << nelsgmesh <<endl;
        //DebugStop();
    }

    //TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> LEMC;
    plasticmat*mohrmat = dynamic_cast<plasticmat *> ( cmesh->MaterialVec() [1] );
    //EAL atrito ;
    REAL atrito = mohrmat->GetPlasticity().fYC.Psi();

    TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> ( cmesh->MaterialVec() [1] );
    TPZAdmChunkVector<TPZElastoPlasticMem>  &mem = pMatWithMem2->GetMemory();

    if ( pMatWithMem2 ) {
        pMatWithMem2->SetUpdateMem ( true );
    }

    int globpt=0;
    for ( int iel=0; iel<nels0; iel++ ) {

        TPZCompEl *celplastic = cmesh->Element ( iel );
        TPZInterpolationSpace *intelplastic = dynamic_cast<TPZInterpolationSpace *> ( celplastic );

        if(celplastic->Material()->Id()<0)
        {
            cout <<"Acessando elemento de contorno na malha plástica"<<endl;
            cout <<"Significa transferir a slolução de forma errada. Elementos não correspontem, até porque a malha RF não tem elementos de contorno."<<endl;
            continue;
            //DebugStop();
        }

        TPZCompEl *celrandom1 = rcmesh1->Element ( iel );
        TPZCompEl *celrandom2 = rcmesh2->Element ( iel );
        TPZInterpolationSpace *intelrandom1 = dynamic_cast<TPZInterpolationSpace *> ( celrandom1 );
        TPZInterpolationSpace *intelrandom2 = dynamic_cast<TPZInterpolationSpace *> ( celrandom2 );

        const TPZIntPoints &intpoints = intelplastic->GetIntegrationRule();

        TPZManVector<REAL,3> point ( 3,0. );

        TPZMaterialData dataplastic,datarandom1,datarandom2;

        datarandom1.fNeedsSol = true;
        datarandom2.fNeedsSol=true;
        intelplastic->InitMaterialData ( dataplastic );
        intelrandom1->InitMaterialData ( datarandom1 );
        intelrandom2->InitMaterialData ( datarandom2 );

        REAL weight=0;
        int nint = intpoints.NPoints();

        for ( long ip =0; ip<nint; ip++ ) {

            intpoints.Point ( ip, point, weight );
            //intpoints1.Point ( ip, point1, weight );

            dataplastic.intLocPtIndex = ip;
            datarandom1.intLocPtIndex=ip;
            datarandom2.intLocPtIndex=ip;

            intelplastic->ComputeRequiredData ( dataplastic, point );
            intelrandom1->ComputeRequiredData ( datarandom1, point );
            intelrandom2->ComputeRequiredData ( datarandom2, point );

            TPZSolVec sol1,sol2;
            ComputeSolution ( celrandom1, datarandom1.phi,sol1);
            ComputeSolution ( celrandom2, datarandom2.phi,sol2);


            int indexplastic = dataplastic.intGlobPtIndex;

            int indexrandom1 = datarandom1.intGlobPtIndex;

            int indexrandom2 = datarandom1.intGlobPtIndex;

            if(indexplastic!=indexrandom1)
            {
//                 std::cout << "index das malhas sao diferentes"<<std::endl;
//                 std::cout << "indexplastic = "<< indexplastic << "indexrandom = "<< indexrandom  << std::endl;
            }

            REAL coesao=sol1[isol][0];
            REAL atrito=sol2[isol][0];
            mem[indexplastic].fPlasticState.fmatprop.Resize ( 2 );
            mem[indexplastic].fPlasticState.fmatprop[0] = coesao;//coesao
            mem[indexplastic].fPlasticState.fmatprop[1] = atrito;//atritointerno

            //dataplastic.intGlobPtIndex = globpt;
            //globpt++;

        }
    }
    pMatWithMem2->SetUpdateMem ( false );

    cmesh->Solution().Zero();
}


TPZGeoMesh * CreateGMeshQuadrado(int ref)
{

    TPZGeoMesh *gmesh  =  new TPZGeoMesh();

    gmesh->SetDimension ( 2 );

  //  string file ="/home/diogo/projects/pz-random/data/teste-coarse.msh";
 //   string file ="/home/diogo/projects/pz-random/data/teste.msh";
   // string file ="/home/diogo/projects/pz-random/data/n1052e1990.msh";
    //string file ="/home/diogo/projects/pz-random/data/n860e1647.msh";
    //string file ="/home/diogo/projects/pz-random/data/n804e1531.msh";//melhor-custo-beneficio
  //  string file ="/home/diogo/projects/pz-random/data/n2978e5671.msh";
    //string file ="/home/diogo/projects/pz-random/data/tri748.msh";
     string file ="/home/diogo/projects/pz-random/data/mesh.msh";
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

        set<int> SETmatRefDir11;
        //SETmatRefDir11.insert(2);
        SETmatRefDir11.insert(2);
        for(int j = 0; j <  4; j++)
		{
			long nel = gmesh->NElements();
			for (long iref = 0; iref < nel; iref++)
			{
				TPZVec<TPZGeoEl*> filhos;
				TPZGeoEl * gelP11 = gmesh->ElementVec()[iref];
                TPZRefPatternTools::RefineUniformIfNeighMat(gelP11, SETmatRefDir11);
            }
        }



   // gmesh->Print(std::cout);
    std::ofstream files ( "teste-mesh.vtk" );
    TPZVTKGeoMesh::PrintGMeshVTK ( gmesh,files,false );
cout << "d" << endl;
    return gmesh;


}


TPZCompMesh * CreateCompMeshMatWithMem(TPZGeoMesh * gmesh,int porder)
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


// TPZGeoMesh *  CreateGMesh ( int ref )
// {
//     const std::string name ( "Slope Problem " );
//
//     TPZGeoMesh *gmesh  =  new TPZGeoMesh();
//
//     gmesh->SetName ( name );
//     gmesh->SetDimension ( 2 );
//
//     TPZVec<REAL> coord ( 2 );
//
//     vector<vector<double>> co= {
//         {0.,0.},{75.,0.},{75.,30.},{45.,30.},{35.,40.},{0.,40.},{0,30},{35.,30.}
//     };
//     vector<vector<int>> meshtopol = {
//         {0,1,2,6},
//         {6,7,4,5},
//         {7,3,4}
//     };
//
//     gmesh->NodeVec().Resize ( co.size() );
//
//     for ( int inode=0; inode<co.size(); inode++ ) {
//         coord[0] = co[inode][0];
//         coord[1] = co[inode][1];
//         gmesh->NodeVec() [inode] = TPZGeoNode ( inode, coord, *gmesh );
//     }
//     TPZVec<long> TopoTri(3);
//     TPZVec<long> TopoQuad(4);
//     TPZVec<long> TopoLine(2);
//
//      for ( int iel=0; iel<meshtopol.size(); iel++ ) {
//             if ( meshtopol[iel].size() ==3 ) {
//                 TopoTri[0] =meshtopol[iel][0];
//                 TopoTri[1] =meshtopol[iel][1];
//                 TopoTri[2] =meshtopol[iel][2];
//                 new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> ( iel, TopoTri, 1,*gmesh );
//
//             }else if(meshtopol[iel].size() ==4){
//
//                 TopoQuad[0] =meshtopol[iel][0];
//                 TopoQuad[1] =meshtopol[iel][1];
//                 TopoQuad[2] =meshtopol[iel][2];
//                 TopoQuad[3] =meshtopol[iel][3];
//                 new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> ( iel, TopoQuad, 1,*gmesh );
//
//             }else{
//                 DebugStop();
//             }
//      }
//
//     int id=0;
//     id++;
//     TopoLine[0] = 0;
//     TopoLine[1] = 1;
//     new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -1, *gmesh );//bottom
//
//     id++;
//     TopoLine[0] = 1;
//     TopoLine[1] = 2;
//     new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -2, *gmesh );//rigth
//
//     id++;
//     TopoLine[0] = 2;
//     TopoLine[1] = 3;
//     new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -3, *gmesh );//top rigth
//
//     id++;
//     TopoLine[0] = 4;
//     TopoLine[1] = 5;
//     new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -4, *gmesh );//top left
//
//     id++;
//     TopoLine[0] = 5;
//     TopoLine[1] = 6;
//     new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -5, *gmesh );//left
//
//         id++;
//     TopoLine[0] = 6;
//     TopoLine[1] = 0;
//     new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -5, *gmesh );//left
//
//
//     id++;
//     TopoLine[0] = 3;
//     TopoLine[1] = 4;
//     new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -6, *gmesh );//ramp
//
//
//     TPZVec<long> TopolPoint(1);
//     id++;
//     TopolPoint[0] = 3;
//     new TPZGeoElRefPattern< pzgeom::TPZGeoPoint> ( id, TopolPoint, 4, *gmesh );
//
//
//     gmesh->BuildConnectivity();
//
//
//         for ( int d = 0; d<4; d++ ) {
//         int nel = gmesh->NElements();
//         TPZManVector<TPZGeoEl *> subels;
//         for ( int iel = 0; iel<nel; iel++ ) {
//             TPZGeoEl *gel = gmesh->ElementVec() [iel];
//             gel->Divide ( subels );
//         }
//     }
// // 		set<int> SETmatRefDir11;
// //         //SETmatRefDir11.insert(2);
// //         SETmatRefDir11.insert(4);
// //         for(int j = 0; j <  4; j++)
// // 		{
// // 			long nel = gmesh->NElements();
// // 			for (long iref = 0; iref < nel; iref++)
// // 			{
// // 				TPZVec<TPZGeoEl*> filhos;
// // 				TPZGeoEl * gelP11 = gmesh->ElementVec()[iref];
// //                 int indexel = gmesh->ElementIndex(gelP11);
// //                 if(gelP11->HasSubElement()|| indexel>5)
// //                 {
// //                     TPZRefPatternTools::RefineUniformIfNeighMat(gelP11, SETmatRefDir11);
// //                 }
// //             }
// //         }
//
//
//
//     std::ofstream files ( "gequadx.vtk" );
//     TPZVTKGeoMesh::PrintGMeshVTK ( gmesh,files,false );
//     return gmesh;
//
// }
// #include "pzgeopoint.h"
// TPZGeoMesh *  CreateGMesh ( int ref )
// {
//     const std::string name ( "Slope Problem " );
//
//     TPZGeoMesh *gmesh  =  new TPZGeoMesh();
//
//     gmesh->SetName ( name );
//     gmesh->SetDimension ( 2 );
//
//     TPZVec<REAL> coord ( 2 );
//
//     vector<vector<double>> co= {
//         {0.,0.},{35.,0.},{45.,0.},{75.,0.},{75.,30.},{45.,30.},{35.,40.},{0.,40.},{0,30},{35.,30.}
//     };
//     vector<vector<int>> meshtopol = {
//         {0,1,9,8},
//         {1,2,5,9},
//         {2,3,4,5},
//         {8,9,6,7},
//         {9,5,6}
//
//     };
//
//     gmesh->NodeVec().Resize ( co.size() );
//
//     for ( int inode=0; inode<co.size(); inode++ ) {
//         coord[0] = co[inode][0];
//         coord[1] = co[inode][1];
//         gmesh->NodeVec() [inode] = TPZGeoNode ( inode, coord, *gmesh );
//     }
//     TPZVec<long> TopoTri(3);
//     TPZVec<long> TopoQuad(4);
//     TPZVec<long> TopoLine(2);
//
//      for ( int iel=0; iel<meshtopol.size(); iel++ ) {
//             if ( meshtopol[iel].size() ==3 ) {
//                 TopoTri[0] =meshtopol[iel][0];
//                 TopoTri[1] =meshtopol[iel][1];
//                 TopoTri[2] =meshtopol[iel][2];
//                 new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> ( iel, TopoTri, 1,*gmesh );
//
//             }else if(meshtopol[iel].size() ==4){
//
//                 TopoQuad[0] =meshtopol[iel][0];
//                 TopoQuad[1] =meshtopol[iel][1];
//                 TopoQuad[2] =meshtopol[iel][2];
//                 TopoQuad[3] =meshtopol[iel][3];
//                 new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> ( iel, TopoQuad, 1,*gmesh );
//
//             }else{
//                 DebugStop();
//             }
//      }
//
//     int id=0;
//     id++;
//     TopoLine[0] = 0;
//     TopoLine[1] = 1;
//     new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -1, *gmesh );//bottom
//
//     id++;
//     TopoLine[0] = 1;
//     TopoLine[1] = 2;
//     new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -1, *gmesh );//bottom
//
//     id++;
//     TopoLine[0] = 2;
//     TopoLine[1] = 3;
//     new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -1, *gmesh );//bottom
//
//     id++;
//     TopoLine[0] = 3;
//     TopoLine[1] = 4;
//     new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -2, *gmesh );//rigth
//
//     id++;
//     TopoLine[0] = 4;
//     TopoLine[1] = 5;
//     new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -3, *gmesh );//top rigth
//
//     id++;
//     TopoLine[0] = 6;
//     TopoLine[1] = 7;
//     new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -4, *gmesh );//top left
//
//     id++;
//     TopoLine[0] = 7;
//     TopoLine[1] = 8;
//     new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -5, *gmesh );//left
//
//         id++;
//     TopoLine[0] = 8;
//     TopoLine[1] = 0;
//     new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -5, *gmesh );//left
//
//
//     id++;
//     TopoLine[0] = 5;
//     TopoLine[1] = 6;
//     new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -6, *gmesh );//ramp
//
//
//     TPZVec<long> TopolPoint(1);
//     id++;
//     TopolPoint[0] = 5;
//     new TPZGeoElRefPattern< pzgeom::TPZGeoPoint> ( id, TopolPoint, 4, *gmesh );
//
//
//     gmesh->BuildConnectivity();
//
//             for ( int d = 0; d<1; d++ ) {
//         int nel = gmesh->NElements();
//         TPZManVector<TPZGeoEl *> subels;
//         for ( int iel = 0; iel<nel; iel++ ) {
//             TPZGeoEl *gel = gmesh->ElementVec() [iel];
//             int indexel = gmesh->ElementIndex(gel);
//             if(gel->HasSubElement()||indexel==4|| indexel>4)
//             {
//                 gel->Divide ( subels );
//             }
//
//         }
//     }
//         for ( int d = 0; d<5; d++ ) {
//         int nel = gmesh->NElements();
//         TPZManVector<TPZGeoEl *> subels;
//         for ( int iel = 0; iel<nel; iel++ ) {
//             TPZGeoEl *gel = gmesh->ElementVec() [iel];
//             gel->Divide ( subels );
//         }
//     }
//
//
//
//     // //         for ( int d = 0; d<1; d++ ) {
// // //         int nel = gmesh->NElements();
// // //         TPZManVector<TPZGeoEl *> subels;
// // //         for ( int iel = 0; iel<nel; iel++ ) {
// // //             TPZGeoEl *gel = gmesh->ElementVec() [iel];
// // //             int indexel = gmesh->ElementIndex(gel);
// // //             if(gel->HasSubElement()||indexel==1||indexel==8||indexel==7|| indexel>11)
// // //             {
// // //                 gel->Divide ( subels );
// // //             }
// // //
// // //         }
// // //     }
//
// 		set<int> SETmatRefDir11;
//         //SETmatRefDir11.insert(2);
//         SETmatRefDir11.insert(4);
//         for(int j = 0; j <  6; j++)
// 		{
// 			long nel = gmesh->NElements();
// 			for (long iref = 0; iref < nel; iref++)
// 			{
// 				TPZVec<TPZGeoEl*> filhos;
// 				TPZGeoEl * gelP11 = gmesh->ElementVec()[iref];
//                 int indexel = gmesh->ElementIndex(gelP11);
//                 if(gelP11->HasSubElement()|| indexel>5)
//                 {
//                     TPZRefPatternTools::RefineUniformIfNeighMat(gelP11, SETmatRefDir11);
//                 }
//             }
//         }
//
//
//
//     std::ofstream files ( "gequadx.vtk" );
//     TPZVTKGeoMesh::PrintGMeshVTK ( gmesh,files,false );
//     return gmesh;
//
// }
// #include "pzgeopoint.h"
// TPZGeoMesh *  CreateGMesh ( int ref )
// {
//     const std::string name ( "Slope Problem " );
//
//     TPZGeoMesh *gmesh  =  new TPZGeoMesh();
//
//     gmesh->SetName ( name );
//     gmesh->SetDimension ( 2 );
//
//     TPZVec<REAL> coord ( 2 );
//
//     vector<vector<double>> co= {
//         {45.,30.},{35.,40.},{0.,40.},{0.,0.},{75.,0.},{75.,30.},{45.,28.},{30.,30.},{30.,40.},{28.,40.},
//         {28.,28.},{26.,40.},{26.,26.},{47.,26.},{47.,30.},{60.,30.},{60.,13.92857},{13.92857,13.92857},
//         {13.92857,40.}
//     };
//     vector<vector<int>> topol = {
//         {8,1,2,9},
//         {11,8,9,10},
//         {13,11,10,12},
//         {18,13,12,19},
//         {4,18,19,3},
//         {5,6,16,17},
//         {17,16,15,14},
//         {14,15,1,7},
//         {7,1,8,11},
//         {13,14,7,11},
//         {18,17,14,13},
//         {4,5,17,18}
//     };
//
//
//     gmesh->NodeVec().Resize ( co.size() );
//
//     for ( int inode=0; inode<co.size(); inode++ ) {
//         coord[0] = co[inode][0];
//         coord[1] = co[inode][1];
//         gmesh->NodeVec() [inode] = TPZGeoNode ( inode, coord, *gmesh );
//     }
//     TPZVec <long> TopoQuad ( 4 );
//     for ( int iel=0; iel<topol.size(); iel++ ) {
//         TopoQuad[0] = topol[iel][0]-1;
//         TopoQuad[1] = topol[iel][1]-1;
//         TopoQuad[2] =	topol[iel][2]-1;
//         TopoQuad[3] = topol[iel][3]-1;
//         new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> ( iel, TopoQuad, 1,*gmesh );
//     }
//     int id = topol.size();
//     TPZVec <long> TopoLine ( 2 );
//     TopoLine[0] = 3;
//     TopoLine[1] = 4;
//     new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -1, *gmesh );//bottom
//
//     id++;
//     TopoLine[0] = 4;
//     TopoLine[1] = 5;
//     new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -2, *gmesh );//rigth
//
//     {
//         id++;
//         TopoLine[0] = 5;
//         TopoLine[1] = 15;
//         new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -3, *gmesh ); // top rigth
//
//         id++;
//         TopoLine[0] = 15;
//         TopoLine[1] = 14;
//         new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -3, *gmesh ); // top rigth
//
//         id++;
//         TopoLine[0] = 14;
//         TopoLine[1] = 0;
//         new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -3, *gmesh ); // top rigth
//     }
//
//     id++;
//     TopoLine[0] = 0;
//     TopoLine[1] = 1;
//     new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -6, *gmesh );//ramp
//
//
//
//     {
//         id++;
//         TopoLine[0] = 1;
//         TopoLine[1] = 8;
//         new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -4, *gmesh ); //top left
//
//         id++;
//         TopoLine[0] = 8;
//         TopoLine[1] = 9;
//         new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -4, *gmesh ); //top left
//
//         id++;
//         TopoLine[0] = 9;
//         TopoLine[1] = 11;
//         new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -4, *gmesh ); //top left
//
//         id++;
//         TopoLine[0] = 11;
//         TopoLine[1] = 18;
//         new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -4, *gmesh ); //top left
//
//         id++;
//         TopoLine[0] = 18;
//         TopoLine[1] = 2;
//         new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -4, *gmesh ); //top left
//     }
//
//
//
//
//     id++;
//     TopoLine[0] = 2;
//     TopoLine[1] = 3;
//     new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -5, *gmesh );//left
//
//
//
//     TPZVec<long> TopolPoint(1);
//     id++;
//     TopolPoint[0] = 0;
//     new TPZGeoElRefPattern< pzgeom::TPZGeoPoint> ( id, TopolPoint, 4, *gmesh );
//
//
// //         id++;
// //         TopoLine[0] = 0;
// //         TopoLine[1] = 7;
// //         new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, 2, *gmesh );
//
//     gmesh->BuildConnectivity();
//
// //     for ( int d = 0; d<1; d++ ) {
// //         int nel = gmesh->NElements();
// //         TPZManVector<TPZGeoEl *> subels;
// //         for ( int iel = 0; iel<nel; iel++ ) {
// //             TPZGeoEl *gel = gmesh->ElementVec() [iel];
// //             int indexel = gmesh->ElementIndex(gel);
// //             if(gel->HasSubElement()||indexel==0 || indexel>11)
// //             {
// //                 gel->Divide ( subels );
// //             }
// //
// //         }
// //     }
// //
// //         for ( int d = 0; d<1; d++ ) {
// //         int nel = gmesh->NElements();
// //         TPZManVector<TPZGeoEl *> subels;
// //         for ( int iel = 0; iel<nel; iel++ ) {
// //             TPZGeoEl *gel = gmesh->ElementVec() [iel];
// //             int indexel = gmesh->ElementIndex(gel);
// //             if(gel->HasSubElement()||indexel==1||indexel==8||indexel==7|| indexel>11)
// //             {
// //                 gel->Divide ( subels );
// //             }
// //
// //         }
// //     }
// //
// //     for ( int d = 0; d<3; d++ ) {
// //         int nel = gmesh->NElements();
// //         TPZManVector<TPZGeoEl *> subels;
// //         for ( int iel = 0; iel<nel; iel++ ) {
// //             TPZGeoEl *gel = gmesh->ElementVec() [iel];
// //             gel->Divide ( subels );
// //         }
// //     }
// //     for ( int d = 0; d<3; d++ ) {
// //         int nel = gmesh->NElements();
// //         TPZManVector<TPZGeoEl *> subels;
// //         for ( int iel = 0; iel<nel; iel++ ) {
// //             TPZGeoEl *gel = gmesh->ElementVec() [iel];
// //             int indexel = gmesh->ElementIndex(gel);
// //             if(gel->HasSubElement()||indexel==0 ||indexel==8|| indexel>12 ||indexel==7)
// //             {
// //                 gel->Divide ( subels );
// //             }
// //
// //         }
// //     }
//         for ( int d = 0; d<3; d++ ) {
//         int nel = gmesh->NElements();
//         TPZManVector<TPZGeoEl *> subels;
//         for ( int iel = 0; iel<nel; iel++ ) {
//             TPZGeoEl *gel = gmesh->ElementVec() [iel];
//             gel->Divide ( subels );
//         }
//     }
// 		set<int> SETmatRefDir11;
//         //SETmatRefDir11.insert(2);
//         SETmatRefDir11.insert(4);
//         for(int j = 0; j < 5; j++)
// 		{
// 			long nel = gmesh->NElements();
// 			for (long iref = 0; iref < nel; iref++)
// 			{
// 				TPZVec<TPZGeoEl*> filhos;
// 				TPZGeoEl * gelP11 = gmesh->ElementVec()[iref];
//                 TPZRefPatternTools::RefineUniformIfNeighMat(gelP11, SETmatRefDir11);
// 			}
// 		}
//
//     std::ofstream files ( "gequad2.vtk" );
//     TPZVTKGeoMesh::PrintGMeshVTK ( gmesh,files,false );
//     return gmesh;
// }
// #include "pzgeopoint.h"
// TPZGeoMesh *  CreateGMesh ( int ref )
// {
//     const std::string name ( "Slope Problem " );
//
//     TPZGeoMesh *gmesh  =  new TPZGeoMesh();
//
//     gmesh->SetName ( name );
//     gmesh->SetDimension ( 2 );
//
//     TPZVec<REAL> coord ( 2 );
//
//     vector<vector<double>> co= {
//         {45.,30.},{35.,40.},{0.,40.},{0.,0.},{75.,0.},{75.,30.},
//         {45.,20.},//7
//         {30.,30.},{30.,40.},
//         {20.,40.},//10
//         {20.,20.},//11
//         {10.,40.},//12
//         {10.,10.},//13
//         {60.,10.},//14
//         {60.,30.}
//     };
//     vector<vector<int>> topol = {
//         {8,1,2,9},
//         {11,8,9,10},
//         {13,11,10,12},
//         {4,13,12,3},
//         {5,6,15,14},
//         {7,14,15,1},
//         {11,7,1,8},
//         {13,14,7,11},
//         {4,5,14,13}
//     };
//
//
//     gmesh->NodeVec().Resize ( co.size() );
//
//     for ( int inode=0; inode<co.size(); inode++ ) {
//         coord[0] = co[inode][0];
//         coord[1] = co[inode][1];
//         gmesh->NodeVec() [inode] = TPZGeoNode ( inode, coord, *gmesh );
//     }
//     TPZVec <long> TopoQuad ( 4 );
//     for ( int iel=0; iel<topol.size(); iel++ ) {
//         TopoQuad[0] = topol[iel][0]-1;
//         TopoQuad[1] = topol[iel][1]-1;
//         TopoQuad[2] =	topol[iel][2]-1;
//         TopoQuad[3] = topol[iel][3]-1;
//         new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> ( iel, TopoQuad, 1,*gmesh );
//     }
//     int id = topol.size();
//     TPZVec <long> TopoLine ( 2 );
//     TopoLine[0] = 3;
//     TopoLine[1] = 4;
//     new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -1, *gmesh );//bottom
//
//     id++;
//     TopoLine[0] = 4;
//     TopoLine[1] = 5;
//     new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -2, *gmesh );//rigth
//
//     {
//         id++;
//         TopoLine[0] = 5;
//         TopoLine[1] = 14;
//         new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -3, *gmesh ); // top rigth
//         id++;
//         TopoLine[0] = 14;
//         TopoLine[1] = 0;
//         new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -3, *gmesh ); // top rigth
//     }
//
//     id++;
//     TopoLine[0] = 0;
//     TopoLine[1] = 1;
//     new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -6, *gmesh );//ramp
//
//
//
//     {
//         id++;
//         TopoLine[0] = 1;
//         TopoLine[1] = 8;
//         new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -4, *gmesh ); //top left
//
//         id++;
//         TopoLine[0] = 8;
//         TopoLine[1] = 9;
//         new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -4, *gmesh ); //top left
//
//         id++;
//         TopoLine[0] = 9;
//         TopoLine[1] = 11;
//         new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -4, *gmesh ); //top left
//         id++;
//         TopoLine[0] = 11;
//         TopoLine[1] = 2;
//         new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -4, *gmesh ); //top left
//     }
//
//
//
//
//     id++;
//     TopoLine[0] = 2;
//     TopoLine[1] = 3;
//     new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -5, *gmesh );//left
//
//
//
//     TPZVec<long> TopolPoint(1);
//     id++;
//     TopolPoint[0] = 0;
//     new TPZGeoElRefPattern< pzgeom::TPZGeoPoint> ( id, TopolPoint, 4, *gmesh );
//
//
//     gmesh->BuildConnectivity();
//
//         for ( int d = 0; d<3; d++ ) {
//         int nel = gmesh->NElements();
//         TPZManVector<TPZGeoEl *> subels;
//         for ( int iel = 0; iel<nel; iel++ ) {
//             TPZGeoEl *gel = gmesh->ElementVec() [iel];
//             gel->Divide ( subels );
//         }
//     }
// 		set<int> SETmatRefDir11;
//         //SETmatRefDir11.insert(2);
//         SETmatRefDir11.insert(4);
//         for(int j = 0; j < 5; j++)
// 		{
// 			long nel = gmesh->NElements();
// 			for (long iref = 0; iref < nel; iref++)
// 			{
// 				TPZVec<TPZGeoEl*> filhos;
// 				TPZGeoEl * gelP11 = gmesh->ElementVec()[iref];
//                 TPZRefPatternTools::RefineUniformIfNeighMat(gelP11, SETmatRefDir11);
// 			}
// 		}
//
//     std::ofstream files ( "gequad2XX3.vtk" );
//     TPZVTKGeoMesh::PrintGMeshVTK ( gmesh,files,false );
//     return gmesh;
// }
// /*
// #include "pzgeopoint.h"
// TPZGeoMesh *  CreateGMesh ( int ref )
// {
//     const std::string name ( "Slope Problem " );
//
//     TPZGeoMesh *gmesh  =  new TPZGeoMesh();
//
//     gmesh->SetName ( name );
//     gmesh->SetDimension ( 2 );
//
//     TPZVec<REAL> coord ( 2 );
//
//     vector<vector<double>> co= {{0., 0.}, {75., 0.}, {75., 30.}, {45., 30.}, {35., 40.}, {0.,40.},
//
//         {35./3., 40.},{2 * 35/3., 40.}, {30., 40.},
//
//         {30., 30.}, {60.,30.},
//
//         {2* 35./3.,2* 35/3.}, {45., 2* 35/3.},
//
//         {35./3., 35/3.}, {60., 35./3.}
//
//     };
//     vector<vector<int>> topol = {
//         {0,  1,  14, 13},
//         {1,  2,  10, 14},
//         {14, 10, 3,  12},
//         {13, 14, 12, 11},
//         {11, 12, 3,  9},
//         {9,  3,  4,  8},
//         {11, 9,  8,  7},
//         {13, 11, 7, 6},
//         {0, 13,  6, 5}
//     };
//
//
//     gmesh->NodeVec().Resize ( co.size() );
//
//     for ( int inode=0; inode<co.size(); inode++ ) {
//         coord[0] = co[inode][0];
//         coord[1] = co[inode][1];
//         gmesh->NodeVec() [inode] = TPZGeoNode ( inode, coord, *gmesh );
//     }
//     TPZVec <long> TopoQuad ( 4 );
//     for ( int iel=0; iel<topol.size(); iel++ ) {
//         TopoQuad[0] = topol[iel][0];
//         TopoQuad[1] = topol[iel][1];
//         TopoQuad[2] =	topol[iel][2];
//         TopoQuad[3] = topol[iel][3];
//         new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> ( iel, TopoQuad, 1,*gmesh );
//     }
//     int id = topol.size();
//     TPZVec <long> TopoLine ( 2 );
//     TopoLine[0] = 0;
//     TopoLine[1] = 1;
//     new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -1, *gmesh );//bottom
//
//     id++;
//     TopoLine[0] = 1;
//     TopoLine[1] = 2;
//     new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -2, *gmesh );//rigth
//
//     {
//         id++;
//         TopoLine[0] = 3;
//         TopoLine[1] = 10;
//         new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -3, *gmesh ); // top rigth
//
//         id++;
//         TopoLine[0] = 10;
//         TopoLine[1] = 2;
//         new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -3, *gmesh ); // top rigth
//     }
//
//     id++;
//     TopoLine[0] = 3;
//     TopoLine[1] = 4;
//     new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -6, *gmesh );//ramp
//
//
//
//     {
//         id++;
//         TopoLine[0] = 5;
//         TopoLine[1] = 6;
//         new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -4, *gmesh ); //top left
//
//         id++;
//         TopoLine[0] = 6;
//         TopoLine[1] = 7;
//         new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -4, *gmesh ); //top left
//
//         id++;
//         TopoLine[0] = 7;
//         TopoLine[1] = 8;
//         new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -4, *gmesh ); //top left
//
//         id++;
//         TopoLine[0] = 8;
//         TopoLine[1] = 4;
//         new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -4, *gmesh ); //top left
//     }
//
//
//
//
//     id++;
//     TopoLine[0] = 0;
//     TopoLine[1] = 5;
//     new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, -5, *gmesh );//left
//
//
//     //refdir
//     id++;
//     TopoLine[0] = 3;
//     TopoLine[1] = 4;
//     new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, 2, *gmesh );
//
//     id++;
//     TopoLine[0] = 3;
//     TopoLine[1] = 9;
//     new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, 3, *gmesh );
//
//     TPZVec<long> TopolPoint(1);
//     id++;
//     TopolPoint[0] = 3;
//     new TPZGeoElRefPattern< pzgeom::TPZGeoPoint> ( id, TopolPoint, 4, *gmesh );
//
//     gmesh->BuildConnectivity();
//
//
//
// //     for ( int d = 0; d<1; d++ ) {
// //         int nel = gmesh->NElements();
// //         TPZManVector<TPZGeoEl *> subels;
// //         for ( int iel = 0; iel<nel; iel++ ) {
// //             TPZGeoEl *gel = gmesh->ElementVec() [iel];
// //             int indexel = gmesh->ElementIndex(gel);
// //             if(gel->HasSubElement()||indexel==5 ||indexel==8|| indexel>9)
// //             {
// //                 gel->Divide ( subels );
// //             }
// //
// //         }
// //     }
//
//     for ( int d = 0; d<3; d++ ) {
//         int nel = gmesh->NElements();
//         TPZManVector<TPZGeoEl *> subels;
//         for ( int iel = 0; iel<nel; iel++ ) {
//             TPZGeoEl *gel = gmesh->ElementVec() [iel];
//             gel->Divide ( subels );
//         }
//     }
//
//
//
// 		set<int> SETmatRefDir11;
//         SETmatRefDir11.insert(4);
//         for(int j = 0; j < 6; j++)
// 		{
// 			long nel = gmesh->NElements();
// 			for (long iref = 0; iref < nel; iref++)
// 			{
// 				TPZVec<TPZGeoEl*> filhos;
// 				TPZGeoEl * gelP11 = gmesh->ElementVec()[iref];
//                 TPZRefPatternTools::RefineUniformIfNeighMat(gelP11, SETmatRefDir11);
// 			}
// 		}
//
// /*
// 		set<int> setmat;
//         setmat.insert(3);
// 		for(int j = 0; j < 5; j++)
// 		{
// 			long nel = gmesh->NElements();
// 			for (long iref = 0; iref < nel; iref++)
// 			{
// 				TPZVec<TPZGeoEl*> filhos;
// 				TPZGeoEl * gelP11 = gmesh->ElementVec()[iref];
//                 TPZRefPatternTools::RefineDirectional(gelP11, setmat);
// 			}
// 		}*/
//
//
// // 				set<int> SETmatRefDir11;
// //         int matRefDir1=3;
// //         int RefDirId=6;
// // 		for(int j = 0; j < 2; j++)
// // 		{
// // 			long nel = gmesh->NElements();
// // 			for (long iref = 0; iref < nel; iref++)
// // 			{
// // 				TPZVec<TPZGeoEl*> filhos;
// // 				TPZGeoEl * gelP11 = gmesh->ElementVec()[iref];
// //                 TPZRefPatternTools::RefineUniformIfNeighMat(gelP11, SETmatRefDir11);
// // 				if(!gelP11) continue;
// // 				SETmatRefDir11.insert(matRefDir1);
// // 				if (RefDirId != 0) {
// // 					int matid= gelP11->MaterialId();
// // 					if(matid==1) TPZRefPatternTools::RefineDirectional(gelP11, SETmatRefDir11,RefDirId);
// // 				}
// // 				else TPZRefPatternTools::RefineDirectional(gelP11, SETmatRefDir11,RefDirId);
// // 			}
// // 		}
// //
// //         matRefDir1=3;
// //         RefDirId=6;
// // 		for(int j = 0; j < 2; j++)
// // 		{
// // 			long nel = gmesh->NElements();
// // 			for (long iref = 0; iref < nel; iref++)
// // 			{
// // 				TPZVec<TPZGeoEl*> filhos;
// // 				TPZGeoEl * gelP11 = gmesh->ElementVec()[iref];
// //                 TPZRefPatternTools::RefineUniformIfNeighMat(gelP11, SETmatRefDir11);
// // 				if(!gelP11) continue;
// // 				SETmatRefDir11.insert(matRefDir1);
// // 				if (RefDirId != 0) {
// // 					int matid= gelP11->MaterialId();
// // 					if(matid==1) TPZRefPatternTools::RefineDirectional(gelP11, SETmatRefDir11,RefDirId);
// // 				}
// // 				else TPZRefPatternTools::RefineDirectional(gelP11, SETmatRefDir11,RefDirId);
// // 			}
// // 		}
//
// 		        std::ofstream files ( "gequadXX.vtk" );
//     TPZVTKGeoMesh::PrintGMeshVTK ( gmesh,files,false );
//     return gmesh;
// }*/
// TPZGeoMesh * CreateGMeshQuadrado(int ref)
// {
//     TPZGeoMesh * gmesh = new TPZGeoMesh();
//
//     gmesh->SetDimension(2);
//
//     gmesh->NodeVec().Resize(4);
//
//     int id=0;
//     TPZVec<REAL> coord(2);
//     coord[0]=-0.5;coord[1]=-0.5;
//     gmesh->NodeVec() [id] = TPZGeoNode ( id, coord, *gmesh );
//
//
//     id++;
//     coord[0]=0.5;coord[1]=-0.5;
//     gmesh->NodeVec() [id] = TPZGeoNode ( id, coord, *gmesh );
//
//     id++;
//     coord[0]=0.5;coord[1]=0.5;
//     gmesh->NodeVec() [id] = TPZGeoNode ( id, coord, *gmesh );
//
//     id++;
//     coord[0]=-0.5;coord[1]=0.5;
//     gmesh->NodeVec() [id] = TPZGeoNode ( id, coord, *gmesh );
//
//     long iel=0;
//     TPZVec<long> topol(4);
//     topol[0]=0;
//     topol[1]=1;
//     topol[2]=2;
//     topol[3]=3;
//     new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> ( iel, topol, 1,*gmesh );
//
//
//
//     iel++;
//     TPZVec<long> topoline(2);
//     topoline[0] = 0;
//     topoline[1] = 1;
//     new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, topoline, -1, *gmesh );//bottom
//     topoline[0] = 1;
//     topoline[1] = 2;
//     new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, topoline, -2, *gmesh );//right
//     topoline[0] = 2;
//     topoline[1] = 3;
//     new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, topoline, -3, *gmesh );//top
//     topoline[0] = 3;
//     topoline[1] = 0;
//     new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, topoline, -4, *gmesh );//left
//
//     gmesh->BuildConnectivity();
//
//     for ( int d = 0; d<ref; d++ ) {
//     int nel = gmesh->NElements();
//         TPZManVector<TPZGeoEl *> subels;
//         for ( int iel = 0; iel<nel; iel++ ) {
//             TPZGeoEl *gel = gmesh->ElementVec() [iel];
//             gel->Divide ( subels );
//         }
//     }
//
//     std::ofstream files ( "ge.vtk" );
//     TPZVTKGeoMesh::PrintGMeshVTK ( gmesh,files,false );
//
//     return gmesh;
//
// }
//
//
// TPZCompMesh * CreateCompMeshMatWithMem(TPZGeoMesh * gmesh,int porder)
// {
//   unsigned int dim  = 2;
//     const std::string name ( "ElastoPlastic COMP MESH Footing Problem " );
//
//     TPZCompMesh * cmesh =  new TPZCompMesh ( gmesh );
//     cmesh->SetName ( name );
//     cmesh->SetDefaultOrder ( porder );
//     cmesh->SetDimModel ( dim );
//
//     // Mohr Coulomb data
// 	REAL mc_cohesion    = 10.;//kpa
//     REAL mc_phi         = 30.*M_PI/180.;
//     REAL mc_psi         = mc_phi;
//
//     /// ElastoPlastic Material using Mohr Coulomb
//     // Elastic predictor
//     TPZElasticResponse ER;
//     REAL nu = 0.49;
//     REAL E = 20000.;
//
//     LEMC plasticstep;
//     //TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> LEMC;
//     //TPZPlasticStepVoigt<TPZMohrCoulombVoigt, TPZElasticResponse> LEMC;
//     ER.SetUp ( E, nu );
//     plasticstep.fER =ER;
//     // LEMC.SetElasticResponse( ER );
//     plasticstep.fYC.SetUp ( mc_phi, mc_psi, mc_cohesion, ER );
//
//     int PlaneStrain = 1;
//
//     plasticmat * material = new plasticmat ( 1,PlaneStrain );
//     //plasticmatcrisfield * material = new plasticmatcrisfield ( 1,PlaneStrain );
//     material->SetPlasticity ( plasticstep );
//
//     material->SetId ( 1 );
//
//     material->SetLoadFactor(1.);
//     material->SetWhichLoadVector(0);//option to compute the total internal force vecor fi=(Bt sigma+ N (b+gradu))
//
//     cmesh->InsertMaterialObject ( material );
//
//     TPZFMatrix<STATE> val1 ( 2,2,0. );
//     TPZFMatrix<STATE>  val2 ( 2,1,0. );
//     int directionadirichlet =3;
//     int newman =1;
//     val2 ( 0,0 ) = 1;
//     val2 ( 1,0 ) = 1;
//     auto * bc_bottom = material->CreateBC ( material, -1,directionadirichlet, val1, val2 );//bottom
//
//
//     val2 ( 0,0 ) = 0;
//     val2 ( 1,0 ) = 0;
//     auto * bc_top = material->CreateBC ( material, -3,newman, val1, val2 );//top rigth
//
//
//
//     cmesh->InsertMaterialObject ( bc_bottom );
//     cmesh->InsertMaterialObject ( bc_top );
//
//
//     //cmesh->InsertMaterialObject ( top );
//     cmesh->SetAllCreateFunctionsContinuousWithMem();
//
//     cmesh->AutoBuild();
//
//     return cmesh;
//
// }

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
    int expansionorder=50;

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
