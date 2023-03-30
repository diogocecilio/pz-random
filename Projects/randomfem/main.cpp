
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

TPZGeoMesh * CreateGMeshQuadrado(int ref);

TPZCompMesh * CreateCompMeshElastMatWithMem(TPZGeoMesh * gmesh,int porder);

TPZCompMesh * CreateCompMeshKL(TPZGeoMesh * gmesh,int porder);

TPZAnalysis * CreateAnalysis(TPZCompMesh *cmesh);

void  CreatePostProcessingMesh ( TPZPostProcAnalysis * PostProcess,TPZCompMesh * cmesh );

void PostProcessVariables ( TPZStack<std::string> &scalNames, TPZStack<std::string> &vecNames );

void SetRandomField(TPZCompMesh * rcmesh,TPZCompMesh * cmesh,int isol);

TPZCompMesh * ComputeField(TPZGeoMesh *gmesh,int porder);

TPZFMatrix<REAL> CreateLogNormalRandomField(TPZFMatrix<REAL> PHI, REAL mean, REAL cov,int samples,string outdata);

void PrintMat ( std::string out,TPZFMatrix<REAL> mat );

void SolveElastic(TPZCompMesh * rcmesh);


int main()
{

    int ref=3;
    int porder=1;
    TPZGeoMesh *gmesh = CreateGMeshQuadrado(ref);

    TPZCompMesh * cmeshfield =  ComputeField(gmesh,porder);

    SolveElastic(cmeshfield);


    return 0;
}

void SolveElastic(TPZCompMesh * rcmesh)
{

    TPZGeoMesh *gmesh = rcmesh->Reference();

    int porder = rcmesh->GetDefaultOrder();

    TPZCompMesh * cmesh = CreateCompMeshElastMatWithMem(gmesh,porder);

    SetRandomField(rcmesh,cmesh,0);

    TPZAnalysis *analysis =  CreateAnalysis(cmesh);
     analysis->Run();
     //analysis->Solution().Print(cout);

     TPZPostProcAnalysis * postprocdeter = new TPZPostProcAnalysis();
     CreatePostProcessingMesh ( postprocdeter, cmesh );

    TPZVec<int> PostProcMatIds ( 1,1 );

    TPZStack<std::string> PostProcVars, scalNames, vecNames;

    PostProcessVariables ( scalNames, vecNames );

    string vtk = "postprocess.vtk";
    postprocdeter->DefineGraphMesh ( 2,scalNames,vecNames,vtk );

    postprocdeter->PostProcess ( 0 );
}


TPZCompMesh * ComputeField(TPZGeoMesh *gmesh,int porder)
{

    TPZCompMesh * cmesh = CreateCompMeshKL(gmesh,porder);

    KLAnalysis * klanal = new KLAnalysis ( cmesh );

    KLMaterial *mat = dynamic_cast<KLMaterial*> (cmesh->MaterialVec()[1]);

    klanal->SetExpansionOrder ( mat->GetExpansioOrder() );

    klanal->Solve();



    TPZFMatrix<REAL> eigenfunctions = klanal->Solution();
    REAL mean=1.;
    REAL cov=0.2;
    int samples=10;
    string outdata = "young.txt";

    TPZFMatrix<REAL> field = CreateLogNormalRandomField(eigenfunctions, mean, cov, samples, outdata);

    cmesh->LoadSolution(field);

    string file1 ="young.vtk";

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


void SetRandomField(TPZCompMesh * rcmesh,TPZCompMesh * cmesh,int isol)
{
    int nels = cmesh->NElements();
    int nels0 = rcmesh->NElements();
    int nelsgmesh = rcmesh->Reference()->NElements();

    if ( nels!=nels0 ) {
        cout << "nels "<< nels << " | nels0 = "<< nels0 << " nelsgmesh = " << nelsgmesh <<endl;
        //DebugStop();
    }

    TPZMatWithMem<TPZElasticMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElasticMem> *> ( cmesh->MaterialVec() [1] );
    TPZAdmChunkVector<TPZElasticMem>  &mem = pMatWithMem2->GetMemory();

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
            DebugStop();
        }

        TPZCompEl *celrandom1 = rcmesh->Element ( iel );
        TPZInterpolationSpace *intelrandom1 = dynamic_cast<TPZInterpolationSpace *> ( celrandom1 );

        const TPZIntPoints &intpoints = intelplastic->GetIntegrationRule();
        const TPZIntPoints &intpoints1 = intelrandom1->GetIntegrationRule();

        TPZManVector<REAL,3> point ( 3,0. );
        TPZManVector<REAL,3> point1 ( 3,0. );

        TPZMaterialData dataplastic,datarandom1,datarandom2;

        datarandom1.fNeedsSol = true;
        datarandom2.fNeedsSol = true;

        intelplastic->InitMaterialData ( dataplastic );
        intelrandom1->InitMaterialData ( datarandom1 );

        REAL weight=0;
        int nint = intpoints.NPoints();

        for ( long ip =0; ip<nint; ip++ ) {

            intpoints.Point ( ip, point, weight );
            //intpoints1.Point ( ip, point1, weight );

            dataplastic.intLocPtIndex = ip;
            datarandom1.intLocPtIndex=ip;

            intelplastic->ComputeRequiredData ( dataplastic, point );
            intelrandom1->ComputeRequiredData ( datarandom1, point );

            int indexplastic = dataplastic.intGlobPtIndex;

            int indexrandom = datarandom1.intGlobPtIndex;

            if(indexplastic!=indexrandom)
            {
//                 std::cout << "index das malhas sao diferentes"<<std::endl;
//                 std::cout << "indexplastic = "<< indexplastic << "indexrandom = "<< indexrandom  << std::endl;
            }

            mem[indexplastic].fE = datarandom1.sol[isol][0];

            //dataplastic.intGlobPtIndex = globpt;
            //globpt++;

        }
    }
    pMatWithMem2->SetUpdateMem ( false );

    cmesh->Solution().Zero();
}


TPZGeoMesh * CreateGMeshQuadrado(int ref)
{
    TPZGeoMesh * gmesh = new TPZGeoMesh();

    gmesh->SetDimension(2);

    gmesh->NodeVec().Resize(4);

    int id=0;
    TPZVec<REAL> coord(2);
    coord[0]=-0.5;coord[1]=-0.5;
    gmesh->NodeVec() [id] = TPZGeoNode ( id, coord, *gmesh );


    id++;
    coord[0]=0.5;coord[1]=-0.5;
    gmesh->NodeVec() [id] = TPZGeoNode ( id, coord, *gmesh );

    id++;
    coord[0]=0.5;coord[1]=0.5;
    gmesh->NodeVec() [id] = TPZGeoNode ( id, coord, *gmesh );

    id++;
    coord[0]=-0.5;coord[1]=0.5;
    gmesh->NodeVec() [id] = TPZGeoNode ( id, coord, *gmesh );

    long iel=0;
    TPZVec<long> topol(4);
    topol[0]=0;
    topol[1]=1;
    topol[2]=2;
    topol[3]=3;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> ( iel, topol, 1,*gmesh );



    iel++;
    TPZVec<long> topoline(2);
    topoline[0] = 0;
    topoline[1] = 1;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, topoline, -1, *gmesh );//bottom
    topoline[0] = 1;
    topoline[1] = 2;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, topoline, -2, *gmesh );//right
    topoline[0] = 2;
    topoline[1] = 3;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, topoline, -3, *gmesh );//top
    topoline[0] = 3;
    topoline[1] = 0;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( iel, topoline, -4, *gmesh );//left

    gmesh->BuildConnectivity();

    for ( int d = 0; d<ref; d++ ) {
    int nel = gmesh->NElements();
        TPZManVector<TPZGeoEl *> subels;
        for ( int iel = 0; iel<nel; iel++ ) {
            TPZGeoEl *gel = gmesh->ElementVec() [iel];
            gel->Divide ( subels );
        }
    }

    std::ofstream files ( "ge.vtk" );
    TPZVTKGeoMesh::PrintGMeshVTK ( gmesh,files,false );

    return gmesh;

}


TPZCompMesh * CreateCompMeshElastMatWithMem(TPZGeoMesh * gmesh,int porder)
{

    int dim = gmesh->Dimension();
    TPZCompMesh * cmesh = new TPZCompMesh ( gmesh );

    int id=1;
    REAL E=1.;
    REAL nu=0.;
    REAL fx=0.;
    REAL fy=0.;
    int planestress=1;

    TPZElasticMem mem;
    mem.fE=8000.;
    mem.fnu=0.;

    TPZElasticityMaterialWithMem<TPZElasticMem> * mat = new TPZElasticityMaterialWithMem<TPZElasticMem> (id,   fx,  fy,  planestress);


    mat->SetDefaultMem(mem);

    cmesh->InsertMaterialObject ( mat );

    TPZFMatrix<STATE> val1 ( 2,2,0. );
    TPZFMatrix<STATE>  val2 ( 2,1,0. );
	int dirichlet =0;
    val2(0,0)=1;
    val2(1,0)=1;
    auto * bc_bottom = mat->CreateBC ( mat, -1,3, val1, val2 );//bottom

    val2(0,0)=0;
    val2(1,0)=-1;
    auto * bc_load = mat->CreateBC ( mat, -3,1, val1, val2 );//bottom

    cmesh->InsertMaterialObject ( bc_bottom );

    cmesh->InsertMaterialObject ( bc_load );

    cmesh->SetDefaultOrder ( porder );
    cmesh->SetDimModel ( dim );
    cmesh->SetAllCreateFunctionsContinuousWithMem();
    cmesh->AutoBuild();

    return cmesh;

}

TPZAnalysis * CreateAnalysis(TPZCompMesh *cmesh)
{

    int numthreads =0;

    TPZAnalysis * analysis =  new TPZAnalysis( cmesh); // Create analysis

    TPZSkylineStructMatrix matskl ( cmesh );

    matskl.SetNumThreads ( numthreads );

    analysis->SetStructuralMatrix ( matskl );

    ///Setting a direct solver
    TPZStepSolver<STATE> step;
    step.SetDirect ( ELDLt );
    analysis->SetSolver ( step );

    return analysis;
}

#include "pzfstrmatrix.h"
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
    scalNames.Push ( "SigmaX" );
    scalNames.Push ( "youngmodulus" );
    vecNames.Push ( "Displacement" );

}

TPZCompMesh * CreateCompMeshKL(TPZGeoMesh * gmesh,int porder)
{

    int id=1;
    int dim = gmesh->Dimension();
    TPZCompMesh * cmesh = new TPZCompMesh ( gmesh );

    REAL lx=1.;
    REAL ly =1.;
    REAL lz=1.;

    int type=3;
    int expansionorder=4;

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
REAL findnodalsol(TPZCompMesh *cmesh) {

    TPZGeoMesh * gmesh =  cmesh->Reference();
    int dim = gmesh->Dimension();
    TPZVec<REAL> xd(dim,0.);
    TPZVec<REAL> mpt(dim,0.);

    xd[0] =-0.5;
    xd[1] = 0.5;
    int id;
    int nels = gmesh->NElements();
    for(int iel=0; iel<nels; iel++) {


        TPZVec<REAL> qsi(dim,0.);
        TPZGeoEl * gel = gmesh->ElementVec()[iel];

        TPZCompEl * cel = gel->Reference();

        if(gel->MaterialId()<0)continue;
        bool check = gel->ComputeXInverse(xd,qsi,1.e-5);

        if(check==true)
        {


// 	   		cout << "elemento encontrado"<<endl;
// 			cout << "index do elemento = " <<gel->Index() <<endl;
            int nnodes=gel->NNodes();
            for(int inode=0; inode<nnodes; inode++)
            {
                TPZVec<REAL> co(2);
                TPZGeoNode* node = gel->NodePtr(inode);
                node->GetCoordinates(co);
                id = node->Id();

                if(fabs(co[0]-xd[0])<1.e-3 && fabs(co[1]-xd[1])<1.e-3)
                {
 					cout << "node id = "<<id <<endl;
 					cout << " Coordinates "<< endl;
 					cout << " x = "<< co[0] << ",  y = " << co[1] << endl;
 					cout << " qsi = "<< qsi << endl;
                    TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *> ( cel );

                    TPZMaterialData data;

                    data.fNeedsSol = true;

                    intel->InitMaterialData ( data );

                    intel->ComputeRequiredData ( data, qsi );

                    return data.sol[0][1];
                }
            }


        }
    }

}
