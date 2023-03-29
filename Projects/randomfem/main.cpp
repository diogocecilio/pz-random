
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

TPZGeoMesh * CreateGMeshQuadrado(int ref);

TPZCompMesh * CreateCompMeshElastMatWithMem(TPZGeoMesh * gmesh,int porder);

TPZCompMesh * CreateCompMeshKL(TPZGeoMesh * gmesh,int porder);

TPZAnalysis * CreateAnalysis(TPZCompMesh *cmesh);

void  CreatePostProcessingMesh ( TPZPostProcAnalysis * PostProcess,TPZCompMesh * cmesh );

void PostProcessVariables ( TPZStack<std::string> &scalNames, TPZStack<std::string> &vecNames );

int main()
{

    //TPZElasticityMaterialWithMem<REAL> * mat = new TPZElasticityMaterialWithMem< REAL >();

    int ref=2;
    int porder=1;


    TPZGeoMesh *gmesh = CreateGMeshQuadrado(ref);

    TPZCompMesh * cmesh = CreateCompMeshElastMatWithMem(gmesh,porder);

    TPZAnalysis *analysis =  CreateAnalysis(cmesh);
//
     analysis->Run();
     analysis->Solution().Print(cout);

     TPZPostProcAnalysis * postprocdeter = new TPZPostProcAnalysis();
     CreatePostProcessingMesh ( postprocdeter, cmesh );

    TPZVec<int> PostProcMatIds ( 1,1 );

    TPZStack<std::string> PostProcVars, scalNames, vecNames;

    PostProcessVariables ( scalNames, vecNames );

    string vtk = "postprocess.vtk";
    postprocdeter->DefineGraphMesh ( 2,scalNames,vecNames,vtk );

    postprocdeter->PostProcess ( 0 );
/*
    TPZVec<std::string> scalarVars ( 3 ),vectorVars ( 1 );
    vectorVars[0] = "Displacement";
    scalarVars[0] = "SigmaX";
    scalarVars[1] = "SigmaY";
    scalarVars[2] = "TauXY";

    analysis->DefineGraphMesh ( 2,scalarVars,vectorVars,vtk);
    analysis->PostProcess ( 0 );*/


//     int ref=2;
//     int porder=1;
//
//
//     TPZGeoMesh *gmesh = CreateGMeshQuadrado(ref);
//
//     TPZCompMesh * cmesh = CreateCompMeshKL(gmesh,porder);
//
//     KLAnalysis * klanal = new KLAnalysis ( cmesh );
//
//     KLMaterial *mat = dynamic_cast<KLMaterial*> (cmesh->MaterialVec()[1]);
//
//     klanal->SetExpansionOrder ( mat->GetExpansioOrder() );
//
//     klanal->Solve();
//
//     cmesh->LoadSolution(klanal->Solution());
//
//     string file1 ="solution.vtk";
//
//     klanal->Post ( file1,2,0 );

    return 0;
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
    REAL E=8000.;
    REAL nu=0.;
    REAL fx=0.;
    REAL fy=-1.;
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

    cmesh->InsertMaterialObject ( bc_bottom );

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
    vecNames.Push ( "Displacement" );

}

TPZCompMesh * CreateCompMeshKL(TPZGeoMesh * gmesh,int porder)
{

    int id=1;
    int dim = gmesh->Dimension();
    TPZCompMesh * cmesh = new TPZCompMesh ( gmesh );

    REAL lx=1.;
    REAL ly =2.;
    REAL lz=1.;

    int type=3;
    int expansionorder=10;

    KLMaterial * klmat = new KLMaterial ( id,lx,ly,lz,dim,type,expansionorder );

    cmesh->SetDefaultOrder ( porder );
    cmesh->SetDimModel ( dim );
    cmesh->InsertMaterialObject ( klmat );
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->AutoBuild();
    return cmesh;


}
