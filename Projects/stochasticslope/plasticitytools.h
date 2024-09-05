#pragma once

#include "pzcmesh.h"
#include "pzdarcymatwithmem.h"
#include "pzgmesh.h"
#include "TPZYCMohrCoulombPV.h"
#include "pzelastoplastic2D.h"
#include "pzelastoplasticmem.h"
#include "TPZPlasticStepPV.h"
#include "pzelastoplasticanalysis.h"
#include <pzskylstrmatrix.h> //symmetric skyline matrix storage
#include <pzstepsolver.h> //for TPZStepSolver
#include "pzfstrmatrix.h"
#include "pzbndcond.h"
#include "darcytools.h"

typedef TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> LEMC;

typedef   TPZMatElastoPlastic2D <LEMC, TPZElastoPlasticMem > plasticmat;

class PlasticityTools
{
protected:

    TPZGeoMesh * fgmesh;

    int fnumthreads;
    REAL fgammaagua;
    REAL fgammasolo;
    REAL fE;
    REAL fnu;
    REAL fc;
    REAL fphi;
    int fporder;
    REAL fntimes=8.;
    bool print;


    TPZVec<REAL> fPlasticDeformSqJ2;

    std::vector<double> fPlasticDeformSqJ22;

public:

    TPZCompMesh * fcmesh;

    PlasticityTools();
    PlasticityTools(TPZGeoMesh * gmesh,REAL E,REAL nu, REAL coes, REAL phi,REAL gammaagua, REAL gammasolo,int  numthreads,int porder);
    PlasticityTools(PlasticityTools & cp);

    ~PlasticityTools();


    TPZCompMesh * CreateCompMesh ( TPZGeoMesh * gmesh,int porder );

    TPZCompMesh * CreateCompMeshPressure ( TPZGeoMesh * gmesh,int porder );

    void LoadingRamp (REAL factor );

    TPZElastoPlasticAnalysis * CreateAnalysis ( );

    REAL GravityIncrease ( );

    REAL ShearRed ( );

    REAL ShearRed ( TPZManVector<TPZCompMesh*,3> vecmesh,int imc);

    REAL Solve(REAL tolfs,int numiterfs,REAL tolres,int numiterres, REAL l,REAL lambda0);

    REAL Solve(REAL tolfs,int numiterfs,REAL tolres,int numiterres, REAL l,REAL lambda0, bool & converge );

    void PostPlasticity(string vtkd);

    void CreatePostProcessingMesh (TPZPostProcAnalysis * PostProcess );

    void PostProcessVariables ( TPZStack<std::string> &scalNames, TPZStack<std::string> &vecNames );

    void TransferSolutionFrom ( TPZManVector<TPZCompMesh*,3> vecmesh,int isol );

    void TransferSolutionFromShearRed ( TPZManVector<TPZCompMesh*,3> vecmesh,int isol ,REAL FS );

    void TransferPermeability();

    void DivideElementsAbove(REAL sqj2, std::set<long> &elindices);

    void PRefineElementsAbove(REAL sqj2, int porder, std::set<long> &elindices);

    void ComputeElementDeformation();

    void ResetMem();

    void ApplyHistory();

    void SetPrint(bool set)
    {
        print=set;
    }

    void SetCP(REAL c, REAL p)
    {
        fc=c;
        fphi=p;
    }
};

PlasticityTools::PlasticityTools()
{

}

PlasticityTools::PlasticityTools(TPZGeoMesh * gmesh,REAL E,REAL nu, REAL coes, REAL phi,REAL gammaagua, REAL gammasolo,int  numthreads,int porder)
{
    fgmesh=gmesh;
    fE=E;
    fnu=nu;
    fc=coes;
    fphi=phi;
    //fcmesh=gmesh->Reference();
    fcmesh=CreateCompMesh(fgmesh,porder);
    fnumthreads=numthreads;
    fporder=porder;
    fgammaagua=gammaagua;
    fgammasolo=gammasolo;
    LoadingRamp(1.);
}

PlasticityTools::PlasticityTools(PlasticityTools & cp):fgmesh(cp.fgmesh),fcmesh(cp.fcmesh),fnumthreads(cp.fnumthreads),
fgammaagua(cp.fgammaagua),fgammasolo(cp.fgammasolo),fporder(cp.fporder)
{

}

PlasticityTools::~PlasticityTools()
{

}

void PlasticityTools::TransferSolutionFrom ( TPZManVector<TPZCompMesh*,3> vecmesh,int isol )
{
    bool debug=false;
    TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> ( fcmesh->MaterialVec() [1] );
    int nels =  fcmesh->NElements();

    TPZMaterial *mat1 = vecmesh[0]->FindMaterial ( 1 );
    if ( !mat1 )
    {
        DebugStop();
    }
    TPZMaterial *mat2 = vecmesh[1]->FindMaterial ( 1 );
    if ( !mat2 )
    {
        DebugStop();
    }

    TPZAdmChunkVector<TPZElastoPlasticMem>  &mem = pMatWithMem2->GetMemory();
    if ( pMatWithMem2 )
    {
        pMatWithMem2->SetUpdateMem ( true );
    }
    else
    {
        DebugStop();
    }

    for ( int iel=0; iel<nels; iel++ )
    {

        TPZCompEl *cel = fcmesh->ElementVec() [iel];
        if ( !cel )
        {
            continue;
        }
        TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *> ( cel );
        if ( !intel )
        {
            continue;
        }
        if ( intel->Material() != pMatWithMem2 )
        {
            continue;
        }

        if ( debug )
        {
            std::stringstream sout;
            TPZGeoEl *gel = cel->Reference();
            gel->Print ( sout );
            REAL co[8][2] = {{-1,-1},{1,-1},{1,1},{-1,1},{0,-1},{1,0},{0,1},{-1,0}};
            for ( int p = 0; p<8; p++ )
            {
                TPZManVector<REAL,3> par ( 2,0. ),x ( 3,0. );
                par[0] = co[p][0];
                par[1] = co[p][1];
                gel->X ( par, x );
                cout << "point " << p << "co " << x << std::endl;
            }
        }

        TPZGeoMesh *gmesh1 = vecmesh[0]->Reference();
        TPZGeoMesh *gmesh2 = vecmesh[1]->Reference();


        const TPZIntPoints &intpoints = intel->GetIntegrationRule();
        int nint = intpoints.NPoints();
        TPZManVector<REAL,3> point ( 3,0. );
        TPZMaterialData data1,data2,data;
        intel->InitMaterialData ( data );
        data.fNeedsSol = false;
        data1.fNeedsSol = true;
        data2.fNeedsSol = true;

        long elementid1 = 0;
        long elementid2 = 0;
        TPZManVector<REAL,3> qsi ( 2,0. );

        for ( long ip =0; ip<nint; ip++ )
        {
            REAL weight;
            intpoints.Point ( ip, point, weight );
            data.intLocPtIndex = ip;
            intel->ComputeRequiredData ( data, point );

            if ( debug )
            {
                int memoryindex = data.intGlobPtIndex;
                std::stringstream sout;
                cout << "Local integration point index " << data.intLocPtIndex << std::endl;
                pMatWithMem2->PrintMem ( sout,memoryindex );
            }

            TPZGeoEl *gel1 = gmesh1->FindElement ( data.x, qsi, elementid1,2 );
            if ( !gel1 )
            {
                DebugStop();
            }
            TPZGeoEl *gel2 = gmesh2->FindElement ( data.x, qsi, elementid2,2 );
            if ( !gel2 )
            {
                DebugStop();
            }
            TPZCompEl *cel1 = gel1->Reference();
            if ( !cel1 )
            {
                DebugStop();
            }
            TPZCompEl *cel2 = gel2->Reference();
            if ( !cel2 )
            {
                DebugStop();
            }

            TPZInterpolationSpace *intel1 = dynamic_cast<TPZInterpolationSpace *> ( cel1 );
            if ( !intel1 )
            {
                DebugStop();
            }

            TPZInterpolationSpace *intel2 = dynamic_cast<TPZInterpolationSpace *> ( cel2 );
            if ( !intel2 )
            {
                DebugStop();
            }

            data1.intLocPtIndex = ip;
            data2.intLocPtIndex = ip;
            data1.fNeedsSol = true;
            intel1->InitMaterialData ( data1 );
            data2.fNeedsSol = true;
            intel1->InitMaterialData ( data2 );

            intel1->ComputeRequiredData ( data1, qsi );
            intel2->ComputeRequiredData ( data2, qsi );

            if ( debug )
            {
                REAL diff = dist ( data1.x,data.x );
                if ( diff > 1.e-6 )
                {
                    std::cout << "Point not found " << data2.x << std::endl;
                    //DebugStop();
                }
            }

            int indexplastic =data.intGlobPtIndex;


            REAL coesao=data1.sol[isol][0];
            REAL atrito=data2.sol[isol][0];

            mem[indexplastic].fPlasticState.fmatprop.Resize ( 3 );
            mem[indexplastic].fPlasticState.fmatprop[0] = coesao;//coesao
            mem[indexplastic].fPlasticState.fmatprop[1] = atrito;//atritointerno

        }
        pMatWithMem2->SetUpdateMem ( false );

    }

}


TPZCompMesh * PlasticityTools::CreateCompMesh ( TPZGeoMesh * gmesh,int porder )
{
    unsigned int dim  = 2;
    const std::string name ( "ElastoPlastic COMP MESH Footing Problem " );

    TPZCompMesh * cmesh =  new TPZCompMesh ( gmesh );
    cmesh->SetName ( name );
    cmesh->SetDefaultOrder ( porder );
    cmesh->SetDimModel ( dim );

    // Mohr Coulomb data
    REAL mc_cohesion    = fc;//kpa
    REAL mc_phi         = fphi;
    REAL mc_psi         = mc_phi;

    /// ElastoPlastic Material using Mohr Coulomb
    // Elastic predictor
    TPZElasticResponse ER;
    REAL nu = fnu;
    REAL E = fE;

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

    material->SetLoadFactor ( 1. );
    material->SetWhichLoadVector ( 0 ); //option to compute the total internal force vecor fi=(Bt sigma+ N (b+gradu))

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



void PlasticityTools::LoadingRamp (  REAL factor )
{
    plasticmat * body= dynamic_cast<plasticmat *> ( fcmesh->FindMaterial ( 1 ) );
    //plasticmatcrisfield * body= dynamic_cast<plasticmatcrisfield *> ( cmesh->FindMaterial ( 1 ) );
    TPZManVector<REAL, 3> force ( 3,0. );
    force[1]= ( fgammaagua-fgammasolo );

    //force[1]=(-gammasolo);
    body->SetLoadFactor ( factor );
    body->SetBodyForce ( force );

}


TPZElastoPlasticAnalysis * PlasticityTools::CreateAnalysis (  )
{

    TPZElastoPlasticAnalysis * analysis =  new TPZElastoPlasticAnalysis ( fcmesh ); // Create analysis

 //   analysis->SetBiCGStab(10,1.e-6);

    TPZSkylineStructMatrix matskl ( fcmesh );
    //TPZFStructMatrix matskl ( cmesh );

    matskl.SetNumThreads ( fnumthreads );

  //  analysis->SetStructuralMatrix ( matskl );

    ///Setting a direct solver
    TPZStepSolver<STATE> step;
    step.SetDirect ( ELDLt );
    //step.SetDirect ( ECholesky );



    long neq = fcmesh->NEquations();
    TPZVec<long> activeEquations;
    analysis->GetActiveEquations(activeEquations);
    TPZEquationFilter filter(neq);
    filter.SetActiveEquations(activeEquations);
    matskl.EquationFilter() = filter;
    analysis->SetStructuralMatrix(matskl);

    analysis->SetSolver ( step );


    return analysis;

}

REAL PlasticityTools::Solve(REAL tolfs,int numiterfs,REAL tolres,int numiterres, REAL l,REAL lambda0 )
{

  //  fcmesh->Solution().Redim(fneq, 1);
    TPZElastoPlasticAnalysis * analysis = CreateAnalysis (  );
    bool converge;
    cout << "Solving incremental plasticity with Arc-length method..."<<endl;
    //REAL fs = analysis->IterativeProcessArcLength ( tolfs,numiterfs,tolres,numiterres,l,lambda0,converge );
    //analysis->IterativeProcessLinearisedArcLength ( );
    REAL ndesirediters=7;
     REAL llimit=l;
    REAL fs = analysis->IterativeProcessHybridArcLength ( tolfs,numiterfs,tolres,numiterres,l,lambda0,converge,ndesirediters,llimit );
   // REAL fs = analysis->IterativeProcessLinearisedArcLength ( tolfs,numiterfs,tolres,numiterres,l,lambda0,converge );
    //REAL fs = analysis->IterativeProcessArcLengthWithLineSearch ( tolfs,numiterfs,tolres,numiterres,l,lambda0,converge );
    if ( converge==false )
    {
        cout << "Arc-length fail to converge. Calling incremental plasticity brutal force."<<endl;
        fs = GravityIncrease (  );
    }
    delete analysis;
    return fs;
}

REAL PlasticityTools::Solve(REAL tolfs,int numiterfs,REAL tolres,int numiterres, REAL l,REAL lambda0, bool & converge )
{

  //  fcmesh->Solution().Redim(fneq, 1);
    TPZElastoPlasticAnalysis * analysis = CreateAnalysis (  );
    cout << "Solving incremental plasticity with Arc-length method..."<<endl;
    //REAL fs = analysis->IterativeProcessArcLength ( tolfs,numiterfs,tolres,numiterres,l,lambda0,converge );
    REAL ndesirediters=7;
     REAL llimit=0.04;
    REAL fs = analysis->IterativeProcessHybridArcLength ( tolfs,numiterfs,tolres,numiterres,l,lambda0,converge,ndesirediters,llimit );
    if ( converge==false )
    {
        cout << "Arc-length fail to converge. Calling incremental plasticity brutal force."<<endl;
        fs = GravityIncrease (  );
    }
    delete analysis;
    return fs;
}

REAL PlasticityTools::GravityIncrease (  )
{

  //  fcmesh->Solution().Redim(fneq, 1);
    REAL FS=0.1,FSmax=1000.,FSmin=0.,tol=0.01;
    int neq = fcmesh->NEquations();
    int maxcount=100;
    TPZFMatrix<REAL> displace ( neq,1 ),displace0 ( neq,1 );

    int counterout = 0;
    REAL factor =1.;
    LoadingRamp ( factor );
    REAL norm = 1000.;
    REAL tol2 = 1.e-2;
    int NumIter = 100;
    bool linesearch = true;
    bool checkconv = false;

    do
    {

        std::cout << "safety factor = " << FS  <<" | Load step = " << counterout << " | Rhs norm = " << norm  << std::endl;
        LoadingRamp ( FS );
        //SetSuportPressure(cmesh,FS);

        TPZElastoPlasticAnalysis  * anal = CreateAnalysis ( );
        chrono::steady_clock sc;
        auto start = sc.now();
        int iters;
        bool conv = anal->IterativeProcess ( cout, tol2, NumIter,linesearch,checkconv,iters );
        //bool conv =anal->IterativeProcess(cout, tol2, NumIter,linesearch,checkconv);
        auto end = sc.now();
        auto time_span = static_cast<chrono::duration<double>> ( end - start );
        //cout << "| total time in iterative process =  " << time_span.count() << std::endl;
        //anal->IterativeProcess ( outnewton, tol2, NumIter);

        norm = Norm ( anal->Rhs() );

        if ( conv==false )
        {
            fcmesh->LoadSolution ( displace0 );
            //cmesh->Solution().Zero();
            FSmax = FS;
            FS = ( FSmin + FSmax ) / 2.;

        }
        else
        {
            // uy+=findnodalsol(cmesh);
            displace0 = anal->Solution();
            FSmin = FS;
            anal->AcceptSolution();
            FS = 1. / ( ( 1. / FSmin + 1. / FSmax ) / 2. );
        }
        // cout << "|asdadadasd =  " << std::endl;
        counterout++;

    }
    while ( ( ( FSmax - FSmin ) / FS > tol && counterout<maxcount ) );
    TPZElastoPlasticAnalysis  * anal = CreateAnalysis (  );
    anal->AcceptSolution();

//    delete anal;
    return FS;
}

void PlasticityTools::ResetMem()
{
        TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> ( fcmesh->MaterialVec()[1]);
        pMatWithMem2->ResetMemory();
        fcmesh->Solution().Zero();

}

REAL PlasticityTools::ShearRed ( )
{

    //fgmesh->ResetReference();
   // fcmesh->LoadReferences();
    LoadingRamp(1.);

    REAL FS=0.8,FSmax=5.,FSmin=0.,tol=0.001;
    int neq = fcmesh->NEquations();

    TPZFMatrix<REAL> displace(neq,1),displace0(neq,1);

    int counterout = 0;

    plasticmat *material = dynamic_cast<plasticmat *> ( fcmesh->MaterialVec() [1] );
    TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> LEMC = material->GetPlasticity();
    TPZElasticResponse ER = LEMC.fER;
    REAL phi0 = LEMC.fYC.Phi();
    REAL cohesion0 = LEMC.fYC.Cohesion();
    REAL phi,psi,c;
    bool conv=false;
    do {

        fcmesh->Solution().Zero();
        std::cout << "FS "<< FS <<  "| step = " << counterout  <<std::endl;
        TPZElastoPlasticAnalysis  * anal = CreateAnalysis();

        REAL norm = 1000.;
        REAL tol2 = 0.01;
        int NumIter = 50;
        bool linesearch = true;
        bool checkconv = false;
        int iters;



        chrono::steady_clock sc;
        auto start = sc.now();

         conv = anal->IterativeProcess ( cout, tol2, NumIter,linesearch,checkconv,iters );
        norm = Norm ( anal->Rhs() );


        //anal->AcceptSolution();


        auto end = sc.now();
        auto time_span = static_cast<chrono::duration<double>> ( end - start );

        std::cout << "safety factor "<< FS <<  "| step = " << counterout <<" | Rhs norm = " << norm  << " | IterativeProcess Time: " << time_span.count() << " seconds !!! " <<std::endl;

        if ( conv==false ) {

            FSmax = FS;
            FS = ( FSmin + FSmax ) / 2.;
        } else {

            FSmin = FS;
			FS = 1. / ( ( 1. / FSmin + 1. / FSmax ) / 2. );
        }

        c=cohesion0/FS;
        std::cout << "coes "<<  c <<std::endl;
        phi=atan ( tan ( phi0 ) /FS );
        psi=phi;
        LEMC.fYC.SetUp ( phi, psi, c, ER );
        material->SetPlasticity ( LEMC );
        counterout++;
        if(( FSmax - FSmin ) / FS < tol)anal->AcceptSolution();
    }  while ( ( FSmax - FSmin ) / FS > tol || conv==false);

        std::cout << "final safety factor "<< FS <<std::endl;
        return ( FSmax + FSmin )/2;
}

REAL PlasticityTools::ShearRed (TPZManVector<TPZCompMesh*,3> vecmesh,int imc )
{

    std::cout << "aqqqqq "<<std::endl;
    LoadingRamp(1.);

    REAL FS=0.8,FSmax=5.,FSmin=0.,tol=0.001;
    int neq = fcmesh->NEquations();

    TPZFMatrix<REAL> displace(neq,1),displace0(neq,1);

    int counterout = 0;

    plasticmat *material = dynamic_cast<plasticmat *> ( fcmesh->MaterialVec() [1] );
    TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> LEMC = material->GetPlasticity();
    TPZElasticResponse ER = LEMC.fER;
    REAL phi0 = LEMC.fYC.Phi();
    REAL cohesion0 = LEMC.fYC.Cohesion();
    REAL phi,psi,c;
    bool conv=false;
    do {

        fcmesh->Solution().Zero();
        std::cout << "FS "<< FS <<  "| step = " << counterout  <<std::endl;
        TPZElastoPlasticAnalysis  * anal = CreateAnalysis();

        TransferSolutionFromShearRed (  vecmesh,imc, FS );
        REAL norm = 1000.;
        REAL tol2 = 0.01;
        int NumIter = 50;
        bool linesearch = true;
        bool checkconv = false;
        int iters;



        chrono::steady_clock sc;
        auto start = sc.now();

         conv = anal->IterativeProcess ( cout, tol2, NumIter,linesearch,checkconv,iters );
        norm = Norm ( anal->Rhs() );


        //anal->AcceptSolution();


        auto end = sc.now();
        auto time_span = static_cast<chrono::duration<double>> ( end - start );

        std::cout << "FS "<< FS <<  "| step = " << counterout <<" | Rhs norm = " << norm  << " | IterativeProcess Time: " << time_span.count() << " seconds !!! " <<std::endl;

        if ( conv==false ) {

            FSmax = FS;
            FS = ( FSmin + FSmax ) / 2.;
        } else {

            FSmin = FS;
			FS = 1. / ( ( 1. / FSmin + 1. / FSmax ) / 2. );
        }

        counterout++;
        if(( FSmax - FSmin ) / FS < tol)anal->AcceptSolution();
    }  while ( ( FSmax - FSmin ) / FS > tol || conv==false);

        return ( FSmax + FSmin )/2;
}

void PlasticityTools::PostPlasticity(string vtkd)
{
    TPZPostProcAnalysis * postprocdeter = new TPZPostProcAnalysis();
    CreatePostProcessingMesh ( postprocdeter );

    TPZVec<int> PostProcMatIds ( 1,1 );

    TPZStack<std::string> PostProcVars, scalNames, vecNames;

    PostProcessVariables ( scalNames, vecNames );

    //string vtkd = "postprocessdeter.vtk";
    postprocdeter->DefineGraphMesh ( 2,scalNames,vecNames,vtkd );

    postprocdeter->PostProcess ( 0 );

    delete postprocdeter;
}
void  PlasticityTools::CreatePostProcessingMesh (TPZPostProcAnalysis * PostProcess )
{
    if ( PostProcess->ReferenceCompMesh() != fcmesh )
    {

        PostProcess->SetCompMesh ( fcmesh );

        TPZVec<int> PostProcMatIds ( 1,1 );
        TPZStack<std::string> PostProcVars, scalNames, vecNames;
        PostProcessVariables ( scalNames, vecNames );

        for ( int i=0; i<scalNames.size(); i++ )
        {
            PostProcVars.Push ( scalNames[i] );
        }
        for ( int i=0; i<vecNames.size(); i++ )
        {
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

void PlasticityTools::PostProcessVariables ( TPZStack<std::string> &scalNames, TPZStack<std::string> &vecNames )
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
    scalNames.Push ( "EEigenFunc" );


}

/// Divide the element using the plastic deformation as threshold
void PlasticityTools::DivideElementsAbove(REAL sqj2, std::set<long> &elindices)
{
        double sum = std::accumulate(fPlasticDeformSqJ22.begin(), fPlasticDeformSqJ22.end(), 0.0);
    double mean = sum / fPlasticDeformSqJ22.size();
    cout << "mean = " <<mean << endl;
    fgmesh->ResetReference();
    fcmesh->LoadReferences();
    TPZManVector<REAL,3> findel(3,0.),qsi(2,0.);


    long nelem = fcmesh->NElements();
    for (long el=0; el<nelem; el++) {
        TPZCompEl *cel = fcmesh->ElementVec()[el];
        if (!cel) {
            continue;
        }

        TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
        if (!intel) {
            DebugStop();
        }
        TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> (cel->Material());
        if (!pMatWithMem2) {
            continue;
        }
        if (fcmesh->ElementSolution()(el,0) < mean*fntimes) {
            continue;
        }
        int porder = intel->GetPreferredOrder();
        TPZStack<long> subels;
        long index = cel->Index();

        intel->Divide(index, subels,0);
        for (int is=0; is<subels.size(); is++) {
            elindices.insert(subels[is]);
            TPZCompEl *subcel = fcmesh->ElementVec()[subels[is]];

            TPZInterpolationSpace *subintel = dynamic_cast<TPZInterpolationSpace *>(subcel);
            if (!subintel) {
                DebugStop();
            }
            subintel->SetPreferredOrder(porder);
        }
    }
    // divide elements with more than one level difference
    bool changed = true;
    while (changed) {
        changed = false;
        std::set<long> eltodivide;
        long nelem = fcmesh->NElements();
        for (long el=0; el<nelem; el++) {
            TPZCompEl *cel = fcmesh->ElementVec()[el];
            if (!cel) {
                continue;
            }
            TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
            if (!intel) {
                DebugStop();
            }
            TPZGeoEl *gel = cel->Reference();
            if (!gel) {
                DebugStop();
            }
            int ns = gel->NSides();
            for (int is=0; is<ns; is++) {
                TPZGeoElSide gelside(gel, is);
                if (gelside.Dimension() != 1) {
                    continue;
                }
                TPZCompElSide big = gelside.LowerLevelCompElementList2(1);
                if (!big) {
                    continue;
                }
                TPZGeoElSide geobig(big.Reference());
                // boundary elements will be refined by AdjustBoundaryElements
                if (geobig.Element()->Dimension() != 2) {
                    continue;
                }
                if (gel->Level()-geobig.Element()->Level() > 1) {
                    eltodivide.insert(big.Element()->Index());
                }
            }
        }
        std::set<long>::iterator it;
        for (it = eltodivide.begin(); it != eltodivide.end(); it++) {
            changed = true;
            long el = *it;
            TPZCompEl *cel = fcmesh->ElementVec()[el];
            if (!cel) {
                continue;
            }
            TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
            if (!intel) {
                DebugStop();
            }
            TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> (cel->Material());
            if (!pMatWithMem2) {
                continue;
            }

            int porder = intel->GetPreferredOrder();
            TPZStack<long> subels;
            long index = cel->Index();
            intel->Divide(index, subels,0);
            for (int is=0; is<subels.size(); is++) {
                elindices.insert(subels[is]);
                TPZCompEl *subcel = fcmesh->ElementVec()[subels[is]];
                TPZInterpolationSpace *subintel = dynamic_cast<TPZInterpolationSpace *>(subcel);
                if (!subintel) {
                    DebugStop();
                }
                subintel->SetPreferredOrder(porder);
            }
        }
    }

    //ApplyHistory(elindices);
    ComputeElementDeformation();
    fcmesh->AdjustBoundaryElements();
   // fcmesh->InitializeBlock();
    fcmesh->Solution().Zero();
   // fneq=fcmesh->NEquations();
    fcmesh->Solution().Resize(0, 0);
    fcmesh->Solution().Redim(fcmesh->NEquations(), 1);
   // fcmesh->LoadReferences();

}

// Get the vector of element plastic deformations
void PlasticityTools::ComputeElementDeformation()
{
    long nelem = fcmesh->NElements();
    fPlasticDeformSqJ2.resize(nelem);
    fPlasticDeformSqJ2.Fill(0.);
    fcmesh->ElementSolution().Redim(nelem, 1);
    TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> (fcmesh->MaterialVec()[1]);
    if (!pMatWithMem2) {
        fPlasticDeformSqJ2.Fill(0.);
    }
    else
    {
        for (long el = 0; el<nelem; el++) {
            TPZCompEl *cel = fcmesh->ElementVec()[el];
            fPlasticDeformSqJ2[el] = 0.;
            if (!cel) {
                continue;
            }
            TPZManVector<long> memindices;
            cel->GetMemoryIndices(memindices);
            int numind = memindices.size();
            REAL sqj2el = 0.;
            for (int ind=0; ind<numind; ind++)
            {
                int memoryindex = memindices[ind];
                if (memoryindex < 0) {
                    continue;
                }
                TPZElastoPlasticMem &mem = pMatWithMem2->MemItem(memindices[ind]);
                TPZTensor<REAL> &plastic = mem.fPlasticState.fEpsP;
                 REAL J2 = plastic.J2();
                 REAL sqj2 = sqrt(J2);
                //REAL val=mem.fPlasticState.fAlpha;
                sqj2el = max(sqj2,sqj2el);
            }
            fPlasticDeformSqJ2[el] = sqj2el;
        }
    }



    int sz=fPlasticDeformSqJ2.size();
    fPlasticDeformSqJ22.resize(sz);
    for(int i=0;i<sz;i++)fPlasticDeformSqJ22[i]=fPlasticDeformSqJ2[i];
    fcmesh->SetElementSolution(0, fPlasticDeformSqJ2);
}

/// Change the polynomial order of element using the plastic deformation as threshold
void PlasticityTools::PRefineElementsAbove(REAL sqj2, int porder, std::set<long> &elindices)
{

    double sum = std::accumulate(fPlasticDeformSqJ22.begin(), fPlasticDeformSqJ22.end(), 0.0);
    double mean = sum / fPlasticDeformSqJ22.size();
    cout << "mean = " <<mean << endl;
    fgmesh->ResetReference();
    fcmesh->LoadReferences();
    long nelem = fcmesh->NElements();
    for (long el=0; el<nelem; el++) {
        TPZCompEl *cel = fcmesh->ElementVec()[el];
        if (!cel) {
            continue;
        }
        TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
        if (!intel) {
            DebugStop();
        }
        TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> (cel->Material());
        if (!pMatWithMem2) {
            continue;
        }

        if (fcmesh->ElementSolution()(el,0) < mean*fntimes) {
            continue;
        }
        TPZStack<long> subels;
        long index = cel->Index();
        elindices.insert(index);
        intel->SetPreferredOrder(porder);
    }
}

void PlasticityTools::TransferSolutionFromShearRed ( TPZManVector<TPZCompMesh*,3> vecmesh,int isol,REAL FS )
{
    bool debug=false;
    TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> ( fcmesh->MaterialVec() [1] );
    int nels =  fcmesh->NElements();

    TPZMaterial *mat1 = vecmesh[0]->FindMaterial ( 1 );
    if ( !mat1 )
    {
        DebugStop();
    }
    TPZMaterial *mat2 = vecmesh[1]->FindMaterial ( 1 );
    if ( !mat2 )
    {
        DebugStop();
    }

    TPZAdmChunkVector<TPZElastoPlasticMem>  &mem = pMatWithMem2->GetMemory();
    if ( pMatWithMem2 )
    {
        pMatWithMem2->SetUpdateMem ( true );
    }
    else
    {
        DebugStop();
    }

    for ( int iel=0; iel<nels; iel++ )
    {

        TPZCompEl *cel = fcmesh->ElementVec() [iel];
        if ( !cel )
        {
            continue;
        }
        TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *> ( cel );
        if ( !intel )
        {
            continue;
        }
        if ( intel->Material() != pMatWithMem2 )
        {
            continue;
        }

        if ( debug )
        {
            std::stringstream sout;
            TPZGeoEl *gel = cel->Reference();
            gel->Print ( sout );
            REAL co[8][2] = {{-1,-1},{1,-1},{1,1},{-1,1},{0,-1},{1,0},{0,1},{-1,0}};
            for ( int p = 0; p<8; p++ )
            {
                TPZManVector<REAL,3> par ( 2,0. ),x ( 3,0. );
                par[0] = co[p][0];
                par[1] = co[p][1];
                gel->X ( par, x );
                cout << "point " << p << "co " << x << std::endl;
            }
        }

        TPZGeoMesh *gmesh1 = vecmesh[0]->Reference();
        TPZGeoMesh *gmesh2 = vecmesh[1]->Reference();


        const TPZIntPoints &intpoints = intel->GetIntegrationRule();
        int nint = intpoints.NPoints();
        TPZManVector<REAL,3> point ( 3,0. );
        TPZMaterialData data1,data2,data;
        intel->InitMaterialData ( data );
        data.fNeedsSol = false;
        data1.fNeedsSol = true;
        data2.fNeedsSol = true;

        long elementid1 = 0;
        long elementid2 = 0;
        TPZManVector<REAL,3> qsi ( 2,0. );

        for ( long ip =0; ip<nint; ip++ )
        {
            REAL weight;
            intpoints.Point ( ip, point, weight );
            data.intLocPtIndex = ip;
            intel->ComputeRequiredData ( data, point );

            if ( debug )
            {
                int memoryindex = data.intGlobPtIndex;
                std::stringstream sout;
                cout << "Local integration point index " << data.intLocPtIndex << std::endl;
                pMatWithMem2->PrintMem ( sout,memoryindex );
            }

            TPZGeoEl *gel1 = gmesh1->FindElement ( data.x, qsi, elementid1,2 );
            if ( !gel1 )
            {
                DebugStop();
            }
            TPZGeoEl *gel2 = gmesh2->FindElement ( data.x, qsi, elementid2,2 );
            if ( !gel2 )
            {
                DebugStop();
            }
            TPZCompEl *cel1 = gel1->Reference();
            if ( !cel1 )
            {
                DebugStop();
            }
            TPZCompEl *cel2 = gel2->Reference();
            if ( !cel2 )
            {
                DebugStop();
            }

            TPZInterpolationSpace *intel1 = dynamic_cast<TPZInterpolationSpace *> ( cel1 );
            if ( !intel1 )
            {
                DebugStop();
            }

            TPZInterpolationSpace *intel2 = dynamic_cast<TPZInterpolationSpace *> ( cel2 );
            if ( !intel2 )
            {
                DebugStop();
            }

            data1.intLocPtIndex = ip;
            data2.intLocPtIndex = ip;
            data1.fNeedsSol = true;
            intel1->InitMaterialData ( data1 );
            data2.fNeedsSol = true;
            intel1->InitMaterialData ( data2 );

            intel1->ComputeRequiredData ( data1, qsi );
            intel2->ComputeRequiredData ( data2, qsi );

            if ( debug )
            {
                REAL diff = dist ( data1.x,data.x );
                if ( diff > 1.e-6 )
                {
                    std::cout << "Point not found " << data2.x << std::endl;
                    //DebugStop();
                }
            }

            int indexplastic =data.intGlobPtIndex;


            REAL coesao=data1.sol[isol][0];
            REAL atrito=data2.sol[isol][0];

            mem[indexplastic].fPlasticState.fmatprop.Resize ( 3 );
            mem[indexplastic].fPlasticState.fmatprop[0]=coesao/FS;
            mem[indexplastic].fPlasticState.fmatprop[1]=atan ( tan ( atrito ) /FS );

        }
        pMatWithMem2->SetUpdateMem ( false );

    }

}


void PlasticityTools::ApplyHistory()
{
    TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> ( fcmesh->MaterialVec() [1] );
    int nels =  fcmesh->NElements();

    for ( int iel=0; iel<nels; iel++ )
    {

        TPZCompEl *cel = fcmesh->ElementVec() [iel];
        if ( !cel )
        {
            continue;
        }
        TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *> ( cel );
        if ( !intel )
        {
            continue;
        }
        if ( intel->Material() != pMatWithMem2 )
        {
            continue;
        }

        // Reset the memory of the integration points of the element
        TPZManVector<long> pointindices;
        cel->GetMemoryIndices(pointindices);
        long npoints = pointindices.size();
        for (long ip = 0; ip<npoints; ip++) {
            long ind = pointindices[ip];
            pMatWithMem2->ResetMemItem(ind);
        }

        }
        pMatWithMem2->SetUpdateMem ( false );

}

// void PlasticityTools::ApplyHistory()
// {
//
//
//
//     std::set<long>::iterator it;
//     for (it=elindices.begin(); it != elindices.end(); it++) {
//         long elindex = *it;
//         TPZCompEl *cel = fcmesh->ElementVec()[elindex];
//         TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> (cel->Material());
//         if (!pMatWithMem2) {
//             DebugStop();
//         }
//         // Reset the memory of the integration points of the element
//         TPZManVector<long> pointindices;
//         cel->GetMemoryIndices(pointindices);
//         long npoints = pointindices.size();
//         for (long ip = 0; ip<npoints; ip++) {
//             long ind = pointindices[ip];
//             pMatWithMem2->ResetMemItem(ind);
//         }
//
//
//
//     }
// }
