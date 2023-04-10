#pragma once

#include "pzcmesh.h"
#include "pzdarcymatwithmem.h"
#include "pzgmesh.h"


class DarcyTools
{
private:

    TPZGeoMesh * fgmesh;
    TPZCompMesh * fcmesh;
    REAL fH;//total slope heigth
    REAL fHw;//water heigth
    REAL fHt;//total heigth
    REAL fGammaW;
    int fporder;
    TPZAutoPointer<TPZFunction<STATE> > fForcingFunction;

public:
    DarcyTools();
    DarcyTools(DarcyTools & copy);
    DarcyTools(TPZGeoMesh* gmesh, REAL H,REAL Hw, REAL Ht,REAL gammaW,int fporder);
    ~DarcyTools();

   // static void Forcing(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);

    TPZCompMesh * CreateCMeshDarcy();
    void SolveDarcyProlem();
    void SetFlux ( TPZCompMesh * cmesh );

    void TransferSolutionFrom ( TPZManVector<TPZCompMesh*,3> vecmesh,int isol );

    void PostProcessVariables ( TPZStack<std::string> &scalNames, TPZStack<std::string> &vecNames );

    void  CreatePostProcessingMesh (TPZPostProcAnalysis * PostProcess );

    void PostDarcy(string vtk);

    void SetForcingFunction(TPZAutoPointer<TPZFunction<STATE> > fp);

    void (*fFunc)(const TPZVec<REAL> &x, TPZVec<REAL> &f);

};

DarcyTools::DarcyTools()
{

}
DarcyTools::~DarcyTools()
{

}

DarcyTools::DarcyTools(DarcyTools & copy):fgmesh(copy.fgmesh),fH(copy.fH),fHw(copy.fHw),fHt(copy.fHt),
fGammaW(copy.fGammaW),fporder(copy.fporder)
{

}

DarcyTools::DarcyTools(TPZGeoMesh* gmesh, REAL H,REAL Hw, REAL Ht,REAL gammaW,int porder)
{
        fgmesh=gmesh;
        fH=H;
        fHw=Hw;
        fHt=Ht;
        fGammaW=gammaW;
        fporder= porder;
        fcmesh = CreateCMeshDarcy( );
}

void DarcyTools::SetForcingFunction(TPZAutoPointer<TPZFunction<STATE> > fp)
{
        fForcingFunction = fp;
}

void Forcing(const TPZVec<REAL> &pt, TPZVec<STATE> &disp)
{
		const auto &x=pt[0];
        const auto &y=pt[1];
        const auto &z=pt[2];
        REAL yy=40-y;
        REAL hw=10.;
        REAL H = 10.;
        REAL gammaaugua=10.;

        if(yy<hw && yy<H)
        {
            disp[0]  =   gammaaugua * yy ;
        }
        else
        {
            disp[0]  = hw*gammaaugua;
        }
}

TPZCompMesh * DarcyTools::CreateCMeshDarcy( )
{
// void (MyClass::*func)(int);
// func = &MyClass::buttonClickedEvent;
 //   void ( DarcyTools::*ForcingBCPressao)(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
//    ForcingBCPressao=&DarcyTools::Forcing;

 //   void ( *ForcingBCPressao)(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
    //ForcingBCPressao=DarcyTools::Forcing;

    int dim=2;
    TPZCompMesh * cmesh = new TPZCompMesh ( fgmesh );
    cmesh->SetDefaultOrder ( fporder );
    cmesh->SetDimModel ( dim );

    //defaut permeability
    TPZDarcyMem mem;
    mem.fpermeability.Resize(3);
    mem.fpermeability[0]=1.;
    mem.fpermeability[1]=1.;
    mem.fpermeability[2]=1.;

    int matid=1;
    TPZDarcyMatWithMem<TPZDarcyMem> * material = new TPZDarcyMatWithMem<TPZDarcyMem> (matid,  dim);

    material->SetDefaultMem(mem);

	material->SetId(1);

    cmesh->InsertMaterialObject ( material );

    //Condition of contours
    TPZFMatrix<STATE>  val1 ( 2,2,0. );
    TPZFMatrix<STATE>  val2 ( 2,1,0. );

  //  this->fFunc(Forcing);

    TPZAutoPointer<TPZFunction<STATE> > pressure = new TPZDummyFunction<STATE>(Forcing);


	//TPZAutoPointer<TPZFunction<STATE> > pressure = new TPZDummyFunction<STATE>((this->*ForcingBCPressao)());

     TPZMaterial * BCond0 = material->CreateBC ( material, -3, 0, val1, val2 );//tr
     BCond0->SetForcingFunction(pressure);

     TPZMaterial * BCond1 = material->CreateBC ( material, -6, 0, val1, val2 );//ramp
     BCond1->SetForcingFunction(pressure);

	//val2(0,0)=-1;
	TPZMaterial * BCond2 = material->CreateBC ( material, -4, 0, val1, val2 );//tl
	BCond2->SetForcingFunction(pressure);


    val2(0,0)=0;
	TPZMaterial * BCond3 = material->CreateBC ( material, -5, 0, val1, val2 );//left
	BCond2->SetForcingFunction(pressure);

    val2(0,0)=10*10;
    TPZMaterial * BCond4 = material->CreateBC ( material, -2, 0, val1, val2 );//right
	BCond2->SetForcingFunction(pressure);

    cmesh->InsertMaterialObject(BCond0);
    cmesh->InsertMaterialObject(BCond1);
	cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    cmesh->InsertMaterialObject(BCond4);

    cmesh->SetAllCreateFunctionsContinuousWithMem();
    cmesh->AutoBuild();

    return cmesh;
}

void DarcyTools::SolveDarcyProlem()
{


    int numthreads = 0;
    TPZAnalysis *analysis =  new TPZAnalysis( fcmesh ); // Create analysis

    TPZSkylineStructMatrix matskl ( fcmesh );
    matskl.SetNumThreads ( numthreads );
    analysis->SetStructuralMatrix ( matskl );

    ///Setting a direct solver
    TPZStepSolver<STATE> step;
    step.SetDirect ( ELDLt );
    analysis->SetSolver ( step );


    analysis->Run();
}

void DarcyTools::SetFlux ( TPZCompMesh * cmesh )
{
    TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> ( cmesh->MaterialVec() [1] );
    int nels =  cmesh->NElements();

    TPZMaterial *mat1 = fcmesh->FindMaterial ( 1 );
    if ( !mat1 )
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

        TPZCompEl *cel = cmesh->ElementVec() [iel];
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

        TPZGeoMesh *gmesh1 = fcmesh->Reference();


        const TPZIntPoints &intpoints = intel->GetIntegrationRule();
        int nint = intpoints.NPoints();
        TPZManVector<REAL,3> point ( 3,0. );
        TPZMaterialData data1,data2,data;
        intel->InitMaterialData ( data );
        data.fNeedsSol = false;
        data1.fNeedsSol = true;
        data2.fNeedsSol = true;

        long elementid1 = 0;

        TPZManVector<REAL,3> qsi ( 2,0. );

        for ( long ip =0; ip<nint; ip++ )
        {
            REAL weight;
            intpoints.Point ( ip, point, weight );
            data.intLocPtIndex = ip;
            data1.intLocPtIndex = ip;
            intel->ComputeRequiredData ( data, point );


            TPZGeoEl *gel1 = gmesh1->FindElement ( data.x, qsi, elementid1,2 );
            if ( !gel1 )
            {
                DebugStop();
            }

            TPZCompEl *cel1 = gel1->Reference();
            if ( !cel1 )
            {
                DebugStop();
            }

            TPZInterpolationSpace *intel1 = dynamic_cast<TPZInterpolationSpace *> ( cel1 );
            if ( !intel1 )
            {
                DebugStop();
            }


            data1.fNeedsSol = true;
            intel1->InitMaterialData ( data1 );
            intel1->ComputeRequiredData ( data1, qsi );

            TPZMaterial *mat1= cel1->Material();

            TPZVec<REAL> flux,pressure;
            int varid =7;//flux
            mat1->Solution(data1,varid,flux);

            varid =1;//pressure
            mat1->Solution(data1,varid,pressure);

            int indexplastic =data.intGlobPtIndex;

            mem[indexplastic].fPlasticState.fflux.Resize ( 3 );
            mem[indexplastic].fPlasticState.fflux[0]=flux[0];
            mem[indexplastic].fPlasticState.fflux[1]=flux[1];
            mem[indexplastic].fPlasticState.fpressure=pressure[0];

        }
        pMatWithMem2->SetUpdateMem ( false );

    }

}

void DarcyTools::TransferSolutionFrom ( TPZManVector<TPZCompMesh*,3> vecmesh,int isol )
{

    TPZMatWithMem<TPZDarcyMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZDarcyMem> *> ( fcmesh->MaterialVec() [1] );
    int nels =  fcmesh->NElements();

    TPZMaterial *mat3 = vecmesh[2]->FindMaterial ( 1 );
    if ( !mat3 )
    {
        DebugStop();
    }

    TPZAdmChunkVector<TPZDarcyMem>  &mem = pMatWithMem2->GetMemory();
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

        TPZGeoMesh *gmesh3 = vecmesh[2]->Reference();


        const TPZIntPoints &intpoints = intel->GetIntegrationRule();
        int nint = intpoints.NPoints();
        TPZManVector<REAL,3> point ( 3,0. );
        TPZMaterialData data3,data;
        intel->InitMaterialData ( data );
        data.fNeedsSol = false;
        data3.fNeedsSol = true;

        long elementid1 = 0;
        long elementid3 = 0;
        TPZManVector<REAL,3> qsi ( 2,0. );

        for ( long ip =0; ip<nint; ip++ )
        {
            REAL weight;
            intpoints.Point ( ip, point, weight );
            data.intLocPtIndex = ip;
            data3.intLocPtIndex = ip;
            intel->ComputeRequiredData ( data, point );

            TPZGeoEl *gel3 = gmesh3->FindElement ( data.x, qsi, elementid1,2 );
            if ( !gel3 )
            {
                DebugStop();
            }
            TPZCompEl *cel3 = gel3->Reference();
            if ( !cel3 )
            {
                DebugStop();
            }

            TPZInterpolationSpace *intel3 = dynamic_cast<TPZInterpolationSpace *> ( cel3 );
            if ( !intel3 )
            {
                DebugStop();
            }


            data3.fNeedsSol = true;
            intel3->InitMaterialData ( data3 );
            intel3->ComputeRequiredData ( data3, qsi );

            int indexplastic =data.intGlobPtIndex;

            REAL perm=data3.sol[isol][0];

            mem[indexplastic].fpermeability.Resize ( 3 );
            mem[indexplastic].fpermeability[0] = perm;//coesao
            mem[indexplastic].fpermeability[1] = perm;//atritointerno
            mem[indexplastic].fpermeability[2] = perm;//atritointerno

        }
        pMatWithMem2->SetUpdateMem ( false );

    }

}

void DarcyTools::PostDarcy(string vtk)
{
    TPZPostProcAnalysis * postprocdeter = new TPZPostProcAnalysis();
    CreatePostProcessingMesh ( postprocdeter );

    TPZVec<int> PostProcMatIds ( 1,1 );

    TPZStack<std::string> PostProcVars, scalNames, vecNames;

    PostProcessVariables ( scalNames, vecNames );

    postprocdeter->DefineGraphMesh ( 2,scalNames,vecNames,vtk );

    postprocdeter->PostProcess ( 0 );
}
void  DarcyTools::CreatePostProcessingMesh (TPZPostProcAnalysis * PostProcess )
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

void DarcyTools::PostProcessVariables ( TPZStack<std::string> &scalNames, TPZStack<std::string> &vecNames )
{
    scalNames.Push ( "Perm" );
    scalNames.Push ( "Solution" );

    vecNames.Push ( "Flux" );
// if (!strcmp("Solution", name.c_str())) return 1;
//     if (!strcmp("Pressure", name.c_str())) return 1;
//     if (!strcmp("Derivative", name.c_str())) return 2;
//     if (!strcmp("GradU", name.c_str())) return 2;
//     if (!strcmp("KDuDx", name.c_str())) return 3;
//     if (!strcmp("KDuDy", name.c_str())) return 4;
//     if (!strcmp("KDuDz", name.c_str())) return 5;
//     if (!strcmp("NormKDu", name.c_str())) return 6;
//     if (!strcmp("MinusKGradU", name.c_str())) return 7;
//     if (!strcmp("Flux", name.c_str())) return 7;
//     if (!strcmp("POrder", name.c_str())) return 8;
//     if (!strcmp("ExactPressure", name.c_str())) return 9;
//     if (!strcmp("ExactSolution", name.c_str())) return 9;
//     if (!strcmp("ExactFlux", name.c_str())) return 10;
//     if (!strcmp("Div", name.c_str())) return 11;
//     if (!strcmp("Divergence", name.c_str())) return 11;
//     if (!strcmp("ExactDiv", name.c_str())) return 12;
//     if (!strcmp("ExactDivergence", name.c_str())) return 12;
//     if (!strcmp("FluxL2", name.c_str())) return 13;
//     if (!strcmp("Perm", name.c_str())) return 14;


}
