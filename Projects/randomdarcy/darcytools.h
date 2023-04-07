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

public:
    DarcyTools();
    DarcyTools(DarcyTools & copy);
    DarcyTools(TPZGeoMesh* gmesh, REAL H,REAL Hw, REAL Ht,REAL gammaW,int fporder);
    ~DarcyTools();

   // static void Forcing(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);

    TPZCompMesh * CreateCMeshDarcy();
    void SolveDarcyProlem();
    void SetFlux ( TPZCompMesh * cmesh );


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
