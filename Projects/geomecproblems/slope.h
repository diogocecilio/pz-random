#include "pzskylstrmatrix.h"
#include "pzfstrmatrix.h"
#include "pzstepsolver.h"

#include "pzmaterial.h"
#include "pzbndcond.h"

#include "TPZVTKGeoMesh.h"
#include "pzgeoelbc.h"

#include "TPZGeoLinear.h"
#include "pzgeopoint.h"
#include "tpzgeoelrefpattern.h"


#include "TPZYCMohrCoulombPV.h"
#include "pzelastoplasticmem.h"
#include "pzelastoplastic2D.h"
#include "TPZPlasticStepPV.h"
#include "pzelastoplasticanalysis.h"

typedef TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> LEMC;

typedef   TPZMatElastoPlastic2D <LEMC, TPZElastoPlasticMem > plasticmat;

TPZGeoMesh *CreateGeoMesh();
TPZCompMesh *CreateCompMesh(TPZGeoMesh *gmesh);
void PostProcessVariables(TPZStack<std::string> &scalNames, TPZStack<std::string> &vecNames);
void  CreatePostProcessingMesh (TPZPostProcAnalysis * PostProcess ,TPZCompMesh * cmesh);
void PostPlasticity(string vtkd,TPZCompMesh * cmesh);
void findid(TPZCompMesh *cmesh,TPZVec<REAL> coord);
void runslope();

void runslope(){


    TPZGeoMesh *gmesh =  CreateGeoMesh();

    TPZCompMesh*cmesh = CreateCompMesh(gmesh);
    TPZVec<REAL> coord(3,0.);
    coord[0]=30.;
    coord[1]=30.;
   // findid(cmesh,coord);

    TPZElastoPlasticAnalysis * analysis = new TPZElastoPlasticAnalysis(cmesh);


    TPZSkylineStructMatrix matskl ( cmesh );

    TPZStepSolver<STATE> step;

    step.SetDirect ( ELDLt );

    analysis->SetStructuralMatrix(matskl);

    analysis->SetSolver ( step );

    REAL tolfs=0.001;
    REAL numiterfs=50;
    REAL numiterres =15;
    REAL tolres =1.e-3;
    REAL l =0.5;
    REAL lambda0=0.1;
    bool converge;

    REAL fs = analysis->IterativeProcessHybridArcLength ( tolfs,numiterfs,tolres,numiterres,l,lambda0,converge );
  //  REAL fs = analysis->IterativeProcessArcLength ( tolfs,numiterfs,tolres,numiterres,l,lambda0,converge );
    string file  = "slope.vtk";
    PostPlasticity(file,cmesh);


}

TPZGeoMesh *CreateGeoMesh() {

    REAL co[11][2] = {
    {0 ,0 },{30, 0},{60,0},
    {60,10},{30,10},{0,10},
    {60,20},{30,20},{0,20},
    {30,30},{0,30}
    };
    long indices[5][4] = {{0,1,4,5},{1,2,3,4},{4,3,6,7},{5,4,7,8},{8,7,9,10}};

    TPZGeoMesh *gmesh = new TPZGeoMesh();


    int dim =2;
    gmesh->SetDimension(dim);

    long nnode = 11;
    long nod;
    for(nod=0; nod<nnode; nod++) {
        long nodind = gmesh->NodeVec().AllocateNewElement();
        TPZVec<REAL> coord(2);
        coord[0] = co[nod][0];
        coord[1] = co[nod][1];
        gmesh->NodeVec()[nodind] = TPZGeoNode(nod,coord,*gmesh);
    }



    long el;
    long nelem = 5;
    for(el=0; el<nelem; el++) {
        TPZVec<long> nodind(4);
        for(nod=0; nod<4; nod++) nodind[nod]=indices[el][nod];
        new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> ( el, nodind, 1,*gmesh );
    }

    TPZVec <long> TopoLine ( 2 );
    long ind;

    TopoLine[0] = 10;
    TopoLine[1] = 8;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( ind, TopoLine, - 1, *gmesh );//left

    TopoLine[0] = 8;
    TopoLine[1] = 5;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( ind, TopoLine, - 1, *gmesh );//left

    TopoLine[0] = 5;
    TopoLine[1] = 0;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( ind, TopoLine, - 1, *gmesh );//left

    TopoLine[0] = 0;
    TopoLine[1] = 1;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( ind, TopoLine, - 2, *gmesh );//bottom

    TopoLine[0] = 1;
    TopoLine[1] = 2;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( ind, TopoLine, - 2, *gmesh );//bottom

    TopoLine[0] = 2;
    TopoLine[1] = 3;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( ind, TopoLine, - 3, *gmesh );//rigth


    TopoLine[0] = 3;
    TopoLine[1] = 6;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( ind, TopoLine, - 3, *gmesh );//rigth



gmesh->BuildConnectivity();

    for(int d=0;d<3;d++) {
    int nel = gmesh->NElements();
    for (int iel=0; iel<nel; iel++) {
        TPZManVector<TPZGeoEl *> subels;
        TPZGeoEl *gel = gmesh->ElementVec()[iel];
        gel->Divide(subels);
    }
	}

    gmesh->BuildConnectivity();
    std::ofstream print("gmesh.txt");
	gmesh->Print(print);
	 std::ofstream files ( "mesh-slope.vtk" );
    TPZVTKGeoMesh::PrintGMeshVTK ( gmesh,files,false );


	return gmesh;
}


TPZCompMesh *CreateCompMesh(TPZGeoMesh *gmesh) {


    unsigned int dim  = 2;
    const std::string name ( "ElastoPlastic COMP MESH Footing Problem " );

    TPZCompMesh * cmesh =  new TPZCompMesh ( gmesh );
    cmesh->SetName ( name );
    cmesh->SetDefaultOrder ( 3 );
    cmesh->SetDimModel ( dim );



      // Mohr Coulomb data
    REAL mc_cohesion    = 50;//kpa
    REAL mc_phi         = 0.001;
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

    REAL factor;
    TPZManVector<REAL, 3> bodyforce ( 3,0. );
    factor=1.;
    bodyforce[1]=-20.;

    material->SetPlasticity ( plasticstep );

    material->SetId ( 1 );

    material->SetWhichLoadVector ( 0 ); //option to compute the total internal force vecor fi=(Bt sigma+ N (b+gradu))

    material->SetLoadFactor ( factor );

    material->SetBodyForce ( bodyforce );

    cmesh->InsertMaterialObject(material);

    TPZFMatrix<REAL> val1(2,2,0.),val2(2,1,0.);

    val2(0,0)=1.;
    auto * bcleft = material->CreateBC(material,-1,3,val1,val2);//left

    val2(0,0)=1.;
    val2(1,0)=1.;
    auto * bcbottom = material->CreateBC(material,-2,3,val1,val2);//bottom

    val2(0,0)=1.;
    val2(1,0)=0;
    auto *  bcrigth = material->CreateBC(material,-3,3,val1,val2);//rigth


	cmesh->InsertMaterialObject(bcleft);

    cmesh->InsertMaterialObject(bcbottom);

    cmesh->InsertMaterialObject(bcrigth);

	cmesh->SetAllCreateFunctionsContinuousWithMem();

    cmesh->AutoBuild();

    std::ofstream print("cmesh.txt");
	cmesh->Print(print);

    return cmesh;
}

void PostProcessVariables ( TPZStack<std::string> &scalNames, TPZStack<std::string> &vecNames )
{




    scalNames.Push ( "PlStrainSqJ2" );
    scalNames.Push ( "PlStrainSqJ2El" );
    scalNames.Push ( "PlAlpha" );

    scalNames.Push ( "DisplacementX" );
    scalNames.Push ( "DisplacementY" );
    scalNames.Push ( "DisplacementZ" );
    vecNames.Push ( "DisplacementTotal" );
    vecNames.Push ( "PrincipalStress" );


}

void PostPlasticity(string vtkd,TPZCompMesh * cmesh)
{
    TPZPostProcAnalysis * postprocdeter = new TPZPostProcAnalysis();
    CreatePostProcessingMesh ( postprocdeter ,cmesh);

    TPZVec<int> PostProcMatIds ( 1,1 );

    TPZStack<std::string> PostProcVars, scalNames, vecNames;

    PostProcessVariables ( scalNames, vecNames );

    //string vtkd = "postprocessdeter.vtk";
    postprocdeter->DefineGraphMesh ( 2,scalNames,vecNames,vtkd );

    postprocdeter->PostProcess ( 0 );

    delete postprocdeter;
}
void  CreatePostProcessingMesh (TPZPostProcAnalysis * PostProcess ,TPZCompMesh * cmesh)
{
    if ( PostProcess->ReferenceCompMesh() != cmesh )
    {

        PostProcess->SetCompMesh ( cmesh );

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

void findid(TPZCompMesh *cmesh,TPZVec<REAL> coord) {

    std::ofstream print("loadvsdisplacement.dat");
    TPZGeoMesh * gmesh =  cmesh->Reference();

    int dim =2;

    int id;
    int nels = gmesh->NElements();
    for(int iel=0; iel<nels; iel++) {


        TPZVec<REAL> qsi(dim,0.);
        TPZGeoEl * gel = gmesh->ElementVec()[iel];

        TPZCompEl * cel = gel->Reference();



        if(gel->MaterialId()<0)continue;
        bool check = gel->ComputeXInverse(coord,qsi,1.e-5);

        if(check==true)
        {

            int nnodes=gel->NNodes();
            for(int inode=0; inode<nnodes; inode++)
            {
                TPZVec<REAL> co(dim);
                TPZGeoNode* node = gel->NodePtr(inode);
                node->GetCoordinates(co);
                id = node->Id();

                long elementid1;
                if(fabs(co[0]-coord[0])<1.e-3 && fabs(co[1]-coord[1])<1.e-3)
                {
//                     TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *> ( cel );
//
//                     TPZMaterialData data;
//
//                     data.fNeedsSol = true;
//
//                     //intel->InitMaterialData ( data );
//
//                     gmesh->FindElement ( coord, qsi, elementid1,2 );
//
//                     //intel->ComputeRequiredData ( data, qsi );
//
//                     cout << " elementid1  = "<<elementid1     << endl;
//  					cout << " data.sol  = "<<data.sol     << endl;
                    //return data.sol[0];
                }
            }


//                         data1.intLocPtIndex = ip;
//             data2.intLocPtIndex = ip;
//             data1.fNeedsSol = true;
//             intel1->InitMaterialData ( data1 );
//             data2.fNeedsSol = true;
//             intel1->InitMaterialData ( data2 );
//
//             intel1->ComputeRequiredData ( data1, qsi );
//             TPZGeoEl *gel1 = gmesh->FindElement ( data.x, qsi, elementid1,2 );
//
//             std::stringstream sout;
//             //TPZGeoEl *gel = cel->Reference();
//             gel->Print ( sout );
//             REAL co[8][2] = {{-1,-1},{1,-1},{1,1},{-1,1},{0,-1},{1,0},{0,1},{-1,0}};
//             for ( int p = 0; p<8; p++ )
//             {
//                 TPZManVector<REAL,3> par ( 2,0. ),x ( 3,0. );
//                 par[0] = co[p][0];
//                 par[1] = co[p][1];
//                 gel->X ( par, x );
//                 cout << "point " << p << "co " << x << std::endl;
//             }


        }
    }

}
