#include "pzvec.h"
#include "pzadmchunk.h"
#include "pzcmesh.h"
#include "pzvec_extras.h"
#include "pzdebug.h"
#include "pzcheckgeom.h"

#include "pzgeoel.h"
#include "pzgnode.h"
#include "pzgeoelside.h"
#include "pzgeoelbc.h"

#include "pzintel.h"
#include "pzcompel.h"

#include "pzmatrix.h"

#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"

#include "pzmaterial.h"
#include "pzbndcond.h"
#include "pzelasmat.h"
#include "pzplaca.h"
#include "pzmat2dlin.h"
#include "pzmathyperelastic.h"
#include "pzmattest3d.h"
#include "pzmatplaca2.h"

#include "pzfunction.h"
#include "TPZVTKGeoMesh.h"

#include "pzlog.h"

#include <pzgeoel.h>
#include "pzgeoelbc.h"
#include "pzfmatrix.h"
#include "pzbstrmatrix.h"
#include <TPZGeoElement.h>
#include "TPZVTKGeoMesh.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZInterfaceEl.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "pzstepsolver.h"
#include "TPZGeoLinear.h"
#include "pzgeopoint.h"
#include "tpzgeoelrefpattern.h"
#include "TPZParFrontStructMatrix.h"
#include <math.h>
#include <cmath>
#include <iostream>
#include <chrono>
#include <cmath>
#include <set>

#include <iostream>
#include <fstream>
#include <string>
#include "pzgmesh.h"
#include "pzstack.h"
#include "TPZVTKGeoMesh.h"
#include "pzanalysis.h"
#include "pzbndcond.h"
#include "pzfunction.h"
#include "TPZDarcyFlow.h"


#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.adaptivity"));
static LoggerPtr loggerconv(Logger::getLogger("pz.adaptivity.conv"));
static LoggerPtr loggerpoint(Logger::getLogger("pz.adaptivity.points"));
#endif

#include <time.h>
#include <stdio.h>
#include <fstream>

TPZGeoMesh *CreateGeoMesh();
TPZCompMesh *CreateMesh(TPZGeoMesh *gmesh);
void SolveSelfWeigthBar(string filename);
void SolveBendingClampedBeam(string filename);

TPZGeoMesh *CreateGeoMeshBending();
TPZCompMesh *CreateMeshBending(TPZGeoMesh *gmesh);
TPZVec<REAL> findnodalsol(TPZCompMesh *cmesh,TPZVec<REAL> coord) ;

// bi-dimensional problem for elasticity
int main() {
    string filename =  "selfweigth.vtk";
    SolveSelfWeigthBar(filename);

    filename =  "clamped.vtk";
    SolveBendingClampedBeam(filename);

}


void SolveSelfWeigthBar(string filename)
{
	// Creating geometric mesh
	TPZGeoMesh *gmesh = CreateGeoMesh();

	// Creating computational mesh (approximation space and materials)
	int p = 3;
    TPZCompEl::SetgOrder(p);
    TPZCompMesh *cmesh = CreateMesh(gmesh);
	// Solving linear equations
	// Initial steps
	TPZAnalysis an (cmesh);
	TPZSkylineStructMatrix strskyl(cmesh);
	an.SetStructuralMatrix(strskyl);
	// Solver (is your choose)
	TPZStepSolver<REAL> *direct = new TPZStepSolver<REAL>;
	direct->SetDirect(ECholesky);
	an.SetSolver(*direct);
	delete direct;
	direct = 0;

	an.Run();

	// Post processing
	TPZManVector<std::string> scalarnames(3), vecnames(1);
	scalarnames[0] = "SigmaX";
	scalarnames[1] = "SigmaY";
	scalarnames[2] = "TauXY";
	vecnames[0] = "displacement";
	//vecnames[1] = "";
	an.DefineGraphMesh(2,scalarnames,vecnames,filename);

	an.PostProcess(0);
}
void SolveBendingClampedBeam(string filename)
{
	// Creating geometric mesh
	TPZGeoMesh *gmesh = CreateGeoMeshBending();

	// Creating computational mesh (approximation space and materials)
	int p = 2;
    TPZCompEl::SetgOrder(p);
    TPZCompMesh *cmesh = CreateMeshBending(gmesh);
	// Solving linear equations
	// Initial steps
	TPZAnalysis an (cmesh);
	TPZSkylineStructMatrix strskyl(cmesh);
	an.SetStructuralMatrix(strskyl);
	// Solver (is your choose)
	TPZStepSolver<REAL> *direct = new TPZStepSolver<REAL>;
	direct->SetDirect(ECholesky);
	an.SetSolver(*direct);
	delete direct;
	direct = 0;

	an.Run();

    //cout << "\n SOLUTION "<< endl;
    TPZFMatrix<REAL> sol = an.Solution();
    TPZVec<REAL> coord(2);
    //coord[0]=0.;
    //coord[1]=0.;
    //TPZVec<REAL> sol = findnodalsol(cmesh,coord);

    cout << "Displacement solution in coord x = "<< coord[0] << " y ="  << coord[1] <<endl;
    cout << "ux = "<< sol(0,0)<< " uy ="  << sol(1,0) <<endl;

	// Post processing
	TPZManVector<std::string> scalarnames(3), vecnames(1);
	scalarnames[0] = "SigmaX";
	scalarnames[1] = "SigmaY";
	scalarnames[2] = "TauXY";
	vecnames[0] = "displacement";
	//vecnames[1] = "";
	an.DefineGraphMesh(2,scalarnames,vecnames,filename);

	an.PostProcess(0);
}
TPZGeoMesh *CreateGeoMeshBending()
{
    REAL co[4][2] = {{0.,0.},{1000,0},{1000,100},{0,100}};
    long indices[1][4] = {{0,1,2,3}};
    TPZGeoEl *elvec[1];
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    gmesh->SetDimension(2);
    long nnode = 4;
    long nod;
    for(nod=0; nod<nnode; nod++) {
        long nodind = gmesh->NodeVec().AllocateNewElement();
        TPZVec<REAL> coord(2);
        coord[0] = co[nod][0];
        coord[1] = co[nod][1];
        gmesh->NodeVec()[nodind] = TPZGeoNode(nod,coord,*gmesh);
    }

    long el;
    long nelem = 1;
    for(el=0; el<nelem; el++) {
        TPZVec<long> nodind(4);
        for(nod=0; nod<4; nod++) nodind[nod]=indices[el][nod];
        //    elvec[el] = new TPZGeoElQ2d(el,nodind,1);
        long index;
        elvec[el] = gmesh->CreateGeoElement(EQuadrilateral,nodind,1,index);
    }

    TPZVec <long> TopoLine ( 2 );

    long index;
    TopoLine[0] = 1;
    TopoLine[1] = 2;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( index, TopoLine, - 1, *gmesh );//clamped in right side

   TPZVec <long> node(1);
   node[0]=0;
   new TPZGeoElRefPattern< pzgeom::TPZGeoPoint> ( index, node, - 2, *gmesh );//load node

    node[0]=1;
   new TPZGeoElRefPattern< pzgeom::TPZGeoPoint> ( index, node, - 3, *gmesh );//bottomrigth node

    gmesh->BuildConnectivity();


    for(int d=0;d<2;d++) {
    int nel = gmesh->NElements();
    for (int iel=0; iel<nel; iel++) {
        TPZManVector<TPZGeoEl *> subels;
        TPZGeoEl *gel = gmesh->ElementVec()[iel];
        gel->Divide(subels);
    }
	}

	return gmesh;
}

TPZCompMesh *CreateMeshBending(TPZGeoMesh *gmesh)
{
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(TPZCompEl::GetgOrder());

    //TPZElasticityMaterial(int id, REAL E, REAL nu, REAL fx, REAL fy, int planestress = 1);
    TPZMaterial * mat = new TPZElasticityMaterial(1,10000.,0.2,0.,0.);//selfweigth

    cmesh->SetDimModel(2);

    TPZFMatrix<REAL> val1(2,2,0.),val2(2,1,0.);
	TPZMaterial *bcload,*bcclamp,*bcnode;

    val2(0,0)=1.;
    val2(1,0)=1.;
    bcclamp = mat->CreateBC(mat,-1,3,val1,val2);//clamped line restrictions

    val2(0,0)=0.;
    val2(1,0)=1.;
    bcnode = mat->CreateBC(mat,-3,3,val1,val2);//bottomrigth node restrictions

    val2(0,0)=0.;
    val2(1,0)=-100.;
    bcload = mat->CreateBC(mat,-2,1,val1,val2);//-100 N in y direction node 4


    cmesh->InsertMaterialObject(mat);
	cmesh->InsertMaterialObject(bcclamp);
    cmesh->InsertMaterialObject(bcnode);
    cmesh->InsertMaterialObject(bcload);


	cmesh->SetAllCreateFunctionsContinuous();

    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();

    return cmesh;
}

TPZGeoMesh *CreateGeoMesh() {

    REAL co[4][2] = {{0.,0.},{1,0},{1,10},{0,10}};
    long indices[1][4] = {{0,1,2,3}};
    TPZGeoEl *elvec[1];
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    long nnode = 4;
    long nod;
    for(nod=0; nod<nnode; nod++) {
        long nodind = gmesh->NodeVec().AllocateNewElement();
        TPZVec<REAL> coord(2);
        coord[0] = co[nod][0];
        coord[1] = co[nod][1];
        gmesh->NodeVec()[nodind] = TPZGeoNode(nod,coord,*gmesh);
    }

    long el;
    long nelem = 1;
    for(el=0; el<nelem; el++) {
        TPZVec<long> nodind(4);
        for(nod=0; nod<4; nod++) nodind[nod]=indices[el][nod];
        //    elvec[el] = new TPZGeoElQ2d(el,nodind,1);
        long index;
        elvec[el] = gmesh->CreateGeoElement(EQuadrilateral,nodind,1,index);
    }

    TPZVec <long> TopoLine ( 2 );

    TopoLine[0] = 0;
    TopoLine[1] = 1;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( 7, TopoLine, - 1, *gmesh );//bottom


    TopoLine[0] = 2;
    TopoLine[1] = 3;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( 8, TopoLine, - 2, *gmesh );//top

   TPZVec <long> node(1);
   node[0]=3;
   new TPZGeoElRefPattern< pzgeom::TPZGeoPoint> ( 9, node, - 3, *gmesh );//top
//     long index;
//     gmesh->CreateGeoElement(EPoint,3,-3,index);
//
//     TPZGeoElBC gbc3(elvec[0],6,-3);

    gmesh->BuildConnectivity();


    for(int d=0;d<2;d++) {
    int nel = gmesh->NElements();
    for (int iel=0; iel<nel; iel++) {
        TPZManVector<TPZGeoEl *> subels;
        TPZGeoEl *gel = gmesh->ElementVec()[iel];
        gel->Divide(subels);
    }
	}

	return gmesh;
}


TPZCompMesh *CreateMesh(TPZGeoMesh *gmesh) {

    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(TPZCompEl::GetgOrder());
//TPZElasticityMaterial(int id, REAL E, REAL nu, REAL fx, REAL fy, int planestress = 1);
    // Creating elasticity material
    TPZMaterial * mat = new TPZElasticityMaterial(1,10000.,0.2,0.,-1.);//selfweigth
   // TPZMaterial * mat = new TPZElasticityMaterial(1,10000.,0.2,0.,0.);
	// Creating four boundary condition
    TPZFMatrix<REAL> val1(2,2,0.),val2(2,1,0.);
	TPZMaterial *bctop,*bcload,*bcpoint;

    val2(1,0)=1.;
    bctop = mat->CreateBC(mat,-2,3,val1,val2);

    val2(0,0)=1.;
    bcpoint = mat->CreateBC(mat,-3,3,val1,val2);

    //val2(1,0)=-1.;
    //bcload = mat->CreateBC(mat,-1,1,val1,val2);
    cmesh->InsertMaterialObject(mat);
	// Inserting boundary conditions into computational mesh
	cmesh->InsertMaterialObject(bctop);

    cmesh->InsertMaterialObject(bcpoint);
    //cmesh->InsertMaterialObject(bcload);

	cmesh->SetAllCreateFunctionsContinuous();

    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();

    return cmesh;
}
TPZVec<REAL> findnodalsol(TPZCompMesh *cmesh,TPZVec<REAL> coord) {

    TPZGeoMesh * gmesh =  cmesh->Reference();
    int dim = gmesh->Dimension()+1;
    TPZVec<REAL> xd(dim,0.);
    TPZVec<REAL> mpt(dim,0.);
    xd[0] =coord[0];
    xd[1] = coord[1];
    int id;
    int nels = gmesh->NElements();
    for(int iel=0; iel<nels; iel++) {


        TPZVec<REAL> qsi(2,0.);
        TPZGeoEl * gel = gmesh->ElementVec()[iel];

        TPZCompEl * cel = gel->Reference();



        if(gel->MaterialId()<0)continue;
        bool check = gel->ComputeXInverse(xd,qsi,1.e-5);

        if(check==true)
        {

            int nnodes=gel->NNodes();
            for(int inode=0; inode<nnodes; inode++)
            {
                TPZVec<REAL> co(dim);
                TPZGeoNode* node = gel->NodePtr(inode);
                node->GetCoordinates(co);
                id = node->Id();

                if(fabs(co[0]-xd[0])<1.e-3 && fabs(co[1]-xd[1])<1.e-3)
                {
                    TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *> ( cel );

                    TPZMaterialData data;

                    data.fNeedsSol = true;

                    intel->InitMaterialData ( data );

                    intel->ComputeRequiredData ( data, qsi );

 					//cout << " data.sol  = "<<data.sol     << endl;
                    return data.sol[0];
                }
            }


        }
    }

}
