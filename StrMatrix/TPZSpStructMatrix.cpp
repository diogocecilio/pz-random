/**
 * @file
 * @brief Contains the implementation of the TPZSpStructMatrix methods. 
 */

#include "TPZSpStructMatrix.h"
#include "pzstrmatrix.h"

#include "pzgeoelbc.h"
#include "pzgmesh.h"
#include "pzcmesh.h"

#include "pzanalysis.h"
#include "pzsolve.h"
#include "pzstepsolver.h"

#include "pzdxmesh.h"
#include <fstream>

#include "pzelmat.h"

#include "pzysmp.h"
#include "pzmetis.h"
#include "pzbndcond.h"
#include "TPZTimer.h"

#include "pzlog.h"
#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.StrMatrix"));
#endif

using namespace std;

TPZStructMatrix * TPZSpStructMatrix::Clone(){
    return new TPZSpStructMatrix(*this);
}
TPZMatrix<STATE> * TPZSpStructMatrix::CreateAssemble(TPZFMatrix<STATE> &rhs,
                                              TPZAutoPointer<TPZGuiInterface> guiInterface){
	
    LOGPZ_DEBUG(logger,"TPZSpStructMatrix::CreateAssemble starting");
	
    long neq = fMesh->NEquations();
    if(fMesh->FatherMesh()) {
		cout << "TPZSpStructMatrix should not be called with CreateAssemble for a substructure mesh\n";
		return new TPZFYsmpMatrix<STATE>(0,0);
    }
    TPZMatrix<STATE> *stiff = Create();//new TPZFYsmpMatrix(neq,neq);
	//    TPZFYsmpMatrix *mat = dynamic_cast<TPZFYsmpMatrix *> (stiff);
    rhs.Redim(neq,1);
    //stiff->Print("Stiffness TPZFYsmpMatrix :: CreateAssemble()");
    TPZTimer before("Assembly of a sparse matrix");
    before.start();
    LOGPZ_DEBUG(logger,"TPZSpStructMatrix::CreateAssemble calling Assemble()");
	Assemble(*stiff,rhs,guiInterface);
    before.stop();
    std::cout << __PRETTY_FUNCTION__ << " " << before << std::endl;
    //    mat->ComputeDiagonal();
    //stiff->Print("Stiffness TPZFYsmpMatrix :: CreateAssemble()");
    LOGPZ_DEBUG(logger,"TPZSpStructMatrix::CreateAssemble exiting");
    return stiff;
}
TPZMatrix<STATE> * TPZSpStructMatrix::Create(){
    long neq = fEquationFilter.NActiveEquations();
	/*    if(fMesh->FatherMesh()) {
	 TPZSubCompMesh *smesh = (TPZSubCompMesh *) fMesh;
	 neq = smesh->NumInternalEquations();
	 }*/
    TPZFYsmpMatrix<STATE> * mat = new TPZFYsmpMatrix<STATE>(neq,neq);
	
    /**
     *Longhin implementation
	 */
    TPZStack<long> elgraph;
    TPZVec<long> elgraphindex;
    //    int nnodes = 0;
    fMesh->ComputeElGraph(elgraph,elgraphindex);
    /**Creates a element graph*/
    TPZMetis metis(elgraphindex.NElements() -1 ,fMesh->NIndependentConnects());
    metis.SetElementGraph(elgraph,elgraphindex);
	
    TPZVec<long> nodegraph;
    TPZVec<long> nodegraphindex;
    /**
     *converts an element graph structure into a node graph structure
     *those vectors have size ZERO !!!
     */
    metis.ConvertGraph(elgraph,elgraphindex,nodegraph,nodegraphindex);
    /**vector sizes*/
    long i;
    long nblock = nodegraphindex.NElements()-1;
    long totalvar = 0;
    long totaleq = 0;
    for(i=0;i<nblock;i++){
		long iblsize = fMesh->Block().Size(i);
		long iblpos = fMesh->Block().Position(i);
        long numactive = fEquationFilter.NumActive(iblpos, iblpos+iblsize);
        if (!numactive) {
            continue;
        }
		totaleq += iblsize;
		long icfirst = nodegraphindex[i];
		long iclast = nodegraphindex[i+1];
		long j;
		//longhin
		totalvar+=iblsize*iblsize;
		for(j=icfirst;j<iclast;j++) {
			long col = nodegraph[j];
			long colsize = fMesh->Block().Size(col);
			long colpos = fMesh->Block().Position(col);
            long numactive = fEquationFilter.NumActive(colpos, colpos+colsize);
            if (!numactive) {
                continue;
            }
			totalvar += iblsize*colsize;
		}
    }
	
    long ieq = 0;
    long pos = 0;
	
    nblock=fMesh->NIndependentConnects();
	
    long * Eq = new long[totaleq+1];
    long * EqCol = new long[totalvar];
    STATE * EqValue = new STATE [totalvar];
    for(i=0;i<nblock;i++){
		long iblsize = fMesh->Block().Size(i);
		long iblpos = fMesh->Block().Position(i);
        long numactive = fEquationFilter.NumActive(iblpos, iblpos+iblsize);
        if (!numactive) {
            continue;
        }
		long ibleq;
		for(ibleq=0; ibleq<iblsize; ibleq++) {
			Eq[ieq] = pos;
			long colsize,colpos,jbleq;
			long diagonalinsert = 0;
			long icfirst = nodegraphindex[i];
			long iclast = nodegraphindex[i+1];
			long j;
			for(j=icfirst;j<iclast;j++) {
				long col = nodegraph[j];
				if(!diagonalinsert && col > i)
				{
					diagonalinsert = 1;
					long colsize = fMesh->Block().Size(i);
					long colpos = fMesh->Block().Position(i);
                    TPZManVector<long> destindices(colsize);
                    for (long i=0; i<colsize; i++) {
                        destindices[i] = colpos+i;
                    }
                    fEquationFilter.Filter(destindices);
					long jbleq;
					for(jbleq=0; jbleq<destindices.size(); jbleq++) {
						//             if(colpos+jbleq == ieq) continue;
						EqCol[pos] = destindices[jbleq];
						EqValue[pos] = 0.;
						//            colpos++;
						pos++;
					}
				}
				colsize = fMesh->Block().Size(col);
				colpos = fMesh->Block().Position(col);
                if (fEquationFilter.NumActive(colpos, colpos+colsize) == 0) {
                    continue;
                }
                TPZManVector<long> destindices(colsize);
                for (long i=0; i<colsize; i++) {
                    destindices[i] = colpos+i;
                }
                fEquationFilter.Filter(destindices);
				for(jbleq=0; jbleq<destindices.size(); jbleq++) {
					EqCol[pos] = destindices[jbleq];
					EqValue[pos] = 0.;
					colpos++;
					pos++;
				}
			}
			if(!diagonalinsert)
			{
				diagonalinsert = 1;
				long colsize = fMesh->Block().Size(i);
				long colpos = fMesh->Block().Position(i);
                TPZManVector<long> destindices(colsize);
                for (long i=0; i<colsize; i++) {
                    destindices[i] = colpos+i;
                }
                fEquationFilter.Filter(destindices);
				long jbleq;
				for(jbleq=0; jbleq<destindices.size(); jbleq++) {
					//             if(colpos+jbleq == ieq) continue;
					EqCol[pos] = destindices[jbleq];
					EqValue[pos] = 0.;
					//            colpos++;
					pos++;
				}
			}
			ieq++;
		}
    }
    Eq[ieq] = pos;
    mat->SetData(Eq,EqCol,EqValue);
    return mat;
}

TPZSpStructMatrix::TPZSpStructMatrix(TPZCompMesh *mesh) : TPZStructMatrix(mesh)
{}

#ifndef STATE_COMPLEX
#include "pzmat2dlin.h"

int TPZSpStructMatrix::main() {
	
	int refine=5;
	int order=5;
	
	TPZGeoMesh gmesh;
	TPZCompMesh cmesh(&gmesh);
	double coordstore[4][3] = {{0.,0.,0.},{1.,0.,0.},{1.,1.,0.},
		{0.,1.,0.}};
	
	int i,j;
	TPZVec<REAL> coord(3,0.);
	for(i=0; i<4; i++) {
		// initializar as coordenadas do no em um vetor
		for (j=0; j<3; j++) coord[j] = coordstore[i][j];
		
		// identificar um espa� no vetor onde podemos armazenar
		// este vetor
		gmesh.NodeVec ().AllocateNewElement ();
		
		// initializar os dados do n�       gmesh.NodeVec ()[nodeindex].Initialize (i,coord,gmesh);
	}
	int el;
	TPZGeoEl *gel;
	for(el=0; el<1; el++) {
		
		// initializar os indices dos nos
		TPZVec<long> indices(4);
		for(i=0; i<4; i++) indices[i] = i;
		// O proprio construtor vai inserir o elemento na malha
		//       gel = new TPZGeoElQ2d(el,indices,1,gmesh);
		long index;
		gel = gmesh.CreateGeoElement(EQuadrilateral,indices,1,index);
	}
	gmesh.BuildConnectivity ();
	
	TPZVec<TPZGeoEl *> subel;
	
	cout << "Refinement ";
	cin >> refine;
	cout << endl;
	DebugStop();
//	UniformRefine(refine,gmesh);
	
	TPZGeoElBC gelbc(gel,4,-4);
	TPZMat2dLin *meumat = new TPZMat2dLin(1);
	TPZFMatrix<STATE> xk(1,1,1.),xc(1,2,0.),xf(1,1,1.);
	meumat->SetMaterial (xk,xc,xf);
	TPZMaterial * meumatptr(meumat);
	cmesh.InsertMaterialObject(meumatptr);
	
	TPZFMatrix<STATE> val1(1,1,0.),val2(1,1,0.);
	TPZMaterial * bnd = meumat->CreateBC (meumatptr,-4,0,val1,val2);
	cmesh.InsertMaterialObject(bnd);
	
	cout << "Interpolation order ";
	cin >> order;
	cout << endl;
	
	//TPZCompEl::gOrder = order;
	cmesh.SetDefaultOrder(order);
	
	cmesh.AutoBuild();
	//	cmesh.AdjustBoundaryElements();
	cmesh.InitializeBlock();
	
	ofstream output("outputPar.dat");
	TPZAnalysis an(&cmesh);
	
	TPZVec<int> numelconnected(cmesh.NEquations(),0);
	TPZSpStructMatrix mat(&cmesh);
	
	an.SetStructuralMatrix(mat);
	
	TPZStepSolver<STATE> sol;
	sol.SetJacobi(100,1.e-5,0);
	
	
	an.SetSolver(sol);
	an.Run(output);
	output.flush();
	cout.flush();
	return 0;
	
}
#endif
