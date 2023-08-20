//$Id: pzelastoplasticanalysis.cpp,v 1.27 2010-11-23 18:58:05 diogo Exp $
#include "pzelastoplasticanalysis.h"
#include "pzcmesh.h"
#include "pzvec.h"
#include "pzmanvector.h"
#include "checkconv.h"
#include "pzstrmatrix.h"
#include "pzelastoplastic.h"
#include "tpzautopointer.h"
#include "pzcompelwithmem.h"
#include "pzelastoplasticmem.h"
#include "pzblockdiag.h"
#include "TPZSpStructMatrix.h"
#include "pzfstrmatrix.h"
#include "pzbdstrmatrix.h"
#include "pzstepsolver.h"
#include "pzmaterial.h"
#include "pzbndcond.h"
#include "pzelastoplastic2D.h"

#include "pzbuildmultiphysicsmesh.h"
#include "TPZPlasticStepPV.h"
#include "TPZPlasticStepPV.h"
#include "TPZElasticResponse.h"
#include "TPZYCMohrCoulombPV.h"
#include "pzelastoplastic2D.h"
#include "pzelastoplastic.h"
#include "TPZMohrCoulombVoigt.h"
#include "TPZPlasticStepVoigt.h"
#include <map>
#include <set>
#include <stdio.h>
#include <fstream>

#include "pzsolve.h"

#include "pzlog.h"

//typedef TPZPlasticStepVoigt<TPZMohrCoulombVoigt,TPZElasticResponse> LEMC;

typedef TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> LEMC;
typedef   TPZMatElastoPlastic2D <LEMC, TPZElastoPlasticMem > plasticmat;

#ifdef LOG4CXX
static LoggerPtr EPAnalysisLogger ( Logger::getLogger ( "pz.analysis.elastoplastic" ) );
static LoggerPtr loggertest ( Logger::getLogger ( "testing" ) );
#endif

using namespace std;


TPZElastoPlasticAnalysis::TPZElastoPlasticAnalysis() : TPZNonLinearAnalysis(), fPrecond ( NULL )
{
    //Mesh()->Solution().Zero(); already performed in the nonlinearanalysis base class
    //fSolution.Zero();
}

TPZElastoPlasticAnalysis::TPZElastoPlasticAnalysis ( TPZCompMesh *mesh,std::ostream &out ) : TPZNonLinearAnalysis ( mesh,out ), fPrecond ( NULL )
{

    int numeq = fCompMesh->NEquations();
    fCumSol.Redim ( numeq,1 );
    fCumSol.Zero();
    fSolution.Redim ( numeq,1 );
    fSolution.Zero();

    LoadSolution();
}

TPZElastoPlasticAnalysis::~TPZElastoPlasticAnalysis()
{
    if ( fPrecond ) delete fPrecond;

#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout << "<<< TPZElastoPlasticAnalysis::~TPZElastoPlasticAnalysis() *** Killing Object\n";
        LOGPZ_INFO ( EPAnalysisLogger,sout.str().c_str() );
    }
#endif
}

REAL TPZElastoPlasticAnalysis::LineSearch ( const TPZFMatrix<REAL> &Wn, const TPZFMatrix<REAL> &DeltaW, TPZFMatrix<REAL> &NextW, REAL RhsNormPrev, REAL &RhsNormResult, int niter, bool & converging )
{

    TPZFMatrix<REAL> Interval = DeltaW;

#ifdef DEBUG
    {
        TPZNonLinearAnalysis::LoadSolution ( Wn );
        AssembleResidual();
        STATE normprev = Norm ( fRhs );
        if ( fabs ( normprev - RhsNormPrev ) > 1.e-6 )
        {
            std::stringstream sout;
            sout << "Norm of Wn " << Norm ( Wn ) << std::endl;
            sout << "Input previous norm " << RhsNormPrev << " Computed Norm " << normprev;
            LOGPZ_ERROR ( EPAnalysisLogger, sout.str() )
        }
    }
#endif
    REAL scalefactor = 1.;
    int iter = 0;
    do
    {
        Interval *= scalefactor;
        NextW = Wn;
        NextW += Interval;
        TPZNonLinearAnalysis::LoadSolution ( NextW );
        AssembleResidual();
        {
            std::cout << "Vertical strain change" << fSolution ( fSolution.Rows()-1,0 ) << std::endl;
            static int count = 0;
            {
                std::stringstream filename,varname;
                filename << "Sol." << count << ".txt";
                varname << "DelSol" << count << " = ";
                ofstream out ( filename.str().c_str() );
                Interval.Print ( varname.str().c_str(),out,EMathematicaInput );
            }
            std::stringstream filename,varname;
            filename << "Rhs." << count << ".txt";
            varname << "Rhs" << count++ << " = ";
            ofstream out ( filename.str().c_str() );
            fRhs.Print ( varname.str().c_str(),out,EMathematicaInput );
        }
        RhsNormResult = Norm ( fRhs );
#ifndef PLASTICITY_CLEAN_OUT
        std::cout << "Scale factor " << scalefactor << " resnorm " << RhsNormResult << std::endl;
#endif
        scalefactor *= 0.5;
        iter++;
    }
    while ( RhsNormResult > RhsNormPrev && iter < niter );
    if ( fabs ( RhsNormResult - RhsNormPrev ) <1.e-6 )
    {
        converging=false;
    }
    else
    {
        converging=true;
    }
    scalefactor *= 2.;
    return scalefactor;

}//void





void TPZElastoPlasticAnalysis::SetUpdateMem ( int update )
{
    if ( !fCompMesh ) return;

    std::map<int, TPZMaterial *> & refMatVec = fCompMesh->MaterialVec();

    std::map<int, TPZMaterial * >::iterator mit;

    TPZMatWithMem<TPZElastoPlasticMem> * pMatWithMem; // defined in file pzelastoplastic.h
    TPZMatWithMem<TPZPoroElastoPlasticMem> * pMatWithMem2; // define in file pzporous.h

    for ( mit=refMatVec.begin(); mit!= refMatVec.end(); mit++ )
    {
        pMatWithMem = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> ( mit->second );
        if ( pMatWithMem != NULL )
        {
            pMatWithMem->SetUpdateMem ( update );
        }
        pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZPoroElastoPlasticMem> *> ( mit->second );
        if ( pMatWithMem2 != NULL )
        {
            pMatWithMem2->SetUpdateMem ( update );
        }
    }

}

#include "pzelasmat.h"

REAL TPZElastoPlasticAnalysis::AcceptSolution ( const int ResetOutputDisplacements )
{

    TPZMaterial *mat = fCompMesh->FindMaterial ( 1 );
    if ( !mat )
    {
        DebugStop();
    }
    TPZElasticityMaterial *elasmat = dynamic_cast<TPZElasticityMaterial *> ( mat );
    if ( elasmat )
    {
        // the material is linear
        return 0.;
    }


    if ( ResetOutputDisplacements )
    {
        fCumSol.Zero();
    }
    else
    {
        fCumSol += fSolution;
    }

#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout << ">>> TTPZElastoPlasticAnalysis::AcceptSolution *** "
             << " with Norm(fCumSol) = " << Norm ( fCumSol );
        LOGPZ_INFO ( EPAnalysisLogger,sout.str().c_str() );
    }
#endif

    this->SetUpdateMem ( true );

    fRhs.Zero();

    AssembleResidual();
    REAL norm = Norm ( fRhs );

    this->SetUpdateMem ( false );

    fSolution.Zero();

    LoadSolution();


    return norm;
}

/** @brief Load the solution into the computable grid, transfering it to the multi physics meshes */
void TPZElastoPlasticAnalysis::LoadSolution()
{
    TPZNonLinearAnalysis::LoadSolution();
    if ( this->IsMultiPhysicsConfiguration() )
    {
        TPZBuildMultiphysicsMesh::TransferFromMultiPhysics ( fMeshVec, fMultiPhysics );
    }

}



void TPZElastoPlasticAnalysis::CheckConv ( std::ostream &out, REAL range )
{

#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout << ">>> TPZElastoPlasticAnalysis::CheckConv() ***"
             << "\nEntering method with parameters:"
             << "\n range = " << range;
        LOGPZ_INFO ( EPAnalysisLogger,sout.str().c_str() );
    }
#endif

    int numeq = fCompMesh->NEquations();

    TPZFMatrix<REAL> rangeMatrix ( numeq, 1, range );

    TPZVec<REAL> coefs ( 1,1. );

    CheckConvergence ( *this,fSolution,rangeMatrix,coefs );

}

void TPZElastoPlasticAnalysis::ComputeTangent ( TPZFMatrix<REAL> &tangent, TPZVec<REAL> &coefs, int icase )
{

    int neq = fCompMesh->NEquations();
    tangent.Redim ( neq,neq );
    TPZFMatrix<REAL> rhs ( neq,1 );
    TPZFStructMatrix substitute ( Mesh() );
    TPZAutoPointer<TPZGuiInterface> guiInterface ( 0 );
    substitute.Assemble ( tangent,rhs,guiInterface );
//	TPZStructMatrix::Assemble(tangent, rhs, *Mesh());
}

int TPZElastoPlasticAnalysis::NumCases()
{
    return 1;
}

void TPZElastoPlasticAnalysis::Residual ( TPZFMatrix<REAL> &residual, int icase )
{
    int neq = fCompMesh->NEquations();
//	TPZFMatrix<REAL> tangent(neq,neq);
    residual.Redim ( neq,1 );
    TPZFStructMatrix substitute ( Mesh() );
    TPZAutoPointer<TPZGuiInterface> guiInterface ( 0 );
    substitute.Assemble ( residual,guiInterface );
//	TPZStructMatrix::Assemble(/*tangent,*/ residual, *Mesh());
    residual *= -1;
}

void TPZElastoPlasticAnalysis::SetPrecond ( TPZMatrixSolver<REAL> &precond )
{
    if ( fPrecond ) delete fPrecond;
    fPrecond = ( TPZMatrixSolver<REAL> * ) precond.Clone();
}

void TPZElastoPlasticAnalysis::UpdatePrecond()
{
    if ( fPrecond )
    {
        TPZMatrix<REAL> * pMatrix = TPZAnalysis::fSolver->Matrix().operator->();
        TPZMatrix<REAL> * pPrecondMat = fPrecond->Matrix().operator->();
        pPrecondMat->Zero();
        TPZBlockDiagonal<REAL> *pBlock = dynamic_cast<TPZBlockDiagonal<REAL> *> ( pPrecondMat );
        pBlock->BuildFromMatrix ( *pMatrix );
    }
}
#include <pzskylstrmatrix.h>
void TPZElastoPlasticAnalysis::SetBiCGStab ( int numiter, REAL tol )
{

    TPZSpStructMatrix StrMatrix(Mesh());
    //TPZSkylineStructMatrix StrMatrix ( Mesh() );
    StrMatrix.SetNumThreads ( 12 );
    this->SetStructuralMatrix ( StrMatrix );
    TPZMatrix<REAL> * mat = StrMatrix.Create();

    TPZBlockDiagonalStructMatrix strBlockDiag ( Mesh() );
    TPZStepSolver<REAL> Pre;
    TPZBlockDiagonal<REAL> * block = new TPZBlockDiagonal<REAL>();

    strBlockDiag.AssembleBlockDiagonal ( *block ); // just to initialize structure
    Pre.SetMatrix ( block );
    // Pre.SetDirect(ELU);
    Pre.SetDirect ( ELDLt );
    TPZStepSolver<REAL> Solver;
    Solver.SetBiCGStab ( numiter, Pre, tol, 0 );
    Solver.SetMatrix ( mat );
    this->SetSolver ( Solver );
    this->SetPrecond ( Pre );

}


void TPZElastoPlasticAnalysis::SetBiCGStab_Jacobi ( int numiter, REAL tol )
{
    TPZSpStructMatrix StrMatrix ( Mesh() );

    this->SetStructuralMatrix ( StrMatrix );
    TPZMatrix<REAL> * mat = StrMatrix.Create();

    TPZBlockDiagonalStructMatrix strBlockDiag ( Mesh() );
    TPZStepSolver<REAL> Pre;
    TPZBlockDiagonal<REAL> * block = new TPZBlockDiagonal<REAL>();

    strBlockDiag.AssembleBlockDiagonal ( *block ); // just to initialize structure
    Pre.SetMatrix ( block );
    //    Pre.SetDirect(ELU);
    Pre.SetDirect ( ELDLt );
    //Pre.SetJacobi(numiter, tol, 0);
    TPZStepSolver<REAL> Solver;
    Solver.SetBiCGStab ( numiter, Pre, tol, 0 );
    Solver.SetMatrix ( mat );
    this->SetSolver ( Solver );
    this->SetPrecond ( Pre );
}

void TPZElastoPlasticAnalysis::SetLU()
{
#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout << ">>> TPZElastoPlasticAnalysis::SetLU() ***\n";
        LOGPZ_INFO ( EPAnalysisLogger,sout.str().c_str() );
    }
#endif

    TPZFStructMatrix StrMatrix ( Mesh() );
    this->SetStructuralMatrix ( StrMatrix );

    TPZMatrix<REAL> * mat = StrMatrix.Create();

    TPZStepSolver<REAL> Solver;
    //Solver.SetDirect(ELU);// ECholesky -> simétrica e positiva definida
    Solver.SetDirect ( ELU );
    Solver.SetMatrix ( mat );

    this->SetSolver ( Solver );
}

void TPZElastoPlasticAnalysis::TransferSolution ( TPZPostProcAnalysis & ppanalysis )
{
    TPZFMatrix<REAL> bkpSolution = fSolution;


    fSolution = fCumSol;
//	 fSolution.Print();
    LoadSolution();//Carrega a solucao convergida no analysis
    //passa o cum sol para o post
    ppanalysis.TransferSolution();//Transfere solucao convergida para o pos processamento


    fSolution = bkpSolution;

    LoadSolution();
}
void TPZElastoPlasticAnalysis::LoadSolution ( TPZFMatrix<STATE> & loadsol )
{
    fSolution = loadsol;
    LoadSolution();
}
REAL TPZElastoPlasticAnalysis::LineSearch ( const TPZFMatrix<STATE> &Wn, TPZFMatrix<STATE> DeltaW, TPZFMatrix<STATE> &NextW, REAL tol, int niter )
{

    REAL error = 2.*tol+1.;
    REAL A = 0.1, B = 2., L = 0, M = 0.;
    TPZFMatrix<STATE> ak, bk, lambdak, muk, Interval;
    REAL NormResLambda = 0., NormResMu = 0.;
    //ak = Wn + 0.1 * DeltaW
    ak = DeltaW;
    ak *= A;
    ak += Wn;
    //bk = Wn + 2. DeltaW
    bk = DeltaW;
    bk *= B;
    bk += Wn;
    //Interval = (bk-ak)
    Interval = bk;
    Interval -= ak;
    int iter = 0;
    int KeptVal = -1; //0 means I have residual(labmda); 1 means I have residual(mu); -1 means I have nothing
    while ( error > tol && iter < niter )
    {
        iter++;
        //cout << "a  " << std::endl;
        if ( KeptVal != 0 )
        {
            L = 0.382* ( B-A )+A;
            //lambdak = ak + 0.382*(bk-ak)
            lambdak = Interval;
            lambdak *= 0.382;
            lambdak += ak;
            //computing residual
            LoadSolution ( lambdak );
#ifdef PZ_LOG
            LOGPZ_DEBUG ( logger,"After LoadSolution" )
#endif
            //		LogWellSolution(*this->Mesh(), 6);
            this->AssembleResidual();
#ifdef PZ_LOG
            LOGPZ_DEBUG ( logger,"After AssembleResidual" )
#endif
            //		LogWellSolution(*this->Mesh(), 6);
            NormResLambda = Norm ( fRhs );
        }

        if ( KeptVal != 1 )
        {
            //muk = ak + 0.618*(bk-ak)
            M = 0.618* ( B-A )+A;
            muk = Interval;
            muk *= 0.618;
            muk += ak;
            LoadSolution ( muk );
            //cout << "b  " << std::endl;
            //muk.Print(std::cout);
            this->AssembleResidual();
            //cout << "c " << std::endl;
            NormResMu = Norm ( fRhs );
        }

        if ( NormResLambda > NormResMu )
        {
            A = L;
            L = M;
            ak = lambdak;
            lambdak = muk;
            NormResLambda = NormResMu;
            KeptVal = 0;
        }
        else
        {
            B = M;
            M = L;
            bk = muk;
            muk = lambdak;
            NormResMu = NormResLambda;
            KeptVal = 1;
        }
        //error = Norm(bk-ak)
        Interval = bk;
        Interval -= ak;
        error = Norm ( Interval );

        //alpha shall be alpha <= 1
        if ( A > 1. && B > 1. ) break;

    }//while

    double ALPHA = 0.5* ( A + B );
    NextW = ak;
    NextW += bk;
    NextW *= 0.5;


#ifdef PZDEBUGLINESEARCH
    //debug: valor do alpha
    TPZFMatrix<REAL> alpha;
    alpha = NextW;
    alpha -= Wn;
    REAL sum = 0.;
    int ncontrib = 0;
    for ( int i = 0; i < alpha.Rows(); i++ )
    {
        if ( DeltaW ( i,0 ) )
        {
            alpha ( i,0 ) = alpha ( i,0 ) /DeltaW ( i,0 );
            sum += alpha ( i,0 );
            ncontrib++;
        }
    }
    //REAL MeanAlpha = sum/ncontrib;
    alphafile << /*MeanAlpha << "\t" <<*/ "ALPHA = " << ALPHA << "\n";
    alphafile.flush();
#endif

    if ( ALPHA > 1. ) //alpha shall be alpha <= 1
    {
        NextW = Wn;
        NextW += DeltaW;
#ifdef PZDEBUGLINESEARCH
        alphafile << "ALPHA LIMIT APPLIED. Alpha = 1.\n";
#endif
        return 1.;
    }

    return ALPHA;

}//void


void TPZElastoPlasticAnalysis::IterativeProcess ( std::ostream &out,REAL tol,int numiter )
{

    int iter = 0;
    REAL error = 1.e10;
    int numeq = fCompMesh->NEquations();
    TPZAutoPointer<TPZMatrix<REAL> > K;
    TPZFMatrix<REAL> rhs,du;
    TPZFMatrix<STATE> prevsol ( fSolution );
    while ( error > tol && iter < numiter )
    {
        cout << "assemble"<< endl;
        Assemble();
        TPZAutoPointer<TPZMatrix<REAL> > K = this->fSolver->Matrix();
        cout << "solve"<< endl;
        int type=0;//llt
        //SolveEigenSparse ( type, K, fRhs, fSolution );
        //rhs=fRhs;
        //SolveEigen ( K, rhs, du );
        //SolveEigenSparse ( type, K, rhs, du );
        Solve();
        //fSolution = du;
        if ( true )
        {
            TPZFMatrix<STATE> nextSol;
            REAL LineSearchTol = 1e-3 * Norm ( fSolution );
            const int niter = 10;
            this->LineSearch ( prevsol, fSolution, nextSol, LineSearchTol, niter );
            fSolution = nextSol;
        }
        //Solve();
        //fSolution = du;
        TPZFMatrix<STATE> sol = fSolution;
        sol += prevsol;


        prevsol -= fSolution;
        REAL normDeltaSol = Norm ( prevsol );
        prevsol = fSolution;
        this->LoadSolution ( fSolution );
        this->AssembleResidual();
        REAL norm = Norm ( fRhs );

        //cout << "Iteracao n : " << ( iter+1 ) << " : normas |Delta(Un)| e |Delta(rhs)| : " << normDeltaSol << " / " << norm << endl;

        error = norm;
        iter++;

    }

}

bool TPZElastoPlasticAnalysis::IterativeProcess ( std::ostream &out,REAL tol,int numiter, bool linesearch, bool checkconv,int &iters )
{

    int iter = 0;
    REAL error = 1.e10,error2=1.e10;
    int numeq = fCompMesh->NEquations();

 //   cout << "number of equations = " << numeq <<endl;

    TPZFMatrix<STATE> prevsol ( fSolution );
    if ( prevsol.Rows() != numeq ) prevsol.Redim ( numeq,1 );

    if ( checkconv )
    {
        TPZVec<REAL> coefs ( 1,1. );
        TPZFMatrix<STATE> range ( numeq,1,1. );
        CheckConvergence ( *this,fSolution,range,coefs );
    }
    bool a=true,b=true,c=true;
    // Assemble();
// 	 TPZAutoPointer<TPZMatrix<REAL> > K = this->fSolver->Matrix();
// 	TPZFMatrix<STATE> rhs =fRhs;
// 	TPZFMatrix<STATE> du;
// 	//chrono::steady_clock sc1;
// 	//auto start = sc1.now();
// 	SolveEigenSparse(0, K, rhs, du );
// 	//auto end = sc1.now();
// 	//auto time_span = static_cast<chrono::duration<double>> ( end - start );
// 	//cout << "| total time taken to solve eigen=  " << time_span.count()<< std::endl;
// 	 fSolution=du;
    AssembleResidual();
  //  std::cout << "asdasdsad " <<std::endl;
    REAL normrhs0 = Norm ( fRhs );
  //  cout << "normrhs0 = " << normrhs0 << endl;
    while ( a && b && c )
    {

//        fSolution.Redim(0,0);
        Assemble();

        chrono::steady_clock sc;
        auto start = sc.now();
        if ( false )
        {
            Eigen::initParallel();
            int n=10;
            //omp_set_num_threads(n);
            setNbThreads ( n );
            //ConjugateGradient = 2 //BiCGSTAB =5
            int type=0;
            type = 1;
            TPZAutoPointer<TPZMatrix<REAL> > K = this->fSolver->Matrix();
            TPZFMatrix<STATE> rhs =fRhs;
            TPZFMatrix<STATE> du;
            //SolveEigen ( K, rhs, du );
            //cout << "solving... "<< std::endl;
            SolveEigenSparse ( type, K, rhs, du );
            fSolution=du;
            auto end = sc.now();
            auto time_span = static_cast<chrono::duration<double>> ( end - start );
           // cout << "| total time taken to solve eigen=  " << time_span.count() << std::endl;
        }
        else
        {
            //cout <<  "sdadas" << endl;
            Solve();
            //cout <<  "aaaaaaaaaa" << endl;
            auto end = sc.now();
            auto time_span = static_cast<chrono::duration<double>> ( end - start );
         //   cout << "| total time taken to solve PZ=  " << time_span.count()<< std::endl;
        }


        //cout << "a  " << std::endl;
        if ( linesearch )
        {
            TPZFMatrix<STATE> nextSol;
            //REAL LineSearchTol = 1e-3 * Norm(fSolution);
            REAL LineSearchTol = 0.001 * Norm ( fSolution );
            const int niter =20;
            this->LineSearch ( prevsol, fSolution, nextSol, LineSearchTol, niter );
            fSolution = nextSol;
        }
        else
        {
            TPZFMatrix<STATE> sol = fSolution;
            sol += prevsol;
        }

        //cout << "b  " << std::endl;
        prevsol -= fSolution;
        //REAL normDeltaSol = Norm(prevsol)/unorm0;
        REAL normu = Norm ( prevsol );

        prevsol = fSolution;
        this->LoadSolution ( fSolution );
        this->AssembleResidual();
        //cout << "c  " << std::endl;
        REAL normf  = Norm ( fRhs );
     //   cout << "Iteracao n : " << ( iter ) << " : normas |Delta(Un)| e |Delta(rhs)/rhs0| : " << normu << " / " << normf/normrhs0 << " | tol = "<<tol << endl;
        a = iter < numiter ;
        b =error2 > tol *1.e-3;
        c= error > tol;

        //if( normDeltaSol>100 || iter >=numiter  || ((normDeltaSol - error2) > 1.e-9 && (NormResLambda - error) > 1.e-9) ) {
        // if((normu>100 || iter >=numiter  ||(normu - error2) > 1.e-3)&& iter>5) {
        // if((normu>1 && iter>5 && (normu - error2) > 1.e-3)|| iter >=numiter) {
//                 if(  ( iter >=numiter || ( iter>2 && normu >1 ) ) ||(normu - error2) > 1.e-3 || fabs((normf - error)) > 1.e-9 ) {
//             cout << "\nDivergent Method\n";
//             return false;
//         }
        if ( ( iter >=numiter || ( iter>2 && normu >1 ) ) )
        {
            //cout << "\nDivergent Method\n";
            return false;
        }
        error = normf;
        error2=normu;
        iter++;
        out.flush();

    }
    iters=iter;
   // cout << "Iteracao n : " << ( iter ) << "Norm ( prevsol ) = "<<Norm ( prevsol ) << "Norm ( fRhs ) = "<<Norm ( fRhs ) << endl;
    return true;
}

void TPZElastoPlasticAnalysis::IterativeProcess ( std::ostream &out,REAL tol,int numiter, bool linesearch, bool checkconv )
{

    int iter = 0;
    REAL error = 1.e10,error2=1.e10;
    int numeq = fCompMesh->NEquations();

    TPZFMatrix<STATE> prevsol ( fSolution );
    if ( prevsol.Rows() != numeq ) prevsol.Redim ( numeq,1 );

    if ( checkconv )
    {
        TPZVec<REAL> coefs ( 1,1. );
        TPZFMatrix<STATE> range ( numeq,1,1. );
        CheckConvergence ( *this,fSolution,range,coefs );
    }
    bool a=true,b=true,c=true;
    while ( a && ( b || c ) && iter<numiter )
    {

//        fSolution.Redim(0,0);
        Assemble();
        Solve();
        if ( linesearch )
        {
            TPZFMatrix<STATE> nextSol;
            REAL LineSearchTol = 1e-3 * Norm ( fSolution );
            const int niter = 10;
            this->LineSearch ( prevsol, fSolution, nextSol, LineSearchTol, niter );
            fSolution = nextSol;
        }
        else
        {
            TPZFMatrix<STATE> sol = fSolution;
            sol += prevsol;
        }

        prevsol -= fSolution;
        REAL normDeltaSol = Norm ( prevsol );
        prevsol = fSolution;
        this->LoadSolution ( fSolution );
        this->AssembleResidual();
        double NormResLambda = Norm ( fRhs );
        double norm = NormResLambda;
        cout << "ttttttttt Iteracao n : " << ( iter+1 ) << " : normas |Delta(Un)| e |Delta(rhs)| : " << normDeltaSol << " / " << NormResLambda << endl;
        a = iter < numiter ;
        b =error2 > tol*1.e-3;
        c= error > tol;

        if ( ( normDeltaSol - error2 ) > 1.e-9 && ( NormResLambda - error ) > 1.e-9 &&  normDeltaSol>0.1 )
        {
            out << "\nDivergent Method\n";
            return;
        }
        error = norm;
        error2=normDeltaSol;
        iter++;
        out.flush();
    }
}


REAL TPZElastoPlasticAnalysis::IterativeProcessArcLength ( REAL tol,int numiter,REAL tol2,int numiter2,REAL l,REAL lambda0,bool &converge )
{

    std::vector<double> fslist;
    REAL diff=1000.;
    int counterout=0;
    plasticmat * material= dynamic_cast<plasticmat *> ( fCompMesh->FindMaterial ( 1 ) );
    REAL lambda=lambda0;
    REAL lambdan=1000;
    TPZFMatrix<REAL> dw,dws,dwb,dww,displace,displacen;
    dw=fSolution;
    dw.Zero();
    cout << " a   " <<endl;
    material->SetWhichLoadVector ( 2 );
    material->SetLoadFactor ( 1. );
    AssembleResidual();
    TPZFMatrix<REAL> FBODY=Rhs();
    dw=Solution();
    displace=Solution();
    dw.Zero();
    displace.Zero();
    REAL normfbody=Norm ( FBODY );
    cout << "normfbody = " << normfbody << endl;
    displacen=displace;
    do//while(counterout<numiter && diff>tol);
    {
        cout << " load step  = " << counterout+1 << " load factor  = " << lambda << " diff = " << diff << " l = " << l <<  endl;
        int counter=0;
        dw.Zero();
        bool conv=true;
        REAL residualrhs=10.;
        REAL normrhsn=10.e12;
        REAL normdu=10.;
        REAL diffnorm=0.;

        do//while( counter<numiter2 && normdu>tol2 );
        {
            material->SetWhichLoadVector ( 0 );
            material->SetLoadFactor ( lambda );
            Assemble();
            this->fSolver->Solve ( fRhs,dws );
            TPZFMatrix<REAL> residual = Rhs();
            normrhsn=residualrhs;
            residualrhs=Norm ( residual );
            diffnorm=residualrhs-normrhsn;
            this->fSolver->Solve ( FBODY,dwb );
            REAL dlamb=0.;
            if ( counter==0 )
            {
                dlamb= computelamda0 ( dwb, dw,l );
            }
            else
            {
                dlamb= computelamda ( dwb, dws,  dw,  l );
            }

            lambda += dlamb;
            lambda0=lambda;
            dww = dwb*dlamb+dws;
            dw += dww;
            displace+=dww;
            normdu=Norm ( dww ) ;
            cout << " counter  = " << counter+1 <<" normrhs = " <<residualrhs << " normdu = " << normdu << " lambda = "<< lambda << " dlamb = "<< dlamb << " l = " << l << endl;
            if ( ( diffnorm>1000 && counter>3 ) || normdu>100. )
            {

                conv=false;
                break;
            }

            counter++;
            LoadSolution ( displace );

//             if(fabs(dlamb)<1.e-5 && residualrhs<1.e-3)
//             {
//                 cout << " fabs(dlamb)<1.e-5 && normdu<1.e-3 : conv=true" <<endl;
//                 conv=true;
//                 break;
//             }
        }
        while ( counter<numiter2 && normdu>tol2 );

        if ( conv==false )
        {
            if ( l<1.e-3 ||normdu>1.e7)
            {
                cout << "  l= "<< l <<  " normdu =  " << normdu <<endl;
                converge =false;
                return 0;
            }
            cout <<" changing initial paramenters... " <<endl;
            fSolution.Zero();
            fRhs.Zero();
            displace.Zero();
            LoadSolution ( displacen );
            l*=0.5;
        }
        else
        {
            diff=fabs ( lambda-lambdan );
            for ( int ifs=0; ifs<fslist.size(); ifs++ )
            {
                if ( fslist[ifs]>lambda )
                {
                    diff=1000.;
                }
            }

            AcceptSolution();
            int ndesi=10;
            l*=REAL ( ndesi ) /counter;
            if ( l>1. )
            {
                l=1.;
            }
            fslist.push_back ( lambda );

            lambdan=lambda;
            displacen=displace;
        }

        if ( counterout>=numiter )
        {
            cout <<"-----Failed to converge! "<< " residualrhs  = " <<residualrhs <<" diffnorm = " <<diffnorm  <<endl;
            converge=false;
            AcceptSolution();
            return lambda;
        }


        counterout++;



    }
    while ( counterout<numiter && diff>tol );

    AcceptSolution();
    cout << " final safety factor = "<<lambda << endl;

    converge =true;
    return lambda;

}

// REAL TPZElastoPlasticAnalysis::IterativeProcessArcLength(REAL tol,int numiter,REAL tol2,int numiter2,REAL l,REAL lambda0){
//
//     REAL diff=1000.;
//     std::ofstream eout("debug.txt");
//     std::ofstream eout2("debug2.nb");
//
//     int counterout=0;
//
//     plasticmat * material= dynamic_cast<plasticmat *> ( fCompMesh->FindMaterial ( 1 ) );
//
//     REAL lambda=lambda0;
//     REAL lambdan=lambda;
//
//
//     TPZFMatrix<REAL> dw,dws,dwb,dww,displace,displacen;
//     dw=fSolution;
//     dw.Zero();
//
//
//     cout << " a   " <<endl;
//     material->SetWhichLoadVector(2);
//     material->SetLoadFactor(1.);
//     AssembleResidual();
//     TPZFMatrix<REAL> FBODY=Rhs();
//     dw=Solution();
//     displace=Solution();
//     dw.Zero();
//     displace.Zero();
//
//     REAL normfbody=Norm(FBODY);
//
//     cout << "normfbody = " << normfbody << endl;
//
//     displacen=displace;
//     do
//     {
//
//         cout << " load step  = " << counterout+1 << " load factor  = " << lambda << " diff = " << diff << " l = " << l <<  endl;
//         int counter=0;
//
//         dw.Zero();
//
//         bool conv=true;
//
//         REAL residualrhs=10.;
//         REAL normrhsn=10.e12;
//         REAL normdu=10.;
//         REAL diffnorm=0.;
//
//         do
//         {
//             chrono::steady_clock sc;
//
//             //auto start = sc.now();
//             material->SetWhichLoadVector(0);
//             material->SetLoadFactor(lambda);
//
//             Assemble();
//
//             //TPZAutoPointer<TPZMatrix<REAL> > KG = this->fSolver->Matrix();
//
//             //fRhs.Print(std::cout);
//
//             //KG->Print(std::cout);
//             //auto end = sc.now();
//             //auto time_span = static_cast<chrono::duration<double>> ( end - start );
//             //cout << "| time to assemble=  " << time_span.count()<< std::endl;
//
//             //start = sc.now();
//             this->fSolver->Solve(fRhs,dws);
//             //end = sc.now();
//             //time_span = static_cast<chrono::duration<double>> ( end - start );
//             //cout << "| time to solve=  " << time_span.count()<< std::endl;
//
//
//             TPZFMatrix<REAL> residual = Rhs();
//
//             normrhsn=residualrhs;
//
//             residualrhs=Norm(residual);
//
//             diffnorm=residualrhs-normrhsn;
//
//             this->fSolver->Solve(FBODY,dwb);
//
//             REAL dlamb=0.;
//             if(counter==0)
//             {
//                 dlamb= computelamda0 ( dwb, dw,l );
//             }else{
//                  dlamb= computelamda ( dwb, dws,  dw,  l );
//             }
//
//             lambda += dlamb;
//
//             dww = dwb*dlamb+dws;
//
//             dw += dww;
//
//             displace+=dww;
//
//              normdu=Norm(dww) ;
//
//             cout << " counter  = " << counter+1 <<" normrhs = " <<residualrhs << " normdu = " << normdu << " lambda = "<< lambda << " dlamb = "<< dlamb << " l = " << l << endl;
//
//             if( (diffnorm>100 && counter>3) || normdu>10.)
//             {
//                 cout <<"-----Failed to converge! "<< " residualrhs  = " <<residualrhs <<" diffnorm = " <<diffnorm  <<endl;
//                 conv=false;
//                 break;
//             }
//
//             counter++;
//
//             LoadSolution(displace);
//
//         }while( counter<numiter2 && normdu>tol2 );
//
//         if(conv==false)
//         {
//
//             if(l<1.e-6)break;
//             lambda=lambdan;
//             fSolution.Zero();
//             fRhs.Zero();
//             displace.Zero();
//             LoadSolution(displacen);
//             l*=0.5;
//         }else{
//             diff=fabs(lambda-lambdan);
//             //LoadSolution(displace);
//             AcceptSolution();
//             REAL delta=1.e-6;
//             //lambda-=1.1*tol;
//             if(lambda<lambdan)
//             {
//                 lambda=lambdan;
//             }else{
//                 lambdan=lambda;
//             }
//             int ndesi=10;
//             l*=REAL(ndesi)/counter;
//             displacen=displace;
//
//             if(l>3.)
//             {
//                 l=3.;
//             }
//         }
//
//
//
//         counterout++;
//
//     }while(counterout<numiter && diff>tol);
//
//     AcceptSolution();
//     cout << "lambda = "<<lambdan << endl;
//
//     return lambdan;
//
// }



#include "pzintel.h"
//#include "pzelctempplus.h"

#include "pzrefpoint.h"
#include "pzgeopoint.h"
#include "pzshapepoint.h"
#include "tpzpoint.h"

#include "pzshapelinear.h"
#include "TPZGeoLinear.h"
#include "TPZRefLinear.h"
#include "tpzline.h"

#include "pzshapetriang.h"
#include "pzreftriangle.h"
#include "pzgeotriangle.h"
#include "tpztriangle.h"

#include "pzrefquad.h"
#include "pzshapequad.h"
#include "pzgeoquad.h"
#include "tpzquadrilateral.h"

#include "pzshapeprism.h"
#include "pzrefprism.h"
#include "pzgeoprism.h"
#include "tpzprism.h"

#include "pzshapetetra.h"
#include "pzreftetrahedra.h"
#include "pzgeotetrahedra.h"
#include "tpztetrahedron.h"

#include "pzshapepiram.h"
#include "pzrefpyram.h"
#include "pzgeopyramid.h"
#include "tpzpyramid.h"

#include "TPZGeoCube.h"
#include "pzshapecube.h"
#include "TPZRefCube.h"
#include "tpzcube.h"

#include "pzelctemp.h"


void TPZElastoPlasticAnalysis::SetAllCreateFunctionsWithMem ( TPZCompMesh *cmesh )
{

    TPZManVector<TCreateFunction,10> functions ( 8 );
    functions[EPoint] = &TPZElastoPlasticAnalysis::CreatePointElWithMem;
    functions[EOned] = TPZElastoPlasticAnalysis::CreateLinearElWithMem;
    functions[EQuadrilateral] = TPZElastoPlasticAnalysis::CreateQuadElWithMem;
    functions[ETriangle] = TPZElastoPlasticAnalysis::CreateTriangElWithMem;
    functions[EPrisma] = TPZElastoPlasticAnalysis::CreatePrismElWithMem;
    functions[ETetraedro] = TPZElastoPlasticAnalysis::CreateTetraElWithMem;
    functions[EPiramide] = TPZElastoPlasticAnalysis::CreatePyramElWithMem;
    functions[ECube] = TPZElastoPlasticAnalysis::CreateCubeElWithMem;
    cmesh->ApproxSpace().SetCreateFunctions ( functions );
}

TPZCompEl * TPZElastoPlasticAnalysis::CreateCubeElWithMem ( TPZGeoEl *gel, TPZCompMesh &mesh, long &index )
{
    return new TPZCompElWithMem< TPZIntelGen< pzshape::TPZShapeCube > > ( mesh,gel,index );
}

TPZCompEl * TPZElastoPlasticAnalysis::CreateLinearElWithMem ( TPZGeoEl *gel, TPZCompMesh &mesh, long &index )
{
    return new TPZCompElWithMem< TPZIntelGen< pzshape::TPZShapeLinear > > ( mesh,gel,index );
}

TPZCompEl * TPZElastoPlasticAnalysis::CreatePointElWithMem ( TPZGeoEl *gel, TPZCompMesh &mesh, long &index )
{
    return new TPZCompElWithMem< TPZIntelGen< pzshape::TPZShapePoint > > ( mesh,gel,index );
}

TPZCompEl * TPZElastoPlasticAnalysis::CreatePrismElWithMem ( TPZGeoEl *gel, TPZCompMesh &mesh, long &index )
{
    return new TPZCompElWithMem< TPZIntelGen< pzshape::TPZShapePrism > > ( mesh,gel,index );
}

TPZCompEl * TPZElastoPlasticAnalysis::CreatePyramElWithMem ( TPZGeoEl *gel, TPZCompMesh &mesh, long &index )
{
    return new TPZCompElWithMem< TPZIntelGen< pzshape::TPZShapePiram > > ( mesh,gel,index );
}

TPZCompEl * TPZElastoPlasticAnalysis::CreateQuadElWithMem ( TPZGeoEl *gel, TPZCompMesh &mesh, long &index )
{
//	return new TPZCompElWithMem< TPZIntelGenPlus<TPZIntelGen< pzshape::TPZShapeQuad > > >(mesh,gel,index);
    return new TPZCompElWithMem< TPZIntelGen< pzshape::TPZShapeQuad > > ( mesh,gel,index );
}

TPZCompEl * TPZElastoPlasticAnalysis::CreateTetraElWithMem ( TPZGeoEl *gel, TPZCompMesh &mesh, long &index )
{
    return new TPZCompElWithMem< TPZIntelGen< pzshape::TPZShapeTetra > > ( mesh,gel,index );
}

TPZCompEl * TPZElastoPlasticAnalysis::CreateTriangElWithMem ( TPZGeoEl *gel, TPZCompMesh &mesh, long &index )
{
    return new TPZCompElWithMem< TPZIntelGen< pzshape::TPZShapeTriang > > ( mesh,gel,index );
}

void TPZElastoPlasticAnalysis::IdentifyEquationsToZero()
{
    fEquationstoZero.clear();
    long nel = fCompMesh->NElements();
    for ( long iel=0; iel<nel; iel++ )
    {
        TPZCompEl *cel = fCompMesh->ElementVec() [iel];
        if ( !cel )
        {
            continue;
        }
        TPZMaterial *mat = cel->Material();
        if ( !mat )
        {
            continue;
        }
        int matid = mat->Id();
        if ( fMaterialIds.find ( matid ) == fMaterialIds.end() )
        {
            continue;
        }
        std::pair<std::multimap<int, int>::iterator,std::multimap<int, int>::iterator> ret;
        ret = fMaterialIds.equal_range ( matid );
        std::multimap<int, int>::iterator it;
        for ( it=ret.first; it != ret.second; it++ )
        {
            int direction = it->second;
            long nc = cel->NConnects();
            for ( long ic=0; ic<nc; ic++ )
            {
                TPZConnect &c = cel->Connect ( ic );
                long seqnum = c.SequenceNumber();
                long pos = fCompMesh->Block().Position ( seqnum );
                int blsize = fCompMesh->Block().Size ( seqnum );
                for ( long i=pos+direction; i<pos+blsize; i+=2 )
                {
                    fEquationstoZero.insert ( i );
                }
            }
        }
    }
#ifdef LOG4CXX
    {
        if ( EPAnalysisLogger->isDebugEnabled() )
        {
            std::stringstream sout;
            sout << "Equations to zero ";
            std::set<long>::iterator it;
            for ( it=fEquationstoZero.begin(); it!= fEquationstoZero.end(); it++ )
            {
                sout << *it << " ";
            }
            LOGPZ_DEBUG ( EPAnalysisLogger, sout.str() )
        }
    }
#endif
}

/// return the vector of active equation indices
void TPZElastoPlasticAnalysis::GetActiveEquations ( TPZVec<long> &activeEquations )
{
    long neq = fCompMesh->NEquations();
    TPZVec<int> equationflag ( neq,1 );
    typedef std::set<long>::iterator setit;
    for ( setit it = fEquationstoZero.begin(); it != fEquationstoZero.end(); it++ )
    {
        equationflag[*it] = 0;
    }
    activeEquations.resize ( neq-fEquationstoZero.size() );
    long count = 0;
    for ( long i=0; i<neq; i++ )
    {
        if ( equationflag[i]==1 )
        {
            activeEquations[count++] = i;
        }
    }
}


