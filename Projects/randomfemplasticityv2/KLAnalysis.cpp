#include "KLAnalysis.h"

#include "pzstrmatrix.h"

#include "tpznodesetcompute.h"
#include "tpzsparseblockdiagonal.h"
#include "pzseqsolver.h"
#include "pzbdstrmatrix.h"
#include "TPZSkylineNSymStructMatrix.h"
#include <Eigen/Eigenvalues>

#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include <sstream>
#include <set>
#include <random>
#include <limits>
//#include "pzfmatrix.h"
#include "pzmatrix.h"

#include "pzvec.h"
#include "pzfmatrix.h"
#include "pzmatrix.h"
#include "pzsolve.h"
#include "fadType.h"
#include "fad.h"

#include <exception>
#include "pzlog.h"
#include <complex>
#include <pzinterpolationspace.h>
#include <tpzintpoints.h>
#include <chrono>
using namespace std;



KLAnalysis::KLAnalysis ( TPZCompMesh *mesh ) : TPZAnalysis ( mesh,true )
{

	OptimizeBandwidth();
    KLStrMatrix *strmat = new KLStrMatrix ( mesh );
    strmat->SetMesh ( mesh );
    fStrMatrix = strmat;
    int numeq = fCompMesh->NEquations();
    fSolution.Redim ( numeq,1 );
    fSolution.Zero();
    LoadSolution();
	

}


KLAnalysis::~KLAnalysis ( void )
{

}


void KLAnalysis::Assemble()
{
    
}

void KLAnalysis::Solve()
{
	std::ofstream outfull("outinfofulltime-p2-h4.txt");
	chrono::steady_clock fulltime;
	auto startfull = fulltime.now();
	std::ofstream out("outinfo-p2-h4.txt");
    TPZFMatrix<REAL> invB,invBC,B,C;
    MatrixXd eigenInvBC,eigenInvB,eigenC,eigenB;
	cout << "Number of Equations =   " << fCompMesh->NEquations() << std::endl;
    cout << "Assembling C  " << std::endl;
    chrono::steady_clock sc;
    auto start = sc.now();

    fStrMatrix->AssembleC ( C );

    //int nrows = C.Rows();
    //int ncols = C.Rows();

    auto end = sc.now();
    auto time_span = static_cast<chrono::duration<double>> ( end - start );
    cout << "| total time taken to assemble C =  " << time_span.count()<< std::endl;

   cout << " Assembling B  " << std::endl;
    start = sc.now();

    fStrMatrix->AssembleB ( B );

    end = sc.now();
    time_span = static_cast<chrono::duration<double>> ( end - start );
    cout << "| total time taken to assemble B =  " << time_span.count()<< std::endl;

    ToEigen ( B,eigenB );
    ToEigen ( C,eigenC );

    eigenInvB = eigenB.inverse();
    eigenInvBC = eigenInvB*eigenC;

    MatrixXd val,vec;

    //GeneralizedEigenSolver<MatrixXd> ces;

    //EigenSolver<MatrixXd> ces;
    ComplexEigenSolver<MatrixXd> ces;

    cout << " Computing Eigensystem  " << std::endl;
    start = sc.now();
    ces.compute ( eigenInvBC );
    //ces.compute(eigenB, eigenC);
    end = sc.now();
    time_span = static_cast<chrono::duration<double>> ( end - start );
   cout << "| total time taken to compute the Eigensystem =  " << time_span.count()<< std::endl;
	ces.info();
	
	auto endfull = fulltime.now();
	time_span = static_cast<chrono::duration<double>> ( endfull - startfull );
    cout << "| total time  =  " << time_span.count()<< std::endl;
	
    int ncols = ces.eigenvectors().cols();
    int M = fExpansionOrder;
    //int M = ncols;
    int nrows = ces.eigenvectors().rows();
    val.resize ( M,1 );
    vec.resize ( nrows,M );

    for ( int icol=0; icol< M; icol++ )
    {
        val ( icol,0 ) =ces.eigenvalues() [nrows-icol-1].real();
        for(int irow=0; irow<nrows; irow++)
        {
            vec ( irow,icol ) =ces.eigenvectors() ( irow,ncols-icol-1 ).real();
        }
    }

    bool print =true;
    if ( print==true )
    {
        cout << "Eigenvalues "<<endl;
        cout << val << endl;
        //cout << "Eigenvectors "<<endl;
        //cout << vec << endl;
        // std::cout << " A "<<endl;
        // std::cout << eigenInvBC <<std::endl;
        //std::cout << " A- V * D * V^(-1) = " << "\n" << (vec * val.asDiagonal() * vec.inverse()) << endl;
    }
    
    cout << "Start to integrate the solution over the domain..  " << std::endl;
    start = sc.now();
    TPZFMatrix<REAL> vecpz(nrows,M);
    fEigenVectors.Resize(M);
    VectorXd intphisqr(M);
    for(int i=0; i<M; i++)fEigenVectors[i].Resize(nrows,1);
    for ( int icol=0; icol< M; icol++ )
    {
        for(int irow=0; irow<nrows; irow++)
        {
            fEigenVectors[icol](irow,0)=vec.col(icol)(irow);
        }

        //fEigenVectors[icol].Print(std::cout);
        LoadSolution(fEigenVectors[icol]);
        int varid=6;
        REAL integral = IntegrateSolution(varid);
        intphisqr(icol)=integral;
        fEigenVectors[icol]*=1./sqrt(integral);
        vec.col(icol)*=1./sqrt(integral);
    }
    end = sc.now();
    time_span = static_cast<chrono::duration<double>> ( end - start );
    cout << "| total time taken to integrate the solution over the domain =  " << time_span.count()<< std::endl;

	//cout << intphisqr << endl;
    for ( int icol=0; icol< M; icol++ )
    {
        for(int irow=0; irow<nrows; irow++)
        {
            vecpz(irow,icol)=sqrt(val(icol))*vec.col(icol)(irow);
        }
    }
    LoadSolution(vecpz);
    REAL totalarea = ComputeTotalArea();
    REAL sumeigenvalues = 0.;
    REAL sumeigenvaluestimessqrvec=0.;
    for ( int i = 0; i < val.size(); i++ ) {
        sumeigenvalues +=  val(i)  ;
        sumeigenvaluestimessqrvec += val(i)* intphisqr(i);
    }


     cout << "total area = "<<totalarea << std::endl;
     cout << "sumeigenvalues= "<<sumeigenvalues << std::endl;
     cout << "sumeigenvaluestimessqrvec= "<<sumeigenvaluestimessqrvec << std::endl;
     //std::cout << " err / (fsig * fsig) = " <<err / (fsig * fsig) << std::endl;
     cout << "mean error 1 = " <<1. - sumeigenvalues/totalarea  << std::endl;
     cout << "mean error 2 = " <<1. - sumeigenvaluestimessqrvec/totalarea  << std::endl;

    FromEigen(val,fEigenValues);


}

// /*
// void KLAnalysis::Solve()
// {
// 
//     TPZFMatrix<REAL> invB,invBC,B,C;
//     MatrixXd eigenInvBC,eigenInvB,eigenC,eigenB;
// 
// 
//     std::cout << "Assembling C  " << std::endl;
//     chrono::steady_clock sc;
//     auto start = sc.now();
// 
//     fStrMatrix->AssembleC ( C );
// 
//     auto end = sc.now();
//     auto time_span = static_cast<chrono::duration<double>> ( end - start );
//     cout << "| total time taken to assemble C =  " << time_span.count()<< std::endl;
// 
// 
//     std::cout << " Assembling B  " << std::endl;
//     start = sc.now();
// 
//     fStrMatrix->AssembleB ( B );
// 
//     end = sc.now();
//     time_span = static_cast<chrono::duration<double>> ( end - start );
//     cout << "| total time taken to assemble B =  " << time_span.count()<< std::endl;
// 
//     ToEigen ( B,eigenB );
//     ToEigen ( C,eigenC );
// 
//     eigenInvB = eigenB.inverse();
//     //std::cout << eigenInvB <<std::endl;
//     eigenInvBC = eigenInvB*eigenC;
// 
//     MatrixXd val,vec;
// 
//     ComplexEigenSolver<MatrixXd> ces;
//     std::cout << " Computing Eigensystem  " << std::endl;
//     start = sc.now();
//     ces.compute ( eigenInvBC );
//     end = sc.now();
//     time_span = static_cast<chrono::duration<double>> ( end - start );
//     cout << "| total time taken to compute the Eigensystem =  " << time_span.count()<< std::endl;
// 
// 
//     int ncols = ces.eigenvectors().cols();
//     int M = fExpansionOrder;
//     //int M = ncols;
//     int nrows = ces.eigenvectors().rows();
//     val.resize ( M,1 );
//     vec.resize ( nrows,M );
// //   for ( int irow=0; irow<nrows; irow++ )
// //   {
// //     val ( irow,0 ) =ces.eigenvalues() [nrows-irow-1].real();
// //     for ( int icol=0; icol<M; icol++ )
// //       {
// //         vec ( irow,icol ) =ces.eigenvectors() ( irow,ncols-icol-1 ).real();
// //       }
// //   }
//     for ( int icol=0; icol< M; icol++ )
//     {
//         val ( icol,0 ) =ces.eigenvalues() [nrows-icol-1].real();
//         for(int irow=0; irow<nrows; irow++)
//         {
//             vec ( irow,icol ) =ces.eigenvectors() ( irow,ncols-icol-1 ).real();
//         }
//     }
// 
//     //std::cout <<"VT * D"   << "\n" << (vec* val.asDiagonal()* vec.inverse() ) << endl;
//     bool print =true;
//     if ( print==true )
//     {
//         std::cout << "Eigenvalues "<<endl;
//         std::cout << val << endl;
//         std::cout << "Eigenvectors "<<endl;
//         std::cout << vec << endl;
//         // std::cout << " A "<<endl;
//         // std::cout << eigenInvBC <<std::endl;
//         //std::cout << " A- V * D * V^(-1) = " << "\n" << (vec * val.asDiagonal() * vec.inverse()) << endl;
//     }
// 
//     std::cout << "Start to integrate the solution over the domain..  " << std::endl;
//     start = sc.now();
//     TPZFMatrix<REAL> vecpz(nrows,M);
//     fEigenVectors.Resize(M);
//     VectorXd intphisqr(M);
//     for(int i=0; i<M; i++)fEigenVectors[i].Resize(nrows,1);
//     for ( int icol=0; icol< M; icol++ )
//     {
//         for(int irow=0; irow<nrows; irow++)
//         {
//             fEigenVectors[icol](irow,0)=vec.col(icol)(irow);
//         }
// 
//         //fEigenVectors[icol].Print(std::cout);
//         LoadSolution(fEigenVectors[icol]);
//         int varid=6;
//         REAL integral = IntegrateSolution(varid);
//         intphisqr(icol)=integral;
//         fEigenVectors[icol]*=1./sqrt(integral);
//         vec.col(icol)*=1./sqrt(integral);
//     }
//     end = sc.now();
//     time_span = static_cast<chrono::duration<double>> ( end - start );
//     cout << "| total time taken to integrate the solution over the domain =  " << time_span.count()<< std::endl;
//     cout << " Integration factors \n";
//     cout<< intphisqr << endl;
// 
// 
// 
// 
//     for ( int icol=0; icol< M; icol++ )
//     {
//         for(int irow=0; irow<nrows; irow++)
//         {
//             vecpz(irow,icol)=fEigenVectors[icol](irow,0);
//         }
//     }
//     LoadSolution(vecpz);
// 
// 
//     /*	% check orthonormality
//     % So that \int_D phi_i(t)phi_j(t) dt = delta_ij
//     AA = ones(d_KL);
//     for i = 1:d_KL
//     for j = 1:d_KL
//         AA(i,j) = trapz(xnod, eigvec(:,i).*eigvec(:,j));
//     end
//     end
//     figure(66); spy(round(AA, 1))*/;
//     cout << " \n Check orthonormality, so that  int_D phi_i(t)phi_j(t) dt = delta_ij "<< std::endl;
//     MatrixXd vecsqr = vec*vec;
//     MatrixXd integra2(M,M);
//     TPZVec<TPZFMatrix<REAL>> vecsqrpz(M);
//     for(int i=0; i<M; i++)vecsqrpz[i].Resize(nrows,1);
// 	
// //cout << "vtv" <<endl;
// //cout << vtv << endl;
//     for ( int icol=0; icol< M; icol++ )
//     {
//         for ( int jcol=0; jcol< M; jcol++ )
//         {
//             for(int irow=0; irow<nrows; irow++)
//             {
//                 vecsqrpz[icol](irow,0)=vec.col(icol)(irow)*vec.col(jcol)(irow);
//             }
//             LoadSolution(vecsqrpz[icol]);
//             int varid=0;
//             REAL integral = IntegrateSolution(varid);
// 			if(integral<1.e-10)integral=0;
//             integra2(icol,jcol)=integral;
//         }
//         //fEigenVectors[icol].Print(std::cout);
//       // LoadSolution(vecsqrpz[icol]);
//       //  int varid=0;
//      //   REAL integral = IntegrateSolution(varid);
//      //   intphisqr(icol)=integral;
// 
//         //vecsqrpz[icol]*=1./integral;
//     }
//     VectorXd vtv = (vec.transpose() * vec);
//     //cout << "vtv" <<endl;
// //cout << vtv << endl;
// 
//     
//     
// 
// IOFormat CommaInitFmt(StreamPrecision, DontAlignCols, ", ", ", ", "", "", " << ", ";");
// IOFormat CleanFmt(3, 0, ", ", "\n", "[", "]");
// IOFormat OctaveFmt(StreamPrecision, 0, ", ", ";\n", "", "", "[", "]");
// IOFormat HeavyFmt(FullPrecision, 0, ", ", ";\n", "[", "]", "[", "]");
//  
// // std::cout << m1 << sep;
// // std::cout << m1.format(CommaInitFmt) << sep;
// // std::cout << m1.format(CleanFmt) << sep;
// // std::cout << m1.format(OctaveFmt) << sep;
// // std::cout << m1.format(HeavyFmt) << sep;
// 
//     cout << integra2 <<endl;
// 	 //   cout << " \n Check orthonormality, so that  int_D phi_i(t)phi_j(t) dt = delta_ij "<< std::endl;
// 	//cout << intphisqr <<endl;
//     //vecsqrpz[0].Print(std::cout);
// 
// 
//     LoadSolution(vecpz);
//     REAL err = 0.;
//     REAL err2=0.;
//     for ( int i = 0; i < val.size(); i++ ) {
//        // err += fabs ( val(i)*intphisqr(i) ) ;
// 		err += fabs ( val(i)* intphisqr(i));
//         err2 += fabs ( val(i) ) ;
//     }
// 
//     REAL totalarea = ComputeTotalArea();
//     std::cout << "total area = "<<totalarea << std::endl;
//     //std::cout << " err / (fsig * fsig) = " <<err / (fsig * fsig) << std::endl;
//     std::cout << "mean error 1 = " <<1. - 1./totalarea *   err2 << std::endl;
//     std::cout << "mean error 2 = " <<1. - 1./totalarea *   err << std::endl;
// 
//     FromEigen(val,fEigenValues);
// //
// // 	vecpz.Print(std::cout);
// 
// }*/
REAL KLAnalysis::IntegrateSolution(int varid)
{


    //int varid=6;//EVECSQR
    //int varid=0;//EVEC
    int nvar=1;
    int nels = fCompMesh->NElements(), iel;
    REAL sum=0.;
    TPZVec<REAL> vecsol ( nvar );


    for ( iel = 0; iel < nels; iel++ )
    {
        TPZCompEl * cel = fCompMesh->Element ( iel );

        vecsol = cel->IntegrateSolution ( varid );
        sum+=vecsol[0];
    }

    return  sum ;

}
void KLAnalysis::SetSolver ( TPZMatrixSolver<STATE> &solver )
{
    if ( fSolver ) delete fSolver;
    fSolver = ( TPZMatrixSolver<STATE> * ) solver.Clone();
}
/// Compute the area of the domain
REAL KLAnalysis::ComputeTotalArea()
{
	int dim = fCompMesh->Dimension();
    REAL area = 0.;
    long nelem = fCompMesh->NElements();
    TPZMaterial *pMatWithMem2 = fCompMesh->MaterialVec()[1];

    if (!pMatWithMem2) {
    }
    else
    {
        for (long el = 0; el<nelem; el++) {
            TPZCompEl *cel = fCompMesh->ElementVec()[el];
            if (!cel) {
                continue;
            }
            TPZGeoEl *gel = cel->Reference();
            if (gel->MaterialId() != 1) {
                continue;
            }
            TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
            if (!intel) {
                DebugStop();
            }
            TPZIntPoints &rule = intel->GetIntegrationRule();
            int np = rule.NPoints();
            for (int ip = 0; ip<np; ip++) {
                TPZManVector<REAL> point(dim,0.);
                REAL weight;
                rule.Point(ip, point, weight);
                TPZFMatrix<REAL> jac,jacinv;
                TPZFMatrix<REAL> axes;
                REAL detjac;
                gel->Jacobian(point, jac, axes, detjac, jacinv);
                area += weight*fabs(detjac);
            }
        }
    }
    return area;
}
