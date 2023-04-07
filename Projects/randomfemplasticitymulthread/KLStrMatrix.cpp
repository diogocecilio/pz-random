#include "KLStrMatrix.h"
#include "pzfstrmatrix.h"
#include "pzfmatrix.h"
#include "pzcmesh.h"
#include "pzsubcmesh.h"
#include <sstream>
#include "pzlog.h"
#include "pzfmatrix.h"
#include "pzmatrix.h"
#include "pzsolve.h"
#ifdef LOG4CXX
static LoggerPtr logger ( Logger::getLogger ( "pz.strmatrix.tpzfstructmatrix" ) );
static LoggerPtr loggerel ( Logger::getLogger ( "pz.strmatrix.element" ) );
#endif


using namespace std;

KLStrMatrix::KLStrMatrix ( TPZCompMesh *mesh ) : TPZStructMatrix ( mesh )
{

}

KLStrMatrix::KLStrMatrix ( TPZAutoPointer<TPZCompMesh> cmesh ) :  TPZStructMatrix ( cmesh )
{

}


TPZMatrix<STATE> * KLStrMatrix::Create()
{

}

TPZMatrix<STATE> * KLStrMatrix::CreateAssemble ( TPZFMatrix<STATE> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface )
{

}

void KLStrMatrix::AssembleC (TPZFMatrix<REAL> &C)
{
  long iel, jel;
  long nelem = fMesh->NElements();
  TPZManVector<TPZManVector<int>> meshtopology;

  TPZElementMatrix temp ( fMesh, TPZElementMatrix::EK );
  TPZAdmChunkVector<TPZCompEl *> &elementvec = fMesh->ElementVec();
  int sz = fMesh->NEquations();
  C.Resize ( sz, sz );

  long count = 0;
  for ( iel = 0; iel < nelem; iel++ )
    {
      //std::cout << "\n iel " << iel << std::endl;
      TPZCompEl *el = elementvec[iel];
      for ( int jel = 0; jel < nelem; jel++ )
        {
          //std::cout << "\n jel " << jel << std::endl;
          TPZCompEl *elj = elementvec[jel];
          TPZElementMatrix ce ( fMesh, TPZElementMatrix::EK );
          el->CalcStiffC (elj, ce );
          //ce.fMat.Print(std::cout);
          int nshape = ce.fMat.Rows();
          TPZManVector<long> SourceIndexIEL, DestinationIndexIEL, SourceIndexJEL, DestinationIndexJEL;
          GetDestIndex ( iel, nshape, SourceIndexIEL, DestinationIndexIEL );
          GetDestIndex ( jel, nshape, SourceIndexJEL, DestinationIndexJEL );

          for ( int irow = 0; irow < DestinationIndexIEL.size(); irow++ )
            {
              for ( int icol = 0; icol < DestinationIndexJEL.size(); icol++ )
                {
                  C ( DestinationIndexIEL[irow], DestinationIndexJEL[icol] ) += ce.fMat ( SourceIndexIEL[irow], SourceIndexJEL[icol] );
                }
            }

        }//jel

    }//iel

 // C.Print ( std::cout );
}
void KLStrMatrix::AssembleB(TPZFMatrix<REAL> &B)
{
  long iel;
  long nelem = fMesh->NElements();
  TPZManVector<TPZManVector<int>> meshtopology;

  TPZElementMatrix temp ( fMesh, TPZElementMatrix::EK );
  TPZAdmChunkVector<TPZCompEl *> &elementvec = fMesh->ElementVec();
  int sz = fMesh->NEquations();
  B.Resize ( sz, sz );

  for ( iel = 0; iel < nelem; iel++ )
    {
      TPZCompEl *el = elementvec[iel];
      TPZElementMatrix be ( fMesh, TPZElementMatrix::EK );
      el->CalcStiffB ( be );
      int nshape = be.fMat.Rows();
      TPZManVector<long> SourceIndexIEL, DestinationIndexIEL;
      GetDestIndex ( iel, nshape, SourceIndexIEL, DestinationIndexIEL );
      for ( int irow = 0; irow < DestinationIndexIEL.size(); irow++ )
          {
            for ( int icol = 0; icol < DestinationIndexIEL.size(); icol++ )
              {
                B ( DestinationIndexIEL[irow], DestinationIndexIEL[icol] ) += be.fMat ( SourceIndexIEL[irow], SourceIndexIEL[icol] );
              }
          }
    }//iel

 // C.Print ( std::cout );
}
// void KLStrMatrix::AssembleB(TPZFMatrix<REAL> &B)
// {
// 
// 
//   long iel;
//   long nelem = fMesh->NElements();
//   TPZElementMatrix ek ( fMesh, TPZElementMatrix::EK ), ef ( fMesh, TPZElementMatrix::EF );
//   int sz = fMesh->NEquations();
//   B.Resize ( sz, sz );
//   TPZAdmChunkVector<TPZCompEl *> &elementvec = fMesh->ElementVec();
// 
//   long count = 0;
//   for ( iel = 0; iel < nelem; iel++ )
//     {
//       TPZCompEl *el = elementvec[iel];
//       ek.Reset();
//       ef.Reset();
//       el->CalcStiffB ( ek );
// 
//       if ( !ek.HasDependency() )
//         {
//           ek.ComputeDestinationIndices();
//           //fEquationFilter.Filter(ek.fSourceIndex, ek.fDestinationIndex);
//           B.AddKel ( ek.fMat, ek.fSourceIndex, ek.fDestinationIndex );
// 
//         }
//       else
//         {
//           // the element has dependent nodes
//           ek.ApplyConstraints();
//           ef.ApplyConstraints();
//           ek.ComputeDestinationIndices();
//           //fEquationFilter.Filter(ek.fSourceIndex, ek.fDestinationIndex);
//           B.AddKel ( ek.fConstrMat, ek.fSourceIndex, ek.fDestinationIndex );
//         }
//     }//fim for iel
// 
// 
//  // B.Print ( std::cout );
// }

void KLStrMatrix::GetDestIndex ( int iel, int nshape, TPZManVector<long> &fSourceIndex, TPZManVector<long> &fDestinationIndex )
{
  fSourceIndex.Resize ( nshape );
  fDestinationIndex.Resize ( nshape );
  long destindex = 0L;
  long fullmatindex = 0L;
  const int numnod = fMesh->ElementVec() [iel]->NConnects();
  for ( int in = 0; in < numnod; in++ )
    {
      const long npindex = fMesh->ElementVec() [iel]->ConnectIndex ( in );
      TPZConnect &np = fMesh->ConnectVec() [npindex];
      long blocknumber = np.SequenceNumber();
      long firsteq = fMesh->Block().Position ( blocknumber );
      int ndf = fMesh->Block().Size ( blocknumber );
      if ( np.HasDependency() || np.IsCondensed() )
        {
          fullmatindex += ndf;
          continue;
        }//for (np)
      for ( int idf = 0; idf < ndf; idf++ )
        {
          fSourceIndex[destindex] = fullmatindex++;
          fDestinationIndex[destindex++] = firsteq + idf;
        }//for idf
    }//for in

}
void KLStrMatrix::Solve()
{
//      TPZFMatrix<REAL> invB,invBC,B,C;
//      MatrixXd eigenInvBC,eigenInvB,eigenC,eigenB;
//      AssembleC(C);
//      AssembleB(B);
//      
// 
//     ToEigen(B,eigenB);
//     ToEigen(C,eigenC);
//     
//     eigenInvB = eigenB.inverse();
//     std::cout << eigenInvB <<std::endl;
//     eigenInvBC = eigenInvB*eigenC;
// 
// 
//     MatrixXd val,vec;
//     
//     ComplexEigenSolver<MatrixXd> ces;
//     ces.compute(eigenInvBC);
//     cout << "The eigenvalues of A are:" << endl << ces.eigenvalues() << endl;
//     cout << "The matrix of eigenvectors, V, is:" << endl << ces.eigenvectors() << endl << endl;
//  
//     
//     double lambda = ces.eigenvalues()[0].real();
//     int ncols = ces.eigenvectors().cols();
//     int nrows = ces.eigenvectors().rows();
//     val.resize(nrows,1);
//     vec.resize(nrows,ncols);
//     for(int i=0;i<nrows;i++)
//     {
//       val(i,0)=ces.eigenvalues()[i].real();
//       for(int j =0; j< ncols;j++)
//       {
//         vec(i,j)=ces.eigenvectors()(i,j).real();
//       }
//       
//     }
// 
//   cout << "The eigenvalues of A are:" << endl << val << endl;
//   cout << "The matrix of eigenvectors, V, is:" << endl << vec << endl << endl;
//   cout << "... and A * v = " << endl << eigenInvBC * vec << endl << endl;
// 
//   cout << "Finally, V * D * V^(-1) = " << endl
//        << vec * val.asDiagonal() * vec.inverse() << endl;
}

TPZStructMatrix * KLStrMatrix::Clone()
{
  return new KLStrMatrix ( *this );
}
