/**
 * @file
 * @brief Contains the TPZFStructMatrix class which implements Full Structural Matrices.
 */

#ifndef KLStrMatrix_H
#define KLStrMatrix_H

#include "pzstrmatrix.h"
#include <Eigen/Core>

using namespace Eigen;
#include <Eigen/Dense>
 
using namespace std;
using namespace Eigen;
class TPZCompMesh;
template<class TVar>
class TPZFMatrix;
template<class TVar>
class TPZMatrix;

/**
 * @brief Implements Full Structural Matrices. \ref structural "Structural Matrix"
 * @ingroup structural
 */
class KLStrMatrix : public TPZStructMatrix {
public:    
	
    KLStrMatrix(TPZCompMesh *);
	
    KLStrMatrix(TPZAutoPointer<TPZCompMesh> );
	
	void AssembleC (TPZFMatrix<REAL> &C);
	
	void AssembleB(TPZFMatrix<REAL> &B);
	
	void Solve();
    
    virtual TPZMatrix<STATE> * Create();
	
    virtual TPZMatrix<STATE> * CreateAssemble(TPZFMatrix<STATE> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
	
    virtual TPZStructMatrix * Clone();
	
	void GetDestIndex(int iel, int nshape, TPZManVector<long> &fSourceIndex, TPZManVector<long> &fDestinationIndex);
	void GetDestIndex (  int iel, TPZElementMatrix *elmat, TPZManVector<long> &fSourceIndex, TPZManVector<long> &fDestinationIndex );
	void SetMesh(TPZCompMesh * cmesh)
	{
		fCompMesh =cmesh;
	}


		 void FromEigen ( MatrixXd eigenmat, TPZFMatrix<REAL>  &pzmat )
    {

		int rows = eigenmat.rows();
		int cols = eigenmat.cols();
		pzmat.Resize(rows,cols);
		for(int irow=0;irow<rows;irow++)
		{
			for(int icol=0;icol<cols;icol++)
			{
				pzmat(irow,icol)=eigenmat(irow,icol);
			}
		}
	
    }

	void ToEigen ( TPZFMatrix<REAL>  pzmat,MatrixXd &eigenmat )
    {
		TPZFMatrix<REAL> intpz(pzmat);
		int rows = pzmat.Rows();
		int cols = pzmat.Cols();
		eigenmat.resize(rows,cols);
		for(int irow=0;irow<rows;irow++)
		{
			for(int icol=0;icol<cols;icol++)
			{
				eigenmat(irow,icol)=intpz(irow,icol);
			}
		}
	
    }
	
private:
	
	TPZFMatrix<REAL> fB;
	TPZFMatrix<REAL> fC;
	

};

#endif //TPZFSTRUCTMATRIX_H
