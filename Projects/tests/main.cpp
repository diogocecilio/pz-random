
#include <iostream>
//#include <stdlib.h>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include "mkl_pardiso.h"
#include "mkl.h"

enum MSystemType {ESymmetric, EHermitian, ENonSymmetric};
enum MStructure {EStructureSymmetric, EStructureNonSymmetric};
enum MProperty {EPositiveDefinite, EIndefinite};
enum ArrayType{EVector, ESparse};

struct ParConfig{

    MSystemType fSystemType = ESymmetric;
    MStructure fStructure = EStructureSymmetric;
    MProperty fProperty = EIndefinite;

    std::vector<long long>  fIA;
    std::vector<long long>  fJA;
    std::vector<double> fA;
    std::vector<double> frhs;

    long long *fPardisoControl[64];

    // adress of the first element of pt;
    void *fHandle[64];
    //long long *fHandle;
    //_MKL_DSS_HANDLE_t fHandle;

    // Array used to pass parameters to Pardiso
    long long fParam[64]; //setsize_64
    //MKL_INT *fParam;

    // Maximum number of factors we will pass to the solver
    long long fMax_num_factors = 1;

    // Factor number we are using
    long long int fMatrix_num = 1;

    // Message level information
    long long fMessageLevel=1;

    // error flag from Pardiso
    long long fError;

    /// permutation vector computed by Pardiso
    std::vector<long long> fPermutation;

    // matrix type, computed based on the structural information and double
    long long fMatrixType;

    /// Compute the matrix type
    //long long MatrixType();

    /// Provides a explanation for the given error
    void Error_check(int error) const;
};

template <class TVar>
void PrintArray(TVar *it,long long size,int n=5, std::string arrayName="vector"/*ArrayType=EVector*/);
void FilliParm(ParConfig &parconf,bool usingCGS);
template <class TVar>
void CheckConsistency(std::vector<TVar> vector, std::string name );

int main(int argc, const char * argv[]) {

    std::string pathtofileA="/home/diogo/Downloads/fMatrix.txt";
    std::string pathtofilerhs="/home/diogo/Downloads/frhs.txt";
    bool usingCGS=true;

    ParConfig parconf;
    long long parCtrl[64];
    for(int i=0;i<64;i++){
       parCtrl[i]=0;
       parconf.fParam[i]=0;

    }

    FilliParm(parconf,usingCGS);

    for (int i = 0; i < 64; i++ )
        {
            parconf.fHandle[i] = 0;
        }

    //Read matrix A
    std::ifstream matrixFile;
    matrixFile.open(pathtofileA);

    std::cout << "Start reading matrix A\n";
    int lineCounter = 0;
    long long ARows = -1;
    std::string dataLine;
    while (std::getline(matrixFile,dataLine)){
        if (lineCounter < 7){
            if(lineCounter == 3){
                std::vector<std::string> cellVec;
                std::stringstream str(dataLine);
                std::string cell;
                while(getline(str,cell,'=')){
                    cellVec.push_back(cell);
                }
                ARows = (long long)atoi(cellVec[1].c_str());
            }
            lineCounter++;
            continue;
        }
        std::vector<std::string> cellVec;
        std::stringstream str(dataLine);
        std::string cell;
        while(getline(str,cell,'\t')){
            cellVec.push_back(cell);
        }

        double elA;
        int eliA,eljA;

        eliA = atoi(cellVec[1].c_str());
        eljA = atoi(cellVec[2].c_str());
        elA = atof(cellVec[3].c_str());

        parconf.fIA.push_back(eliA);
        parconf.fJA.push_back(eljA);
        parconf.fA.push_back(elA);

        lineCounter++;
    }
    std::cout << "Stop reading matrix A\n";
    //Read rhs
    std::cout << "Start reading matrix b\n";
    std::ifstream rhsFile;
    rhsFile.open(pathtofilerhs);
    lineCounter=0;
    while(std::getline(rhsFile,dataLine)){
        if(lineCounter < 1){
            lineCounter++;
            continue;
        }
        std::stringstream str(dataLine);
        double elrhs;
        std::string cell;
        if(getline(str,cell))
            elrhs = atof(cell.c_str());
        parconf.frhs.push_back(elrhs);
        lineCounter++;
    }

    std::cout << "Stop reading matrix b\n";

    parconf.fIA.resize(ARows+1);

    if(0){
      int rhsnorm1=0;
    for (int ii=0; ii<parconf.frhs.size();ii++){
      rhsnorm1+=abs(parconf.frhs[ii]);
    }
    std::cout << "rhsnorm1=\t" << rhsnorm1 << "\tparconf.frhs.size()=\t" << parconf.frhs.size() << std::endl;

    int fAnorm1=0;
    for (int ii=0; ii<parconf.fA.size();ii++){
      fAnorm1+=abs(parconf.fA[ii]);
    }

      CheckConsistency<double>( parconf.frhs,  "rhs" );
      CheckConsistency<long long>( parconf.fIA,  "IA" );
      CheckConsistency<double>( parconf.fA,  "A" );
      CheckConsistency<long long>( parconf.fJA,  "JA" );

      std::cout << "parconf.fA.size():\t"<< parconf.fA.size() << std::endl;
      std::cout << "parconf.fIA.size():\t"<< parconf.fIA.size() << std::endl;
      std::cout << "parconf.fJA.size():\t"<< parconf.fJA.size() << std::endl;

      double JAnorm1=0;
      double IAnorm1=0;
      double IAmax=0;
      double JAmax=0;
      for(long long i=0; i< parconf.fJA.size(); i++){
          long long currentJA = abs(parconf.fJA[i]);
          JAnorm1+=currentJA;
          if(currentJA>JAmax)
            JAmax=currentJA;
          if(i<parconf.fIA.size()){
            long long currentIA = abs(parconf.fIA[i]);
            IAnorm1+=abs(parconf.fIA[i]);
            if(currentIA>IAmax)
              IAmax=currentIA;
            }
      }
      std::cout << "JAnorm1=\t" << JAnorm1 << std::endl;
      std::cout << "IAnorm1=\t" << IAnorm1 << std::endl;
      std::cout << "JAmax=\t" << JAmax << std::endl;
      std::cout << "IAmax=\t" << IAmax << std::endl;
      return 0;
    }

    std::cout << "linecounter=\t" << lineCounter <<std::endl;
    if(1){
        std::cout << "parconf.frhs.size():\t" << parconf.frhs.size() << "\n" << "ARows:\t" << ARows << std::endl;
        if(parconf.frhs.size() != ARows){
            std::cout << "Error" <<std::endl;
            return -1;
        }
    }

    if(parconf.fSystemType==ESymmetric ){
        if(parconf.fProperty == EIndefinite)
            parconf.fMatrixType =-2;
        else {
            parconf.fMatrixType =2;
        }
    }

    { //Decomposition
        long long n = 0;
        double bval = 0., xval = 0.;
        double *a,*b = &bval, *x = &xval;
        long long *ia,*ja;

        if (parconf.fA.size()==0) {
            return -1;
        }

        a = &(parconf.fA[0]);
        ia = (long long *) &(parconf.fIA[0]);
        ja = (long long *) &(parconf.fJA[0]);
        n = ARows;

        for(int i=0;i<n;i++){
            parconf.fPermutation.push_back(0);
        }
        long long *perm = 0,nrhs = 1;
        long long Error = 0;

        parconf.fPermutation.resize(n);
        perm = &(parconf.fPermutation[0]);


        /// analyse and factor the equations
        const long long int phase = 12;
        for (long long i=0; i<n; i++) {
           parconf.fPermutation[i] = i;
        }
        perm = &(parconf.fPermutation[0]);



        PrintArray<long long>(parconf.fParam,64,16, "fParam");

        std::cout << "pardiso_64" <<std::endl;
        pardiso_64(parconf.fHandle,  &parconf.fMax_num_factors, &parconf.fMatrix_num, &parconf.fMatrixType, &phase, &n, &parconf.fA[0],&parconf.fIA[0], &parconf.fJA[0], perm /*&idumgit*/, &nrhs, &parconf.fParam[0] /*iparm*/, &parconf.fMessageLevel,&parconf.frhs[0], x, &Error);

        if (Error) {
            parconf.Error_check(int(Error));
            return -1;
        }
        if (Error) {
            std::cout << __PRETTY_FUNCTION__ << " error code " << Error << std::endl;
            return -1;
        }
    }

    double *x;
    double sol[ARows];
    for(int ii=0; ii<ARows; ii++)
        sol[ii]=0;

    x = &sol[0];//&sol(0,0);
    mkl_set_num_threads(1);
    { //Solve
        long long n=0;
        double *a,*b;
        long long *ia,*ja;

        a = &(parconf.fA[0]);
        ia = (long long *) &(parconf.fIA[0]);
        ja = (long long *) &(parconf.fJA[0]);
        n = ARows;

        long long *perm,nrhs;
        long long Error = 0;
        nrhs = 1; //rhs.Cols();
        b = &(parconf.frhs)[0]; //&rhs(0,0);

        perm = &(parconf.fPermutation)[0]; //&fPermutation[0];
        /// forward and backward substitution
        long long phase = 33;

        pardiso_64(parconf.fHandle, &parconf.fMax_num_factors, &parconf.fMatrix_num, &parconf.fMatrixType, &phase, &n, &parconf.fA[0], &parconf.fIA[0], &parconf.fJA[0], perm /*&idumgit*/,
                    &nrhs, &parconf.fParam[0], &parconf.fMessageLevel, &parconf.frhs[0], x, &Error);

        if(parconf.fParam[19]>150){
            std::cout << "Pardiso:: Number of iterations " << parconf.fParam[19] << " > 150, calling numerical factorization... " << std::endl;
            phase = 23;
            pardiso_64 (parconf.fHandle,  &parconf.fMax_num_factors, &parconf.fMatrix_num, &parconf.fMatrixType, &phase, &n, a, ia, ja, perm,
                        &nrhs, &parconf.fParam[0], &parconf.fMessageLevel, b, x, &Error);
        }

        int rest = abs(parconf.fParam[19]%10); // CG/CGS error report
        if(parconf.fParam[19] <= 0){
            switch (rest) {
                case 1:{
                    std::cout << "Pardiso:: fluctuations of the residuum are too large. " << std::endl;
                }
                    break;

                case 2:{
                    std::cout << "Pardiso:: Slow convergence - Main matrix and matrix for preconditioner differ a lot. " << std::endl;
                }
                    break;

                case 4:{
                    std::cout << "Pardiso:: perturbed pivots caused iterative refinement. " << std::endl;
                }
                    break;

                case 5:{
                    std::cout << "Pardiso:: factorization is too fast for this matrix. It is better to use the factorization method with iparm[3] = 0 " << std::endl;
                    parconf.fParam[3] = 0;
                }
                    break;
                case 6:{
                    std::cout << "Pardiso:: There is not a diagnostig. " << std::endl;
                }
                    break;
                default:
                    break;
            }

        }


        if (Error<0) {
            parconf.Error_check(int(Error));
            std::cout << "Pardiso:: Calling a numerical factorization. \n";
            phase = 23;
            pardiso_64 (parconf.fHandle,  &parconf.fMax_num_factors, &parconf.fMatrix_num, &parconf.fMatrixType, &phase, &n, &parconf.fA[0], &parconf.fIA[0], &parconf.fJA[0], perm,
                        &nrhs, &parconf.fParam[0], &parconf.fMessageLevel, &parconf.frhs[0], x, &Error);
        }

        if (Error) {
            parconf.Error_check(int(Error));
            return -1; //DebugStop();
        }
    }

    std::vector<double> residue(ARows,0);
    std::vector<double> Ax(ARows,0);

    for(int i=0; i< ARows; i++){
        int nonZeroCols=parconf.fIA[i+1]-parconf.fIA[i];
        for(int icol=0; icol < nonZeroCols; icol++){
            int currentCol=parconf.fJA[parconf.fIA[i]+icol];
            int symCol = i;
            int symRow = currentCol;
            Ax[i] += parconf.fA[parconf.fIA[i]+icol]*x[currentCol];
            if(currentCol != i){
                Ax[symRow] += parconf.fA[parconf.fIA[i]+icol]*x[symCol];
            }
        }
    }

    double norm1=0;
    for(int i=0; i< ARows; i++){
        residue[i]=Ax[i]-parconf.frhs[i];
        norm1+=abs(residue[i]);
    }
    std::cout << "Norm1(Ax-b)=\t" << norm1 << std::endl;

    return 0;
}

void FilliParm(ParConfig &parconf, bool usingCGS){
    parconf.fParam[0]=1;
    parconf.fParam[1]=2;
    parconf.fParam[9]=8;
    parconf.fParam[17]=-1;
    parconf.fParam[18]=-1;
    parconf.fParam[20]=1;
    parconf.fParam[34] = 1;
    // Do not use OOC
    parconf.fParam[59] = 0;

    /// analyse and factor the equations
    // LU preconditioned CGS (10*L+K) where K={1:CGS,2:CG} and L=10^-L stopping threshold
    if (parconf.fProperty == EIndefinite) {
        parconf.fParam[4] = 1;
        if(parconf.fSystemType == ESymmetric){ // The factorization is always computed as required by phase.
            parconf.fParam[3 ] = 10*6+2;
        }else{ // CGS iteration replaces the computation of LU. The preconditioner is LU that was computed at a previous step (the first step or last step with a failure) in a sequence of solutions needed for identical sparsity patterns.
            parconf.fParam[3 ] = 10*6+1;
            parconf.fParam[10] = 1;
            parconf.fParam[12] = 1;
        }
    }else{

        if(parconf.fSystemType == ESymmetric){ // CGS iteration for symmetric positive definite matrices replaces the computation of LLT. The preconditioner is LLT that was computed at a previous step (the first step or last step with a failure) in a sequence of solutions needed for identical sparsity patterns.
            parconf.fParam[3 ] = 10*6+2;
        }else{
            parconf.fParam[3 ] = 10*6+1;
            parconf.fParam[10] = 1;
            parconf.fParam[12] = 1;
        }
    }

    if(!usingCGS){
        parconf.fParam[4] = 0;
        parconf.fParam[3] = 0;
        parconf.fParam[10] = 1;
        parconf.fParam[12] = 1;
    }
}

void ParConfig::Error_check(int error) const {

    switch (error) {
        case -1:
            std::cout << "Pardiso:: Input inconsistent." << std::endl;
            break;
        case -2:
            std::cout << "Pardiso:: Not enough memory." << std::endl;
            break;
        case -3:
            std::cout << "Pardiso:: Reordering problem." << std::endl;
            break;
        case -4:
            std::cout << "Pardiso:: Zero pivot, numerical fact. or iterative refinement problem. " << std::endl;
            break;
        case -5:
            std::cout << "Pardiso:: Unclassified (internal) error. " << std::endl;
            break;
        case -6:
            std::cout << "Pardiso:: Preordering failed (matrix types 11, 13 only). " << std::endl;
            break;
        case -7:
            std::cout << "Pardiso:: Diagonal matrix problem. " << std::endl;
            break;
        case -8:
            std::cout << "Pardiso:: 32-bit integer overflow problem. " << std::endl;
            break;
        default:
            std::cout << "Pardiso:: There is not a explanation. " << std::endl;
            break;
    }

}

template <class TVar>
void PrintArray(TVar *it,long long size, int n, std::string arrayName /*ArrayType*/ ){
    std::cout << arrayName << "[" << size << "]=\n\t\{";
    for(int ii=0;ii<size; ii++){
        if(ii%n==0 && ii!=0){
            std::cout<< "\n\t ";
        }
        std::cout << *(it+ii) <<",\t";
    }
    std::cout << "}\n";
    return;
}

template <class TVar>
void CheckConsistency(std::vector<TVar> vector, std::string name ){
    TVar *rhsPointer=&vector[0];
    for(int ii=0; ii<vector.size()-2; ii++){
       rhsPointer++;
    }
    std::cout << "\n "<< name <<" vector:\t" <<  vector[vector.size()-2] << "\t" << vector[vector.size()-1] << "\n";
    std::cout << name << " Pointer:\t" <<  *rhsPointer << "\t";
    rhsPointer++;
    std::cout << *rhsPointer << "\n";
}

//void PrintVector(std::vector<TVar> vec, long long size, std::string arrayName, ArrayType ){
//    std::cout << arrayName << "[" << size << "]=\n\t\{";
//    for(int ii=0;ii<size; ii++){
//        if(ii%5==0 && ii!=0){
//            std::cout<< "\n\t ";
//        }
//        std::cout << *(it+ii) <<",\t";
//    }
//    std::cout << "}\n";
//    return;
//}


// #include "pzfmatrix.h"
// #define EIGEN_USE_MKL_ALL
// #include <Eigen/Core>
//
// #include <Eigen/Sparse>
// #include <unsupported/Eigen/SparseExtra>
// #include <iostream>
// #include <Eigen/OrderingMethods>
// #include <Eigen/PardisoSupport>
// #include <Eigen/Sparse>
// #include <unsupported/Eigen/SparseExtra>
// #include <iostream>
// #include <Eigen/OrderingMethods>
// #include <Eigen/PardisoSupport>
//
// #include "pzvec.h"
// #include "pzmatrix.h"
//
// #include <pzskylmat.h>
// #include <pzstepsolver.h>
// #include <mkl.h>
// #include <stdio.h>
// #include <stdlib.h>
//
// #include <iostream>
// #include <chrono>
// #include "pzblockdiag.h"
// #include "tpzsparseblockdiagonal.h"
// using namespace std;
// using namespace Eigen;
// using namespace std;
// using std::chrono::high_resolution_clock;
// using std::chrono::duration_cast;
// using std::chrono::duration;
// using std::chrono::milliseconds;
//
// void FillMat ( TPZMatrix<REAL> &mat,REAL max, REAL min );
//
// void SolveSkyLine ( TPZSkylMatrix<REAL>  *skyline,  TPZFMatrix<REAL> rhs );
//
// void SolveSkyLine1 ( TPZSkylMatrix<REAL>  *skyline,  TPZFMatrix<REAL> rhs );
//
// void CopyFromSky ( TPZSkylMatrix<REAL> skyline,TPZFMatrix<REAL> & mat );
//
// void CopyFromFull ( TPZFMatrix<REAL> fullmat,TPZSkylMatrix<REAL> & skymat );
//
// void ReadMatrixEigen ( SparseMatrix<double> &A,VectorXd & b );
//
// void SolveEigen(MatrixXd mat,VectorXd rhs);
//
// void FromEigen ( MatrixXd eigenmat, TPZFMatrix<REAL>  &pzmat );
//
// void ToEigen ( TPZFMatrix<REAL>  pzmat,MatrixXd &eigenmat );
//
// void ToEigen ( TPZFMatrix<REAL>  pzmat,VectorXd &eigenmat );
//
// void CreateMatrix(TPZFMatrix<REAL> &fullmat,TPZSkylMatrix<REAL>  &skylinemat,  TPZFMatrix<REAL> &rhs,int sz);
//
// void ToSparseMatrixEigen (MatrixXd source, SparseMatrix<double> &A );
//
// int main()
// {
//     TPZFMatrix<REAL> fullmat,rhs;
//     TPZSkylMatrix<REAL>  skylinemat;
//     MatrixXd eigenmat;
//     VectorXd eigenrhs;
//     int sz=20000;
//     CreateMatrix(fullmat,skylinemat,rhs,sz);
//     //fullmat.Print ( "Matriz cheia ",cout );
//     //skylinemat.Print ( "Matriz skyline ",cout );
//     ToEigen(fullmat,eigenmat);
//     ToEigen(rhs,eigenrhs);
//     SolveEigen(eigenmat,eigenrhs);
//     SolveSkyLine1 ( &skylinemat,  rhs );
//    // std::ofstream sout("outmat.dat");
//    // fullmat.Print("",sout,EFormatted);
//
//     return 0;
// }
// void SolveEigen(MatrixXd mat,VectorXd rhs)
// {
//      SparseMatrix<double> A;
//     ToSparseMatrixEigen (mat,A );
//     SimplicialLLT< SparseMatrix<double> > solver;
//
//     std::cout << "start solving with Eigen SimplicialLLT"<< endl ;
//     auto t1 = high_resolution_clock::now();
//
//     VectorXd x = solver.compute ( A ).solve ( rhs );
//
//     auto t2 = high_resolution_clock::now();
//     auto ms_int = duration_cast<milliseconds> ( t2 - t1 );
//     std::cout << "tempo total  = "<<ms_int.count() << " ms\n";
//
//     //cout << x << endl;
// }
//
// void ToSparseMatrixEigen (MatrixXd source, SparseMatrix<double> &A )
// {
//     int m = source.rows();
//     A.resize(m,m);
//     for ( int i = 0; i < m; i++ )
//     {
//         for ( int j = 0; j < m; j++ )
//         {
//             if ( fabs ( source ( i,j ) ) !=0. )
//             {
//                 A.coeffRef ( i, j ) =source(i,j);
//             }
//         }
//
//     }
// }
//
// void FromEigen ( MatrixXd eigenmat, TPZFMatrix<REAL>  &pzmat )
// {
//
//     int rows = eigenmat.rows();
//     int cols = eigenmat.cols();
//     pzmat.Resize ( rows,cols );
//     for ( int irow=0; irow<rows; irow++ )
//     {
//         for ( int icol=0; icol<cols; icol++ )
//         {
//             pzmat ( irow,icol ) =eigenmat ( irow,icol );
//         }
//     }
//
// }
//
// void ToEigen ( TPZFMatrix<REAL>  pzmat,MatrixXd &eigenmat )
// {
//     TPZFMatrix<REAL> intpz ( pzmat );
//     int rows = pzmat.Rows();
//     int cols = pzmat.Cols();
//     eigenmat.resize ( rows,cols );
//     for ( int irow=0; irow<rows; irow++ )
//     {
//         for ( int icol=0; icol<cols; icol++ )
//         {
//             eigenmat ( irow,icol ) =intpz ( irow,icol );
//         }
//     }
//
// }
// void ToEigen ( TPZFMatrix<REAL>  pzmat,VectorXd &eigenmat )
// {
//     TPZFMatrix<REAL> intpz ( pzmat );
//     int rows = pzmat.Rows();
//     eigenmat.resize ( rows );
//     for ( int irow=0; irow<rows; irow++ )
//     {
//         eigenmat ( irow ) =intpz ( irow,0 );
//     }
//
// }
//
// void CreateMatrix(TPZFMatrix<REAL> &fullmat,TPZSkylMatrix<REAL>  &skylinemat,  TPZFMatrix<REAL> &rhs,int sz)
// {
//    // TPZBlockDiagonal<REAL> blockmat(sz,sz);
//
//     int64_t neq=sz;
//     skylinemat.Resize ( neq,neq );
//     skylinemat.AutoFill ( neq,neq,1 );
//    // blockmat.AutoFill(neq,neq,1);
//     CopyFromSky ( skylinemat,fullmat );
//
//     rhs.Resize(neq,1);
//     rhs.AutoFill();
//    // fullmat.Print ( "Matriz cheia ",cout );
//     //skylinemat.Print ( "Matriz skyline ",cout );
//
// }
//
// void CopyFromSky ( TPZSkylMatrix<REAL> skyline,TPZFMatrix<REAL> & mat )
// {
//     int rows=skyline.Rows();
//     int cols=skyline.Cols();
//     mat.Resize ( rows,cols );
//     for ( int irow=0; irow<rows; irow++ )
//     {
//         for ( int icol=0; icol<cols; icol++ )
//         {
//             mat.PutVal ( irow,icol,skyline.GetVal ( irow,icol ) );
//         }
//     }
//
//
// }
// void CopyFromFull ( TPZFMatrix<REAL> fullmat,TPZSkylMatrix<REAL> & skymat )
// {
//     int rows=fullmat.Rows();
//     int cols=fullmat.Cols();
//     skymat.Resize ( rows,cols );
//     for ( int irow=0; irow<skymat.Rows(); irow++ )
//     {
//         for ( int icol=0; icol<skymat.Cols(); icol++ )
//         {
//             skymat.PutVal ( irow,icol,fullmat.GetVal ( irow,icol ) );
//         }
//     }
//
//
// }
// void SolveSkyLine ( TPZSkylMatrix<REAL>  *skyline,  TPZFMatrix<REAL> rhs )
// {
//
//     TPZFMatrix<REAL> sol;
//     TPZStepSolver<REAL> step ( skyline );
//     TPZStepSolver<REAL> precond ( step );
//     int64_t numiterpre =2;
//     int64_t numiter = 5;
//     double overrelax = 1.1;
//     double tol = 1e-8;
//     precond.SetSSOR ( numiterpre,overrelax,tol,0 );
//     step.SetCG ( numiter,precond,tol,0 );
//
//     std::cout << "start solving with SetSSOR and SetCG"<< endl ;
//     auto t1 = high_resolution_clock::now();
//     step.Solve ( rhs,sol );
//     auto t2 = high_resolution_clock::now();
//     auto ms_int = duration_cast<milliseconds> ( t2 - t1 );
//
//     std::cout << "tempo total  = "<<ms_int.count() << " ms\n";
//
//    // sol.Print(std::cout);
//
// }
//
// void SolveSkyLine1 ( TPZSkylMatrix<REAL>  *skyline,  TPZFMatrix<REAL> rhs )
// {
//
//     TPZFMatrix<REAL> sol;
//     TPZStepSolver<REAL> step ( skyline );
//     step.SetDirect ( ECholesky );
//
//     std::cout << "start solving direct ECholesky"<< endl ;
//     auto t1 = high_resolution_clock::now();
//     step.Solve ( rhs,sol );
//     auto t2 = high_resolution_clock::now();
//     auto ms_int = duration_cast<milliseconds> ( t2 - t1 );
//
//     std::cout << "tempo total  = "<<ms_int.count() << " ms\n";
//
//     //sol.Print(std::cout);
//
// }
//
