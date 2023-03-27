/**
 * @file
 */


#include "TPZPlasticStepPV.h"
#include "TPZElasticResponse.h"

#include "pzsandlerextPV.h"
#include "TPZYCMohrCoulombPV.h"

#include "TPZElasticResponse.h"

#include "pzlog.h"

//#ifdef LOG4CXX
//static LoggerPtr logger(Logger::getLogger("pz.material.TPZPlasticStepPV"));
//#endif
#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("plasticity.poroelastoplastic2"));
#endif

#ifdef LOG4CXX
static LoggerPtr logger2(Logger::getLogger("plasticity.poroelastoplastic"));
#endif

template <class YC_t, class ER_t>
void TPZPlasticStepPV<YC_t, ER_t>::ApplyStrainComputeSigma(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma)
{
	TPZTensor<REAL>::TPZDecomposed DecompSig; // It may be SigTr or SigPr Decomposition, dependes on the part of this method
	TPZTensor<REAL> sigtr;
	if(fN.fmatprop.size()!=0 && fN.fmatprop[0]>1.e-3)
 	{
	 fYC.SetLocalMatState(fN);
	}
	//
	TPZTensor<REAL> epsTr,epsPN,epsElaNp1;
	epsPN = fN.fEpsP;
	epsTr = epsTotal;
	epsTr -= epsPN; // Porque soh tem implementado o operator -=
	
	// Compute and Decomposition of SigTrial
	fER.Compute(epsTr, sigtr); // sigma = lambda Tr(E)I + 2 mu E
	sigtr.EigenSystem(DecompSig);
	TPZManVector<REAL,3> sigtrvec(DecompSig.fEigenvalues), sigprvec(3,0.);
	
	// ReturMap in the principal values
	STATE nextalpha = -6378.;
	fYC.ProjectSigma(sigtrvec, fN.fAlpha, sigprvec, nextalpha);
	fN.fAlpha = nextalpha;
#ifdef LOG4CXX_KEEP
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout << "Sig Trial " << sigtrvec << "\nSig Project " << sigprvec << std::endl;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif

	// Reconstruction of sigmaprTensor
	DecompSig.fEigenvalues = sigprvec; // CHANGING THE EIGENVALUES FOR THE ONES OF SIGMAPR
	sigma = TPZTensor<REAL>(DecompSig);
	
	fER.ComputeDeformation(sigma,epsElaNp1);
	fN.fEpsT = epsTotal;
	epsPN = epsTotal;
	epsPN -= epsElaNp1; // Transforma epsPN em epsPNp1
	fN.fEpsP = epsPN; 
}


// template <class YC_t, class ER_t>
// void TPZPlasticStepPV<YC_t, ER_t>::ApplyStrainComputeDep(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma, TPZFMatrix<REAL> &Dep)
// {
//     TPZTensor<REAL>::TPZDecomposed DecompSig,DecompEps; // It may be SigTr or SigPr Decomposition, dependes on the part of this method
//     TPZTensor<REAL> sigtr;
//
// 	if(fN.fmatprop.size()!=0 && fN.fmatprop[0]>1.e-3)
//  	{
// 	 fYC.SetLocalMatState(fN);
// 	}
//
// 	//fYC.fPhi = cohesion;
//     //
//     TPZTensor<REAL> epsTr,epsPN,epsElaNp1;
//     epsPN = fN.fEpsP;
//     epsTr = epsTotal;
//     epsTr -= epsPN; // Porque soh tem implementado o operator -=
//
//     // Compute and Decomposition of SigTrial
//     fER.Compute(epsTr, sigtr); // sigma = lambda Tr(E)I + 2 mu E
//     epsTr.EigenSystem(DecompEps);
//
//     //
//     sigtr.EigenSystem(DecompSig);
//     TPZManVector<REAL,3> sigtrvec(DecompSig.fEigenvalues), sigprvec(3,0.);
//
//     // Pegando os autovetores
//     TPZManVector<TPZFNMatrix<3>,3> epsegveFromProj(3);
//     TPZManVector<TPZFNMatrix<9,REAL>, 3 > EigenvecMat(3);
// 	for (int i = 0; i < 3; i++)
// 	{
// 		EigenvecMat[i] = DecompSig.fEigenvectors[i];
// 		epsegveFromProj[i].Resize(3,1);
//         STATE maxvecnorm = 0;
//         int maxvecindex = 0;
// 		for	(int k = 0 ; k < 3 ; k++){
//             STATE vecnorm=0.;
//             for (int j=0; j<3; j++) {
//                 vecnorm += EigenvecMat[i](j,k)*EigenvecMat[i](j,k);
//             }
//             if (vecnorm> maxvecnorm) {
//                 maxvecindex = k;
//                 maxvecnorm = vecnorm;
//             }
// 		}
//         for (int j=0; j<3; j++) {
//             epsegveFromProj[i](j,0) = EigenvecMat[i](j,maxvecindex);
//         }
// 	}
//     for (int i = 0; i < 3; i++) {
//         REAL normvec = 0.;
//         normvec = NormVecOfMat(epsegveFromProj[i]);
//         for (int j = 0; j < 3; j++) {
//             epsegveFromProj[i](j,0) /= normvec;
//         }
//     }
// #ifdef LOG4CXX
//     if(logger->isDebugEnabled())
//     {
//         std::stringstream sout;
//         sout << "\nEigenvectors:\n " << std::endl;
//         for (int i = 0 ; i < 3; i++) {
//             sout << "i = " << i << std::endl;
//             epsegveFromProj[i].Print("epsegveFromProj",sout);
//             EigenvecMat[i].Print("EigenvecMat",sout);
//             sout << std::endl;
//         }
//         LOGPZ_DEBUG(logger, sout.str())
//     }
// #endif
//
//
//
//
//     // ReturMap in the principal values
//     STATE nextalpha = -6378.;
//     STATE printPlastic = fN.Alpha();
//     TPZFNMatrix<9> GradSigma(3,3,0.);
// 	TPZTensor<REAL> sigtrtensor;
//     fYC.ProjectSigmaDep(sigtrvec,sigtrtensor, fN.fAlpha, sigprvec, nextalpha, GradSigma);
//     //GradSigma.Print("Grad");
//     fN.fAlpha = nextalpha;
//
//
//
// #ifdef LOG4CXX
//     if(logger->isDebugEnabled())
//     {
//         std::stringstream sout;
//         sout << "Sig Trial " << sigtrvec << "\nSig Project " << sigprvec << std::endl;
//         GradSigma.Print("GradSigma", sout,EMathematicaInput);
//         LOGPZ_DEBUG(logger, sout.str())
//     }
// #endif
//
//     // Aqui calculo minha matriz tangente ------------------------------------
//     // Criando matriz tangente
//     TPZFNMatrix<36> dSigDe(6,6,0.);
//
//     //Montando a matriz tangente
//     int kival[] = {0,0,0,1,1,2};
//     int kjval[] = {0,1,2,1,2,2};
//     REAL G = fER.G();
//     REAL lambda = fER.Lambda();
//     int ki, kj;
//     for (int k = 0; k < 6; k++)
//     {
//
//         ki = kival[k];
//         kj = kjval[k];
//         for (int i = 0; i < 3; i++)
//         {
//             for (int j = 0; j < 3; j++)
//             {
//                 for (int l = 0; l < 6; l++)
//                 {
//                     REAL temp = 2 * G * EigenvecMat[j](kj,ki); // * EigenvecMat[j](j,ki);
//
//                     if (ki == kj)
//                     {
//                         temp += lambda;
//                     }
//                     else {
//                         temp *= 2.;
//                     }
//
//                     temp *= GradSigma(i,j);
//                     dSigDe(l,k) += temp * DecompSig.fEigenvectors[i][l];
//                 }///l
//             }///j
//         }///i
//     }///k
//
//     REAL deigensig = 0., deigeneps = 0.;
//     TPZFNMatrix<36> RotCorrection(6,6,0.);
//     // Correcao do giro rigido
//     for (int i = 0; i < 2; i++) {
//         for (int j = i+1; j<3 ; j++) {
//             deigeneps = DecompEps.fEigenvalues[i]  - DecompEps.fEigenvalues[j];
//             deigensig = sigprvec[i] - sigprvec[j];
//             TPZFNMatrix<9,REAL> tempMat(3,3,0.);
//             REAL factor = 0.;
//             if (!IsZero(deigeneps)) {
//                 factor = deigensig / deigeneps;
//             }
//             else {
//                 factor = fER.G() * ( GradSigma(i,i) - GradSigma(i,j) - GradSigma(j,i) + GradSigma(j,j) );
//             }
//             tempMat = ProdT(epsegveFromProj[i],epsegveFromProj[j]) + ProdT(epsegveFromProj[j],epsegveFromProj[i]);
//             for (int k = 0 ; k < 6 ; k++){
//                 ki = kival[k];
//                 kj = kjval[k];
//                 TPZFNMatrix<9> ColCorr(3,3,0.);
//                 TPZFNMatrix<6> ColCorrV(6,1,0.);
//                 if (ki == kj) {
//                     REAL tempval = (epsegveFromProj[j](ki,0) * epsegveFromProj[i](kj,0) );
//                     ColCorr = tempval * factor * tempMat;
//                 }
//                 else {
//                     ColCorr = (epsegveFromProj[j](ki,0) * epsegveFromProj[i](kj,0) + epsegveFromProj[j](kj,0) * epsegveFromProj[i](ki,0) ) * factor * tempMat;
//                 }
//                 ColCorrV = FromMatToVoight(ColCorr);
//                 for (int l = 0; l < 6; l++) {
//                     RotCorrection(l,k) += ColCorrV(l,0);
//                 }
//             }
//         }
//     }
//
//     dSigDe += RotCorrection;
//
// #ifdef LOG4CXX
//     {
//         if(logger->isDebugEnabled())
//         {
//             std::stringstream str;
//             str << "\n**********************MATRIZ TANGENTE**********************" << endl;
//             dSigDe.Print("Matriz Tangente:",str);
//             str << "\n**********************CORRECAO GIRO**********************" << endl;
//             RotCorrection.Print("GiroCorrection",str);
//             LOGPZ_DEBUG(logger,str.str())
//         }
//     }
// #endif
//
//     // Reconstruction of sigmaprTensor
//
//     DecompSig.fEigenvalues = sigprvec; // CHANGING THE EIGENVALUES FOR THE ONES OF SIGMAPR
//     sigma = TPZTensor<REAL>(DecompSig);
//
//
//
//
//
//     fER.ComputeDeformation(sigma,epsElaNp1);
//     fN.fEpsT = epsTotal;
//     epsPN = epsTotal;
//     epsPN -= epsElaNp1; // Transforma epsPN em epsPNp1
//     fN.fEpsP = epsPN;
//     Dep = dSigDe;
//
//
// #ifdef LOG4CXX
//     if(logger2->isDebugEnabled())
//     {
//         if(fabs(printPlastic-fN.fAlpha)>1.e-4)
//         {
//             std::stringstream sout;
//             TPZVec<STATE> phi;
//             TPZTensor<STATE> epsElastic(fN.fEpsT);
//             epsElastic-=fN.fEpsP;
//             Phi(epsElastic, phi);
//             sout << " \n phi = [";
//             for (int i=0;i<phi.size();i++)
//             {
//                 sout << phi[i] <<" ";
//             }
//
//             sout << " ] "<<endl;
//
//             sout << " \n eigenvalues Sigma = [";
//             for (int i=0;i<3;i++)
//             {
//                 sout << DecompSig.fEigenvalues[i] <<" ";
//             }
//
//             sout << " ] "<<endl;
//
//
//
//             LOGPZ_DEBUG(logger2, sout.str())
//         }
//     }
// #endif
// }


// template <class YC_t, class ER_t>
// void TPZPlasticStepPV<YC_t, ER_t>::ApplyStrainComputeDep(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma, TPZFMatrix<REAL> &Dep)
// {
//     TPZTensor<REAL>::TPZDecomposed DecompSig,DecompEps; // It may be SigTr or SigPr Decomposition, dependes on the part of this method
//     TPZTensor<REAL> sigtr;
//
// 	if(fN.fmatprop.size()!=0 && fN.fmatprop[0]>1.e-3)
//  	{
// 	 fYC.SetLocalMatState(fN);
// 	}
//
// 	//fYC.fPhi = cohesion;
//     //
//     TPZTensor<REAL> epsTr,epsPN,epsElaNp1;
//     epsPN = fN.fEpsP;
//     epsTr = epsTotal;
//     epsTr -= epsPN; // Porque soh tem implementado o operator -=
//
//     // Compute and Decomposition of SigTrial
//     fER.Compute(epsTr, sigtr); // sigma = lambda Tr(E)I + 2 mu E
//     epsTr.EigenSystem(DecompEps);
//
//     //
//     sigtr.EigenSystem(DecompSig);
//     TPZManVector<REAL,3> sigtrvec(DecompSig.fEigenvalues), sigprvec(3,0.);
//
//     // Pegando os autovetores
//     TPZManVector<TPZFNMatrix<3>,3> epsegveFromProj(3);
//     TPZManVector<TPZFNMatrix<9,REAL>, 3 > EigenvecMat(3);
// 	for (int i = 0; i < 3; i++)
// 	{
// 		EigenvecMat[i] = DecompSig.fEigenvectors[i];
// 		epsegveFromProj[i].Resize(3,1);
//         STATE maxvecnorm = 0;
//         int maxvecindex = 0;
// 		for	(int k = 0 ; k < 3 ; k++){
//             STATE vecnorm=0.;
//             for (int j=0; j<3; j++) {
//                 vecnorm += EigenvecMat[i](j,k)*EigenvecMat[i](j,k);
//             }
//             if (vecnorm> maxvecnorm) {
//                 maxvecindex = k;
//                 maxvecnorm = vecnorm;
//             }
// 		}
//         for (int j=0; j<3; j++) {
//             epsegveFromProj[i](j,0) = EigenvecMat[i](j,maxvecindex);
//         }
// 	}
//     for (int i = 0; i < 3; i++) {
//         REAL normvec = 0.;
//         normvec = NormVecOfMat(epsegveFromProj[i]);
//         for (int j = 0; j < 3; j++) {
//             epsegveFromProj[i](j,0) /= normvec;
//         }
//     }
// #ifdef LOG4CXX
//     if(logger->isDebugEnabled())
//     {
//         std::stringstream sout;
//         sout << "\nEigenvectors:\n " << std::endl;
//         for (int i = 0 ; i < 3; i++) {
//             sout << "i = " << i << std::endl;
//             epsegveFromProj[i].Print("epsegveFromProj",sout);
//             EigenvecMat[i].Print("EigenvecMat",sout);
//             sout << std::endl;
//         }
//         LOGPZ_DEBUG(logger, sout.str())
//     }
// #endif
//
//
//
//
//     // ReturMap in the principal values
//     STATE nextalpha = -6378.;
//     STATE printPlastic = fN.Alpha();
//     TPZFNMatrix<9> GradSigma ;
//     fYC.ProjectSigmaDep(sigtrvec, sigtr,fN.fAlpha, sigprvec, nextalpha, GradSigma);
//     //GradSigma.Print("Grad");
// //     cout << "sig proj " << endl;
// //     for(int i=0;i<3;i++)std::cout<< sigprvec[i] << std::endl;
//     fN.fAlpha = nextalpha;
//     Dep.Resize(6,6);
//
//     Dep= GradSigma;
//         // Reconstruction of sigmaprTensor
//
//     DecompSig.fEigenvalues = sigprvec; // CHANGING THE EIGENVALUES FOR THE ONES OF SIGMAPR
//     sigma = TPZTensor<REAL>(DecompSig);
//
//     fER.ComputeDeformation(sigma,epsElaNp1);
//     fN.fEpsT = epsTotal;
//     epsPN = epsTotal;
//     epsPN -= epsElaNp1; // Transforma epsPN em epsPNp1
//     fN.fEpsP = epsPN;
// //
// }

template <class YC_t, class ER_t>
void TPZPlasticStepPV<YC_t, ER_t>::ApplyStrainComputeDep(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma, TPZFMatrix<REAL> &Dep)
{
    TPZTensor<REAL>::TPZDecomposed DecompSig,DecompEps; // It may be SigTr or SigPr Decomposition, dependes on the part of this method
    TPZTensor<REAL> sigtr;

	if(fN.fmatprop.size()!=0 && fN.fmatprop[0]>1.e-3)
 	{
	 fYC.SetLocalMatState(fN);
	}

	//fYC.fPhi = cohesion;
    //
    TPZTensor<REAL> epsTr,epsPN,epsElaNp1;
    epsPN = fN.fEpsP;
    epsTr = epsTotal;
    epsTr -= epsPN; // Porque soh tem implementado o operator -=

    // Compute and Decomposition of SigTrial
    fER.Compute(epsTr, sigtr); // sigma = lambda Tr(E)I + 2 mu E
    epsTr.EigenSystem(DecompEps);

    //
    sigtr.EigenSystem(DecompSig);
    TPZManVector<REAL,3> sigtrvec(DecompSig.fEigenvalues), sigprvec(3,0.);

    // Pegando os autovetores
    TPZManVector<TPZFNMatrix<3>,3> epsegveFromProj(3);
    TPZManVector<TPZFNMatrix<9,REAL>, 3 > EigenvecMat(3);
	for (int i = 0; i < 3; i++)
	{
		EigenvecMat[i] = DecompSig.fEigenvectors[i];
		epsegveFromProj[i].Resize(3,1);
        STATE maxvecnorm = 0;
        int maxvecindex = 0;
		for	(int k = 0 ; k < 3 ; k++){
            STATE vecnorm=0.;
            for (int j=0; j<3; j++) {
                vecnorm += EigenvecMat[i](j,k)*EigenvecMat[i](j,k);
            }
            if (vecnorm> maxvecnorm) {
                maxvecindex = k;
                maxvecnorm = vecnorm;
            }
		}
        for (int j=0; j<3; j++) {
            epsegveFromProj[i](j,0) = EigenvecMat[i](j,maxvecindex);
        }
	}
    for (int i = 0; i < 3; i++) {
        REAL normvec = 0.;
        normvec = NormVecOfMat(epsegveFromProj[i]);
        for (int j = 0; j < 3; j++) {
            epsegveFromProj[i](j,0) /= normvec;
        }
    }
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout << "\nEigenvectors:\n " << std::endl;
        for (int i = 0 ; i < 3; i++) {
            sout << "i = " << i << std::endl;
            epsegveFromProj[i].Print("epsegveFromProj",sout);
            EigenvecMat[i].Print("EigenvecMat",sout);
            sout << std::endl;
        }
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif




    // ReturMap in the principal values
    STATE nextalpha = -6378.;
    STATE printPlastic = fN.Alpha();
    TPZFNMatrix<9> GradSigma(3,3,0.);
	TPZTensor<REAL> sigtrtensor;
    fYC.ProjectSigmaDep(sigtrvec,sigtrtensor, fN.fAlpha, sigprvec, nextalpha, GradSigma);
    //GradSigma.Print("Grad");
    fN.fAlpha = nextalpha;



#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout << "Sig Trial " << sigtrvec << "\nSig Project " << sigprvec << std::endl;
        GradSigma.Print("GradSigma", sout,EMathematicaInput);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif

    // Aqui calculo minha matriz tangente ------------------------------------
    // Criando matriz tangente
    TPZFNMatrix<36> dSigDe(6,6,0.);

    //Montando a matriz tangente
    int kival[] = {0,0,0,1,1,2};
    int kjval[] = {0,1,2,1,2,2};
    REAL G = fER.G();
    REAL lambda = fER.Lambda();
    int ki, kj;
    for (int k = 0; k < 6; k++)
    {

        ki = kival[k];
        kj = kjval[k];
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                for (int l = 0; l < 6; l++)
                {
                    REAL temp = 2 * G * EigenvecMat[j](kj,ki); // * EigenvecMat[j](j,ki);

                    if (ki == kj)
                    {
                        temp += lambda;
                    }
                    else {
                        temp *= 2.;
                    }

                    temp *= GradSigma(i,j);
                    dSigDe(l,k) += temp * DecompSig.fEigenvectors[i][l];
                }///l
            }///j
        }///i
    }///k

    REAL deigensig = 0., deigeneps = 0.;
    TPZFNMatrix<36> RotCorrection(6,6,0.);
    // Correcao do giro rigido
    for (int i = 0; i < 2; i++) {
        for (int j = i+1; j<3 ; j++) {
            deigeneps = DecompEps.fEigenvalues[i]  - DecompEps.fEigenvalues[j];
            deigensig = sigprvec[i] - sigprvec[j];
            TPZFNMatrix<9,REAL> tempMat(3,3,0.);
            REAL factor = 0.;
            if (!IsZero(deigeneps)) {
                factor = deigensig / deigeneps;
            }
            else {
                factor = fER.G() * ( GradSigma(i,i) - GradSigma(i,j) - GradSigma(j,i) + GradSigma(j,j) );
            }
            tempMat = ProdT(epsegveFromProj[i],epsegveFromProj[j]) + ProdT(epsegveFromProj[j],epsegveFromProj[i]);
            for (int k = 0 ; k < 6 ; k++){
                ki = kival[k];
                kj = kjval[k];
                TPZFNMatrix<9> ColCorr(3,3,0.);
                TPZFNMatrix<6> ColCorrV(6,1,0.);
                if (ki == kj) {
                    REAL tempval = (epsegveFromProj[j](ki,0) * epsegveFromProj[i](kj,0) );
                    ColCorr = tempval * factor * tempMat;
                }
                else {
                    ColCorr = (epsegveFromProj[j](ki,0) * epsegveFromProj[i](kj,0) + epsegveFromProj[j](kj,0) * epsegveFromProj[i](ki,0) ) * factor * tempMat;
                }
                ColCorrV = FromMatToVoight(ColCorr);
                for (int l = 0; l < 6; l++) {
                    RotCorrection(l,k) += ColCorrV(l,0);
                }
            }
        }
    }

    dSigDe += RotCorrection;

#ifdef LOG4CXX
    {
        if(logger->isDebugEnabled())
        {
            std::stringstream str;
            str << "\n**********************MATRIZ TANGENTE**********************" << endl;
            dSigDe.Print("Matriz Tangente:",str);
            str << "\n**********************CORRECAO GIRO**********************" << endl;
            RotCorrection.Print("GiroCorrection",str);
            LOGPZ_DEBUG(logger,str.str())
        }
    }
#endif

    // Reconstruction of sigmaprTensor

    DecompSig.fEigenvalues = sigprvec; // CHANGING THE EIGENVALUES FOR THE ONES OF SIGMAPR
    sigma = TPZTensor<REAL>(DecompSig);





    fER.ComputeDeformation(sigma,epsElaNp1);
    fN.fEpsT = epsTotal;
    epsPN = epsTotal;
    epsPN -= epsElaNp1; // Transforma epsPN em epsPNp1
    fN.fEpsP = epsPN;
    Dep = dSigDe;


#ifdef LOG4CXX
    if(logger2->isDebugEnabled())
    {
        if(fabs(printPlastic-fN.fAlpha)>1.e-4)
        {
            std::stringstream sout;
            TPZVec<STATE> phi;
            TPZTensor<STATE> epsElastic(fN.fEpsT);
            epsElastic-=fN.fEpsP;
            Phi(epsElastic, phi);
            sout << " \n phi = [";
            for (int i=0;i<phi.size();i++)
            {
                sout << phi[i] <<" ";
            }

            sout << " ] "<<endl;

            sout << " \n eigenvalues Sigma = [";
            for (int i=0;i<3;i++)
            {
                sout << DecompSig.fEigenvalues[i] <<" ";
            }

            sout << " ] "<<endl;



            LOGPZ_DEBUG(logger2, sout.str())
        }
    }
#endif
}
// template <class YC_t, class ER_t>
// void TPZPlasticStepPV<YC_t, ER_t>::ApplyStrainComputeDep(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma, TPZFMatrix<REAL> &Dep) {
//
//     TPZTensor<REAL>::TPZDecomposed sig_eigen_system, eps_eigen_system;
//     TPZTensor<REAL> sigtr;
//
//     TPZTensor<REAL> epsTr, epsPN, epsElaNp1;
//    //epsPN = fN.m_eps_p;
// 	epsPN = fN.fEpsP;
//     epsTr = epsTotal;
//     epsTr -= epsPN; // Porque soh tem implementado o operator -=
//
//     // Compute and Decomposition of SigTrial
//     fER.Compute(epsTr, sigtr); // sigma = lambda Tr(E)I + 2 mu E
//     epsTr.EigenSystem(eps_eigen_system);
//     //epsTr.ComputeEigenvectors(eps_eigen_system);
//     sigtr.EigenSystem(sig_eigen_system);
//     //sigtr.ComputeEigenvectors(sig_eigen_system);
//
//     TPZManVector<REAL, 3> sigtrvec(sig_eigen_system.fEigenvalues), sigprvec(3, 0.);
//
// #ifdef PZ_LOG
//     if (logger.isDebugEnabled()) {
//         std::stringstream sout;
// 		sig_eigen_system.Print(sout);
//         LOGPZ_DEBUG(logger, sout.str())
//     }
//     STATE printPlastic = fN.VolHardening();
// #endif
//
//     // ReturnMap in the principal values
//     STATE nextalpha = 0;
//     TPZFNMatrix<9> GradSigma(3, 3, 0.);
//   //  fYC.ProjectSigmaDep(sigtrvec, fN.fAlpha, sigprvec, nextalpha, GradSigma);
// 		TPZTensor<REAL> sigtrtensor;
//     fYC.ProjectSigmaDep(sigtrvec,sigtrtensor, fN.fAlpha, sigprvec, nextalpha, GradSigma);
//     fN.fAlpha = nextalpha;
//  	//cout << "Sig Trial " << sigtrvec << "\nSig Project " << sigprvec << std::endl;
// #ifdef PZ_LOG
//     if (logger.isDebugEnabled()) {
//         std::stringstream sout;
//         sout << "Sig Trial " << sigtrvec << "\nSig Project " << sigprvec << std::endl;
//         GradSigma.Print("GradSigma", sout, EMathematicaInput);
//         LOGPZ_DEBUG(logger, sout.str())
//     }
// #endif
//
//     // Reconstruction of sigmaprTensor
//     sig_eigen_system.fEigenvalues = sigprvec; // updating the projected values used inside TangentOperator method.
//     sigma = TPZTensor<REAL>(sig_eigen_system);
//     TangentOperator(GradSigma, eps_eigen_system, sig_eigen_system, Dep);
//
//
//     fER.ComputeDeformation(sigma, epsElaNp1);
//     fN.fEpsT = epsTotal;
//     epsPN = epsTotal;
//     epsPN -= epsElaNp1; // Transforma epsPN em epsPNp1
//     fN.fEpsP = epsPN;
//
//
// #ifdef PZ_LOG
//     if (logger2.isDebugEnabled()) {
//         if (fabs(printPlastic - fN.m_hardening) > 1.e-4) {
//             std::stringstream sout;
//             TPZVec<STATE> phi;
//             TPZTensor<STATE> epsElastic(fN.m_eps_t);
//             epsElastic -= fN.m_eps_p;
//             Phi(epsElastic, phi);
//             sout << " \n phi = [";
//             for (int i = 0; i < phi.size(); i++) {
//                 sout << phi[i] << " ";
//             }
//
//             sout << " ] " << std::endl;
//
//             sout << " \n eigenvalues Sigma = [";
//             for (int i = 0; i < 3; i++) {
//                 sout << sig_eigen_system.fEigenvalues[i] << " ";
//             }
//
//             sout << " ] " << std::endl;
//
//
//
//             LOGPZ_DEBUG(logger2, sout.str())
//         }
//     }
// #endif
// }
//template <class YC_t, class ER_t>
//void TPZPlasticStepPV<YC_t, ER_t>::
template <class YC_t, class ER_t>
void TPZPlasticStepPV<YC_t, ER_t>::TangentOperator(TPZFMatrix<REAL> & gradient,TPZTensor<REAL>::TPZDecomposed & eps_eigen_system, TPZTensor<REAL>::TPZDecomposed & sig_eigen_system, TPZFMatrix<REAL> & Tangent){
	
	Tangent.Resize(6,6);
	  //  Pegando os autovetores
    TPZManVector<TPZFNMatrix<3>,3> epsegveFromProj(3);
    TPZManVector<TPZFNMatrix<9,REAL>, 3 > EigenvecMat(3);
	for (int i = 0; i < 3; i++)
	{
		EigenvecMat[i] = eps_eigen_system.fEigenvectors[i];
		epsegveFromProj[i].Resize(3,1);
        STATE maxvecnorm = 0;
        int maxvecindex = 0;
		for	(int k = 0 ; k < 3 ; k++){
            STATE vecnorm=0.;
            for (int j=0; j<3; j++) {
                vecnorm += EigenvecMat[i](j,k)*EigenvecMat[i](j,k);
            }
            if (vecnorm> maxvecnorm) {
                maxvecindex = k;
                maxvecnorm = vecnorm;
            }
		}
        for (int j=0; j<3; j++) {
            epsegveFromProj[i](j,0) = EigenvecMat[i](j,maxvecindex);
        }
	}
    for (int i = 0; i < 3; i++) {
        REAL normvec = 0.;
        normvec = NormVecOfMat(epsegveFromProj[i]);
        for (int j = 0; j < 3; j++) {
            epsegveFromProj[i](j,0) /= normvec;
        }
    }
    
    
//         // Aqui calculo minha matriz tangente ------------------------------------
    // Criando matriz tangente
    TPZFNMatrix<36> dSigDe(6,6,0.);
    
    //Montando a matriz tangente
    int kival[] = {0,0,0,1,1,2};
    int kjval[] = {0,1,2,1,2,2};
    REAL G = fER.G();
    REAL lambda = fER.Lambda();
    int ki, kj;
    for (int k = 0; k < 6; k++)
    {
        
        ki = kival[k];
        kj = kjval[k];
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                for (int l = 0; l < 6; l++)
                {
                    REAL temp = 2 * G * EigenvecMat[j](kj,ki); // * EigenvecMat[j](j,ki);
                   // REAL temp = 2 * G * EigenvecMat[j](kj,ki)* EigenvecMat[j](j,ki);
                    
                    if (ki == kj)
                    {
                        temp += lambda;
                    }
                    else {
                        temp *= 2.;
                    }
                    
                    temp *= gradient(i,j);
                    Tangent(l,k) += temp * eps_eigen_system.fEigenvectors[i][l];
                }///l
            }///j
        }///i
    }///k
    //Tangent.Print("tanzaquiii");
    REAL deigensig = 0., deigeneps = 0.;
    TPZFNMatrix<36> RotCorrection(6,6,0.);
    // Correcao do giro rigido
    for (int i = 0; i < 2; i++) {
        for (int j = i+1; j<3 ; j++) {
            deigeneps = eps_eigen_system.fEigenvalues[i]  - eps_eigen_system.fEigenvalues[j];
            deigensig = sig_eigen_system.fEigenvalues[i]- sig_eigen_system.fEigenvalues[j];
            TPZFNMatrix<9,REAL> tempMat(3,3,0.);
            REAL factor = 0.;
            if (!IsZero(deigeneps)) {
                factor = deigensig / deigeneps;
            }
            else {
                factor = fER.G() * ( gradient(i,i) - gradient(i,j) - gradient(j,i) + gradient(j,j) );
            }
            //cout << "epsegveFromProj[0] "<<endl;
			//cout << epsegveFromProj[0] <<endl;
            tempMat = ProdT(epsegveFromProj[i],epsegveFromProj[j]) + ProdT(epsegveFromProj[j],epsegveFromProj[i]);
            for (int k = 0 ; k < 6 ; k++){
                ki = kival[k];
                kj = kjval[k];
                TPZFNMatrix<9> ColCorr(3,3,0.);
                TPZFNMatrix<6> ColCorrV(6,1,0.);
                if (ki == kj) {
                    REAL tempval = (epsegveFromProj[j](ki,0) * epsegveFromProj[i](kj,0) );
                    ColCorr = tempval * factor * tempMat;
                }
                else {
                    ColCorr = (epsegveFromProj[j](ki,0) * epsegveFromProj[i](kj,0) + epsegveFromProj[j](kj,0) * epsegveFromProj[i](ki,0) ) * factor * tempMat;
                }
                ColCorrV = FromMatToVoight(ColCorr);
				//ColCorrV.Print("colCorrV");
                for (int l = 0; l < 6; l++) {
                    RotCorrection(l,k) += ColCorrV(l,0);
                }
            }
        }
    }
    //RotCorrection.Print("root");
    Tangent += RotCorrection;
    
    
    
    
  /*  
	
    Tangent.Resize(6,6);
	//TPZFNMatrix<36> dSigDe(6,6,0.);
	cout << "gradient " << gradient ;
	cout << "\neps_eigen_system";
	eps_eigen_system.Print(cout);
    //Montando a matriz tangente
    //unsigned int kival[] = {0, 0, 0, 1, 1, 2};
    //unsigned int kjval[] = {0, 1, 2, 1, 2, 2};
	int kival[] = {0,0,0,1,1,2};
    int kjval[] = {0,1,2,1,2,2};
    REAL G = fER.G();
    REAL lambda = fER.Lambda();

	

    // Coluna da matriz tangente
    for (unsigned int k = 0; k < 6; ++k) {
        const unsigned int ki = kival[k];
        const unsigned int kj = kjval[k];
        for (unsigned int i = 0; i < 3; ++i) {
            for (unsigned int j = 0; j < 3; ++j) {
				//REAL temp = 2 * G * EigenvecMat[j](kj,ki); // * EigenvecMat[j](j,ki);
				TPZManVector<REAL,6> v1,v2;
			v1[0]=eps_eigen_system.fEigenvectors[i].XX();
			v1[1]=eps_eigen_system.fEigenvectors[i].XY();
			v1[2]=eps_eigen_system.fEigenvectors[i].XZ();
			v1[3]=eps_eigen_system.fEigenvectors[i].YY();
			v1[4]=eps_eigen_system.fEigenvectors[i].YZ();
			v1[5]=eps_eigen_system.fEigenvectors[i].ZZ();
			
			v2[0]=eps_eigen_system.fEigenvectors[j].XX();
			v2[1]=eps_eigen_system.fEigenvectors[j].XY();
			v2[2]=eps_eigen_system.fEigenvectors[j].XZ();
			v2[3]=eps_eigen_system.fEigenvectors[j].YY();
			v2[4]=eps_eigen_system.fEigenvectors[j].YZ();
			v2[5]=eps_eigen_system.fEigenvectors[j].ZZ();
				
                REAL temp = 2 * G * v1[kj] * v2[ki];
                if (ki == kj) {
                    temp += lambda;
                } else {
                    temp *= 2.;
                }
                for (int l = 0; l < 6; ++l) {
                    const unsigned int li = kival[l];
                    const unsigned int lj = kjval[l];
					
					//dSigDe(l,k) += temp * DecompSig.fEigenvectors[i][l];
					//Tangent(l, k) += gradient(i, j) * temp * eps_eigen_system.fEigenvectors[i][l];
                    Tangent(l, k) += temp * gradient(i, j) * v1[li] * v2[lj];
                }/// l
            }///j
        }///i
    }///k
    cout << "\naqui";
	Tangent.Print("tanzqui");
    REAL deigensig = 0., deigeneps = 0.;
    TPZFNMatrix<9, REAL> tempMat(3, 3, 0.);
    TPZFNMatrix<9, REAL> temp_mat(3, 3, 0.),temp_mat2(3, 3, 0.);
//    TPZFNMatrix<9> ColCorr(3, 3, 0.);
    TPZFNMatrix<6> ColCorrV(6, 1, 0.);
    
    // Correction of the eigenvectors variation
    for (unsigned int i = 0; i < 2; ++i) {
        for (unsigned int j = i + 1; j < 3; ++j) {
            deigeneps = eps_eigen_system.fEigenvalues[i] - eps_eigen_system.fEigenvalues[j];
            deigensig = sig_eigen_system.fEigenvalues[i] - sig_eigen_system.fEigenvalues[j];
    
            REAL factor = 0.;
            if (!IsZero(deigeneps)) {
                factor = deigensig / deigeneps;
            } else {
                factor = fER.G() * (gradient(i, i) - gradient(i, j) - gradient(j, i) + gradient(j, j)); // expression C.20
            }
            
			TPZManVector<REAL,3> v1,v2;
			cout << "eps_eigen_system.fEigenvectors[0] "<<endl;
			eps_eigen_system.fEigenvectors[0].Print(cout);
			//<< endl;
			v1[0]=eps_eigen_system.fEigenvectors[i].XX();
			v1[1]=eps_eigen_system.fEigenvectors[i].XY();
			v1[2]=eps_eigen_system.fEigenvectors[i].YY();
			
			v2[0]=eps_eigen_system.fEigenvectors[j].XX();
			v2[1]=eps_eigen_system.fEigenvectors[j].XY();
			v2[2]=eps_eigen_system.fEigenvectors[j].YY();
			
            ProdT(v1, v2,temp_mat);

            for (unsigned int it = 0; it < 3; ++it) {
                for (unsigned int jt = 0; jt < 3; ++jt) {
                    tempMat(it,jt) += temp_mat(it,jt);
                }
            }
            
			ProdT(v2,v1,tempMat);

            for (unsigned int it = 0; it < 3; ++it) {
                for (unsigned int jt = 0; jt < 3; ++jt) {
                    tempMat(it,jt) += temp_mat2(it,jt);
                }
            }
            
            // expression C.14
            for (unsigned int k = 0; k < 6; ++k) {
                const unsigned int ki = kival[k];
                const unsigned int kj = kjval[k];
                if (ki == kj) {
                    temp_mat = (eps_eigen_system.fEigenvectors[j][ki] * eps_eigen_system.fEigenvectors[i][kj]) * factor * tempMat;
                } else {
                    temp_mat = (eps_eigen_system.fEigenvectors[j][ki] * eps_eigen_system.fEigenvectors[i][kj] + eps_eigen_system.fEigenvectors[j][kj] * eps_eigen_system.fEigenvectors[i][ki]) * factor * tempMat;
                }
                ColCorrV = FromMatToVoight(temp_mat);
                for (int l = 0; l < 6; l++) {
                    Tangent(l, k) += ColCorrV(l, 0);
                }
            }
        } // j
    } // i
    //Tangent=dSigDe;*/
    
	
	
}



template <class YC_t, class ER_t>
void TPZPlasticStepPV<YC_t, ER_t>::TaylorCheck(TPZTensor<REAL> &EpsIni, TPZTensor<REAL> &deps, REAL kprev, TPZVec<REAL> &conv)
{
	TPZTensor<REAL> eps1,eps2, SigmaTemp,Sigma1,Sigma2;
	TPZFNMatrix <36> dSigDe(6,6,0.);
	TPZStack<REAL> coef;
	
	fN.fEpsP.Scale(0.);
	fN.fEpsT.Scale(0.);
	fN.fAlpha = kprev;
	this->ApplyStrainComputeDep(EpsIni,SigmaTemp,dSigDe);
#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout << "EpsIni " << EpsIni << "\nSigmaTemp " << SigmaTemp << "\ndSidDe " << dSigDe << std::endl;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
	fN.fEpsP.Scale(0.);
	fN.fEpsT.Scale(0.);
	fN.fAlpha = kprev;
	
	REAL scale = 1.;
	REAL alphatable[] = {0.1,0.2,0.3,0.4,0.5,0.6};
	for (int i = 0; i < 6; i++) {
		alphatable[i] *= scale;
	}
	for (int ia = 0 ; ia < 5; ia++) {
		REAL alpha1 = alphatable[0];
		REAL alpha2 = alphatable[ia+1];
		eps1.Scale(0.);
		eps2.Scale(0.);
		eps1 = EpsIni;
		eps2 = EpsIni;
		eps1.Add(deps, alpha1);
		eps2.Add(deps, alpha2);
		
		fN.fEpsT = EpsIni;
		this->ApplyStrainComputeSigma(eps1,Sigma1);
		fN.fEpsP.Scale(0.);
		fN.fEpsT.Scale(0.);
		fN.fAlpha = kprev;
		
		fN.fEpsT = EpsIni;
		this->ApplyStrainComputeSigma(eps2,Sigma2);
		fN.fEpsP.Scale(0.);
		fN.fEpsT.Scale(0.);
		fN.fAlpha = kprev;
		
		TPZFNMatrix <6> deps1(6,1,0.),deps2(6,1,0.);
		TPZFNMatrix <9> depsMat(3,3,0.);
		depsMat = deps;
		deps1 = FromMatToVoight(depsMat);
		deps2 = FromMatToVoight(depsMat);
		
		TPZFNMatrix <6> tanmult1(6,1,0.), tanmult2(6,1,0.);
		dSigDe.Multiply(deps1, tanmult1);
		dSigDe.Multiply(deps2, tanmult2);
		
		for (int i = 0 ; i < 6; i++) {
			tanmult1(i,0) *= alpha1;
			tanmult2(i,0) *= alpha2;
		}
		
		TPZFNMatrix <9> SigMatTemp33(3,3,0.);
		TPZFNMatrix <6> sigprMat(6,1,0.),sigpr1Mat(6,1,0.),sigpr2Mat(6,1,0.);
		SigMatTemp33 = SigmaTemp;
		sigprMat = FromMatToVoight(SigMatTemp33);
		SigMatTemp33 = Sigma1;
		sigpr1Mat = FromMatToVoight(SigMatTemp33);
		SigMatTemp33 = Sigma2;
		sigpr2Mat = FromMatToVoight(SigMatTemp33);
		
		TPZFNMatrix<6> error1(6,1,0.), error2(6,1,0.);
#ifdef LOG4CXX
        if (logger->isDebugEnabled()) {
            std::stringstream sout;
            sigprMat.Print("sigprMat",sout);
            sigpr1Mat.Print("sigpr1Mat",sout);
            tanmult1.Print("tanmult1",sout);
            sigpr2Mat.Print("sigpr2Mat",sout);
            tanmult2.Print("tanmult2",sout);
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
		for (int i = 0 ; i < 6; i++) {
			error1(i,0) = sigpr1Mat(i,0) - sigprMat(i,0) - tanmult1(i,0);
			error2(i,0) = sigpr2Mat(i,0) - sigprMat(i,0) - tanmult2(i,0);
		}
#ifdef LOG4CXX
        if (logger->isDebugEnabled()) {
            std::stringstream sout;
            error1.Print("error1:",sout);
            error2.Print("error2:",sout);
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
		REAL n;
		REAL norm1, norm2;
		norm1 = NormVecOfMat(error1);
		norm2 = NormVecOfMat(error2);
		n = ( log(norm1) - log(norm2) ) / ( log(alpha1) - log(alpha2) );
		coef.push_back(n);
	}
	conv = coef;
	std::cout << "coef = " << coef << std::endl;
}

template <class YC_t, class ER_t>
REAL TPZPlasticStepPV<YC_t, ER_t>::ComputeNFromTaylorCheck(REAL alpha1, REAL alpha2, TPZFMatrix<REAL> &error1Mat, TPZFMatrix<REAL> &error2Mat)
{
	REAL norm1, norm2, n;
	norm1 = NormVecOfMat(error1Mat);
	norm2 = NormVecOfMat(error2Mat);
	n = log(norm1/norm2) / log(alpha1/alpha2);
	return n;
}

REAL NormVecOfMat(TPZFNMatrix <9> mat)
{
	REAL norm = 0.;
	for (int i = 0; i < mat.Rows(); i++) {
		norm += mat(i,0) * mat(i,0);
	}
	norm = sqrt(norm);
	return norm;
}

REAL InnerVecOfMat(TPZFMatrix<REAL> &m1,TPZFMatrix<REAL> &m2)
{
	REAL dot = 0.;
	for (int i = 0; i < m1.Rows(); i++) {
		dot += m1(i,0) * m2(i,0);
	}
	return dot;
}

TPZFMatrix<REAL> ProdT(TPZFMatrix<REAL> &m1,TPZFMatrix<REAL> &m2)
{
	TPZFMatrix<REAL> mat(3,3,0.);
	for (int i = 0; i < 3; i++) {
		for (int j = 0 ; j < 3; j++) {
			mat(i,j) = m1(i,0) * m2(j,0);
		}
	}
	return mat;
}
void ProdT(TPZManVector<REAL,3> &v1, TPZManVector<REAL,3> &v2, TPZFMatrix<REAL> & mat) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            mat(i, j) = v1[i] * v2[j];
        }
    }
}
TPZFNMatrix <6> FromMatToVoight(TPZFNMatrix <9> mat)
{
	TPZFNMatrix <6> voi(6,1,0.);
	int k = 0;
	for (int i = 0; i < 3; i++) {
		for (int j = i ; j < 3; j++) {
			voi(k++,0) = mat(i,j);
		}
	}
	return voi;	
}
TPZManVector<REAL,3>  FromToManVec(TPZFNMatrix <9> mat)
{

}



template <class YC_t, class ER_t>
void TPZPlasticStepPV<YC_t, ER_t>::ApplyStrain(const TPZTensor<REAL> &epsTotal)
{

    std::cout<< " \n this method is not implemented in PlasticStepPV. ";
    DebugStop();
    
}

template <class YC_t, class ER_t>
void TPZPlasticStepPV<YC_t, ER_t>::ApplyLoad(const TPZTensor<REAL> & GivenStress, TPZTensor<REAL> &epsTotal)
{
    
    TPZPlasticState<STATE> prevstate=GetState();
    epsTotal=prevstate.fEpsP;
    TPZTensor<STATE> GuessStress,Diff,Diff2,deps;
    TPZFNMatrix<36,STATE> Dep(6,6);
    TPZFNMatrix<6,STATE> GuessStressFN(6,1),DiffFN(6,1);
    
    ApplyStrainComputeDep(epsTotal, GuessStress, Dep);
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        Dep.Print("Dep = ",sout,EMathematicaInput);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    Diff=GivenStress;
    Diff-=GuessStress;
    
    STATE norm = Norm(Diff),normprev;
    STATE tol = 1.e-7;
    int counter = 0;
    TPZVec<STATE> conv;
    
    
    while (norm>tol && counter<30)
    {
        CopyFromTensorToFMatrix(Diff,DiffFN);
        std::list<long> singular;
        Dep.Solve_LU(&DiffFN,singular);
        CopyFromFMatrixToTensor(DiffFN,Diff);
        TPZTensor<STATE> epsprev(epsTotal);
        normprev=norm;
        STATE scale=1.;
        int counter2=0;
        do{
            for(int i=0;i<6;i++)epsTotal.fData[i]=epsprev.fData[i]+scale*Diff.fData[i];
            
            ApplyStrainComputeDep(epsTotal, GuessStress,Dep);
#ifdef LOG4CXX
            if (logger->isDebugEnabled()) {
                std::stringstream sout;
                Dep.Print("Dep = ",sout,EMathematicaInput);
                LOGPZ_DEBUG(logger,sout.str())
            }
#endif
            
            fN=prevstate;
            Diff2=GivenStress;
            Diff2-=GuessStress;
            CopyFromTensorToFMatrix(Diff2, DiffFN);
            Dep.Solve_LU(&DiffFN, singular);
            norm=Norm(Diff2);
            //scale*=0.5;
            counter2++;
        }while (norm>=normprev && counter2<30);
        Diff=Diff2;
        counter++;
        
    }
    ApplyStrainComputeDep(epsTotal, GuessStress,Dep);
 
}

template <class YC_t, class ER_t>
TPZPlasticState<STATE>  TPZPlasticStepPV<YC_t, ER_t>::GetState() const
{
    return fN;
}

template <class YC_t, class ER_t>
void TPZPlasticStepPV<YC_t, ER_t>::Phi(const TPZTensor<STATE> &eps, TPZVec<REAL> &phi) const
{
    TPZTensor<STATE> sigma;
    fER.Compute(eps, sigma);
    TPZTensor<STATE>::TPZDecomposed DecSig;
    sigma.EigenSystem(DecSig);
    TPZVec<STATE> sigvec(3);
    sigvec[0]=DecSig.fEigenvalues[0];
    sigvec[1]=DecSig.fEigenvalues[1];
    sigvec[2]=DecSig.fEigenvalues[2];
    fYC.Phi(sigvec,fN.Alpha(),phi);
}


template <class YC_t, class ER_t>
void TPZPlasticStepPV<YC_t, ER_t>::SetState(const TPZPlasticState<REAL> &state)
{
    fN=state;
}



template <class YC_t, class ER_t>
void TPZPlasticStepPV<YC_t, ER_t>::Read(TPZStream &buf)
{
	fYC.Read(buf);
	fER.Read(buf);
	buf.Read(&fResTol);
	buf.Read(&fMaxNewton);
	fN.Read(buf);
}

/** @brief Object which represents the yield criterium */
//YC_t fYC;

/** @brief Object representing the elastic response */
//ER_t fER;

/** @brief Residual tolerance accepted in the plastic loop processes */
//REAL fResTol;

/** @brief Maximum number of Newton interations allowed in the nonlinear solvers */
//int fMaxNewton;	// COLOCAR = 30 (sugestao do erick!)




/** @brief Plastic State Variables (EpsT, EpsP, Alpha) at the current time step */
//TPZPlasticState<STATE> fN;



template <class YC_t, class ER_t>
void TPZPlasticStepPV<YC_t, ER_t>::Write(TPZStream &buf) const
{
    fYC.Write(buf);
    fER.Write(buf);
    buf.Write(&fResTol);
    buf.Write(&fMaxNewton);
    fN.Write(buf);
    
}

template <class YC_t, class ER_t>
void TPZPlasticStepPV<YC_t, ER_t>::CopyFromFMatrixToTensor(TPZFMatrix<STATE> FNM,TPZTensor<STATE> &copy)
{
    FNM.Resize(6,1);
    copy.XX()=FNM(0,0);
    copy.XY()=FNM(1,0);
    copy.XZ()=FNM(2,0);
    copy.YY()=FNM(3,0);
    copy.YZ()=FNM(4,0);
    copy.ZZ()=FNM(5,0);
}

template <class YC_t, class ER_t>
void TPZPlasticStepPV<YC_t, ER_t>::CopyFromTensorToFMatrix(TPZTensor<STATE> tensor,TPZFMatrix<STATE> &copy)
{
    
    copy(0,0)=tensor.XX();
    copy(1,0)=tensor.XY();
    copy(2,0)=tensor.XZ();
    copy(3,0)=tensor.YY();
    copy(4,0)=tensor.YZ();
    copy(5,0)=tensor.ZZ();
}
// template <class YC_t, class ER_t>
// 	REAL TPZPlasticStepPV<YC_t, ER_t>::NormVecOfMat(TPZFNMatrix <9> mat)
// {
//     REAL norm = 0.;
//     for (int i = 0; i < mat.Rows(); i++) {
//         norm += mat(i, 0) * mat(i, 0);
//     }
//     norm = sqrt(norm);
//     return norm;
// }
// template <class YC_t, class ER_t>
// REAL TPZPlasticStepPV<YC_t, ER_t>::InnerVecOfMat(TPZFMatrix<REAL> &m1,TPZFMatrix<REAL> &m2)
// {
//     REAL dot = 0.;
//     for (int i = 0; i < m1.Rows(); i++) {
//         dot += m1(i, 0) * m2(i, 0);
//     }
//     return dot;
// }
// template <class YC_t, class ER_t>
// TPZFMatrix<REAL> TPZPlasticStepPV<YC_t, ER_t>::ProdT(TPZManVector<REAL,3> &v1, TPZManVector<REAL,3> &v2) {
//     TPZFMatrix<REAL> mat(3, 3, 0.);
//     for (int i = 0; i < 3; i++) {
//         for (int j = 0; j < 3; j++) {
//             mat(i, j) = v1[i] * v2[j];
//         }
//     }
//     return mat;
// }
// template <class YC_t, class ER_t>
// void TPZPlasticStepPV<YC_t, ER_t>::ProdT(TPZManVector<REAL,3> &v1, TPZManVector<REAL,3> &v2, TPZFMatrix<REAL> & mat) {
//     for (int i = 0; i < 3; i++) {
//         for (int j = 0; j < 3; j++) {
//             mat(i, j) = v1[i] * v2[j];
//         }
//     }
// }
// template <class YC_t, class ER_t>
// TPZFNMatrix <6> TPZPlasticStepPV<YC_t, ER_t>::FromMatToVoight(TPZFNMatrix <9> mat)
// {
//     TPZFNMatrix <6> voi(6, 1, 0.);
//     int k = 0;
//     for (int i = 0; i < 3; i++) {
//         for (int j = i; j < 3; j++) {
//             voi(k++, 0) = mat(i, j);
//         }
//     }
//     return voi;
// }

template class TPZPlasticStepPV<TPZSandlerExtended, TPZElasticResponse>;
template class TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>;
//template class TPZPlasticStepPV<TPZMohrCoulombVoigt, TPZElasticResponse>;

/*
 // Correcao do giro rigido
 for (int i = 0; i < 2; i++) {
 for (int j = i+1; j<3 ; j++) {
 deigeneps = DecompEps.fEigenvalues[i]  - DecompEps.fEigenvalues[j];
 deigensig = sigprvec[i] - sigprvec[j];
 TPZFNMatrix<9,REAL> tempMat(3,1,0.);
 depsMat.Multiply(epsegveFromProj[i], tempMat);
 REAL deij = InnerVecOfMat(tempMat,epsegveFromProj[j]);
 REAL factor = 0.;
 if (!IsZero(deigeneps)) {
 factor = deigensig * deij / deigeneps;
 }
 else {
 factor = fER.G() * ( GradSigma(i,i) - GradSigma(i,j) - GradSigma(j,i) + GradSigma(j,j) ) * deij;
 }
 std::cout << "factor = " << factor << std::endl;
 std::cout << "G = " << fER.G() << std::endl;
 GradSigma.Print("GradSigma");
 tempMat.Redim(3, 3);
 tempMat = ProdT(epsegveFromProj[i],epsegveFromProj[j]) + ProdT(epsegveFromProj[j],epsegveFromProj[i]);
 factorMat += tempMat * factor;
 
 }
 }
*/
