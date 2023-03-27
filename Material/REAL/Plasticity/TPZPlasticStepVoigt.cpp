#include "TPZPlasticStepVoigt.h"
#include "TPZElasticResponse.h"


#include "TPZMohrCoulombVoigt.h"
#include "TPZElasticResponse.h"


template <class YC_t, class ER_t>
void TPZPlasticStepVoigt<YC_t, ER_t>::ApplyStrainComputeSigma(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma)
{
    std::cout<< " \n this method is not implemented in TPZPlasticStepVoigt. ";
    DebugStop();
}

template <class YC_t, class ER_t>
void TPZPlasticStepVoigt<YC_t, ER_t>::ApplyStrainComputeDep(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma, TPZFMatrix<REAL> &Dep)
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
    TPZFMatrix<REAL> tempeps,tempsig,sigprojtemp,epstemp;
    //FromTensorToMatVoigt(epsTr,tempeps);
    epsTr.FromTensorToNRmatrix(tempeps);

    //fER.ComputeStress(epsTr,sigtr);
    fER.Compute(epsTr,sigtr);

    // ReturMap in the principal values
    STATE nextalpha = -0.;
    STATE printPlastic = fN.Alpha();

    fYC.ProjectSigmaDep(sigtr, sigma,Dep, nextalpha);


    //cout << " ------------------------------------------Dep " << endl;
    //Dep.Print(std::cout);

    sigma.FromTensorToNRmatrix(sigprojtemp);
    //FromTensorToMatVoigt(sigma,sigprojtemp);

    fN.fAlpha = nextalpha;

    fER.ComputeDeformation(sigma,epsElaNp1);
    //fER.ComputeStrain(sigma,epsElaNp1);

    fN.fEpsT = epsTotal;
    epsPN = epsTotal;
    epsPN -= epsElaNp1; // Transforma epsPN em epsPNp1
    fN.fEpsP = epsPN;
}

template <class YC_t, class ER_t>
void TPZPlasticStepVoigt<YC_t, ER_t>::TangentOperator(TPZFMatrix<REAL> & gradient,TPZTensor<REAL>::TPZDecomposed & eps_eigen_system, TPZTensor<REAL>::TPZDecomposed & sig_eigen_system, TPZFMatrix<REAL> & Tangent){
	
    std::cout<< " \n this method is not implemented in TPZPlasticStepVoigt. ";
    DebugStop();
	
}


template <class YC_t, class ER_t>
void TPZPlasticStepVoigt<YC_t, ER_t>::TaylorCheck(TPZTensor<REAL> &EpsIni, TPZTensor<REAL> &deps, REAL kprev, TPZVec<REAL> &conv)
{
    std::cout<< " \n this method is not implemented in TPZPlasticStepVoigt. ";
    DebugStop();
}


template <class YC_t, class ER_t>
void TPZPlasticStepVoigt<YC_t, ER_t>::ApplyStrain(const TPZTensor<REAL> &epsTotal)
{

    std::cout<< " \n this method is not implemented in TPZPlasticStepVoigt. ";
    DebugStop();
    
}

template <class YC_t, class ER_t>
void TPZPlasticStepVoigt<YC_t, ER_t>::ApplyLoad(const TPZTensor<REAL> & GivenStress, TPZTensor<REAL> &epsTotal)
{
     std::cout<< " \n this method is not implemented in TPZPlasticStepVoigt. ";
    DebugStop();
}

template <class YC_t, class ER_t>
TPZPlasticState<STATE>  TPZPlasticStepVoigt<YC_t, ER_t>::GetState() const
{
    return fN;
}

template <class YC_t, class ER_t>
void TPZPlasticStepVoigt<YC_t, ER_t>::Phi(const TPZTensor<STATE> &eps, TPZVec<REAL> &phi) const
{
    std::cout<< " \n this method is not implemented in TPZPlasticStepVoigt. ";
    DebugStop();
}


template <class YC_t, class ER_t>
void TPZPlasticStepVoigt<YC_t, ER_t>::SetState(const TPZPlasticState<REAL> &state)
{
    fN=state;
}



template <class YC_t, class ER_t>
void TPZPlasticStepVoigt<YC_t, ER_t>::Read(TPZStream &buf)
{
	fYC.Read(buf);
	fER.Read(buf);
	buf.Read(&fResTol);
	buf.Read(&fMaxNewton);
	fN.Read(buf);
}



template <class YC_t, class ER_t>
void TPZPlasticStepVoigt<YC_t, ER_t>::Write(TPZStream &buf) const
{
    fYC.Write(buf);
    fER.Write(buf);
    buf.Write(&fResTol);
    buf.Write(&fMaxNewton);
    fN.Write(buf);
    
}

// template <class YC_t, class ER_t>
// void TPZPlasticStepVoigt<YC_t, ER_t>::CopyFromFMatrixToTensor(TPZFMatrix<STATE> FNM,TPZTensor<STATE> &copy)
// {
//     //FNM.Resize(6,1);
//     copy.XX()=FNM(0,0);
//     copy.XY()=FNM(1,0);
//     copy.XZ()=FNM(2,0);
//     copy.YY()=FNM(3,0);
//     copy.YZ()=FNM(4,0);
//     copy.ZZ()=FNM(5,0);
// }
//
// template <class YC_t, class ER_t>
// void TPZPlasticStepVoigt<YC_t, ER_t>::CopyFromTensorToFMatrix(TPZTensor<STATE> tensor,TPZFMatrix<STATE> &copy)
// {
//     copy.Resize(6,1);
//     copy(0,0)=tensor.XX();
//     copy(1,0)=tensor.XY();
//     copy(2,0)=tensor.XZ();
//     copy(3,0)=tensor.YY();
//     copy(4,0)=tensor.YZ();
//     copy(5,0)=tensor.ZZ();
// }


template class TPZPlasticStepVoigt<TPZMohrCoulombVoigt, TPZElasticResponse>;
