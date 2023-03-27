/*
 *  TPZMohrCoulomb.h
 *  FEMPZ
 *
 *  Created by Nathan Shauer on 5/4/13.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef TPZYCMOHRCOULOMBPV_H
#define TPZYCMOHRCOULOMBPV_H

#include "pzlog.h"
#include "TPZTensor.h"
#include "pzvec_extras.h"
#include "pzsave.h"
#include "TPZPlasticState.h"
#include "TPZElasticResponse.h"

#ifdef LOG4CXX
static LoggerPtr loggerMohrCoulombPV ( Logger::getLogger ( "pz.plasticity.mohrcoulombpv" ) );
#endif

class TPZYCMohrCoulombPV
{

public:
    enum {NYield=3};

private:
    REAL fPhi;
    REAL fPsi;
    REAL fc;
    TPZElasticResponse fER;
    REAL fa;
    REAL fda;
    REAL fd2a;
    REAL ftheta;


protected:
    REAL fEpsPlasticBar;

public:

    /// structure which contains the decision tree of the return map
    // we can only expect a consistent tangent matrix if the decision tree remains the same
    struct TComputeSequence {
        TComputeSequence() : fWhichPlane ( ENoPlane ), fGamma ( 0 )
        {

        }

        TComputeSequence ( const TComputeSequence &copy ) : fWhichPlane ( copy.fWhichPlane ), fGamma ( copy.fGamma )
        {

        }

        TComputeSequence &operator= ( const TComputeSequence &copy )
        {
            fWhichPlane = copy.fWhichPlane;
            fGamma = copy.fGamma;
            return *this;
        }

        enum MPlane {ENoPlane, EElastic, EMainPlane, ERightEdge, ELeftEdge, EApex };

        MPlane fWhichPlane;

        TPZManVector<REAL> fGamma;
    };

public:

    /**
     * @brief empty constructor
     */
    TPZYCMohrCoulombPV();

    /**
     * @brief Constructor seting yc parameters
     */
    TPZYCMohrCoulombPV ( REAL Phi, REAL Psi, REAL c, TPZElasticResponse &ER );

    /**
     * @brief Copy Constructor
     */
    TPZYCMohrCoulombPV ( const TPZYCMohrCoulombPV &cp );

    /**
     * @brief Sets up the data
     */
    void SetUp ( REAL Phi, REAL Psi, REAL c, TPZElasticResponse &ER )
    {
        fPhi = Phi;
        fPsi = Psi;
        fc = c;
        fER = ER;

    }

    void SetLocalMatState ( TPZPlasticState<REAL> & state )
    {
        //if ( fc<1.e-3 ) DebugStop();
        fc =   state.fmatprop[0];
        fPhi = state.fmatprop[1];
        fPsi = state.fmatprop[1];
	//	std::cout << "fc = "<< fc <<endl;
    }

    /**
     * @brief Operator =
     */
    TPZYCMohrCoulombPV & operator= ( const TPZYCMohrCoulombPV &cp );


    void Read ( TPZStream &buf );

    void Write ( TPZStream &buf ) const;


    /**
     * @brief Sets epsbar
     */
    void SetEpsBar ( REAL &epsbar )
    {
        fEpsPlasticBar = epsbar;
    }

    /**
     * @brief Print Method
     */
    void Print ( std::ostream &out ) const
    {
        out << "TPZYCMohrCoulombPV\n";
        out << "Still have to implement the print" << std::endl;
    }

    /**
     * @brief Calculates the value c(epsp) and its derivative
     */
    template <class T>
    void PlasticityFunction ( const T epsp, T &c, T &H ) const;

    /**
     * @brief sigma = lambda Tr(E)I + 2 mu E
     */
    template<class T>
    TPZVec<T> SigmaElastPV ( const TPZVec<T> &deform ) const;

    /**
     * @brief Calcula o valor da funcao criteiro de plastificacao
     */
    template<class T>
    T PhiPlane ( const TPZVec<T> &sigma ) const;

    /**
     * @brief Implements the return map in the plane of the surface
     */
    template<class T>
    bool ReturnMapPlane ( const TPZVec<T> &sigma_trial, TPZVec<T> &sigma_projected,
                          TComputeSequence &memory, REAL &epsbarnew ) const;

    /**
     * @brief Computes dsigmapr/dsigmatr for the ReturnMapPlane
     */
    void ComputePlaneTangent ( TPZMatrix<REAL> &tang, REAL &epsbarp ) const;

    /**
     * @brief Implements the return map in the left edge of the surface
     */
    template<class T>
    bool ReturnMapLeftEdge ( const TPZVec<T> &sigma_trial, TPZVec<T> &sigma_projected,
                             TComputeSequence &memory, REAL &epsbarnew ) const;

    /**
     * @brief Computes dsigmapr/dsigmatr for the ReturnMapLeftEdge
     */
    void ComputeLeftEdgeTangent ( TPZMatrix<REAL> &tang, REAL &epsbarp ) const;

    /**
     * @brief Implements the return map in the right edge of the surface
     */
    template<class T>
    bool ReturnMapRightEdge ( const TPZVec<T> &sigma_trial, TPZVec<T> &sigma_projected,
                              TComputeSequence &memory, REAL &epsbarnew ) const;

    /**
     * @brief Computes dsigmapr/dsigmatr for the ReturnMapRightEdge
     */
    void ComputeRightEdgeTangent ( TPZMatrix<REAL> &tang, REAL &epsbarp ) const;

    /**
     * @brief Implements the return map in the apex
     */
    template<class T>
    bool ReturnMapApex ( const TPZVec<T> &sigma_trial, TPZVec<T> &sigma_projected,
                         TComputeSequence &memory, REAL &epsbarnew ) const;

    /**
     * @brief Computes dsigmapr/dsigmatr for the ReturnMapApex
     */
    void ComputeApexTangent ( TPZMatrix<REAL> &tang, REAL &epsbarp ) const;

    /**
     * @brief Choses the correct projection and returns projected sigma and new epspbar
     */
    void ProjectSigma ( const TPZVec<STATE> &sigma_trial, STATE eprev, TPZVec<STATE> &sigma, STATE &eproj );

    /**
     * @brief Choses the correct projection and returns projected sigma, new epspbar and tangent matrix
     */
    void ProjectSigmaDep ( const TPZVec<STATE> &sigmatrial, TPZTensor<REAL>sigtrtensor,STATE kprev, TPZVec<STATE> &sigmaproj, STATE &kproj, TPZFMatrix<STATE> &tang );

    /**
     * @brief Calculates the value of phi based on eps
     */
    void Phi ( TPZVec<STATE> sigvec,STATE alpha,TPZVec<STATE> &phi ) const;

    STATE Phi()
    {
        return fPhi;
    }
    STATE Psi()
    {
        return fPsi;
    }
    STATE Cohesion()
    {
        return fc;
    }
    STATE E()
    {
        return fER.E();
    }
    STATE Poisson()
    {
        return fER.Poisson();
    }

    REAL PhiInvars ( TPZTensor<REAL> tensor )
    {
        REAL I1=tensor.I1();
        REAL J2=tensor.J2();
        REAL thetaval = theta ( tensor );
        REAL atheta = fa;//A(tensor,thetaval);
        REAL thetaint =theta ( tensor );
        REAL s1=2.*sqrt ( J2 ) /sqrt ( 3. ) * ( sin ( thetaint+2.*M_PI/3. ) )+I1/3.;
        //cout << "s1= "<< s1<< endl;
        return 1./3. *I1*sin ( fPhi )+sqrt ( J2 ) *atheta-fc*cos ( fPhi );
    }

    REAL theta ( TPZTensor<REAL> tensor )
    {
        REAL J2=tensor.J2();
        REAL J3=tensor.J3();
        REAL val = -3*sqrt ( 3. ) *J3/ ( 2.*pow ( J2,1.5 ) );
        if ( val>1. ) {
            val=1.;
        }
        if ( val<-1. ) {
            val=-1.;
        }
        return 1/3.*asin ( val );
    }
    void ComputeConsistentPlaneTangent ( TPZTensor<REAL> &trialstress,TPZFMatrix<REAL> & Dep );
    void ComputeConsistentEdgeTangent ( TPZTensor<REAL> &trialstress,TPZFMatrix<REAL> & Dep, int &m_type );


//     TPZFMatrix<REAL> GetElasticMatrix()
//     {
//         TPZFMatrix<REAL> C ( 6, 6, 0. );
//         REAL G = fER.G();
//         REAL K = fER.K();
//         C ( _XX_,_XX_ ) = ( 4 * G ) / 3 + K;
//         C ( _XX_,_YY_ ) = - ( ( 2 * G ) / 3 ) + K;
//         C ( _XX_,_ZZ_ ) = - ( ( 2 * G ) / 3 ) + K;
//         C ( _XX_,_XZ_ ) = 0.;
//         C ( _XY_,_YZ_ ) = 0.;
//         C ( _XX_,_XY_ ) = 0.;
//
//
//         C ( _YY_,_XX_ ) = - ( ( 2 * G ) / 3 ) + K;
//         C ( _YY_,_YY_ ) = ( 4 * G ) / 3 + K;
//         C ( _YY_,_ZZ_ ) = - ( ( 2 * G ) / 3 ) + K;
//         C ( _YY_,_XZ_ ) = 0.;
//         C ( _YY_,_YZ_ ) = 0.;
//         C ( _YY_,_XY_ ) = 0.;
//
//
//         C ( _ZZ_,_XX_ ) = - ( ( 2 * G ) / 3 ) + K;
//         C ( _ZZ_,_YY_ ) = - ( ( 2 * G ) / 3 ) + K;
//         C ( _ZZ_,_ZZ_ ) = ( 4 * G ) / 3 + K;
//         C ( _ZZ_,_XZ_ ) = 0.;
//         C ( _ZZ_,_YZ_ ) = 0.;
//         C ( _ZZ_,_XY_ ) = 0.;
//
//
//         C ( _XZ_,_XX_ ) = 0;
//         C ( _XZ_,_YY_ ) = 0;
//         C ( _XZ_,_ZZ_ ) = 0;
//         C ( _XZ_,_XZ_ ) = G ;
//         C ( _XZ_,_YZ_ ) = 0.;
//         C ( _XZ_,_XY_ ) = 0.;
//
//
//         C ( _YZ_,_XX_ ) = 0;
//         C ( _YZ_,_YY_ ) = 0;
//         C ( _YZ_,_ZZ_ ) = 0;
//         C ( _YZ_,_XZ_ ) = 0.;
//         C ( _YZ_,_YZ_ ) = G ;
//         C ( _YZ_,_XY_ ) = 0.;
//
//
//         C ( _XY_,_XX_ ) = 0;
//         C ( _XY_,_YY_ ) = 0;
//         C ( _XY_,_ZZ_ ) = 0;
//         C ( _XY_,_XZ_ ) = 0.;
//         C ( _XY_,_YZ_ ) = 0.;
//         C ( _XY_,_XY_ ) = G ;
//
//
//         return C;
//     }


    TPZFMatrix<REAL> avec ( TPZTensor<REAL> tensor )
    {
        //Eq. 14 (Crisfield Eng. Comput., 1987)
        REAL c1=C1(),c2=C2 ( tensor ),c3=C3 ( tensor );
        TPZFMatrix<REAL> a1 ( 6,1 ),a2 ( 6,1 ),a3 ( 6,1 ),a;
        tensor.dI1().FromTensorToNRmatrix ( a1 );
        tensor.dJ2().FromTensorToNRmatrix ( a2 );
        tensor.dJ3().FromTensorToNRmatrix ( a3 );

		//tensor.Print(cout);
		//cout << "c1" << c1 << endl;
		//cout << "c2" << c2 << endl;
		//cout << "c3" << c3 << endl;
		//a1.Print("a1");
		//a2.Print("a2");
		//a3.Print("a3");
        a1*=c1;
        a2*=c2;
        a3*=c3;
        a=a1;
        a+=a2;
        a+=a3;

        return a;
    }

    REAL C1()
    {
        //return sin ( fPhi ) /3.;
        if ( fabs ( ftheta ) <29.99*M_PI/180 ) {
            return sin ( fPhi ) /3.;
        } else {
            return 0;
        }

    }
    REAL C2 ( TPZTensor<REAL> tensor )
    {
        REAL thetaval = ftheta;
        REAL J2=tensor.J2();
        REAL da=fda;//dA(tensor,thetaval);

        REAL a=fa;//A(tensor,thetaval);
        //return 1/2. * pow(J2,-0.5)*(a-tan(3*thetaval)*da);
        //cout << "thetaval = "<< thetaval << endl;
        if ( fabs ( ftheta ) <29.99*M_PI/180 ) {
            return 1/2. * pow ( J2,-0.5 ) * ( a-tan ( 3*thetaval ) *da );
        } else {
            return 1/2. * pow ( J2,-0.5 ) *a;
        }
    }
    REAL C3 ( TPZTensor<REAL> tensor )
    {
        REAL thetaval = ftheta;
        REAL J2=tensor.J2();
        REAL da=fda;//dA(tensor,thetaval);
        REAL a=fa;//A(tensor,thetaval);
        //return -sqrt(3.)*da/(2.*J2*cos(3*thetaval));
        if ( fabs ( ftheta ) <29.99*M_PI/180 ) {
            return -sqrt ( 3. ) *da/ ( 2.*J2*cos ( 3*thetaval ) );
        } else {
            return 0;
        }
    }

    REAL C23 ( TPZTensor<REAL> tensor )
    {
        REAL thetaval = ftheta;
        REAL J2=tensor.J2();
        REAL da=fda;//dA(tensor,thetaval);
        REAL a=fa;//A(tensor,thetaval);

        REAL c4 = C4 ( tensor );
        return ( 1./2.*tan ( 3.*thetaval ) *c4+da ) *sqrt ( 3. ) / ( 2*J2*J2*cos ( 3.*thetaval ) );
    }
    REAL C22 ( TPZTensor<REAL> tensor )
    {
        REAL thetaval = ftheta;
        REAL J2=tensor.J2();
        REAL da=fda;//dA(tensor,thetaval);

        REAL a=fa;//A(tensor,thetaval);
        REAL c4 = C4 ( tensor );

        return - ( a-pow ( tan ( 3*thetaval ),2 ) *c4-3.*tan ( 3.*thetaval ) *da ) / ( 4.*pow ( J2,3./2. ) );
    }
    REAL C33 ( TPZTensor<REAL> tensor )
    {
        REAL J2=tensor.J2();
        REAL thetaval = ftheta;
        REAL c4 = C4 ( tensor );
        return 3.*c4/ ( 4.*pow ( J2,5/2. ) *cos ( 3.*thetaval ) *cos ( 3.*thetaval ) );
    }
    REAL C4 ( TPZTensor<REAL> tensor )
    {
        REAL thetaval = ftheta;
        REAL da=fda;//dA(tensor,thetaval);
        REAL d2a=fd2a;//d2A(tensor,thetaval);
        if ( fabs ( ftheta ) <29.99*M_PI/180 ) {
            return d2a+3.*tan ( 3.*thetaval ) *da;
        } else {
            return d2a;
        }
    }


    TPZFMatrix<REAL>  dAdsig ( TPZTensor<REAL> tensor )
    {
        //Eq. 14 (Crisfield Eng. Comput., 1987)
        REAL c1=C1(),c2=C2 ( tensor ),c3=C3 ( tensor ),c4=C4 ( tensor ),c23=C23 ( tensor ),c22=C22 ( tensor ),c33=C33 ( tensor );
		//cout << "c1 " << c1<<endl;
		//cout << "c2 "<<c2<<endl;
		//cout << "c3 "<< c3 << endl;
		//cout << "c4 "<< c4<<endl;
		//cout << "c23 "<< c23<<endl;
		//cout << "c22 "<< c22 <<endl;
		//cout << "c33 "<< c33 <<endl;
        REAL c32 =c23;
        TPZFMatrix<REAL> da2dsig = d2J2d2sig();
        TPZFMatrix<REAL> da3dsig = d2J3d2sig ( tensor );
				//cout << "da2dsig"<<endl;
        //da2dsig.Print(cout);
				//		cout << "da3dsig"<<endl;
       // da3dsig.Print(cout);
        TPZFMatrix<REAL> dadsig ( 6,6 ),a2 ( 6,1 ),a3 ( 6,1 ),a2t,a3t,temp1,temp2,temp3,temp4;
        //cout << "dj2"<<endl;
        //tensor.dJ2().Print(cout);
        tensor.dJ2().FromTensorToNRmatrix ( a2 );
        //cout << "dj2"<<endl;
       // a2.Print(cout);
        tensor.dJ3().FromTensorToNRmatrix ( a3 );
		//cout << "a3"<<endl;
        //a3.Print(cout);
        a2.Transpose ( &a2t );
        a3.Transpose ( &a3t ); //tranposto é deitado
        da2dsig*=c2;
        da3dsig*=c3;

        a2.Multiply ( a2t,temp1 );
        temp1*=c22;

        a2.Multiply ( a3t,temp2 );
        temp2*=c23;

        a3.Multiply ( a2t,temp3 );
        temp3*=c32;

        a3.Multiply ( a3t,temp4 );
        temp4*=c33;

        dadsig=da2dsig;
        dadsig+=da3dsig;
        dadsig+=temp1;
        dadsig+=temp2;
        dadsig+=temp3;
        dadsig+=temp4;

        return dadsig;
    }


    TPZFMatrix<REAL> d2J2d2sig()
    {
        TPZFMatrix<REAL>  P ( 6, 6 );
        P ( _XX_,_XX_ ) = 2. / 3.;
        P ( _XX_,_YY_ ) = -1. / 3.;
        P ( _XX_,_ZZ_ ) = -1. / 3.;
        P ( _XX_,_XZ_ ) = 0.;
        P ( _XX_,_YZ_ ) = 0.;
        P ( _XX_,_XY_ ) = 0.;


        P ( _YY_,_XX_ ) = -1. / 3.;
        P ( _YY_,_YY_ ) = 2. / 3.;
        P ( _YY_,_ZZ_ ) = -1. / 3.;
        P ( _YY_,_XZ_ ) = 0.;
        P ( _YY_,_YZ_ ) = 0.;
        P ( _YY_,_XY_ ) = 0.;


        P ( _ZZ_,_XX_ ) = -1. / 3.;
        P ( _ZZ_,_YY_ ) = -1. / 3.;
        P ( _ZZ_,_ZZ_ ) = 2. / 3.;
        P ( _ZZ_,_XZ_ ) = 0.;
        P ( _ZZ_,_YZ_ ) = 0.;
        P ( _ZZ_,_XY_ ) = 0.;


        P ( _XZ_,_XX_ ) = 0.;
        P ( _XZ_,_YY_ ) = 0.;
        P ( _XZ_,_ZZ_ ) = 0.;
        P ( _XZ_,_XZ_ ) = 2.;
        P ( _XZ_,_YZ_ ) = 0.;
        P ( _XZ_,_XY_ ) = 0.;


        P ( _YZ_,_XX_ ) = 0.;
        P ( _YZ_,_YY_ ) = 0.;
        P ( _YZ_,_ZZ_ ) = 0.;
        P ( _YZ_,_XZ_ ) = 0.;
        P ( _YZ_,_YZ_ ) = 2.;
        P ( _YZ_,_XY_ ) = 0.;


        P ( _XY_,_XX_ ) = 0.;
        P ( _XY_,_YY_ ) = 0.;
        P ( _XY_,_ZZ_ ) = 0.;
        P ( _XY_,_XZ_ ) = 0.;
        P ( _XY_,_YZ_ ) = 0.;
        P ( _XY_,_XY_ ) = 2.;

        return P;
    }

    TPZFMatrix<REAL>   d2J3d2sig ( TPZTensor<REAL> tensor )
    {

        TPZFMatrix<REAL>  P ( 6, 6 ),P2 ( 6,6 );
        REAL sxx =tensor.XX();
        REAL syy=tensor.YY();
        REAL szz=tensor.ZZ();
        REAL sxy=tensor.XY();
        REAL syz=tensor.YZ();
        REAL sxz=tensor.XZ();

        P ( 0,0 ) = ( 4*sxx - 2* ( syy + szz ) ) /9.;
        P ( 0,1 ) = ( -2*sxx - 2*syy + 4*szz ) /9.;
        P ( 0,2 ) = ( -2*sxx + 4*syy - 2*szz ) /9.;
        P ( 0,3 ) = ( 2*sxz ) /3.,P ( 0,4 ) = ( -4*syz ) /3.;
        P ( 0,5 ) = ( 2*sxy ) /3.;

        P2 ( _XX_,_XX_ ) = P ( 0,0 );
        P2 ( _XX_,_YY_ ) = P ( 0,1 );
        P2 ( _XX_,_ZZ_ ) = P ( 0,2 );
        P2 ( _XX_,_XZ_ ) = P ( 0,3 );
        P2 ( _XX_,_YZ_ ) = P ( 0,4 );
        P2 ( _XX_,_XY_ ) = P ( 0,5 );

        P ( 1,0 ) = ( -2*sxx - 2*syy + 4*szz ) /9.;
        P ( 1,1 ) = ( -2*sxx + 4*syy - 2*szz ) /9.;
        P ( 1,2 ) = ( 4*sxx - 2*syy - 2*szz ) /9.;
        P ( 1,3 ) = ( -4*sxz ) /3.,P ( 1,4 ) = ( 2*syz ) /3.;
        P ( 1,5 ) = ( 2*sxy ) /3.;

        P2 ( _YY_,_XX_ ) = P ( 1,0 );
        P2 ( _YY_,_YY_ ) = P ( 1,1 );
        P2 ( _YY_,_ZZ_ ) = P ( 1,2 );
        P2 ( _YY_,_XZ_ ) = P ( 1,3 );
        P2 ( _YY_,_YZ_ ) = P ( 1,4 );
        P2 ( _YY_,_XY_ ) = P ( 1,5 );

        P ( 2,0 ) = ( -2*sxx + 4*syy - 2*szz ) /9.;
        P ( 2,1 ) = ( 4*sxx - 2*syy - 2*szz ) /9.;
        P ( 2,2 ) = ( -2*sxx - 2*syy + 4*szz ) /9.;
        P ( 2,3 ) = ( 2*sxz ) /3.,P ( 2,4 ) = ( 2*syz ) /3.;
        P ( 2,5 ) = ( -4*sxy ) /3.;

        P2 ( _ZZ_,_XX_ ) = P ( 2,0 );
        P2 ( _ZZ_,_YY_ ) = P ( 2,1 );
        P2 ( _ZZ_,_ZZ_ ) = P ( 2,2 );
        P2 ( _ZZ_,_XZ_ ) = P ( 2,3 );
        P2 ( _ZZ_,_YZ_ ) = P ( 2,4 );
        P2 ( _ZZ_,_XY_ ) = P ( 2,5 );

        P ( 3,0 ) = ( 2*sxz ) /3.;
        P ( 3,1 ) = ( -4*sxz ) /3.;
        P ( 3,2 ) = ( 2*sxz ) /3.;
        P ( 3,3 ) = ( 2* ( sxx - 2*syy + szz ) ) /3.,P ( 3,4 ) =2*sxy;
        P ( 3,5 ) =2*syz;

        P2 ( _XZ_,_XX_ ) = P ( 3,0 );
        P2 ( _XZ_,_YY_ ) = P ( 3,1 );
        P2 ( _XZ_,_ZZ_ ) = P ( 3,2 );
        P2 ( _XZ_,_XZ_ ) = P ( 3,3 );
        P2 ( _XZ_,_YZ_ ) = P ( 3,4 );
        P2 ( _XZ_,_XY_ ) = P ( 3,5 );

        P ( 4,0 ) = ( -4*syz ) /3.;
        P ( 4,1 ) = ( 2*syz ) /3.;
        P ( 4,2 ) = ( 2*syz ) /3.;
        P ( 4,3 ) =2*sxy,P ( 4,4 ) = ( 2* ( -2*sxx + syy + szz ) ) /3.;
        P ( 4,5 ) =2*sxz;

        P2 ( _YZ_,_XX_ ) = P ( 4,0 );
        P2 ( _YZ_,_YY_ ) = P ( 4,1 );
        P2 ( _YZ_,_ZZ_ ) = P ( 4,2 );
        P2 ( _YZ_,_XZ_ ) = P ( 4,3 );
        P2 ( _YZ_,_YZ_ ) = P ( 4,4 );
        P2 ( _YZ_,_XY_ ) = P ( 4,5 );

        P ( 5,0 ) = ( 2*sxy ) /3.;
        P ( 5,1 ) = ( 2*sxy ) /3.;
        P ( 5,2 ) = ( -4*sxy ) /3.;
        P ( 5,3 ) =2*syz,P ( 5,4 ) =2*sxz;
        P ( 5,5 ) = ( 2* ( sxx + syy - 2*szz ) ) /3.;

        P2 ( _XY_,_XX_ ) = P ( 5,0 );
        P2 ( _XY_,_YY_ ) = P ( 5,1 );
        P2 ( _XY_,_ZZ_ ) = P ( 5,2 );
        P2 ( _XY_,_XZ_ ) = P ( 5,3 );
        P2 ( _XY_,_YZ_ ) = P ( 5,4 );
        P2 ( _XY_,_XY_ ) = P ( 5,5 );


        return P2;
    }
    
    
    

};


#endif //TPZYCMOHRCOULOMBPV_H
