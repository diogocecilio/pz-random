/*
 *  TPZMohrCoulomb.h
 *  FEMPZ
 *
 *  Created by Nathan Shauer on 5/4/13.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef TPZMohrCoulombVoigt_H
#define TPZMohrCoulombVoigt_H

#include "pzlog.h"
#include "TPZTensor.h"
#include "pzvec_extras.h"
#include "pzsave.h"
#include "TPZPlasticState.h"
#include "TPZElasticResponse.h"


class TPZMohrCoulombVoigt
{

public:
    enum {NYield=3};

private:
    REAL fPhi;
    REAL fPsi;
    REAL fc;
    TPZElasticResponse fER;


protected:
    REAL fEpsPlasticBar;

public:

    /// structure which contains the decision tree of the return map
    // we can only expect a consistent tangent matrix if the decision tree remains the same


public:

    /**
     * @brief empty constructor
     */
    TPZMohrCoulombVoigt();

    /**
     * @brief Constructor seting yc parameters
     */
    TPZMohrCoulombVoigt ( REAL Phi, REAL Psi, REAL c, TPZElasticResponse &ER );

    /**
     * @brief Copy Constructor
     */
    TPZMohrCoulombVoigt ( const TPZMohrCoulombVoigt &cp );

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
        if ( fc<1.e-3 ) DebugStop();
        fc =   state.fmatprop[0];
        fPhi = state.fmatprop[1];
        fPsi = state.fmatprop[1];
	//	std::cout << "fc = "<< fc <<endl;
    }

    /**
     * @brief Operator =
     */
    TPZMohrCoulombVoigt & operator= ( const TPZMohrCoulombVoigt &cp );


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

    REAL ReturnMapPlane (  TPZTensor<REAL> &sigma_trial, TPZTensor<REAL> &sigma_proj,TPZFMatrix<REAL>&dep, REAL &epsbarnew ) ;


    bool ReturnMapLeftEdge (TPZTensor<REAL> &sigma_trial, TPZTensor<REAL> &sigma_proj,TPZFMatrix<REAL>&dep, REAL &epsbarnew) ;


    void ComputeLeftEdgeTangent ( TPZTensor<REAL> &sigma_trial, TPZTensor<REAL> &sigma_proj,TPZFMatrix<REAL>&dep, REAL &epsbarnew) ;


    bool ReturnMapRightEdge ( TPZTensor<REAL> &sigma_trial, TPZTensor<REAL> &sigma_proj,TPZFMatrix<REAL>&dep, REAL &epsbarnew ) ;


    bool ReturnMapApex ( TPZTensor<REAL> &sigma_trial, TPZTensor<REAL> &sigma_proj,TPZFMatrix<REAL>&dep, REAL &epsbarnew ) ;


    void ProjectSigmaDep ( TPZTensor<REAL> &sigma_trial, TPZTensor<REAL> &sigma_proj,TPZFMatrix<REAL>&dep, REAL &epsbarnew);

    void ProjectSigmaDep2( TPZTensor<REAL> &sigma_trial, TPZTensor<REAL> &sigma_proj,TPZFMatrix<REAL>&dep, REAL &epsbarnew );

    bool CheckOrder(TPZManVector<REAL,3>  eigenvalues);

    TPZManVector<REAL,3> EigenVal(TPZTensor<REAL> & tensor);

    inline void FlowVector(TPZTensor<REAL> &sigma,REAL & a, REAL & dadt, REAL & d2Adt,TPZTensor<REAL> &flowvec)
    {

        TPZTensor<REAL> di1 = sigma.dI1();
        TPZTensor<REAL> dj2 = sigma.dJ2();
        TPZTensor<REAL> dj3 = sigma.dJ3();

        REAL c1=C1();
        REAL c2=C2 ( sigma, a, dadt);
        REAL c3 = C3 ( sigma,dadt );
        di1*=c1;
        dj2*=c2;
        dj3*=c3;


        flowvec=di1;
        flowvec+=dj2;
        flowvec+=dj3;
    }

    inline void Yield(TPZTensor<REAL> &sigma,REAL &a,REAL &f)
    {
         REAL I1 = sigma.I1();
         REAL J2 = sigma.J2();
         f = -(fc*cos(fPhi)) + a*sqrt(J2) + 0.3333333333333333*I1*sin(fPhi);
    }


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



    REAL theta(TPZTensor<REAL> tensor)
    {
        REAL J2 =tensor.J2();

        if( J2 < 0.000000001)
        {
           J2=0.000000001;
        }

        REAL denom=pow(J2,1.5);
        REAL val = (-2.598076211353316*tensor.J3())/denom;

        REAL valone=0.9999999999999999999999999999999999;

        if(val>valone )
        {
            return 1./3.*asin(valone);

        }

        if(val<-valone){

            return 1./3.*asin(-valone);
        }

        return 1./3.*asin(val);

    }

    TPZTensor<REAL> FlowVectorMain(TPZTensor<REAL> & tensor)
    {
         TPZTensor<REAL> flowvec;
         REAL dadt = dAdt(tensor);
         REAL d2adt = d2Adt(tensor);
         REAL a=A(tensor);
         FlowVector(tensor,a,dadt,d2adt,flowvec);
         return flowvec;
    }

    TPZTensor<REAL> FlowVectorLeft(TPZTensor<REAL> & tensor)
    {
         TPZTensor<REAL> flowvec;
         REAL dadt = dA3dt(tensor);
         REAL d2adt = d2A3dt(tensor);
         REAL a=A3(tensor);
         FlowVector(tensor,a,dadt,d2adt,flowvec);
         return flowvec;
    }

    TPZTensor<REAL> FlowVectorRight(TPZTensor<REAL> & tensor)
    {
         TPZTensor<REAL> flowvec;
         REAL dadt = dA2dt(tensor);
         REAL d2adt = d2A2dt(tensor);
         REAL a=A2(tensor);
         FlowVector(tensor,a,dadt,d2adt,flowvec);
         return flowvec;
    }

    REAL FMain(TPZTensor<REAL> & tensor)
    {
         REAL I1 = tensor.I1();
         REAL J2 = tensor.J2();
         REAL a =A(tensor);
         return (-(fc*cos(fPhi)) + a*sqrt(J2) + 0.3333333333333333*I1*sin(fPhi));
    }
    REAL FLeft(TPZTensor<REAL> & tensor)
    {
         REAL I1 = tensor.I1();
         REAL J2 = tensor.J2();
         REAL a =A3(tensor);
         return (-(fc*cos(fPhi)) + a*sqrt(J2) + 0.3333333333333333*I1*sin(fPhi));
    }
    REAL FRight(TPZTensor<REAL> & tensor)
    {
         REAL I1 = tensor.I1();
         REAL J2 = tensor.J2();
         REAL a =A2(tensor);
         return (-(fc*cos(fPhi)) + a*sqrt(J2) + 0.3333333333333333*I1*sin(fPhi));
    }


    REAL A(TPZTensor<REAL> tensor)
    {

        REAL t= theta(tensor);
        return cos(t) - 0.5773502691896258*sin(fPhi)*sin(t);
    }

    REAL dAdt(TPZTensor<REAL> tensor)
    {

        REAL t= theta(tensor);
        return -0.5773502691896258*cos(t)*sin(fPhi) - sin(t);
    }

    REAL d2Adt(TPZTensor<REAL> tensor)
    {

        REAL t= theta(tensor);
        return -cos(t) + 0.5773502691896258*sin(fPhi)*sin(t);
    }

    REAL A2(TPZTensor<REAL> tensor)
    {

        REAL t= theta(tensor);
        return 0.5*cos(t)*(1 + sin(fPhi)) + 0.2886751345948129*(-3. + sin(fPhi))*sin(t);
    }

    REAL dA2dt(TPZTensor<REAL> tensor)
    {
        REAL t= theta(tensor);
        return 0.2886751345948129*cos(t)*(-3. + sin(fPhi)) - 0.5*(1. + sin(fPhi))*sin(t);

    }

    REAL d2A2dt(TPZTensor<REAL> tensor)
    {
        REAL t= theta(tensor);
        //-0.5*Cos(t)*(1. + Sin(phi)) - 0.2886751345948129*(-3. + Sin(phi))*Sin(t)
        return -0.5*cos(t)*(1. + sin(fPhi)) - 0.2886751345948129*(-3. + sin(fPhi))*sin(t);
    }


    REAL A3(TPZTensor<REAL> tensor)
    {
        REAL t= theta(tensor);
        return 0.5*cos(t)*(1. - sin(fPhi)) + 0.2886751345948129*(3. + sin(fPhi))*sin(t);
    }

    REAL dA3dt(TPZTensor<REAL> tensor)
    {
        REAL t= theta(tensor);
        return 0.2886751345948129*cos(t)*(3. + sin(fPhi)) - 0.5*(1. - sin(fPhi))*sin(t);
    }

    REAL d2A3dt(TPZTensor<REAL> tensor)
    {
        REAL t= theta(tensor);
        return -0.5*cos(t)*(1. - sin(fPhi)) - 0.2886751345948129*(3. + sin(fPhi))*sin(t);
    }


    REAL C1()
    {
        return sin ( fPhi ) /3.;
    }
//     REAL C2 ( TPZTensor<REAL> tensor, REAL (*funcA)(TPZTensor<REAL>),REAL (*dadt)(TPZTensor<REAL>)  )
//     {
//         REAL J2=tensor.J2();
//         REAL a=funcA(tensor);
//         return (0.5*(a - dadt(tensor)*tan(3.*theta(tensor))))/pow(J2,0.5);
//
//     }

    REAL C2 ( TPZTensor<REAL> tensor, REAL funcA,REAL dadt  )
    {
        REAL thetaval = theta(tensor);
        REAL J2=tensor.J2();
        REAL a=funcA;
        REAL powj2=pow(J2,0.5);
        REAL tantheta=tan(3.*thetaval);
        REAL num = (a - dadt*tantheta);
        REAL c2=(0.5*num)/powj2;
        return c2;

    }


    REAL C3 ( TPZTensor<REAL> tensor,REAL dadt)
    {
        REAL thetaval = theta(tensor);
        REAL J2=tensor.J2();
        REAL da=dadt;//dA(tensor,thetaval);
        return -sqrt ( 3. ) *da/ ( 2.*J2*cos ( 3*thetaval ) );
    }



    REAL C4 ( TPZTensor<REAL> tensor, REAL dadt,REAL da2dt )
    {
        REAL thetaval = theta(tensor);
        REAL da=dadt;
        REAL d2a=da2dt;
        return d2a+3.*tan ( 3.*thetaval ) *da;

    }



    REAL C23 ( TPZTensor<REAL> tensor, REAL funcA,REAL dadt,REAL da2dt )
    {
        REAL thetaval = theta(tensor);
        REAL J2=tensor.J2();
        REAL da=dadt;//dA(tensor,thetaval);
        REAL a=funcA;//A(tensor,thetaval);

        REAL c4 = C4 ( tensor ,dadt, da2dt);
        return ( 1./2.*tan ( 3.*thetaval ) *c4+da ) *sqrt ( 3. ) / ( 2*J2*J2*cos ( 3.*thetaval ) );
    }



    REAL C22 ( TPZTensor<REAL> tensor, REAL funcA,REAL dadt,REAL da2dt)
    {
        REAL thetaval = theta(tensor);
        REAL J2=tensor.J2();
        REAL da=dadt;//dA(tensor,thetaval);
        REAL a=funcA;//A(tensor,thetaval);
        REAL c4 = C4 ( tensor ,dadt, da2dt);

        return - ( a-pow ( tan ( 3*thetaval ),2 ) *c4-3.*tan ( 3.*thetaval ) *da ) / ( 4.*pow ( J2,3./2. ) );
    }



    REAL C33 ( TPZTensor<REAL> tensor, REAL funcA,REAL dadt,REAL da2dt )
    {
        REAL J2=tensor.J2();
        REAL thetaval = theta(tensor);
        REAL c4 = C4 ( tensor ,dadt, da2dt);
        return 3.*c4/ ( 4.*pow ( J2,5/2. ) *cos ( 3.*thetaval ) *cos ( 3.*thetaval ) );
    }


    void ComputeConsistentPlaneTangent ( TPZTensor<REAL> &trialstress,TPZFMatrix<REAL> & Dep );
    void ComputeConsistentEdgeTangent ( TPZTensor<REAL> &trialstress,TPZFMatrix<REAL> & Dep, int &m_type );


 TPZFMatrix<REAL>  dAdsig ( TPZTensor<REAL> tensor, REAL funcA,REAL dadt,REAL da2dt )
    {


        REAL c1=C1();
        REAL c2=C2 ( tensor,funcA,dadt );
        REAL c3 = C3 ( tensor,dadt );
        REAL c4=C4 (  tensor, dadt,da2dt);
        REAL c23=C23 ( tensor,funcA,dadt,da2dt );
        REAL c22=C22 ( tensor,funcA,dadt,da2dt);
        REAL c33=C33 ( tensor,funcA,dadt,da2dt );
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

// #define _XX_ 0
// #define _XY_ 1
// #define _XZ_ 2
// #define _YY_ 3
// #define _YZ_ 4
// #define _ZZ_ 5

};


#endif //TPZYCMOHRCOULOMBPV_H
