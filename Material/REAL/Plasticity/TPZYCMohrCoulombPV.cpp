#include "TPZYCMohrCoulombPV.h"

typedef TFad<3,REAL> fadtype;

TPZYCMohrCoulombPV::TPZYCMohrCoulombPV() : fPhi ( 0. ),fPsi ( 0. ), fc ( 0. ), fEpsPlasticBar ( 0. ), fER()
{

}

TPZYCMohrCoulombPV::TPZYCMohrCoulombPV ( REAL Phi, REAL Psi, REAL c, TPZElasticResponse &ER ) : fPhi ( Phi ), fPsi ( Psi ), fc ( c ), fEpsPlasticBar ( 0. ), fER ( ER )
{

}

TPZYCMohrCoulombPV::TPZYCMohrCoulombPV ( const TPZYCMohrCoulombPV &cp ) // : 	fPhi(cp.fPhi), fPsi(cp.fPsi), fc(cp.fc), fEpsPlasticBar(cp.fEpsPlasticBar), fER(cp.fER)
{
    TPZYCMohrCoulombPV::operator= ( cp );
}

TPZYCMohrCoulombPV & TPZYCMohrCoulombPV::operator= ( const TPZYCMohrCoulombPV &cp )
{
    fPhi = cp.fPhi;
    fPsi = cp.fPsi;
    fc = cp.fc;
    fEpsPlasticBar = cp.fEpsPlasticBar;
    fER = cp.fER;
    return *this;
}


void TPZYCMohrCoulombPV::Read ( TPZStream &buf )
{
    buf.Read ( &fPhi );
    buf.Read ( &fPsi );
    buf.Read ( &fc );
    buf.Read ( &fEpsPlasticBar );
    fER.Read ( buf );

}

void TPZYCMohrCoulombPV::Write ( TPZStream &buf ) const
{
    buf.Write ( &fPhi );
    buf.Write ( &fPsi );
    buf.Write ( &fc );
    buf.Write ( &fEpsPlasticBar );
    fER.Write ( buf );
}


void TPZYCMohrCoulombPV::Phi ( TPZVec<STATE> sigvec,STATE alpha,TPZVec<STATE> &phi ) const
{
    phi.resize ( 3 );
    for ( int i=0; i<3; i++ ) phi[i]=0;
    phi[0]= PhiPlane ( sigvec );
}

template <class T>
void TPZYCMohrCoulombPV::PlasticityFunction ( const T epsp, T &c, T &H ) const
{
    c = fc; // c(epsp)
    H = 0.; // dc(epsp)/depsp
}

template<class T>
TPZVec<T> TPZYCMohrCoulombPV::SigmaElastPV ( const TPZVec<T> &deform ) const
{
    T trace = deform[0]+deform[1]+deform[2];
    TPZVec<T> sigma ( 3,0. );
    sigma = trace * fER.Lambda() + 2 * fER.G() * deform[0];
    sigma = trace * fER.Lambda() + 2 * fER.G() * deform[1];
    sigma = trace * fER.Lambda() + 2 * fER.G() * deform[2];

    return sigma;
}

template<class T>
T TPZYCMohrCoulombPV::PhiPlane ( const TPZVec<T> &sigma ) const
{
    const REAL sinphi = sin ( fPhi );
    const REAL cosphi = cos ( fPhi );
    T c,H;
    PlasticityFunction ( T ( fEpsPlasticBar ),c, H );

    return sigma[0]-sigma[2]+ ( sigma[0]+sigma[2] ) *sinphi-2.*c*cosphi;
	//return (sigma[0]-sigma[2])/2+ ( sigma[0]+sigma[2] )/2 *sinphi-c*cosphi;
}

template<class T>
bool TPZYCMohrCoulombPV::ReturnMapPlane ( const TPZVec<T> &sigma_trial, TPZVec<T> &sigma_projected,
        TComputeSequence &memory, REAL &epsbarnew ) const
{
    sigma_projected = sigma_trial;
		//cout <<sigma_trial[0]<<endl;
	//cout <<sigma_trial[1]<<endl;
	//cout <<sigma_trial[2]<<endl;
    TPZManVector<T,3> eigenvalues = sigma_projected;
    const REAL sinphi = sin ( fPhi );
    const REAL sinpsi = sin ( fPsi );
    const REAL cosphi = cos ( fPhi );
    const REAL sinphi2 = sinphi*sinphi;
    const REAL cosphi2 = 1.-sinphi2;
    const REAL constA = 4.* fER.G() * ( 1.+ sinphi*sinpsi/3. ) + 4.*fER.K() * sinphi*sinpsi;
    T sigmay,H;
    T epsbar = T ( fEpsPlasticBar+memory.fGamma[0]*2.*cosphi ); //diogo aqui
    //cout << " fPhi = " << fPhi << endl;
    //cout << " cos = " << fc << endl;
    PlasticityFunction ( epsbar,sigmay, H );
    T phi = eigenvalues[0]-eigenvalues[2]+ ( eigenvalues[0]+eigenvalues[2] ) *sinphi-2.*sigmay*cosphi;
    T gamma = memory.fGamma[0];
    REAL phival = shapeFAD::val ( phi );
    REAL tolerance = 1.e-8;
    do {
        T denom = -constA- T ( 4.*cosphi2 ) *H;
        //            T d = T(-4.*G()*(1.+sinphi*sinpsi/3.)-4.*K()*sinphi*sinpsi)-T(4.*cosphi2)*H;
        T deriv_gamma = -phi/denom;
        gamma += deriv_gamma;
        epsbar = T ( fEpsPlasticBar )+gamma*T ( 2.*cosphi ); //errado esta inicializando toda vez. diogo aqui
        PlasticityFunction ( epsbar, sigmay, H );
//		if (shapeFAD::val(H) < 0.) {
//			DebugStop();
//		}

        phi = eigenvalues[0]-eigenvalues[2]+ ( eigenvalues[0]+eigenvalues[2] ) *sinphi-2.*sigmay*cosphi - constA*gamma;
        phival = shapeFAD::val ( phi );

    } while ( abs ( phival ) > tolerance );

    memory.fGamma[0] = shapeFAD::val ( gamma );
    eigenvalues[0] -= T ( 2.*fER.G() * ( 1+sinpsi/3. )+2.*fER.K() *sinpsi ) *gamma;
    eigenvalues[1] += T ( ( 4.*fER.G() /3. - fER.K() *2. ) *sinpsi ) *gamma;
    eigenvalues[2] += T ( 2.*fER.G() * ( 1-sinpsi/3. )-2.*fER.K() *sinpsi ) *gamma;
	//cout <<" \n gamma" <<endl;
	//cout <<gamma<<endl;
    sigma_projected = eigenvalues;
	//cout <<eigenvalues[0]<<endl;
	//cout <<eigenvalues[1]<<endl;
	//cout <<eigenvalues[2]<<endl;
    epsbarnew = shapeFAD::val ( epsbar );

#ifdef DEBUG
    phi = eigenvalues[0]-eigenvalues[2]+ ( eigenvalues[0]+eigenvalues[2] ) *sinphi-2.*sigmay*cosphi;
#endif
    return ( shapeFAD::val ( eigenvalues[0] ) >shapeFAD::val ( eigenvalues[1] ) && shapeFAD::val ( eigenvalues[1] ) > shapeFAD::val ( eigenvalues[2] ) );
}

void TPZYCMohrCoulombPV::ComputePlaneTangent ( TPZMatrix<REAL> &tang, REAL &epsbarp ) const
{
    const REAL sinphi = sin ( fPhi );
    const REAL sinpsi = sin ( fPsi );
    const REAL cosphi = cos ( fPhi );
    const REAL cosphi2 = cosphi*cosphi;
    const REAL G = fER.G(), K = fER.K();
    const REAL c1 = 2.*G* ( 1.+1./3.*sinpsi ) + 2.*K*sinpsi;
    const REAL c2 = ( 4.*G/3.-2.*K ) *sinpsi;
    const REAL c3 = 2.*G* ( 1.-1./3.*sinpsi ) - 2.*K*sinpsi;
    const REAL constA = 4.* G * ( 1.+ sinphi*sinpsi/3. ) + 4.* K * sinphi*sinpsi;
    REAL epsbar = epsbarp;
    REAL c, H;
    PlasticityFunction ( epsbar, c, H );
    const REAL denom = constA + 4 * cosphi2*H;
    const REAL dGds1 = ( 1+sinphi ) /denom; // Derivate of gamma with respect to Sigma1tr
    const REAL dGds2 = 0.; // Only created to remember that PHIfunc doesnt depend in it
    const REAL dGds3 = ( -1+sinphi ) /denom;
    tang.Redim ( 3, 3 );

    // First column
    tang ( 0,0 ) = 1.-c1*dGds1;
    tang ( 1,0 ) = c2*dGds1;
    tang ( 2,0 ) = c3*dGds1;

    // Second column
    tang ( 1,1 ) = 1.; // The others are 0

    // Third column
    tang ( 0,2 ) = -c1*dGds3;
    tang ( 1,2 ) = c2*dGds3;
    tang ( 2,2 ) = 1.+c3*dGds3;
}

template<class T>
bool TPZYCMohrCoulombPV::ReturnMapLeftEdge ( const TPZVec<T> &sigma_trial, TPZVec<T> &sigma_projected,
        TComputeSequence &memory, REAL &epsbarnew ) const
{

    sigma_projected = sigma_trial;
    TPZManVector<T,3> eigenvalues = sigma_projected;
    const REAL sinphi = sin ( fPhi );
    const REAL sinpsi = sin ( fPsi );
    const REAL cosphi = cos ( fPhi );
    const REAL sinphi2 = sinphi*sinphi;
    const REAL cosphi2 = 1.-sinphi2;
    TPZManVector<T,2> gamma ( 2,0. ),phi ( 2,0. ),sigma_bar ( 2,0. ),ab ( 2,0. );
    gamma[0] = memory.fGamma[0];
    gamma[1] = memory.fGamma[1];
    TPZManVector<REAL,2> phival ( 2,0. );
    TPZManVector<TPZManVector<T,2>,2> d ( 2 ),dinverse ( 2 );
    for ( int i = 0; i < 2; i++ ) {
        d[i].Resize ( 2,0. );
        dinverse[i].Resize ( 2,0. );
    }
    sigma_bar[0] = eigenvalues[0]-eigenvalues[2]+ ( eigenvalues[0]+eigenvalues[2] ) *T ( sinphi );
    sigma_bar[1] = eigenvalues[1]-eigenvalues[2]+ ( eigenvalues[1]+eigenvalues[2] ) *T ( sinphi );
    T sigmay,H;
    T epsbar = T ( fEpsPlasticBar ) + ( gamma[0]+gamma[1] ) *T ( 2.*cosphi ); //diogo aqui
    PlasticityFunction ( epsbar,sigmay, H );
    phi[0] = sigma_bar[0] - T ( 2.*cosphi ) *sigmay;
    phi[1] = sigma_bar[1] - T ( 2.*cosphi ) *sigmay;
    ab[0] = T ( 4.*fER.G() * ( 1+sinphi*sinpsi/3. )+4.*fER.K() *sinphi*sinpsi );
    ab[1] = T ( 2.*fER.G() * ( 1.-sinphi-sinpsi-sinphi*sinpsi/3. )+4.*fER.K() *sinphi*sinpsi );
    T residual =1;
    REAL tolerance = 1.e-8;
    do {
        d[0][0] = -ab[0]-T ( 4.*cosphi2 ) *H;
        d[1][0] = -ab[1]-T ( 4.*cosphi2 ) *H;
        d[0][1] = -ab[1]-T ( 4.*cosphi2 ) *H;
        d[1][1] = -ab[0]-T ( 4.*cosphi2 ) *H;
        T detd = d[0][0]*d[1][1]-d[0][1]*d[1][0];
        dinverse[0][0] = d[1][1]/detd;
        dinverse[1][0] = -d[1][0]/detd;
        dinverse[0][1]= -d[0][1]/detd;
        dinverse[1][1] = d[0][0]/detd;
        gamma[0] -= ( dinverse[0][0]*phi[0]+dinverse[0][1]*phi[1] );
        gamma[1] -= ( dinverse[1][0]*phi[0]+dinverse[1][1]*phi[1] );
        epsbar = T ( fEpsPlasticBar )+ ( gamma[0]+gamma[1] ) *T ( 2.*cosphi );
        PlasticityFunction ( epsbar, sigmay, H );
        phi[0] = sigma_bar[0] - ab[0]*gamma[0] - ab[1]*gamma[1] - T ( 2.*cosphi ) *sigmay;
        phi[1] = sigma_bar[1] - ab[1]*gamma[0] - ab[0]*gamma[1] - T ( 2.*cosphi ) *sigmay;
        phival[0] = shapeFAD::val ( phi[0] );
        phival[1] = shapeFAD::val ( phi[1] );
        residual= ( fabs ( phival[0] )+fabs ( phival[1] ) );
    } while ( residual>tolerance );

    memory.fGamma[0] = shapeFAD::val ( gamma[0] );
    memory.fGamma[1] = shapeFAD::val ( gamma[1] );
    eigenvalues[0] += - T ( 2.*fER.G() * ( 1+sinpsi/3. )+2.*fER.K() *sinpsi ) *gamma[0]+T ( ( 4.*fER.G() /3.-2.*fER.K() ) *sinpsi ) *gamma[1];
    eigenvalues[1] += T ( ( 4.*fER.G() /3.- fER.K() *2. ) *sinpsi ) *gamma[0]-T ( 2.*fER.G() * ( 1.+sinpsi/3. )+2.*fER.K() *sinpsi ) *gamma[1];
    eigenvalues[2] += T ( 2.*fER.G() * ( 1-sinpsi/3. )-2.*fER.K() *sinpsi ) * ( gamma[0]+gamma[1] );
    sigma_projected = eigenvalues;
    epsbarnew	= shapeFAD::val ( epsbar );

    return ( shapeFAD::val ( eigenvalues[0] ) >=shapeFAD::val ( eigenvalues[1] ) && shapeFAD::val ( eigenvalues[1] ) >= shapeFAD::val ( eigenvalues[2] ) );
}

/**
 * @brief Computes dsigmapr/dsigmatr for the ReturnMapLeftEdge
 */
void TPZYCMohrCoulombPV::ComputeLeftEdgeTangent ( TPZMatrix<REAL> &tang, REAL &epsbarp ) const
{
    const REAL sinphi = sin ( fPhi );
    const REAL sinpsi = sin ( fPsi );
    const REAL cosphi = cos ( fPhi );
    const REAL cosphi2 = cosphi*cosphi;
    const REAL G = fER.G(), K = fER.K();
    const REAL c1 = 2.*G* ( 1.+1./3.*sinpsi ) + 2.*K*sinpsi;
    const REAL c2 = ( 4.*G/3.-2.*K ) *sinpsi;
    const REAL c3 = 2.*G* ( 1.-1./3.*sinpsi ) - 2.*K*sinpsi;
    const REAL constA = 4.* G * ( 1.+ sinphi*sinpsi/3. ) + 4.* K * sinphi*sinpsi;
    const REAL constB = 2.*G* ( 1-sinphi-sinpsi-1./3.*sinphi*sinpsi ) + 4.*K*sinphi*sinpsi;
    REAL epsbar = epsbarp;
    REAL c, H;
    PlasticityFunction ( epsbar, c, H );
    const REAL cos2H4 = 4.*cosphi2*H;
    const REAL denom = ( constA-constB ) * ( constA+constB+8.*cosphi2*H );
    const REAL dGads1 = ( cos2H4* ( 1.+sinphi )+constA* ( 1.+sinphi ) ) /denom; // Derivative of DgammaA with respect to Sigma1tr
    const REAL dGads2 = ( -cos2H4* ( 1.+sinphi )-constB* ( 1.+sinphi ) ) /denom; // Derivative of DgammaA with respect to Sigma2tr
    const REAL dGads3 = ( constA* ( -1.+sinphi )-constB* ( -1.+sinphi ) ) /denom; // Derivative of DgammaA with respect to Sigma3tr
    const REAL dGbds1 = ( -cos2H4* ( 1.+sinphi )-constB* ( 1.+sinphi ) ) /denom; // Derivative of DgammaB with respect to Sigma1tr
    const REAL dGbds2 = ( cos2H4* ( 1.+sinphi )+constA* ( 1.+sinphi ) ) /denom; // Derivative of DgammaA with respect to Sigma2tr
    const REAL dGbds3 = ( constA* ( -1.+sinphi )-constB* ( -1.+sinphi ) ) /denom; // Derivative of DgammaA with respect to Sigma3tr

    tang.Redim ( 3, 3 );

    // First column
    tang ( 0,0 ) = 1.-c1*dGads1+c2*dGbds1;
    tang ( 1,0 ) = c2*dGads1-c1*dGbds1;
    tang ( 2,0 ) = c3* ( dGads1+dGbds1 );

    // Second column
    tang ( 0,1 ) = -c1*dGads2+c2*dGbds2;
    tang ( 1,1 ) = 1.+c2*dGads2-c1*dGbds2;
    tang ( 2,1 ) = c3* ( dGads2+dGbds2 );

    // Third column
    tang ( 0,2 ) = -c1*dGads3+c2*dGbds3;
    tang ( 1,2 ) = c2*dGads3-c1*dGbds3;
    tang ( 2,2 ) = 1.+c3* ( dGads3+dGbds3 );
}


/**
 * @brief Implements the return map in the right edge of the surface
 */
template<class T>
bool TPZYCMohrCoulombPV::ReturnMapRightEdge ( const TPZVec<T> &sigma_trial, TPZVec<T> &sigma_projected,
        TComputeSequence &memory, REAL &epsbarnew ) const
{
    sigma_projected = sigma_trial;
    TPZManVector<T,3> eigenvalues = sigma_projected;
    const REAL sinphi = sin ( fPhi );
    const REAL sinpsi = sin ( fPsi );
    const REAL cosphi = cos ( fPhi );
    const REAL sinphi2 = sinphi*sinphi;
    const REAL cosphi2 = 1.-sinphi2;
    const REAL KV = fER.K();
    const REAL GV = fER.G();
    TPZManVector<T,2> gamma ( 2,0. ),phi ( 2,0. ),sigma_bar ( 2,0. ),ab ( 2,0. );
    gamma[0] = memory.fGamma[0];
    gamma[1] = memory.fGamma[1];
    TPZManVector<REAL,2> phival ( 2,0. );
    TPZManVector<TPZManVector<T,2>,2> d ( 2 ),dinverse ( 2 );
    for ( int i = 0; i < 2; i++ ) {
        d[i].Resize ( 2,0. );
        dinverse[i].Resize ( 2,0. );
    }
    sigma_bar[0] = eigenvalues[0]-eigenvalues[2]+ ( eigenvalues[0]+eigenvalues[2] ) *T ( sinphi );
    sigma_bar[1] = eigenvalues[0]-eigenvalues[1]+ ( eigenvalues[0]+eigenvalues[1] ) *T ( sinphi );
    T sigmay,H;
    T epsbar = T ( fEpsPlasticBar )+ ( gamma[0]+gamma[1] ) *T ( 2.*cosphi );
    PlasticityFunction ( epsbar,sigmay, H );
    phi[0] = sigma_bar[0] - T ( 2.*cosphi ) *sigmay;
    phi[1] = sigma_bar[1] - T ( 2.*cosphi ) *sigmay;
    ab[0] = T ( 4.*GV* ( 1+sinphi*sinpsi/3. )+4.*KV*sinphi*sinpsi );
    ab[1] = T ( 2.*GV* ( 1.+sinphi+sinpsi-sinphi*sinpsi/3. )+4.*KV*sinphi*sinpsi );

#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout << "phi = " << phi << std::endl;
        LOGPZ_DEBUG ( loggerMohrCoulombPV, sout.str() )
    }
#endif

    REAL tolerance = 1.e-8;
    int iter = 0;
    T residual =1;
    do {
#ifdef LOG4CXX
        {
            std::stringstream sout;
            sout << "epsbar = " << epsbar << std::endl;
            sout << "sigmay = " << sigmay << std::endl;
            sout << "H = " << H << std::endl;
            LOGPZ_DEBUG ( loggerMohrCoulombPV, sout.str() )
        }
#endif
        d[0][0] = -ab[0]-T ( 4.*cosphi2 ) *H;
        d[1][0] = -ab[1]-T ( 4.*cosphi2 ) *H;
        d[0][1] = -ab[1]-T ( 4.*cosphi2 ) *H;
        d[1][1] = -ab[0]-T ( 4.*cosphi2 ) *H;
        T detd = d[0][0]*d[1][1]-d[0][1]*d[1][0];
        dinverse[0][0] = d[1][1]/detd;
        dinverse[1][0] = -d[1][0]/detd;
        dinverse[0][1]= -d[0][1]/detd;
        dinverse[1][1] = d[0][0]/detd;
        gamma[0] -= ( dinverse[0][0]*phi[0]+dinverse[0][1]*phi[1] );
        gamma[1] -= ( dinverse[1][0]*phi[0]+dinverse[1][1]*phi[1] );
        epsbar = T ( fEpsPlasticBar )+ ( gamma[0]+gamma[1] ) *T ( 2.*cosphi );
        PlasticityFunction ( epsbar, sigmay, H );
//		if (shapeFAD::val(H) < 0.) {
//			DebugStop();
//		}
        iter++;
        phi[0] = sigma_bar[0] - ab[0]*gamma[0] - ab[1]*gamma[1] - T ( 2.*cosphi ) *sigmay;
        phi[1] = sigma_bar[1] - ab[1]*gamma[0] - ab[0]*gamma[1] - T ( 2.*cosphi ) *sigmay;
        phival[0] = shapeFAD::val ( phi[0] );
        phival[1] = shapeFAD::val ( phi[1] );
#ifdef LOG4CXX
        {
            std::stringstream sout;
            sout << "iter = " << iter << " phi = " << phival << std::endl;
            LOGPZ_DEBUG ( loggerMohrCoulombPV, sout.str() )
        }
#endif
        residual= ( fabs ( phival[0] )+fabs ( phival[1] ) );
    } while ( residual>tolerance );

    memory.fGamma[0] = shapeFAD::val ( gamma[0] );
    memory.fGamma[1] = shapeFAD::val ( gamma[1] );

    eigenvalues[0] -= T ( 2.*GV* ( 1+sinpsi/3. )+2.*KV*sinpsi ) * ( gamma[0]+gamma[1] );
    eigenvalues[1] += T ( ( 4.*GV/3.- KV*2. ) *sinpsi ) *gamma[0]+T ( 2.*GV* ( 1.-sinpsi/3. )-2.*KV*sinpsi ) *gamma[1];
    eigenvalues[2] += T ( 2.*GV* ( 1-sinpsi/3. )-2.*KV*sinpsi ) *gamma[0]+T ( ( 4.*GV/3.-2.*KV ) *sinpsi ) *gamma[1];
    sigma_projected = eigenvalues;
    epsbarnew = shapeFAD::val ( epsbar );

    return ( shapeFAD::val ( eigenvalues[0] ) >=shapeFAD::val ( eigenvalues[1] ) && shapeFAD::val ( eigenvalues[1] ) >= shapeFAD::val ( eigenvalues[2] ) );
}

/**
 * @brief Computes dsigmapr/dsigmatr for the ReturnMapRightEdge
 */
void TPZYCMohrCoulombPV::ComputeRightEdgeTangent ( TPZMatrix<REAL> &tang, REAL &epsbarp ) const
{

    const REAL sinphi = sin ( fPhi );
    const REAL sinpsi = sin ( fPsi );
    const REAL cosphi = cos ( fPhi );
    const REAL cosphi2 = cosphi*cosphi;
    const REAL G = fER.G(), K = fER.K();
    const REAL c1 = 2.*G* ( 1.+1./3.*sinpsi ) + 2.*K*sinpsi;
    const REAL c2 = ( 4.*G/3.-2.*K ) *sinpsi;
    const REAL c3 = 2.*G* ( 1.-1./3.*sinpsi ) - 2.*K*sinpsi;
    const REAL constA = 4.* G * ( 1.+ sinphi*sinpsi/3. ) + 4.* K * sinphi*sinpsi;
    const REAL constB = 2.*G* ( 1+sinphi+sinpsi-1./3.*sinphi*sinpsi ) + 4.*K*sinphi*sinpsi;
    REAL epsbar = epsbarp;
    REAL c, H;
    PlasticityFunction ( epsbar, c, H );
    const REAL cos2H4 = 4.*cosphi2*H;
    const REAL denom = ( constA-constB ) * ( constA+constB+8.*cosphi2*H );
    const REAL dGads1 = ( -constB* ( 1+sinphi )+constA* ( 1+sinphi ) ) /denom; // Derivative of DgammaA with respect to Sigma1tr
    const REAL dGads2 = ( -cos2H4* ( -1.+sinphi )-constB* ( -1.+sinphi ) ) /denom; // Derivative of DgammaA with respect to Sigma2tr
    const REAL dGads3 = ( cos2H4* ( -1.+sinphi )+constA* ( -1.+sinphi ) ) /denom; // Derivative of DgammaA with respect to Sigma3tr
    const REAL dGbds1 = ( constA* ( 1.+sinphi )-constB* ( 1.+sinphi ) ) /denom; // Derivative of DgammaB with respect to Sigma1tr
    const REAL dGbds2 = ( cos2H4* ( -1.+sinphi )+constA* ( -1.+sinphi ) ) /denom; // Derivative of DgammaA with respect to Sigma2tr
    const REAL dGbds3 = ( -cos2H4* ( -1.+sinphi )-constB* ( -1.+sinphi ) ) /denom; // Derivative of DgammaA with respect to Sigma3tr

    tang.Redim ( 3, 3 );

    // First column
    tang ( 0,0 ) = 1.-c1* ( dGads1+dGbds1 );
    tang ( 1,0 ) = c2*dGads1+c3*dGbds1;
    tang ( 2,0 ) = c3*dGads1+c2*dGbds1;

    // Second column
    tang ( 0,1 ) = -c1* ( dGads2+dGbds2 );
    tang ( 1,1 ) = 1.+c2*dGads2+c3*dGbds2;
    tang ( 2,1 ) = c3*dGads2+c2*dGbds2;

    // Third column
    tang ( 0,2 ) = -c1* ( dGads3+dGbds3 );
    tang ( 1,2 ) = c2*dGads3+c3*dGbds3;
    tang ( 2,2 ) = 1.+c3*dGads3+c2*dGbds3;
}

template<class T>
bool TPZYCMohrCoulombPV::ReturnMapApex ( const TPZVec<T> &sigmatrial, TPZVec<T> &sigma_projected,
        TComputeSequence &memory, REAL &epsbarnew ) const
{
    const REAL K = fER.K(), G = fER.G();
    const REAL sinphi = sin ( fPhi );
    const REAL sinpsi = sin ( fPsi );
    const REAL cosphi = cos ( fPhi );
    const REAL cotphi = 1./tan ( fPhi );
    T ptrnp1 = 0.;
    for ( int i = 0; i < 3; i++ ) {
        ptrnp1 += T ( sigmatrial[i] );
    }
    ptrnp1 /= 3.;
    T DEpsPV = 0.;
    T epsbarnp1 = T ( fEpsPlasticBar );
    T c,H;
    PlasticityFunction ( epsbarnp1, c, H );
    //std::cout << "ReturnMap do Apex: c = " << c << "\tH = " << H << "\tptrnp1 = " << ptrnp1 << std::endl;
    T alpha = cos ( fPhi ) /sin ( fPsi );
    REAL tol = 1.e-8;

    T res = c*cotphi-ptrnp1;
    T pnp1;

    for ( int i = 0; i < 30; i++ ) {
        const T d = H*T ( cosphi*cotphi ) /T ( sinpsi ) + T ( K );
        DEpsPV -= res/d;

        epsbarnp1 = T ( fEpsPlasticBar )+T ( alpha ) *DEpsPV;
        pnp1 = ptrnp1 - T ( K ) * DEpsPV;
        PlasticityFunction ( epsbarnp1, c, H );
        res = c*cotphi - pnp1;
        if ( fabs ( res ) < tol ) break;
    }
    epsbarnew = shapeFAD::val ( epsbarnp1 );

    for ( int i = 0; i < 3; i++ ) {
        sigma_projected[i] = pnp1;
    }
    return true; // If it is in this ReturnMap it surely is this type of ReturnMap (ProjectSigma manages this)
}

void TPZYCMohrCoulombPV::ComputeApexTangent ( TPZMatrix<REAL> &tang, REAL &epsbarp ) const
{
    REAL c,H;
    const REAL cosphi = cos ( fPhi );
    const REAL sinpsi	= sin ( fPsi );
    const REAL cotphi = 1./tan ( fPhi );
    const REAL K = fER.K();
    const REAL alpha = cosphi/sinpsi;
    this->PlasticityFunction ( epsbarp, c, H );
    const REAL num = H*alpha*cotphi/K;
    const REAL denom = 1. + num;
    const REAL dpdptr = num/denom;
    const REAL dsigdsigtr = dpdptr/3.;
    tang.Redim ( 3, 3 );
    for ( int i = 0; i < 3; i++ ) {
        for ( int j = 0; j < 3 ; j++ ) {
            tang ( i,j ) = dsigdsigtr;
        }
    }
}

void TPZYCMohrCoulombPV::ProjectSigma ( const TPZVec<STATE> &sigma_trial, STATE eprev, TPZVec<STATE> &sigma,STATE &eproj )
{
    this->SetEpsBar ( eprev );
    REAL epsbartemp = eprev; // it will be defined by the correct returnmap
    TComputeSequence memory;
    REAL phi = PhiPlane<REAL> ( sigma_trial );
    if ( phi <= 0. ) {
		//cout << "elastic"<<endl;
        memory.fWhichPlane = TComputeSequence::EElastic;
        memory.fGamma.Resize ( 0 );
        sigma = sigma_trial;
        return;
    }
    TPZVec<REAL> sigma_projected;
    memory.fGamma.Resize ( 1 );
    memory.fGamma[0] = 0.;
    if ( this->ReturnMapPlane<REAL> ( sigma_trial, sigma_projected, memory, epsbartemp ) ) {
        eproj = epsbartemp;
        this->SetEpsBar ( eproj );
        sigma = sigma_projected;
        memory.fWhichPlane = TComputeSequence::EMainPlane;
		//cout << "plane"<<endl;
    } else {
        memory.fGamma.Resize ( 2 );
        memory.fGamma[0] = 0.;
        memory.fGamma[1] = 0.;
        bool IsEdge = false;

        const REAL sinpsi = sin ( fPsi );
        REAL val = ( 1-sinpsi ) *sigma_trial[0]-2.*sigma_trial[1]+ ( 1+sinpsi ) *sigma_trial[2];
        if ( val > 0. ) {
			//cout << "ERightEdge"<<endl;
            IsEdge = this->ReturnMapRightEdge<REAL> ( sigma_trial, sigma_projected, memory,epsbartemp );
            memory.fWhichPlane = TComputeSequence::ERightEdge;
        } else {
			//cout << "ELeftEdge"<<endl;
            IsEdge = this->ReturnMapLeftEdge<REAL> ( sigma_trial, sigma_projected, memory,epsbartemp );
            memory.fWhichPlane = TComputeSequence::ELeftEdge;
        }
        if ( !IsEdge ) {
			//cout << "apex"<<endl;
            this->ReturnMapApex ( sigma_trial, sigma_projected, memory,epsbartemp );
            memory.fWhichPlane = TComputeSequence::EApex;
        }

        eproj = epsbartemp;
        this->SetEpsBar ( eproj );
        sigma = sigma_projected;
    }
}
/*
void TPZYCMohrCoulombPV::ProjectSigmaDep(const TPZVec<STATE> &sigma_trial,TPZTensor<REAL>sigtrtensor, STATE eprev, TPZVec<STATE> &sigma, STATE &eproj, TPZFMatrix<STATE> &GradSigma)
{
	this->SetEpsBar(eprev);
	REAL epsbartemp = -6738.; // it will be defined by the correct returnmap
	TComputeSequence memory;
	REAL phi = PhiPlane<REAL>(sigma_trial);

//     cout << "trial tensor" << endl;
//  	cout << sigtrtensor.XX()  << endl;
//  	cout << sigtrtensor.YY()  << endl;
//     cout << sigtrtensor.ZZ()  << endl;
//     cout << sigtrtensor.XZ()  << endl;
//     cout << sigtrtensor.YZ()  << endl;
//     cout << sigtrtensor.XY()  << endl;


	if (phi <= 0.) {
		memory.fWhichPlane = TComputeSequence::EElastic;
		//GradSigma.Identity();
		memory.fGamma.Resize(0);
		eproj = eprev;
		sigma = sigma_trial;
        GradSigma=GetElasticMatrix();
		return;
	}
	TPZManVector<REAL, 3> strial = sigma_trial;
           // cout << "trial" << endl;
  	//cout << sigma_trial[0]  << endl;
  	//cout << sigma_trial[1]  << endl;
    // cout << sigma_trial[2]  << endl;
	TPZManVector<REAL,3> sigma_projected;
	memory.fGamma.Resize(1);
	memory.fGamma[0] = 0.;
	if (this->ReturnMapPlane<REAL>(sigma_trial, sigma_projected, memory, epsbartemp)) {
		eproj = epsbartemp;
		//this->ComputePlaneTangent(GradSigma, epsbartemp);
        ComputeConsistentPlaneTangent (sigtrtensor,GradSigma );
		sigma = sigma_projected;
		memory.fWhichPlane = TComputeSequence::EMainPlane;
		  //          cout << "sigma_projected" << endl;
  	///cout << sigma_projected[0]  << endl;
  	//cout << sigma_projected[1]  << endl;
    // cout << sigma_projected[2]  << endl;
	}
	else {
		memory.fGamma.Resize(2);
		memory.fGamma[0] = 0.;
		memory.fGamma[1] = 0.;
		bool IsEdge = false, IsRight = false;

		const REAL sinpsi = sin(fPsi);
		REAL val = (1-sinpsi)*sigma_trial[0]-2.*sigma_trial[1]+(1+sinpsi)*sigma_trial[2];
		if (val > 0.) {
			IsEdge = this->ReturnMapRightEdge<REAL>(sigma_trial, sigma_projected, memory,epsbartemp);
			memory.fWhichPlane = TComputeSequence::ERightEdge;
			IsRight = true;
            int type=1;
            ComputeConsistentEdgeTangent ( sigtrtensor,GradSigma, type);
		}
		else {
			IsEdge = this->ReturnMapLeftEdge<REAL>(sigma_trial, sigma_projected, memory,epsbartemp);
			memory.fWhichPlane = TComputeSequence::ELeftEdge;
			IsRight = false;
            int type=0;
            ComputeConsistentEdgeTangent ( sigtrtensor,GradSigma, type);
		}
		if (!IsEdge) {
			this->ReturnMapApex(sigma_trial, sigma_projected, memory,epsbartemp);
			memory.fWhichPlane = TComputeSequence::EApex;
            //this->ComputeApexTangent(GradSigma, epsbartemp);
            GradSigma.Resize(6,6);
            GradSigma.Zero();


		}
// 		if (IsEdge && IsRight) this->ComputeRightEdgeTangent(GradSigma, epsbartemp);
// 		else if (IsEdge && !IsRight) this->ComputeLeftEdgeTangent(GradSigma, epsbartemp);
// 		else this->ComputeApexTangent(GradSigma, epsbartemp);

		eproj = epsbartemp;
		this->SetEpsBar(eproj);
		sigma = sigma_projected;
	}
}*/

void TPZYCMohrCoulombPV::ProjectSigmaDep ( const TPZVec<STATE> &sigma_trial,TPZTensor<REAL>sigtrtensor, STATE eprev, TPZVec<STATE> &sigma, STATE &eproj, TPZFMatrix<STATE> &GradSigma )
{
    this->SetEpsBar ( eprev );
    REAL epsbartemp = -6738.; // it will be defined by the correct returnmap
    TComputeSequence memory;
    REAL phi = PhiPlane<REAL> ( sigma_trial );

    if ( phi <= 0. ) {
        memory.fWhichPlane = TComputeSequence::EElastic;
        GradSigma.Identity();
        memory.fGamma.Resize ( 0 );
        eproj = eprev;
        sigma = sigma_trial;
        //GradSigma=GetElasticMatrix();
		//cout << "elastic"<<endl;
        return;
    }
    TPZManVector<REAL, 3> strial = sigma_trial;

    TPZManVector<REAL,3> sigma_projected;
    memory.fGamma.Resize ( 1 );
    memory.fGamma[0] = 0.;
    if ( this->ReturnMapPlane<REAL> ( sigma_trial, sigma_projected, memory, epsbartemp ) ) {
        eproj = epsbartemp;
        this->ComputePlaneTangent ( GradSigma, epsbartemp );
        //ComputeConsistentPlaneTangent (sigtrtensor,GradSigma );
        sigma = sigma_projected;
        memory.fWhichPlane = TComputeSequence::EMainPlane;
		//cout << "plane"<<endl;
    } else {
        memory.fGamma.Resize ( 2 );
        memory.fGamma[0] = 0.;
        memory.fGamma[1] = 0.;
        bool IsEdge = false, IsRight = false;

        const REAL sinpsi = sin ( fPsi );
        REAL val = ( 1-sinpsi ) *sigma_trial[0]-2.*sigma_trial[1]+ ( 1+sinpsi ) *sigma_trial[2];
        if ( val > 0. ) {
			//cout << "ERightEdge"<<endl;
            IsEdge = this->ReturnMapRightEdge<REAL> ( sigma_trial, sigma_projected, memory,epsbartemp );
            memory.fWhichPlane = TComputeSequence::ERightEdge;
            IsRight = true;
            int type=1;
            //ComputeConsistentEdgeTangent ( sigtrtensor,GradSigma, type);
        } else {
			//cout << "ELeftEdge"<<endl;
            IsEdge = this->ReturnMapLeftEdge<REAL> ( sigma_trial, sigma_projected, memory,epsbartemp );
            memory.fWhichPlane = TComputeSequence::ELeftEdge;
            IsRight = false;
            int type=0;

        }
        if ( !IsEdge ) {
			//cout << "EApex"<<endl;
            this->ReturnMapApex ( sigma_trial, sigma_projected, memory,epsbartemp );
            memory.fWhichPlane = TComputeSequence::EApex;
            this->ComputeApexTangent ( GradSigma, epsbartemp );
            GradSigma.Resize ( 6,6 );
            GradSigma.Zero();
        }
        if ( IsEdge && IsRight ) this->ComputeRightEdgeTangent ( GradSigma, epsbartemp );
        else if ( IsEdge && !IsRight ) this->ComputeLeftEdgeTangent ( GradSigma, epsbartemp );
        else this->ComputeApexTangent ( GradSigma, epsbartemp );

        eproj = epsbartemp;
        this->SetEpsBar ( eproj );
        sigma = sigma_projected;
    }
}

void TPZYCMohrCoulombPV::ComputeConsistentPlaneTangent ( TPZTensor<REAL> & trialstress,TPZFMatrix<REAL> & Dep )
{

	REAL sig1,sig2,sig3;
	//cout << "\n trialstress ";
	//trialstress.Print(cout);
	REAL I1=trialstress.I1();
    REAL J2=trialstress.J2();
    REAL J3=trialstress.J3();
    REAL val = -3*sqrt ( 3. ) *J3/ ( 2.*pow ( J2,1.5 ) );
    if ( val>1. ) {
        val=1.;
    }
    if ( val<-1. ) {
        val=-1.;
    }
    ftheta= 1/3.*asin ( val );

    fa = cos ( ftheta )-1./sqrt ( 3. ) *sin ( ftheta ) *sin ( fPhi );

    fda=- ( ( cos ( ftheta ) *sin ( fPhi ) ) /sqrt ( 3 ) ) - sin ( ftheta );

    fd2a = -cos ( ftheta ) + ( sin ( fPhi ) *sin ( ftheta ) ) /sqrt ( 3 );

	//cout << "\n ftheta " << ftheta << endl;
	//cout << "\n fa " << fa << endl;
    //Compute Deltalambda
    TPZFMatrix<REAL> E= fER.GetElasticMatrixReal();
    TPZFMatrix<REAL> atemp,denom;
    REAL phiinvars=PhiInvars ( trialstress );
	//REAL phiinvars = 2*(1./3. * I1 *sin ( fPhi )+sqrt ( J2 ) *fa-fc*cos ( fPhi ));
	//cout << "\n phi " << phiinvars << endl;
    //Doub phiinvars=Yields(trialstress)[0];
    TPZFMatrix<REAL> a,at;
    a=avec ( trialstress );
	//cout << "\n a "  << endl;
	//a.Print("a");
    //a =N(trialstress)[0];
    a.Transpose ( &at );
    at.Multiply ( E,atemp );
    atemp.Multiply ( a,denom );
    if ( denom.Rows() >1 ) DebugStop();
    REAL Deltalambda = phiinvars/denom ( 0,0 );
	//cout << "denom ( 0,0 )"<< denom ( 0,0 ) << endl;
    //cout << "Deltalambda"<< Deltalambda << endl;
    TPZFMatrix<REAL> dadsig = dAdsig ( trialstress ),Etmodular ( E ),InverseTemp,Q ( 6,6 );
    Q.Identity ( );

	//dadsig.Print("dadsig");
	
	//E.Print("E");
	
    TPZFMatrix<REAL>  Ea,atEt,atEa,Et,atE,num;
    E.Transpose ( &Et );
    E.Multiply ( a,Ea );
    at.Multiply ( Et,atEt );
    at.Multiply ( E,atE );
    atE.Multiply ( a,atEa );
    Ea.Multiply ( atEt,num );
    if ( atEa.Rows() >1 ) DebugStop();
    num*=1/atEa ( 0,0 );
    Etmodular-=num;

    TPZFMatrix<REAL> temp;

    dadsig.Multiply ( E,temp );
    temp*=Deltalambda;
	//Etmodular.Print("Etmodular");
    Q-=temp;
	//Q.Print("Q");
    Etmodular.Multiply ( Q,Dep );


}

void TPZYCMohrCoulombPV::ComputeConsistentEdgeTangent ( TPZTensor<REAL> &trialstress,TPZFMatrix<REAL> & Dep, int &m_type )
{

    REAL J2=trialstress.J2();
    REAL J3=trialstress.J3();
    REAL val = -3*sqrt ( 3. ) *J3/ ( 2.*pow ( J2,1.5 ) );
    if ( val>1. ) {
        val=1.;
    }
    if ( val<-1. ) {
        val=-1.;
    }
    ftheta= 1/3.*asin ( val );

    fa = cos ( ftheta )-1./sqrt ( 3. ) *sin ( ftheta ) *sin ( fPhi );

    fda=- ( ( cos ( ftheta ) *sin ( fPhi ) ) /sqrt ( 3 ) ) - sin ( ftheta );

    fd2a = -cos ( ftheta ) + ( sin ( fPhi ) *sin ( ftheta ) ) /sqrt ( 3 );


    //cout << "theta "<<endl;
    //cout <<ftheta*180/M_PI<<endl;
    // NRvector<NRmatrix<Doub>> nvecs = N ( trialstress );
    // NRvector<Doub> phis = Yields ( trialstress );
    TPZFMatrix<REAL>E = fER.GetElasticMatrixReal();
    TPZFMatrix<REAL> a,b,at,bt;

    a=avec ( trialstress );

    REAL ca,cb,cd,faa,f2a;
    faa=PhiInvars ( trialstress );
    //dadsig para thetax
    //a=nvecs[0];
    //faa=phis[0];

    TPZFMatrix<REAL> dbdsig,dadsig;
    TPZFMatrix<REAL> Ea,Eb,atE,btE,tempa,tempb,tempd;
    dadsig = dAdsig ( trialstress );
    if ( m_type==1 ) { //rigth f1 e f2
        // b=nvecs[1];
        //   f2a=phis[1];
        //dbdsig para theta30
        //ftheta = 30*M_PI/180.;
        fa= ( cos ( ftheta ) * ( 1 + sin ( fPhi ) ) ) /2. + ( ( -3 + sin ( fPhi ) ) *sin ( ftheta ) ) / ( 2.*sqrt ( 3 ) );
        fda= ( cos ( ftheta ) * ( -3 + sin ( fPhi ) ) ) / ( 2.*sqrt ( 3 ) ) - ( ( 1 + sin ( fPhi ) ) *sin ( ftheta ) ) /2.;
        fd2a=- ( cos ( ftheta ) * ( 1 + sin ( fPhi ) ) ) /2. - ( ( -3 + sin ( fPhi ) ) *sin ( ftheta ) ) / ( 2.*sqrt ( 3 ) );

    } else { //left f1 e f3
        //b=nvecs[2];
        //f2a=phis[2];
        //ftheta = -29.99*M_PI/180.;
        fa= ( cos ( ftheta ) * ( 1 - sin ( fPhi ) ) ) /2. + ( ( 3 + sin ( fPhi ) ) *sin ( ftheta ) ) / ( 2.*sqrt ( 3 ) );
        fda= ( cos ( ftheta ) * ( 3 + sin ( fPhi ) ) ) / ( 2.*sqrt ( 3 ) ) - ( ( 1 - sin ( fPhi ) ) *sin ( ftheta ) ) /2.;
        fd2a=- ( cos ( ftheta ) * ( 1 - sin ( fPhi ) ) ) /2. - ( ( 3 + sin ( fPhi ) ) *sin ( ftheta ) ) / ( 2.*sqrt ( 3 ) );
    }

    dbdsig = dAdsig ( trialstress );
    TPZFMatrix<REAL> symcheck;

    b=avec ( trialstress );
    f2a=PhiInvars ( trialstress );

    a.Transpose ( &at );
    b.Transpose ( &bt );

    at.Multiply ( E,atE );
    bt.Multiply ( E,btE );

    atE.Multiply ( a,tempa );
    btE.Multiply ( b,tempb );
    atE.Multiply ( b,tempd );
    REAL q = ( tempa ( 0,0 ) *tempb ( 0,0 )-pow ( tempd ( 0,0 ),2. ) );
    REAL DeltalambdaA = ( tempb ( 0,0 ) *faa-tempd ( 0,0 ) *f2a ) / q;
    REAL DeltalambdaB = ( tempa ( 0,0 ) *f2a-tempd ( 0,0 ) *faa ) / q;

    E.Multiply ( a,Ea );
    E.Multiply ( b,Eb );

    TPZFMatrix<REAL>temp1,temp2,temp3,temp4,Et2 ( E );

    Ea.Multiply ( atE,temp1 );
    temp1*=tempb ( 0,0 );

    Ea.Multiply ( btE,temp2 );
    temp2*=tempd ( 0,0 );

    Eb.Multiply ( atE,temp3 );
    temp3*=tempd ( 0,0 );

    Eb.Multiply ( btE,temp4 );
    temp4*=tempa ( 0,0 );

    temp1-=temp2;
    temp1-=temp3;
    temp1+=temp4;
    temp1*=1./q;
    Et2-=temp1;
    //Dep=Et2;
    //dadsig.Print();
    //dbdsig.Print();
    //DebugStop();

    TPZFMatrix<REAL> T ( 6,6 ),dadsigE,dbdsigE;
    T.Identity();


    dadsig.Multiply ( E,dadsigE );
    //dadsigE.Mult(a,symcheck);
    //cout << "symcheck" << symcheck[0][0] <<endl;

    dbdsig.Multiply ( E,dbdsigE );
    //dbdsigE.Mult(b,symcheck);
    //cout << "symcheck B" << symcheck[0][0] <<endl;


    dadsigE*=DeltalambdaA;
    dbdsigE*=DeltalambdaB;

    T-=dadsigE;
    T-=dbdsigE;
    //dbdsig.Print();
    Et2.Multiply ( T,Dep );

//   Dep.Print("DEP");
}



template void TPZYCMohrCoulombPV::PlasticityFunction<REAL> ( const REAL epsp, REAL &sigmay, REAL &H ) const;
template void TPZYCMohrCoulombPV::PlasticityFunction<fadtype> ( const fadtype epsp, fadtype &sigmay, fadtype &H ) const;

template TPZVec<REAL> TPZYCMohrCoulombPV::SigmaElastPV<REAL> ( const TPZVec<REAL> &deform ) const;
template TPZVec<fadtype> TPZYCMohrCoulombPV::SigmaElastPV<fadtype> ( const TPZVec<fadtype> &deform ) const;

template REAL TPZYCMohrCoulombPV::PhiPlane<REAL> ( const TPZVec<REAL> &sigma ) const;
template fadtype TPZYCMohrCoulombPV::PhiPlane<fadtype> ( const TPZVec<fadtype> &sigma ) const;

template bool TPZYCMohrCoulombPV::ReturnMapPlane<REAL> ( const TPZVec<REAL> &sigma_trial, TPZVec<REAL> &sigma_projected,
        TComputeSequence &memory, REAL &epsbarnew ) const;
template bool TPZYCMohrCoulombPV::ReturnMapPlane<fadtype> ( const TPZVec<fadtype> &sigma_trial, TPZVec<fadtype> &sigma_projected,
        TComputeSequence &memory, REAL &epsbarnew ) const;

template bool TPZYCMohrCoulombPV::ReturnMapLeftEdge<REAL> ( const TPZVec<REAL> &sigma_trial, TPZVec<REAL> &sigma_projected,
        TComputeSequence &memory, REAL &epsbarnew ) const;
template bool TPZYCMohrCoulombPV::ReturnMapLeftEdge<fadtype> ( const TPZVec<fadtype> &sigma_trial, TPZVec<fadtype> &sigma_projected,
        TComputeSequence &memory, REAL &epsbarnew ) const;

template bool TPZYCMohrCoulombPV::ReturnMapRightEdge<REAL> ( const TPZVec<REAL> &sigma_trial, TPZVec<REAL> &sigma_projected,
        TComputeSequence &memory, REAL &epsbarnew ) const;
template bool TPZYCMohrCoulombPV::ReturnMapRightEdge<fadtype> ( const TPZVec<fadtype> &sigma_trial, TPZVec<fadtype> &sigma_projected,
        TComputeSequence &memory, REAL &epsbarnew ) const;

template bool TPZYCMohrCoulombPV::ReturnMapApex<REAL> ( const TPZVec<REAL> &sigma_trial, TPZVec<REAL> &sigma_projected,
        TComputeSequence &memory, REAL &epsbarnew ) const;
template bool TPZYCMohrCoulombPV::ReturnMapApex<fadtype> ( const TPZVec<fadtype> &sigma_trial, TPZVec<fadtype> &sigma_projected,
        TComputeSequence &memory, REAL &epsbarnew ) const;
