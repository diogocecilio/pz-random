#include "TPZMohrCoulombVoigt.h"
#include "TPZTensor.h"

TPZMohrCoulombVoigt::TPZMohrCoulombVoigt() : fPhi ( 0. ),fPsi ( 0. ), fc ( 0. ), fEpsPlasticBar ( 0. ), fER()
{

}

TPZMohrCoulombVoigt::TPZMohrCoulombVoigt ( REAL Phi, REAL Psi, REAL c, TPZElasticResponse &ER ) : fPhi ( Phi ), fPsi ( Psi ), fc ( c ), fEpsPlasticBar ( 0. ), fER ( ER )
{

}

TPZMohrCoulombVoigt::TPZMohrCoulombVoigt ( const TPZMohrCoulombVoigt &cp ) // : 	fPhi(cp.fPhi), fPsi(cp.fPsi), fc(cp.fc), fEpsPlasticBar(cp.fEpsPlasticBar), fER(cp.fER)
{
    TPZMohrCoulombVoigt::operator= ( cp );
}

TPZMohrCoulombVoigt & TPZMohrCoulombVoigt::operator= ( const TPZMohrCoulombVoigt &cp )
{
    fPhi = cp.fPhi;
    fPsi = cp.fPsi;
    fc = cp.fc;
    fEpsPlasticBar = cp.fEpsPlasticBar;
    fER = cp.fER;
    return *this;
}


void TPZMohrCoulombVoigt::Read ( TPZStream &buf )
{
    buf.Read ( &fPhi );
    buf.Read ( &fPsi );
    buf.Read ( &fc );
    buf.Read ( &fEpsPlasticBar );
    fER.Read ( buf );

}

void TPZMohrCoulombVoigt::Write ( TPZStream &buf ) const
{
    buf.Write ( &fPhi );
    buf.Write ( &fPsi );
    buf.Write ( &fc );
    buf.Write ( &fEpsPlasticBar );
    fER.Write ( buf );
}

template<class T>
T TPZMohrCoulombVoigt::PhiPlane ( const TPZVec<T> &sigma ) const
{

}


REAL TPZMohrCoulombVoigt::ReturnMapPlane (  TPZTensor<REAL> &sigma_trial, TPZTensor<REAL> &sigma_proj,TPZFMatrix<REAL>&dep, REAL &epsbarnew)
{

    TPZTensor<REAL> avec;
     REAL a =A(sigma_trial);
     REAL I1 = sigma_trial.I1();
     REAL J2 = sigma_trial.J2();
     REAL dadt = dAdt(sigma_trial);
     REAL d2adt = d2Adt(sigma_trial);

     REAL c1=C1();
     REAL c2=C2 ( sigma_trial, a, dadt);
     REAL c3 = C3 ( sigma_trial,dadt );
     REAL c4=C4 (  sigma_trial, dadt,d2adt);
     TPZTensor<REAL> di1 = sigma_trial.dI1();
     TPZTensor<REAL> dj2 = sigma_trial.dJ2();
     TPZTensor<REAL> dj3 = sigma_trial.dJ3();

     REAL f1 = -(fc*cos(fPhi)) + a*sqrt(J2) + 0.3333333333333333*I1*sin(fPhi);

     di1*=c1;
     dj2*=c2;
     dj3*=c3;

     //TPZTensor<REAL> avec;
     avec=di1;
     avec+=dj2;
     avec+=dj3;

   //  avec.Print(std::cout);

//
      REAL dlvoigt=f1;
      TPZFMatrix<REAL> cmat = fER.GetElasticMatrixReal();
     ///     TPZFMatrix<REAL> cmat;
  ///  fER.ElasticMat(cmat);
//
    //  cmat.Print(std::cout);

      TPZFMatrix<REAL> amat,amatT,temp,temp2,tempT;

      avec.FromTensorToNRmatrix(amat);
      //FromTensorToMatVoigt(avec,amat);

      amat.Transpose(&amatT);

      cmat.Multiply(amat,temp);
      temp.Transpose(&tempT);
      //amat.Print(std::cout);
      //temp.Print(std::cout);
      tempT.Multiply(amat,temp2);
     //cmat.Multiply(avec,

      dlvoigt *= 1./temp2.Get(0,0);
      //sigtr - dlvoigt cmat.a;
      //FromMatToTensor(temp,sigma_proj);
      sigma_proj.CopyFrom(temp);

      sigma_proj*=-dlvoigt;
      sigma_proj+=sigma_trial;

      //sigma_proj.Print(std::cout);

      TPZFMatrix<REAL> Q1(6,6,0.),dadsig;
      Q1(0,0)=1.;Q1(1,1)=1.;Q1(2,2)=1.;Q1(3,3)=1.;Q1(4,4)=1.;Q1(5,5)=1.;

      dadsig =  dAdsig(  sigma_trial,  a, dadt, d2adt );

      TPZFMatrix<REAL> tempx,tempy,tempz,ect;

      dadsig.Multiply(cmat,tempx);
      tempx*=dlvoigt;
      Q1-=tempx;
      tempx.Zero();

      amat.Multiply(amatT,tempx);

      cmat.Multiply(tempx,tempy);

      tempy.Multiply(cmat,tempz);

      tempz*=1./temp2.Get(0,0);

      ect=cmat;

      ect-=tempz;

      ect.Multiply(Q1,dep);

     // dep.Print(std::cout);

      return dlvoigt;
}

/*
bool TPZMohrCoulombVoigt::ReturnMapLeftEdge ( TPZTensor<REAL> &sigma_trial, TPZTensor<REAL> &sigma_proj,TPZFMatrix<REAL>&dep, REAL &epsbarnew)
{
    TPZTensor<REAL> avec,bvec;
    TPZFMatrix<REAL> cmat = fER.GetElasticMatrixReal();
     REAL I1 = sigma_trial.I1();
     REAL J2 = sigma_trial.J2();

     REAL a =A(sigma_trial);
     REAL dadt = dAdt(sigma_trial);
     REAL d2adt = d2Adt(sigma_trial);

     FlowVector(sigma_trial,a,dadt,d2adt,avec);

     REAL a3 =A3(sigma_trial);
     REAL da3dt = dA3dt(sigma_trial);
     REAL d2a3dt = d2A3dt(sigma_trial);

     FlowVector(sigma_trial,a3,da3dt,d2a3dt,bvec);


     TPZFMatrix<REAL> temp,temp2,temp3,temp4;


     avec.FromTensorToNRmatrix(temp);
     //FromTensorToMatVoigt(avec,temp);
     cmat.Multiply(temp,temp2);
     REAL ax = Dot(temp,temp2);



     //FromTensorToMatVoigt(bvec,temp);
     bvec.FromTensorToNRmatrix(temp);
     cmat.Multiply(temp,temp2);
     REAL bx =Dot(temp,temp2);


     avec.FromTensorToNRmatrix(temp);
     //FromTensorToMatVoigt(avec,temp);
     //FromTensorToMatVoigt(bvec,temp2);
     bvec.FromTensorToNRmatrix(temp2);
     cmat.Multiply(temp2,temp3);
     REAL dx = Dot(temp,temp3);

     REAL f11,f22;

     Yield(sigma_trial,a,f11);
     Yield(sigma_trial,a3,f22);


   REAL q = ax*bx - dx*dx;
   REAL  dl11 = (bx*f11 - dx*f22)/(ax* bx - dx*dx);
   REAL dl22 = (ax*f22 - dx*f11)/(ax*bx - dx*dx);


   avec.FromTensorToNRmatrix(temp);
   //FromTensorToMatVoigt(avec,temp);
   bvec.FromTensorToNRmatrix(temp2);
   //FromTensorToMatVoigt(bvec,temp2);
   cmat.Multiply(temp,temp3);
   cmat.Multiply(temp2,temp4);
   temp3*=-dl11;
   temp4*=-dl22;

   temp3+=temp4;

   sigma_proj.CopyFrom(temp3);
   //FromMatToTensor(temp3,sigma_proj);
   sigma_proj+=sigma_trial;

   TPZFMatrix<REAL> dadsig1 = dAdsig(sigma_trial,  a, dadt, d2adt);
   TPZFMatrix<REAL> dadsig2 = dAdsig(sigma_trial,  a3, da3dt, d2a3dt);

   //Outer[Times, avec, avec]

    TPZFMatrix<REAL> avectemp,avectempT,bvectemp,bvectempT,aaT,abT,baT,bbT;

    avec.FromTensorToNRmatrix(avectemp);
    //FromTensorToMatVoigt(avec,avectemp);
    bvec.FromTensorToNRmatrix(bvectemp);
    //FromTensorToMatVoigt(bvec,bvectemp);

    avectemp.Transpose(&avectempT);
    bvectemp.Transpose(&bvectempT);

    //Outer[Times, avec, avec]
    avectemp.Multiply(avectempT,aaT);
    //Outer[Times, avec, bvec]
    avectemp.Multiply(bvectempT,abT);
    //Outer[Times, bvec, avec]
    bvectemp.Multiply(avectempT,baT);
    //Outer[Times, bvec, bvec]
    bvectemp.Multiply(bvectempT,bbT);

    TPZFMatrix<REAL> cmataaT,cmataaTcmat,cmatabT,cmatabTcmat,cmatbaT,cmatbaTcmat,cmatbbT,cmatbbTcmat,tempfinal,et2;

    //cmat.Outer[Times, avec, avec]
    cmat.Multiply(aaT,cmataaT);
    //cmat.Outer[Times, avec, avec].cmat
    cmataaT.Multiply(cmat,cmataaTcmat);

    cmataaTcmat*=bx;

    cmat.Multiply(abT,cmatabT);
    cmatabT.Multiply(cmat,cmatabTcmat);

    cmatabTcmat*=dx;


    cmat.Multiply(baT,cmatbaT);
    cmatbaT.Multiply(cmat,cmatbaTcmat);

    cmatbaTcmat*=dx;

    cmat.Multiply(bbT,cmatbbT);
    cmatbbT.Multiply(cmat,cmatbbTcmat);

    cmatbbTcmat*=ax;

    tempfinal=cmataaTcmat;

    tempfinal-=cmatabTcmat;
    tempfinal-=cmatbaTcmat;
    tempfinal+=cmatbbTcmat;

    tempfinal*=1./q;

    et2=cmat;

    et2-=tempfinal;

    TPZFMatrix<REAL> T(6,6,0.);
    T(0,0)=1.;T(1,1)=1.;T(2,2)=1.;T(3,3)=1.;T(4,4)=1.;T(5,5)=1.;

    TPZFMatrix<REAL> partea,parteb;
    dadsig1.Multiply(cmat,partea);
    dadsig2.Multiply(cmat,parteb);
    partea*=-dl11;
    parteb*=-dl22;
    T+=partea;
    T+=parteb;

    et2.Multiply(T,dep);

    return true;

}*/

bool TPZMohrCoulombVoigt::ReturnMapLeftEdge ( TPZTensor<REAL> &sigma_trial, TPZTensor<REAL> &sigma_proj,TPZFMatrix<REAL>&dep, REAL &epsbarnew)
{
    TPZTensor<REAL> avec,bvec;
    TPZFMatrix<REAL> cmat = fER.GetElasticMatrixReal();
      //  TPZFMatrix<REAL> cmat;
    //fER.ElasticMat(cmat);


     REAL a =A(sigma_trial);
     REAL dadt = dAdt(sigma_trial);
     REAL d2adt = d2Adt(sigma_trial);

     avec=FlowVectorMain(sigma_trial);

     REAL a2 =A3(sigma_trial);
     REAL da2dt = dA3dt(sigma_trial);
     REAL d2a2dt = d2A3dt(sigma_trial);

     bvec=FlowVectorLeft(sigma_trial);

     TPZFMatrix<REAL> temp,temp2,temp3,temp4;

     avec.FromTensorToNRmatrix(temp);
     //FromTensorToMatVoigt(avec,temp);
     cmat.Multiply(temp,temp2);
     REAL ax = Dot(temp,temp2);


     bvec.FromTensorToNRmatrix(temp);
     //FromTensorToMatVoigt(bvec,temp);
     cmat.Multiply(temp,temp2);
     REAL bx =Dot(temp,temp2);


     avec.FromTensorToNRmatrix(temp);
     bvec.FromTensorToNRmatrix(temp2);

     cmat.Multiply(temp2,temp3);
     REAL dx = Dot(temp,temp3);

     REAL f11,f22;

     f11= FMain(sigma_trial);
     f22=FLeft(sigma_trial);

   REAL q = ax*bx - dx*dx;
   REAL  dl11 = (bx*f11 - dx*f22)/(ax* bx - dx*dx);
   REAL dl22 = (ax*f22 - dx*f11)/(ax*bx - dx*dx);
//cout << "bb "  <<endl;

    avec.FromTensorToNRmatrix(temp);
    bvec.FromTensorToNRmatrix(temp2);
   //FromTensorToMatVoigt(avec,temp);
   //FromTensorToMatVoigt(bvec,temp2);
   cmat.Multiply(temp,temp3);
   cmat.Multiply(temp2,temp4);
   temp3*=-dl11;
   temp4*=-dl22;

   temp3+=temp4;

   sigma_proj.CopyFrom(temp3);
   //FromMatToTensor(temp3,sigma_proj);
   sigma_proj+=sigma_trial;

   TPZFMatrix<REAL> dadsig1 = dAdsig(sigma_trial,  a, dadt, d2adt);
   TPZFMatrix<REAL> dadsig2 = dAdsig(sigma_trial,  a2, da2dt, d2a2dt);

   //Outer[Times, avec, avec]

    TPZFMatrix<REAL> avectemp,avectempT,bvectemp,bvectempT,aaT,abT,baT,bbT;

    avec.FromTensorToNRmatrix(avectemp);
    bvec.FromTensorToNRmatrix(bvectemp);
    //FromTensorToMatVoigt(avec,avectemp);
    //FromTensorToMatVoigt(bvec,bvectemp);

    avectemp.Transpose(&avectempT);
    bvectemp.Transpose(&bvectempT);

    //Outer[Times, avec, avec]
    avectemp.Multiply(avectempT,aaT);
    //Outer[Times, avec, bvec]
    avectemp.Multiply(bvectempT,abT);
    //Outer[Times, bvec, avec]
    bvectemp.Multiply(avectempT,baT);
    //Outer[Times, bvec, bvec]
    bvectemp.Multiply(bvectempT,bbT);

    TPZFMatrix<REAL> cmataaT,cmataaTcmat,cmatabT,cmatabTcmat,cmatbaT,cmatbaTcmat,cmatbbT,cmatbbTcmat,tempfinal,et2;

    //cmat.Outer[Times, avec, avec]
    cmat.Multiply(aaT,cmataaT);
    //cmat.Outer[Times, avec, avec].cmat
    cmataaT.Multiply(cmat,cmataaTcmat);

    cmataaTcmat*=bx;

    cmat.Multiply(abT,cmatabT);
    cmatabT.Multiply(cmat,cmatabTcmat);

    cmatabTcmat*=dx;


    cmat.Multiply(baT,cmatbaT);
    cmatbaT.Multiply(cmat,cmatbaTcmat);

    cmatbaTcmat*=dx;

    cmat.Multiply(bbT,cmatbbT);
    cmatbbT.Multiply(cmat,cmatbbTcmat);

    cmatbbTcmat*=ax;

    tempfinal=cmataaTcmat;

    tempfinal-=cmatabTcmat;
    tempfinal-=cmatbaTcmat;
    tempfinal+=cmatbbTcmat;

    tempfinal*=1./q;

    et2=cmat;

    et2-=tempfinal;

    TPZFMatrix<REAL> T(6,6,0.);
    T(0,0)=1.;T(1,1)=1.;T(2,2)=1.;T(3,3)=1.;T(4,4)=1.;T(5,5)=1.;

    TPZFMatrix<REAL> partea,parteb;
    dadsig1.Multiply(cmat,partea);
    dadsig2.Multiply(cmat,parteb);
    partea*=-dl11;
    parteb*=-dl22;
    T+=partea;
    T+=parteb;

    et2.Multiply(T,dep);
    //cout << "cc "  <<endl;

    return true;

}
bool TPZMohrCoulombVoigt::ReturnMapRightEdge ( TPZTensor<REAL> &sigma_trial, TPZTensor<REAL> &sigma_proj,TPZFMatrix<REAL>&dep, REAL &epsbarnew)
{
    //cout << "aa "  <<endl;
    TPZTensor<REAL> avec,bvec;
    TPZFMatrix<REAL> cmat = fER.GetElasticMatrixReal();
    //TPZFMatrix<REAL> cmat;
    //fER.ElasticMat(cmat);
     REAL I1 = sigma_trial.I1();
     REAL J2 = sigma_trial.J2();

     REAL a =A(sigma_trial);
     REAL dadt = dAdt(sigma_trial);
     REAL d2adt = d2Adt(sigma_trial);

     FlowVector(sigma_trial,a,dadt,d2adt,avec);

     REAL a2 =A2(sigma_trial);
     REAL da2dt = dA2dt(sigma_trial);
     REAL d2a2dt = d2A2dt(sigma_trial);

     FlowVector(sigma_trial,a2,da2dt,d2a2dt,bvec);


     TPZFMatrix<REAL> temp,temp2,temp3,temp4;


     avec.FromTensorToNRmatrix(temp);
     //FromTensorToMatVoigt(avec,temp);
     cmat.Multiply(temp,temp2);
     REAL ax = Dot(temp,temp2);


     bvec.FromTensorToNRmatrix(temp);
     //FromTensorToMatVoigt(bvec,temp);
     cmat.Multiply(temp,temp2);
     REAL bx =Dot(temp,temp2);


     avec.FromTensorToNRmatrix(temp);
     bvec.FromTensorToNRmatrix(temp2);
     //FromTensorToMatVoigt(avec,temp);
     //FromTensorToMatVoigt(bvec,temp2);
     cmat.Multiply(temp2,temp3);
     REAL dx = Dot(temp,temp3);

     REAL f11,f22;

     Yield(sigma_trial,a,f11);
     Yield(sigma_trial,a2,f22);

   REAL q = ax*bx - dx*dx;
   REAL  dl11 = (bx*f11 - dx*f22)/(ax* bx - dx*dx);
   REAL dl22 = (ax*f22 - dx*f11)/(ax*bx - dx*dx);
//cout << "bb "  <<endl;

    avec.FromTensorToNRmatrix(temp);
    bvec.FromTensorToNRmatrix(temp2);
   //FromTensorToMatVoigt(avec,temp);
   //FromTensorToMatVoigt(bvec,temp2);
   cmat.Multiply(temp,temp3);
   cmat.Multiply(temp2,temp4);
   temp3*=-dl11;
   temp4*=-dl22;

   temp3+=temp4;

   sigma_proj.CopyFrom(temp3);
   //FromMatToTensor(temp3,sigma_proj);
   sigma_proj+=sigma_trial;

   TPZFMatrix<REAL> dadsig1 = dAdsig(sigma_trial,  a, dadt, d2adt);
   TPZFMatrix<REAL> dadsig2 = dAdsig(sigma_trial,  a2, da2dt, d2a2dt);

   //Outer[Times, avec, avec]

    TPZFMatrix<REAL> avectemp,avectempT,bvectemp,bvectempT,aaT,abT,baT,bbT;

    avec.FromTensorToNRmatrix(avectemp);
    bvec.FromTensorToNRmatrix(bvectemp);
    //FromTensorToMatVoigt(avec,avectemp);
    //FromTensorToMatVoigt(bvec,bvectemp);

    avectemp.Transpose(&avectempT);
    bvectemp.Transpose(&bvectempT);

    //Outer[Times, avec, avec]
    avectemp.Multiply(avectempT,aaT);
    //Outer[Times, avec, bvec]
    avectemp.Multiply(bvectempT,abT);
    //Outer[Times, bvec, avec]
    bvectemp.Multiply(avectempT,baT);
    //Outer[Times, bvec, bvec]
    bvectemp.Multiply(bvectempT,bbT);

    TPZFMatrix<REAL> cmataaT,cmataaTcmat,cmatabT,cmatabTcmat,cmatbaT,cmatbaTcmat,cmatbbT,cmatbbTcmat,tempfinal,et2;

    //cmat.Outer[Times, avec, avec]
    cmat.Multiply(aaT,cmataaT);
    //cmat.Outer[Times, avec, avec].cmat
    cmataaT.Multiply(cmat,cmataaTcmat);

    cmataaTcmat*=bx;

    cmat.Multiply(abT,cmatabT);
    cmatabT.Multiply(cmat,cmatabTcmat);

    cmatabTcmat*=dx;


    cmat.Multiply(baT,cmatbaT);
    cmatbaT.Multiply(cmat,cmatbaTcmat);

    cmatbaTcmat*=dx;

    cmat.Multiply(bbT,cmatbbT);
    cmatbbT.Multiply(cmat,cmatbbTcmat);

    cmatbbTcmat*=ax;

    tempfinal=cmataaTcmat;

    tempfinal-=cmatabTcmat;
    tempfinal-=cmatbaTcmat;
    tempfinal+=cmatbbTcmat;

    tempfinal*=1./q;

    et2=cmat;

    et2-=tempfinal;

    TPZFMatrix<REAL> T(6,6,0.);
    T(0,0)=1.;T(1,1)=1.;T(2,2)=1.;T(3,3)=1.;T(4,4)=1.;T(5,5)=1.;

    TPZFMatrix<REAL> partea,parteb;
    dadsig1.Multiply(cmat,partea);
    dadsig2.Multiply(cmat,parteb);
    partea*=-dl11;
    parteb*=-dl22;
    T+=partea;
    T+=parteb;

    et2.Multiply(T,dep);
    //cout << "cc "  <<endl;

    return true;
}

bool TPZMohrCoulombVoigt::ReturnMapApex ( TPZTensor<REAL> &sigma_trial, TPZTensor<REAL> &sigma_proj,TPZFMatrix<REAL>&dep, REAL &epsbarnew )
{
    REAL ccotanphi = fc*cos(fPhi)/sin(fPhi);
    sigma_proj.XX()=ccotanphi;sigma_proj.YY()=ccotanphi;sigma_proj.ZZ()=ccotanphi;
    sigma_proj.XZ()=0.;sigma_proj.XY()=0.;sigma_proj.YZ()=0.;
    TPZFMatrix<REAL> dumb(6,6,0.);
    TPZFMatrix<REAL> cmat=fER.GetElasticMatrixReal();
    //TPZFMatrix<REAL> cmat;
    //fER.ElasticMat(cmat);
    REAL young = fER.E();
    cmat*=0.;
    dep=cmat;
    return true;

// TPZFMatrix<REAL> cmat=fER.GetElasticMatrixReal();
//     const REAL K = fER.K(), G = fER.G();
//     const REAL sinphi = sin ( fPhi );
//     const REAL sinpsi = sin ( fPsi );
//     const REAL cosphi = cos ( fPhi );
//     const REAL cotphi = 1./tan ( fPhi );
//     TPZManVector<REAL,3>  sigmatrial= EigenVal(sigma_trial);
//     REAL ptrnp1 = 0.;
//     for ( int i = 0; i < 3; i++ ) {
//         ptrnp1 +=  sigmatrial[i] ;
//     }
//     ptrnp1 /= 3.;
//     REAL DEpsPV = 0.;
//     REAL epsbarnp1 = epsbarnew;
//     REAL c=fc,H=0.;
//
//     REAL alpha = cos ( fPhi ) /sin ( fPsi );
//
//     REAL res = c*cotphi-ptrnp1;
//     REAL pnp1;
//
//
//     const REAL d =  K ;
//     DEpsPV -= res/d;
//
//     pnp1 = ptrnp1 - K  * DEpsPV;
//
//
//     sigma_proj.XX()=pnp1;sigma_proj.YY()=pnp1;sigma_proj.ZZ()=pnp1;
//     sigma_proj.XZ()=0.;sigma_proj.XY()=0.;sigma_proj.YZ()=0.;
//
//     cmat*=0.0000001;
//     dep=cmat;
//
//
//
//     return true; // If it is in this ReturnMap it surely is this type of ReturnMap (ProjectSigma manages this)

}

void TPZMohrCoulombVoigt::ProjectSigmaDep2 ( TPZTensor<REAL> &sigma_trial, TPZTensor<REAL> &sigma_proj,TPZFMatrix<REAL>&dep, REAL &epsbarnew )
{
    bool flag;

    dep.Resize(6,6);
    REAL f1 = FMain(sigma_trial);

    if(f1>0)
    {

        TPZTensor<REAL> avec,bvec,bvec1,bvec2,ab;

        avec = FlowVectorMain(sigma_trial);

        REAL dlamb=0.;
        REAL epsbarnew=0.;

        dlamb = ReturnMapPlane ( sigma_trial, sigma_proj,dep,  epsbarnew);

        ab = FlowVectorMain(sigma_proj);

         REAL beta1;

        beta1= acos( avec.Dot(ab)/(avec.Norm()*ab.Norm()) )*180/M_PI;

        //cout << "beta1 = " << beta1 <<endl;

        REAL f1proj = FMain(sigma_proj);

        if(fabs(f1proj)<0.001)
        {
            //cout<< "projetou na f1 = " << f1proj << endl;
            return;
        }
        else
        {
            REAL J2 = sigma_trial.J2();
            REAL D2 = -tan(3.* theta(sigma_trial))/(2.*J2);
            REAL D3 = -sqrt(3.)/( 2.*pow(J2,1.5)* cos(3.* theta(sigma_trial ) ));

            TPZTensor<REAL> dj2 = sigma_trial.dJ2();
            TPZTensor<REAL> dj3 = sigma_trial.dJ3();

            TPZFMatrix<REAL> dj2mat,dj3mat,cc;

            dj2.FromTensorToNRmatrix(dj2mat);
            dj3.FromTensorToNRmatrix(dj3mat);

            cc=D2*dj2mat+D3*dj3mat;

            TPZFMatrix<REAL>  temp,amat,temp2,cmat,avecmat,ccT;

            cc.Transpose(&ccT);

            avec.FromTensorToNRmatrix(avecmat);

            cmat=fER.GetElasticMatrixReal();
            //fER.ElasticMat(cmat);

            cmat.Multiply(avecmat,temp);

            ccT.Multiply(temp,temp2);

            temp2*=-dlamb;

            REAL dt = temp2.Get(0,0)*180./M_PI;
            if(dt<0)
            {
                flag= ReturnMapRightEdge(sigma_trial, sigma_proj,dep, epsbarnew);
                if(beta1>89.99)
                {
                   flag= ReturnMapApex(sigma_trial, sigma_proj,dep, epsbarnew);
                }
            }
            else
            {
                flag= ReturnMapLeftEdge(sigma_trial, sigma_proj,dep, epsbarnew);
                if(beta1>89.99)
                {
                   flag= ReturnMapApex(sigma_trial, sigma_proj,dep, epsbarnew);
                }
            }
        }
    }
    else
    {
        sigma_proj=sigma_trial;
        //cout << "k = "  <<endl;
        dep=fER.GetElasticMatrixReal();
        //fER.ElasticMat(dep);
        //cout << "l = "  <<endl;
        epsbarnew=0;

    }

}


void TPZMohrCoulombVoigt::ProjectSigmaDep( TPZTensor<REAL> &sigma_trial, TPZTensor<REAL> &sigma_proj,TPZFMatrix<REAL>&dep, REAL &epsbarnew )
{
    epsbarnew=0.;

    //cout << "sigma_trial"<< sigma_trial <<endl;
    REAL phi = FMain( sigma_trial );
    if ( phi <= 0. ) {
        //cout << "Elastic"<<endl;
        sigma_proj=sigma_trial;
        dep=fER.GetElasticMatrixReal();
        //fER.ElasticMat(dep);
        epsbarnew=0;
        return;
    }

    TPZManVector<REAL,3>  eigenvaltrial= EigenVal(sigma_trial);
    REAL sinpsi = sin(fPhi);
    REAL val = ( 1-sinpsi ) *eigenvaltrial[0]-2.*eigenvaltrial[1]+ ( 1+sinpsi ) *eigenvaltrial[2];
    bool rigth;
    if(val>=0)
    {
        rigth=true;
    }else{
        rigth=false;

    }
    REAL dlamb = ReturnMapPlane ( sigma_trial, sigma_proj,dep,  epsbarnew);

    REAL f1proj = FMain(sigma_proj);
    if(fabs(f1proj)<0.0001 )
    {
        //cout<< "projetou na f1 = " << f1proj << endl;
       // cout << "Main"<<endl;
//        dep=fER.GetElasticMatrixReal();
        return;

    }

//     REAL sz=sigma_trial.ZZ();
//     REAL s,diff1,diff2;
//     diff1=fabs(sz-eigenvaltrial[0]);
//     diff2=fabs(sz-eigenvaltrial[2]);
//     if(diff1<diff2)
//     {
//         s=-1.;
//     }else{
//         s=1.;
//     }
//     REAL perturb=sz+0.1*s*(eigenvaltrial[0]-eigenvaltrial[2]);
//     sigma_trial.ZZ()=perturb;

//     REAL sz2=sigma_proj.ZZ();
//     diff1=fabs(sz-eigenvaltrial[0]);
//     diff2=fabs(sz-eigenvaltrial[2]);
//     if(diff1<diff2)
//     {
//         s=-1.;
//     }else{
//         s=1.;
//     }
//     TPZManVector<REAL,3>  eigenvalproj= EigenVal(sigma_proj);
//     perturb=sz2+0.1*s*(eigenvalproj[0]-eigenvalproj[2]);
//     sigma_proj.ZZ()=perturb;

    if ( rigth) {
       //cout << "ERightEdge"<<endl;
        bool flag= ReturnMapRightEdge(sigma_trial, sigma_proj,dep, epsbarnew);

    } else {
       // cout << "ELeftEdge"<<endl;
        bool flag= ReturnMapLeftEdge(sigma_trial, sigma_proj,dep, epsbarnew);

    }


//dep=fER.GetElasticMatrixReal();
    TPZTensor<REAL> aa,ab;
    aa = FlowVectorMain(sigma_trial);
    ab = FlowVectorMain(sigma_proj);
    REAL beta= acos( aa.Dot(ab)/(aa.Norm()*ab.Norm()) )*180/M_PI;
    if ( beta>89.99) {
       // cout << "APEX"<<endl;
        //cout << "prval"<< prval << endl;
        bool flag= ReturnMapApex(sigma_trial, sigma_proj,dep, epsbarnew);

    }


}

 bool TPZMohrCoulombVoigt::CheckOrder(TPZManVector<REAL,3>  eigenvalues)
 {
    return (eigenvalues[0]  >= eigenvalues[1]  &&  eigenvalues[1] >=  eigenvalues[2]);
 }
TPZManVector<REAL,3>  TPZMohrCoulombVoigt::EigenVal(TPZTensor<REAL> & tensor)
 {
    TPZTensor<REAL>::TPZDecomposed DecompSig;
	tensor.EigenSystem(DecompSig);
	TPZManVector<REAL,3> eigenvalues(DecompSig.fEigenvalues);
    return eigenvalues;
 }
