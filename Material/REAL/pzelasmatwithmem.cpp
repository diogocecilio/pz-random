
#include "pzelasmatwithmem.h"
#include "pzbndcond.h"

template <class TMEM>
TPZElasticityMaterialWithMem<TMEM>::TPZElasticityMaterialWithMem(): TPZMatWithMem<TMEM>()
{

//     TMEM memory;
//     this->SetDefaultMem(memory);
    ff[0]	= 0.; // X component of the body force
	ff[1]	= 0.; // Y component of the body force
	ff[2] = 0.; // Z component of the body force - not used for this class

	fPlaneStress = 0;
}

template <class TMEM>
TPZElasticityMaterialWithMem<TMEM>::TPZElasticityMaterialWithMem(int id): TPZMatWithMem<TMEM>(id)
{
//         TMEM memory;
//     this->SetDefaultMem(memory);
    ff[0]	= 0.; // X component of the body force
	ff[1]	= 0.; // Y component of the body force
	ff[2] = 0.; // Z component of the body force - not used for this class

	fPlaneStress = 0;
}

template <class TMEM>
TPZElasticityMaterialWithMem<TMEM>::TPZElasticityMaterialWithMem(int id, REAL fx, REAL fy, int planestress): TPZMatWithMem<TMEM>(id)
{

//     TMEM memory;
//     this->SetDefaultMem(memory);

	ff[0]	= fx; // X component of the body force
	ff[1]	= fy; // Y component of the body force
	ff[2] = 0.; // Z component of the body force - not used for this class
	fPlaneStress = planestress;

}

template <class TMEM>
int TPZElasticityMaterialWithMem<TMEM>::NStateVariables() {
	return 2;
}

template <class TMEM>
void TPZElasticityMaterialWithMem<TMEM>::Print(std::ostream &out) {

	out << "\tF   = " << ff[0] << ' ' << ff[1]   << std::endl;

}

template <class TMEM>
void TPZElasticityMaterialWithMem<TMEM>::Contribute(TPZMaterialData &data, REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef)
{

	TPZFMatrix<REAL> &dphi = data.dphix;
	TPZFMatrix<REAL> &phi = data.phi;
	TPZFMatrix<REAL> &axes=data.axes;

	int phc,phr,dphc,dphr,efr,efc,ekr,ekc;
	phc = phi.Cols();
	phr = phi.Rows();
	dphc = dphi.Cols();
	dphr = dphi.Rows();
	efr = ef.Rows();
	efc = ef.Cols();
	ekr = ek.Rows();
	ekc = ek.Cols();
	if(phc != 1 || dphr != 2 || phr != dphc){
		PZError << "\nTPZElasticityMaterial.contr, inconsistent input data : \n" <<
		"phi.Cols() = " << phi.Cols() << " dphi.Cols() = " << dphi.Cols() <<
		" phi.Rows = " << phi.Rows() << " dphi.Rows = " <<
		dphi.Rows() << "\nek.Rows() = " << ek.Rows() << " ek.Cols() = "
	    << ek.Cols() <<
		"\nef.Rows() = " << ef.Rows() << " ef.Cols() = "
	    << ef.Cols() << "\n";
		return;
		//		PZError.show();
	}

    int index = data.intGlobPtIndex;


    TMEM localmem;

    localmem = this->MemItem(index);


	REAL nu=localmem.fnu;
    REAL E= localmem.fE;

	REAL fEover1MinNu2 = E/(1-nu*nu);  //G = E/2(1-nu);
	REAL fEover21PlusNu = E/(2.*(1+nu));//E/(1-nu)

    if (localmem.fE<1.e-6)
    {
        std::cout << "modulo de elasticidade nulo"<<std::endl;
        DebugStop(); //deve inicializar o qsi pelo SetDefaultMemory
    }


	TPZFNMatrix<4,STATE> du(2,2);
	/*
	 * Plain strain materials values
	 */
	REAL nu1 = 1. - nu;//(1-nu)
	REAL nu2 = (1.-2.*nu)/2.;
	REAL F = E/((1.+nu)*(1.-2.*nu));

for( int in = 0; in < phr; in++ ) {
		du(0,0) = dphi(0,in)*axes(0,0)+dphi(1,in)*axes(1,0);//dvx
		du(1,0) = dphi(0,in)*axes(0,1)+dphi(1,in)*axes(1,1);//dvy

        for (int col = 0; col < efc; col++)
        {
					ef(2*in, col) += weight * (ff[0]*phi(in,0) - du(0,0)*0 - du(1,0)*0);  // direcao x
					ef(2*in+1, col) += weight * (ff[1]*phi(in,0) - du(0,0)*0 - du(1,0)*0);// direcao y <<<----
        }
		for( int jn = 0; jn < phr; jn++ ) {
			du(0,1) = dphi(0,jn)*axes(0,0)+dphi(1,jn)*axes(1,0);//dux
			du(1,1) = dphi(0,jn)*axes(0,1)+dphi(1,jn)*axes(1,1);//duy


			if (fPlaneStress != 1){
				/* Plane Strain State */
				ek(2*in,2*jn) += weight * (
										   nu1 * du(0,0)*du(0,1)+ nu2 * du(1,0)*du(1,1)
										   ) * F;

				ek(2*in,2*jn+1) += weight * (
											 nu*du(0,0)*du(1,1)+ nu2*du(1,0)*du(0,1)
											 ) * F;

				ek(2*in+1,2*jn) += weight * (
											 nu*du(1,0)*du(0,1)+ nu2*du(0,0)*du(1,1)
											 ) * F;

				ek(2*in+1,2*jn+1) += weight * (
											   nu1*du(1,0)*du(1,1)+ nu2*du(0,0)*du(0,1)
											   ) * F;
			}
			else{
				/* Plain stress state */
				ek(2*in,2*jn) += weight * (
										   fEover1MinNu2 * du(0,0)*du(0,1)+ fEover21PlusNu * du(1,0)*du(1,1)
										   );

				ek(2*in,2*jn+1) += weight * (
											 fEover1MinNu2*nu*du(0,0)*du(1,1)+ fEover21PlusNu*du(1,0)*du(0,1)
											 );

				ek(2*in+1,2*jn) += weight * (
											 fEover1MinNu2*nu*du(1,0)*du(0,1)+ fEover21PlusNu*du(0,0)*du(1,1)
											 );

				ek(2*in+1,2*jn+1) += weight * (
											   fEover1MinNu2*du(1,0)*du(1,1)+ fEover21PlusNu*du(0,0)*du(0,1)
											   );
			}
		}
	}

}

template <class TMEM>
void TPZElasticityMaterialWithMem<TMEM>::ContributeBC(TPZMaterialData &data,REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
{

	TPZFMatrix<REAL> &phi = data.phi;
     int dim = Dimension();

	const REAL BIGNUMBER  = TPZMaterial::gBigNumber;

	int phr = phi.Rows();
	short in,jn;

    if (ef.Cols() != bc.NumLoadCases()) {
        DebugStop();
    }

//		In general when the problem is  needed to stablish any convention for ContributeBC implementations

//     REAL v2[2];
// 	v2[0] = bc.Val2()(0,0);
// 	v2[1] = bc.Val2()(1,0);
    int nstate = NStateVariables();

    TPZFMatrix<STATE> &v1 = bc.Val1();


    switch (bc.Type()) {
        case 0 :			// Dirichlet condition
        {
            for(in = 0 ; in < phr; in++) {
                for (int il = 0; il<0; il++)
                {
                    REAL v2[2];
                    v2[0] = bc.Val2(il)(0,0);
                    v2[1] = bc.Val2(il)(1,0);
                    ef(2*in,il)   += BIGNUMBER * v2[0] * phi(in,0) * weight;        // forced v2 displacement
                    ef(2*in+1,il) += BIGNUMBER * v2[1] * phi(in,0) * weight;        // forced v2 displacement
                }
                for (jn = 0 ; jn < phi.Rows(); jn++)
                {
                    ek(2*in,2*jn)     += BIGNUMBER * phi(in,0) *phi(jn,0) * weight;
                    ek(2*in+1,2*jn+1) += BIGNUMBER * phi(in,0) *phi(jn,0) * weight;
                }
            }
        }
            break;

        case 1 :		// Neumann condition
        {
            for (in = 0; in < phr; in++)
            {
                for (int il = 0; il <0; il++)
                {
                    TPZFNMatrix<2,STATE> v2 = bc.Val2(il);
                    ef(2*in,il) += v2(0,0) * phi(in,0) * weight;        // force in x direction
                    ef(2*in+1,il) +=  v2(1,0) * phi(in,0) * weight;      // force in y direction
                }
            }
        }
            break;

        case 2 :		// Mixed Condition
        {
            for(in = 0 ; in < phi.Rows(); in++)
            {
                for (int il = 0; il <0; il++)
                {
                    TPZFNMatrix<2,STATE> v2 = bc.Val2(il);
                    ef(2*in,il) += v2(0,0) * phi(in,0) * weight;        // force in x direction
                    ef(2*in+1,il) += v2(1,0) * phi(in,0) * weight;      // forced in y direction
                }

                for (jn = 0 ; jn < phi.Rows(); jn++) {
                    ek(2*in,2*jn) += bc.Val1()(0,0) * phi(in,0) *
                    phi(jn,0) * weight;         // peso de contorno => integral de contorno
                    ek(2*in+1,2*jn) += bc.Val1()(1,0) * phi(in,0) *
                    phi(jn,0) * weight;
                    ek(2*in+1,2*jn+1) += bc.Val1()(1,1) * phi(in,0) *
                    phi(jn,0) * weight;
                    ek(2*in,2*jn+1) += bc.Val1()(0,1) * phi(in,0) *
                    phi(jn,0) * weight;
                }
            }   // este caso pode reproduzir o caso 0 quando o deslocamento


        case 3: // Directional Null Dirichlet - displacement is set to null in the non-null vector component direction

    for(in = 0 ; in < phr; in++) {
      //ef(2*in+0,0) += BIGNUMBER *   bc.Val2()(0,0) * phi(in,0) * weight;
      //ef(2*in+1,0) += BIGNUMBER *   bc.Val2()(1,0) * phi(in,0) * weight;
      for (jn = 0 ; jn < phr; jn++) {
        ek(2*in+0,2*jn+0) += BIGNUMBER * phi(in,0) * phi(jn,0) * weight * bc.Val2()(0,0);
        ek(2*in+1,2*jn+1) += BIGNUMBER * phi(in,0) * phi(jn,0) * weight * bc.Val2()(1,0);
      }//jn
    }//in
            break;


        case 4: // stressField Neumann condition
            {
                REAL v2[2];
                for(in = 0; in < dim; in ++)
                {
                    v2[in] =  ( v1(in,0) * data.normal[0] +
                                v1(in,1) * data.normal[1]);
                }
                // The normal vector points towards the neighbour. The negative sign is there to
                // reflect the outward normal vector.
                for(in = 0 ; in < phi.Rows(); in++) {
                    ef(nstate*in+0,0) += v2[0] * phi(in,0) * weight;
                    ef(nstate*in+1,0) += v2[1] * phi(in,0) * weight;
                    //	cout << "normal:" << data.normal[0] << ' ' << data.normal[1] << ' ' << data.normal[2] << endl;
                    //	cout << "val2:  " << v2[0]  << endl;
                }
            }
            break;

        case 5://PRESSAO DEVE SER POSTA NA POSICAO 0 DO VETOR v2
            {
                TPZFNMatrix<2,STATE> res(2,1,0.);
                for(in = 0 ; in < phi.Rows(); in++)
                {
                    for (int il=0; il<0; il++)
                    {
                        ef(nstate*in+0,0) += (bc.Val2(il)(0,0)*data.normal[0]) * phi(in,0) * weight ;
                        ef(nstate*in+1,0) += (bc.Val2(il)(0,0)*data.normal[1]) * phi(in,0) * weight ;
                    }
                    for(jn=0; jn<phi.Rows(); jn++)
                    {
                        for(int idf=0; idf<2; idf++) for(int jdf=0; jdf<2; jdf++)
                        {
                            ek(nstate*in+idf,nstate*jn+jdf) += bc.Val1()(idf,jdf)*data.normal[idf]*data.normal[jdf]*phi(in,0)*phi(jn,0)*weight;
                            //BUG FALTA COLOCAR VAL2
                            //                        DebugStop();
                        }
                    }

                }
            }
            break;

        case 6://PRESSAO DEVE SER POSTA NA POSICAO 0 DO VETOR v2
            {
                TPZFNMatrix<2,STATE> res(2,1,0.);
                for(in = 0 ; in < phi.Rows(); in++)
                {
                    for (int il=0; il<0; il++)
                    {
                        ef(nstate*in+0,0) += (bc.Val2(il)(0,0)*data.normal[0]) * phi(in,0) * weight ;
                        ef(nstate*in+1,0) += (bc.Val2(il)(0,0)*data.normal[1]) * phi(in,0) * weight ;
                    }
                    for(jn=0; jn<phi.Rows(); jn++)
                    {
                        for(int idf=0; idf<2; idf++) for(int jdf=0; jdf<2; jdf++)
                        {
                            ek(nstate*in+idf,nstate*jn+jdf) += bc.Val1()(idf,jdf)*phi(in,0)*phi(jn,0)*weight;
                            //BUG FALTA COLOCAR VAL2
                            //                        DebugStop();
                        }
                    }

                }

            }
            break;

        }      // �nulo introduzindo o BIGNUMBER pelos valores da condi�o
    } // 1 Val1 : a leitura �00 01 10 11
}



template <class TMEM>
int TPZElasticityMaterialWithMem<TMEM>::VariableIndex(const std::string &name){


	if(!strcmp("displacement",name.c_str()))     return 9;
	if(!strcmp("Displacement",name.c_str()))     return 9;
	if(!strcmp("DisplacementMem",name.c_str()))     return 9;
	if(!strcmp("Pressure",name.c_str()))         return 1;
	if(!strcmp("MaxStress",name.c_str()))        return 2;
	if(!strcmp("PrincipalStress1",name.c_str())) return 3;
	if(!strcmp("PrincipalStress2",name.c_str())) return 4;
	if(!strcmp("SigmaX",name.c_str()))           return 5;
	if(!strcmp("SigmaY",name.c_str()))           return 6;
	if(!strcmp("TauXY",name.c_str()))            return 8;//Cedric
	if(!strcmp("Strain",name.c_str()))           return 11;//Philippe
	if(!strcmp("SigmaZ",name.c_str()))           return 12;//Philippe

	if(!strcmp("sig_x",name.c_str()))            return 5;
	if(!strcmp("sig_y",name.c_str()))            return 6;
	if(!strcmp("tau_xy",name.c_str()))           return 8;//Cedric
	if(!strcmp("Displacement6",name.c_str()))    return 7;
	if(!strcmp("Stress",name.c_str()))           return 10;
	if(!strcmp("Flux",name.c_str()))           return 10;
    if(!strcmp("J2",name.c_str()))           return 20;
    if(!strcmp("I1",name.c_str()))           return 21;
    if(!strcmp("J2Stress",name.c_str()))           return 20;
    if(!strcmp("I1Stress",name.c_str()))           return 21;
    if(!strcmp("Alpha",name.c_str()))        return 22;
    if(!strcmp("PlasticSqJ2",name.c_str()))        return 22;
    if(!strcmp("PlasticSqJ2El",name.c_str()))        return 22;
    if(!strcmp("YieldSurface",name.c_str()))        return 27;
    if(!strcmp("NormalStress",name.c_str()))        return 23;
    if(!strcmp("ShearStress",name.c_str()))        return 24;
    if(!strcmp("NormalStrain",name.c_str()))        return 25;
    if(!strcmp("ShearStrain",name.c_str()))        return 26;



	//   cout << "TPZElasticityMaterial::VariableIndex Error\n";
	return TPZMaterial::VariableIndex(name);
}

template <class TMEM>
int TPZElasticityMaterialWithMem<TMEM>::NSolutionVariables(int var){

	switch(var) {
		case 0:
			return 2;
		case 1:
		case 2:
			return 1;
		case 3:
		case 4:
			return 2;
		case 5:
		case 6:
		case 8:
			return 1;
		case 7:
			return 6;
		case 9:
			return 3;
		case 10 : //Stress Tensor
			return 3;
        case 11 : //Strain Tensor
            return 3;
            // SigZ
        case 12:
            return 1;
        case 20:
            return 1;
        case 21:
            return 1;
        case 22:
            return 1;
        case 23:
        case 24:
        case 25:
        case 26:
        case 27:
            return 3;
		default:
			return TPZMaterial::NSolutionVariables(var);
	}
}

template<class TMEM>
void TPZElasticityMaterialWithMem<TMEM>::Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout)
{
    int numbersol = data.dsol.size();

    if (numbersol != 1) {
        std::cout << "Solution is empty"<<std::endl;
        DebugStop();
    }
    int ipos = 0;

    TPZVec<STATE> &Sol = data.sol[ipos];
    TPZFMatrix<STATE> &DSol = data.dsol[ipos];
    TPZFMatrix<REAL> &axes = data.axes;
    TPZFNMatrix<4,STATE> DSolxy(2,2);

	REAL epsx;
	REAL epsy;
	REAL epsxy;
    REAL epsz = 0.;
	REAL SigX;
	REAL SigY;
    REAL SigZ;
	REAL TauXY,aux,Sig1,Sig2,angle;

    // dudx - dudy
	DSolxy(0,0) = DSol(0,0)*axes(0,0)+DSol(1,0)*axes(1,0);
	DSolxy(1,0) = DSol(0,0)*axes(0,1)+DSol(1,0)*axes(1,1);
	// dvdx - dvdy
	DSolxy(0,1) = DSol(0,1)*axes(0,0)+DSol(1,1)*axes(1,0);
	DSolxy(1,1) = DSol(0,1)*axes(0,1)+DSol(1,1)*axes(1,1);

    epsx = DSolxy(0,0);// du/dx
    epsy = DSolxy(1,1);// dv/dy
    epsxy = 0.5*(DSolxy(1,0)+DSolxy(0,1));

    TMEM localmem;

    int index = data.intGlobPtIndex;
    localmem = this->MemItem(index);


	REAL nu=localmem.fnu;
    REAL E= localmem.fE;

	REAL fEover1MinNu2 = E/(1-nu*nu);  //G = E/2(1-nu);
	REAL fEover21PlusNu = E/(2.*(1+nu));//E/(1-nu)

    if (localmem.fE<1.e-6)
    {
        std::cout << "modulo de elasticidade nulo"<<std::endl;
        DebugStop(); //deve inicializar o qsi pelo SetDefaultMemory
    }
    REAL lambda=(nu*E)/((1.+nu)*(1.-2.*nu)) ;
    REAL mu =E/(2.*(1.+nu));

    if (this->fPlaneStress == 1) {

        epsz = -lambda*(epsx+epsy)/(lambda+2.*mu);
    }
    else {
        epsz = 0.;
    }
    TauXY = 2*mu*epsxy;
#ifdef DEBUG
    REAL TauXY2 = E*epsxy/(1.+nu);
    #ifdef REALfloat
    if (fabs(TauXY-TauXY2) > 1.e-10) {
        DebugStop();
    }
    #else
    if (fabs(TauXY-TauXY2) > 1.e-6) {
        DebugStop();
    }
	#endif
#endif
    if (this->fPlaneStress == 1){
        SigX = fEover1MinNu2*(epsx+nu*epsy);
        SigY = fEover1MinNu2*(nu*epsx+epsy);
        SigZ = 0;
    }
    else
    {
        SigX = E/((1.-2.*nu)*(1.+nu))*((1.-nu)*epsx+nu*epsy);
        SigY = E/((1.-2.*nu)*(1.+nu))*(nu*epsx+(1.-nu)*epsy);
        SigZ = 0+lambda*(epsx+epsy);
    }

	switch(var) {
		case 0:
			//numvar = 2;
			Solout[0] = Sol[0];
			Solout[1] = Sol[1];
			break;
		case 7:
			//numvar = 6;
			Solout[0] = Sol[0];
			Solout[1] = Sol[1];
			Solout[2] = 0.;
			Solout[3] = 0.;
			Solout[4] = 0.;
			Solout[5] = 0.;
			break;
		case 1:
		case 2:
		case 3:
		case 4:
		case 5:
		case 6:
		case 8:
		case 10:

			//numvar = 1;
			Solout[0] = SigX+SigY+SigZ;
            // Pressure variable
			if(var == 1) {
				Solout[0] = SigX+SigY+SigZ;
				return;
			}
            // TauXY variable
			if(var == 8) {
				Solout[0] = TauXY;
				return;
			}
			if(var ==5) {
				Solout[0] = SigX;
				return;
			}
			if(var == 6) {
				Solout[0] = SigY;
				return;
			}
			aux = sqrt(0.25*(SigX-SigY)*(SigX-SigY)
					   +(TauXY)*(TauXY));
			// Philippe 13/5/99
			//         if(abs(Tau) < 1.e-10 && abs(SigY-SigX) < 1.e-10) angle = 0.;
			if(fabs(TauXY) < 1.e-10 && fabs(SigY-SigX) < 1.e-10) angle = 0.;
			else angle = atan2(2*TauXY,SigY-SigX)/2.;
			Sig1 = 0.5*(SigX+SigY)+aux;
			Sig2 = 0.5*(SigX+SigY)-aux;
			if(var == 3 ){
				//numvar = 2;
				Solout[0] = Sig1*cos(angle);
				Solout[1] = Sig1*sin(angle);
				return;
			}
			if(var == 4 ) {
				//numvar = 2;
				Solout[0] = -Sig2*sin(angle);
				Solout[1] = Sig2*cos(angle);
				return;
			}
			if(var == 2) {
				REAL sigmax;
				sigmax = (fabs(Sig1) < fabs(Sig2))? fabs(Sig2) : fabs(Sig1);
				Solout[0] = sigmax;
				return;
			}
			if (var ==10)
			{
				Solout[0] = SigX;
				Solout[1] = SigY;
				Solout[2] = TauXY;
				return;
			}
			std::cout << "Very critical error TPZElasticityMaterial::Solution\n";
			exit(-1);
			//         Solout[0] /= 0.;
			break;
		case 9:
			Solout[0] = Sol[0];
			Solout[1] = Sol[1];
			Solout[2] = 0.;
			break;
        case 11:
            Solout[0] = epsx;
			Solout[1] = epsy;
			Solout[2] = epsxy;
            break;
        case 12:
            Solout[0] = SigZ;
            break;

        case 20:
        {

           REAL J2 = (pow(SigX + SigY,2) - (3*(-pow(SigX,2) - pow(SigY,2) + pow(SigX + SigY,2) - 2*pow(TauXY,2)))/2.)/2.;

            Solout[0]=J2;
            break;
        }
        case 21:
        {
            REAL I1 = SigX+SigY;
            Solout[0]=I1;
            break;
        }
        case 22:
            Solout[0] = 0.;
            break;
        case 23:
            // normal stress
            Solout[0] = SigX;
            Solout[1] = SigY;
            Solout[2] = SigZ;
            break;
        case 24:
            // shear stress
            Solout[0] = TauXY;
            Solout[1] = 0.;
            Solout[2] = 0.;
            break;
        case 25:
            Solout[0] = epsx;
            Solout[1] = epsy;
            Solout[2] = epsz;
            break;
        case 26:
            Solout[0] = epsxy;
            Solout[1] = 0.;
            Solout[2] = 0.;
            break;
        case 27:
            Solout[0] = 0.;
            Solout[1] = 0.;
            Solout[2] = 0.;
            break;
		default:
			TPZMaterial::Solution(Sol,DSol,axes,var,Solout);
			break;
	}
}



template <class TMEM>
void TPZElasticityMaterialWithMem<TMEM>::FillDataRequirements(TPZMaterialData &data)
{
    data.fNeedsSol = true;
    data.fNeedsNormal = false;

}

template <class TMEM>
void TPZElasticityMaterialWithMem<TMEM>::FillBoundaryConditionDataRequirement(int type,TPZMaterialData &data)
{
    data.fNeedsSol = false;
    data.fNeedsNormal = false;
    if (type == 4 || type == 5 || type == 6) {
        data.fNeedsNormal = true;
    }
}

template class TPZElasticityMaterialWithMem<TPZElasticMem>;
