
#include "pzdarcymatwithmem.h"
#include "pzbndcond.h"
#include "pzaxestools.h"
#include <pzbndcond.h>

// template <class TMEM>
// TPZDarcyMatWithMem<TMEM>::TPZDarcyMatWithMem(): TPZMatWithMem<TMEM>(-1),fDim(-1)
// {
//
//
// }
//
// template <class TMEM>
// TPZDarcyMatWithMem<TMEM>::TPZDarcyMatWithMem(int id): TPZMatWithMem<TMEM>(id),fDim(-1)
// {
//
// }

template <class TMEM>
TPZDarcyMatWithMem<TMEM>::TPZDarcyMatWithMem(int id,int dim): TPZMatWithMem<TMEM>(id),fDim(dim)
{

}

template <class TMEM>
void TPZDarcyMatWithMem<TMEM>::Print(std::ostream &out) {
    out << "Material Name: " << this->Name() << "\n";
    out << "Material Id: " << this->Id() << "\n";
    out << "Dimension: " << this->Dimension() << "\n\n";

}

template <class TMEM>
void TPZDarcyMatWithMem<TMEM>::TPZDarcyMatWithMem::SetDimension(int dim) {

    fDim = dim;
}
template <class TMEM>
void TPZDarcyMatWithMem<TMEM>::UpdateMem(TPZMaterialData &data, TPZManVector<REAL,3> permeability)
{
    int index = data.intGlobPtIndex;
	this->MemItem(index).fpermeability=permeability;

}

template <class TMEM>
void TPZDarcyMatWithMem<TMEM>::Contribute(TPZMaterialData &data, REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef)
{
   const TPZFMatrix<REAL> &phi = data.phi;
   const TPZFMatrix<REAL> &dphi = data.dphix;

   auto phr = dphi.Cols();

    int index = data.intGlobPtIndex;

    TMEM localmem;

    localmem = this->MemItem(index);
	//STATE perm = GetPermeability(data.x);
	//std:: cout << "perm "<< perm<<std::endl;
	TPZManVector<REAL,3> perm = localmem.fpermeability;

    STATE source_term = 0;

    // Darcy's equation
    for (int in = 0; in < phr; in++) {
        ef(in, 0) -= weight * source_term * (phi.Get(in, 0));
        for (int jn = 0; jn < phr; jn++) {
            for (int kd = 0; kd < fDim; kd++) {
                ek(in, jn) += weight * (dphi.Get(kd, in) * perm[kd]* dphi.Get(kd, jn));
            }
        }
    }

}

template <class TMEM>
void TPZDarcyMatWithMem<TMEM>::ContributeBC(TPZMaterialData &data,REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
{
  const TPZFMatrix<REAL> &phi = data.phi;
    const TPZFMatrix<REAL> &axes = data.axes;
    int phr = phi.Rows();
    int in, jn;

    STATE v2 = bc.Val2()(0,0);


    if (bc.HasForcingFunction()) {
        TPZManVector<STATE, 1> rhs_val(1);
        TPZFNMatrix<1, STATE> mat_val(1, 1);
        bc.ForcingFunction()->Execute(data.x, rhs_val);
        v2 = rhs_val[0];
    }

    switch (bc.Type()) {
        case 0 : // Dirichlet condition
            for (in = 0; in < phr; in++) {
                //std::cout << "pt "<<data.x << std::endl;
                ef(in, 0) += (STATE) (TPZMaterial::gBigNumber * phi.Get(in, 0) * weight) * v2;
                for (jn = 0; jn < phr; jn++) {
                    ek(in, jn) += TPZMaterial::gBigNumber * phi.Get(in, 0) * phi.Get(jn, 0) * weight;
                }
            }
            break;
        case 1 : // Neumann condition
            for (in = 0; in < phi.Rows(); in++) {
                ef(in, 0) += v2 * (STATE) (phi.Get(in, 0) * weight);
            }
            break;
        case 2 : // Robin condition
            for (in = 0; in < phi.Rows(); in++) {
                ef(in, 0) += v2 * (STATE) (phi.Get(in, 0) * weight);
                for (jn = 0; jn < phi.Rows(); jn++) {
                    ek(in, jn) += bc.Val1()(0, 0) * (STATE) (phi.Get(in, 0) * phi.Get(jn, 0) * weight);
                }
            }
            break;
        default:
            PZError << __PRETTY_FUNCTION__
                    << "\nBoundary condition type not implemented. Please use one of the following:\n"
                    << "\t 0: Dirichlet\n"
                    << "\t 1: Neumann\n"
                    << "\t 2: Robin\n";
            DebugStop();
    }
}



template <class TMEM>
int TPZDarcyMatWithMem<TMEM>::VariableIndex(const std::string &name){
if (!strcmp("Solution", name.c_str())) return 1;
    if (!strcmp("Pressure", name.c_str())) return 1;
    if (!strcmp("Derivative", name.c_str())) return 2;
    if (!strcmp("GradU", name.c_str())) return 2;
    if (!strcmp("KDuDx", name.c_str())) return 3;
    if (!strcmp("KDuDy", name.c_str())) return 4;
    if (!strcmp("KDuDz", name.c_str())) return 5;
    if (!strcmp("NormKDu", name.c_str())) return 6;
    if (!strcmp("MinusKGradU", name.c_str())) return 7;
    if (!strcmp("Flux", name.c_str())) return 7;
    if (!strcmp("POrder", name.c_str())) return 8;
    if (!strcmp("ExactPressure", name.c_str())) return 9;
    if (!strcmp("ExactSolution", name.c_str())) return 9;
    if (!strcmp("ExactFlux", name.c_str())) return 10;
    if (!strcmp("Div", name.c_str())) return 11;
    if (!strcmp("Divergence", name.c_str())) return 11;
    if (!strcmp("ExactDiv", name.c_str())) return 12;
    if (!strcmp("ExactDivergence", name.c_str())) return 12;
    if (!strcmp("FluxL2", name.c_str())) return 13;
    if (!strcmp("Perm", name.c_str())) return 14;
    return TPZMaterial::VariableIndex(name);

}

template <class TMEM>
int TPZDarcyMatWithMem<TMEM>::NSolutionVariables(int var){

    if (var == 1) return 1;      // Solution/Pressure
    if (var == 2) return fDim;   // Derivative/GradU
    if (var == 3) return 1;      // KDuDx;
    if (var == 4) return 1;      // KDuDy;
    if (var == 5) return 1;      // KDuDz;
    if (var == 6) return 1;      // NormKDu;
    if (var == 7) return fDim;   // MinusKGradU/Flux;
    if (var == 8) return 1;      // POrder
    if (var == 9) return 1;      // ExactPressure/ExactSolution
    if (var == 10) return fDim;  // ExactFlux
    if (var == 11) return 1;     // Div/Divergence
    if (var == 12) return 1;     // ExactDiv/ExactDivergence
    if (var == 13) return fDim;  // FluxL2
    if (var == 14) return 1;  // FluxL2


    return TPZMaterial::NSolutionVariables(var);
}

template<class TMEM>
void TPZDarcyMatWithMem<TMEM>::Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout)
{

    TMEM localmem;

    int index = data.intGlobPtIndex;
    localmem = this->MemItem(index);

	TPZManVector<REAL,3> perm=localmem.fpermeability;

 switch (var) {
        case 1: {
            Solout.resize(1);
            // Solution/Pressure
            Solout[0] = data.sol[0][0];
            return;
        }
        case 2: {
            // Derivative/GradU
            TPZFNMatrix<9, STATE> dsoldx;
            TPZAxesTools<STATE>::Axes2XYZ(data.dsol[0], dsoldx, data.axes);
            for (int id = 0; id < fDim; id++) {
                Solout[id] = dsoldx(id, 0);
            }
            return;
        }
        case 3: {
            // KDuDx;
            TPZFNMatrix<9, STATE> dsoldx;
            TPZAxesTools<STATE>::Axes2XYZ(data.dsol[0], dsoldx, data.axes);
            //const STATE perm = GetPermeability(data.x);

            Solout[0] = perm[0] * dsoldx(0, 0);
            return;
        }
        case 4: {
            // KDuDy;
            TPZFNMatrix<9, STATE> dsoldx;
            TPZAxesTools<STATE>::Axes2XYZ(data.dsol[0], dsoldx, data.axes);
            //const STATE perm = GetPermeability(data.x);

            Solout[0] = perm[1] * dsoldx(1, 0);
            return;
        }
        case 5: {
            // KDuDz;
            TPZFNMatrix<9, STATE> dsoldx;
            TPZAxesTools<STATE>::Axes2XYZ(data.dsol[0], dsoldx, data.axes);
            //const STATE perm = GetPermeability(data.x);

            Solout[0] = perm[2] * dsoldx(2, 0);
            return;
        }
        case 6: {
            // NormKDu;
            TPZFNMatrix<9, STATE> dsoldx;
            TPZAxesTools<STATE>::Axes2XYZ(data.dsol[0], dsoldx, data.axes);

            //const STATE perm = GetPermeability(data.x);

            STATE res = 0;
            for (int id = 0; id < fDim; id++) {
                res += perm[id] * dsoldx(id, 0) * perm[id] * dsoldx(id, 0);
            }
            Solout[0] = sqrt(res);
            return;
        }
        case 7: {
            // MinusKGradU/Flux;
            TPZFNMatrix<9, STATE> dsoldx;
            TPZAxesTools<STATE>::Axes2XYZ(data.dsol[0], dsoldx, data.axes);
            //const STATE perm = GetPermeability(data.x);

            Solout.resize(fDim);
            for (int id = 0; id < fDim; id++) {
                Solout[id] =  perm[id] * dsoldx(id, 0);
            }
            return;
        }
        case 8: {
            // POrder
            Solout[0] = data.p;
            return;
        }
        case 9: {
            // ExactPressure/ExactSolution
            TPZVec<STATE> exact_pressure(1);
            TPZFMatrix<STATE> exact_flux(fDim, 1);
            //fExactSol(data.x, exact_pressure, exact_flux);
            Solout[0] = exact_pressure[0];
            return;
        }
        case 10: {
            // ExactFlux
            TPZVec<STATE> exact_pressure(1);
            TPZFMatrix<STATE> exact_flux(fDim, 1);
            //fExactSol(data.x, exact_pressure, exact_flux);
            for (int id = 0; id < fDim; id++) {
                Solout[id] = exact_flux[id];
            }
            return;
        }
        case 11: {

            Solout[0] = 0;
            return;
        }
        case 12: {
            // ExactDiv/ExactDivergence
            TPZVec<STATE> exact_pressure(1);
            TPZFMatrix<STATE> exact_flux(fDim, 1);
            //fExactSol(data.x, exact_pressure, exact_flux);
            STATE res = 0;
            for (int id = 0; id < fDim; id++) {
                res += exact_flux(id, 0);
            }
            Solout[0] = res;

            return;
        }
            case 14: {
            Solout[0] = perm[0];

            return;
        }

    default: {
            PZError << __PRETTY_FUNCTION__ << "\n Post-processing variable index not implemented!\n";
            DebugStop();
        }
    }
}



template <class TMEM>
void TPZDarcyMatWithMem<TMEM>::FillDataRequirements(TPZMaterialData &data)
{
    data.fNeedsSol = true;
    data.fNeedsNormal = false;

}

template <class TMEM>
void TPZDarcyMatWithMem<TMEM>::FillBoundaryConditionDataRequirement(int type,TPZMaterialData &data)
{
    data.fNeedsSol = false;
    data.fNeedsNormal = false;
    if (type == 4 || type == 5 || type == 6) {
        data.fNeedsNormal = true;
    }
}

template class TPZDarcyMatWithMem<TPZDarcyMem>;
