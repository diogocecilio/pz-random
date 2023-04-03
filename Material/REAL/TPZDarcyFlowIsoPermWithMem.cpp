//
// Created by Gustavo Batistela on 5/13/21.
//

#include "TPZDarcyFlowIsoPerm.h"
#include "pzaxestools.h"
#include <pzbndcond.h>

//class TPZIsotropicPermeabilityBC;

TPZDarcyFlowIsoPerm::TPZDarcyFlowIsoPerm() :  fDim(-1) {}

TPZDarcyFlowIsoPerm::TPZDarcyFlowIsoPerm(int id, int dim) : fDim(dim) {}

void TPZDarcyFlowIsoPerm::SetDimension(int dim) {
    if (dim > 3 || dim < 1) DebugStop();
    fDim = dim;
}

void TPZDarcyFlowIsoPerm::Contribute( TPZMaterialData &data, STATE weight, TPZFMatrix<STATE> &ek,
                              TPZFMatrix<STATE> &ef) {

    const TPZFMatrix<REAL> &phi = data.phi;
    const TPZFMatrix<REAL> &dphi = data.dphix;
    const TPZVec<REAL> &x = data.x;
    const TPZFMatrix<REAL> &axes = data.axes;
    const TPZFMatrix<REAL> &jacinv = data.jacinv;
    auto phr = dphi.Cols();


	//STATE perm = GetPermeability(data.x);

	//std:: cout << "perm "<< perm<<std::endl;
	TPZManVector<REAL,3> perm = fConstPermeability;

	
    STATE source_term = 0;
    if (this->HasForcingFunction()) {
        TPZManVector<STATE, 1> res(1);
        ForcingFunction()->Execute(x, res);
        source_term = -res[0];
    }
	

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

void TPZDarcyFlowIsoPerm::ContributeBC( TPZMaterialData &data, STATE weight, TPZFMatrix<STATE> &ek,
                                TPZFMatrix<STATE> &ef, TPZBndCond &bc) {

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
                ef(in, 0) += (STATE) (gBigNumber * phi.Get(in, 0) * weight) * v2;
                for (jn = 0; jn < phr; jn++) {
                    ek(in, jn) += gBigNumber * phi.Get(in, 0) * phi.Get(jn, 0) * weight;
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

int TPZDarcyFlowIsoPerm::VariableIndex(const std::string &name)  {

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

    return TPZMaterial::VariableIndex(name);
}

int TPZDarcyFlowIsoPerm::NSolutionVariables(int var)  {

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


    return TPZMaterial::NSolutionVariables(var);
}

void TPZDarcyFlowIsoPerm::Solution( TPZMaterialData &data, int var, TPZVec<STATE> &solOut) {

    switch (var) {
        case 1: {
            solOut.resize(1);
            // Solution/Pressure
            solOut[0] = data.sol[0][0];
            return;
        }
        case 2: {
            // Derivative/GradU
            TPZFNMatrix<9, STATE> dsoldx;
            TPZAxesTools<STATE>::Axes2XYZ(data.dsol[0], dsoldx, data.axes);
            for (int id = 0; id < fDim; id++) {
                solOut[id] = dsoldx(id, 0);
            }
            return;
        }
        case 3: {
            // KDuDx;
            TPZFNMatrix<9, STATE> dsoldx;
            TPZAxesTools<STATE>::Axes2XYZ(data.dsol[0], dsoldx, data.axes);
            //const STATE perm = GetPermeability(data.x);
			TPZManVector<REAL,3>  perm = fConstPermeability;
            solOut[0] = perm[0] * dsoldx(0, 0);
            return;
        }
        case 4: {
            // KDuDy;
            TPZFNMatrix<9, STATE> dsoldx;
            TPZAxesTools<STATE>::Axes2XYZ(data.dsol[0], dsoldx, data.axes);
            //const STATE perm = GetPermeability(data.x);
			TPZManVector<REAL,3>  perm = fConstPermeability;
            solOut[0] = perm[1] * dsoldx(1, 0);
            return;
        }
        case 5: {
            // KDuDz;
            TPZFNMatrix<9, STATE> dsoldx;
            TPZAxesTools<STATE>::Axes2XYZ(data.dsol[0], dsoldx, data.axes);
            //const STATE perm = GetPermeability(data.x);
			TPZManVector<REAL,3>  perm = fConstPermeability;
            solOut[0] = perm[2] * dsoldx(2, 0);
            return;
        }
        case 6: {
            // NormKDu;
            TPZFNMatrix<9, STATE> dsoldx;
            TPZAxesTools<STATE>::Axes2XYZ(data.dsol[0], dsoldx, data.axes);

            //const STATE perm = GetPermeability(data.x);
			TPZManVector<REAL,3> perm = fConstPermeability;
            STATE res = 0;
            for (int id = 0; id < fDim; id++) {
                res += perm[id] * dsoldx(id, 0) * perm[id] * dsoldx(id, 0);
            }
            solOut[0] = sqrt(res);
            return;
        }
        case 7: {
            // MinusKGradU/Flux;
            TPZFNMatrix<9, STATE> dsoldx;
            TPZAxesTools<STATE>::Axes2XYZ(data.dsol[0], dsoldx, data.axes);
            //const STATE perm = GetPermeability(data.x);
			TPZManVector<REAL,3> perm = fConstPermeability;
            solOut.resize(fDim);
            for (int id = 0; id < fDim; id++) {
                solOut[id] =  perm[id] * dsoldx(id, 0);
            }
            return;
        }
        case 8: {
            // POrder
            solOut[0] = data.p;
            return;
        }
        case 9: {
            // ExactPressure/ExactSolution
            TPZVec<STATE> exact_pressure(1);
            TPZFMatrix<STATE> exact_flux(fDim, 1);
            //fExactSol(data.x, exact_pressure, exact_flux);
            solOut[0] = exact_pressure[0];
            return;
        }
        case 10: {
            // ExactFlux
            TPZVec<STATE> exact_pressure(1);
            TPZFMatrix<STATE> exact_flux(fDim, 1);
            //fExactSol(data.x, exact_pressure, exact_flux);
            for (int id = 0; id < fDim; id++) {
                solOut[id] = exact_flux[id];
            }
            return;
        }
        case 11: {

            solOut[0] = 0;
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
            solOut[0] = res;

            return;
        }

    default: {
            PZError << __PRETTY_FUNCTION__ << "\n Post-processing variable index not implemented!\n";
            DebugStop();
        }
    }
}

void TPZDarcyFlowIsoPerm::GetSolDimensions(uint64_t &u_len, uint64_t &du_row, uint64_t &du_col) const {
    u_len=1;
    du_row=fDim;
    du_col=1;
}

void TPZDarcyFlowIsoPerm::Errors( TPZMaterialData &data,
                          TPZVec<REAL> &errors) {
    const TPZVec<REAL> &x = data.x;
    const TPZVec<STATE> &sol = data.sol[0];
    const TPZFMatrix<STATE> &dsol = data.dsol[0];
    const TPZFMatrix<REAL> &axes = data.axes;

#ifdef PZDEBUG
    if(!this->HasExactSol()){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"\nThe material has no associated exact solution. Aborting...\n";
        DebugStop();
    }
#endif
    
    errors.Resize(NEvalErrors(), 0.);

    TPZVec<STATE> exact_pressure(1, 0);
    TPZFMatrix<STATE> exact_flux(fDim, 1, 0);
    //fExactSol(x, exact_pressure, exact_flux);

    TPZFNMatrix<3,STATE> gradu(3,1);
    TPZAxesTools<STATE>::Axes2XYZ(dsol,gradu,axes);

    // errors[1] - L2 norm error
    REAL diff = fabs(sol[0] - exact_pressure[0]);
    errors[1] = diff * diff;

    // errors[2] - H1 semi-norm: |H1| = K*(grad[u] - grad[u_exact])
     //STATE perm = GetPermeability(data.x);
	TPZManVector<REAL,3>  perm = fConstPermeability;

    TPZVec<REAL> graduDiff(fDim, 0);
    for (int id = 0; id < fDim; id++) {
        graduDiff[id] += fabs(gradu(id) - exact_flux(id, 0));
    }
    diff = 0;
    for (int id = 0; id < fDim; id++) {
        diff += perm[id] * graduDiff[id];
    }
    errors[2] += abs(diff * diff);

    // errors[0] - H1/Energy norm
    errors[0] = errors[1] + errors[2];

    // TODO confirm with Phil is the following norms are correct
    // errors[3] - L2 norm of the x-component of the flux
    // errors[4] - L2 norm of the y-component of the flux, if applicable
    // errors[5] - L2 norm of the z-component of the flux, if applicable
    TPZFNMatrix<9, STATE> dsoldx;
    TPZAxesTools<STATE>::Axes2XYZ(dsol, dsoldx, axes);
    TPZManVector<STATE, 3> flux_sol(fDim, 0);
    for (int id = 0; id < fDim; id++) {
        flux_sol[id] = - perm[id] * dsoldx(id, 0);
    }

    for (int id = 0; id < fDim; id++) {
        diff = fabs(exact_flux[id] - flux_sol[id]);
        errors[3 + id] = diff * diff;
    }
}

int TPZDarcyFlowIsoPerm::ClassId() const {
    return -999993;
}

// TPZMaterial *TPZDarcyFlow::NewMaterial() const {
//     return new TPZDarcyFlow(*this);
// }

void TPZDarcyFlowIsoPerm::Print(std::ostream &out) const {
    out << "Material Name: " << this->Name() << "\n";
    out << "Material Id: " << this->Id() << "\n";
    out << "Dimension: " << this->Dimension() << "\n\n";
}

void TPZDarcyFlowIsoPerm::FillDataRequirements(TPZMaterialData &data) const {
    data.SetAllRequirements(false);
}


void TPZDarcyFlowIsoPerm::FillBoundaryConditionDataRequirements(int type, TPZMaterialData &data) const {

    data.SetAllRequirements(false);
    if (type == 50) {
        data.fNeedsSol = true;
    }
    if (type == 3) {
        data.fNeedsNormal = true;
    }
}
