//
// Created by Gustavo Batistela on 5/13/21.
//

#include "TPZIsotropicPermeability.h"


void TPZIsotropicPermeability::SetConstantPermeability( STATE constant) {
    fConstantPermeability = constant;
}

void TPZIsotropicPermeability::SetPermeabilityFunction(PermeabilityFunctionType &perm_function) {
    fPermeabilityFunction = perm_function;
}

STATE TPZIsotropicPermeability::GetPermeability( TPZVec<REAL> &coord) {
    return fPermeabilityFunction ? fPermeabilityFunction(coord) : fConstantPermeability;
}

int TPZIsotropicPermeability::ClassId() const {
    return -9999991;
}
