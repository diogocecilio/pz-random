#ifndef ELASMAT_WM
#define ELASMAT_WM

#include <iostream>
#include "pzmatwithmem.h"
#include "pzelasticmem.h"

template<class TMEM = TPZElasticMem>
class TPZElasticityMaterialWithMem : public TPZMatWithMem<TMEM>{

public:

	TPZElasticityMaterialWithMem();

	TPZElasticityMaterialWithMem(int id);

	TPZElasticityMaterialWithMem(int id, REAL fx, REAL fy, int planestress=1);

	virtual void Contribute(TPZMaterialData &data, REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef);

	virtual void ContributeBC(TPZMaterialData &data,REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);

	virtual void FillDataRequirements(TPZMaterialData &data);

    virtual void FillBoundaryConditionDataRequirement(int type, TPZMaterialData &data);

	void UpdateMem(TPZMaterialData &data, REAL E, REAL nu);

	int Dimension() const { return 2;}

	virtual  int NStateVariables();

	void Print(std::ostream &out);

	int VariableIndex(const std::string &name);

	int NSolutionVariables(int var);

	void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout);

	REAL ff[3];

	int fPlaneStress;

};

#endif
