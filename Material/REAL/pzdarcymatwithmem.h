#ifndef TPZDarcyMatWithMem_WM
#define TPZDarcyMatWithMem_WM

#include <iostream>
#include "pzmatwithmem.h"
#include "pzdarcymatwithmem.h"
#include "pzdarcymem.h"
template<class TMEM = TPZDarcyMem>
class TPZDarcyMatWithMem : public TPZMatWithMem<TMEM>{

public:

	//TPZDarcyMatWithMem();

	//TPZDarcyMatWithMem(int id);

	TPZDarcyMatWithMem(int id, int dim);

	virtual std::string Name()   const{ return "TPZDarcyMatmem"; }

	virtual int Dimension() const override { return this->fDim; }

	virtual int NStateVariables(){ return 1;}

	virtual void SetDimension(int dim);

	virtual void Contribute(TPZMaterialData &data, REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef);

	virtual void ContributeBC(TPZMaterialData &data,REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);

	int VariableIndex(const std::string &name);

	int NSolutionVariables(int var);

	void Solution( TPZMaterialData &data, int var, TPZVec<STATE> &solOut) ;

	virtual void FillDataRequirements(TPZMaterialData &data);

    virtual void FillBoundaryConditionDataRequirement(int type, TPZMaterialData &data);

	void UpdateMem(TPZMaterialData &data, TPZManVector<REAL,3> permeability);

	void Print(std::ostream &out);

	int fDim;

};

#endif
