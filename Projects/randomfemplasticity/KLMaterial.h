/**
 * @file pzmaterial.h
 * @brief Header file for abstract class KLMaterial.\n
 * It implements the weak statement of the differential equation within the PZ environment.
 */

#ifndef KLMaterial_H
#define KLMaterial_H

#include "pzreal.h"
#include "pzvec.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzadmchunk.h"
#include "tpzautopointer.h"
#include "pzsave.h"
#include "pzmaterialdata.h"
#include "pzfunction.h"
#include "pzcompel.h"
#include "pzmaterial.h"
#include <iostream>
#include <string>

class TPZBndCond;

class TPZMaterialData;
class TPZIntPoints;

/**
 * @ingroup material
 * @brief This abstract class defines the behaviour which each derived class needs to implement
 */
/**
 * Classes derived from the KLMaterial class implement the weak statement of the differential equation
 * within the PZ environment \n
 * It is noteworthy to observe that this definition does not depend on the definition of the interpolation space \n
 * KLMaterial objects also need to implement the interface for post processing the results
 */
class  KLMaterial : public TPZMaterial
{
	 enum SOLUTIONVARS{ENone = -1,
	  // Strain
	  EVEC = 0,
	  EVEC1 = 1,
	  EVEC2 = 2,
	  EVEC3 = 3,
	  EVEC4 = 4,
	  EVEC5 = 5,
	  EVECSQR =6
};
private:
    int fId;
    
protected:
   

	int fExpasionOrder;
	
	REAL fLx;
	
	REAL fLy;
	
	REAL fLz;
	
	int ftype;
	
	int fDim;
    
public:
    /** @brief Big number to penalization method, used for Dirichlet conditions */
    static REAL gBigNumber;
    
    /** @brief Creates a material object and inserts it in the vector of material pointers of the mesh. */
	/** Upon return vectorindex contains the index of the material object within the vector */
    KLMaterial(int id,REAL Lx,REAL Ly,REAL Lz,int dim, int type,int expansionorder);
	
	KLMaterial(int id);
    
    /** @brief Default constructor */
    KLMaterial();
    
    /** @brief Creates a material object based on the referred object and inserts it in the vector of material pointers of the mesh. */
	/** Upon return vectorindex contains the index of the material object within the vector */
    KLMaterial(const KLMaterial &mat);
    /** @brief Default destructor */
    virtual ~KLMaterial();
    
    /** 
	 * @brief Fill material data parameter with necessary requirements for the
	 * @since April 10, 2007
	 */
	/** 
	 * Contribute method. Here, in base class, all requirements are considered as necessary. 
	 * Each derived class may optimize performance by selecting only the necessary data.
     */
    virtual void FillDataRequirements(TPZMaterialData &data)
	{
		    data.fNeedsSol = true;
			data.fNeedsNormal = false;
	}
	
	/** 
	 * @brief Fill material data parameter with necessary requirements for the
	 * Contribute method. Here, in base class, all requirements are considered as necessary. 
	 * Each derived class may optimize performance by selecting only the necessary data.
     */
	virtual void FillDataRequirements(TPZVec<TPZMaterialData > &datavec)
	{
		FillDataRequirements( datavec);
	}
    
    /** @brief This method defines which parameters need to be initialized in order to compute the contribution of the boundary condition */
    virtual void FillBoundaryConditionDataRequirement(int type,TPZMaterialData &data)
    {
        // default is no specific data requirements
        if(type == 50)
        {
            data.fNeedsSol = true;
        }
    }
    
    /** @brief This method defines which parameters need to be initialized in order to compute the contribution of the boundary condition */
    virtual void FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec)
    {
        // default is no specific data requirements
        int nref = datavec.size();
        if(type == 50)
        {
            for(int iref = 0; iref<nref; iref++){
                datavec[iref].fNeedsSol = true;
            }
        }
    }

    /** @brief Returns the name of the material */
    virtual std::string Name() { return "KLMaterial"; }
    
    /** @brief Returns the integrable dimension of the material */
	int Dimension() const { return fDim;}
    
    
    void SetDim(int dim)
	{
		fDim=dim;
	}
//     int Id() const 
//     { 
// 		return fId; 
// 		
// 	}
//     void SetId(int id) {
// 		
//         if(id == 0) 
// 		{
//             std::cout << "\n*** Material Id can't be ZERO! ***\n";
//             std::cout << "*** This Will Be a Disaster!!! ***\n";
//             DebugStop();
//         }
//         fId = id; 
// 		
// 	}
    
    /** @brief Returns the number of state variables associated with the material */
	int NStateVariables(){ return 1;}
    
    /** @brief Returns the number of components which form the flux function */
    virtual int NFluxes() {return 0;}
    
    /** @brief returns the number of load cases for this material object */
    int NumLoadCases()
    {
        return fNumLoadCases;
    }
    
    /** @brief returns the minimum number of load cases for this material */
    virtual int MinimumNumberofLoadCases()
    {
        return 1;
    }
    
    /** @brief changes the number of load cases for this material */
    void SetNumLoadCases(int numloadcases)
    {
        if(numloadcases <= 0)
        {
            std::cout << __PRETTY_FUNCTION__ << " numloadcases " << numloadcases << " cannot be less or equal to zero\n";
            DebugStop();
        }
        fNumLoadCases = numloadcases;
    }
    
    /** @brief indicates which variable should be post processed */
    void SetPostProcessIndex(int index)
    {
#ifdef DEBUG
        if (index < 0 || index >= fNumLoadCases)
        {
            DebugStop();
        }
#endif
        fPostProcIndex = index;
    }
	
    /** @brief Prints out the data associated with the material */
    virtual void Print(std::ostream &out)
	{
		out << "KLMaterial" << std::endl;
		out << "Lx = " << fLx << std::endl;
		out << "Ly = " << fLy <<std::endl;
		out << "Lz = " << fLz <<std::endl;
		out << "dim = " << fDim <<std::endl;
		out << "Type = " << ftype <<std::endl;
	}
    
    /** @brief Returns the variable index associated with the name */
    virtual int VariableIndex(const std::string &name);
    
    /** 
	 * @brief Returns the number of variables associated with the variable indexed by var. 
	 * @param var Index variable into the solution, is obtained by calling VariableIndex
	 */
    virtual int NSolutionVariables(int var);
    
	virtual void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout);
    
    /** @brief Computes the value of the flux function to be used by ZZ error estimator */
    virtual void Flux(TPZVec<REAL> &x, TPZVec<STATE> &Sol,
                      TPZFMatrix<STATE> &DSol, TPZFMatrix<REAL> &axes,
                      TPZVec<STATE> &flux) {}
    
    virtual void ContributeB(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek);
	
	virtual void ContributeC(TPZMaterialData &data1,TPZMaterialData &data2, REAL w1,REAL w2, TPZFMatrix<STATE> &ek);
    
	REAL AutocorrelationFunc ( TPZManVector<REAL,3>  x1, TPZManVector<REAL,3>  x2 );

	virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);


	virtual void ContributeBC(TPZMaterialData &data,REAL weight,
							  TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);

    void SetForcingFunction(TPZAutoPointer<TPZFunction<STATE> > fp)
    {
			fForcingFunction = fp;
    }
    void SetForcingFunction(void (*fp)(const TPZVec<REAL> &loc, TPZVec<STATE> &result) )
		{
				if(fp) fForcingFunction = new TPZDummyFunction<STATE>(fp);
				else fForcingFunction = NULL;
		}

	/** @brief Returns a procedure as source function for the material */
	TPZAutoPointer<TPZFunction<STATE> > &ForcingFunction() {
		return fForcingFunction;
	}
	
    /** 
	 * @brief Sets a procedure as exact solution for the problem
	 * @param fp pointer of exact solution function
	 */
	void SetForcingFunctionExact(TPZAutoPointer<TPZFunction<STATE> > fp)
	{
		fForcingFunctionExact = fp;
	}
	
    /** 
	 * @brief Sets a procedure as source function for the material.
	 * @param fp pointer of the forces function
	 * @note Parameter loc corresponds to the coordinate of the point where the source function is applied
	 * @note Parameter result contains the forces resulting
	 */
    void SetTimeDependentForcingFunction(TPZAutoPointer<TPZFunction<STATE> > fp)
    {
		fTimeDependentForcingFunction = fp;
    }
	
    /** 
	 * @brief Sets a procedure as exact solution for the problem
	 * @param fp pointer of exact solution function
	 */
	void SetTimeDependentFunctionExact(TPZAutoPointer<TPZFunction<STATE> > fp)
	{
		fTimedependentFunctionExact = fp;
	}
	
    /** 
     * @brief Sets a procedure as variable boundary condition
     * @param fp pointer of exact solution function
     */
    void SetfBCForcingFunction(TPZAutoPointer<TPZFunction<STATE> > fp)
    {
        fBCForcingFunction = fp;
    }
    
    /** 
     * @brief Sets a procedure as variable boundary condition
     * @param fp pointer of exact solution function
     */
    void SetTimedependentBCForcingFunction(TPZAutoPointer<TPZFunction<STATE> > fp)
    {
        fTimedependentBCForcingFunction = fp;
    }    
        
    
    virtual int HasForcingFunction() {return (fForcingFunction != 0);}
	virtual int HasfForcingFunctionExact() {return (fForcingFunctionExact != 0);}
    virtual int HasffBCForcingFunction() {return (fBCForcingFunction != 0);}
    virtual int HasfTimedependentBCForcingFunction() {return (fTimedependentBCForcingFunction != 0);}    
    
//     /** @brief Gets the order of the integration rule necessary to integrate an element with polinomial order p */
//     virtual int IntegrationRuleOrder(int elPMaxOrder) const
//     {
// 		return IntegrationRuleOrder( elPMaxOrder);
// 	}
// 	
// 	/** @brief Gets the order of the integration rule necessary to integrate an element multiphysic */
//     virtual int IntegrationRuleOrder(TPZVec<int> &elPMaxOrder) const
//     {
// 		IntegrationRuleOrder(elPMaxOrder);
// 	}
	
    virtual void Errors(TPZMaterialData &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors)
    {
        TPZManVector<STATE,3> flux;
        Flux(data.x, data.sol[0], data.dsol[0], data.axes, flux);
        Errors(data.x, data.sol[0], data.dsol[0], data.axes, flux, u_exact, du_exact, errors );
    }
    /**
	 * @brief Computes the error due to the difference between the interpolated flux \n
	 * and the flux computed based on the derivative of the solution
	 */
    virtual void Errors(TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol,
                        TPZFMatrix<REAL> &axes, TPZVec<STATE> &flux,
                        TPZVec<STATE> &uexact, TPZFMatrix<STATE> &duexact,
                        TPZVec<REAL> &val) {
        PZError << __PRETTY_FUNCTION__ << std::endl;
        PZError << "Method not implemented! Error comparison not available. Please, implement it." << std::endl;
    }
	virtual	void ErrorsHdiv(TPZMaterialData &data, TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values) {
		PZError << __PRETTY_FUNCTION__ << std::endl;
		PZError << "Nao sei o q fazer." << std::endl;
		
	}
    /** @brief Returns the number of norm errors. Default is 3: energy, L2 and H1. */
    virtual int NEvalErrors() {return 3;}
    
    /** @brief To create another material of the same type*/
    //virtual KLMaterial * NewMaterial();
    
    /** @brief Reads data of the material from a istream (file data)*/
    virtual void SetData(std::istream &data)
	{
		SetData(data);
	}
    
    /** @brief Creates a copy of the material object and put it in the vector which is passed on */
    virtual void Clone(std::map<int, KLMaterial * > &matvec)
	{
		Clone(matvec);
	}
    
    /** @brief To return a numerical flux type to apply over the interfaces of the elements */
    virtual int FluxType() { return 2; }
	
    virtual void ContributeErrors(TPZMaterialData &data,
                                  REAL weight,
                                  TPZVec<REAL> &nk,
                                  int &errorid){
        PZError << "Error at " << __PRETTY_FUNCTION__ << " - Method not implemented\n";
    }
    
    /**
     * @brief Computes square of residual of the differential equation at one integration point.
     * @param X is the point coordinate (x,y,z)
     * @param sol is the solution vector
     * @param dsol is the solution derivative with respect to x,y,z as computed in TPZShapeDisc::Shape2DFull
     */    
    virtual REAL ComputeSquareResidual(TPZVec<REAL>& X, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol){
        PZError << "Error at " << __PRETTY_FUNCTION__ << " - Method not implemented\n";
        return -1.;
    }
    
    /**
     * @brief Pushes a new entry in the context of materials with memory,
     * returning its index at the internal storage stack.
	 */
	/** To be implemented only in the proper materials. */
    virtual int PushMemItem(int sourceIndex)
	{
		PushMemItem(sourceIndex);
	}
    
    /** @brief Frees an entry in the material with memory internal history storage */
    virtual void FreeMemItem(int index)
	{
		FreeMemItem(index);
	}
    
    /** @brief Sets fLinearContext attribute */
    void SetLinearContext(bool IsLinear)
	{
		SetLinearContext(IsLinear);
	}
    
    /** @brief Returns fLinearContext attribute */
    bool GetLinearContext() const {
        return fLinearContext;
    }
    
    /** @{
     * @name Save and Load methods
     */
    
    /** @brief Unique identifier for serialization purposes */
    virtual int ClassId() const
    {
		return 9112;
	}
    
    /** @brief Saves the element data to a stream */
    virtual void Write(TPZStream &buf, int withclassid)
	{
		Write(buf,withclassid);
	}
    
    /** @brief Reads the element data from a stream */
    virtual void Read(TPZStream &buf, void *context)
	{
		Read(buf, context);
	}
    

	int GetExpansioOrder()
	{
		return fExpasionOrder;
	}
};

/** @brief Extern variable - Vector of force values */
extern TPZVec< void(*) (const TPZVec<REAL> &, TPZVec<STATE>& ) > GFORCINGVEC;

#endif

