/**
 * @file
 */

#ifndef PZELASTOPLASTIC_H
#define PZELASTOPLASTIC_H


#include "pzmaterial.h"
#include "pzmatwithmem.h"
#include "pzelastoplasticmem.h"
#include "pzporoelastoplasticmem.h"

  /**
   * Implements an elastoplastic material and uses the memory feature to store the damage variables
   * This material works only together with the Plasticity Library.
   */

//typedef TPZMatWithMem<TPZElastoPlasticMem> BASE_MATWITHMEM;

template <class T, class TMEM = TPZElastoPlasticMem>
class  TPZMatElastoPlastic : public TPZMatWithMem<TMEM>
{
   public:

      enum SOLUTIONVARS{ENone = -1,
	         
		  
	  EStrainVol = 0,
	  EStrainXX = 1,
	  EStrainYY = 2,
	  EStrainZZ = 3,
	  EStrainXY = 4,
	  EStrainXZ = 5,
	  EStrainYZ = 6,
	  // Elastic Strain
	  EElStrainVol = 7,
	  EElStrainXX = 8,
	  EElStrainYY = 9,
	  EElStrainZZ = 10,
	  EElStrainXY = 11,
	  EElStrainXZ = 12,
	  EElStrainYZ = 13,
	  // Plastic Strain
	  EPlStrainVol = 14,
	  EPlStrainXX = 15,
	  EPlStrainYY = 16,
	  EPlStrainZZ = 17,
	  EPlStrainXY = 18,
	  EPlStrainXZ = 19,
	  EPlStrainYZ = 20,
	  EPlStrainSqJ2 = 21,
	  EPlStrainSqJ2El = 22,
	  EPlAlpha = 23,
	  // Displacement
	  EDisplacementX = 24,
	  EDisplacementY = 25,
	  EDisplacementZ = 26,
	  EDisplacementTotal = 27,
	  // Total Stress
	  ETotStressI1 = 28,
	  ETotStressJ2 = 29,
	  ETotStressXX = 30,
	  ETotStressYY = 31,
	  ETotStressZZ = 32,
	  ETotStressXY = 33,
	  ETotStressXZ = 34,
	  ETotStressYZ = 35,
	  ETotStress1 = 36,
	  ETotStress2 = 37,
	  ETotStress3 = 38,
	  // Effective stress
	  EEffStressI1 = 39,
	  EEffStressJ2 = 40,
	  EEffStressXX = 41,
	  EEffStressYY = 42,
	  EEffStressZZ = 43,
	  EEffStressXY = 44,
	  EEffStressXZ = 45,
	  EEffStressYZ = 46,
	  EEffStress1 = 47,
	  EEffStress2 = 48,
	  EEffStress3 = 49,
	  // Yield Surface
	  EYieldSurface1 = 50,
	  EYieldSurface2 = 51,
	  EYieldSurface3 = 52,
	  // Simulation
	  EPOrder = 53,
	  ENSteps = 54,
	  // Pore pressure
	  EPorePressure = 55,
	  // Material
	  EMatPorosity = 56,
	  EMatE = 57,
	  EMatPoisson = 58,
	  ECohes =59,
	  EFric = 60,
      EFluxX=61,
      EFluxY=62,
      EFlux=63,
      EPressure=64,
      EPrincipalStress=65,
      EEnergy=66
};
		
		
      /**
       * Default constructor
       */
      TPZMatElastoPlastic();		
		
      /** Creates a material object and inserts it in the vector of
       *  material pointers of the mesh. Upon return vectorindex
       *  contains the index of the material object within the
       *  vector
       */
      TPZMatElastoPlastic(int id);

      /** Creates a material object based on the referred object and
       *  inserts it in the vector of material pointers of the mesh.
       *  Upon return vectorindex contains the index of the material
       *  object within the vector
       */
      TPZMatElastoPlastic(const TPZMatElastoPlastic &mat);

      virtual ~TPZMatElastoPlastic();

	  /** Sets the plasticity model already with proper parameters */
      void SetPlasticity(T & plasticity);
	  
	  T GetPlasticity()
	  {
		return fPlasticity;  
	  }
	
	  /** Sets the material bulk density */
      void SetBulkDensity(REAL & RhoB);
	
      /** returns the name of the material*/
      virtual std::string Name();

      /**returns the integrable dimension of the material*/
      virtual int Dimension() const { return 3; }

      /** returns the number of state variables associated with the material*/
      virtual int NStateVariables() { return 3; }

      /** print out the data associated with the material*/
      virtual void Print(std::ostream &out, const int memory);

    /** print out the data associated with the material*/
    virtual void Print(std::ostream &out);
    
      /**returns the variable index associated with the name*/
      virtual int VariableIndex(const std::string &name);

      /** returns the number of variables associated with the variable
	  indexed by var.  var is obtained by calling VariableIndex*/
      virtual int NSolutionVariables(int var);

      /**returns the solution associated with the var index based on
       * the finite element approximation*/
      virtual void Solution(TPZMaterialData &data, int var, TPZVec<REAL> &Solout);
	
      /** Return the number of components which form the flux function
       * Method not implemented.
       */
      virtual int NFluxes() 
	  {
         PZError << "TPZMatElastoPlastic::NFluxes() - Method not implemented\n";
         return 0;
      }

      /** Compute the value of the flux function to be used by ZZ error estimator.
       * Method not implemented.
       */
      virtual void Flux(TPZVec<REAL> &x, TPZVec<REAL> &Sol, TPZFMatrix<REAL> &DSol, TPZFMatrix<REAL> &axes, TPZVec<REAL> &flux)
	  {
         PZError << "TPZMatElastoPlastic::Flux - Method not implemented\n";
      }

      /** Evaluate error between approximate (FEM) and exact solutions.
	   *  Method not implemented
       */
      virtual void Errors(TPZVec<REAL> &x,TPZVec<REAL> &u, TPZFMatrix<REAL> &dudx,
                          TPZFMatrix<REAL> &axes, TPZVec<REAL> &flux,
                          TPZVec<REAL> &u_exact,TPZFMatrix<REAL> &du_exact,TPZVec<REAL> &values);
      /**
       * Returns the number of norm errors: 3 (Semi H1, L2 and H1)
	   * Method not implemented
       */
      virtual int NEvalErrors() {return 3;}

      /**
       * It computes a contribution to the stiffness matrix and load vector at one integration point.
       * @param data [in] stores all input data
       * @param weight [in] is the weight of the integration rule
       * @param ek [out] is the stiffness matrix
       * @param ef [out] is the load vector
       */
      virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef);

      /**
       * It computes a contribution to the stiffness matrix and load vector at one BC integration point.
       * @param data [in] stores all input data
       * @param weight [in] is the weight of the integration rule
       * @param ek [out] is the stiffness matrix
       * @param ef [out] is the load vector
       * @param bc [in] is the boundary condition material
       */
      virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef, TPZBndCond &bc);

      /**
       * It computes a contribution to the residual vector at one integration point.
       * @param data [in] stores all input data
       * @param weight [in] is the weight of the integration rule
       * @param ef [out] is the residual vector
       */
      virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ef);

      /**
       * It computes a contribution to the stiffness matrix and load vector at one BC integration point.
       * @param data [in] stores all input data
       * @param weight [in] is the weight of the integration rule
       * @param ef [out] is the load vector
       * @param bc [in] is the boundary condition material
       */
      virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ef, TPZBndCond &bc);

	  /** Evaluates the Strain vector based on an available DSol (solution derivatives set) vector.
	   * @param data [in] 
	   * @param Strain [out]
	   */
	  void ComputeStrainVector(TPZMaterialData & data, TPZFMatrix<REAL> &Strain);
	
	  /** Evaluates the Strain vector based on an available DSol (solution derivatives set) vector.
	   * @param DeltaStrain [out] 
	   * @param data [in]
	   */
	  void ComputeDeltaStrainVector(TPZMaterialData & data, TPZFMatrix<REAL> &DeltaStrain);
	
	  /** Evaluates the Stress vector based on an available DSol (solution derivatives set) vector.
	   * @param data [in]
	   * @param Stress [out] 
	   */
	  void ComputeStressVector(TPZMaterialData & data, TPZFMatrix<REAL> &Stress);
	
    /** Method that checks DEP consistency
     *  @param data [in]
     *  @param DeltaStrain [in]
     */
    void CheckConvergence(TPZMaterialData & data,TPZFMatrix<REAL> & DeltaStrain);
    
	  /** Calls the plasticity template aggregate applyStrainComputeDep method
	   *  @param data [in]
	   *  @param DeltaStrain [in]
	   *  @param Stress [out]
	   *  @param Dep [out]
	   */
	 void ApplyDeltaStrainComputeDep(TPZMaterialData & data, TPZFMatrix<REAL> & DeltaStrain, 
												TPZFMatrix<REAL> & Stress, TPZFMatrix<REAL> & Dep);
	
	  /** Calls the plasticity template aggregate applyStrain method
	   *  @param data [in]
	   *  @param DeltaStrain [in]
	   *  @param Stress [out]
	   */
	  void ApplyDeltaStrain(TPZMaterialData & data, TPZFMatrix<REAL> & DeltaStrain, 
												TPZFMatrix<REAL> & Stress);

	  /** Applies the tensor in vectorial form to an internally stored direction
	   * @param vectorTensor [in]
	   * @param Out [out]
	   */
      void ApplyDirection(TPZFMatrix<REAL> &vectorTensor, TPZVec<REAL> &Out);
	
	  /** Converts the stress vector onto a symmetric stress tensor
	   * @param vectorTensor [in]
	   * @param Tensor [out]
	   */
	  void vectorToTensor(const TPZFMatrix<REAL> & vectorTensor, TPZFMatrix<REAL> & Tensor);
	
	  /** Evaluates the eigenvalues of the tensor vectorTensor (in compact vectorial form)
	   * @param vectorTensor [in] compact vectorial form of a symmetrical tensor
	   * @param ev [out] evaluated eigenvalues
	   */
	  void EigenValues(TPZFMatrix<REAL> & vectorTensor, TPZVec<REAL> & ev);

	  /** Evaluates the eigenvectors of the tensor vectorTensor (in compact vectorial form)
	   * @param vectorTensor [in] compact vectorial form of a symmetrical tensor
	   * @param Solout [out] evaluated eigenvector
	   * @param direction [in] selected direction (0 to 2)
	   */
      void EigenVectors(TPZFMatrix<REAL> &vectorTensor, TPZVec< REAL > &Solout, int direction);
	
      /**To create another material of the same type*/
      virtual TPZMaterial * NewMaterial();

      /**
       * Unique identifier for serialization purposes
       */
      virtual int ClassId() const;

      /**
       * Save the element data to a stream
       */
      virtual void Write(TPZStream &buf, int withclassid);

      /**
       * Read the element data from a stream
       */
      virtual void Read(TPZStream &buf, void *context);
	
	  /**
	   * Sets the tolerance value for post-processing purposes
	   */
	  void SetTol(const REAL & tol);

	/**
	 * Sets the SetBulkDensity of the material
	 */
	void SetBulkDensity(const REAL & bulk);
	
	

	  /**
	   * Defining what parameters the material needs. In particular this material needs the
	   * evaluation of normal vector for the sake of boundary conditions
	   */
		virtual void FillDataRequirements(TPZMaterialData &data);

        /**
         * This method defines which parameters need to be initialized in order to compute the contribution of the boundary condition
         */
		virtual void FillBoundaryConditionDataRequirement(int type,TPZMaterialData &data);


		void SetBodyForce(TPZManVector<REAL, 3> Force)
		{
			fForce = Force;
		}
		void SetLoadFactor(REAL factor)
		{
			ffactor = factor;
		}
		void SetWhichLoadVector(int option)
		{
			fwhichinternalforce = option;
		}
		
protected:

	/**
       * gravity acceleration
       */
      TPZManVector<REAL, 3> fForce;

	  REAL ffactor;


	  //0 ef = (Bt sigma + (b + gradu))
	  //1 ef = Bt sigma
	  //2 ef = b + gradu
	  int fwhichinternalforce;
	
	  /**
	   * bulk density of rock
	   */
	//  REAL fRhoB; 	
 
	  /**
	   * Post Processing direction
	   */
      TPZManVector<REAL,3> fPostProcessDirection;
	
	  /**
	   * Elastoplastic material object instantiation
	   * this instantiation avoids several instantiations
	   * for each use of this object
	   */
	
	  T fPlasticity;
	
	  /**
	   * Tolerance for posto-processing purposes
	   */
	  REAL fTol;
	
};

#endif
