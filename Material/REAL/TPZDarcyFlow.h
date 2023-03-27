//
// Created by Gustavo Batistela on 5/13/21.
//

#ifndef TPZDARCYFLOW_H
#define TPZDARCYFLOW_H

#include "pzmaterialdata.h"
#include "pzmaterial.h"
#include "pzfmatrix.h"
#include "pzbndcond.h"
#include "pzlog.h"
#include "TPZIsotropicPermeability.h"

/**
 * @ingroup material
 * @brief This class implements an H1-conforming approximation for the Darcy flow equation for isotropic materials.
 *
 * The Darcy flow equation is given by: \f[- \nabla \cdot K \nabla u = f,\f] where \f$u\f$ is the pressure field
 * to be solved, \f$K\f$ is the permeability tensor and \f$f\f$ is the source term.
 *
 * @see TPZMixedDarcyFlow For an approximation using the mixed method.
 * @see TPZIsotropicPermeability For setting the permeability field.
 */

class TPZDarcyFlow : public  TPZMaterial {

    // type alias to improve constructor readability

public:
    /**
     * @brief Default constructor
     */
    TPZDarcyFlow();

    /**
	 * @brief Class constructor
	 * @param [in] id material id
	 * @param [in] dim problem dimension
	 */
    TPZDarcyFlow(int id, int dim);

    /**
	 * @brief Returns a 'std::string' with the name of the material
	 */
    virtual std::string Name()   const{ return "TPZDarcyFlow"; }

    /**
	 * @brief Returns the problem dimension
	 */
    virtual int Dimension() const override { return this->fDim; }

    /**
	 * @brief Returns the number of state variables
	 */
    virtual int NStateVariables(){ return 1;}

    /**
	 * @brief Returns the number of errors to be evaluated
     *
     * Returns the number of errors to be evaluated, that is, the number of error norms associated
     * with the problem.
     */
    int NEvalErrors() const  { return 6; }

    /**
     * @brief Sets problem dimension
     */
    virtual void SetDimension(int dim);

    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point
     * @param[in] data stores all input data
     * @param[in] weight is the weight of the integration rule
     * @param[out] ek is the element matrix
     * @param[out] ef is the rhs vector
     */
     void Contribute( TPZMaterialData &data, STATE weight, TPZFMatrix<STATE> &ek,
                    TPZFMatrix<STATE> &ef) ;

    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point
     * @param[in] data stores all input data
     * @param[in] weight is the weight of the integration rule
     * @param[out] ek is the element matrix
     * @param[out] ef is the rhs vector
     * @param[in] bc is the boundary condition material
     */
    void ContributeBC( TPZMaterialData &data, STATE weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef,
                      TPZBndCond &bc) ;

    /**
     * @brief Returns an integer associated with a post-processing variable name
     * @param [in] name string containing the name of the post-processing variable. Ex: "Pressure".
     */
     int VariableIndex(const std::string &name)  ;

    /**
     * @brief Returns an integer with the dimension of a post-processing variable
     * @param [in] var index of the post-processing variable, according to TPZDarcyFlow::VariableIndex method.
     */
     int NSolutionVariables(int var)  ;

    /**
     * @brief Returns the solution associated with the var index based on the
     * finite element approximation at a point
     * @param [in] data material data associated with a given integration point
     * @param [in] var index of the variable to be calculated
     * @param [out] solOut vector to store the solution
     */
    void Solution( TPZMaterialData &data, int var, TPZVec<STATE> &solOut) ;

    /**
     * @brief Get the dimensions of the solution for each state variable
     *
     * This will be used for initializing the corresponding TPZMaterialData
     * @param [out] u_len solution dimension
     * @param [out] du_row number of rows of the derivative
     * @param [out] du_col number of columns of the derivative
     */
    void GetSolDimensions(uint64_t &u_len, uint64_t &du_row, uint64_t &du_col) const ;

    /**
     * @brief Calculates the approximation error at a point
     * @param [in] data material data of the integration point
     * @param [out] errors calculated errors
     */
    void Errors( TPZMaterialData &data, TPZVec<REAL> &errors) ;

    /*
     * @brief fill requirements for volumetric contribute
     */
    void FillDataRequirements(TPZMaterialData &data) const ;

    /*
     * @brief fill requirements for boundary contribute
     */
    void FillBoundaryConditionDataRequirements(int type, TPZMaterialData &data) const ;

    /**
     * @brief Returns an unique class identifier
     */
    [[nodiscard]] int ClassId() const ;

    /**
     * @brief Creates another material of the same type
     */
   // [[nodiscard]] TPZMaterial *NewMaterial() const ;

    /**
     * @brief Prints data associated with the material.
     */
    void Print(std::ostream & out) const ;
	
	void SetConstantPermeability(REAL perm)
	{
		fConstPermeability =  perm;
	}

protected:
    /**
     * @brief Problem dimension
     */
    int fDim;
	
	REAL fConstPermeability;
};

#endif //TPZDARCYFLOW_H
