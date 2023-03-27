/**
 * @file
 * @brief Contains the TPZExtendGridDimension class which generates a three dimensional mesh as an extension of a two dimensional mesh.
 */

/**
 * This class extends a dimension mesh n = 1,2 
 * for a dimension mesh n+1 = 2,3 respectively, 
 * the entrance data are the fine mesh and the thickness of the new mesh
 * restrictions:
 * - for now the mesh it must be plain
 * - the nodes of the mesh must be given with three co-ordinated and counter-clock sense
 */

#ifndef TPZEXTENDGRID_HPP
#define TPZEXTENDGRID_HPP

class TPZGeoMesh;
#include "pzstack.h"
#include <fstream>
#include <iostream>
#include "pzmatrix.h"

#include "tpzautopointer.h"


/**
 * @ingroup pre
 * @brief Generates a three dimensional mesh as an extension of a two dimensional mesh. \ref pre "Getting Data"
 */
class TPZExtendGridDimension {
	
	/**
	 * @brief Thickness of the mesh (+  or -)
	 */
	REAL fThickness;
	
	/**
	 * @brief Name of the fine mesh to be extended
	 */
	std::ifstream fFineFileMesh;
	
	/**
	 * @brief Fine geometric mesh generated by the NeoPZ
	 */
	TPZAutoPointer<TPZGeoMesh> fFineGeoMesh;
    
	/**
	 * @brief Vector of n surfaces to be connected the first one correpsonds to the base
	 */
    TPZVec<TPZAutoPointer<TPZGeoMesh> > fSurfaces;
	
public:
	/** @brief Constructor using filename with gmesh data and thickness */
	TPZExtendGridDimension(char *geofile,REAL thickness);
	/** @brief Constructor using geometric mesh one- or two- dimensional and thickness */
	TPZExtendGridDimension(TPZAutoPointer<TPZGeoMesh> &finegeomesh,REAL thickness);
	TPZExtendGridDimension(TPZGeoMesh* finegeomesh,REAL thickness);

	/** @brief Destructor default */
	~TPZExtendGridDimension(){};
    
	/**
	 * @brief It reads the mesh since the archive of entrance finemesh, or since the fFineGeoMesh
	 * passed in the constructor, and returns extended mesh.
	 */
    TPZGeoMesh* ExtendedMesh();
    
    /**
     * @brief Apply transformation to a given geomesh
     */
    static  void DeformMesh(TPZFMatrix<STATE> &Tr, TPZGeoMesh * GeoSurface);
    
	/**
	 * @brief It reads the mesh since the archive of entrance finemesh, or since the fFineGeoMesh
	 * passed in the constructor, and returns extended mesh. \n
	 * The extension is from n (=1,2) dimensional mesh to (n+1) (=2,3) dimensional.
	 * @param naumentedlayers Numbers of layers to be incremented.
	 * @param matidbottom Material id to bottom boundary surface after to extrude process.
	 * @param matidtop Material id to top boundary surface after to extrude process.
	 */
	TPZGeoMesh* ExtendedMesh(int naumentedlayers,int matidbottom=0,int matidtop=0);

	/**
	 * @brief Prints the generated mesh
	 */
	void PrintGeneratedMesh(std::ostream &out = std::cout);

};

#endif
