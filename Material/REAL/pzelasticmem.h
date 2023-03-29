//$Id: pzelastoplasticmem.h,v 1.7 2009-10-04 05:44:22 erick Exp $

#ifndef PZELASTICMEM_H
#define PZELASTICMEM_H

#include "pzfmatrix.h"
  /**
   * This class defines the material memory that should be stored at each integration point
   * for the purposes of an elastoplastic material.
   */

class TPZElasticMem
{
public:
	TPZElasticMem();
	
	TPZElasticMem(const TPZElasticMem & source);
	
	const TPZElasticMem & operator=(const TPZElasticMem & source);
	
	virtual ~TPZElasticMem();


	const std::string Name()const;

	const int ClassId()const;

    void Write(TPZStream &buf, int withclassid);

    void Read(TPZStream &buf, void *context);

	virtual void Print(std::ostream &out = std::cout)const;


    friend std::ostream& operator<<( std::ostream& Out, const TPZElasticMem & s )
	{
		s.Print(Out);
		return Out;
	}


	REAL fE;

	REAL fnu;

	
};



#endif
