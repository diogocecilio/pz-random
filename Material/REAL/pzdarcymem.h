//$Id: pzelastoplasticmem.cpp,v 1.6 2023-04-06 00:55:14 diogo  $

#ifndef TPZDarcyMem_H
#define TPZDarcyMem_H

#include "pzfmatrix.h"
  /**
   * This class defines the material memory that should be stored at each integration point
   * for the purposes of an elastoplastic material.
   */

class TPZDarcyMem
{
public:
	TPZDarcyMem();
	
	TPZDarcyMem(const TPZDarcyMem & source);
	
	const TPZDarcyMem & operator=(const TPZDarcyMem & source);
	
	virtual ~TPZDarcyMem();


	const std::string Name()const;

	const int ClassId()const;

    void Write(TPZStream &buf, int withclassid);

    void Read(TPZStream &buf, void *context);

	virtual void Print(std::ostream &out = std::cout)const;


    friend std::ostream& operator<<( std::ostream& Out, const TPZDarcyMem & s )
	{
		s.Print(Out);
		return Out;
	}

	TPZManVector<REAL,3> fpermeability;
	
};



#endif
