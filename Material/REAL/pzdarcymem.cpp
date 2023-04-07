//$Id: pzelastoplasticmem.cpp,v 1.6 2023-04-06 00:55:14 diogo Exp $

#include "pzdarcymem.h"


TPZDarcyMem::TPZDarcyMem():  fpermeability(0){ }
	
TPZDarcyMem::TPZDarcyMem(const TPZDarcyMem & source):
                          fpermeability(source.fpermeability) { }


TPZDarcyMem::~TPZDarcyMem(){ }


const TPZDarcyMem & TPZDarcyMem::operator=(const TPZDarcyMem & source)
{
	fpermeability = source.fpermeability;
	return *this;
}

const std::string TPZDarcyMem::Name()const
{
		return "TPZDarcyMem";
}

const int TPZDarcyMem::ClassId()const
{
		return 763;
}

void TPZDarcyMem::Write(TPZStream &buf, int withclassid)
{
	buf.Write(&fpermeability[0],3);


}

void TPZDarcyMem::Read(TPZStream &buf, void *context)
{
	buf.Read(&fpermeability[0],3);
}

void TPZDarcyMem::Print(std::ostream &out)const
{
	out << "\TPZDarcyMem";
	out << "\fpermeability = " << fpermeability;

}
