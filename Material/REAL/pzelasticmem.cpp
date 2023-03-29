//$Id: pzelastoplasticmem.cpp,v 1.6 2009-06-22 00:55:14 erick Exp $

#include "pzelasticmem.h"


TPZElasticMem::TPZElasticMem():  fE(0),fnu(0.){ }
	
TPZElasticMem::TPZElasticMem(const TPZElasticMem & source):
                          fE(source.fE), fnu(source.fnu) { }


TPZElasticMem::~TPZElasticMem(){ }


const TPZElasticMem & TPZElasticMem::operator=(const TPZElasticMem & source)
{
	fE = source.fE;
	fnu = source.fnu;
	return *this;
}

const std::string TPZElasticMem::Name()const
{
		return "TPZElasticMem";
}

const int TPZElasticMem::ClassId()const
{
		return 762;
}

void TPZElasticMem::Write(TPZStream &buf, int withclassid)
{
	buf.Write(&fE,1);
	buf.Write(&fnu,1);
}

void TPZElasticMem::Read(TPZStream &buf, void *context)
{
	buf.Read(&fE,1);
	buf.Read(&fnu,1);
}

void TPZElasticMem::Print(std::ostream &out)const
{

	out << "\ElasticMem";
	out << "\fE = " << fE;
	out << "\fnu = " << fnu;

}
