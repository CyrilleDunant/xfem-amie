//
// C++ Implementation: mechanical analytic homogenization
//
// Description:
//
//
// Author:  Alain Giorla, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "properties_base.h"
#include "converter.h"
#include "elastic_homogenization.h"

using namespace Mu ;

Properties::Properties()
{
	ptag = TAG_NULL ;
	p = 0. ;
}
Properties::Properties(double v)
{
	ptag = TAG_UNIVERSAL ;
	p = v ;
}
Properties::Properties(Tag t, double v)
{
	ptag = t ;
	p = v ;
}
Properties::Properties(const Properties & prop)
{
	ptag = prop.tag() ;
	p = prop.val() ;
}

void Properties::print()
{
	switch(ptag)
	{
	case TAG_NULL:
		std::cout << "NULL PROPERTIES" << std::endl ;
		break;
	case TAG_UNIVERSAL:
		std::cout << "UNDEFINED PROPERTIES = " ;
		break;
	case TAG_VOLUME:
		std::cout << "VOLUME = " ;
		break;			
	case TAG_VOLUME_FRACTION: 
		std::cout << "VOLUME FRACTION = " ; 
		break;			
	case TAG_VOLUME_TOTAL:
		std::cout << "VOLUME TOTAL = " ; 
		break;			
	case TAG_MASS:
		std::cout << "MASS = " ;
		break;			
	case TAG_MASS_FRACTION:
		std::cout << "MASS FRACTION = " ; 
		break;			
	case TAG_MASS_TOTAL:
		std::cout << "MASS TOTAL = " ; 
		break;			
	case TAG_DENSITY:
		std::cout << "DENSITY = " ; 
		break;			
	case TAG_YOUNG_MODULUS:
		std::cout << "YOUNG MODULUS = " ; 
		break;			
	case TAG_POISSON_RATIO:
		std::cout << "POISSON RATIO = " ; 
		break;			
	case TAG_BULK_MODULUS:
		std::cout << "BULK MODULUS = " ; 
		break;			
	case TAG_SHEAR_MODULUS:
		std::cout << "SHEAR MODULUS = " ; 
		break;			
	case TAG_EXPANSION_COEFFICIENT: 
		std::cout << "EXPANSION COEFFICIENT = " ; 
		break;			
	case TAG_CRACK_DENSITY: 
		std::cout << "CRACK DENSITY = " ;
		break;			
	case TAG_ELLIPSE_A:
		std::cout << "ELLIPSE MAJOR RADIUS = " ; 
		break;			
	case TAG_ELLIPSE_B:
		std::cout << "ELLIPSE MINOR RADIUS = " ;
		break;			
	case TAG_AREA:
		std::cout << "AREA = " ;
		break;			
	case TAG_PERIMETER:
		std::cout << "PERIMETER = " ;
		break;			
	case TAG_ELLIPSE_FIRST_COMPLETE_INTEGRAL:
		std::cout << "ELLIPTIC FIRST COMPLETE INTEGRAL = " ;
		break;			
	case TAG_ELLIPSE_SECOND_COMPLETE_INTEGRAL:
		std::cout << "ELLIPTIC SECOND COMPLETE INTEGRAL = " ; 
		break;			
	case TAG_CIRCLE_RADIUS:
		std::cout << "CIRCLE RADIUS = " ;
		break;			

	}
	if(!isNull())
		std::cout << p << std::endl ;
}












Material::Material()
{

}
Material::Material(PredefinedMaterial mat)
{
	switch(mat)
	{
	case MAT_DUMMY:
		push_back(Properties(TAG_YOUNG_MODULUS,1.)) ;
		push_back(Properties(TAG_POISSON_RATIO,0.2)) ;
		push_back(Properties(TAG_DENSITY,1.)) ;
		break;
	case MAT_AGGREGATE:
		push_back(Properties(TAG_YOUNG_MODULUS,60.)) ;
		push_back(Properties(TAG_POISSON_RATIO,0.2)) ;
		push_back(Properties(TAG_DENSITY,2.7)) ;
		break;
	case MAT_CEMENT:
		push_back(Properties(TAG_YOUNG_MODULUS,14.)) ;
		push_back(Properties(TAG_POISSON_RATIO,0.2)) ;
		push_back(Properties(TAG_DENSITY,1.5)) ;
		break;
	}
//	this->print() ;
}
Material::Material(const Properties & p)
{
	push_back(p) ;
}
Material::Material(const std::vector<Properties> & p)
{
	for(size_t i = 0 ; i < p.size() ; i++)
		push_back(p[i]) ;
}
Material::Material(const Matrix & cauchy)
{
	double E = 0. ;
	double nu = 0. ;
	double k = 0. ;
	double mu = 0. ;

	double A = cauchy[0][0] ;
	double B = cauchy[0][1] ;

	if(cauchy.size() > 9)
	{
		nu = B / (A+B) ;
		E = (1.+nu)*(1.-2.*nu)*B/nu ;
	} else {
		nu = B / A ;
		E = A * (1.-nu*nu) ;
	}

	GeneralConverter conv(TAG_UNIVERSAL) ;

	k = conv.getBulkModulus(E,nu) ;
	mu = conv.getShearModulus(E,nu) ;

	if(conv.isOK() || conv.status() == STATUS_RESET)
	{
		push_back(Properties(TAG_YOUNG_MODULUS,E)) ;
		push_back(Properties(TAG_POISSON_RATIO,nu)) ;
		push_back(Properties(TAG_BULK_MODULUS,k)) ;
		push_back(Properties(TAG_SHEAR_MODULUS,mu)) ;
	}

}

std::vector<size_t> Material::getIndex(Tag t) const
{
	std::vector<size_t> index ;
	for(size_t i = 0 ; i < size() ; i++)
	{
		if((*this)[i].is(t))
			index.push_back(i) ;
	}
	return index ;
}

size_t Material::getIndex(Tag t, size_t i) const
{
	std::vector<size_t> index = getIndex(t) ;
	if(index.size() == 0)
		return -1 ;
	if(i+1 == 0 || i+1 > index.size())
		return index[index.size()-1] ;

	return index[i] ;
}

double Material::val(Tag t, size_t i) const
{
	size_t j = getIndex(t,i) ;
	if(j+1 == 0)
		return -1 ;
	return (*this)[j].val() ;
}

bool Material::replace(Properties p)
{
	if(isSet(p.tag()))
		return false ;

	for(size_t i = 0 ; i < size() ; i++)
	{
		if((*this)[i].is(p.tag()))
			(*this)[i].kill() ;
	}
	push_back(p) ;
	return true ;
}

bool Material::replaceForce(Properties p)
{
	if(replace(p))
		return true ;

	for(size_t i = 0 ; i < size() ; i++)
	{
		if((*this)[i].is(p.tag()))
			(*this)[i].kill() ;
	}
	push_back(p) ;
	return false ;
	
}

bool Material::isSet(Tag t) const
{
	for(size_t i = 0 ; i < tagset.size() ; i++)
	{
		if(tagset[i] == t)
			return true ;
	}
	return false ;
}

bool Material::set(Tag t)
{
	if(isSet(t))
		return true ;

	if(getIndex(t).size() == 0)
		return false ;

	tagset.push_back(t) ;
	return true ;
}

bool Material::set(Tag t, size_t i)
{
	if((*this)[i].is(t))
	{
		set(t) ;
		for(int j = 0 ; j < size() ; j++)
		{
			if((*this)[j].is(t) && (j != i))
				(*this)[j].kill() ; 
		}
		return true ;
	}
	return false ;
}

bool Material::combine(Material m, std::vector<Tag> compare, Tag combine)
{
	bool success = true ;
	for(size_t i = 0 ; i < compare.size() && success ; i++)
	{
		size_t thisi = this->getIndex(compare[i],-1) ;
		size_t otheri = m.getIndex(compare[i],-1) ;

		if((thisi+1)*(otheri+1) == 0)
			success = false ;

		if(std::abs((*this)[thisi].val() - m[otheri].val()) > 1e-6)
			success = false ;

	}
	if(success)
	{
		size_t thisi = this->getIndex(combine, -1) ;
		size_t otheri = m.getIndex(combine, -1) ;
		if((thisi+1)*(otheri+1) == 0)
			success = false ;
		else
			(*this)[thisi].set((*this)[thisi].val() + m[otheri].val()) ;
	}

	return success ;
}

bool Material::merge(Material m)
{
	bool success = true ;
	for(size_t i = 0 ; i < m.size() ; i++)
		success &= replace(m[i]) ;
	return success ;
}

bool Material::findMissing(std::vector<Tag> t)
{
	bool found = true ;
	for(size_t i = 0 ; i < t.size() ; i++)
		found = found && findMissing(t[i]) ;
	return found ;
}

bool Material::findMissing(Tag t)
{
//	Properties dummy(t,0.) ;
//	dummy.print() ;

	if(getIndex(t).size() > 0)
		return true ;

	GeneralConverter conv(t) ;
	std::vector<Properties> prop = conv.homogenize(*this) ;
//	prop[0].print() ;
//	conv.print() ;
	if(conv.isOK())
	{
//		prop[0].print() ;
		this->push_back(prop[0]) ;
	}
	return conv.isOK() ;
}

void Material::print()
{
	for(size_t i = 0 ; i < size() ; i++)
		(*this)[i].print() ;
}


