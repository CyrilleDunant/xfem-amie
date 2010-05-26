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
#include "scheme_base.h"
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

void Properties::set(std::string s)
{
	set(TAG_NULL) ;
	if(s.compare("UNIVERSAL") == 0)
		set(TAG_UNIVERSAL) ; 
	if(s.compare("VOLUME") == 0)
		set(TAG_VOLUME) ;
	if(s.compare("VOLUME_FRACTION") == 0)
		set(TAG_VOLUME_FRACTION) ;
	if(s.compare("VOLUME_TOTAL") == 0)
		set(TAG_VOLUME_TOTAL) ;
	if(s.compare("MASS") == 0)
		set(TAG_MASS) ;
	if(s.compare("MASS_FRACTION") == 0)
		set(TAG_MASS_FRACTION) ;
	if(s.compare("MASS_TOTAL") == 0)
		set(TAG_MASS_TOTAL) ;
	if(s.compare("DENSITY") == 0)
		set(TAG_DENSITY) ;
	if(s.compare("YOUNG_MODULUS") == 0)
		set(TAG_YOUNG_MODULUS) ;
	if(s.compare("POISSON_RATIO") == 0)
		set(TAG_POISSON_RATIO) ;
	if(s.compare("BULK_MODULUS") == 0)
		set(TAG_BULK_MODULUS) ;	
	if(s.compare("SHEAR_MODULUS") == 0)
		set(TAG_SHEAR_MODULUS) ;
	if(s.compare("LAME_COEFFICIENT") == 0)
		set(TAG_LAME_COEFFICIENT) ;
	if(s.compare("STRAIN") == 0)
		set(TAG_STRAIN) ;
	if(s.compare("STRESS") == 0)
		set(TAG_STRESS) ;
	if(s.compare("EXPANSION_COEFFICIENT") == 0)
		set(TAG_EXPANSION_COEFFICIENT) ; 
	if(s.compare("CRACK_DENSITY") == 0)
		set(TAG_CRACK_DENSITY) ;
	if(s.compare("CRACK_LENGTH") == 0)
		set(TAG_CRACK_LENGTH) ;
	if(s.compare("ELLIPSE_A") == 0)
		set(TAG_ELLIPSE_A) ;
	if(s.compare("ELLIPSE_B") == 0)
		set(TAG_ELLIPSE_B) ;
	if(s.compare("AREA") == 0)
		set(TAG_AREA) ;
	if(s.compare("PERIMETER") == 0)
		set(TAG_PERIMETER) ;
	if(s.compare("ELLIPSE_FIRST_COMPLETE_INTEGRAL") == 0)
		set(TAG_ELLIPSE_FIRST_COMPLETE_INTEGRAL) ;
	if(s.compare("ELLIPSE_SECOND_COMPLETE_INTEGRAL") == 0)
		set(TAG_ELLIPSE_SECOND_COMPLETE_INTEGRAL) ;
	if(s.compare("CIRCLE_RADIUS") == 0)
		set(TAG_CIRCLE_RADIUS) ;
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

void Properties::print(std::string indent)
{
	std::cout << indent ;
	print() ;
}











Material::Material()
{
	name = "MAT" ;

}
Material::Material(std::string n)
{
	name = n ;
}

Material::Material(PredefinedMaterial mat)
{
	name = "MAT" ;

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
	name = "MAT" ;

	push_back(p) ;
}
Material::Material(const std::vector<Properties> & p)
{
	name = "MAT" ;

	for(size_t i = 0 ; i < p.size() ; i++)
		push_back(p[i]) ;
}
Material::Material(const Matrix & cauchy)
{
	name = "MAT" ;

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

std::vector<int> Material::getIndex(Tag t) const
{
	std::vector<int> index ;
	for(size_t i = 0 ; i < size() ; i++)
	{
		if((*this)[i].is(t))
			index.push_back((int) i) ;
	}
	return index ;
}

int Material::getIndex(Tag t, int i) const
{
	std::vector<int> index = getIndex(t) ;
	if(index.size() == 0)
		return -1 ;
	if(i+1 == 0 || i+1 > index.size())
		return index[index.size()-1] ;

	return index[i] ;
}

double Material::val(Tag t, int i) const
{
	int j = getIndex(t,i) ;
	if(j+1 == 0)
		return 0 ;
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

bool Material::set(Tag t, int i)
{
	if((*this)[i].is(t))
	{
		set(t) ;
		for(size_t j = 0 ; j < size() ; j++)
		{
			if((*this)[j].is(t) && (j != (size_t) i))
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
		int thisi = this->getIndex(compare[i],-1) ;
		int otheri = m.getIndex(compare[i],-1) ;

		if((thisi+1)*(otheri+1) == 0)
			success = false ;

		if(std::abs((*this)[thisi].val() - m[otheri].val()) > 1e-6)
			success = false ;

	}
	if(success)
	{
		int thisi = this->getIndex(combine, -1) ;
		int otheri = m.getIndex(combine, -1) ;
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
	std::cout << "---- " << name << " ----" << std::endl ;
	for(size_t i = 0 ; i < size() ; i++)
		(*this)[i].print() ;
	for(size_t i = 0 ; i < phases.size() ; i++)
		phases[i].print("|") ;
}

void Material::print(std::string indent)
{
	std::cout << indent << "---- " << name << " ----" << std::endl ;
	for(size_t i = 0 ; i < size() ; i++)
		(*this)[i].print(indent) ;
	for(size_t i = 0 ; i < phases.size() ; i++)
		phases[i].print(indent.append(indent)) ;
}

Material Material::operator*(std::string s)
{
	double f = 0. ;
	bool go_on = true ;
	int i = 0 ;
	double dec = 1. ;
	Properties prop ;
	while(go_on)
	{
		char c = s.at(i) ;
		bool valid = false ;

		if(c == '0' ||
		   c == '1' ||
		   c == '2' ||
		   c == '3' ||
		   c == '4' ||
		   c == '5' ||
		   c == '6' ||
		   c == '7' ||
		   c == '8' ||
		   c == '9')
		{
			if(dec < 1.)
			{
				f += dec * atoi(&c) ;
				dec *= 0.1 ;
			} else {
				f = f*10 + atoi(&c) ;
			}
			valid = true ;
		}

		if(c == '.')
		{
			dec = 0.1 ;
			valid = true ;
		}

		if(c == '_')
		{
			std::string tag = s.substr(i+1) ;
			prop.set(tag) ;
			prop.set(f) ;
			valid = false ;
		}

		go_on = valid ;
		i++ ;
		if(i+1 > (int) s.length())
			go_on = false ;

	}

	replace(prop) ;

	return (*this) ;

}

Material Material::operator+(Material m)
{
	phases.push_back(m) ;
	return (*this) ;
}

bool Material::build(Scheme * s)
{
	merge(s->homogenize(phases)) ;
	s->print() ;
	return s->isOK() ;
}
