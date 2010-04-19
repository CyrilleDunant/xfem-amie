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
#include "elastic_homogenization.h"

using namespace Mu ;


size_t Mu::standardNVal(Mu::PropertiesType p)
{
	switch(p)
	{
		case(VOID_PROP):
			return 0 ;
		case(ABSTRACT):
			return -1 ;
		case(FRACTION):
			return 1 ;
		case(HOOKE):
			return 2 ;
		case(BULK_SHEAR):
			return 2 ;
		case(EXPANSION):
			return 1 ;
	}
	return 0 ;
}

bool Mu::conversionPossible(Mu::PropertiesType p1, Mu::PropertiesType p2)
{
	switch(p2)
	{
		case VOID_PROP:
			return true ;
		case ABSTRACT:
			return true ;
		case FRACTION:
			return (p1 == ABSTRACT || p1 == FRACTION) ;
		case HOOKE:
			return (p1 == ABSTRACT || p1 == HOOKE || p1 == BULK_SHEAR) ;
		case BULK_SHEAR:
			return (p1 == ABSTRACT || p1 == HOOKE || p1 == BULK_SHEAR) ;
		case EXPANSION:
			return (p1 == ABSTRACT || p1 == EXPANSION) ;
	}
	return false ;
}

Properties::Properties()
{
	pType = VOID_PROP ;
	nVal = 0 ;
}

Properties::Properties(double v)
{
	pType = ABSTRACT ;
	values.resize(1) ;
	values[0] = v ;
	nVal = 1 ;
}

Properties::Properties(const std::pair<double, double> & v)
{
	pType = ABSTRACT ;
	values.resize(2) ;
	values[0] = (v.first) ;
	values[1] = (v.second) ;
	nVal = 2 ;
}

Properties::Properties(const std::vector<double> & v)
{
	pType = ABSTRACT ;
	values.resize(v.size()) ;
	for(size_t i = 0 ; i < v.size() ; i++)
		values[i]=(v[i]) ;
	nVal = values.size() ;
}

Properties::Properties(const Vector & v)
{
	pType = ABSTRACT ;
	values.resize(v.size()) ;
	for(size_t i = 0 ; i < v.size() ; i++)
		values[i]=(v[i]) ;
	nVal = values.size() ;
}

Properties::Properties(PropertiesType p, double v)
{
	pType = p ;
	nVal = standardNVal(p) ;
	if(nVal != 0)
	{
		if(nVal+1 > 0)
		{
			values.resize(nVal) ;
		} else {
			values.resize(1) ;
		}
		values[0] = v ;
	} else {
		values.resize(0) ;
	}
}

Properties::Properties(PropertiesType p, const std::pair<double, double> & v)
{
	pType = p ;
	nVal = standardNVal(p) ;
	if(nVal != 0)
	{
		if(nVal+1 > 0)
		{
			values.resize(nVal) ;
		} else {
			values.resize(2) ;
		}
		values[0] = (v.first) ;
		if(values.size()>1)
			values[1] = (v.second) ;
	} else {
		values.resize(0) ;
	}
}

Properties::Properties(PropertiesType p, const std::vector<double> & v)
{
	pType = p ;
	nVal = standardNVal(p) ;
	if(nVal+1 == 0)
		nVal = v.size() ;
	values.resize(nVal) ;
	for(size_t i = 0 ; i < std::min(v.size(),nVal) ; i++)
		values[i] = (v[i]) ;
	nVal = values.size() ;
}

Properties::Properties(PropertiesType p, const Vector & v)
{
	pType = p ;
	nVal = standardNVal(p) ;
	if(nVal+1 == 0)
		nVal = v.size() ;
	values.resize(nVal) ;
	for(size_t i = 0 ; i < std::min(v.size(),nVal) ; i++)
		values[i] = (v[i]) ;
	nVal = values.size() ;
}

Properties::Properties(PropertiesType p, const Matrix & m)
{
	pType = p ;
	nVal = standardNVal(p) ;
	if(nVal+1 == 0)
	{
		nVal = m.array().size() ;
	}
	values.resize(nVal) ;
	for(size_t i = 0 ; i < values.size() ; i++)
		values[i] = (m.array())[i] ;
	if(p == HOOKE || p == BULK_SHEAR)
	{
		double E = 0 ;
		double nu = 0 ;
		if(m.size() == 9)
		{
			nu = m[0][0]/m[0][1] ;
			E = m[0][0] * (1-nu*nu) ;
		}
		if(m.size() == 36) ;
		{
			nu = m[0][1] / (m[0][0]+m[0][1]) ;
			E = m[0][0] * ((1.+nu)*(1.-2.*nu)) / (1.-nu) ;
		}
		switch(p)
		{
			case HOOKE:
			{
				values[0] = E ;
				values[1] = nu ;
				break;
			}
			case BULK_SHEAR:
			{
				values[0] = E / (3*(1-2*nu)) ;
				values[1] = E / (2*(1+nu)) ;
				break;
			}
		}
	}
}


Properties::Properties(const Properties & p)
{
	pType = p.getPropertiesType() ;
	values.resize(p.getValues().size()) ;
	for(size_t i = 0 ; i < p.getNVal() ; i++)
		values[i]=(p.getValue(i)) ;
	nVal = values.size() ;
}

std::pair<bool, Properties> Properties::convert(PropertiesType p_out)
{
	switch(p_out)
	{
		case VOID_PROP:
			return std::make_pair(true, Properties()) ;
		case ABSTRACT:
			return std::make_pair(true, Properties(values)) ;
		case FRACTION:
		{
			switch(pType)
			{
				case VOID_PROP:
					return std::make_pair(false, Properties()) ;
				case ABSTRACT:
				{
					if(nVal > 0)
						return std::make_pair(true, Properties(FRACTION, values[0])) ;
					else
						return std::make_pair(false, Properties()) ;
				}
				case FRACTION:
					return std::make_pair(true, Properties(*this)) ;
				case HOOKE:
					return std::make_pair(false, Properties()) ;
				case BULK_SHEAR:
					return std::make_pair(false, Properties()) ;
				case EXPANSION:
					return std::make_pair(false, Properties()) ;
			}
		}
		case HOOKE:
		{
			switch(pType)
			{
				case VOID_PROP:
					return std::make_pair(false, Properties()) ;
				case ABSTRACT:
				{
					if(nVal > 1)
						return std::make_pair(true, Properties(HOOKE, std::make_pair(values[0],values[1]))) ;
					else
						return std::make_pair(false, Properties()) ;
				}
				case FRACTION:
					return std::make_pair(false, Properties()) ;
				case HOOKE:
					return std::make_pair(true, Properties(*this)) ;
				case BULK_SHEAR:
				{
					double k = values[0] ;
					double mu = values[1] ;
					double E = 9*k*mu / (3*k+mu) ;
					double nu = (3*k-2*mu) / (6*k+2*mu) ;
					return std::make_pair(true, Properties(HOOKE, std::make_pair(E,nu))) ;
				}	
				case EXPANSION:
					return std::make_pair(false, Properties()) ;
			}
		}
		case BULK_SHEAR:
		{
			switch(pType)
			{
				case VOID_PROP:
					return std::make_pair(false, Properties()) ;
				case ABSTRACT:
				{
					if(nVal > 1)
						return std::make_pair(true, Properties(BULK_SHEAR, std::make_pair(values[0],values[1]))) ;
					else
						return std::make_pair(false, Properties()) ;
				}
				case FRACTION:
					return std::make_pair(false, Properties()) ;
				case HOOKE:
				{
					double E = values[0] ;
					double nu = values[1] ;
					double k = E / (3*(1-2*nu)) ;
					double mu = E / (2*(1+nu)) ;
					return std::make_pair(true, Properties(BULK_SHEAR, std::make_pair(k,mu))) ;
				}
				case BULK_SHEAR:
					return std::make_pair(true, Properties(*this)) ;
				case EXPANSION:
					return std::make_pair(false, Properties()) ;
			}
		}
	}
}

std::pair<bool,Matrix> Properties::getCauchyGreen(SpaceDimensionality dim) const
{
	if(pType == HOOKE || pType == BULK_SHEAR)
	{
		Matrix cg = cauchyGreen(std::make_pair(values[0],values[1]),pType == HOOKE, dim) ;
		return std::make_pair(true,cg) ;
	}
	Matrix cg = cauchyGreen(std::make_pair(1e-9,0.2),true,dim) ;
	return std::make_pair(false,cg) ;
}


void Properties::print()
{
	switch(pType)
	{
		case VOID_PROP:
		{
			std::cout << "VOID PROPERTIES" << std::endl ;
			break ;
		}
		case ABSTRACT:
		{
			std::cout << "ABSTRACT" << std::endl ;
			break ;
		}
		case FRACTION:
		{
			std::cout << "FRACTION" << std::endl ;
			break ;
		}
		case HOOKE:
		{
			std::cout << "HOOKE" << std::endl ;
			break ;
		}
		case BULK_SHEAR:
		{
			std::cout << "BULK SHEAR" << std::endl ;
			break ;
		}
		case EXPANSION:
		{
			std::cout << "EXPANSION COEFFICIENT" << std::endl ;
			break;
		}
	}
	for(size_t i = 0 ; i < nVal ; i++)
		std::cout << values[i] << ";" ; 
	
	std::cout << std::endl ;
}

Material::Material()
{
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

size_t Material::getFirstIndex(PropertiesType p) const
{
	for(size_t i = 0 ; i < this->size() ; i++)
	{
		if((*this)[i].getPropertiesType() == p)
			return i ;
	}
	return -1 ;
}


size_t Material::getHooke()
{
	size_t i = getFirstIndex(HOOKE) ;
	if((i+1) > 0)
		return i ;
	i = getFirstIndex(BULK_SHEAR) ;
	if((i+1) > 0)
	{
		Properties bulk((*this)[i]) ;
		std::pair<bool, Properties> hooke = bulk.convert(HOOKE) ;
		if(hooke.first)
		{
			this->push_back(hooke.second) ;
			return (this->size() -1) ;
		} else {
			return -1 ;
		}
	} else {
		return -1 ;
	}
}

size_t Material::getBulkShear()
{
	size_t i = getFirstIndex(BULK_SHEAR) ;
	if((i+1) > 0)
		return i ;
	i = getFirstIndex(HOOKE) ;
	if((i+1) > 0)
	{
		Properties hooke((*this)[i]) ;
		std::pair<bool, Properties> bulk = hooke.convert(BULK_SHEAR) ;
		if(bulk.first)
		{
			this->push_back(bulk.second) ;
			return (this->size() -1) ;
		} else {
			return -1 ;
		}
	} else {
		return -1 ;
	}
}

bool Material::equals(Material m,PropertiesType p) const
{
	size_t it = this->getFirstIndex(p) ;
	size_t im = m.getFirstIndex(p) ;
//	std::cout << it << ";" << im << std::endl ;
	if(it+1 > 0 && im+1 > 0)
	{
//		std::cout << "here" << std::endl ;
		if((*this)[it].getNVal() == m[im].getNVal())
		{
//			std::cout << "here" << std::endl ;
			bool v = true ;
			for(size_t i = 0 ; i < m[im].getNVal() ; i++)
				v &= (std::abs((*this)[it].getValue(i) - m[im].getValue(i)) < 1e-6) ;
			return v ;
		}
	}
	return false ;
}

void Material::combine(Material m, PropertiesType p)
{
	for(size_t i = 0 ; i < this->size() ; i++)
	{
		if((*this)[i].getPropertiesType() == p)
		{
			size_t im = getFirstIndex(p) ;
			for(size_t j = 0 ; j < (*this)[i].getNVal() ; j++)
				(*this)[i].setValue(j,m[im].getValue(j)+(*this)[i].getValue(i)) ;
		}
	}
}














/*

SimpleMaterial::SimpleMaterial(double E, double nu)
{
	Young_Poisson = std::make_pair(E,nu) ;
}

SimpleMaterial::SimpleMaterial(std::pair<double,double> E_nu)
{
	Young_Poisson = E_nu ;
}

SimpleMaterial::SimpleMaterial(SimpleMaterial & mat)
{
	Young_Poisson = mat.getEnu() ;
}

SimpleMaterial::SimpleMaterial(SimpleMaterial * mat)
{
	Young_Poisson = mat->getEnu() ;
}

SimpleMaterial::SimpleMaterial(HomogenizationScheme scheme, std::pair<double,SimpleMaterial *> inclusions, SimpleMaterial * matrix)
{
	std::pair<double,double> k_mu_inc = inclusions.second->getkmu() ;
	double f_inc = inclusions.first ;
	double k_inc = k_mu_inc.first ;
	double mu_inc = k_mu_inc.second ;

	std::pair<double,double> k_mu_mat = matrix->getkmu() ;
	double f_mat = 1 - f_inc ;
	double k_mat = k_mu_mat.first ;
	double mu_mat = k_mu_mat.second ;

	double k_hom = k_mat ;
	double mu_hom = mu_mat ;

	switch(scheme)
	{
		case DILUTED:
		{
			k_hom += f_inc * (k_inc - k_mat) * (3*k_mat + 4*mu_mat) / (3*k_inc + 4*mu_mat) ;
			mu_hom += f_inc * (mu_inc - mu_mat) * (5*mu_mat*(3*k_mat + 4*mu_mat)) / (mu_mat*(9*k_mat + 8*mu_mat) + 6*mu_inc*(k_mat + 2*mu_mat)) ;
			break;
		}
		case INCREMENTAL:
		{
			double alpha = 1e-6 ;
			double dalpha = 1e-6 ;
			SimpleMaterial * incremental = new SimpleMaterial(matrix) ;
			while(alpha < f_inc)
			{
				incremental = new SimpleMaterial(DILUTED, std::make_pair(dalpha,inclusions.second), incremental) ;
				alpha += dalpha ;
			}
			incremental = new SimpleMaterial(DILUTED, std::make_pair(f_inc - (alpha-dalpha), inclusions.second), incremental) ;
			std::pair<double,double> k_mu_hom = incremental->getkmu() ;
			k_hom = k_mu_hom.first ;
			mu_hom = k_mu_hom.second ;
			break;
		}
		case MORI_TANAKA:
		{
			double Sk = (4*mu_mat) / (3*k_mat) ;
			double Smu = (3 + 2*Sk) / (2 + 3*Sk) ;
			double K = k_inc / k_mat ;
			double Mu = mu_inc / mu_mat ;
			k_hom *= (K/Sk + 1 + f_inc*(K-1)) ;
			k_hom /= (1 + (K + f_inc*(1-K))/Sk) ;
			mu_hom *= (Mu/Smu + 1 + f_inc*(Mu-1)) ;
			mu_hom /= (1 + (Mu + f_inc*(1-Mu))/Smu) ;
			break;
		}
		case SELF_CONSISTENT:
		{
			double k_min = std::min(k_inc,k_mat) ;
			double mu_min = std::min(mu_inc,mu_mat) ;
			SimpleMaterial * sc_mat = new SimpleMaterial(kmu2Enu(std::make_pair(k_min,mu_min))) ;
			std::vector<std::pair<double, SimpleMaterial *> > sc_inc ;
			sc_inc.push_back(std::make_pair(f_inc,inclusions.second)) ;
			sc_inc.push_back(std::make_pair(f_mat,matrix)) ;
			SimpleMaterial * next_sc_mat = new SimpleMaterial(MORI_TANAKA,sc_inc,sc_mat) ;
			while(next_sc_mat->relativeDifference(sc_mat) > 1e-4)
			{
				sc_mat = new SimpleMaterial(next_sc_mat) ;
				next_sc_mat = new SimpleMaterial(MORI_TANAKA,sc_inc,sc_mat) ;
			}
			std::pair<double,double> k_mu_hom = next_sc_mat->getkmu() ;
			k_hom = k_mu_hom.first ;
			mu_hom = k_mu_hom.second ;
			break;
		}
	}
	Young_Poisson = kmu2Enu(std::make_pair(k_hom,mu_hom)) ;
}

SimpleMaterial::SimpleMaterial(HomogenizationScheme scheme, std::vector<std::pair<double,SimpleMaterial *> > inclusions, SimpleMaterial * matrix)
{
	SimpleMaterial * homogenized = new SimpleMaterial(matrix) ;
	switch(scheme)
	{
		case DILUTED:
		{
			for(size_t i = 0 ; i < inclusions.size() ; i++)
				homogenized = new SimpleMaterial(DILUTED, inclusions[i], homogenized) ;
			break;
		}
		case INCREMENTAL:
		{
			for(size_t i = 0 ; i < inclusions.size() ; i++)
				homogenized = new SimpleMaterial(INCREMENTAL, inclusions[i], homogenized) ;
			break;
		}
		case MORI_TANAKA:
		{
			double f_mat = 1 ;
			for(size_t i = 0 ; i < inclusions.size() ; i++)
				f_mat -= inclusions[i].first ;
			if(f_mat < 1e-9)
				f_mat = 1e-9 ;
			homogenized->multiply(f_mat) ;
			for(size_t i = 0 ; i < inclusions.size() ; i++)
			{
				SimpleMaterial * mt = new SimpleMaterial(MORI_TANAKA, inclusions[i], matrix) ;
				mt->multiply(inclusions[i].first) ;
				homogenized->add(mt) ;
			}
			break;
		}
		case SELF_CONSISTENT:
		{
			double k_min = matrix->getkmu().first ;
			double mu_min = matrix->getkmu().second ;
			double f_mat = 1 ;
			for(size_t i = 0 ; i < inclusions.size() ; i++)
			{
				k_min = std::max(k_min, inclusions[i].second->getkmu().first) ;
				mu_min = std::max(mu_min, inclusions[i].second->getkmu().second) ;
				f_mat -= inclusions[i].first ;
			}
			SimpleMaterial * sc_mat = new SimpleMaterial(kmu2Enu(std::make_pair(k_min,mu_min))) ;
			std::vector<std::pair<double, SimpleMaterial *> > sc_inc ;
			for(size_t i = 0 ; i < inclusions.size() ; i++)
				sc_inc.push_back(inclusions[i]) ;
			sc_inc.push_back(std::make_pair(f_mat,matrix)) ;
			homogenized = new SimpleMaterial(MORI_TANAKA,sc_inc,sc_mat) ;
			while(homogenized->relativeDifference(sc_mat) > 1e-4)
			{
				sc_mat = new SimpleMaterial(homogenized) ;
				homogenized = new SimpleMaterial(MORI_TANAKA,sc_inc,sc_mat) ;
			}
			break;
		}
	}
	Young_Poisson = homogenized->getEnu() ;
}

SimpleMaterial::SimpleMaterial(XMLTree * xml)
{
	std::pair<bool,double> E ;
	E.first = false ;
	std::pair<bool,double> nu ;
	nu.first = false ;
	if(xml->match("material"))
	{
		if(xml->nChildren() > 1)
		{
			if(xml->getChild(0)->match("young"))
				E = xml->getChild(0)->buildDouble() ;
			if(xml->getChild(1)->match("poisson"))
				nu = xml->getChild(0)->buildDouble() ;
		}
	}
	if(E.first && nu.first)
		Young_Poisson = std::make_pair(E.second,nu.second) ;
	else
		Young_Poisson = kmu2Enu(std::make_pair(0.,0.)) ;
}

std::pair<double, double> SimpleMaterial::getkmu() 
{
	std::pair<double, double> kmu = Enu2kmu(Young_Poisson) ; 
	return kmu ;
}


VOID_PROP SimpleMaterial::add(SimpleMaterial * mat)
{
	double k = this->getkmu().first + mat->getkmu().first ;
	double mu = this->getkmu().second + mat->getkmu().second ;
	Young_Poisson = kmu2Enu(std::make_pair(k,mu)) ;
}
VOID_PROP SimpleMaterial::add(SimpleMaterial * mat, double f)
{
	double k = this->getkmu().first + f * mat->getkmu().first ;
	double mu = this->getkmu().second + f * mat->getkmu().second ;
	Young_Poisson = kmu2Enu(std::make_pair(k,mu)) ;
}

VOID_PROP SimpleMaterial::multiply(double f) 
{
	double k = f * this->getkmu().first ;
	double mu = f * this->getkmu().second ;
	Young_Poisson = kmu2Enu(std::make_pair(k,mu)) ;
}

double SimpleMaterial::relativeDifference(SimpleMaterial * mat)
{
	double diff = (Young_Poisson.first - mat->getEnu().first) / Young_Poisson.first ;
	if(diff < 0)
		diff *= -1 ;
	
	return diff ;
}

XMLTree * SimpleMaterial::toXML()
{
	XMLTree * mat = new XMLTree("material") ;
	mat->addChild(new XMLTree("young",Young_Poisson.first)) ;
	mat->addChild(new XMLTree("poisson",Young_Poisson.second)) ;
	return mat ;
}


std::pair<double,double> Enu2kmu(std::pair<double,double> E_nu)
{
	double E = E_nu.first ;
	double nu = E_nu.second ;
	double k = E / (3*(1-2*nu)) ;
	double mu = E / (2*(1+nu)) ;
	return std::make_pair(k,mu) ;
} ;

std::pair<double,double> kmu2Enu(std::pair<double,double> k_mu)
{
	double k = k_mu.first ;
	double mu = k_mu.second ;
	double E = 9*k*mu / (3*k+mu) ;
	double nu = (3*k-2*mu) / (6*k+2*mu) ;
	return std::make_pair(E,nu) ;
} ;*/




