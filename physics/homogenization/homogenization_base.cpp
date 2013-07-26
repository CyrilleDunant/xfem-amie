#include "homogenization_base.h"
#include "scheme_template.h"

using namespace Mu ;

Properties::Properties()
{
	tag = P_nullptr_PROPERTIES ;
	value = 0. ;
}

Properties::Properties(PType t, double v)
{
	tag = t ;
	value = v ;
}

Properties::Properties(const Properties & p)
{
	tag = p.type() ;
	value = p.val() ;
}

PType Properties::type() const
{
	return tag ;
}

bool Properties::is(PType t) const
{
	return tag == t ;
}

bool Properties::isExtensible() const
{
	return is(P_VOLUME) || is(P_MASS) ;
}

double Properties::val() const
{
	return value ;
}

double Properties::val(PType t) const
{
	if(is(t))
	{
		return value ;
	}
	else
	{
		return 0. ;
	}
}

void Properties::set(double v)
{
	value = v ;
}

bool Properties::equals(Properties p, bool extensible) const
{
	if(extensible && isExtensible())
	{
		return true ;
	}
	else
	{
		if(p.is(tag))
		{
			return p.val() == value ;
		}
		else
		{
			return false ;
		}
	}
}






Material::Material()
{
	prop.clear() ;
	phases.clear() ;
}

Material::Material(std::vector<Properties> p)
{
	prop.clear() ;
	for(size_t i = 0 ; i < p.size() ; i++)
	{
		prop.push_back(p[i]) ;
	}
	phases.clear() ;
}

Material::Material(std::vector<Material> p)
{
	prop.clear() ;
	phases.clear() ;
	for(size_t i = 0 ; i < p.size() ; i++)
	{
		phases.push_back(p[i]) ;
	}
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

	setProperties(P_YOUNG_MODULUS, E) ;
	setProperties(P_POISSON_RATIO, nu) ;

	searchProperties(P_BULK_MODULUS, SEARCH_TRY_ANYTHING, false) ;
	searchProperties(P_SHEAR_MODULUS, SEARCH_TRY_ANYTHING, false) ;
}


void Material::addProperties(Properties p)
{
	prop.push_back(p) ;
}

void Material::addProperties(PType t, double v) 
{
	prop.push_back(Properties(t,v)) ;
}

void Material::aloneProperties(bool last)
{
	int i = 0 ;
	while(i < (int) prop.size())
	{
		aloneProperties(prop[i].type(),last) ;
		i++ ;
	}
}

void Material::aloneProperties(PType t, bool last)
{
	while(indexOfProperties(t, true) != indexOfProperties(t, false))
	{
		removeProperties(indexOfProperties(t,!last)) ;
	}
}

void Material::clearProperties()
{
	prop.clear() ;
}

Properties Material::getProperties(int i) const
{
	if(i < 0 || i > (int) prop.size())
	{
		return Properties(P_BAD_INDEX, 0.) ;
	}
	else
	{
		return prop[i] ;
	}
}

Properties Material::getProperties(PType t) const 
{
	return getProperties(indexOfProperties(t)) ;
}

bool Material::hasProperties(PType t) const 
{
	return indexOfProperties(t) > -1 ;
}

int Material::indexOfProperties(PType t, bool first) const 
{
	if(first)
	{
		for(int i = 0 ; i < (int) prop.size() ; i++)
		{
			if(prop[i].is(t))
			{
				return i ;
			}
		}
	}
	else
	{
		for(int i = (int) prop.size() ; i > 0 ; i--)
		{
			if(prop[i-1].is(t))
			{
				return i-1 ;
			}
		}
	}
	return -1 ;
}

void Material::removeProperties(int i) 
{
	if(!(i < 0 || i > (int) prop.size()))
	{
		std::vector<Properties> tmp ;
		for(int j = 0 ; j < (int) prop.size() ; j++)
		{
			if(j != i)
			{
				tmp.push_back(prop[j]) ;
			}
		}
		prop.clear() ;
		for(int j = 0 ; j < (int) tmp.size() ; j++)
		{
			prop.push_back(tmp[j]) ;
		}
	}
}

void Material::removeProperties(PType t) 
{
	removeProperties(indexOfProperties(t)) ;
}

void Material::setProperties(PType t, double v) 
{
	int i = indexOfProperties(t) ;
	if(i != -1)
	{
		prop[i].set(v) ;
	}
	else
	{
		prop.push_back(Properties(t, v)) ;
	}
}

void Material::setProperties(Properties p) 
{
	int i = indexOfProperties(p.type()) ;
	if(i != -1)
	{
		prop[i].set(p.val()) ;
	}
	else
	{
		prop.push_back(p) ;
	}
}

void Material::setProperties(Material m)
{
	for(int i = 0 ; i < m.sizeProperties() ; i++)
	{
		setProperties(m.getProperties(i)) ;
	}
}

int Material::sizeProperties() const
{
	return (int) prop.size() ;
}

double Material::valProperties(int i) const
{
	if(!(i < 0 || i > (int) prop.size()))
	{
		return prop[i].val() ;
	}
	else
	{
		return 0. ;
	}
}

double Material::valProperties(PType t) const
{
	return valProperties(indexOfProperties(t)) ;
}

void Material::addPhase(Material m, bool merge) 
{
	phases.push_back(m) ;
	if(merge)
	{
		mergePhase() ;
	}
}

void Material::clearPhase()
{
	while((int) phases.size() > 0)
	{
		removePhase(0) ;
	}
}

Material Material::getPhase(int i) const
{
	if(i < 0 || i > (int) phases.size()-1)
	{
		return Material() ;
	}
	else
	{
		return phases[i] ;
	}
}

void Material::mergePhase()
{
	for(int i = 0 ; i < (int) phases.size()-1 ; i++)
	{
		for(int j = i+1 ; j < (int) phases.size() ; j++)
		{
			if(phases[i].equals(phases[j], true))
			{
				for(int k = 0 ; k < phases[i].sizeProperties() ; k++)
				{
					if(phases[i].getProperties(k).isExtensible())
					{
						Properties extension = phases[i].getProperties(k) ;
						phases[i].setProperties(extension.type(),extension.val() + phases[j].valProperties(extension.type())) ;
					}
				}
				removePhase(j) ;
			}
		}
	}
}

void Material::removePhase(int i) 
{
	if(!(i < 0 || i > (int) phases.size()))
	{
		std::vector<Material> tmp ;
		for(int j = 0 ; j < (int) phases.size() ; j++)
		{
			if(j != i)
			{
				tmp.push_back(phases[j]) ;
			}
		}
		phases.clear() ;
		for(int j = 0 ; j < (int) tmp.size() ; j++)
		{
			phases.push_back(tmp[j]) ;
		}
	}
}

int Material::sizePhase() const
{
	return (int) phases.size() ;
}

double Material::valPhase(int i, PType t) const
{
	return getPhase(i).valProperties(t) ;
}

std::vector<double> Material::valPhase(PType t) const 
{
	std::vector<double> val ;
	for(int i = 0 ; i < (int) phases.size() ; i++)
	{
		val.push_back(valPhase(i,t)) ;
	}
	return val ;
}

bool Material::equals(Material m, bool extensible) const
{
	for(int i = 0 ; i < (int) prop.size() ; i++)
	{
		if(!prop[i].equals(m.getProperties(prop[i].type()), extensible))
		{
			return false ;
		}
	}
	return true ;
}

Properties Material::searchProperties(PType t, SearchPattern p, bool next, Material father)
{
	if(!hasProperties(t) && p != SEARCH_TRY_NOTHING)
	{
	
		SearchPattern pattern = SEARCH_TRY_NOTHING ;
		if(next)
		{
			pattern = SEARCH_TRY_ANYTHING ;
		}

		switch(t)
		{
		case P_BULK_MODULUS:
			if(hasProperties(P_POISSON_RATIO) && hasProperties(P_YOUNG_MODULUS))
			{
				double E = valProperties(P_YOUNG_MODULUS) ;
				double nu = valProperties(P_POISSON_RATIO) ;
				setProperties(P_BULK_MODULUS,E/(3.*(1.-2.*nu))) ;
			}
			break ;

		case P_POISSON_RATIO:
			if(hasProperties(P_BULK_MODULUS) && hasProperties(P_SHEAR_MODULUS))
			{
				double k = valProperties(P_BULK_MODULUS) ;
				double mu = valProperties(P_SHEAR_MODULUS) ;
				setProperties(P_POISSON_RATIO,(3.*k-2.*mu)/(6.*k+2.*mu)) ;
			}
			break ;

		case P_SHEAR_MODULUS:
			if(hasProperties(P_POISSON_RATIO) && hasProperties(P_YOUNG_MODULUS))
			{
				double E = valProperties(P_YOUNG_MODULUS) ;
				double nu = valProperties(P_POISSON_RATIO) ;
				setProperties(P_SHEAR_MODULUS,E/(2.*(1.+nu))) ;
			}
			break ;
		
		case P_YOUNG_MODULUS:
			if(hasProperties(P_BULK_MODULUS) && hasProperties(P_SHEAR_MODULUS))
			{
				double k = valProperties(P_BULK_MODULUS) ;
				double mu = valProperties(P_SHEAR_MODULUS) ;
				setProperties(P_YOUNG_MODULUS,(9.*k*mu)/(3.*k+mu)) ;
			}
			break ;

		case P_EXPANSION_COEFFICIENT:
			setProperties(P_EXPANSION_COEFFICIENT, 0.) ;
			break ;

		case P_VOLUME:
			switch(p)
			{
			case SEARCH_TRY_ANYTHING:
				searchProperties(P_VOLUME, SEARCH_VOLUME_FROM_MASS_DENSITY, false) ;
				searchProperties(P_VOLUME, SEARCH_VOLUME_FROM_FATHER_VOLUME, false, father) ;
				searchProperties(P_VOLUME, SEARCH_VOLUME_FROM_CHILDREN_VOLUME, false) ;
				break ;
			
			case SEARCH_VOLUME_FROM_MASS_DENSITY:
				if(searchProperties(P_MASS, pattern, false).is(P_MASS) 
					&& searchProperties(P_DENSITY, pattern, false).is(P_DENSITY))
				{
					double mass = valProperties(P_MASS) ;
					double d = valProperties(P_DENSITY) ;
					setProperties(P_VOLUME, mass/d) ;
				}
				break ;
			
			case SEARCH_VOLUME_FROM_FATHER_VOLUME:
				if(searchProperties(P_VOLUME_FRACTION, pattern, false).is(P_VOLUME_FRACTION) 
					&& father.searchProperties(P_VOLUME, pattern, false).is(P_VOLUME))
				{
					double f = valProperties(P_VOLUME_FRACTION) ;
					double tot = father.valProperties(P_VOLUME) ;
					setProperties(P_VOLUME, f*tot) ;
				}
				break ;
				
			case SEARCH_VOLUME_FROM_CHILDREN_VOLUME:
				bool has = (sizePhase() > 0) ;
				if(has)
				{
					double v = 0 ;
					int i = 0 ;
					while(has && i < sizePhase())
					{
						has &= phases[i].searchProperties(P_VOLUME, pattern, next).is(P_VOLUME) ;
						v += phases[i].valProperties(P_VOLUME) ;
						i++ ;
					}
					if(has)
					{
						setProperties(P_VOLUME, v) ;
					}
				}
				break ;
			
			}
			break ;

		case P_VOLUME_FRACTION:
			if(searchProperties(P_VOLUME, pattern, false).is(P_VOLUME)
				 && father.searchProperties(P_VOLUME, SEARCH_VOLUME_FROM_CHILDREN_VOLUME, true).is(P_VOLUME))
			{
				double tot = father.valProperties(P_VOLUME) ;
				double vol = valProperties(P_VOLUME) ;
				setProperties(P_VOLUME_FRACTION, vol/tot) ;
			}
			break ;

		case P_MASS:
			switch(p)
			{
			case SEARCH_TRY_ANYTHING:
				searchProperties(P_MASS, SEARCH_MASS_FROM_VOLUME_DENSITY, false) ;
				searchProperties(P_MASS, SEARCH_MASS_FROM_FATHER_MASS, false, father) ;
				searchProperties(P_MASS, SEARCH_MASS_FROM_CHILDREN_MASS, false) ;
				break ;
			
			case SEARCH_MASS_FROM_VOLUME_DENSITY:
				if(searchProperties(P_VOLUME, pattern, false).is(P_VOLUME)
					 && searchProperties(P_DENSITY, pattern, false).is(P_DENSITY))
				{
					double vol = valProperties(P_VOLUME) ;
					double d = valProperties(P_DENSITY) ;
					setProperties(P_MASS, vol*d) ;
				}
				break ;
			
			case SEARCH_MASS_FROM_FATHER_MASS:
				if(searchProperties(P_MASS_FRACTION, pattern, false).is(P_MASS_FRACTION)
					 && father.searchProperties(P_MASS, pattern, false).is(P_MASS))
				{
					double f = valProperties(P_MASS_FRACTION) ;
					double tot = father.valProperties(P_MASS) ;
					setProperties(P_MASS, f*tot) ;
				}
				break ;
				
			case SEARCH_MASS_FROM_CHILDREN_MASS:
				bool has = (sizePhase() > 0) ;
				if(has)
				{
					double m = 0 ;
					int i = 0 ;
					while(has && i < sizePhase())
					{
						has &= phases[i].searchProperties(P_MASS, pattern, false).is(P_MASS) ;
						m += phases[i].valProperties(P_MASS) ;
						i++ ;
					}
					if(has)
					{
						setProperties(P_MASS, m) ;
					}
				}
				break ;
			
			}
			break ;

		case P_MASS_FRACTION:
			if(searchProperties(P_MASS, pattern, false).is(P_MASS)
				 && father.searchProperties(P_MASS, SEARCH_MASS_FROM_CHILDREN_MASS, true).is(P_MASS))
			{
				double tot = father.valProperties(P_MASS) ;
				double mass = valProperties(P_MASS) ;
				setProperties(P_MASS_FRACTION, mass/tot) ;
			}
			break ;
		
		case P_DENSITY:
			if(searchProperties(P_MASS, pattern, false).is(P_MASS)
				 && searchProperties(P_VOLUME, pattern, false).is(P_VOLUME))
			{
				double mass = valProperties(P_MASS) ;
				double vol = valProperties(P_VOLUME) ;
				setProperties(P_DENSITY, mass/vol) ;
			}
			break ;
		
		}

	}
	return getProperties(t) ;
}

bool Material::homogenize(HomogenizationScheme s)
{
	SchemeTemplate scheme(s) ;
	if(scheme.cast(*this))
	{
		for(int i = 0 ; i < scheme.sizeResult() ; i++)
		{
			setProperties(scheme.getResult(i)) ;
		}
		return true ;
	}
	else
	{
		return false ;
	}
}

Matrix Material::cauchyGreen(std::pair<double,double> prop, bool hooke, SpaceDimensionality dim, planeType pt) 
{
	double E = prop.first ;
	double nu = prop.second ;
	if(!hooke)
	{
		double k = prop.first ;
		double mu = prop.second ;
		E = 9*k*mu / (3*k+mu) ;
		nu = (3*k-2*mu) / (6*k+2*mu) ;
	}
	switch(dim)
	{
		case SPACE_TWO_DIMENSIONAL:
		{
			Matrix cg(3,3) ;

			if(pt == PLANE_STRESS)
			{
			cg[0][0] = 1. ; cg[0][1] = nu ; cg[0][2] = 0 ;
			cg[1][0] = nu ; cg[1][1] = 1. ; cg[1][2] = 0 ;
			cg[2][0] = 0 ; cg[2][1] = 0 ; cg[2][2] = (1.-nu) ;
			cg *= E/(1.-nu*nu) ;
			}
			else
			{
			cg[0][0] = 1.-nu ; cg[0][1] = nu ; cg[0][2] = 0 ;
			cg[1][0] = nu ; cg[1][1] = 1.-nu ; cg[1][2] = 0 ;
			cg[2][0] = 0 ; cg[2][1] = 0 ; cg[2][2] = (1-2.*nu) ;
			cg *= E/((1.+nu)*(1.-2.*nu)) ;
			}
			return cg ;
		}
		case SPACE_THREE_DIMENSIONAL:
		{
			Matrix cgg(6,6) ;
			cgg[0][0] = 1. - nu ; cgg[0][1] = nu ; cgg[0][2] = nu ;
			cgg[1][0] = nu ; cgg[1][1] = 1. - nu ; cgg[1][2] = nu ;
			cgg[2][0] = nu ; cgg[2][1] = nu ; cgg[2][2] = 1. - nu ;
			cgg[3][3] = (0.5 - nu)*.5 ;
			cgg[4][4] = (0.5 - nu)*.5 ;
			cgg[5][5] = (0.5 - nu)*.5 ;
			cgg *= E/((1.+nu)*(1.-2.*nu)) ;
			return cgg ;
		}
	}
	return Matrix(0,0) ;
}

Matrix Material::cauchyGreen(double p1, double p2, bool hooke, SpaceDimensionality dim, planeType pt) 
{
	double E = p1 ;
	double nu = p2 ;
	if(!hooke)
	{
		double k = p1 ;
		double mu = p2 ;
		E = 9*k*mu / (3*k+mu) ;
		nu = (3*k-2*mu) / (6*k+2*mu) ;
	}
	switch(dim)
	{
		case SPACE_TWO_DIMENSIONAL:
		{
			Matrix cg(3,3) ;

			if(pt == PLANE_STRESS)
			{
			cg[0][0] = 1. ; cg[0][1] = nu ; cg[0][2] = 0 ;
			cg[1][0] = nu ; cg[1][1] = 1. ; cg[1][2] = 0 ;
			cg[2][0] = 0 ; cg[2][1] = 0 ; cg[2][2] = (1.-nu) ;
			cg *= E/(1.-nu*nu) ;
			}
			else
			{
			cg[0][0] = 1.-nu ; cg[0][1] = nu ; cg[0][2] = 0 ;
			cg[1][0] = nu ; cg[1][1] = 1.-nu ; cg[1][2] = 0 ;
			cg[2][0] = 0 ; cg[2][1] = 0 ; cg[2][2] = (1-2.*nu) ;
			cg *= E/((1.+nu)*(1.-2.*nu)) ;
			}
			return cg ;
		}
		case SPACE_THREE_DIMENSIONAL:
		{
			Matrix cgg(6,6) ;
			cgg[0][0] = 1. - nu ; cgg[0][1] = nu ; cgg[0][2] = nu ;
			cgg[1][0] = nu ; cgg[1][1] = 1. - nu ; cgg[1][2] = nu ;
			cgg[2][0] = nu ; cgg[2][1] = nu ; cgg[2][2] = 1. - nu ;
			cgg[3][3] = 0.5 - nu ;
			cgg[4][4] = 0.5 - nu ;
			cgg[5][5] = 0.5 - nu ;
			cgg *= E/((1.+nu)*(1.-2.*nu)) ;
			return cgg ;
		}
	}
	return Matrix(0,0) ;
}

Matrix Material::orthothropicCauchyGreen(double E_1, double E_2, double G,  double nu, planeType pt) 
{

	Matrix cg(3,3) ;

	if(pt == PLANE_STRESS)
	{
		if(E_1 > POINT_TOLERANCE_2D && E_2 > POINT_TOLERANCE_2D)
		{
// 			Matrix A(2,2) ;
// 			A[0][0] = E_1/std::max(E_1, E_2) ; A[0][1] = -E_2/std::max(E_1, E_2) ;
// 			A[1][0] = 0 ; A[1][1] = 1. ;
// 			Vector b(2) ; b[0] = 0 ; b[1] = nu ;
// 			b = inverse2x2Matrix(A)*b ;
// 			double nu_12 = b[1];
// 			double nu_21 = b[0] ;
			
			//nu_12*nu_21 = nu*nu ;
			double nu_21 = (nu/E_1)*sqrt(E_1*E_2) ;
			double nu_12 = (nu/E_2)*sqrt(E_1*E_2) ;

			double gamma = 1./(1.-nu_12*nu_21) ;

			cg[0][0] = E_1*gamma ; cg[0][1] = nu_21*E_1*gamma ;
			cg[1][0] = cg[0][1] ;  cg[1][1] = E_2*gamma ;
			
			G = E_1*E_2/(E_1*(1.-nu_12*nu_21)+E_2*(1.-nu_21*nu_12)) ;
			
			cg[2][2] = G ;
		}
		else if(E_1 > POINT_TOLERANCE_2D)
		{
			cg[0][0] = E_1 ;
			cg[2][2] = 0 ;

		}
		else if(E_2 > POINT_TOLERANCE_2D)
		{
			cg[1][1] = E_2 ;
			cg[2][2] = 0 ;
		}
		else
			cg.array() = 0 ;
	}
	else if (pt == PLANE_STRAIN)
	{
		if(E_1 > POINT_TOLERANCE_2D && E_2 > POINT_TOLERANCE_2D)
		{
			Matrix A(2,2) ;

			double nu21 = (nu/E_1)*sqrt(E_1*E_2) ;
			double nu12 = (nu/E_2)*sqrt(E_1*E_2) ;
			double nu23 = nu ;
			double nu32 = nu ;
			double nu13 = nu ;
			double nu31 = nu ;
			double nupe = 1.-nu21*nu12-nu23*nu32-nu13*nu31-nu12*nu23*nu31-nu21*nu32*nu13 ;

			cg[0][0] = E_1*(1.-nu32*nu23)/nupe ; cg[0][1] = (nu21-nu23*nu31)*E_1/nupe;
			cg[1][0] = cg[0][1] ;  cg[1][1] = E_2*(1.-nu32*nu23)/nupe ;
			cg[2][2] = E_1*E_2/(E_2*(1+nu12)*(1-2.*nu12)+E_1*(1+nu21)*(1-2.*nu21)) ;
		}
		else if(E_1 > POINT_TOLERANCE_2D)
		{
			cg[0][0] = E_1 ;
			cg[2][2] = G ;
		}
		 else if(E_2 > POINT_TOLERANCE_2D)
		{
			cg[1][1] = E_2 ;
			cg[2][2] = G ;
		}
		else
			cg.array() = 0 ;
	}
	return cg ;
}

Matrix Material::orthothropicCauchyGreen(double E_1, double E_2, double E_3, double G_1, double G_2, double G_3,  double nu) 
{

	Matrix cg(6,6) ;
	double nu_12 = nu ;
	double nu_13 = nu ;
	double nu_23 = nu ;
	double nu_21 = nu_12*E_2/E_1 ;
	double nu_31 = nu_13*E_3/E_1 ;
	double nu_32 = nu_23*E_3/E_2 ;
	double gamma = 1./(1.-nu_12*nu_21-nu_23*nu_32-nu_13*nu_31-2.*nu_12*nu_32*nu_13) ;
	cg[0][0] = E_1*(1.-nu_23*nu_32)*gamma ;
	cg[1][1] = E_2*(1.-nu_13*nu_31)*gamma ;
	cg[2][2] = E_3*(1.-nu_12*nu_21)*gamma ;
	cg[0][1] = E_1*(nu_21+nu_31*nu_23)*gamma ; cg[1][0] = cg[0][1] ;
	cg[0][2] = E_1*(nu_31+nu_21*nu_32)*gamma ; cg[2][0] = cg[0][2] ;
	cg[1][2] = E_2*(nu_32+nu_12*nu_31)*gamma ; cg[2][1] = cg[1][2] ;
	cg[3][3] = G_1*0.5 ;
	cg[4][4] = G_2*0.5 ;
	cg[5][5] = G_3*0.5 ;
	
	return cg ;
}


