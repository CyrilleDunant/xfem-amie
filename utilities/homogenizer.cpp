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

#include "homogenizer.h"

namespace Mu
{



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


void SimpleMaterial::add(SimpleMaterial * mat)
{
	double k = this->getkmu().first + mat->getkmu().first ;
	double mu = this->getkmu().second + mat->getkmu().second ;
	Young_Poisson = kmu2Enu(std::make_pair(k,mu)) ;
}
void SimpleMaterial::add(SimpleMaterial * mat, double f)
{
	double k = this->getkmu().first + f * mat->getkmu().first ;
	double mu = this->getkmu().second + f * mat->getkmu().second ;
	Young_Poisson = kmu2Enu(std::make_pair(k,mu)) ;
}

void SimpleMaterial::multiply(double f) 
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
} ;

}


