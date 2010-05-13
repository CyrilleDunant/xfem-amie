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

#include "scheme_base.h"
#include "../../geometry/geometry_base.h"

using namespace Mu ;



Scheme::Scheme(int n)
{
	s = STATUS_RESET ;
	p = n ;
}

Scheme::Scheme(int n, Tag in)
{
	s = STATUS_RESET ;
	p = n ;
	input.push_back(in) ;
	output.push_back(in) ;
}

Scheme::Scheme(int n, Tag in, Tag out)
{
	s = STATUS_RESET ;
	p = n ;
	input.push_back(in) ;
	output.push_back(out) ;
}

Scheme::Scheme(int n, std::vector<Tag> & in)
{
	s = STATUS_RESET ;
	p = n ;
	for(size_t i = 0 ; i < in.size() ; i++)
	{
		input.push_back(in[i]) ;
		output.push_back(in[i]) ;
	}
}

Scheme::Scheme(int n, std::vector< Tag >& in, std::vector< Tag >& out)
{
	s = STATUS_RESET ;
	p = n ;
	for(size_t i = 0 ; i < in.size() ; i++)
		input.push_back(in[i]) ;
	for(size_t j = 0 ; j < out.size() ; j++)
		output.push_back(out[j]) ;
}

std::vector<Properties> Scheme::homogenize(const std::vector<Material> & mat)
{
	s = STATUS_RESET ;

	std::vector<Properties> hom ;

	size_t n = p ;
	if(n+1 == 0)
		n = mat.size() ;

	if(mat.size() < n)
	{
		s = STATUS_MATERIAL_NOT_FOUND ;
		return hom ;
	}

	Matrix data(n,input.size()) ;

	for(int i = 0 ; i < n ; i++)
	{
		for(int j = 0 ; j < input.size() ; j++)
		{
			size_t k = mat[i].getIndex(input[j],-1) ;
			if(k+1 == 0)
			{
				s = STATUS_PROPERTIES_NOT_FOUND ;
				return hom ;
			}
			data[i][j] = (mat[i])[k].val() ;
		}
	}

	Vector processed = this->process(data) ;

	if(s == STATUS_RESET)
	{
		if(processed.size() != output.size())
		{
			s = STATUS_BAD_HOMOGENIZATION ;
		} else {
			s = STATUS_OK ;
		}
	}

	size_t nout = output.size() ;
	if(processed.size() < output.size())
		nout = processed.size() ;

	for(size_t i = 0 ; i < nout ; i++)
		hom.push_back(Properties(output[i],processed[i])) ;

	return hom ;
}

std::vector<Properties> Scheme::homogenize(const Material & mat)
{
	std::vector<Material> m ;
	m.push_back(mat) ;
	return this->homogenize(m) ;
}

std::vector<Properties> Scheme::homogenize(const Material & m1, const Material & m2)
{
	std::vector<Material> m ;
	m.push_back(m1) ;
	m.push_back(m2) ;
	return this->homogenize(m) ;
}

Vector Scheme::process(const Matrix & data)
{
	std::cout << "no!" << std::endl ;
	Vector processed(output.size()) ;
	return processed ;
}

bool Scheme::check(bool r)
{
	bool b = isOK() ;
	if(r)
		reset() ;
	return b ;
}

bool Scheme::equalsZero(double x)
{
	if(std::abs(x) < POINT_TOLERANCE)
	{
		s = STATUS_BAD_HOMOGENIZATION ;
		return true ;
	}
	return false ;
}

bool Scheme::lessThanZero(double x)
{
	if(x < 0.)
	{
		s = STATUS_BAD_HOMOGENIZATION ;
		return true ;
	}
	return false ;
}

double Scheme::simpleDivision(double num, double denom)
{
	if(equalsZero(denom))
		return 0. ;
	return num / denom ;
}

double Scheme::simpleSquareRoot(double square)
{
	if(lessThanZero(square))
		return 0. ;
	return std::sqrt(square) ;
}

void Scheme::print()
{
	switch(s)
	{
	case STATUS_RESET:
		std::cout << "STATUS RESET" << std::endl ;
		break;
	case STATUS_OK:
		std::cout << "STATUS OK" << std::endl ;
		break;
	case STATUS_MATERIAL_NOT_FOUND:
		std::cout << "MATERIAL_NOT_FOUND" << std::endl ;
		break;
	case STATUS_PROPERTIES_NOT_FOUND:
		std::cout << "PROPERTIES_NOT_FOUND" << std::endl ;
		break;
	case STATUS_BAD_HOMOGENIZATION:
		std::cout << "BAD_HOMOGENIZATION" << std::endl ;
		break;
	}
}

MeanScheme::MeanScheme(bool v, bool p, Tag t) : Scheme(-1)
{
	parallel = p ;

	if(v)
	{
		input.push_back(TAG_VOLUME_FRACTION) ;
	} else {
		input.push_back(TAG_MASS_FRACTION) ;
	}

	input.push_back(t) ;
	output.push_back(t) ;

}

MeanScheme::MeanScheme(bool v, bool p, std::vector<Tag> t) : Scheme(-1)
{
	parallel = p ;

	if(v)
	{
		input.push_back(TAG_VOLUME_FRACTION) ;
	} else {
		input.push_back(TAG_MASS_FRACTION) ;
	}

	for(size_t i = 0 ; i < t.size() ; i++)
	{
		input.push_back(t[i]) ;
		output.push_back(t[i]) ;
	}

}


Vector MeanScheme::process(const Matrix & data)
{
	if(parallel)
		return processParallel(data) ;

	Vector processed(output.size()) ;
	double ftot = 0. ;

	for(size_t i = 0 ; i < data.size() ; i++)
	{
		ftot += data[i][0] ;
		for(size_t j = 0 ; j < processed.size() ; j++)
			processed[j] += simpleDivision(data[i][0], data[i][j+1]) ;
	}

	for(size_t j = 0 ; j < processed.size() ; j++)
		processed[j] = simpleDivision(ftot,processed[j]) ;

	return processed ;
}


Vector MeanScheme::processParallel(const Matrix & data)
{
	Vector processed(output.size()) ;
	double ftot = 0. ;

	for(size_t i = 0 ; i < data.numRows() ; i++)
	{
		ftot += data[i][0] ;
		for(size_t j = 0 ; j < processed.size() ; j++)
			processed[j] += data[i][0] * data[i][j+1] ;
	}

	for(size_t j = 0 ; j < processed.size() ; j++)
		processed[j] = simpleDivision(processed[j],ftot) ;

	return processed ;
}










