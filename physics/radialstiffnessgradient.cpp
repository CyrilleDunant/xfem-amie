//
// C++ Implementation: radialstiffnessgradient
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "radialstiffnessgradient.h"

using namespace Mu ;

RadialStiffnessGradient::RadialStiffnessGradient(double E_int, double nu_int, double rint, double E_ext, double nu_ext, double rext, Point c) : LinearForm(Matrix(3,3), false, true, 2), paramAlt(3,3), r_ext(rext), r_int(rint), dr(r_ext-r_int), centre(c)
{

	param[0][0] = E_int/(1-nu_int*nu_int) ; param[0][1] =E_int/(1-nu_int*nu_int)*nu_int ; param[0][2] = 0 ;
	param[1][0] = E_int/(1-nu_int*nu_int)*nu_int ; param[1][1] = E_int/(1-nu_int*nu_int) ; param[1][2] = 0 ; 
	param[2][0] = 0 ; param[2][1] = 0 ; param[2][2] = E_int/(1-nu_int*nu_int)*(1.-nu_int)/2. ; 
	
	paramAlt[0][0] = E_ext/(1-nu_ext*nu_ext) ;        paramAlt[0][1] = E_ext/(1-nu_ext*nu_ext)*nu_ext ; paramAlt[0][2] = 0 ;
	paramAlt[1][0] = E_ext/(1-nu_ext*nu_ext)*nu_ext ; paramAlt[1][1] = E_ext/(1-nu_ext*nu_ext) ;        paramAlt[1][2] = 0 ; 
	paramAlt[2][0] = 0 ; paramAlt[2][1] = 0 ;         paramAlt[2][2] = E_ext/(1-nu_ext*nu_ext)*(1.-nu_ext)/2. ; 
	this->space_d = true ;
}

RadialStiffnessGradient::~RadialStiffnessGradient() { } ;

void RadialStiffnessGradient::transform(const Function & x, const Function & y)
{
	r = f_sqrt(((x-centre.x)^2)+((y-centre.y)^2)) ;
}

Matrix RadialStiffnessGradient::apply(const Function & p_i, const Function & p_j, const IntegrableEntity *e) const
{
	VirtualMachine vm ;
	
	FunctionMatrix C(3,3) ;
	std::vector<Variable> v ;
	v.push_back(XI);
	v.push_back(ETA);

	for(size_t i = 0 ; i < 3 ; i++)
	{
		for(size_t j = 0 ; j < 3 ; j++)
		{
// 				double val =  (paramAlt[i][j] - param[i][j])/dr ;
			C[i][j] = ((r-r_int)/dr)*paramAlt[i][j] - ((r-r_ext)/dr)*param[i][j];
		}
	}
	
	return vm.ieval(Gradient(p_i) * C * Gradient(p_j, true), e,v) ;
}

bool RadialStiffnessGradient::fractured() const
{
	return false ;
}

Matrix RadialStiffnessGradient::getTensor(const Point & p) const
{
	VirtualMachine vm ;
	
	FunctionMatrix C(3,3) ;
	
	for(size_t i = 0 ; i < 3 ; i++)
	{
		for(size_t j = 0 ; j < 3 ; j++)
		{
// 				double val =  (paramAlt[i][j] - param[i][j])/dr ;
			C[i][j] = ((r-r_int)/dr)*paramAlt[i][j] - ((r-r_ext)/dr)*param[i][j];
		}
	}
	
	return vm.eval(C, p.x, p.y) ;
}

Matrix RadialStiffnessGradient::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const
{
	
	FunctionMatrix C(3,3) ;
	std::vector<Variable> v ;
	v.push_back(XI);
	v.push_back(ETA);

	for(size_t i = 0 ; i < 3 ; i++)
	{
		for(size_t j = 0 ; j < 3 ; j++)
		{
			C[i][j] = ((r-r_int)/dr)*paramAlt[i][j] - ((r-r_ext)/dr)*param[i][j];
		}
	}
	
	return VirtualMachine().ieval(Gradient(p_i) * C * Gradient(p_j, true), gp, Jinv,v) ;
}

Form * RadialStiffnessGradient::getCopy() const 
{
	return new RadialStiffnessGradient(*this) ;
}

Vector RadialStiffnessGradient::getForces(const ElementState & s, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const 
{
	return Vector(0) ;
}
