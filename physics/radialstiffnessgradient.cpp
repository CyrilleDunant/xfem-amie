//
// C++ Implementation: radialstiffnessgradient
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "radialstiffnessgradient.h"

using namespace Mu ;

RadialStiffnessGradient::RadialStiffnessGradient(double E_int, double nu_int, double rint, double E_ext, double nu_ext, double rext, Point c) : LinearForm(Matrix(3,3), false, true, 2), paramAlt(3,3), r_ext(rext), r_int(rint), dr(r_ext-r_int), centre(c), dfunc()
{

	param[0][0] = E_int/(1-nu_int*nu_int) ; param[0][1] =E_int/(1-nu_int*nu_int)*nu_int ; param[0][2] = 0 ;
	param[1][0] = E_int/(1-nu_int*nu_int)*nu_int ; param[1][1] = E_int/(1-nu_int*nu_int) ; param[1][2] = 0 ; 
	param[2][0] = 0 ; param[2][1] = 0 ; param[2][2] = E_int/(1-nu_int*nu_int)*(1.-nu_int)/2. ; 
	
	paramAlt[0][0] = E_ext/(1-nu_ext*nu_ext) ;        paramAlt[0][1] = E_ext/(1-nu_ext*nu_ext)*nu_ext ; paramAlt[0][2] = 0 ;
	paramAlt[1][0] = E_ext/(1-nu_ext*nu_ext)*nu_ext ; paramAlt[1][1] = E_ext/(1-nu_ext*nu_ext) ;        paramAlt[1][2] = 0 ; 
	paramAlt[2][0] = 0 ; paramAlt[2][1] = 0 ;         paramAlt[2][2] = E_ext/(1-nu_ext*nu_ext)*(1.-nu_ext)/2. ; 
	this->space_d = true ;

	criterion = nullptr ;
	frac = false ;
	change  = false ;
	previousDamage = 0 ;
	damage = 0 ;
}

RadialStiffnessGradient::~RadialStiffnessGradient() 
{
	delete criterion ;
}

void RadialStiffnessGradient::transform(const Function & x, const Function & y)
{
	r = f_sqrt(((x-centre.x)^2)+((y-centre.y)^2)) -r_int;
	r.compile() ;
}

void RadialStiffnessGradient::setFractureCriterion(FractureCriterion * crit)
{
	criterion = crit ;
}

void RadialStiffnessGradient::step(double timestep, ElementState & currentState, double maxscore) 
{
	if(criterion == nullptr)
		return ;

	change = false ;

	if(!frac  )
	{
		dfunc.step(currentState, maxscore) ;
		previousDamage = damage ;
		
		Vector state = dfunc.getState() ;
		damage = 0 ;
		for(size_t i = 0 ; i < state.size() ; i++)
			damage += state[i] ;
		change = true ;//std::abs(damage-previousDamage) > 1e-12 ;
		if(damage > .9)
		{
			frac = true ;
// 			damage = .9999 ;
// 			param[0][1] = 0 ;param[0][1] = 0 ;
// 			param[2][2] *= 0.0001 ;
// 			this->type = VOID_BEHAVIOUR ;
		}
	}

}

bool RadialStiffnessGradient::changed() const
{
	return change ;
} 

bool RadialStiffnessGradient::fractured() const
{
	return frac;
}

Matrix RadialStiffnessGradient::getTensor(const Point & p, IntegrableEntity * e, int g) const
{
	VirtualMachine vm ;
	
	FunctionMatrix C(3,3) ;
	Matrix dE ((paramAlt-param)/dr) ;
	Matrix d(3,3) ;
	for(size_t i = 0 ; i < 3 ; i++)
	{
		for(size_t j = 0 ; j < 3 ; j++)
		{
// 				double val =  (paramAlt[i][j] - param[i][j])/dr ;
// 			C[i][j] = ((r-r_int)/dr)*paramAlt[i][j] - ((r-r_ext)/dr)*param[i][j];
			C[i][j] =  r*dE[i][j] + param[i][j];
		}
	}
	
	return dfunc.apply(vm.eval(C, p.x, p.y)) ;
}

void RadialStiffnessGradient::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{
	
	FunctionMatrix C(3,3) ;
	std::vector<Variable> v ;
	v.push_back(XI);
	v.push_back(ETA);

	Matrix dE ((paramAlt-param)/dr) ;
	Matrix d(3,3) ;
	d.array() = 1 ;
	d = dfunc.apply(d) ;
	for(size_t i = 0 ; i < 3 ; i++)
	{
		for(size_t j = 0 ; j < 3 ; j++)
		{
			C[i][j] =  r*dE[i][j]*d[i][j] + param[i][j]*d[i][j];
		}
	}
	
	vm->ieval(Gradient(p_i) * C * Gradient(p_j, true), gp, Jinv,v,ret) ;
}

Form * RadialStiffnessGradient::getCopy() const 
{
	return new RadialStiffnessGradient(*this) ;
}
