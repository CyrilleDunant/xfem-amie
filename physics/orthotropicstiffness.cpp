//
// C++ Implementation: stiffness
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "orthotropicstiffness.h"
#include "../mesher/delaunay.h"
#include "fracturecriteria/vonmises.h"
#include "fracturecriteria/mohrcoulomb.h"
#include "fracturecriteria/nonlocalvonmises.h"
#include <valarray>


using namespace Mu ;

OrthothropicStiffness::OrthothropicStiffness(double E_1, double E_2, double G,  double nu, double angle) : LinearForm(Material::orthothropicCauchyGreen(E_1, E_2, G,  nu), true, false, 2) , E_1(E_1),
		E_2(E_2),
		E_3(0),
		G_1(G),
		G_2(0),
		G_3(0), 
		nu(nu),
		angle(angle)
{
	v.push_back(XI);
	v.push_back(ETA);
	change = false ;
} ;

OrthothropicStiffness::OrthothropicStiffness(double E_1, double E_2, double E_3, double G_1, double G_2, double G_3,  double nu, double angle) : LinearForm(Material::orthothropicCauchyGreen(E_1, E_2, E_3, G_1, G_2, G_3,  nu), true, false, 3) , E_1(E_1),
		E_2(E_2),
		E_3(E_3),
		G_1(G_1),
		G_2(G_2),
		G_3(G_3), 
		nu(nu),
		angle(angle)
{
	v.push_back(XI);
	v.push_back(ETA);
	v.push_back(ZETA);
	change = false ;
} ;

OrthothropicStiffness::~OrthothropicStiffness() { } ;

void OrthothropicStiffness::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{
	Matrix transform(3,3) ;
	transform[0][0] = cos(angle)*cos(angle) ;      transform[0][1] = sin(angle)*sin(angle) ;     transform[0][2] = 2.*sin(angle)*cos(angle) ;
	transform[1][0] = sin(angle)*sin(angle) ;   transform[1][1] = cos(angle)*cos(angle) ;     transform[1][2] = -2.*sin(angle)*cos(angle) ;
	transform[2][0] = -sin(angle)*cos(angle) ;  transform[2][1] = sin(angle)*cos(angle) ; transform[2][2] = cos(angle)*cos(angle)-sin(angle)*sin(angle) ; 
	Matrix effective = transform*(param*transform.transpose()) ;
	vm->ieval(Gradient(p_i) * effective * Gradient(p_j, true), gp, Jinv,v, ret) ;

}

void OrthothropicStiffness::step(double timestep, ElementState & currentState)
{
	change = false ;
	if(timestep > POINT_TOLERANCE_2D)
	{
		angle +=.01 ;
		change = true ;
		currentState.getParent()->behaviourUpdated = true ;
	}
}

bool OrthothropicStiffness::fractured() const
{
	return false ;
}

Form * OrthothropicStiffness::getCopy() const 
{
	if(v.size() == 2)
		return new OrthothropicStiffness(E_1, E_2, G_1,  nu, angle) ;
	
	return new OrthothropicStiffness(E_1, E_2, E_3, G_1, G_2, G_3,  nu, angle) ;
}




