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

OrthothropicStiffness::OrthothropicStiffness(double E_1, double E_2, double nu_12,  double nu_21, double angle, bool poissondefined) : LinearForm(Material::orthothropicCauchyGreen(E_1, E_2, E_1*E_2/(E_1*(1.+nu_12)+E_2*(1.+nu_21)),  (nu_12+nu_21)*.5), true, false, 2) , E_1(E_1),
E_2(E_2),
E_3(0),
G_1(E_1*E_2/(E_1*(1.+nu_12)+E_2*(1.+nu_21))),
G_2(0),
G_3(0), 
nu((nu_12+nu_21)*.5),
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
	if(v.size() == 2)
	{
		Matrix effective(3,3); 
		getTensor(effective) ;
		vm->ieval(Gradient(p_i) * effective * Gradient(p_j, true), gp, Jinv,v, ret) ;
	}
	else
	{
		Matrix effective(6,6); 
		getTensor(effective) ;
		vm->ieval(Gradient(p_i) * effective * Gradient(p_j, true), gp, Jinv,v, ret) ;
	}
	

}

void OrthothropicStiffness::getTensor(Matrix & m) const
{
	if(v.size() == 2)
	{
		Matrix transform(3,3) ;
		transform[0][0] =  cos(angle)*cos(angle) ;  transform[0][1] = sin(angle)*sin(angle) ; transform[0][2] =  2.*sin(angle)*cos(angle) ;
		transform[1][0] =  sin(angle)*sin(angle) ;  transform[1][1] = cos(angle)*cos(angle) ; transform[1][2] = -2.*sin(angle)*cos(angle) ;
		transform[2][0] = -sin(angle)*cos(angle) ;  transform[2][1] = sin(angle)*cos(angle) ; transform[2][2] =     cos(angle)*cos(angle) - sin(angle)*sin(angle) ; 
		m = transform*(param*transform.transpose()) ;
	}
	else
	{
		Matrix transform(6,6) ;
		//this rotates only along the x axis
		double lx = cos(angle) ;
		double ly = sin(angle) ;
		double lz = 0 ;
		double rx = -sin(angle) ;
		double ry =  cos(angle) ;
		double rz = 0 ;
		double tx = 0 ;
		double ty = 0 ;
		double tz = 0 ;
		transform[0][0] = lx*lx ; transform[0][1] = ly*ly ; transform[0][2] = lz*lz ; transform[0][3] = lx*ly ; transform[0][4] = lx*lz ; transform[0][5] = 2.*ly*lz ;
		transform[1][0] = rx*rx ; transform[1][1] = ry*ry ; transform[1][2] = rz*rz ; transform[1][3] = rx*ry ; transform[1][4] = rx*rz ; transform[1][5] = 2.*ry*rz ; 
		transform[2][0] = tx*tx ; transform[2][1] = ty*ty ; transform[2][2] = tz*tz ; transform[2][3] = tx*ty ; transform[2][4] = tx*tz ; transform[2][5] = 2.*ty*tz ;
		transform[3][0] = lx*rx ; transform[3][1] = ly*ry ; transform[3][5] = lz*rz ; transform[3][3] = lx*ry+ly*rx ; transform[3][4] = lz*rx+lx*rz ; transform[3][5] = ly*rz+lz*ry ;
		transform[4][0] = lx*tx ; transform[4][1] = ly*ty ; transform[4][2] = lz*tz ; transform[4][3] = tx*ly+ly*lx ; transform[4][4] = tz*lx+tx*lz ; transform[4][5] = ty*lz+tz*ly ;
		transform[5][0] = rx*tx ; transform[5][1] = ry*ty ; transform[5][2] = rz*tz ; transform[5][3] = rx*ty+ry*tx ; transform[5][4] = rz*tx+rx*tz ; transform[5][5] = ry*tz+rz*ty ;
		m = transform*(param*transform.transpose()) ;
	}
}

Matrix OrthothropicStiffness::getTensor(const Point & p)
{
	if(v.size() == 2)
	{
		Matrix transform(3,3) ;
		transform[0][0] =  cos(angle)*cos(angle) ;  transform[0][1] = sin(angle)*sin(angle) ; transform[0][2] =  2.*sin(angle)*cos(angle) ;
		transform[1][0] =  sin(angle)*sin(angle) ;  transform[1][1] = cos(angle)*cos(angle) ; transform[1][2] = -2.*sin(angle)*cos(angle) ;
		transform[2][0] = -sin(angle)*cos(angle) ;  transform[2][1] = sin(angle)*cos(angle) ; transform[2][2] =     cos(angle)*cos(angle) - sin(angle)*sin(angle) ; 
		return transform*(param*transform.transpose()) ;
	}
	else
	{
		Matrix transform(6,6) ;
		//this rotates only along the x axis
		double lx = cos(angle) ;
		double ly = sin(angle) ;
		double lz = 0 ;
		double rx = -sin(angle) ;
		double ry = cos(angle) ;
		double rz = 0 ;
		double tx = 0 ;
		double ty = 0 ;
		double tz = 0 ;
		transform[0][0] = lx*lx ; transform[0][1] = ly*ly ; transform[0][2] = lz*lz ; transform[0][3] = lx*ly ; transform[0][4] = lx*lz ; transform[0][5] = 2.*ly*lz ;
		transform[1][0] = rx*rx ; transform[1][1] = ry*ry ; transform[1][2] = rz*rz ; transform[1][3] = rx*ry ; transform[1][4] = rx*rz ; transform[1][5] = 2.*ry*rz ; 
		transform[2][0] = tx*tx ; transform[2][1] = ty*ty ; transform[2][2] = tz*tz ; transform[2][3] = tx*ty ; transform[2][4] = tx*tz ; transform[2][5] = 2.*ty*tz ;
		transform[3][0] = lx*rx ; transform[3][1] = ly*ry ; transform[3][5] = lz*rz ; transform[3][3] = lx*ry+ly*rx ; transform[3][4] = lz*rx+lx*rz ; transform[3][5] = ly*rz+lz*ry ;
		transform[4][0] = lx*tx ; transform[4][1] = ly*ty ; transform[4][2] = lz*tz ; transform[4][3] = tx*ly+ly*lx ; transform[4][4] = tz*lx+tx*lz ; transform[4][5] = ty*lz+tz*ly ;
		transform[5][0] = rx*tx ; transform[5][1] = ry*ty ; transform[5][2] = rz*tz ; transform[5][3] = rx*ty+ry*tx ; transform[5][4] = rz*tx+rx*tz ; transform[5][5] = ry*tz+rz*ty ;
		return transform*(param*transform.transpose()) ;
	}
}

void OrthothropicStiffness::step(double timestep, ElementState & currentState)
{
	change = false ;
	if(timestep > POINT_TOLERANCE_2D)
	{
		angle +=.1 ;
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




