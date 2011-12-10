//
// C++ Implementation: kelvinvoight
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2008-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
//
#include "kelvinvoight.h"
#include "homogenization/composite.h"
#include "../utilities/matrixops.h"
#include "../features/boundarycondition.h"

using namespace Mu ;

KelvinVoight::KelvinVoight(const Matrix & rig, const Matrix & e) : LinearForm(rig, false, false, rig.numRows()/3+1), eta(e)
{
	rig.print() ;
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() > 9)
		v.push_back(ZETA);
	v.push_back(TIME_VARIABLE);
} ;

KelvinVoight::~KelvinVoight() { } ;

void KelvinVoight::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{

	Matrix temp(ret) ;
	Matrix temp0(ret) ;
	Matrix temp1(ret) ;
	
	vm->ieval(GradientDot(p_i) * param * Gradient(p_j, true), gp, Jinv,v,temp0) ;
	vm->ieval(Gradient(p_i) * param * GradientDot(p_j, true), gp, Jinv,v,temp1);

	vm->ieval(Gradient(p_j) * param * Gradient(p_i, true), gp, Jinv,v,ret);
	
	if(std::abs((temp0+temp1).array()).max() < 1e-6*std::abs(ret.array()).max())
		ret *= 0 ;
	else 
		ret /= vm->ieval(p_j, gp) ;//[81000 825000]
	

	vm->ieval(GradientDot(p_i) * eta * GradientDot(p_j, true), gp, Jinv,v,temp);

// 	temp.print() ;
	
	ret = ret + temp ;
}

bool KelvinVoight::fractured() const
{
	return false ;
}

bool KelvinVoight::changed() const
{
	return false ;
} 

Form * KelvinVoight::getCopy() const 
{
	return new KelvinVoight(*this) ;
}

IncrementalKelvinVoight::IncrementalKelvinVoight(const Matrix &rig, const Matrix &e, double dt) : LinearForm(rig, false, false, rig.numRows()/3+1), eta(e), stiff(rig), tau(dt)
{
    v.push_back(XI);
    v.push_back(ETA);
    if(param.size() > 9)
	    v.push_back(ZETA);

    Matrix I = Composite::I4(rig) ;
    Matrix K = eta ;
    if(K.size()==36)
	invert6x6Matrix(K) ;
    else
	invert3x3Matrix(K) ;
    K *= stiff ;
    K *= -dt ; // K = - dt*C/E

    Matrix L = exp(K) ;
    if(L.size()==36)
	invert6x6Matrix(L) ;
    else
	invert3x3Matrix(L) ;

    Matrix J = I - L ; // J = (1 - exp( -dt*C/E ) )
    N.resize(J.numRows(), J.numCols()) ;
    N = I - L ;// in fact, N = J

    if(K.size()==36)
	invert6x6Matrix(K) ;
    else
	invert3x3Matrix(K) ;
    // K = -E/(dt*C)
    J *= K ; // J = -E/(dt*C) * (1 - exp( -dt*C/E ) )
    J += I ; // J = (1 -E/(dt*C) * (1 - exp( -dt*C/E ) ) )

    if(J.size()==36)
	invert6x6Matrix(J) ;
    else
	invert3x3Matrix(J) ;

    param *= J ;

    phi.resize(param.numCols()) ; // initialize phi = 0 ;
    up.resize(param.numCols());
}

IncrementalKelvinVoight::~IncrementalKelvinVoight()
{

}

void IncrementalKelvinVoight::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{
	vm->ieval(Gradient(p_i) * param * Gradient(p_j, true), gp, Jinv,v,ret) ;
}


Form * IncrementalKelvinVoight::getCopy() const
{
    return new IncrementalKelvinVoight(*this) ;
}

Vector IncrementalKelvinVoight::getImposedStress(const Point &p) const
{
    Matrix S = param * N ;
    return S * phi ;
}

std::vector<BoundaryCondition * > IncrementalKelvinVoight::getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const
{
    Vector phi_n(v.size()) ;
    int n = (id - id%v.size())/v.size() ; // n is point id
    for(int i = 0 ; i < param.numCols() ; i++)
	phi_n[i] = phi[n*param.numCols()+i] ;

    Matrix S = param * N ;

    Vector f = VirtualMachine().ieval(Gradient(p_i) * (S * phi_n), gp, Jinv,v) ;

	std::vector<BoundaryCondition * > ret ;
	if(f.size() == 2)
	{
		ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_XI, dynamic_cast<ElementarySurface *>(s.getParent()), id, f[0]));
		ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_ETA, dynamic_cast<ElementarySurface *>(s.getParent()), id, f[1]));
	}
	if(f.size() == 3)
	{
		ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_XI, dynamic_cast<ElementaryVolume *>(s.getParent()), id, f[0]));
		ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_ETA, dynamic_cast<ElementaryVolume *>(s.getParent()), id, f[1]));
		ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_ZETA, dynamic_cast<ElementaryVolume *>(s.getParent()), id, f[2]));
	}
	return ret ;
}

void IncrementalKelvinVoight::step(double timestep, ElementState &s)
{
    PointArray bc = s.getParent()->getBoundingPoints() ;
    Vector sigma = s.getStress(bc) ;
    Vector epsilon = s.getStrain(bc) ;

    Matrix compliance = stiff ;
    if(compliance.size()==36)
	invert6x6Matrix(compliance) ;
    else
	invert3x3Matrix(compliance) ;

    for(size_t i = 0 ; i < bc.size() ; i++)
    {
	size_t id = bc[i]->id ;
	if(!up[id])
	{
	    Vector stress(param.numCols()) ;
	    Vector strain(param.numCols()) ;
	    Vector phi_n(param.numCols()) ;

	    for(size_t j = 0 ; j < param.numCols() ; j++)
	    {
		stress[j] = sigma[i*param.numCols() + j] ;
		strain[j] = epsilon[i*param.numCols() + j] ;
		phi_n[j] = phi[id*param.numCols() + j] ;
	    }

	    stress = compliance * stress ;
	    phi_n = phi_n - stress + strain ;

	    up[i] = true ;
	}
    }

    bool done = true ;
    int n = 0 ;
    while(done && n < up.size())
    {
	done = up[n] ;
	n++ ;
    }
    if(done && n==up.size())
	up.resize(up.size());

}

void IncrementalKelvinVoight::resize(size_t num_points)
{
    std::cerr << "----" << std::endl ;
    std::cerr << num_points << std::endl ;
    phi.resize(num_points*param.numCols()) ;
    up.resize(num_points) ;
}
