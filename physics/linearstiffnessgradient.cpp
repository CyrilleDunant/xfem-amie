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

#include "linearstiffnessgradient.h"
#include "../polynomial/vm_function_extra.h"

using namespace Amie ;

LinearStiffnessGradient::LinearStiffnessGradient(double E_int, double nu_int, double E_ext, double nu_ext, Point l, Point r) : LinearForm(Matrix(3,3), false, true, 2), paramAlt(3,3), left(l), right(r)
{

    v.push_back(XI);
    v.push_back(ETA);
    param[0][0] = E_int/(1.-nu_int*nu_int) ;
    param[0][1] =E_int/(1.-nu_int*nu_int)*nu_int ;
    param[0][2] = 0 ;
    param[1][0] = E_int/(1.-nu_int*nu_int)*nu_int ;
    param[1][1] = E_int/(1.-nu_int*nu_int) ;
    param[1][2] = 0 ;
    param[2][0] = 0 ;
    param[2][1] = 0 ;
    param[2][2] = E_int/(1-nu_int*nu_int)*(1.-nu_int)/2. ;

    paramAlt[0][0] = E_ext/(1.-nu_ext*nu_ext) ;
    paramAlt[0][1] = E_ext/(1.-nu_ext*nu_ext)*nu_ext ;
    paramAlt[0][2] = 0 ;
    paramAlt[1][0] = E_ext/(1.-nu_ext*nu_ext)*nu_ext ;
    paramAlt[1][1] = E_ext/(1.-nu_ext*nu_ext) ;
    paramAlt[1][2] = 0 ;
    paramAlt[2][0] = 0 ;
    paramAlt[2][1] = 0 ;
    paramAlt[2][2] = E_ext/(1.-nu_ext*nu_ext)*(1.-nu_ext)/2. ;
    this->space_d = true ;
}

LinearStiffnessGradient::~LinearStiffnessGradient() { }

void LinearStiffnessGradient::transform(const ElementarySurface * e)
{
    Function l(left, e) ;
    Function r(right, e) ;
    double t = dist(left, right) ;
    t *=t ;
    s = .5 + f_sqrt(2.-(2./t)*l*l) ;
// 	s.setTransform(e);
}

bool LinearStiffnessGradient::fractured() const
{
    return false ;
}

Matrix LinearStiffnessGradient::getTensor(const Point & p, IntegrableEntity * e, int g) const
{
    VirtualMachine vm ;

    FunctionMatrix C(3,3) ;

    for(size_t i = 0 ; i < 3 ; i++)
    {
        for(size_t j = 0 ; j < 3 ; j++)
        {
            C[i][j] = s*paramAlt[i][j] - (s-1.)*param[i][j];
        }
    }
    return vm.eval(C, p.getX(), p.getY()) ;
}

void LinearStiffnessGradient::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{

    FunctionMatrix C(3,3) ;

    for(size_t i = 0 ; i < 3 ; i++)
    {
        for(size_t j = 0 ; j < 3 ; j++)
        {
            C[i][j] = s*paramAlt[i][j] - (s-1.)*param[i][j];
// 			C[i][j].compile() ;
        }
    }

    vm->ieval(Gradient(p_i) * C * Gradient(p_j, true), gp, Jinv,v,ret) ;
}

Form * LinearStiffnessGradient::getCopy() const
{

    LinearStiffnessGradient * copy =  new LinearStiffnessGradient(*this) ;

    return copy ;
}

