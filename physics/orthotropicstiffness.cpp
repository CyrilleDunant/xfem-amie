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
#include "../utilities/tensor.h"
#include "fracturecriteria/vonmises.h"
#include "fracturecriteria/mohrcoulomb.h"
#include "fracturecriteria/nonlocalvonmises.h"
#include <valarray>


using namespace Amie ;

void OrthotropicStiffness::setAngle(double angle)
{
//     if(std::abs(angle) < POINT_TOLERANCE)
//     {
//         if(v.size() == 2)
//             transform = identity(3) ;
//         else
//             transform = identity(6) ;
//         transformt = transform.transpose() ;
//         transformset = true ;
//         return ;
//     }
    
    if(v.size() == 2)
    {
//         double c = cos(angle) ;
//         double s = sin(angle) ;
//         transform[0][0] =  c*c ;
//         transform[0][1] = s*s ;
//         transform[0][2] =  2.*s*c ;
//         transform[1][0] =  s*s ;
//         transform[1][1] = c*c ;
//         transform[1][2] = -2.*s*c ;
//         transform[2][0] = -s*c ;
//         transform[2][1] = s*c ;
//         transform[2][2] = c*c - s*s ;
//         transformt = transform.transpose() ;
        transformset = true ;
    }
    else
    {
//         double lx = cos(angle) ;
//         double ly = 0 ;
//         double lz = 0 ;
//         double rx = cos(angle) ;
//         double ry = 0 ;
//         double rz = 0 ;
//         double tx = cos(angle) ;
//         double ty = 0 ;
//         double tz = 0 ;
//         transform[0][0] = lx*lx ;
//         transform[0][1] = ly*ly ;
//         transform[0][2] = lz*lz ;
//         transform[0][3] = lx*ly ;
//         transform[0][4] = lx*lz ;
//         transform[0][5] = ly*lz ;
//         transform[1][0] = rx*rx ;
//         transform[1][1] = ry*ry ;
//         transform[1][2] = rz*rz ;
//         transform[1][3] = rx*ry ;
//         transform[1][4] = rx*rz ;
//         transform[1][5] = ry*rz ;
//         transform[2][0] = tx*tx ;
//         transform[2][1] = ty*ty ;
//         transform[2][2] = tz*tz ;
//         transform[2][3] = tx*ty ;
//         transform[2][4] = tx*tz ;
//         transform[2][5] = ty*tz ;
//         transform[3][0] = 2.*lx*rx ;
//         transform[3][1] = 2.*ly*ry ;
//         transform[3][5] = 2.*lz*rz ;
//         transform[3][3] = lx*ry+ly*rx ;
//         transform[3][4] = lz*rx+lx*rz ;
//         transform[3][5] = ly*rz+lz*ry ;
//         transform[4][0] = 2.*lx*tx ;
//         transform[4][1] = 2.*ly*ty ;
//         transform[4][2] = 2.*lz*tz ;
//         transform[4][3] = tx*ly+ly*lx ;
//         transform[4][4] = tz*lx+tx*lz ;
//         transform[4][5] = ty*lz+tz*ly ;
//         transform[5][0] = 2.*rx*tx ;
//         transform[5][1] = 2.*ry*ty ;
//         transform[5][2] = 2.*rz*tz ;
//         transform[5][3] = rx*ty+ry*tx ;
//         transform[5][4] = rz*tx+rx*tz ;
//         transform[5][5] = ry*tz+rz*ty ;
//         transformt = transform.transpose() ;
        transformset = true ;
    }
    
    alpha = angle ;
    param = paramBase ;
    Tensor::rotate4thOrderTensor2D( param, alpha) ;
}

void OrthotropicStiffness::setStiffness(double E_1, double E_2, double G, double nu)
{
    this->E_1 = E_1 ;
    this->E_2 = E_2 ;
    paramBase = Tensor::orthotropicCauchyGreen(E_1, E_2, G,  nu, PLANE_STRAIN) ;
    param = paramBase ;
//     Tensor::rotate4thOrderTensor2D( param, alpha) ;
    if(transformset)
        Tensor::rotate4thOrderTensor2D( param, alpha, POINT_TOLERANCE) ; ;
}

void OrthotropicStiffness::setStiffness(double E_1, double E_2, double E_3, double G_1, double G_2, double G_3, double nu)
{
    this->E_1 = E_1 ;
    this->E_2 = E_2 ;
    this->E_3 = E_3 ;
    paramBase = Tensor::orthotropicCauchyGreen(E_1, E_2,E_3, G_1, G_2, G_3,  nu) ;
    if(transformset)
    {
        param = paramBase ;
	Tensor::rotate4thOrderTensor2D( param, alpha, POINT_TOLERANCE) ; ;
    }
}

OrthotropicStiffness::OrthotropicStiffness(double E_1, double E_2, double G,  double nu, double angle) : LinearForm(Tensor::orthotropicCauchyGreen(E_1, E_2, G,  nu, PLANE_STRAIN), true, false, 2) , E_1(E_1), E_2(E_2) ,
    transform(3,3),
    transformt(3,3)
{
    transformset = false ;
    setStiffness(E_1, E_2, G, nu) ;
    v.push_back(XI);
    v.push_back(ETA);
    paramBase = param ;
    setAngle(angle);
} 

OrthotropicStiffness::OrthotropicStiffness(const OrthotropicStiffness * source) : LinearForm(source->param, true, false, source->getNumberOfDegreesOfFreedom()) , transform(source->transform),
    transformt(source->transformt), paramBase(source->paramBase),transformset(source->transformset),v(source->v)
{
} 

// OrthotropicStiffness::OrthotropicStiffness(double E_1, double E_2, double nu_12,  double nu_21, double angle, bool poissondefined) : LinearForm(Tensor::orthotropicCauchyGreen(E_1, E_2, E_1*E_2/(E_1*(1.+nu_12)+E_2*(1.+nu_21)),  (nu_12+nu_21)*.5, PLANE_STRAIN), true, false, 2) , E_1(E_1), E_2(E_2),
//     transform(3,3),
//     transformt(3,3)
// {
//     transformset = false ;
//     v.push_back(XI);
//     v.push_back(ETA);
//     paramBase = param ;
//     setAngle(angle);
//     
// }

OrthotropicStiffness::OrthotropicStiffness(double E_1, double E_2, double E_3, double G_1, double G_2, double G_3,  double nu, double angle) : LinearForm(Tensor::orthotropicCauchyGreen(E_1, E_2, E_3, G_1, G_2, G_3,  nu), true, false, 3),E_1(E_1), E_2(E_2),E_3(E_3),
    transform(6,6),
    transformt(6,6)
{
    transformset = false ;
    setStiffness(E_1, E_2, E_3, G_1, G_2, G_3,  nu) ;
    v.push_back(XI);
    v.push_back(ETA);
    v.push_back(ZETA);
    paramBase = param ;
    setAngle(angle);
}

OrthotropicStiffness::~OrthotropicStiffness() { } 

void OrthotropicStiffness::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{
    vm->ieval(Gradient(p_i) * param * Gradient(p_j, true), gp, Jinv,v, ret) ;

}

Matrix OrthotropicStiffness::getTensor(const Point& p, IntegrableEntity* e , int g ) const
{
    return param ;
}

bool OrthotropicStiffness::fractured() const
{
    return false ;
}

Form * OrthotropicStiffness::getCopy() const
{
    if(v.size() == 2)
    {
        OrthotropicStiffness * copy = new OrthotropicStiffness(this) ;

        return copy ;
    }

    OrthotropicStiffness * copy =  new OrthotropicStiffness(this) ;

    return copy ;
}




