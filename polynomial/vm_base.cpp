// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
//
// Copyright: See COPYING file that comes with this distribution



#include "vm_base.h"
#include "../delaunay.h"

using namespace Mu ;

VirtualMachine::VirtualMachine() :stack() { } ;

double VirtualMachine::eval(const Function &f, const double x, const double y, const double z, const double t, const double u, const double v, const double w) 
{
	stack.reset() ;
	Context context(stack, x, y, z, t, u, v, w) ;
	
	size_t size = f.size() ;
	for(size_t i = 0 ; i < size  ; i++)
		f.tokenEval(i,context) ;

// 	if(std::abs(stack[0]) < 1e-6)
// 		return 0 ;
		
	return  stack[0] ;
}

double VirtualMachine::eval(const std::vector<RefCountedToken> &f, const double x, const double y, const double z, const double t, const double u, const double v, const double w) 
 {
 	stack.reset() ;
	Context context(stack, x, y, z, t, u, v, w) ;
	size_t size = f.size() ;
	for(size_t i = 0 ; i < size  ; i++)
		f[i]->eval(context) ;

// 	if(std::abs(stack[0]) < 1e-6)
// 		return 0 ;
		
	return  stack[0] ;
 }

double VirtualMachine::eval(const Function &f, const Point & p, const double t, const Point & p_) 
{
	return eval(f, p.x, p.y, p.z,t,p_.x, p_.y, p_.z) ;
}


double VirtualMachine::eval(const Function &f, const Point *p, const double t, const Point * p_) 
{
	if(p_)
		return eval(f, p->x, p->y, p->z, t, p_->x, p_->y, p_->z) ;
	return eval(f, p->x, p->y, p->z, t) ;
}


void VirtualMachine::print(const Function &f) const
{
	std::cout << f.size() << " tokens" << std::endl ;
	for(size_t i = 0 ; i < f.size()  ; i++)
	{
		std::cout << f.getToken(i)->print() << " " << std::flush ;
	}
	
	std::cout << std::endl ;
}

void VirtualMachine::print(const std::vector<RefCountedToken> &f) const
{
	std::cout << f.size() << " tokens" << std::endl ;
	for(size_t i = 0 ; i < f.size()  ; i++)
	{
		std::cout << f[i]->print() << " " << std::flush ;
	}
	
	std::cout << std::endl ;
}

Matrix VirtualMachine::eval(const FunctionMatrix &f, const double x, const double y, const double z, const double t, const double u, const double v, const double w) 
{
	Matrix ret(f.numRows(), f.numCols()) ;
	
	for(size_t i = 0 ;  i < f.numRows() ; i++)
	{
		for(size_t j = 0 ;  j < f.numCols() ; j++)
		{
			ret[i][j] = eval(f[i][j], x, y, z, t, u, v, w) ;
		}
	}
	
	return ret ;
	
}

Matrix VirtualMachine::eval(const FunctionMatrix &f, const Point & p, const double t, const Point & p_) 
{
	return eval(f, p.x, p.y, p.z, t, p_.x, p_.y, p_.z) ;
}

Matrix VirtualMachine::eval(const FunctionMatrix &f, const Point * p, const double t, const Point *p_) 
{
	if(p_)
		return eval(f, p->x, p->y, p->z, t, p_->x, p_->y, p_->z) ;
	return eval(f, p->x, p->y, p->z, t) ;
}

double VirtualMachine::deval(const Function &f, const Point&p,  const double x, const double y, const double z , const double eps , bool normed)
{
	double d = deval(f,XI,x,y,z,0,0,0,0,eps)*p.x + deval(f,ETA,x,y,z,0,0,0,0,eps)*p.y + deval(f,ZETA,x,y,z,0,0,0,0,eps)*p.z ;
	if(normed)
		return d ;
	else
		return d/p.norm() ;
}

Matrix VirtualMachine::deval(const FunctionMatrix &f,const Point&p,  const double x, const double y, const double z , const double eps, bool normed )
{
	Matrix ret(f.numRows(), f.numCols()) ;
	
	for(size_t i = 0 ; i < f.numRows() ; i++)
	{
		for(size_t j = 0 ; j < f.numCols() ; j++)
		{
			ret[i][j] = deval(f[i][j], p, x, y, z, eps, normed) ;
		}
	}
	
	return ret ;
}

double VirtualMachine::deval(const Function &f,  const Point&p, const Point x, const double eps, bool normed )
{
	return deval(f, p, x.x, x.y, x.z, eps, normed) ;
}

Matrix VirtualMachine::deval(const FunctionMatrix &f,  const Point&p, const Point x, const double eps, bool normed )
{
	return deval(f, p, x.x, x.y, x.z, eps, normed) ;
}

double VirtualMachine::deval(const Function &f, const Variable v_,  const double x, const double y , const double z,  const double t, const double u, const double v, const double w, const double eps) 
{
	if(f.isDifferentiable())
	{
		return eval(f.d(v_), x, y, z, t, u, v, w) ;
	}
	else
	{
		switch(v_)
		{
		case ONE : 
			{
				return 0 ;
			}
		case XI : 
			{
// 				std::cout << "called" << std::endl ;
				//f(x) = (f(x-2*eps) - 8*f(x-eps) + 8*f(x+eps) - f(x+eps))/(12*eps)
// 				return ( eval(f, x+eps, y, z, t, u, v, w) - eval(f, x-eps, y, z, t, u, v, w))/(2.*eps) ;
				return (eval(f, x-2.*eps, y, z, t, u, v, w) -8.*eval(f, x-eps, y, z, t, u, v, w) + 8.*eval(f, x+eps, y, z, t, u, v, w) - eval(f, x-2.*eps, y, z, t, u, v, w))/(eps*12.) ;
			}
		case ETA:
			{
// 				return ( eval(f, x, y+eps, z, t, u, v, w) - eval(f, x, y-eps, z, t, u, v, w))/(2.*eps) ;
				return (eval(f, x, y-2.*eps, z, t, u, v, w) -8.*eval(f, x, y-eps, z, t, u, v, w) + 8.*eval(f, x, y+eps, z, t, u, v, w) - eval(f, x, y-2.*eps, z, t, u, v, w))/(eps*12.) ;
			}
		case ZETA:
			{
// 				return ( eval(f, x, y, z+eps, t, u, v, w) - eval(f, x, y, z-eps, t, u, v, w))/(2.*eps) ;
				return (eval(f, x, y, z-2.*eps, t, u, v, w) -8.*eval(f, x, y, z-eps, t, u, v, w) + 8.*eval(f, x, y, z+eps, t, u, v, w) - eval(f, x, y, z-2.*eps, t, u, v, w))/(eps*12.) ;
			}
		case TIME_VARIABLE : 
			{
				return ( eval(f, x, y, z, t+eps, u, v, w) - eval(f, x, y, z, t-eps, u, v, w))/(2.*eps) ;
				//f(x) = (f(x-2*eps) - 8*f(x-eps) + 8*f(x+eps) - f(x+eps))/(12*eps)
// 				return (eval(f, x, y, z, t-2.*eps, u, v, w) -8.*eval(f, x, y, z, t-eps, u, v, w) + 8.*eval(f, x, y, z, t+eps, u, v, w) - eval(f, x, y, z, t-2.*eps, u, v, w))/(eps*12.) ;
			}
		case U_VARIABLE:
			{
				return ( eval(f, x, y, z, t, u+eps, v, w) - eval(f, x, y, z, t, u-eps, v, w))/(2.*eps) ;
// 				return (eval(f, x, y, z, t, u-2.*eps, v, w) -8.*eval(f, x, y, z, t, u-eps, v, w) + 8.*eval(f, x, y, z, t, u+eps, v, w) - eval(f, x, y, z, t, u-2.*eps, v, w))/(eps*12.) ;
			}
		case V_VARIABLE:
			{
				return ( eval(f, x, y, z, t, u, v+eps, w) - eval(f, x, y, z, t, u, v-eps, w))/(2.*eps) ;
// 				return (eval(f, x, y, z,t, u, v-2.*eps, w) -8.*eval(f, x, y, z, t, u, v-eps, w) + 8.*eval(f, x, y, z, t, u, v+eps, w) - eval(f, x, y, z, t, u, v-2.*eps, w))/(eps*12.) ;
			}
		case W_VARIABLE:
			{
				return ( eval(f, x, y, z, t, u, v, w+eps) - eval(f, x, y, z, t, u, v, w-eps))/(2.*eps) ;
// 				return (eval(f, x, y, z, t, u, v, w-2.*eps) -8.*eval(f, x, y, z, t, u, v, w-eps) + 8.*eval(f, x, y, z, t, u, v, w+eps) - eval(f, x, y, z, t, u, v, w-2.*eps))/(eps*12.) ;
			}
		}
	}
	
	return 0 ;
}

Matrix VirtualMachine::deval(const FunctionMatrix &f, const Variable v_,  const double x, const double y, const double z, const double t, const double u, const double v, const double w, const double eps)
{
	Matrix ret(f.numRows(), f.numCols()) ;
	
	for(size_t i = 0 ; i < f.numRows() ; i++)
	{
		for(size_t j = 0 ; j < f.numCols() ; j++)
		{
			ret[i][j] = deval(f[i][j], v_, x, y, z,t,u,v,w,  eps) ;
		}
	}
	
	return ret ;
}

double VirtualMachine::deval(const Function &f, const Variable v, const Point p, const double t, const Point p_, const double eps)
{
	return deval(f, v, p.x, p.y, p.z, t, p_.x, p_.y, p_.z, eps) ;
}

Matrix VirtualMachine::deval(const FunctionMatrix &f, const Variable v, const Point p, const double t, const Point p_, const double eps)
{
	return deval(f, v, p.x, p.y, p.z, t, p_.x, p_.y, p_.z, eps) ;
}

double VirtualMachine::ieval(const Function &f, const std::valarray< std::pair<Point, double> > &gp)
{
	double ret = 0 ;
	for(size_t i = 0 ; i < gp.size() ; i++)
		ret += eval(f, gp[i].first)*gp[i].second ;
	
	return ret ;
}

double VirtualMachine::ieval(const Function &f, const IntegrableEntity *e)
{

	std::valarray< std::pair<Point, double> > gp = e->getGaussPoints() ;
	double ret = 0 ;
	for(size_t i = 0 ; i < gp.size() ; i++)
		ret += eval(f, gp[i].first)*gp[i].second ;
	
	return ret ;
}

double VirtualMachine::ieval(const Function &f, const std::vector<std::pair<Segment *, IntegrableEntity *> > & gamma)
{
	double ret = 0;
	for(size_t i = 0 ; i < gamma.size() ; i++)
	{
		std::vector< std::pair<Point, double> > gaussPoints = gamma[i].first->getGaussPoints() ;
		for(size_t j = 0 ; j <  gaussPoints.size() ; j++)
		{
			ret += eval(f, gamma[i].second->inLocalCoordinates(gaussPoints[j].first))*gaussPoints[j].second ;
		}
	}
	
	return ret ;
}


double VirtualMachine::ieval(Vector &f, const IntegrableEntity *e)
{
	std::valarray< std::pair<Point, double> > gp = e->getGaussPoints() ;
	double ret = 0 ;
	for(size_t i = 0 ; i < gp.size() ; i++)
		ret += f[i]*gp[i].second ;
	
	return ret ;
}

Matrix VirtualMachine::ieval(const FunctionMatrix &f, const IntegrableEntity *e)
{

	Matrix ret(f.numRows(), f.numCols()) ;
		
	for(size_t i = 0 ; i < f.numRows() ; i++)
	{
		for(size_t j = 0 ; j < f.numCols() ; j++)
		{
			ret[i][j] = ieval(f[i][j],e) ;
		}
	}
		
	return ret ;
}

Matrix VirtualMachine::ieval(const FMtMtFM &f, const IntegrableEntity *e)
{
	Matrix a = ieval(f.first, e) ;
	
	Matrix c = ieval(f.third, e) ;
	
	return a*(f.second)*c ;
}


Matrix VirtualMachine::geval(const Function &f, const IntegrableEntity *e, const double x, const double y, const double z, bool transpose)
{
	Matrix Jinv = e->getInverseJacobianMatrix(Point(x,y,z)) ;
	return geval(f, Jinv, x, y, z, transpose) ;
}

Matrix VirtualMachine::geval(const Gradient &f, const IntegrableEntity *e, const double x, const double  y, const double z)
{
	return geval(f.f, e, x, y, z, f.transpose) ;
}

Matrix VirtualMachine::geval(const Function &f, const Matrix & m, const Point& p, bool transpose )
{
	return geval(f, m, p.x, p.y, p.z, transpose) ;
}

Matrix VirtualMachine::geval(const Function &f, const Matrix & m, const double x, const double y, const double  z, bool transpose )
{
	if(m.size() == 9 )
	{
		if(transpose)
		{
			double dxi = deval(f, XI, x,y,z) ;
			double deta = deval(f, ETA, x,y,z) ;
			double dzeta = deval(f, ZETA, x,y,z) ;
			Matrix ret(6,3) ;
			ret[0][0] = dxi*m[0][0] + deta*m[0][1] + dzeta*m[0][2];
			ret[1][1] = dxi*m[1][0] + deta*m[1][1] + dzeta*m[1][2];
			ret[2][2] = dxi*m[2][0] + deta*m[2][1] + dzeta*m[2][2];
			ret[3][1] = ret[2][2] ;
			ret[3][2] = ret[1][1] ;
			ret[4][0] = ret[2][2] ;
			ret[4][2] = ret[0][0] ;
			ret[5][0] = ret[1][1] ;
			ret[5][1] = ret[0][0] ;
			
			return ret ;
		}
		else
		{
			Matrix ret(3,6) ;
			double dxi = deval(f, XI, x,y,z) ;
			double deta = deval(f, ETA, x,y,z) ;
			double dzeta = deval(f, ZETA, x,y,z) ;
			ret[0][0] = dxi*m[0][0] + deta*m[0][1] + dzeta*m[0][2];
			ret[1][1] = dxi*m[1][0] + deta*m[1][1] + dzeta*m[1][2];
			ret[2][2] = dxi*m[2][0] + deta*m[2][1] + dzeta*m[2][2];
			ret[1][3] = ret[2][2] ;
			ret[2][3] = ret[1][1] ;
			ret[0][4] = ret[2][2] ;
			ret[2][4] = ret[0][0] ;
			ret[0][5] = ret[1][1] ;
			ret[1][5] = ret[0][0] ;
			
			return ret ;
		}
	}
	else if (m.size() == 4)
	{
		if(transpose)
		{
			Matrix ret(3,2) ;
			double dxi = deval(f, XI, x,y,z) ;
			double deta = deval(f, ETA, x,y,z) ;
			ret[0][0] = dxi*m[0][0] + deta*m[0][1] ;
			ret[1][1] = dxi*m[1][0] + deta*m[1][1] ;
			ret[2][0] = ret[1][1] ;
			ret[2][1] = ret[0][0] ;
			
			return ret ;
		}
		else
		{
			Matrix ret(2,3) ;
			double dxi = deval(f, XI, x,y,z) ;
			double deta = deval(f, ETA, x,y,z) ;
			ret[0][0] = dxi*m[0][0] + deta*m[0][1] ;
			ret[1][1] = dxi*m[1][0] + deta*m[1][1] ;
			ret[0][2] = ret[1][1] ;
			ret[1][2] = ret[0][0] ;
			
			return ret ;
		}
	}
	return Matrix(0,0) ;
}

Matrix VirtualMachine::gveval(const Function &f, const IntegrableEntity *e, const double x, const double y, const double z, bool transpose)
{
	Matrix Jinv = e->getInverseJacobianMatrix(Point(x,y,z)) ;
	return gveval(f, Jinv, x, y, z, transpose) ;
}

Matrix VirtualMachine::gveval(const VectorGradient &f, const IntegrableEntity *e, const double x, const double  y, const double z)
{
	return gveval(f.f, e, x, y, z, f.transpose) ;
}

Matrix VirtualMachine::gveval(const Function &f, const Matrix & m, const double x, const double y, const double  z, bool transpose )
{
	if(m.size() == 9 )
	{
		if(transpose)
		{
			Matrix ret(3,1) ;
			ret[0][0] = deval(f, XI, x,y,z)*m[0][0] + deval(f, ETA, x,y,z)*m[0][1] + deval(f, ZETA, x,y,z)*m[0][2];
			ret[1][0] = deval(f, XI, x,y,z)*m[1][0] + deval(f, ETA, x,y,z)*m[1][1] + deval(f, ZETA, x,y,z)*m[1][2];
			ret[2][0] = deval(f, XI, x,y,z)*m[2][0] + deval(f, ETA, x,y,z)*m[2][1] + deval(f, ZETA, x,y,z)*m[2][2];

			
			return ret ;
		}
		else
		{
			Matrix ret(1,3) ;
			ret[0][0] = deval(f, XI, x,y,z)*m[0][0] + deval(f, ETA, x,y,z)*m[0][1] + deval(f, ZETA, x,y,z)*m[0][2];
			ret[0][1] = deval(f, XI, x,y,z)*m[1][0] + deval(f, ETA, x,y,z)*m[1][1] + deval(f, ZETA, x,y,z)*m[1][2];
			ret[0][2] = deval(f, XI, x,y,z)*m[2][0] + deval(f, ETA, x,y,z)*m[2][1] + deval(f, ZETA, x,y,z)*m[2][2];
	
			return ret ;
		}
	}
	else
	{
		if(transpose)
		{
			Matrix ret(2,1) ;
			ret[0][0] = deval(f, XI, x,y)*m[0][0] + deval(f, ETA, x,y)*m[0][1] ;
			ret[1][0] = deval(f, XI, x,y)*m[1][0] + deval(f, ETA, x,y)*m[1][1] ;

			
			return ret ;
		}
		else
		{
			Matrix ret(1,2) ;
			ret[0][0] = deval(f, XI, x,y)*m[0][0] + deval(f, ETA, x,y)*m[0][1] ;
			ret[0][1] = deval(f, XI, x,y)*m[1][0] + deval(f, ETA, x,y)*m[1][1] ;
			
			return ret ;
		}
	}
}

Matrix VirtualMachine::gveval(const VectorGradient &f, const Matrix & m, const double x, const double y , const double z)
{
	return gveval(f.f, m, x, y, z, f.transpose) ;
}


Matrix VirtualMachine::ieval(const GtM &f, const IntegrableEntity *e)
{

	std::valarray< std::pair<Point, double> > gp = e->getGaussPoints() ;
	Matrix B = geval(f.first.f, e, gp[0].first.x, gp[0].first.y, gp[0].first.z, f.first.transpose)*gp[0].second ;
	for(size_t i = 1 ; i < gp.size() ; i++)
	{
		B += geval(f.first.f, e, gp[i].first.x, gp[i].first.y, gp[i].first.z, f.first.transpose)*gp[i].second ;
	}
	
	return B*f.second ;
}

Matrix VirtualMachine::ieval(const VGtM &f, const IntegrableEntity *e)
{

	std::valarray< std::pair<Point, double> > gp = e->getGaussPoints() ;
	Matrix B = gveval(f.first.f, e, gp[0].first.x, gp[0].first.y, gp[0].first.z, f.first.transpose)*gp[0].second ;
	for(size_t i = 1 ; i < gp.size() ; i++)
	{
		B += geval(f.first.f, e, gp[i].first.x, gp[i].first.y, gp[i].first.z, f.first.transpose)*gp[i].second ;
	}
	
	return B*f.second ;
}



Matrix VirtualMachine::ieval(const GtMtG &f, const IntegrableEntity *e)
{

	std::valarray< std::pair<Point, double> > gp = e->getGaussPoints() ;
	std::valarray<Matrix> Jinv(gp.size()) ;
	for(size_t i = 0 ; i < gp.size() ; i++)
	{
		Jinv[i] = e->getInverseJacobianMatrix(gp[i].first) ;
	}
		
		return ieval(f, gp, Jinv) ;

}

double VirtualMachine::ieval(const VGtMtVG &f, const IntegrableEntity *e)
{

	std::valarray< std::pair<Point, double> > gp = e->getGaussPoints() ;
	std::valarray<Matrix> Jinv(gp.size()) ;
	for(size_t i = 0 ; i < gp.size() ; i++)
	{
		Jinv[i] = e->getInverseJacobianMatrix(gp[i].first) ;
	}
		
	return ieval(f, gp, Jinv) ;

}

Matrix VirtualMachine::ieval(const GtMtG &f, const std::valarray< std::pair<Point, double> > &gp, const std::valarray<Matrix> &Jinv)
{

	Matrix B = geval(f.first.f, Jinv[0], gp[0].first.x, gp[0].first.y, gp[0].first.z, f.first.transpose);
	Matrix B_ = geval(f.third.f, Jinv[0], gp[0].first.x, gp[0].first.y, gp[0].first.z, f.third.transpose) ;
	
	Matrix r_  = B*f.second*B_ ;
	Matrix ret = r_ * gp[0].second;
	
	for(size_t i = 1 ; i < gp.size() ; i++)
	{
		B  = geval(f.first.f, Jinv[i], gp[i].first.x, gp[i].first.y, gp[i].first.z, f.first.transpose) ;
		B_ = geval(f.third.f, Jinv[i], gp[i].first.x, gp[i].first.y, gp[i].first.z, f.third.transpose) ;
		Matrix r_  = B*f.second*B_ ;
		ret += r_ * gp[i].second ;
	}
	
	return ret ;
}

double VirtualMachine::ieval(const VGtMtVG &f, const std::valarray< std::pair<Point, double> > &gp, const std::valarray<Matrix> &Jinv)
{

	Matrix B = gveval(f.first.f, Jinv[0], gp[0].first.x, gp[0].first.y, gp[0].first.z, f.first.transpose);
	Matrix B_ = gveval(f.third.f, Jinv[0], gp[0].first.x, gp[0].first.y, gp[0].first.z, f.third.transpose) ;
	
	Matrix r_  = B*f.second*B_ ;
	double ret = r_[0][0] * gp[0].second;
	
	for(size_t i = 1 ; i < gp.size() ; i++)
	{
		B  = gveval(f.first.f, Jinv[i], gp[i].first.x, gp[i].first.y, gp[i].first.z, f.first.transpose) ;
		B_ = gveval(f.third.f, Jinv[i], gp[i].first.x, gp[i].first.y, gp[i].first.z, f.third.transpose) ;
		Matrix r_  = B*f.second*B_ ;
		ret += r_[0][0] * gp[i].second ;
	}
	
	return ret ;
}

Vector VirtualMachine::ieval(const GtV &f, const std::valarray< std::pair<Point, double> > &gp, const std::valarray<Matrix> &Jinv)
{

	Matrix B = geval(f.first.f, Jinv[0], gp[0].first.x, gp[0].first.y, gp[0].first.z, f.first.transpose);
	
	Vector temp(double(0),f.second.size()) ;
	for(size_t i = 0 ; i < temp.size() ; i++)
		temp[i] = f.second[i] ;
	
	Vector r_  = B*temp ;
	Vector ret = r_ * gp[0].second;
	
	for(size_t i = 1 ; i < gp.size() ; i++)
	{
		B  = geval(f.first.f, Jinv[i], gp[i].first.x, gp[i].first.y, gp[i].first.z, f.first.transpose) ;
		for(size_t j = 0 ; j < temp.size() ; j++)
		{
			temp[j] = f.second[j] ;
		}
		Vector r_  = B*temp ;
		ret += r_ * gp[i].second ;
	}
	
	return ret ;
}

double VirtualMachine::ieval(const VGtV &f, const std::valarray< std::pair<Point, double> > &gp, const std::valarray<Matrix> &Jinv)
{

	Matrix B = gveval(f.first.f, Jinv[0], gp[0].first.x, gp[0].first.y, gp[0].first.z, f.first.transpose);
	
	Vector temp(B.numCols()) ;
	for(size_t i = 0 ; i < temp.size() ; i++)
		temp[i] = f.second[i] ;
	
	Vector r_  = B*temp ;
	double ret = r_[0] * gp[0].second;
	
	for(size_t i = 1 ; i < gp.size() ; i++)
	{
		B  = gveval(f.first.f, Jinv[i], gp[i].first.x, gp[i].first.y, gp[i].first.z, f.first.transpose) ;
		for(size_t j = 0 ; j < temp.size() ; j++)
			temp[j] = f.second[i*temp.size()+j] ;
		Vector r_  = B*temp ;
		ret += r_[0] * gp[i].second ;
	}
	
	return ret ;
}

std::vector<Point> VirtualMachine::allHints(const Function &f0, const Function &f1,const IntegrableEntity *e )
{
	std::vector<Point> hints(0) ;
	
	if (e!= NULL && f0.hasIntegrationHint() || f1.hasIntegrationHint())
	{
		hints.push_back(e->getBoundingPoint(0)) ;
		hints.push_back(e->getBoundingPoint(e->getBoundingPoints().size()/3)); 
		hints.push_back(e->getBoundingPoint(2*e->getBoundingPoints().size()/3)) ;
// 		for(size_t i = 0 ; i < t->getBoundingPoints()->size() ; i++)
// 		{
// 			if(i%t->getOrder() != 0)
// 				hints.push_back(Point(*t->getBoundingPoint(i))) ;
// 		}
		if(f0.hasIntegrationHint())
		{
			for(size_t i = 0 ; i<  f0.getIntegrationHint().size() ; i++)
			{
				hints.push_back(f0.getIntegrationHint(i)) ;
			}
		}
		if(f1.hasIntegrationHint())
		{
			for(size_t i = 0 ; i<  f1.getIntegrationHint().size() ; i++)
			{
				hints.push_back(f1.getIntegrationHint(i)) ;
			}
		}
	}
	
	return hints ;
}

std::vector<Point> VirtualMachine::allHints(const Function &f0,const IntegrableEntity *e )
{
	std::vector<Point> hints(0) ;
	

	if (e!= NULL && f0.hasIntegrationHint() )
	{
		hints.push_back(e->getBoundingPoint(0)) ;
		hints.push_back(e->getBoundingPoint(e->getBoundingPoints().size()/3)); 
		hints.push_back(e->getBoundingPoint(2*e->getBoundingPoints().size()/3)) ;
		for(size_t i = 0 ; i < e->getBoundingPoints().size() ; i++)
		{
			if(i%e->getOrder() != 0)
				hints.push_back(Point(e->getBoundingPoint(i))) ;
		}
		if(f0.hasIntegrationHint())
		{
			for(size_t i = 0 ; i<  f0.getIntegrationHint().size() ; i++)
			{
				hints.push_back(f0.getIntegrationHint(i)) ;
			}
		}
	}
	
	return hints ;
}

Matrix VirtualMachine::ieval(const GtFMtG &f, const IntegrableEntity *e)
{

		std::valarray< std::pair<Point, double> > gp = e->getGaussPoints() ;
	std::valarray<Matrix> Jinv(gp.size()) ;
		for(size_t i = 0 ; i < gp.size() ; i++)
		{
			Jinv[i] = e->getInverseJacobianMatrix(gp[i].first) ;
		}
	
		return ieval(f, gp, Jinv) ;
	
}

Matrix VirtualMachine::ieval(const GtFMtG &f, const std::valarray< std::pair<Point, double> > &gp, const std::valarray<Matrix> &Jinv)
{

	Matrix B = geval(f.first.f, Jinv[0], gp[0].first.x, gp[0].first.y, gp[0].first.z, f.first.transpose);
	Matrix B_ = geval(f.third.f, Jinv[0], gp[0].first.x, gp[0].first.y, gp[0].first.z, f.third.transpose) ;
	
	Matrix r_  = B*eval(f.second,gp[0].first )*B_ ;
	Matrix ret = r_ * gp[0].second;
	
	for(size_t i = 1 ; i < gp.size() ; i++)
	{
		B  = geval(f.first.f, Jinv[i], gp[i].first.x, gp[i].first.y, gp[i].first.z, f.first.transpose) ;
		B_ = geval(f.third.f, Jinv[i], gp[i].first.x, gp[i].first.y, gp[i].first.z, f.third.transpose) ;
		Matrix r_  = B*eval(f.second,gp[i].first )*B_ ;
		ret += r_ * gp[i].second ;
	}
	
	return ret ;
}

double VirtualMachine::ieval(const Differential & d, const std::valarray< std::pair<Point, double> > &gp, const std::valarray<Matrix> &Jinv)
{
	switch(d.v)
	{
	case XI:
		{
			double ret = 0 ;
			for(size_t i = 0 ; i < gp.size() ; i++)
			{
				for(size_t j = 0 ; j < Jinv[i].numCols() ; j++)
				{
					ret += deval(d.f, (Variable)j, gp[i].first)*Jinv[i][0][j]*gp[i].second ;
				}
			}
			
			return ret ;
		}
	case ETA:
		{
			double ret = 0 ;
			for(size_t i = 0 ; i < gp.size() ; i++)
			{
				for(size_t j = 0 ; j < Jinv[i].numCols() ; j++)
				{
					ret += deval(d.f, (Variable)j, gp[i].first)*Jinv[i][1][j]*gp[i].second ;
				}
			}
			
			return ret ;
		}
	case ZETA:
		{
			double ret = 0 ;
			for(size_t i = 0 ; i < gp.size() ; i++)
			{
				for(size_t j = 0 ; j < Jinv[i].numCols() ; j++)
				{
					ret += deval(d.f, (Variable)j, gp[i].first)*Jinv[i][2][j]*gp[i].second ;
				}
			}
			
			return ret ;
		}
	case TIME_VARIABLE:
		{
			double ret = 0 ;
			for(size_t i = 0 ; i < gp.size() ; i++)
			{
				ret += deval(d.f, TIME_VARIABLE, gp[i].first)*gp[i].second ;
			}
			
			return ret ;
		}
	default:
		{
			std::cout << "differential operator not implemented" << std::endl ;
			return 0 ;
		}
	}
}


double VirtualMachine::ieval(const DtF & d, const std::valarray< std::pair<Point, double> > &gp, const std::valarray<Matrix> &Jinv)
{
	switch(d.d.v)
	{
	case XI:
		{
			double ret = 0 ;
			for(size_t i = 0 ; i < gp.size() ; i++)
			{
				for(size_t j = 0 ; j < Jinv[i].numCols() ; j++)
				{
					ret += deval(d.d.f, (Variable)j, gp[i].first)*
						   deval(d.f, (Variable)j, gp[i].first)*Jinv[i][0][j]*gp[i].second ;
				}
			}
			
			return ret ;
		}
	case ETA:
		{
			double ret = 0 ;
			for(size_t i = 0 ; i < gp.size() ; i++)
			{
				for(size_t j = 0 ; j < Jinv[i].numCols() ; j++)
				{
					ret += deval(d.d.f, (Variable)j, gp[i].first)*
						deval(d.f, (Variable)j, gp[i].first)*Jinv[i][1][j]*gp[i].second ;
				}
			}
			
			return ret ;
		}
	case ZETA:
		{
			double ret = 0 ;
			for(size_t i = 0 ; i < gp.size() ; i++)
			{
				for(size_t j = 0 ; j < Jinv[i].numCols() ; j++)
				{
					ret += deval(d.d.f, (Variable)j, gp[i].first)*
						deval(d.f, (Variable)j, gp[i].first)*Jinv[i][2][j]*gp[i].second ;
				}
			}
			
			return ret ;
		}
	case TIME_VARIABLE:
		{
			double ret = 0 ;
			for(size_t i = 0 ; i < gp.size() ; i++)
			{
				ret += deval(d.d.f, TIME_VARIABLE, gp[i].first)*
					deval(d.f, TIME_VARIABLE, gp[i].first)*gp[i].second ;
			}
			
			return ret ;
		}
	default:
		{
			std::cout << "differential operator not implemented" << std::endl ;
			return 0 ;
		}
	}
}

double VirtualMachine::ieval(const Differential & d, const IntegrableEntity *e)
{
	std::valarray< std::pair<Point, double> > gp = e->getGaussPoints() ;
			
	switch(d.v)
	{
	case XI:
		{
			double ret = 0 ;
			
			for(size_t i = 0 ; i < gp.size() ; i++)
			{
				Matrix Jinv = e->getInverseJacobianMatrix(gp[i].first);
				for(size_t j = 0 ; j < Jinv.numCols() ; j++)
				{
					ret += deval(d.f, (Variable)j, gp[i].first)*Jinv[0][j]*gp[i].second ;
				}
			}
			
			return ret ;
		}
	case ETA:
		{
			double ret = 0 ;
			
			for(size_t i = 0 ; i < gp.size() ; i++)
			{
				Matrix Jinv = e->getInverseJacobianMatrix(gp[i].first);
				for(size_t j = 0 ; j < Jinv.numCols() ; j++)
				{
					ret += deval(d.f, (Variable)j, gp[i].first)*Jinv[1][j]*gp[i].second ;
				}
			}
			
			return ret ;
		}
	case ZETA:
		{
			double ret = 0 ;
			
			for(size_t i = 0 ; i < gp.size() ; i++)
			{
				Matrix Jinv = e->getInverseJacobianMatrix(gp[i].first);
				for(size_t j = 0 ; j < Jinv.numCols() ; j++)
				{
					ret += deval(d.f, (Variable)j, gp[i].first)*Jinv[2][j]*gp[i].second ;
				}
			}
			
			return ret ;
		}
	case TIME_VARIABLE:
		{
			double ret = 0 ;
			for(size_t i = 0 ; i < gp.size() ; i++)
			{
				ret += deval(d.f, TIME_VARIABLE, gp[i].first)*gp[i].second ;
			}
			
			return ret ;
		}
	default:
		{
			std::cout << "differential operator not implemented" << std::endl ;
			return 0 ;
		}
	}
}

double VirtualMachine::ieval(const DtF & d, const IntegrableEntity *e)
{
	std::valarray< std::pair<Point, double> > gp = e->getGaussPoints() ;
	
	switch(d.d.v)
	{
	case XI:
		{
			double ret = 0 ;
			
			for(size_t i = 0 ; i < gp.size() ; i++)
			{
				Matrix Jinv = e->getInverseJacobianMatrix(gp[i].first);
				for(size_t j = 0 ; j < Jinv.numCols() ; j++)
				{
					ret += deval(d.d.f, (Variable)j, gp[i].first)*
						   deval(d.f, (Variable)j, gp[i].first)*Jinv[0][j]*gp[i].second ;
				}
			}
			
			return ret ;
		}
	case ETA:
		{
			double ret = 0 ;
			
			for(size_t i = 0 ; i < gp.size() ; i++)
			{
				Matrix Jinv = e->getInverseJacobianMatrix(gp[i].first);
				for(size_t j = 0 ; j < Jinv.numCols() ; j++)
				{
					ret += deval(d.d.f, (Variable)j, gp[i].first)*
						   deval(d.f, (Variable)j, gp[i].first)*Jinv[1][j]*gp[i].second ;
				}
			}
			
			return ret ;
		}
	case ZETA:
		{
			double ret = 0 ;
			
			for(size_t i = 0 ; i < gp.size() ; i++)
			{
				Matrix Jinv = e->getInverseJacobianMatrix(gp[i].first);
				for(size_t j = 0 ; j < Jinv.numCols() ; j++)
				{
					ret += deval(d.d.f, (Variable)j, gp[i].first)*
						   deval(d.f, (Variable)j, gp[i].first)*Jinv[2][j]*gp[i].second ;
				}
			}
			
			return ret ;
		}
	case TIME_VARIABLE:
		{
			double ret = 0 ;
			for(size_t i = 0 ; i < gp.size() ; i++)
			{
				ret += deval(d.f, TIME_VARIABLE, gp[i].first)*
					   deval(d.f, TIME_VARIABLE, gp[i].first)*gp[i].second ;
			}
			
			return ret ;
		}
	default:
		{
			std::cout << "differential operator not implemented" << std::endl ;
			return 0 ;
		}
	}
}


GtFMtG GtFM::operator*(const Mu::Gradient & f) const
{
	return GtFMtG(this->first, this->second, f) ;
}

Mu::GtFM operator *(const Mu::Gradient g, const Mu::FunctionMatrix m)
{
	return Mu::GtFM(g, m) ;
}
