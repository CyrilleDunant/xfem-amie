// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution



#include "vm_base.h"
#include "../mesher/delaunay.h"
#include <limits>


using namespace Mu ;

VirtualMachine::VirtualMachine(){ } ;

double VirtualMachine::eval(const Function &f, const double x, const double y, const double z, const double t, const double u, const double v, const double w) 
{
	stack.memory.reset() ;
	stack.set(x, y, z, t, u, v, w) ;
// 	Context context(stack, x, y, z, t, u, v, w) ;
	
	for(size_t i = 0 ; i < f.size()  ; ++i)
		f.tokenEval(i,stack) ;

// 	if(std::abs(stack[0]) < 1e-6)
// 		return 0 ;
		
	return  *stack.memory.top_pos ;
}

double VirtualMachine::eval(const Function &f, std::vector<std::pair<std::string, double> > vars, const double x, const double y, const double z, const double t, const double u, const double v, const double w) 
{
	stack.memory.reset() ;
	for(size_t i = 0 ; i < vars.size() ; i++)
		stack.memory.setNamedVariable(vars[i].first, vars[i].second);
	
	stack.set(x, y, z, t, u, v, w) ;
// 	Context context(stack, x, y, z, t, u, v, w) ;
	
	for(size_t i = 0 ; i < f.size()  ; ++i)
		f.tokenEval(i,stack) ;

// 	if(std::abs(stack[0]) < 1e-6)
// 		return 0 ;
		
	return  *stack.memory.top_pos ;
}

Vector VirtualMachine::eval(const Function &f, const GaussPointArray &gp)
{
	if(f.precalculated(gp))
		return f.getPrecalculatedValue(gp) ;
	
	Vector ret(gp.gaussPoints.size()) ;
	for(size_t i = 0 ; i < ret.size() ; i++)
		ret[i] = eval(f, gp.gaussPoints[i].first.x,
				 gp.gaussPoints[i].first.y,
				 gp.gaussPoints[i].first.z,
				 gp.gaussPoints[i].first.t) ;

	return ret ;
} 

Vector VirtualMachine::oeval(const Function &f, const GaussPointArray &gp, const Point & offset)
{
	Vector ret(gp.gaussPoints.size()) ;
	for(size_t i = 0 ; i < ret.size() ; i++)
		ret[i] = eval(f, gp.gaussPoints[i].first.x+offset.x,
				 gp.gaussPoints[i].first.y+offset.y,
				 gp.gaussPoints[i].first.z+offset.z,
				 gp.gaussPoints[i].first.t+offset.t) ;

	return ret ;
}

double VirtualMachine::eval(const std::vector<RefCountedToken> &f, const double & x, const double & y, const double & z, const double &t, const double &u, const double &v, const double &w) 
{
	stack.memory.reset() ;
	stack.set(x, y, z, t, u, v, w) ;

	size_t size = f.size() ;
	for(size_t i = 0 ; i < size  ; ++i)
		f[i]->eval(stack) ;

// 	if(std::abs(stack[0]) < 1e-6)
// 		return 0 ;
		
	 return  *stack.memory.top_pos ;
 }

double VirtualMachine::eval(const Function &f, const Point & p, const Point & p_) 
{
	return eval(f, p.x, p.y, p.z,p.t,p_.x, p_.y, p_.z) ;
}


double VirtualMachine::eval(const Function &f, const Point *p, const Point * p_) 
{
	if(p_)
		return eval(f, p->x, p->y, p->z, p->t, p_->x, p_->y, p_->z) ;
	return eval(f, p->x, p->y, p->z, p->t) ;
}


void VirtualMachine::print(const Function &f) const
{
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

Matrix VirtualMachine::eval(const FunctionMatrix &f, const Point & p, const Point & p_) 
{
	return eval(f, p.x, p.y, p.z, p.t, p_.x, p_.y, p_.z) ;
}

Matrix VirtualMachine::eval(const FunctionMatrix &f, const Point * p, const Point *p_) 
{
	if(p_)
		return eval(f, p->x, p->y, p->z, p->t, p_->x, p_->y, p_->z) ;
	return eval(f, p->x, p->y, p->z, p->t) ;
}

double VirtualMachine::deval(const Function &f, const Point&p,  const double x, const double y, const double z, const double t , const double eps , bool normed)
{
	double d = deval(f,XI,x,y,z,t,0,0,0,eps)*p.x + deval(f,ETA,x,y,z,t,0,0,0,eps)*p.y + deval(f,ZETA,x,y,z,t,0,0,0,eps)*p.z ;
	if(normed)
		return d ;
	else
		return d/p.norm() ;
}

Matrix VirtualMachine::deval(const FunctionMatrix &f,const Point&p,  const double x, const double y, const double z, const double t , const double eps, bool normed )
{
	Matrix ret(f.numRows(), f.numCols()) ;
	
	for(size_t i = 0 ; i < f.numRows() ; i++)
	{
		for(size_t j = 0 ; j < f.numCols() ; j++)
		{
			ret[i][j] = deval(f[i][j], p, x, y, z,t, eps, normed) ;
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

Vector VirtualMachine::deval(const Function&f, const Variable v_, const GaussPointArray & gp, const double eps)
{
	if(f.precalculated(gp, v_))
		return f.getPrecalculatedValue(gp, v_) ;
	
	Point poffset ;
	Point noffset ;
	
	Vector ret(double(0), gp.gaussPoints.size()) ;
	double h = sqrt(std::numeric_limits<double>::epsilon()) ;
	switch(v_)
	{
	case ONE : 
		{
			return ret ;
		}
	case XI : 
		{
			volatile double temp = poffset.x+h ;
			h = temp -poffset.x ;
			poffset.x = h ;
			noffset.x = -h ;
			break ;
		}
	case ETA:
		{
			volatile double temp = poffset.y+h ;
			h = temp -poffset.y ;
			poffset.y = h ;
			noffset.y = -h ;
			break ;
		}
	case ZETA:
		{
			volatile double temp = poffset.z+h ;
			h = temp -poffset.z ;
			poffset.z = h ;
			noffset.z = -h ;
			break ;
		}
	case TIME_VARIABLE : 
		{
			volatile double temp = poffset.t+h ;
			h = temp -poffset.t ;
			poffset.t = h ;
			noffset.t = -h ;
			break ;
		}
	default:
		return ret ;
	}
	
	return (oeval(f, gp, poffset)-oeval(f, gp, noffset))/(2.*h) ;
}

Vector VirtualMachine::ddeval(const Function&f, const Variable v_0, const Variable v_1, const GaussPointArray & gp, const double eps)
{
	Vector ret(gp.gaussPoints.size());
	for(size_t i = 0 ; i < ret.size() ; i++)
	{
		ret[i] = ddeval(f, v_0, v_1, gp.gaussPoints[i].first.x, gp.gaussPoints[i].first.y, gp.gaussPoints[i].first.z, gp.gaussPoints[i].first.t, 0, 0, 0, eps ) ;
	}
	
	return ret ;
}

double VirtualMachine::ddeval(const Function &f, const Variable v_0, const Variable v_1,  const double x, const double y , const double z, const double t, const double u, const double v, const double w, const double eps) 
{
	switch(v_0)
	{
		case ONE : 
		{
			switch(v_1)
			{
			case ONE : 
				{
					return 0 ;
				}
			case XI : 
				{
					return ( eval(f, x+eps, y, z, t, u, v, w) - eval(f, x-eps, y, z, t, u, v, w))/(2.*eps) ;
				}
			case ETA:
				{
					return ( eval(f, x, y+eps, z, t, u, v, w) - eval(f, x, y-eps, z, t, u, v, w))/(2.*eps) ;
				}
			case ZETA:
				{
					return ( eval(f, x, y, z+eps, t, u, v, w) - eval(f, x, y, z-eps, t, u, v, w))/(2.*eps) ;
				}
			case TIME_VARIABLE : 
				{
					return (eval(f, x, y, z, t+eps, u, v, w) - eval(f, x, y, z, t-eps, u, v, w))/(2.*eps) ;
				}
			case U_VARIABLE:
				{
					return ( eval(f, x, y, z, t, u+eps, v, w) - eval(f, x, y, z, t, u-eps, v, w))/(2.*eps) ;
				}
			case V_VARIABLE:
				{
					return ( eval(f, x, y, z, t, u, v+eps, w) - eval(f, x, y, z, t, u, v-eps, w))/(2.*eps) ;
				}
			case W_VARIABLE:
				{
					return ( eval(f, x, y, z, t, u, v, w+eps) - eval(f, x, y, z, t, u, v, w-eps))/(2.*eps) ;
				}
			}
		}
		case XI : 
		{
			switch(v_1)
			{
			case ONE : 
				{
					return ( eval(f, x+eps, y, z, t, u, v, w) - eval(f, x-eps, y, z, t, u, v, w))/(2.*eps) ;
				}
			case XI : 
				{
					return ( eval(f, x+eps+eps, y, z, t, u, v, w) - 2.*eval(f, x, y, z, t, u, v, w) + eval(f, x-eps-eps, y, z, t, u, v, w))/(4.*eps*eps) ;
				}
			case ETA:
				{
					return ( eval(f, x+eps, y+eps, z, t, u, v, w) - eval(f, x-eps, y+eps, z, t, u, v, w)- eval(f, x+eps, y-eps, z, t, u, v, w) + eval(f, x-eps, y-eps, z, t, u, v, w))/(4.*eps*eps) ;
				}
			case ZETA:
				{
					return ( eval(f, x+eps, y, z+eps, t, u, v, w) - eval(f, x-eps, y, z+eps, t, u, v, w)- eval(f, x+eps, y, z-eps, t, u, v, w) + eval(f, x-eps, y, z-eps, t, u, v, w))/(4.*eps*eps) ;
				}
			case TIME_VARIABLE : 
				{
					return( eval(f, x+1e-4, y, z, t+1e-4, u, v, w) - eval(f, x-1e-4, y, z, t+1e-4, u, v, w)- eval(f, x+1e-4, y, z, t-1e-4, u, v, w) + eval(f, x-1e-4, y, z, t-1e-4, u, v, w))/(4.*eps*eps) ;
				}
			case U_VARIABLE:
				{
					return ( eval(f, x+eps, y, z, t, u+eps, v, w) - eval(f, x-eps, y, z, t, u+eps, v, w)- eval(f, x+eps, y, z, t, u-eps, v, w) + eval(f, x-eps, y, z, t, u-eps, v, w))/(4.*eps*eps) ;
				}
			case V_VARIABLE:
				{
					return ( eval(f, x+eps, y, z, t, u, v+eps, w) - eval(f, x-eps, y, z, t, u, v+eps, w)- eval(f, x+eps, y, z, t, u, v-eps, w) + eval(f, x-eps, y, z, t, u, v-eps, w))/(4.*eps*eps) ;
				}
			case W_VARIABLE:
				{
					return ( eval(f, x+eps, y, z, t, u, v, w+eps) - eval(f, x-eps, y, z, t, u, v, w+eps)- eval(f, x+eps, y, z, t, u, v, w-eps) + eval(f, x-eps, y, z, t, u, v, w-eps))/(4.*eps*eps) ;
				}
			}
		}
		case ETA:
		{
			switch(v_1)
			{
			case ONE : 
				{
					return ( eval(f, x, y+eps, z, t, u, v, w) - eval(f, x, y-eps, z, t, u, v, w))/(2.*eps) ;
				}
			case XI : 
				{
					return ( eval(f, x+eps, y+eps, z, t, u, v, w) - eval(f, x-eps, y+eps, z, t, u, v, w)- eval(f, x+eps, y-eps, z, t, u, v, w) + eval(f, x-eps, y-eps, z, t, u, v, w))/(4.*eps*eps) ;
				}
			case ETA:
				{
					return ( eval(f, x, y+eps+eps, z, t, u, v, w) - 2.*eval(f, x, y, z, t, u, v, w) + eval(f, x, y-eps-eps, z, t, u, v, w))/(4.*eps*eps) ;

				}
			case ZETA:
				{
					return ( eval(f, x, y+eps, z+eps, t, u, v, w) - eval(f, x, y-eps, z+eps, t, u, v, w)- eval(f, x, y+eps, z-eps, t, u, v, w) + eval(f, x, y-eps, z-eps, t, u, v, w))/(4.*eps*eps) ;
				}
			case TIME_VARIABLE : 
				{
					return ( eval(f, x, y+1e-3, z, t+1e-3, u, v, w) - eval(f, x, y-1e-3, z, t+1e-3, u, v, w)- eval(f, x, y+1e-3, z, t-1e-3, u, v, w) + eval(f, x, y-1e-3, z, t-1e-3, u, v, w))/(4.*1e-3*1e-3) ;
				}
			case U_VARIABLE:
				{
					return ( eval(f, x, y+eps, z, t, u+eps, v, w) - eval(f, x, y-eps, z, t, u+eps, v, w)- eval(f, x, y+eps, z, t, u-eps, v, w) + eval(f, x, y-eps, z, t, u-eps, v, w))/(4.*eps*eps) ;
				}
			case V_VARIABLE:
				{
					return ( eval(f, x, y+eps, z, t, u, v+eps, w) - eval(f, x, y-eps, z, t, u, v+eps, w)- eval(f, x, y+eps, z, t, u, v-eps, w) + eval(f, x, y-eps, z, t, u, v-eps, w))/(4.*eps*eps) ;
				}
			case W_VARIABLE:
				{
					return ( eval(f, x, y+eps, z, t, u, v, w+eps) - eval(f, x, y-eps, z, t, u, v, w+eps)- eval(f, x, y+eps, z, t, u, v, w-eps) + eval(f, x, y-eps, z, t, u, v, w-eps))/(4.*eps*eps) ;
				}
			}
		}
		case ZETA:
		{
			switch(v_1)
			{
			case ONE : 
				{
					return ( eval(f, x, y, z+eps, t, u, v, w) - eval(f, x, y, z-eps, t, u, v, w))/(2.*eps) ;
				}
			case XI : 
				{
					return ( eval(f, x+eps, y+eps, z, t, u, v, w) - eval(f, x-eps, y, z+eps, t, u, v, w)- eval(f, x+eps, y, z-eps, t, u, v, w) + eval(f, x-eps, y, z-eps, t, u, v, w))/(4.*eps*eps) ;
				}
			case ETA:
				{

					return ( eval(f, x, y+eps, z+eps, t, u, v, w) - eval(f, x, y-eps, z+eps, t, u, v, w)- eval(f, x, y+eps, z-eps, t, u, v, w) + eval(f, x, y-eps, z-eps, t, u, v, w))/(4.*eps*eps) ;
				}
			case ZETA:
				{
					return ( eval(f, x, y, z+eps+eps, t, u, v, w) - 2.*eval(f, x, y, z, t, u, v, w) + eval(f, x, y, z-eps-eps, t, u, v, w))/(4.*eps*eps) ;
				}
			case TIME_VARIABLE : 
				{
					return( eval(f, x, y, z+eps, t+eps, u, v, w) - eval(f, x, y, z-eps, t+eps, u, v, w)- eval(f, x, y, z+eps, t-eps, u, v, w) + eval(f, x, y, z-eps, t-eps, u, v, w))/(4.*eps*eps) ;
				}
			case U_VARIABLE:
				{
					return ( eval(f, x, y, z+eps, t, u+eps, v, w) - eval(f, x, y, z-eps, t, u+eps, v, w)- eval(f, x, y, z+eps, t, u-eps, v, w) + eval(f, x, y, z-eps, t, u-eps, v, w))/(4.*eps*eps) ;
				}
			case V_VARIABLE:
				{
					return ( eval(f, x, y, z+eps, t, u, v+eps, w) - eval(f, x, y, z-eps, t, u, v+eps, w)- eval(f, x, y, z+eps, t, u, v-eps, w) + eval(f, x, y, z-eps, t, u, v-eps, w))/(4.*eps*eps) ;
				}
			case W_VARIABLE:
				{
					return ( eval(f, x, y, z+eps, t, u, v, w+eps) - eval(f, x, y, z-eps, t, u, v, w+eps)- eval(f, x, y, z+eps, t, u, v, w-eps) + eval(f, x, y, z-eps, t, u, v, w-eps))/(4.*eps*eps) ;
				}
			}
		}
		case TIME_VARIABLE : 
		{
			switch(v_1)
			{
			case ONE : 
				{
					return ( eval(f, x, y, z, t+eps, u, v, w) - eval(f, x, y, z, t-eps, u, v, w))/(2.*eps) ;
				}
			case XI : 
				{
					return ( eval(f, x+eps, y, z, t+eps, u, v, w) - eval(f, x-eps, y, z, t+eps, u, v, w)- eval(f, x+eps, y, z, t-eps, u, v, w) + eval(f, x-eps, y, z, t-eps, u, v, w))/(4.*eps*eps) ;
				}
			case ETA:
				{
					return( eval(f, x, y+eps, z, t+eps, u, v, w) - eval(f, x, y-eps, z, t+eps, u, v, w)- eval(f, x, y+eps, z, t-eps, u, v, w) + eval(f, x-eps, y, z, t-eps, u, v, w))/(4.*eps*eps) ;
				}
			case ZETA:
				{
					return ( eval(f, x, y, z+eps, t+eps, u, v, w) - eval(f, x, y, z+eps, t-eps, u, v, w)- eval(f, x, y, z-eps, t+eps, u, v, w) + eval(f, x, y, z-eps, t-eps, u, v, w))/(4.*eps*eps) ;
				}
			case TIME_VARIABLE : 
				{
					return ( eval(f, x, y, z, t+eps+eps, u, v, w) - 2.*eval(f, x, y, z, t, u, v, w) + eval(f, x, y, z, t-eps-eps, u, v, w))/(4.*eps*eps) ;

				}
			case U_VARIABLE:
				{
					return ( eval(f, x, y, z, t+eps, u+eps, v, w) - eval(f, x, y, z, t-eps, u+eps, v, w)- eval(f, x, y, z, t+eps, u-eps, v, w) + eval(f, x, y, z, t-eps, u-eps, v, w))/(4.*eps*eps) ;
				}
			case V_VARIABLE:
				{
					return ( eval(f, x, y, z, t+eps, u, v+eps, w) - eval(f, x, y, z, t-eps, u, v+eps, w)- eval(f, x, y, z, t+eps, u, v-eps, w) + eval(f, x, y, z, t-eps, u, v-eps, w))/(4.*eps*eps) ;
				}
			case W_VARIABLE:
				{
					return ( eval(f, x, y, z, t+eps, u, v, w+eps) - eval(f, x, y, z, t-eps, u, v, w+eps)- eval(f, x, y, z, t+eps, u, v, w-eps) + eval(f, x, y, z, t-eps, u, v, w-eps))/(4.*eps*eps) ;
				}
			}
		}
		case U_VARIABLE:
		{
			std::cerr << "not implemented (ddeval)" << std::endl ;
			switch(v_1)
			{
			case ONE : 
				{
					return ( eval(f, x, y+eps, z, t, u, v, w) - eval(f, x, y-eps, z, t, u, v, w))/(2.*eps) ;
				}
			case XI : 
				{
					return ( eval(f, x+eps, y+eps, z, t, u, v, w) - eval(f, x-eps, y+eps, z, t, u, v, w)- eval(f, x+eps, y-eps, z, t, u, v, w) + eval(f, x-eps, y-eps, z, t, u, v, w))/(4.*eps*eps) ;
				}
			case ETA:
				{
					return ( eval(f, x, y+eps+eps, z, t, u, v, w) - 2.*eval(f, x, y, z, t, u, v, w) + eval(f, x, y-eps-eps, z, t, u, v, w))/(4.*eps*eps) ;

				}
			case ZETA:
				{
					return ( eval(f, x, y+eps, z+eps, t, u, v, w) - eval(f, x, y-eps, z+eps, t, u, v, w)- eval(f, x, y+eps, z-eps, t, u, v, w) + eval(f, x-eps, y, z-eps, t, u, v, w))/(4.*eps*eps) ;
				}
			case TIME_VARIABLE : 
				{
					return( eval(f, x, y+eps, z, t+eps, u, v, w) - eval(f, x, y-eps, z, t+eps, u, v, w)- eval(f, x, y+eps, z, t-eps, u, v, w) + eval(f, x-eps, y, z, t-eps, u, v, w))/(4.*eps*eps) ;
				}
			case U_VARIABLE:
				{
					return ( eval(f, x, y+eps, z, t, u+eps, v, w) - eval(f, x, y-eps, z, t, u+eps, v, w)- eval(f, x, y+eps, z, t, u-eps, v, w) + eval(f, x, y-eps, z, t, u-eps, v, w))/(4.*eps*eps) ;
				}
			case V_VARIABLE:
				{
					return ( eval(f, x, y+eps, z, t, u, v+eps, w) - eval(f, x, y-eps, z, t, u, v+eps, w)- eval(f, x, y+eps, z, t, u, v-eps, w) + eval(f, x, y-eps, z, t, u, v-eps, w))/(4.*eps*eps) ;
				}
			case W_VARIABLE:
				{
					return ( eval(f, x, y+eps, z, t, u, v, w+eps) - eval(f, x, y-eps, z, t, u, v, w+eps)- eval(f, x, y+eps, z, t, u, v, w-eps) + eval(f, x, y-eps, z, t, u, v, w-eps))/(4.*eps*eps) ;
				}
			}
		}
		case V_VARIABLE:
		{
			std::cerr << "not implemented (ddeval)" << std::endl ;
			
			switch(v_1)
			{
			case ONE : 
				{
					return ( eval(f, x, y+eps, z, t, u, v, w) - eval(f, x, y-eps, z, t, u, v, w))/(2.*eps) ;
				}
			case XI : 
				{
					return ( eval(f, x+eps, y+eps, z, t, u, v, w) - eval(f, x-eps, y+eps, z, t, u, v, w)- eval(f, x+eps, y-eps, z, t, u, v, w) + eval(f, x-eps, y-eps, z, t, u, v, w))/(4.*eps*eps) ;
				}
			case ETA:
				{
					return ( eval(f, x, y+eps+eps, z, t, u, v, w) - 2.*eval(f, x, y, z, t, u, v, w) + eval(f, x, y-eps-eps, z, t, u, v, w))/(4.*eps*eps) ;

				}
			case ZETA:
				{
					return ( eval(f, x, y+eps, z+eps, t, u, v, w) - eval(f, x, y-eps, z+eps, t, u, v, w)- eval(f, x, y+eps, z-eps, t, u, v, w) + eval(f, x-eps, y, z-eps, t, u, v, w))/(4.*eps*eps) ;
				}
			case TIME_VARIABLE : 
				{
					return( eval(f, x, y+eps, z, t+eps, u, v, w) - eval(f, x, y-eps, z, t+eps, u, v, w)- eval(f, x, y+eps, z, t-eps, u, v, w) + eval(f, x-eps, y, z, t-eps, u, v, w))/(4.*eps*eps) ;
				}
			case U_VARIABLE:
				{
					return ( eval(f, x, y+eps, z, t, u+eps, v, w) - eval(f, x, y-eps, z, t, u+eps, v, w)- eval(f, x, y+eps, z, t, u-eps, v, w) + eval(f, x, y-eps, z, t, u-eps, v, w))/(4.*eps*eps) ;
				}
			case V_VARIABLE:
				{
					return ( eval(f, x, y+eps, z, t, u, v+eps, w) - eval(f, x, y-eps, z, t, u, v+eps, w)- eval(f, x, y+eps, z, t, u, v-eps, w) + eval(f, x, y-eps, z, t, u, v-eps, w))/(4.*eps*eps) ;
				}
			case W_VARIABLE:
				{
					return ( eval(f, x, y+eps, z, t, u, v, w+eps) - eval(f, x, y-eps, z, t, u, v, w+eps)- eval(f, x, y+eps, z, t, u, v, w-eps) + eval(f, x, y-eps, z, t, u, v, w-eps))/(4.*eps*eps) ;
				}
			}
		}
		case W_VARIABLE:
		{
			std::cerr << "not implemented (ddeval)" << std::endl ;
			
			switch(v_1)
			{
			case ONE : 
				{
					return ( eval(f, x, y+eps, z, t, u, v, w) - eval(f, x, y-eps, z, t, u, v, w))/(2.*eps) ;
				}
			case XI : 
				{
					return ( eval(f, x+eps, y+eps, z, t, u, v, w) - eval(f, x-eps, y+eps, z, t, u, v, w)- eval(f, x+eps, y-eps, z, t, u, v, w) + eval(f, x-eps, y-eps, z, t, u, v, w))/(4.*eps*eps) ;
				}
			case ETA:
				{
					return ( eval(f, x, y+eps+eps, z, t, u, v, w) - 2.*eval(f, x, y, z, t, u, v, w) + eval(f, x, y-eps-eps, z, t, u, v, w))/(4.*eps*eps) ;
				}
			case ZETA:
				{
					return ( eval(f, x, y+eps, z+eps, t, u, v, w) - eval(f, x, y-eps, z+eps, t, u, v, w)- eval(f, x, y+eps, z-eps, t, u, v, w) + eval(f, x-eps, y, z-eps, t, u, v, w))/(4.*eps*eps) ;
				}
			case TIME_VARIABLE : 
				{
					return( eval(f, x, y+eps, z, t+eps, u, v, w) - eval(f, x, y-eps, z, t+eps, u, v, w)- eval(f, x, y+eps, z, t-eps, u, v, w) + eval(f, x-eps, y, z, t-eps, u, v, w))/(4.*eps*eps) ;
				}
			case U_VARIABLE:
				{
					return ( eval(f, x, y+eps, z, t, u+eps, v, w) - eval(f, x, y-eps, z, t, u+eps, v, w)- eval(f, x, y+eps, z, t, u-eps, v, w) + eval(f, x, y-eps, z, t, u-eps, v, w))/(4.*eps*eps) ;
				}
			case V_VARIABLE:
				{
					return ( eval(f, x, y+eps, z, t, u, v+eps, w) - eval(f, x, y-eps, z, t, u, v+eps, w)- eval(f, x, y+eps, z, t, u, v-eps, w) + eval(f, x, y-eps, z, t, u, v-eps, w))/(4.*eps*eps) ;
				}
			case W_VARIABLE:
				{
					return ( eval(f, x, y+eps, z, t, u, v, w+eps) - eval(f, x, y-eps, z, t, u, v, w+eps)- eval(f, x, y+eps, z, t, u, v, w-eps) + eval(f, x, y-eps, z, t, u, v, w-eps))/(4.*eps*eps) ;
				}
			}
		}
	}
	return 0 ;
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
				double h = std::max(eval(f, x, y, z, t, u, v, w), 1.)*eps ;
				volatile double temp = x+h ;
				h = temp - x ;
				return .5*( eval(f, x+h, y, z, t, u, v, w) - eval(f, x-h, y, z, t, u, v, w))/h ;
			}
		case ETA:
			{
				double h = std::max(eval(f, x, y, z, t, u, v, w), 1.)*eps ;
				volatile double temp = y+h ;
				h = temp - y ;
				return .5*( eval(f, x, y+h, z, t, u, v, w) - eval(f, x, y-h, z, t, u, v, w))/h ;
			}
		case ZETA:
			{
				double h = std::max(eval(f, x, y, z, t, u, v, w), 1.)*eps ;
				volatile double temp = z+h ;
				h = temp - z ;
				return .5*( eval(f, x, y, z+h, t, u, v, w) - eval(f, x, y, z-h, t, u, v, w))/h ;
			}
		case TIME_VARIABLE : 
			{
				double h = std::max(eval(f, x, y, z, t, u, v, w), 1.)*eps ;
				volatile double temp = t+h ;
				h = temp - t ;
				return .5*(eval(f, x, y, z, t+h, u, v, w) - eval(f, x, y, z, t-h, u, v, w))/(h) ;
			}
		case U_VARIABLE:
			{
				double h = std::max(eval(f, x, y, z, t, u, v, w), 1.)*eps ;
				volatile double temp = u+h ;
				h = temp -u ;
				return .5*( eval(f, x, y, z, t, u+h, v, w) - eval(f, x, y, z, t, u-h, v, w))/(h) ;
			}
		case V_VARIABLE:
			{
				double h = std::max(eval(f, x, y, z, t, u, v, w), 1.)*eps ;
				volatile double temp = v+h ;
				h = temp -v ;
				return .5*( eval(f, x, y, z, t, u, v+h, w) - eval(f, x, y, z, t, u, v-h, w))/(h) ;
			}
		case W_VARIABLE:
			{
				double h = std::max(eval(f, x, y, z, t, u, v, w), 1.)*eps ;
				volatile double temp = w+h ;
				h = temp -w ;
				return .5*( eval(f, x, y, z, t, u, v, w+h) - eval(f, x, y, z, t, u, v, w-h))/(h) ;
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

double VirtualMachine::deval(const Function &f, const Variable v, const Point p, const Point p_, const double eps)
{
	return deval(f, v, p.x, p.y, p.z, p.t, p_.x, p_.y, p_.z, eps) ;
}

Matrix VirtualMachine::deval(const FunctionMatrix &f, const Variable v, const Point p, const Point p_, const double eps)
{
	return deval(f, v, p.x, p.y, p.z, p.t, p_.x, p_.y, p_.z, eps) ;
}

double VirtualMachine::ieval(const Function &f, const GaussPointArray &gp)
{
	double ret = 0 ;
	for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
		ret += eval(f, gp.gaussPoints[i].first)*gp.gaussPoints[i].second ;
	
	return ret ;
}

double VirtualMachine::ieval(const Function &f, IntegrableEntity *e)
{

	GaussPointArray gp = e->getGaussPoints() ;
	double ret = 0 ;
	for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
		ret += eval(f, gp.gaussPoints[i].first)*gp.gaussPoints[i].second ;
	
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


double VirtualMachine::ieval(Vector &f, IntegrableEntity *e)
{
	GaussPointArray gp = e->getGaussPoints() ;
	double ret = 0 ;
	for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
		ret += f[i]*gp.gaussPoints[i].second ;
	
	return ret ;
}

Matrix VirtualMachine::ieval(const FunctionMatrix &f, IntegrableEntity *e)
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

Matrix VirtualMachine::ieval(const FMtMtFM &f, IntegrableEntity *e)
{
	Matrix a(ieval(f.first, e)) ;
	
	Matrix c(ieval(f.third, e)) ;
	
	return a*(f.second)*c ;
}


Matrix VirtualMachine::geval(const Function &f, IntegrableEntity *e, const std::vector<Variable> & vars, const double x, const double y, const double z, const double t, bool transpose)
{
	Matrix Jinv ; e->getInverseJacobianMatrix(Point(x,y,z,t), Jinv) ;
	return geval(f, Jinv, vars, x, y, z,t, transpose) ;
}

Matrix VirtualMachine::geval(const Gradient &f, IntegrableEntity *e, const std::vector<Variable> & vars, const double x, const double  y, const double z, const double t)
{
	return geval(f.f, e,vars, x, y, z,t, f.transpose) ;
}

Matrix VirtualMachine::geval(const Function &f, const Matrix & m, const std::vector<Variable> & vars, const Point& p, bool transpose )
{
	return geval(f, m,vars, p.x, p.y, p.z, p.t, transpose) ;
}

std::valarray<Matrix> VirtualMachine::geval(const Function &f, const std::valarray<Matrix> & m, const std::vector<Variable> & var, const GaussPointArray & gp, bool transpose )
{
	if(var.size() == 3 )
	{
		if(var[2] == TIME_VARIABLE)
		{
			if(transpose)
			{
				std::valarray<Matrix> ret(Matrix(3,2), gp.gaussPoints.size()) ;
				
				Vector dxi = deval(f, var[0], gp) ;
				Vector deta = deval(f, var[1], gp) ;
				
				for(size_t i = 0 ; i < ret.size() ; i++)
				{
					ret[i][0][0] = dxi[i]*m[i][0][0] + deta[i]*m[i][0][1] ;
					ret[i][1][1] = dxi[i]*m[i][1][0] + deta[i]*m[i][1][1] ;
					ret[i][2][0] = ret[i][1][1] ;
					ret[i][2][1] = ret[i][0][0] ;
				}
				
				return ret ;
			}
			else
			{
				std::valarray<Matrix> ret(Matrix(2,3), gp.gaussPoints.size()) ;

				Vector dxi = deval(f, var[0], gp) ;
				Vector deta = deval(f, var[1], gp) ;
				
				for(size_t i = 0 ; i < ret.size() ; i++)
				{
					ret[i][0][0] = dxi[i]*m[i][0][0] + deta[i]*m[i][0][1] ;
					ret[i][1][1] = dxi[i]*m[i][1][0] + deta[i]*m[i][1][1] ;
					ret[i][0][2] = ret[i][1][1] ;
					ret[i][1][2] = ret[i][0][0] ;
				}
				return ret ;
			}
		}
		else
		{
			if(transpose)
			{
				std::valarray<Matrix> ret(Matrix(6,3), gp.gaussPoints.size()) ;
				
				Vector dxi = deval(f, var[0], gp) ;
				Vector deta = deval(f, var[1], gp) ;
				Vector dzeta = deval(f, var[2],gp) ;

				for(size_t i = 0 ; i < ret.size() ; i++)
				{
					ret[i][0][0] = dxi[i]*m[i][0][0] + deta[i]*m[i][0][1] + dzeta[i]*m[i][0][2] ;
					ret[i][1][1] = dxi[i]*m[i][1][0] + deta[i]*m[i][1][1] + dzeta[i]*m[i][1][2] ;
					ret[i][2][2] = dxi[i]*m[i][2][0] + deta[i]*m[i][2][1] + dzeta[i]*m[i][2][2] ;
					ret[i][3][1] = ret[i][2][2] ;
					ret[i][3][2] = ret[i][1][1] ;
					ret[i][4][0] = ret[i][2][2] ;
					ret[i][4][2] = ret[i][0][0] ;
					ret[i][5][0] = ret[i][1][1] ;
					ret[i][5][1] = ret[i][0][0] ;
				}
				
				return ret ;
			}
			else
			{
				std::valarray<Matrix> ret(Matrix(3,6), gp.gaussPoints.size()) ;
				
				Vector dxi = deval(f, var[0],gp) ;
				Vector deta = deval(f, var[1],gp) ;
				Vector dzeta = deval(f, var[2], gp) ;
				
				for(size_t i = 0 ; i < ret.size() ; i++)
				{
					ret[i][0][0] = dxi[i]*m[i][0][0] + deta[i]*m[i][0][1] + dzeta[i]*m[i][0][2];
					ret[i][1][1] = dxi[i]*m[i][1][0] + deta[i]*m[i][1][1] + dzeta[i]*m[i][1][2];
					ret[i][2][2] = dxi[i]*m[i][2][0] + deta[i]*m[i][2][1] + dzeta[i]*m[i][2][2];
					ret[i][1][3] = ret[i][2][2] ;
					ret[i][2][3] = ret[i][1][1] ;
					ret[i][0][4] = ret[i][2][2] ;
					ret[i][2][4] = ret[i][0][0] ;
					ret[i][0][5] = ret[i][1][1] ;
					ret[i][1][5] = ret[i][0][0] ;
				}
				
				return ret ;
			}
		}
	}
	else if (var.size() == 2)
	{
		if(transpose)
		{
			
			std::valarray<Matrix> ret(Matrix(3,2), gp.gaussPoints.size()) ;
			
			Vector dxi = deval(f, var[0], gp) ;
			Vector deta = deval(f, var[1], gp) ;
			
			for(size_t i = 0 ; i < ret.size() ; i++)
			{
				ret[i][0][0] = dxi[i]*m[i][0][0] + deta[i]*m[i][0][1] ;
				ret[i][1][1] = dxi[i]*m[i][1][0] + deta[i]*m[i][1][1] ;
				ret[i][2][0] = ret[i][1][1] ;
				ret[i][2][1] = ret[i][0][0] ;
			}
			
			return ret ;
		}
		else
		{
			std::valarray<Matrix> ret(Matrix(2,3), gp.gaussPoints.size()) ;
			
			Vector dxi = deval(f, var[0], gp) ;
			Vector deta = deval(f, var[1], gp) ;
			
			for(size_t i = 0 ; i < ret.size() ; i++)
			{
				ret[i][0][0] = dxi[i]*m[i][0][0] + deta[i]*m[i][0][1] ;
				ret[i][1][1] = dxi[i]*m[i][1][0] + deta[i]*m[i][1][1] ;
				ret[i][0][2] = ret[i][1][1] ;
				ret[i][1][2] = ret[i][0][0] ;
			}
			
			return ret ;
		}
	}
	else if(var.size() == 4 )
	{
		if(transpose)
		{
			std::valarray<Matrix> ret(Matrix(6,3), gp.gaussPoints.size()) ;
			
			Vector dxi = deval(f, var[0], gp) ;
			Vector deta = deval(f, var[1], gp) ;
			Vector dzeta = deval(f, var[2], gp) ;
			Vector dteta = deval(f, var[4], gp) ;
			
			for(size_t i = 0 ; i < ret.size() ; i++)
			{
				ret[i][0][0] = dxi[i]*m[i][0][0] + deta[i]*m[i][0][1] + dzeta[i]*m[i][0][2] + dteta[i]*m[i][0][3];
				ret[i][1][1] = dxi[i]*m[i][1][0] + deta[i]*m[i][1][1] + dzeta[i]*m[i][1][2] + dteta[i]*m[i][1][3];
				ret[i][2][2] = dxi[i]*m[i][2][0] + deta[i]*m[i][2][1] + dzeta[i]*m[i][2][2] + dteta[i]*m[i][2][3];
				ret[i][3][1] = ret[i][2][2] ;
				ret[i][3][2] = ret[i][1][1] ;
				ret[i][4][0] = ret[i][2][2] ;
				ret[i][4][2] = ret[i][0][0] ;
				ret[i][5][0] = ret[i][1][1] ;
				ret[i][5][1] = ret[i][0][0] ;
			}
			
			return ret ;
		}
		else
		{
			std::valarray<Matrix> ret(Matrix(3,6), gp.gaussPoints.size()) ;
			
			Vector dxi = deval(f, var[0], gp) ;
			Vector deta = deval(f, var[1], gp) ;
			Vector dzeta = deval(f, var[2], gp) ;
			Vector dteta = deval(f, var[4], gp) ;
			
			for(size_t i = 0 ; i < ret.size() ; i++)
			{
				ret[i][0][0] = dxi[i]*m[i][0][0] + deta[i]*m[i][0][1] + dzeta[i]*m[i][0][2]+ dteta[i]*m[i][0][3];
				ret[i][1][1] = dxi[i]*m[i][1][0] + deta[i]*m[i][1][1] + dzeta[i]*m[i][1][2]+ dteta[i]*m[i][0][3];
				ret[i][2][2] = dxi[i]*m[i][2][0] + deta[i]*m[i][2][1] + dzeta[i]*m[i][2][2]+ dteta[i]*m[i][0][3];
				ret[i][1][3] = ret[i][2][2] ;
				ret[i][2][3] = ret[i][1][1] ;
				ret[i][0][4] = ret[i][2][2] ;
				ret[i][2][4] = ret[i][0][0] ;
				ret[i][0][5] = ret[i][1][1] ;
				ret[i][1][5] = ret[i][0][0] ;
			}
			
			return ret ;
		}
	}
	
	return std::valarray<Matrix>(Matrix(0,0), gp.gaussPoints.size()) ; ;
}

void VirtualMachine::geval(const Function &f, const std::valarray<Matrix> & m, const std::vector<Variable> & var, const GaussPointArray & gp, bool transpose, std::valarray<Matrix> &ret)
{
	if(var.size() == 3 )
	{
		if(var[2] == TIME_VARIABLE)
		{
			if(transpose)
			{
				
				Vector dxi = deval(f, var[0], gp) ;
				Vector deta = deval(f, var[1], gp) ;
				
				for(size_t i = 0 ; i < ret.size() ; i++)
				{
					ret[i][0][0] = dxi[i]*m[i][0][0] + deta[i]*m[i][0][1] ;
					ret[i][1][1] = dxi[i]*m[i][1][0] + deta[i]*m[i][1][1] ;
					ret[i][2][0] = ret[i][1][1] ;
					ret[i][2][1] = ret[i][0][0] ;
				}
				
			}
			else
			{
				
				Vector dxi = deval(f, var[0], gp) ;
				Vector deta = deval(f, var[1], gp) ;
				
				for(size_t i = 0 ; i < ret.size() ; i++)
				{
					ret[i][0][0] = dxi[i]*m[i][0][0] + deta[i]*m[i][0][1] ;
					ret[i][1][1] = dxi[i]*m[i][1][0] + deta[i]*m[i][1][1] ;
					ret[i][0][2] = ret[i][1][1] ;
					ret[i][1][2] = ret[i][0][0] ;
				}
				ret[0].print() ;
			}
		}
		else
		{
			if(transpose)
			{
				
				Vector dxi = deval(f, var[0], gp) ;
				Vector deta = deval(f, var[1], gp) ;
				Vector dzeta = deval(f, var[2],gp) ;
				
				for(size_t i = 0 ; i < ret.size() ; i++)
				{
					ret[i][0][0] = dxi[i]*m[i][0][0] + deta[i]*m[i][0][1] + dzeta[i]*m[i][0][2] ;
					ret[i][1][1] = dxi[i]*m[i][1][0] + deta[i]*m[i][1][1] + dzeta[i]*m[i][1][2] ;
					ret[i][2][2] = dxi[i]*m[i][2][0] + deta[i]*m[i][2][1] + dzeta[i]*m[i][2][2] ;
					ret[i][3][1] = ret[i][2][2] ;
					ret[i][3][2] = ret[i][1][1] ;
					ret[i][4][0] = ret[i][2][2] ;
					ret[i][4][2] = ret[i][0][0] ;
					ret[i][5][0] = ret[i][1][1] ;
					ret[i][5][1] = ret[i][0][0] ;
				}
				
			}
			else
			{
				
				Vector dxi = deval(f, var[0],gp) ;
				Vector deta = deval(f, var[1],gp) ;
				Vector dzeta = deval(f, var[2], gp) ;
				
				for(size_t i = 0 ; i < ret.size() ; i++)
				{
					ret[i][0][0] = dxi[i]*m[i][0][0] + deta[i]*m[i][0][1] + dzeta[i]*m[i][0][2];
					ret[i][1][1] = dxi[i]*m[i][1][0] + deta[i]*m[i][1][1] + dzeta[i]*m[i][1][2];
					ret[i][2][2] = dxi[i]*m[i][2][0] + deta[i]*m[i][2][1] + dzeta[i]*m[i][2][2];
					ret[i][1][3] = ret[i][2][2] ;
					ret[i][2][3] = ret[i][1][1] ;
					ret[i][0][4] = ret[i][2][2] ;
					ret[i][2][4] = ret[i][0][0] ;
					ret[i][0][5] = ret[i][1][1] ;
					ret[i][1][5] = ret[i][0][0] ;
				}
				
			}
		}
	}
	else if (var.size() == 2)
	{
		if(transpose)
		{
			
			
			Vector dxi = deval(f, var[0], gp) ;
			Vector deta = deval(f, var[1], gp) ;
			
			for(size_t i = 0 ; i < ret.size() ; i++)
			{
				ret[i][0][0] = dxi[i]*m[i][0][0] + deta[i]*m[i][0][1] ;
				ret[i][1][1] = dxi[i]*m[i][1][0] + deta[i]*m[i][1][1] ;
				ret[i][2][0] = ret[i][1][1] ;
				ret[i][2][1] = ret[i][0][0] ;
			}
			
		}
		else
		{
			
			Vector dxi = deval(f, var[0], gp) ;
			Vector deta = deval(f, var[1], gp) ;
			
			for(size_t i = 0 ; i < ret.size() ; i++)
			{
				ret[i][0][0] = dxi[i]*m[i][0][0] + deta[i]*m[i][0][1] ;
				ret[i][1][1] = dxi[i]*m[i][1][0] + deta[i]*m[i][1][1] ;
				ret[i][0][2] = ret[i][1][1] ;
				ret[i][1][2] = ret[i][0][0] ;
			}
			
		}
	}
	else if(var.size() == 4 )
	{
		if(transpose)
		{
			
			Vector dxi = deval(f, var[0], gp) ;
			Vector deta = deval(f, var[1], gp) ;
			Vector dzeta = deval(f, var[2], gp) ;
			Vector dteta = deval(f, var[4], gp) ;
			
			for(size_t i = 0 ; i < ret.size() ; i++)
			{
				ret[i][0][0] = dxi[i]*m[i][0][0] + deta[i]*m[i][0][1] + dzeta[i]*m[i][0][2] + dteta[i]*m[i][0][3];
				ret[i][1][1] = dxi[i]*m[i][1][0] + deta[i]*m[i][1][1] + dzeta[i]*m[i][1][2] + dteta[i]*m[i][1][3];
				ret[i][2][2] = dxi[i]*m[i][2][0] + deta[i]*m[i][2][1] + dzeta[i]*m[i][2][2] + dteta[i]*m[i][2][3];
				ret[i][3][1] = ret[i][2][2] ;
				ret[i][3][2] = ret[i][1][1] ;
				ret[i][4][0] = ret[i][2][2] ;
				ret[i][4][2] = ret[i][0][0] ;
				ret[i][5][0] = ret[i][1][1] ;
				ret[i][5][1] = ret[i][0][0] ;
			}
			
		}
		else
		{
			
			Vector dxi = deval(f, var[0], gp) ;
			Vector deta = deval(f, var[1], gp) ;
			Vector dzeta = deval(f, var[2], gp) ;
			Vector dteta = deval(f, var[4], gp) ;
			
			for(size_t i = 0 ; i < ret.size() ; i++)
			{
				ret[i][0][0] = dxi[i]*m[i][0][0] + deta[i]*m[i][0][1] + dzeta[i]*m[i][0][2]+ dteta[i]*m[i][0][3];
				ret[i][1][1] = dxi[i]*m[i][1][0] + deta[i]*m[i][1][1] + dzeta[i]*m[i][1][2]+ dteta[i]*m[i][0][3];
				ret[i][2][2] = dxi[i]*m[i][2][0] + deta[i]*m[i][2][1] + dzeta[i]*m[i][2][2]+ dteta[i]*m[i][0][3];
				ret[i][1][3] = ret[i][2][2] ;
				ret[i][2][3] = ret[i][1][1] ;
				ret[i][0][4] = ret[i][2][2] ;
				ret[i][2][4] = ret[i][0][0] ;
				ret[i][0][5] = ret[i][1][1] ;
				ret[i][1][5] = ret[i][0][0] ;
			}
			
		}
	}
	
}

std::valarray<Matrix> VirtualMachine::gdeval(const Function &f, const std::valarray<Matrix> & m, const std::vector<Variable> & var, const GaussPointArray & gp, bool transpose )
{
	if(var.size() == 3 )
	{
		if(var[2] == TIME_VARIABLE)
		{
			if(transpose)
			{
				std::valarray<Matrix> ret(Matrix(3,2), gp.gaussPoints.size()) ;
				
				Vector dxi = ddeval(f, var[0],TIME_VARIABLE, gp) ;
				Vector deta = ddeval(f, var[1],TIME_VARIABLE, gp) ;

				for(size_t i = 0 ; i < ret.size() ; i++)
				{
					ret[i][0][0] = dxi[i]*m[i][0][0] + deta[i]*m[i][0][1] ;
					ret[i][1][1] = dxi[i]*m[i][1][0] + deta[i]*m[i][1][1] ;
					ret[i][2][0] = ret[i][1][1] ;
					ret[i][2][1] = ret[i][0][0] ;
				}
				
				return ret ;
			}
			else
			{
				std::valarray<Matrix> ret(Matrix(2,3), gp.gaussPoints.size()) ;
				
				Vector dxi = ddeval(f, var[0],TIME_VARIABLE, gp) ;
				Vector deta = ddeval(f, var[1],TIME_VARIABLE, gp) ;
				
				for(size_t i = 0 ; i < ret.size() ; i++)
				{
					ret[i][0][0] = dxi[i]*m[i][0][0] + deta[i]*m[i][0][1] ;
					ret[i][1][1] = dxi[i]*m[i][1][0] + deta[i]*m[i][1][1] ;
					ret[i][0][2] = ret[i][1][1] ;
					ret[i][1][2] = ret[i][0][0] ;
				}
				
				return ret ;
			}
		}
		else
		{
			if(transpose)
			{
				std::valarray<Matrix> ret(Matrix(6,3), gp.gaussPoints.size()) ;
				
				Vector dxi = ddeval(f, var[0],TIME_VARIABLE, gp) ;
				Vector deta = ddeval(f, var[1],TIME_VARIABLE, gp) ;
				Vector dzeta = ddeval(f, var[2],TIME_VARIABLE,gp) ;
				
				for(size_t i = 0 ; i < ret.size() ; i++)
				{
					ret[i][0][0] = dxi[i]*m[i][0][0] + deta[i]*m[i][0][1] + dzeta[i]*m[i][0][2] ;
					ret[i][1][1] = dxi[i]*m[i][1][0] + deta[i]*m[i][1][1] + dzeta[i]*m[i][1][2] ;
					ret[i][2][2] = dxi[i]*m[i][2][0] + deta[i]*m[i][2][1] + dzeta[i]*m[i][2][2] ;
					ret[i][3][1] = ret[i][2][2] ;
					ret[i][3][2] = ret[i][1][1] ;
					ret[i][4][0] = ret[i][2][2] ;
					ret[i][4][2] = ret[i][0][0] ;
					ret[i][5][0] = ret[i][1][1] ;
					ret[i][5][1] = ret[i][0][0] ;
				}
				
				return ret ;
			}
			else
			{
				std::valarray<Matrix> ret(Matrix(3,6), gp.gaussPoints.size()) ;
				
				Vector dxi = ddeval(f, var[0],TIME_VARIABLE,gp) ;
				Vector deta = ddeval(f, var[1],TIME_VARIABLE,gp) ;
				Vector dzeta = ddeval(f, var[2],TIME_VARIABLE, gp) ;
				
				for(size_t i = 0 ; i < ret.size() ; i++)
				{
					ret[i][0][0] = dxi[i]*m[i][0][0] + deta[i]*m[i][0][1] + dzeta[i]*m[i][0][2];
					ret[i][1][1] = dxi[i]*m[i][1][0] + deta[i]*m[i][1][1] + dzeta[i]*m[i][1][2];
					ret[i][2][2] = dxi[i]*m[i][2][0] + deta[i]*m[i][2][1] + dzeta[i]*m[i][2][2];
					ret[i][1][3] = ret[i][2][2] ;
					ret[i][2][3] = ret[i][1][1] ;
					ret[i][0][4] = ret[i][2][2] ;
					ret[i][2][4] = ret[i][0][0] ;
					ret[i][0][5] = ret[i][1][1] ;
					ret[i][1][5] = ret[i][0][0] ;
				}
				
				return ret ;
			}
		}
	}
	else if (var.size() == 2)
	{
		if(transpose)
		{
			
			std::valarray<Matrix> ret(Matrix(3,2), gp.gaussPoints.size()) ;
			
			Vector dxi = ddeval(f, var[0],TIME_VARIABLE, gp) ;
			Vector deta = ddeval(f, var[1],TIME_VARIABLE, gp) ;
			
			for(size_t i = 0 ; i < ret.size() ; i++)
			{
				ret[i][0][0] = dxi[i]*m[i][0][0] + deta[i]*m[i][0][1] ;
				ret[i][1][1] = dxi[i]*m[i][1][0] + deta[i]*m[i][1][1] ;
				ret[i][2][0] = ret[i][1][1] ;
				ret[i][2][1] = ret[i][0][0] ;
			}
			
			return ret ;
		}
		else
		{
			std::valarray<Matrix> ret(Matrix(2,3), gp.gaussPoints.size()) ;
			
			Vector dxi = ddeval(f, var[0],TIME_VARIABLE, gp) ;
			Vector deta = ddeval(f, var[1],TIME_VARIABLE, gp) ;
			
			for(size_t i = 0 ; i < ret.size() ; i++)
			{
				ret[i][0][0] = dxi[i]*m[i][0][0] + deta[i]*m[i][0][1] ;
				ret[i][1][1] = dxi[i]*m[i][1][0] + deta[i]*m[i][1][1] ;
				ret[i][0][2] = ret[i][1][1] ;
				ret[i][1][2] = ret[i][0][0] ;
			}
			
			return ret ;
		}
	}
	else if(var.size() == 4 )
	{
		if(transpose)
		{
			std::valarray<Matrix> ret(Matrix(6,3), gp.gaussPoints.size()) ;
			
			Vector dxi = ddeval(f, var[0],TIME_VARIABLE, gp) ;
			Vector deta = ddeval(f, var[1],TIME_VARIABLE, gp) ;
			Vector dzeta = ddeval(f, var[2],TIME_VARIABLE, gp) ;
			Vector dteta = ddeval(f, var[3],TIME_VARIABLE, gp) ;
			
			for(size_t i = 0 ; i < ret.size() ; i++)
			{
				ret[i][0][0] = dxi[i]*m[i][0][0] + deta[i]*m[i][0][1] + dzeta[i]*m[i][0][2] + dteta[i]*m[i][0][3];
				ret[i][1][1] = dxi[i]*m[i][1][0] + deta[i]*m[i][1][1] + dzeta[i]*m[i][1][2] + dteta[i]*m[i][1][3];
				ret[i][2][2] = dxi[i]*m[i][2][0] + deta[i]*m[i][2][1] + dzeta[i]*m[i][2][2] + dteta[i]*m[i][2][3];
				ret[i][3][1] = ret[i][2][2] ;
				ret[i][3][2] = ret[i][1][1] ;
				ret[i][4][0] = ret[i][2][2] ;
				ret[i][4][2] = ret[i][0][0] ;
				ret[i][5][0] = ret[i][1][1] ;
				ret[i][5][1] = ret[i][0][0] ;
			}
			
			return ret ;
		}
		else
		{
			std::valarray<Matrix> ret(Matrix(3,6), gp.gaussPoints.size()) ;
			
			Vector dxi = ddeval(f, var[0],TIME_VARIABLE, gp) ;
			Vector deta = ddeval(f, var[1],TIME_VARIABLE, gp) ;
			Vector dzeta = ddeval(f, var[2],TIME_VARIABLE, gp) ;
			Vector dteta = ddeval(f, var[3],TIME_VARIABLE, gp) ;
			
			for(size_t i = 0 ; i < ret.size() ; i++)
			{
				ret[i][0][0] = dxi[i]*m[i][0][0] + deta[i]*m[i][0][1] + dzeta[i]*m[i][0][2]+ dteta[i]*m[i][0][3];
				ret[i][1][1] = dxi[i]*m[i][1][0] + deta[i]*m[i][1][1] + dzeta[i]*m[i][1][2]+ dteta[i]*m[i][0][3];
				ret[i][2][2] = dxi[i]*m[i][2][0] + deta[i]*m[i][2][1] + dzeta[i]*m[i][2][2]+ dteta[i]*m[i][0][3];
				ret[i][1][3] = ret[i][2][2] ;
				ret[i][2][3] = ret[i][1][1] ;
				ret[i][0][4] = ret[i][2][2] ;
				ret[i][2][4] = ret[i][0][0] ;
				ret[i][0][5] = ret[i][1][1] ;
				ret[i][1][5] = ret[i][0][0] ;
			}
			
			return ret ;
		}
	}
	
	return std::valarray<Matrix>(Matrix(0,0), gp.gaussPoints.size()) ; ;
}


Matrix VirtualMachine::geval(const Function &f, const Matrix & m, const std::vector<Variable> & var, const double x, const double y, const double  z, const double  t, bool transpose )
{
	if(var.size() == 3 )
	{
		if(var[2] == TIME_VARIABLE)
		{
			if(transpose)
			{
				double dxi = deval(f, var[0], x,y,z,t) ;
				double deta = deval(f, var[1], x,y,z,t) ;

				Matrix ret(3,2) ;
				ret[0][0] = dxi*m[0][0] + deta*m[0][1] ;
				ret[1][1] = dxi*m[1][0] + deta*m[1][1] ;
				ret[2][0] = ret[1][1] ;
				ret[2][1] = ret[0][0] ;
				
				return ret ;
			}
			else
			{
				Matrix ret(2,3) ;
				double dxi = deval(f, var[0], x,y,z) ;
				double deta = deval(f, var[1], x,y,z) ;

				ret[0][0] = dxi*m[0][0] + deta*m[0][1] ;
				ret[1][1] = dxi*m[1][0] + deta*m[1][1] ;
				ret[0][2] = ret[1][1] ;
				ret[1][2] = ret[0][0] ;
				
				return ret ;
			}
		}
		else
		{
			if(transpose)
			{
				double dxi = deval(f, var[0], x,y,z,t) ;
				double deta = deval(f, var[1], x,y,z,t) ;
				double dzeta = deval(f, var[2], x,y,z,t) ;
				Matrix ret(6,3) ;
				ret[0][0] = dxi*m[0][0] + deta*m[0][1] + dzeta*m[0][2] ;
				ret[1][1] = dxi*m[1][0] + deta*m[1][1] + dzeta*m[1][2] ;
				ret[2][2] = dxi*m[2][0] + deta*m[2][1] + dzeta*m[2][2] ;
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
				double dxi = deval(f, var[0], x,y,z,t) ;
				double deta = deval(f, var[1], x,y,z,t) ;
				double dzeta = deval(f, var[2], x,y,z,t) ;
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
	}
	else if (var.size() == 2)
	{
		if(transpose)
		{
			
			Matrix ret(3,2) ;
			double dxi = deval(f, var[0], x,y,z) ;
			double deta = deval(f, var[1], x,y,z) ;
			ret[0][0] = dxi*m[0][0] + deta*m[0][1] ;
			ret[1][1] = dxi*m[1][0] + deta*m[1][1] ;
			ret[2][0] = ret[1][1] ;
			ret[2][1] = ret[0][0] ;
			
			return ret ;
		}
		else
		{
			Matrix ret(2,3) ;
			double dxi = deval(f, var[0], x,y,z) ;
			double deta = deval(f, var[1], x,y,z) ;
			ret[0][0] = dxi*m[0][0] + deta*m[0][1] ;
			ret[1][1] = dxi*m[1][0] + deta*m[1][1] ;
			ret[0][2] = ret[1][1] ;
			ret[1][2] = ret[0][0] ;
			
			return ret ;
		}
	}
	else if(var.size() == 4 )
	{
		if(transpose)
		{
			double dxi = deval(f, var[0], x,y,z,t) ;
			double deta = deval(f, var[1], x,y,z,t) ;
			double dzeta = deval(f, var[2], x,y,z,t) ;
			double dteta = deval(f, var[4], x,y,z,t) ;
			Matrix ret(6,3) ;
			ret[0][0] = dxi*m[0][0] + deta*m[0][1] + dzeta*m[0][2] + dteta*m[0][3];
			ret[1][1] = dxi*m[1][0] + deta*m[1][1] + dzeta*m[1][2] + dteta*m[1][3];
			ret[2][2] = dxi*m[2][0] + deta*m[2][1] + dzeta*m[2][2] + dteta*m[2][3];
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
			double dxi = deval(f, var[0], x,y,z) ;
			double deta = deval(f, var[1], x,y,z) ;
			double dzeta = deval(f, var[2], x,y,z) ;
			double dteta = deval(f, var[4], x,y,z,t) ;
			ret[0][0] = dxi*m[0][0] + deta*m[0][1] + dzeta*m[0][2]+ dteta*m[0][3];
			ret[1][1] = dxi*m[1][0] + deta*m[1][1] + dzeta*m[1][2]+ dteta*m[0][3];
			ret[2][2] = dxi*m[2][0] + deta*m[2][1] + dzeta*m[2][2]+ dteta*m[0][3];
			ret[1][3] = ret[2][2] ;
			ret[2][3] = ret[1][1] ;
			ret[0][4] = ret[2][2] ;
			ret[2][4] = ret[0][0] ;
			ret[0][5] = ret[1][1] ;
			ret[1][5] = ret[0][0] ;
			
			return ret ;
		}
	}
	return Matrix(0,0) ;
}

void VirtualMachine::geval(const Function &f, const Matrix & m, const std::vector<Variable> & var, const double x, const double y, const double  z, const double  t, bool transpose, Matrix & ret )
{
	if(var.size() == 3 )
	{
		if(var[2] == TIME_VARIABLE)
		{
			if(transpose)
			{
				double dxi = deval(f, var[0], x,y,z,t) ;
				double deta = deval(f, var[1], x,y,z,t) ;
				if(ret.isNull() || ret.numRows() !=3 ||ret.numCols() != 2)
					ret.resize(3,2) ;
				ret[0][0] = dxi*m[0][0] + deta*m[0][1] ;
				ret[1][1] = dxi*m[1][0] + deta*m[1][1] ;
				ret[2][0] = ret[1][1] ;
				ret[2][1] = ret[0][0] ;
				
			}
			else
			{
				double dxi = deval(f, var[0], x,y,z,t) ;
				double deta = deval(f, var[1], x,y,z,t) ;
				if(ret.isNull() || ret.numRows() != 2 ||ret.numCols() != 3)
					ret.resize(2,3) ;
				ret[0][0] = dxi*m[0][0] + deta*m[0][1] ;
				ret[1][1] = dxi*m[1][0] + deta*m[1][1] ;
				ret[0][2] = ret[1][1] ;
				ret[1][2] = ret[0][0] ;
			}
		}
		else
		{
			if(transpose)
			{
				double dxi = deval(f, var[0], x,y,z,t) ;
				double deta = deval(f, var[1], x,y,z,t) ;
				double dzeta = deval(f, var[2], x,y,z,t) ;
				
				if(ret.isNull()|| ret.numRows() != 6 ||ret.numCols() != 3)
					ret.resize(6,3) ;
				ret[0][0] = dxi*m[0][0] + deta*m[0][1] + dzeta*m[0][2] ;
				ret[1][1] = dxi*m[1][0] + deta*m[1][1] + dzeta*m[1][2] ;
				ret[2][2] = dxi*m[2][0] + deta*m[2][1] + dzeta*m[2][2] ;
				ret[3][1] = ret[2][2] ;
				ret[3][2] = ret[1][1] ;
				ret[4][0] = ret[2][2] ;
				ret[4][2] = ret[0][0] ;
				ret[5][0] = ret[1][1] ;
				ret[5][1] = ret[0][0] ;
				
			}
			else
			{
				double dxi = deval(f, var[0], x,y,z,t) ;
				double deta = deval(f, var[1], x,y,z,t) ;
				double dzeta = deval(f, var[2], x,y,z,t) ;
				if(ret.isNull()|| ret.numRows() != 3 ||ret.numCols() != 6)
					ret.resize(3,6) ;
				ret[0][0] = dxi*m[0][0] + deta*m[0][1] + dzeta*m[0][2];
				ret[1][1] = dxi*m[1][0] + deta*m[1][1] + dzeta*m[1][2];
				ret[2][2] = dxi*m[2][0] + deta*m[2][1] + dzeta*m[2][2];
				ret[1][3] = ret[2][2] ;
				ret[2][3] = ret[1][1] ;
				ret[0][4] = ret[2][2] ;
				ret[2][4] = ret[0][0] ;
				ret[0][5] = ret[1][1] ;
				ret[1][5] = ret[0][0] ;
				
			}
		}
	}
	else if (var.size() == 2)
	{
		if(transpose)
		{
			if(ret.isNull()|| ret.numRows() != 3 ||ret.numCols() != 2)
				ret.resize(3,2) ;
			double dxi = deval(f, var[0], x,y,z) ;
			double deta = deval(f, var[1], x,y,z) ;
			ret[0][0] = dxi*m[0][0] + deta*m[0][1] ;
			ret[1][1] = dxi*m[1][0] + deta*m[1][1] ;
			ret[2][0] = ret[1][1] ;
			ret[2][1] = ret[0][0] ;
		
			
		}
		else
		{
			if(ret.isNull()|| ret.numRows() != 2 ||ret.numCols() != 3)
				ret.resize(2,3) ;
			double dxi = deval(f, var[0], x,y,z) ;
			double deta = deval(f, var[1], x,y,z) ;
			ret[0][0] = dxi*m[0][0] + deta*m[0][1] ;
			ret[1][1] = dxi*m[1][0] + deta*m[1][1] ;
			ret[0][2] = ret[1][1] ;
			ret[1][2] = ret[0][0] ;	
			
		}
	}
	else if(var.size() == 4 )
	{
		if(transpose)
		{
			if(ret.isNull()|| ret.numRows() != 6 ||ret.numCols() != 3)
				ret.resize(6,3) ;
			double dxi = deval(f, var[0], x,y,z,t) ;
			double deta = deval(f, var[1], x,y,z,t) ;
			double dzeta = deval(f, var[2], x,y,z,t) ;
			double dteta = deval(f, var[4], x,y,z,t) ;
			ret[0][0] = dxi*m[0][0] + deta*m[0][1] + dzeta*m[0][2] + dteta*m[0][3];
			ret[1][1] = dxi*m[1][0] + deta*m[1][1] + dzeta*m[1][2] + dteta*m[1][3];
			ret[2][2] = dxi*m[2][0] + deta*m[2][1] + dzeta*m[2][2] + dteta*m[2][3];
			ret[3][1] = ret[2][2] ;
			ret[3][2] = ret[1][1] ;
			ret[4][0] = ret[2][2] ;
			ret[4][2] = ret[0][0] ;
			ret[5][0] = ret[1][1] ;
			ret[5][1] = ret[0][0] ;
			
		}
		else
		{
			if(ret.isNull()|| ret.numRows() != 3 ||ret.numCols() != 6)
				ret.resize(3,6) ;
			double dxi = deval(f, var[0], x,y,z) ;
			double deta = deval(f, var[1], x,y,z) ;
			double dzeta = deval(f, var[2], x,y,z) ;
			double dteta = deval(f, var[4], x,y,z,t) ;
			ret[0][0] = dxi*m[0][0] + deta*m[0][1] + dzeta*m[0][2]+ dteta*m[0][3];
			ret[1][1] = dxi*m[1][0] + deta*m[1][1] + dzeta*m[1][2]+ dteta*m[0][3];
			ret[2][2] = dxi*m[2][0] + deta*m[2][1] + dzeta*m[2][2]+ dteta*m[0][3];
			ret[1][3] = ret[2][2] ;
			ret[2][3] = ret[1][1] ;
			ret[0][4] = ret[2][2] ;
			ret[2][4] = ret[0][0] ;
			ret[0][5] = ret[1][1] ;
			ret[1][5] = ret[0][0] ;
			
		}
	}
}


void VirtualMachine::geval(const Function &f, const Matrix & m, const std::vector<Variable> & var, const Point & p, bool transpose, Matrix & ret )
{
	geval(f, m, var, p.x, p.y, p.z, p.t, transpose, ret ) ;
}

Matrix VirtualMachine::gdeval(const Function &f, const Matrix & m, const std::vector<Variable> & var, const double x, const double y, const double  z, const double  t, bool transpose )
{
	if(var.size() == 3 )
	{
		if(var[2] == TIME_VARIABLE)
		{
			if(transpose)
			{
				double dxi = ddeval(f, var[0],TIME_VARIABLE, x,y,z,t) ;
				double deta = ddeval(f, var[1],TIME_VARIABLE, x,y,z,t) ;
				
				Matrix ret(3,2) ;
				ret[0][0] = dxi*m[0][0] + deta*m[0][1] ;
				ret[1][1] = dxi*m[1][0] + deta*m[1][1] ;
				ret[2][0] = ret[1][1] ;
				ret[2][1] = ret[0][0] ;

				return ret ;
			}
			else
			{
				Matrix ret(2,3) ;
				double dxi = ddeval(f, var[0],TIME_VARIABLE, x,y,z) ;
				double deta = ddeval(f, var[1],TIME_VARIABLE, x,y,z) ;
				
				ret[0][0] = dxi*m[0][0] + deta*m[0][1] ;
				ret[1][1] = dxi*m[1][0] + deta*m[1][1] ;
				ret[2][0] = ret[1][1] ;
				ret[2][1] = ret[0][0] ;
				return ret ;
			}
		}
		else
		{
			if(transpose)
			{
				double dxi = ddeval(f, var[0],TIME_VARIABLE, x,y,z,t) ;
				double deta = ddeval(f, var[1],TIME_VARIABLE, x,y,z,t) ;
				double dzeta = ddeval(f, var[2],TIME_VARIABLE, x,y,z,t) ;
				Matrix ret(6,3) ;
				ret[0][0] = dxi*m[0][0] + deta*m[0][1] + dzeta*m[0][2] ;
				ret[1][1] = dxi*m[1][0] + deta*m[1][1] + dzeta*m[1][2] ;
				ret[2][2] = dxi*m[2][0] + deta*m[2][1] + dzeta*m[2][2] ;
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
				double dxi = ddeval(f, var[0],TIME_VARIABLE, x,y,z,t) ;
				double deta = ddeval(f, var[1],TIME_VARIABLE, x,y,z,t) ;
				double dzeta = ddeval(f, var[2],TIME_VARIABLE, x,y,z,t) ;
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
	}
	else if (var.size() == 2)
	{
		if(transpose)
		{
			
			Matrix ret(3,2) ;
			double dxi = ddeval(f, var[0],TIME_VARIABLE, x,y,z) ;
			double deta = ddeval(f, var[1],TIME_VARIABLE, x,y,z) ;
			ret[0][0] = dxi*m[0][0] + deta*m[0][1] ;
			ret[1][1] = dxi*m[1][0] + deta*m[1][1] ;
			ret[2][0] = ret[1][1] ;
			ret[2][1] = ret[0][0] ;
			
			return ret ;
		}
		else
		{
			Matrix ret(2,3) ;
			double dxi = ddeval(f, var[0],TIME_VARIABLE, x,y,z) ;
			double deta = ddeval(f, var[1],TIME_VARIABLE, x,y,z) ;
			ret[0][0] = dxi*m[0][0] + deta*m[0][1] ;
			ret[1][1] = dxi*m[1][0] + deta*m[1][1] ;
			ret[0][2] = ret[1][1] ;
			ret[1][2] = ret[0][0] ;
			
			return ret ;
		}
	}
	else if(var.size() == 4 )
	{
		if(transpose)
		{
			double dxi = ddeval(f, var[0],TIME_VARIABLE, x,y,z,t) ;
			double deta = ddeval(f, var[1],TIME_VARIABLE, x,y,z,t) ;
			double dzeta = ddeval(f, var[2],TIME_VARIABLE, x,y,z,t) ;
			double dteta = ddeval(f, var[4],TIME_VARIABLE, x,y,z,t) ;
			Matrix ret(6,3) ;
			ret[0][0] = dxi*m[0][0] + deta*m[0][1] + dzeta*m[0][2] + dteta*m[0][3];
			ret[1][1] = dxi*m[1][0] + deta*m[1][1] + dzeta*m[1][2] + dteta*m[1][3];
			ret[2][2] = dxi*m[2][0] + deta*m[2][1] + dzeta*m[2][2] + dteta*m[2][3];
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
			double dxi = ddeval(f, var[0],TIME_VARIABLE, x,y,z) ;
			double deta = ddeval(f, var[1],TIME_VARIABLE, x,y,z) ;
			double dzeta = ddeval(f, var[2],TIME_VARIABLE, x,y,z) ;
			double dteta = ddeval(f, var[4],TIME_VARIABLE, x,y,z,t) ;
			ret[0][0] = dxi*m[0][0] + deta*m[0][1] + dzeta*m[0][2]+ dteta*m[0][3];
			ret[1][1] = dxi*m[1][0] + deta*m[1][1] + dzeta*m[1][2]+ dteta*m[0][3];
			ret[2][2] = dxi*m[2][0] + deta*m[2][1] + dzeta*m[2][2]+ dteta*m[0][3];
			ret[1][3] = ret[2][2] ;
			ret[2][3] = ret[1][1] ;
			ret[0][4] = ret[2][2] ;
			ret[2][4] = ret[0][0] ;
			ret[0][5] = ret[1][1] ;
			ret[1][5] = ret[0][0] ;
			
			return ret ;
		}
	}
	return Matrix(0,0) ;
}

Matrix VirtualMachine::gveval(const Function &f, IntegrableEntity *e, const std::vector<Variable> & var, const double x, const double y, const double z, const double t, bool transpose)
{
	Matrix Jinv ;e->getInverseJacobianMatrix(Point(x,y,z,t), Jinv) ;
	return gveval(f, Jinv, var, x, y, z,t , transpose) ;
}

Matrix VirtualMachine::gveval(const VectorGradient &f, IntegrableEntity *e, const std::vector<Variable> & var, const double x, const double  y, const double z, const double t)
{
	return gveval(f.f, e,var, x, y, z,t, f.transpose) ;
}

Matrix VirtualMachine::gveval(const Function &f, const Matrix & m, const std::vector<Variable> & var, const double x, const double y, const double  z, const double  t, bool transpose )
{
	if(var.size() == 3 )
	{
		if(var[2] == ZETA)
		{
			if(transpose)
			{
				Matrix ret(3,1) ;
				double fx = deval(f, var[0], x,y,z,t) ;
				double fy = deval(f, var[1], x,y,z,t) ;
				double fz = deval(f, var[2], x,y,z,t) ;
				ret[0][0] = fx*m[0][0] + fy*m[0][1] + fz*m[0][2];
				ret[1][0] = fx*m[1][0] + fy*m[1][1] + fz*m[1][2];
				ret[2][0] = fx*m[2][0] + fy*m[2][1] + fz*m[2][2];
	
				
				return ret ;
			}
			else
			{
				Matrix ret(1,3) ;
				double fx = deval(f, var[0], x,y,z,t) ;
				double fy = deval(f, var[1], x,y,z,t) ;
				double fz = deval(f, var[2], x,y,z,t) ;
				ret[0][0] = fx*m[0][0] + fy*m[0][1] + fz*m[0][2];
				ret[0][1] = fx*m[1][0] + fy*m[1][1] + fz*m[1][2];
				ret[0][2] = fx*m[2][0] + fy*m[2][1] + fz*m[2][2];
				
		
				return ret ;
			}
		}
		else
		{
			if(transpose)
			{
				Matrix ret(2,1) ;
				double fx = deval(f, var[0], x,y,z,t) ;
				double fy = deval(f, var[1], x,y,z,t) ;
				ret[0][0] = fx*m[0][0] + fy*m[0][1] ;
				ret[1][0] = fx*m[1][0] + fy*m[1][1] ;

				
				
				return ret ;
			}
			else
			{
				Matrix ret(1,2) ;
				double fx = deval(f, var[0], x,y,z,t) ;
				double fy = deval(f, var[1], x,y,z,t) ;
				ret[0][0] = fx*m[0][0] + fy*m[0][1] ;
				ret[0][1] = fx*m[1][0] + fy*m[1][1] ;
				
				return ret ;
			}
		}
	}
	else if(var.size() == 2)
	{
		if(transpose)
		{
			Matrix ret(2,1) ;
			double fx = deval(f, var[0], x,y,z,t) ;
			double fy = deval(f, var[1], x,y,z,t) ;
			ret[0][0] = fx*m[0][0] + fy*m[0][1] ;
			ret[1][0] = fx*m[1][0] + fy*m[1][1] ;

			
			return ret ;
		}
		else
		{
			Matrix ret(1,2) ;
			double fx = deval(f, var[0], x,y,z,t) ;
			double fy = deval(f, var[1], x,y,z,t) ;
			ret[0][0] = fx*m[0][0] + fy*m[0][1] ;
			ret[0][1] = fx*m[1][0] + fy*m[1][1] ;
			
			return ret ;
		}
	}
	else if(var.size() == 4)
	{
		if(transpose)
		{
			Matrix ret(3,1) ;
			double fx = deval(f, var[0], x,y,z,t) ;
			double fy = deval(f, var[1], x,y,z,t) ;
			double fz = deval(f, var[2], x,y,z,t) ;
			ret[0][0] = fx*m[0][0] + fy*m[0][1] + fz*m[0][2];
			ret[1][0] = fx*m[1][0] + fy*m[1][1] + fz*m[1][2];
			ret[2][0] = fx*m[2][0] + fy*m[2][1] + fz*m[2][2];

			return ret ;
		}
		else
		{
			Matrix ret(1,3) ;
			double fx = deval(f, var[0], x,y,z,t) ;
			double fy = deval(f, var[1], x,y,z,t) ;
			double fz = deval(f, var[2], x,y,z,t) ;
			ret[0][0] = fx*m[0][0] + fy*m[0][1] + fz*m[0][2];
			ret[0][1] = fx*m[1][0] + fy*m[1][1] + fz*m[1][2];
			ret[0][2] = fx*m[2][0] + fy*m[2][1] + fz*m[2][2];

			return ret ;
		}
	}
	
	return Matrix(1,1) ;
}

Matrix VirtualMachine::gveval(const VectorGradient &f, const Matrix & m, const std::vector<Variable> & var, const double x, const double y , const double z, const double t)
{
	return gveval(f.f, m,var, x, y, z, f.transpose) ;
}


Matrix VirtualMachine::ieval(const GtM &f, IntegrableEntity *e, const std::vector<Variable> & var)
{

	GaussPointArray gp = e->getGaussPoints() ;
	Matrix B (geval(f.first.f, e,var, gp.gaussPoints[0].first.x, gp.gaussPoints[0].first.y, gp.gaussPoints[0].first.z, gp.gaussPoints[0].first.t, f.first.transpose)*gp.gaussPoints[0].second) ;
	for(size_t i = 1 ; i < gp.gaussPoints.size() ; i++)
	{
		B += geval(f.first.f, e,var, gp.gaussPoints[i].first.x, gp.gaussPoints[i].first.y, gp.gaussPoints[i].first.z, gp.gaussPoints[i].first.t, f.first.transpose)*gp.gaussPoints[i].second ;
	}
	
	return B*f.second ;
}

Matrix VirtualMachine::ieval(const VGtM &f, IntegrableEntity *e, const std::vector<Variable> & var)
{

	GaussPointArray gp = e->getGaussPoints() ;
	Matrix B(gveval(f.first.f, e,var, gp.gaussPoints[0].first.x, gp.gaussPoints[0].first.y, gp.gaussPoints[0].first.z, gp.gaussPoints[0].first.t, f.first.transpose)*gp.gaussPoints[0].second) ;
	for(size_t i = 1 ; i < gp.gaussPoints.size() ; i++)
	{
		B += geval(f.first.f, e,var, gp.gaussPoints[i].first.x, gp.gaussPoints[i].first.y, gp.gaussPoints[i].first.z, gp.gaussPoints[i].first.t, f.first.transpose)*gp.gaussPoints[i].second ;
	}
	
	return B*f.second ;
}



Matrix VirtualMachine::ieval(const GtMtG &f, IntegrableEntity *e, const std::vector<Variable> & var)
{

	GaussPointArray gp = e->getGaussPoints() ;
	std::valarray<Matrix> Jinv( Matrix(), gp.gaussPoints.size()) ;
	for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
	{
		e->getInverseJacobianMatrix(gp.gaussPoints[i].first, Jinv[i]) ;
	}
		
	return ieval(f, gp, Jinv,var) ;

}

double VirtualMachine::ieval(const VGtMtVG &f, IntegrableEntity *e, const std::vector<Variable> & var)
{

	GaussPointArray gp = e->getGaussPoints() ;
	std::valarray<Matrix> Jinv(Matrix(), gp.gaussPoints.size()) ;
	for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
	{
		e->getInverseJacobianMatrix(gp.gaussPoints[i].first, Jinv[i]) ;
	}
		
	return ieval(f, gp, Jinv,var) ;

}

Matrix VirtualMachine::ieval(const GtMtG &f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & var)
{

	geval(f.third.f, Jinv[0], var, gp.gaussPoints[0].first, f.third.transpose,B_) ;
	geval(f.first.f, Jinv[0], var, gp.gaussPoints[0].first, f.first.transpose,B) ;
	
	Matrix ret(B*f.second*B_) ;
	ret *= gp.gaussPoints[0].second ;
	
	for(size_t i = 1 ; i  < gp.gaussPoints.size() ; i++)
	{
		geval(f.third.f, Jinv[i], var, gp.gaussPoints[i].first, f.third.transpose, B_) ;
		geval(f.first.f, Jinv[i], var, gp.gaussPoints[i].first, f.first.transpose, B) ;
		matrix_matrix_matrix_multiply_and_add(B, f.second, B_, gp.gaussPoints[i].second, ret) ;
	}

	return ret ;
}

void VirtualMachine::ieval(const GtMtG &f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & var, Matrix & ret)
{
	geval(f.third.f, Jinv[0], var, gp.gaussPoints[0].first, f.third.transpose,B_) ;
	geval(f.first.f, Jinv[0], var, gp.gaussPoints[0].first, f.first.transpose,B) ;
	
	ret = B*f.second*B_ ;
	ret *= gp.gaussPoints[0].second ;
	
	for(size_t i = 1 ; i  < gp.gaussPoints.size() ; i++)
	{
		geval(f.third.f, Jinv[i], var, gp.gaussPoints[i].first, f.third.transpose, B_) ;
		geval(f.first.f, Jinv[i], var, gp.gaussPoints[i].first, f.first.transpose, B) ;
		matrix_matrix_matrix_multiply_and_add(B, f.second, B_, gp.gaussPoints[i].second, ret) ;
	}


	
}


Matrix VirtualMachine::ieval( const DtGtMtG & d, IntegrableEntity *e, const std::vector<Variable> & var)
{
	GaussPointArray gp = e->getGaussPoints() ;
	std::valarray<Matrix> Jinv(Matrix(), gp.gaussPoints.size()) ;
	for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
	{
		e->getInverseJacobianMatrix(gp.gaussPoints[i].first, Jinv[i]) ;
	}
	
	return ieval(d, gp, Jinv, var);
}

Matrix VirtualMachine::ieval(const DtGtMtG & d, const GaussPointArray &gp_, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars)
{
	GaussPointArray gp_a(gp_);
	gp_a.id = -1 ;
	GaussPointArray gp_b(gp_);
	gp_b.id = -1 ;
	switch(d.first.v)
	{
		case XI :
		{
			for(size_t i = 0 ; i < gp_a.gaussPoints.size() ; i++)
			{
				gp_a.gaussPoints[i].first.x += default_derivation_delta ;
			}
			
			Matrix x_a(ieval(d.second, gp_a, Jinv, vars)) ;
			
			for(size_t i = 0 ; i < gp_b.gaussPoints.size() ; i++)
			{
				gp_b.gaussPoints[i].first.x -= default_derivation_delta ;
			}
			
			Matrix x_b(ieval(d.second, gp_b, Jinv, vars)) ;
			
			return (x_a-x_b)/(2.*default_derivation_delta) ;
		}
	case ETA :
		{
			for(size_t i = 0 ; i < gp_a.gaussPoints.size() ; i++)
			{
				gp_a.gaussPoints[i].first.y += default_derivation_delta ;
			}
			
			Matrix x_a (ieval(d.second, gp_a, Jinv, vars)) ;
			
			for(size_t i = 0 ; i < gp_b.gaussPoints.size() ; i++)
			{
				gp_b.gaussPoints[i].first.y -= default_derivation_delta ;
			}
			
			Matrix x_b (ieval(d.second, gp_b, Jinv, vars)) ;
			
			return (x_a-x_b)/(2.*default_derivation_delta) ;
		}
	case ZETA :
		{
			for(size_t i = 0 ; i < gp_a.gaussPoints.size() ; i++)
			{
				gp_a.gaussPoints[i].first.z += default_derivation_delta ;
			}
			
			Matrix x_a(ieval(d.second, gp_a, Jinv, vars)) ;
			
			for(size_t i = 0 ; i < gp_b.gaussPoints.size() ; i++)
			{
				gp_b.gaussPoints[i].first.z -= default_derivation_delta ;
			}
			
			Matrix x_b (ieval(d.second, gp_b, Jinv, vars)) ;
			
			return (x_a-x_b)/(2.*default_derivation_delta) ;
		}
	case TIME_VARIABLE :
		{
			for(size_t i = 0 ; i < gp_a.gaussPoints.size() ; i++)
			{
				gp_a.gaussPoints[i].first.t += default_derivation_delta ;
			}
			
			Matrix x_a(ieval(d.second, gp_a, Jinv, vars)) ;
			x_a.print() ;
			
			for(size_t i = 0 ; i < gp_b.gaussPoints.size() ; i++)
			{
				gp_b.gaussPoints[i].first.t -= default_derivation_delta ;
			}
			
			Matrix x_b (ieval(d.second, gp_b, Jinv, vars)) ;
			x_b.print() ;

			return (x_a-x_b)/(2.*default_derivation_delta) ;
		}
	default:
		{
			std::cerr << "operator not implemented" << std::endl ;
			return Matrix() ;
		}
	}
}

Matrix VirtualMachine::ieval(const GDtMtGD & f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars)
{
	std::valarray<Matrix> B = gdeval(f.first.f, Jinv, vars, gp, f.first.transpose) ;
	std::valarray<Matrix> B_ = gdeval(f.third.f, Jinv, vars, gp, f.third.transpose) ;

	Matrix ret( (Matrix)(B[0]*f.second*B_[0])*gp.gaussPoints[0].second) ;

	for(size_t i = 1 ; i  < gp.gaussPoints.size() ; i++)
	{
		ret += (Matrix)(B[i]*f.second*B_[i])*gp.gaussPoints[i].second ;
	}
	
	return ret ;
}

void VirtualMachine::ieval(const GDtMtGD & f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars, Matrix & ret)
{
	std::valarray<Matrix> B = gdeval(f.first.f, Jinv, vars, gp, f.first.transpose) ;
	std::valarray<Matrix> B_ = gdeval(f.third.f, Jinv, vars, gp, f.third.transpose) ;

	ret = (B[0]*f.second*B_[0]);
	ret *= gp.gaussPoints[0].second ;

	for(size_t i = 1 ; i  < gp.gaussPoints.size() ; i++)
	{
		ret += B[i]*f.second*(B_[i]*gp.gaussPoints[i].second) ;
	}
}

Matrix VirtualMachine::ieval(const GDtMtG & f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars)
{
	std::valarray<Matrix> B = gdeval(f.first.f, Jinv, vars, gp, f.first.transpose) ;
	std::valarray<Matrix> B_ = geval(f.third.f, Jinv, vars, gp, f.third.transpose) ;

	Matrix ret( (Matrix)(B[0]*f.second*B_[0])*gp.gaussPoints[0].second) ;

	for(size_t i = 1 ; i  < gp.gaussPoints.size() ; i++)
	{
		ret += (Matrix)(B[i]*f.second*B_[i])*gp.gaussPoints[i].second ;
	}
	
	return ret ;
}

Matrix VirtualMachine::ieval(const GDtMtG & f, IntegrableEntity * e , const std::vector<Variable> & vars)
{
	GaussPointArray gp = e->getGaussPoints() ;
	std::valarray<Matrix> Jinv(Matrix(), gp.gaussPoints.size()) ;
	for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
	{
		e->getInverseJacobianMatrix(gp.gaussPoints[i].first, Jinv[i]) ;
	}
	
	std::valarray<Matrix> B = gdeval(f.first.f, Jinv, vars, gp, f.first.transpose) ;
	std::valarray<Matrix> B_ = geval(f.third.f, Jinv, vars, gp, f.third.transpose) ;

	Matrix ret( (Matrix)(B[0]*f.second*B_[0])*gp.gaussPoints[0].second) ;

	for(size_t i = 1 ; i  < gp.gaussPoints.size() ; i++)
	{
		ret += (Matrix)(B[i]*f.second*B_[i])*gp.gaussPoints[i].second ;
	}
	
	return ret ;
}

void VirtualMachine::ieval(const GDtMtG & f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars, Matrix & ret)
{
	std::valarray<Matrix> B = gdeval(f.first.f, Jinv, vars, gp, f.first.transpose) ;
	std::valarray<Matrix> B_ = geval(f.third.f, Jinv, vars, gp, f.third.transpose) ;

	ret = (B[0]*f.second*B_[0]);
	ret *= gp.gaussPoints[0].second ;

	for(size_t i = 1 ; i  < gp.gaussPoints.size() ; i++)
	{
		ret += B[i]*f.second*(B_[i]*gp.gaussPoints[i].second) ;
	}
}

Matrix VirtualMachine::ieval(const GtMtGD & f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars)
{
	std::valarray<Matrix> B = geval(f.first.f, Jinv, vars, gp, f.first.transpose) ;
	std::valarray<Matrix> B_ = gdeval(f.third.f, Jinv, vars, gp, f.third.transpose) ;

	Matrix ret( (Matrix)(B[0]*f.second*B_[0])*gp.gaussPoints[0].second) ;

	for(size_t i = 1 ; i  < gp.gaussPoints.size() ; i++)
	{
		ret += (Matrix)(B[i]*f.second*B_[i])*gp.gaussPoints[i].second ;
	}
	
	return ret ;
}


Matrix VirtualMachine::ieval(const GtMtGD & f, IntegrableEntity * e , const std::vector<Variable> & vars)
{
	GaussPointArray gp = e->getGaussPoints() ;
	std::valarray<Matrix> Jinv(Matrix(), gp.gaussPoints.size()) ;
	for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
	{
		e->getInverseJacobianMatrix(gp.gaussPoints[i].first, Jinv[i]) ;
	}
	std::valarray<Matrix> B = geval(f.first.f, Jinv, vars, gp, f.first.transpose) ;
	std::valarray<Matrix> B_ = gdeval(f.third.f, Jinv, vars, gp, f.third.transpose) ;

	Matrix ret( (Matrix)(B[0]*f.second*B_[0])*gp.gaussPoints[0].second) ;

	for(size_t i = 1 ; i  < gp.gaussPoints.size() ; i++)
	{
		ret += (Matrix)(B[i]*f.second*B_[i])*gp.gaussPoints[i].second ;
	}
	
	return ret ;
}

void VirtualMachine::ieval(const GtMtGD & f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars, Matrix & ret)
{
	std::valarray<Matrix> B = geval(f.first.f, Jinv, vars, gp, f.first.transpose) ;
	std::valarray<Matrix> B_ = gdeval(f.third.f, Jinv, vars, gp, f.third.transpose) ;

	ret = (B[0]*f.second*B_[0]);
	ret *= gp.gaussPoints[0].second ;

	for(size_t i = 1 ; i  < gp.gaussPoints.size() ; i++)
	{
		ret += B[i]*f.second*(B_[i]*gp.gaussPoints[i].second) ;
	}
}


double VirtualMachine::ieval(const VGtMtVG &f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & var)
{

	Matrix B (gveval(f.first.f, Jinv[0],var, gp.gaussPoints[0].first.x, gp.gaussPoints[0].first.y, gp.gaussPoints[0].first.z, gp.gaussPoints[0].first.t, f.first.transpose));
	Matrix B_(gveval(f.third.f, Jinv[0],var, gp.gaussPoints[0].first.x, gp.gaussPoints[0].first.y, gp.gaussPoints[0].first.z, gp.gaussPoints[0].first.t, f.third.transpose)) ;

	Matrix r_  (B*f.second*B_) ;
	double ret = r_[0][0] * gp.gaussPoints[0].second;
	
	for(size_t i = 1 ; i < gp.gaussPoints.size() ; i++)
	{
		B  = gveval(f.first.f, Jinv[i],var, gp.gaussPoints[i].first.x, gp.gaussPoints[i].first.y, gp.gaussPoints[i].first.z, gp.gaussPoints[i].first.t, f.first.transpose) ;
		B_ = gveval(f.third.f, Jinv[i],var, gp.gaussPoints[i].first.x, gp.gaussPoints[i].first.y, gp.gaussPoints[i].first.z, gp.gaussPoints[i].first.t, f.third.transpose) ;
		Matrix r_ (B*f.second*B_) ;
		ret += r_[0][0] * gp.gaussPoints[i].second ;
	}
	
	return ret ;
}

Vector VirtualMachine::ieval(const GtV &f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & var)
{
	if(!Jinv.size())
	{
		int size = var.size() - (var[var.size()-1] == TIME_VARIABLE) ;
		return Vector(double(0), size) ;
	}
	Matrix M ;
	geval(f.first.f, Jinv[0],var, gp.gaussPoints[0].first.x, gp.gaussPoints[0].first.y, gp.gaussPoints[0].first.z, gp.gaussPoints[0].first.t, f.first.transpose, M);
	Vector temp(double(0),f.second.size()) ;
	for(size_t i = 0 ; i < temp.size() ; i++)
		temp[i] = f.second[i] ;
	
	Vector r_  = M*temp ;
	Vector ret = r_ * gp.gaussPoints[0].second;
	for(size_t i = 1 ; i < gp.gaussPoints.size() ; i++)
	{
		  geval(f.first.f, Jinv[i],var, gp.gaussPoints[i].first.x, gp.gaussPoints[i].first.y, gp.gaussPoints[i].first.z, gp.gaussPoints[i].first.t, f.first.transpose, M) ;
		
		for(size_t j = 0 ; j < temp.size() ; j++)
		{
			temp[j] = f.second[j] ;
		}
		r_  = M*temp ;
		ret += r_ * gp.gaussPoints[i].second ;
	}
	
	return ret ;
}

double VirtualMachine::ieval(const VGtV &f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & var)
{

	Matrix B (gveval(f.first.f, Jinv[0],var, gp.gaussPoints[0].first.x, gp.gaussPoints[0].first.y, gp.gaussPoints[0].first.z, gp.gaussPoints[0].first.t, f.first.transpose));
	
	Vector temp(B.numCols()) ;
	for(size_t i = 0 ; i < temp.size() ; i++)
		temp[i] = f.second[i] ;
	
	Vector r_  = B*temp ;
	double ret = r_[0] * gp.gaussPoints[0].second;
	
	for(size_t i = 1 ; i < gp.gaussPoints.size() ; i++)
	{
		B  = gveval(f.first.f, Jinv[i],var, gp.gaussPoints[i].first.x, gp.gaussPoints[i].first.y, gp.gaussPoints[i].first.z, gp.gaussPoints[i].first.t, f.first.transpose) ;
		for(size_t j = 0 ; j < temp.size() ; j++)
			temp[j] = f.second[i*temp.size()+j] ;
		Vector r_  = B*temp ;
		ret += r_[0] * gp.gaussPoints[i].second ;
	}
	
	return ret ;
}

std::vector<Point> VirtualMachine::allHints(const Function &f0, const Function &f1,IntegrableEntity *e )
{
	std::vector<Point> hints(0) ;
	
	if ((e!= NULL && f0.hasIntegrationHint()) || f1.hasIntegrationHint())
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

std::vector<Point> VirtualMachine::allHints(const Function &f0,IntegrableEntity *e )
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

Matrix VirtualMachine::ieval(const GtFMtG &f, IntegrableEntity *e, const std::vector<Variable> & var)
{

	GaussPointArray gp = e->getGaussPoints() ;
	std::valarray<Matrix> Jinv(Matrix(), gp.gaussPoints.size()) ;
	for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
	{
		e->getInverseJacobianMatrix(gp.gaussPoints[i].first, Jinv[i]) ;
	}
	
	return ieval(f, gp, Jinv,var) ;
	
}

Matrix VirtualMachine::ieval(const GtFMtG &f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & var)
{
	Matrix B(geval(f.first.f, Jinv[0],var, gp.gaussPoints[0].first.x, gp.gaussPoints[0].first.y, gp.gaussPoints[0].first.z, gp.gaussPoints[0].first.t, f.first.transpose));
	Matrix B_(geval(f.third.f, Jinv[0],var, gp.gaussPoints[0].first.x, gp.gaussPoints[0].first.y, gp.gaussPoints[0].first.z, gp.gaussPoints[0].first.t, f.third.transpose)) ;
	
	Matrix r_ (B*eval(f.second,gp.gaussPoints[0].first )*B_) ;
	Matrix ret(r_ * gp.gaussPoints[0].second);
	
	for(size_t i = 1 ; i < gp.gaussPoints.size() ; i++)
	{
		B  = geval(f.first.f, Jinv[i],var, gp.gaussPoints[i].first.x, gp.gaussPoints[i].first.y, gp.gaussPoints[i].first.z, gp.gaussPoints[i].first.t, f.first.transpose) ;
		B_ = geval(f.third.f, Jinv[i],var, gp.gaussPoints[i].first.x, gp.gaussPoints[i].first.y, gp.gaussPoints[i].first.z, gp.gaussPoints[i].first.t, f.third.transpose) ;
		Matrix r_ (B*eval(f.second,gp.gaussPoints[i].first )*B_) ;
		ret += r_ * gp.gaussPoints[i].second ;
	}
	
	return ret ;
}

void VirtualMachine::ieval(const GtFMtG &f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & var, Matrix &ret)
{

	geval(f.third.f, Jinv[0], var, gp.gaussPoints[0].first, f.third.transpose,B_) ;
	geval(f.first.f, Jinv[0], var, gp.gaussPoints[0].first, f.first.transpose,B) ;
	
	ret = B*eval(f.second, gp.gaussPoints[0].first)*B_ ;
	ret *= gp.gaussPoints[0].second ;
	
	for(size_t i = 1 ; i  < gp.gaussPoints.size() ; i++)
	{
		geval(f.third.f, Jinv[i], var, gp.gaussPoints[i].first, f.third.transpose, B_) ;
		geval(f.first.f, Jinv[i], var, gp.gaussPoints[i].first, f.first.transpose, B) ;
		matrix_matrix_matrix_multiply_and_add(B, eval(f.second, gp.gaussPoints[0].first), B_, gp.gaussPoints[i].second, ret) ;
	}
}


double VirtualMachine::ieval(const Differential & d, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & var)
{

	double ret = 0 ;
	size_t var_line = 0 ;
	for(size_t i = 0 ; i < var.size() ; i++)
	{
		if(var[i] == d.v)
		{
			var_line = i ;
			break ;
		}
	}
	
	for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
	{
		for(size_t j = 0 ; j < Jinv[i].numCols() ; j++)
		{
			ret += deval(d.f, var[j], gp.gaussPoints[i].first)*Jinv[i][var_line][j]*gp.gaussPoints[i].second ;
		}
	}
	
	return ret ;
}


double VirtualMachine::ieval(const DtF & d, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & var)
{
	double ret = 0 ;
	size_t var_line = 0 ;
	for(size_t i = 0 ; i < var.size() ; i++)
	{
		if(var[i] == d.d.v)
		{
			var_line = i ;
			break ;
		}
	}
	for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
	{
		double fx = eval(d.f, gp.gaussPoints[i].first) ;
		for(size_t j = 0 ; j < Jinv[i].numCols() ; j++)
		{
			ret += deval(d.d.f, var[j], gp.gaussPoints[i].first)*fx*Jinv[i][var_line][j]*gp.gaussPoints[i].second ;
		}
	}
	
	return ret ;
}

double VirtualMachine::ieval(const DtD & d, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & var)
{
	double ret = 0 ;
	size_t var_line = 0 ;
	for(size_t i = 0 ; i < var.size() ; i++)
	{
		if(var[i] == d.d.v)
		{
			var_line = i ;
			break ;
		}
	}
	for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
	{
// 		print(d.d.f) ;
// 		print(d.f.f) ;
		for(size_t j = 0 ; j < Jinv[i].numCols() ; j++)
		{

			ret += deval(d.d.f, var[j], gp.gaussPoints[i].first)
				*deval(d.f.f, var[j], gp.gaussPoints[i].first)
				*Jinv[i][var_line][j]*gp.gaussPoints[i].second ;
// 				std::cout << deval(d.d.f, var[j], gp.gaussPoints[i].first) << "   " << deval(d.f.f, var[j], gp.gaussPoints[i].first) << std::endl ;
		}
	}
		
	return ret ;

}

double VirtualMachine::ieval(const Differential & d, IntegrableEntity *e, const std::vector<Variable> & var)
{
	GaussPointArray gp = e->getGaussPoints() ;
	double ret = 0 ;
	size_t var_line = 0 ;
	
	for(size_t i = 0 ; i < var.size() ; i++)
	{
		if(var[i] == d.v)
		{
			var_line = i ;
			break ;
		}
	}
	
	for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
	{
		Matrix Jinv ; e->getInverseJacobianMatrix(gp.gaussPoints[i].first, Jinv);
		for(size_t j = 0 ; j < Jinv.numCols() ; j++)
		{
			ret += deval(d.f, var[j], gp.gaussPoints[i].first)*Jinv[var_line][j]*gp.gaussPoints[i].second ;
		}
	}
	
	return ret ;
}

double VirtualMachine::ieval(const DtF & d, IntegrableEntity *e, const std::vector<Variable> & var)
{
	GaussPointArray gp = e->getGaussPoints() ;
	
	double ret = 0 ;
	size_t var_line = 0 ;
	
	for(size_t i = 0 ; i < var.size() ; i++)
	{
		if(var[i] == d.d.v)
		{
			var_line = i ;
			break ;
		}
	}
	
	for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
	{
		Matrix Jinv ; e->getInverseJacobianMatrix(gp.gaussPoints[i].first, Jinv);
		for(size_t j = 0 ; j < Jinv.numCols() ; j++)
		{
			ret += deval(d.d.f, var[j], gp.gaussPoints[i].first)*
				eval(d.f, gp.gaussPoints[i].first)*Jinv[var_line][j]*gp.gaussPoints[i].second ;
		}
	}
	
	return ret ;
	
}

double VirtualMachine::ieval(const DtD & d, IntegrableEntity *e, const std::vector<Variable> & var)
{
	GaussPointArray gp = e->getGaussPoints() ;
	
	double ret = 0 ;
	size_t var_line = 0 ;
	
	for(size_t i = 0 ; i < var.size() ; i++)
	{
		if(var[i] == d.d.v)
		{
			var_line = i ;
			break ;
		}
	}
	
	for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
	{
		Matrix Jinv ; e->getInverseJacobianMatrix(gp.gaussPoints[i].first, Jinv);
		for(size_t j = 0 ; j < Jinv.numCols() ; j++)
		{
			ret += deval(d.d.f, var[j], gp.gaussPoints[i].first)*
				deval(d.f.f, var[j], gp.gaussPoints[i].first)*Jinv[var_line][j]*gp.gaussPoints[i].second ;
		}
	}
	
	return ret ;
	
}


