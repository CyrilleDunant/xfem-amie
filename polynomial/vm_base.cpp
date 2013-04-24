// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution



#include "vm_base.h"
#include "../mesher/delaunay.h"
#include <limits>
#include <iomanip>


using namespace Mu ;

VirtualMachine::VirtualMachine(){ } ;


double VirtualMachine::eval(const Function &f, const double x, const double y, const double z, const double t, const double u, const double v, const double w) 
{
	size_t size = f.byteCodeSize ;
	stack.memory.heap[1] = x ;
	stack.memory.heap[2] = y ;
	stack.memory.heap[3] = z ;
	stack.memory.heap[4] = t ;
	stack.memory.heap[5] = u ;
	stack.memory.heap[6] = v ;
	stack.memory.heap[7] = w ;
	
	std::reverse_copy(&f.values[0],&f.values[f.constNumber],&stack.memory.heap[FUNCTION_LENGTH-f.constNumber] ) ;
// 	for(size_t i = 0 ; i < f.constNumber  ; ++i)
// 	{
// 		stack.memory.heap[511-i] = f.values[i] ;
// 	}
	  
	for(size_t i = 0 ; i < size  ; ++i)
	{
#define REG_A stack.memory.heap[f.adress_a[i*4]]
#define REG_B stack.memory.heap[f.adress_a[i*4+1]]
#define REG_C stack.memory.heap[f.adress_a[i*4+2]]
#define REG_0 stack.memory.heap[0]
		if(f.byteCode[i] >= TOKEN_OPERATION_PLUS && f.byteCode[i] <= TOKEN_OPERATION_INPLACE_POWER)
		{
			if(f.use_temp[i] == NO_TEMPORARY)
			{
				if (f.byteCode[i] ==  TOKEN_OPERATION_TIMES)
				{
					REG_C = REG_A*REG_B ;
				}
				else if (f.byteCode[i] ==  TOKEN_OPERATION_INPLACE_PLUS)
				{
					REG_A += REG_B ;
				}
				else if (f.byteCode[i] ==  TOKEN_OPERATION_POWER)
				{
					double val = REG_A ;
					int pow = REG_B-1 ;
					REG_C = val ;
					for(int j = 0 ; j < pow ; ++j)
					{
						REG_C *= val;
					}
				}
				else if (f.byteCode[i] ==  TOKEN_OPERATION_INPLACE_TIMES)
				{
					REG_A *= REG_B ;
				}
				else if (f.byteCode[i] == TOKEN_OPERATION_PLUS)
				{
					REG_C = REG_A+REG_B ;
				}
				else if (f.byteCode[i] == TOKEN_OPERATION_MINUS)
				{
					REG_C = REG_A-REG_B ;
				}
				else if (f.byteCode[i] ==  TOKEN_OPERATION_DIVIDES)
				{
					REG_C = REG_A/REG_B ;
				}
				else if (f.byteCode[i] ==  TOKEN_OPERATION_INPLACE_MINUS)
				{
					REG_A -= REG_B ;
				}
				else 
				{
					REG_A /= REG_B ;
				}
			}
			else if (f.use_temp[i] == SET_TEMPORARY)
			{
				if (f.byteCode[i] ==  TOKEN_OPERATION_TIMES)
				{
					REG_0 = REG_A*REG_B ;
				}
				else if (f.byteCode[i] ==  TOKEN_OPERATION_INPLACE_PLUS)
				{
					REG_0 += REG_B ;
				}
				else if (f.byteCode[i] ==  TOKEN_OPERATION_POWER)
				{
					double val = REG_A ;
					int pow = REG_B-1 ;
					REG_0 = val ;
					for(int j = 0 ; j < pow ; ++j)
					{
						REG_0 *= val;
					}
				}
				else if (f.byteCode[i] == TOKEN_OPERATION_PLUS)
				{
					REG_0 = REG_A+REG_B ;
				}
				else if (f.byteCode[i] ==  TOKEN_OPERATION_INPLACE_TIMES)
				{
					REG_0 *= REG_B ;
				}
				else if (f.byteCode[i] == TOKEN_OPERATION_MINUS)
				{
					REG_0 = REG_A-REG_B ;
				}
				else if (f.byteCode[i] ==  TOKEN_OPERATION_DIVIDES)
				{
					REG_0 = REG_A/REG_B ;
				}
				else if (f.byteCode[i] ==  TOKEN_OPERATION_INPLACE_MINUS)
				{
					REG_0 -= REG_B ;
				}
				else 
				{
					REG_0 /= REG_B ;
				}
			}
			else if(f.use_temp[i] == GET_TEMPORARY_A)
			{
				if (f.byteCode[i] ==  TOKEN_OPERATION_TIMES)
				{
					REG_C = REG_0*REG_B ;
				}
				else if (f.byteCode[i] == TOKEN_OPERATION_PLUS)
				{
					REG_C = REG_0+REG_B ;
				}
				else if (f.byteCode[i] ==  TOKEN_OPERATION_MINUS)
				{
					REG_C = REG_0-REG_B ;
				}
				else if (f.byteCode[i] ==  TOKEN_OPERATION_DIVIDES)
				{
					REG_C = REG_0/REG_B ;
				}
				else
				{
					double val = REG_0 ;
					int pow = REG_B-1 ;
					REG_C = val ;
					for(int j = 0 ; j < pow ; ++j)
					{
						REG_C *= val;
					}
				}
			}
			else if(f.use_temp[i] == GET_TEMPORARY_B)
			{
				 if (f.byteCode[i] == TOKEN_OPERATION_TIMES)
				{
					REG_C = REG_A*REG_0 ;
				}
				else if (f.byteCode[i] == TOKEN_OPERATION_PLUS)
				{
					REG_C = REG_A+REG_0 ;
				}
				else if (f.byteCode[i] == TOKEN_OPERATION_INPLACE_PLUS)
				{
					REG_A += REG_0 ;
				}
				else if (f.byteCode[i] == TOKEN_OPERATION_POWER)
				{
					double val = REG_A ;
					int pow = REG_0-1 ;
					REG_C = val ;
					for( int j = 0 ; j < pow ; ++j)
					{
						REG_C *= val;
					}
				}
				else if (f.byteCode[i] == TOKEN_OPERATION_MINUS)
				{
					REG_C = REG_A-REG_0 ;
				}
				else if (f.byteCode[i] == TOKEN_OPERATION_DIVIDES)
				{
					REG_C = REG_A/REG_0 ;
				}
				else if (f.byteCode[i] == TOKEN_OPERATION_INPLACE_MINUS)
				{
					REG_A -= REG_0 ;
				}
				else if (f.byteCode[i] == TOKEN_OPERATION_INPLACE_TIMES)
				{
					REG_A *= REG_0 ;
				}
				else
				{
					REG_A /= REG_0 ;
				}
			}
			else if(f.use_temp[i] == SET_GET_TEMPORARY_A)
			{
				if (f.byteCode[i] == TOKEN_OPERATION_TIMES)
				{
					REG_0 *= REG_B ;
				}
				else if (f.byteCode[i] == TOKEN_OPERATION_INPLACE_PLUS)
				{
					REG_0 += REG_B ;
				}
				else if (f.byteCode[i] == TOKEN_OPERATION_POWER)
				{
					double val = REG_0 ;
					int pow = REG_B-1 ;
					REG_0 = val ;
					for( int j = 0 ; j < pow ; ++j)
					{
						REG_0 *= val;
					}
				}
				else if (f.byteCode[i] == TOKEN_OPERATION_INPLACE_TIMES)
				{
					REG_0 *= REG_B ;
				}
				else if (f.byteCode[i] == TOKEN_OPERATION_PLUS)
				{
					REG_0 += REG_B ;
				}
				else if (f.byteCode[i] == TOKEN_OPERATION_MINUS)
				{
					REG_0 -= REG_B ;
				}
				else if (f.byteCode[i] == TOKEN_OPERATION_DIVIDES)
				{
					REG_0 /= REG_B ;
				}
				else if (f.byteCode[i] == TOKEN_OPERATION_INPLACE_MINUS)
				{
					REG_0 -= REG_B ;
				}
				else
				{
					REG_0 /= REG_B ;
				}
			}
			else 
			{
				if (f.byteCode[i] == TOKEN_OPERATION_TIMES)
				{
					REG_0 *= REG_A ;
				}
				else if (f.byteCode[i] == TOKEN_OPERATION_POWER)
				{
					double val = REG_A ;
					int pow = REG_0-1 ;
					REG_0 = val ;
					for( int j = 0 ; j < pow ; ++j)
					{
						REG_0 *= val;
					}
				}
				else if (f.byteCode[i] == TOKEN_OPERATION_PLUS)
				{
					REG_0 += REG_A ;
				}
				else if (f.byteCode[i] == TOKEN_OPERATION_INPLACE_TIMES)
				{
					REG_A *= REG_0 ;
				}
				else if (f.byteCode[i] == TOKEN_OPERATION_MINUS)
				{
					REG_0 = REG_A-REG_0 ;
				}
				else if (f.byteCode[i] == TOKEN_OPERATION_DIVIDES)
				{
					REG_0 = REG_A/REG_0 ;
				}
				else if (f.byteCode[i] == TOKEN_OPERATION_INPLACE_MINUS)
				{
					REG_A -= REG_0 ;
				}
				else
				{
					REG_A /= REG_0 ;
				}
			}
			continue ;
		}

		if(f.byteCode[i] == TOKEN_OPERATION_CONSTANT ||
			f.byteCode[i] == TOKEN_OPERATION_X ||
			f.byteCode[i] == TOKEN_OPERATION_Y ||
			f.byteCode[i] == TOKEN_OPERATION_Z ||
			f.byteCode[i] == TOKEN_OPERATION_T ||
			f.byteCode[i] == TOKEN_OPERATION_U ||
			f.byteCode[i] == TOKEN_OPERATION_V ||
			f.byteCode[i] == TOKEN_OPERATION_W 
		)
		{
			REG_C=REG_A ;
			continue ;
		}
		
		switch(f.byteCode[i])
		{
			case TOKEN_OPERATION_GEO_OPERATION:
			{
				f.geo_op[i]->eval(&REG_A, &REG_B, &REG_C) ;
				break ;
			}
			case TOKEN_OPERATION_COS:
			{
				REG_C = cos(REG_A) ;
				break ;
			}
			case TOKEN_OPERATION_ABS:
			{
				REG_C = std::abs(REG_A) ;
				break ;
			}
			case TOKEN_OPERATION_TAN:
			{
				REG_C = tan(REG_A) ;
				break ;
			}
			case TOKEN_OPERATION_SIN:
			{
				REG_C = sin(REG_A) ;
				break ;
			}
			case TOKEN_OPERATION_EXP:
			{
				REG_C = exp(REG_A) ;
				break ;
			}
			case TOKEN_OPERATION_SIGN:
			{
				REG_C = sign(REG_A) ;
				break ;
			}
			case TOKEN_OPERATION_POSITIVITY:
			{
				double s0 = sign(REG_A) ;
				s0 >= 0 ? REG_C = 1 : REG_C = 0 ;
				break ;
			}
			case TOKEN_OPERATION_NEGATIVITY:
			{
				double s0 = sign(REG_A) ;
				s0 < 0 ? REG_C = 1 : REG_C = 0 ;
				break ;
			}
			case TOKEN_OPERATION_LOG:
			{
				REG_C = log(REG_A) ;
				break ;
			}
			case TOKEN_OPERATION_COSH:
			{
				REG_C = cosh(REG_A) ;
				break ;
			}
			case TOKEN_OPERATION_SINH:
			{
				REG_C = sinh(REG_A) ;
				break ;
			}
			case TOKEN_OPERATION_TANH:
			{
				REG_C = tanh(REG_A) ;
				break ;
			}
			case TOKEN_OPERATION_SQRT:
			{
				REG_C = sqrt(REG_A) ;
				break ;
			}
			case TOKEN_OPERATION_ATAN2:
			{
				REG_C = atan2(REG_B, REG_A) ;

				break ;
			}
		}

	}

	return  stack.memory.heap[8] ;
}

Vector VirtualMachine::eval(const Function &f, const GaussPointArray &gp)
{
	if(f.precalculated(gp))
		return f.getPrecalculatedValue(gp) ;
	
	Vector ret(gp.gaussPoints.size()) ;
	for(size_t i = 0 ; i < ret.size() ; i++)
	{
		ret[i] = eval(f, gp.gaussPoints[i].first.x,
				 gp.gaussPoints[i].first.y,
				 gp.gaussPoints[i].first.z,
				 gp.gaussPoints[i].first.t) ;
	}

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
	std::cout << f.size() << " instructions:\n" << std::endl ;
	for(size_t i = 0 ; i < f.size()  ; i++)
	{
		bool done = false ;
		switch(f.byteCode[i])
		{
			case TOKEN_OPERATION_CONSTANT:
			{
				std::cout << "sto " << f.values[i*4] << "@" << f.adress_a[i*4]+100 << "  "<< std::endl ;
				done = true ;
				break ;
			}
			case TOKEN_OPERATION_X:
			{
				std::cout << "sto " << "x" << "@" << f.adress_a[i*4]+100 << "  "<< std::endl ;
				done = true ;
				break ;
			}
			case TOKEN_OPERATION_Y:
			{
				std::cout << "sto " << "y" << "@" << f.adress_a[i*4]+100 << "  "<< std::endl ;
				done = true ;
				break ;
			}
			case TOKEN_OPERATION_Z:
			{
				std::cout << "sto " << "z" << "@" << f.adress_a[i*4]+100 << "  "<< std::endl ;
				done = true ;
				break ;
			}
			case TOKEN_OPERATION_T:
			{
				std::cout << "sto " << "t" << "@" << f.adress_a[i*4]+100 << "  "<< std::endl ;
				done = true ;
				break ;
			}
		}
	
		if(done)
			continue ;
		
		switch(f.byteCode[i])
		{
			case TOKEN_OPERATION_PLUS:
			{
				if(f.use_temp[i] == NO_TEMPORARY)
					std::cout << "add " << "@" << f.adress_a[i*4]+100 <<"    @" << f.adress_a[i*4+1]+100 << "    : @" << f.adress_a[i*4+2]+100 << "  "<< std::endl ;
				if(f.use_temp[i] == SET_TEMPORARY)
					std::cout << "add " << "@" << f.adress_a[i*4]+100 <<"    @" << f.adress_a[i*4+1]+100 << "    : @tmp" << "  "<< std::endl ;
				if(f.use_temp[i] == GET_TEMPORARY_A)
					std::cout << "add " << "@tmp" <<"    @" << f.adress_a[i*4+1]+100 << "    : @" << f.adress_a[i*4+2]+100 << "  "<< std::endl ;
				if(f.use_temp[i] == GET_TEMPORARY_B)
					std::cout << "add " << "@" << f.adress_a[i*4]+100 <<"    @tmp" << "    : @" << f.adress_a[i*4+2]+100 << "  "<< std::endl ;
				if(f.use_temp[i] == SET_GET_TEMPORARY_A)
					std::cout << "add " << "@tmp" <<"    @" << f.adress_a[i*4+1]+100 << "    : @tmp" << "  "<< std::endl ;
				if(f.use_temp[i] == SET_GET_TEMPORARY_B)
					std::cout << "add " << "@" << f.adress_a[i*4]+100 <<"    @tmp" << "    : @tmp" << "  "<< std::endl ;
				done = true ;
				break ;
			}
			case TOKEN_OPERATION_MINUS:
			{
				if(f.use_temp[i] == NO_TEMPORARY)
					std::cout << "sub " << "@" << f.adress_a[i*4]+100 <<"    @" << f.adress_a[i*4+1]+100 << "    : @" << f.adress_a[i*4+2]+100 << "  "<< std::endl ;
				if(f.use_temp[i] == SET_TEMPORARY)
					std::cout << "sub " << "@" << f.adress_a[i*4]+100 <<"    @" << f.adress_a[i*4+1]+100 << "    : @tmp" << "  "<< std::endl ;
				if(f.use_temp[i] == GET_TEMPORARY_A)
					std::cout << "sub " << "@tmp" <<"    @" << f.adress_a[i*4+1]+100 << "    : @" << f.adress_a[i*4+2]+100 << "  "<< std::endl ;
				if(f.use_temp[i] == GET_TEMPORARY_B)
					std::cout << "sub " << "@" << f.adress_a[i*4]+100 <<"    @tmp" << "    : @" << f.adress_a[i*4+2]+100 << "  "<< std::endl ;
				if(f.use_temp[i] == SET_GET_TEMPORARY_A)
					std::cout << "sub " << "@tmp" <<"    @" << f.adress_a[i*4+1]+100 << "    : @tmp" << "  "<< std::endl ;
				if(f.use_temp[i] == SET_GET_TEMPORARY_B)
					std::cout << "sub " << "@" << f.adress_a[i*4]+100 <<" @tmp" << "    : @tmp" << "  "<< std::endl ;
				done = true ;
				break ;
			}
			case TOKEN_OPERATION_TIMES:
			{
				if(f.use_temp[i] == NO_TEMPORARY)
					std::cout << "mul " << "@" << f.adress_a[i*4]+100 <<"    @" << f.adress_a[i*4+1]+100 << "    : @" << f.adress_a[i*4+2]+100 << "  "<< std::endl ;
				if(f.use_temp[i] == SET_TEMPORARY)
					std::cout << "mul " << "@" << f.adress_a[i*4]+100 <<"    @" << f.adress_a[i*4+1]+100 << "    : @tmp" << "  "<< std::endl ;
				if(f.use_temp[i] == GET_TEMPORARY_A)
					std::cout << "mul " << "@tmp" <<"    @" << f.adress_a[i*4+1]+100 << "    : @" << f.adress_a[i*4+2]+100 << "  "<< std::endl ;
				if(f.use_temp[i] == GET_TEMPORARY_B)
					std::cout << "mul " << "@" << f.adress_a[i*4]+100 <<"     @tmp" << "    : @" << f.adress_a[i*4+2]+100 << "  "<< std::endl ;
				if(f.use_temp[i] == SET_GET_TEMPORARY_A)
					std::cout << "mul " << "@tmp" <<"    @" << f.adress_a[i*4+1]+100 << "    : @tmp" << "  "<< std::endl ;
				if(f.use_temp[i] == SET_GET_TEMPORARY_B)
					std::cout << "mul " << "@" << f.adress_a[i*4]+100 <<"    @tmp" << "    : @tmp" << "  "<< std::endl ;
				done = true ;
				break ;
			}
			case TOKEN_OPERATION_DIVIDES:
			{
				if(f.use_temp[i] == NO_TEMPORARY)
					std::cout << "div " << "@" << f.adress_a[i*4]+100 <<"    @" << f.adress_a[i*4+1]+100 << "    : @" << f.adress_a[i*4+2]+100 << "  "<< std::endl ;
				if(f.use_temp[i] == SET_TEMPORARY)
					std::cout << "div " << "@" << f.adress_a[i*4]+100 <<"    @" << f.adress_a[i*4+1]+100 << "    : @tmp" << "  "<< std::endl ;
				if(f.use_temp[i] == GET_TEMPORARY_A)
					std::cout << "div " << "@tmp" <<"    @" << f.adress_a[i*4+1]+100 << "    : @" << f.adress_a[i*4+2]+100 << "  "<< std::endl ;
				if(f.use_temp[i] == GET_TEMPORARY_B)
					std::cout << "div " << "@" << f.adress_a[i*4]+100 <<"    @tmp" << "    : @" << f.adress_a[i*4+2]+100 << "  "<< std::endl ;
				if(f.use_temp[i] == SET_GET_TEMPORARY_A)
					std::cout << "div " << "@tmp" <<"    @" << f.adress_a[i*4+1]+100 << "    : @tmp" << "  "<< std::endl ;
				if(f.use_temp[i] == SET_GET_TEMPORARY_B)
					std::cout << "div " << "@" << f.adress_a[i*4]+100 <<"    @tmp" << "    : @tmp" << "  "<< std::endl ;
				done = true ;
				break ;
			}
			case TOKEN_OPERATION_POWER:
			{
				if(f.use_temp[i] == NO_TEMPORARY)
					std::cout << "pow " << "@" << f.adress_a[i*4]+100 <<"    @" << f.adress_a[i*4+1]+100 << "    : @" << f.adress_a[i*4+2]+100 << "  "<< std::endl ;
				if(f.use_temp[i] == SET_TEMPORARY)
					std::cout << "pow " << "@" << f.adress_a[i*4]+100 <<"    @" << f.adress_a[i*4+1]+100 << "    : @tmp" << "  "<< std::endl ;
				if(f.use_temp[i] == GET_TEMPORARY_A)
					std::cout << "pow " << "@tmp" <<"    @" << f.adress_a[i*4+1]+100 << "    : @" << f.adress_a[i*4+2]+100 << "  "<< std::endl ;
				if(f.use_temp[i] == GET_TEMPORARY_B)
					std::cout << "pow " << "@" << f.adress_a[i*4]+100 <<"     @tmp" << "    : @" << f.adress_a[i*4+2]+100 << "  "<< std::endl ;
				if(f.use_temp[i] == SET_GET_TEMPORARY_A)
					std::cout << "pow " << "@tmp" <<"    @" << f.adress_a[i*4+1]+100 << "    : @tmp" << "  "<< std::endl ;
				if(f.use_temp[i] == SET_GET_TEMPORARY_B)
					std::cout << "pow " << "@" << f.adress_a[i*4]+100 <<"    @tmp" << "    : @tmp" << "  "<< std::endl ;
				done = true ;
				break ;
			}
			case TOKEN_OPERATION_INPLACE_PLUS:
			{
				if(f.use_temp[i] == NO_TEMPORARY)
					std::cout << "add " << "@" << f.adress_a[i*4+1]+100 <<"     : @" << f.adress_a[i*4]+100 << "  "<< std::endl ;
				if(f.use_temp[i] == SET_TEMPORARY)
					std::cout << "add " << "@" << f.adress_a[i*4+1]+100 <<"     : @tmp" << "  "<< std::endl ;
				if(f.use_temp[i] == GET_TEMPORARY_A)
					std::cout << "add " << "@" << f.adress_a[i*4+1]+100 <<"     : @tmp"  << "  "<< std::endl ;
				if(f.use_temp[i] == GET_TEMPORARY_B)
					std::cout << "add " << "@tmp" <<"            : @" << f.adress_a[i*4]+100 << "  "<< std::endl ;
				if(f.use_temp[i] == SET_GET_TEMPORARY_A)
					std::cout << "add " << "@" << f.adress_a[i*4+1]+100 <<"     : @tmp"  << "  "<< std::endl ;
				if(f.use_temp[i] == SET_GET_TEMPORARY_B)
					std::cout << "add " << "@tmp" <<"     : @tmp"  << "  "<< std::endl ;
				done = true ;
				break ;
			}
			case TOKEN_OPERATION_INPLACE_MINUS:
			{
				if(f.use_temp[i] == NO_TEMPORARY)
					std::cout << "sub " << "@" << f.adress_a[i*4+1]+100 <<"     : @" << f.adress_a[i*4]+100 << "  "<< std::endl ;
				if(f.use_temp[i] == SET_TEMPORARY)
					std::cout << "sub " << "@" << f.adress_a[i*4+1]+100 <<"     : @tmp" << "  "<< std::endl ;
				if(f.use_temp[i] == GET_TEMPORARY_A)
					std::cout << "sub " << "@" << f.adress_a[i*4+1]+100 <<"     : @tmp"  << "  "<< std::endl ;
				if(f.use_temp[i] == GET_TEMPORARY_B)
					std::cout << "sub " << "@tmp" <<"            : @" << f.adress_a[i*4]+100 << "  "<< std::endl ;
				if(f.use_temp[i] == SET_GET_TEMPORARY_A)
					std::cout << "sub " << "@" << f.adress_a[i*4+1]+100 <<"     : @tmp"  << "  "<< std::endl ;
				if(f.use_temp[i] == SET_GET_TEMPORARY_B)
					std::cout << "sub " << "@tmp" <<"     : @tmp"  << "  "<< std::endl ;
				done = true ;
				break ;
			}
			case TOKEN_OPERATION_INPLACE_TIMES:
			{
				if(f.use_temp[i] == NO_TEMPORARY)
					std::cout << "mul " << " " << f.adress_a[i*4+1]+100 <<"     : @" << f.adress_a[i*4]+100 << "  "<< std::endl ;
				if(f.use_temp[i] == SET_TEMPORARY)
					std::cout << "mul " << "@" << f.adress_a[i*4+1]+100 <<"     : @tmp" << "  "<< std::endl ;
				if(f.use_temp[i] == GET_TEMPORARY_A)
					std::cout << "mul " << "@" << f.adress_a[i*4+1]+100 <<"     : @tmp"  << "  "<< std::endl ;
				if(f.use_temp[i] == GET_TEMPORARY_B)
					std::cout << "mul " << "@tmp" <<"            : @" << f.adress_a[i*4]+100 << "  "<< std::endl ;
				if(f.use_temp[i] == SET_GET_TEMPORARY_A)
					std::cout << "mul " << "@" << f.adress_a[i*4+1]+100 <<"     : @tmp"  << "  "<< std::endl ;
				if(f.use_temp[i] == SET_GET_TEMPORARY_B)
					std::cout << "mul " << "@tmp" <<"     : @tmp"  << "  "<< std::endl ;
				done = true ;
				break ;
			}
			case TOKEN_OPERATION_INPLACE_DIVIDES:
			{
				if(f.use_temp[i] == NO_TEMPORARY)
					std::cout << "div " << "@" << f.adress_a[i*4+1]+100 <<"     : @" << f.adress_a[i*4]+100 << "  "<< std::endl ;
				if(f.use_temp[i] == SET_TEMPORARY)
					std::cout << "div " << "@" << f.adress_a[i*4+1]+100 <<"     : @tmp" << "  "<< std::endl ;
				if(f.use_temp[i] == GET_TEMPORARY_A)
					std::cout << "div " << "@" << f.adress_a[i*4+1]+100 <<"     : @tmp"  << "  "<< std::endl ;
				if(f.use_temp[i] == GET_TEMPORARY_B)
					std::cout << "div " << "@tmp" <<"            : @" << f.adress_a[i*4]+100 << "  "<< std::endl ;
				if(f.use_temp[i] == SET_GET_TEMPORARY_A)
					std::cout << "div " << "@" << f.adress_a[i*4+1]+100 <<"     : @tmp"  << "  "<< std::endl ;
				if(f.use_temp[i] == SET_GET_TEMPORARY_B)
					std::cout << "div " << "@tmp" <<"     : @tmp"  << "  "<< std::endl ;
				done = true ;
				break ;
			}
			case TOKEN_OPERATION_INPLACE_POWER:
			{
				if(f.use_temp[i] == NO_TEMPORARY)
					std::cout << "pow " << "@" << f.adress_a[i*4+1]+100 <<"     : @" << f.adress_a[i*4]+100 << "  "<< std::endl ;
				if(f.use_temp[i] == SET_TEMPORARY)
					std::cout << "pow " << "@" << f.adress_a[i*4+1]+100 <<"     : @tmp" << "  "<< std::endl ;
				if(f.use_temp[i] == GET_TEMPORARY_A)
					std::cout << "pow " << "@" << f.adress_a[i*4+1]+100 <<"     : @tmp"  << "  "<< std::endl ;
				if(f.use_temp[i] == GET_TEMPORARY_B)
					std::cout << "pow " << "@tmp" <<"            : @" << f.adress_a[i*4]+100 << "  "<< std::endl ;
				if(f.use_temp[i] == SET_GET_TEMPORARY_A)
					std::cout << "pow " << "@" << f.adress_a[i*4+1]+100 <<"     : @tmp"  << "  "<< std::endl ;
				if(f.use_temp[i] == SET_GET_TEMPORARY_B)
					std::cout << "pow " << "@tmp" <<"     : @tmp"  << "  "<< std::endl ;
				done = true ;
				break ;
			}
		}
	
		if(done)
			continue ;
		
		switch(f.byteCode[i])
		{

			case TOKEN_OPERATION_U:
			{
				std::cout << "sto " << "u" << "@" << f.adress_a[i*4]+100 << "  "<< std::endl ;
				done = true ;
// 				std::cout << "u" << "  " << std::flush ;
				break ;
			}
			case TOKEN_OPERATION_V:
			{
				std::cout << "sto " << "v" << "@" << f.adress_a[i*4]+100 << "  "<< std::endl ;
				done = true ;
				break ;
			}
			case TOKEN_OPERATION_W:
			{
				std::cout << "sto " << "w" << "@" << f.adress_a[i*4]+100 << "  "<< std::endl ;
				done = true ;
				break ;
			}
		}
	
		if(done)
			continue ;
		
		switch(f.byteCode[i])
		{
			case TOKEN_OPERATION_GEO_OPERATION:
			{
				std::cout << "geo " << "@" << f.adress_a[i*4]+100 << "    @" << f.adress_a[i*4+1]+100 <<"    : @" << f.adress_a[i*4+2]+100 << "  "<< std::endl ;
				break ;
			}
			case TOKEN_OPERATION_COS:
			{
				std::cout << "cos " <<  "@" << f.adress_a[i*4]+100 << "            : @"<< f.adress_a[i*4+2]+100 << std::endl ;
				break ;
			}
			case TOKEN_OPERATION_ABS:
			{
				std::cout << "abs " <<  "@" << f.adress_a[i*4]+100 << "            : @"<< f.adress_a[i*4+2]+100 << std::endl ;
				break ;
			}
			case TOKEN_OPERATION_TAN:
			{
				std::cout << "tan " <<  "@" << f.adress_a[i*4]+100 << "            : @"<< f.adress_a[i*4+2]+100 << std::endl ;
				break ;
			}
			case TOKEN_OPERATION_SIN:
			{
				std::cout << "sin " << "@" << f.adress_a[i*4]+100 << "            : @"<< f.adress_a[i*4+2]+100 << std::endl ;
				break ;
			}
			case TOKEN_OPERATION_EXP:
			{
				std::cout << "exp "  << "@" << f.adress_a[i*4]+100 << "            : @"<< f.adress_a[i*4+2]+100 << std::endl ;
				break ;
			}
			case TOKEN_OPERATION_SIGN:
			{
				std::cout << "sgn " << "@" << f.adress_a[i*4]+100 << "            : @"<< f.adress_a[i*4+2]+100 << std::endl ;
				break ;
			}
			case TOKEN_OPERATION_POSITIVITY:
			{
				std::cout << "pos " << "@" << f.adress_a[i*4]+100 << "            : @"<< f.adress_a[i*4+2]+100 << std::endl ;
				break ;
			}
			case TOKEN_OPERATION_NEGATIVITY:
			{
				std::cout << "neg " <<"@" << f.adress_a[i*4]+100 << "            : @"<< f.adress_a[i*4+2]+100 << std::endl ;
				break ;
			}
			case TOKEN_OPERATION_LOG:
			{
				std::cout << "log " <<  "@" << f.adress_a[i*4]+100 << "            : @"<< f.adress_a[i*4+2]+100 << std::endl ;
				break ;
			}
			case TOKEN_OPERATION_COSH:
			{
				std::cout << "cosh " <<  "@" << f.adress_a[i*4]+100 << "            : @"<< f.adress_a[i*4+2]+100 << std::endl ;
				break ;
			}
			case TOKEN_OPERATION_SINH:
			{
				std::cout << "sinh " <<  "@" << f.adress_a[i*4]+100 << "            : @"<< f.adress_a[i*4+2]+100 << std::endl ;
				break ;
			}
			case TOKEN_OPERATION_TANH:
			{
				std::cout << "tanh " <<  "@" << f.adress_a[i*4]+100 << "            : @"<< f.adress_a[i*4+2]+100 << std::endl ;
				break ;
			}
			case TOKEN_OPERATION_SQRT:
			{
				std::cout << "sqrt " << "@" << f.adress_a[i*4]+100 << "            : @"<< f.adress_a[i*4+2]+100 << std::endl ;
				break ;
			}
			case TOKEN_OPERATION_ATAN2:
			{
				std::cout << "sto " << " atan2 " << "@" << f.adress_a[i*4+1]+100 << "    @" << f.adress_a[i*4]+100 << "    : @"<< f.adress_a[i*4+2]+100 << std::endl ;

				break ;
			}
		}
	}
	
	std::cout << "    return @108\n" << std::endl ;
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

	if(f.isDifferentiable(v_))
	{
		return eval( f.d(v_), gp) ;
	}
	
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

double VirtualMachine::ddeval(const Function &f, const Variable v_0, const Variable v_1,  const Point & p, const double eps) 
{
	return ddeval(f, v_0, v_1, p.x,p.y,p.z,p.t, 0., 0., 0., eps) ;
}

double VirtualMachine::ddeval(const Function &f, const Variable v_0, const Variable v_1,  const double x, const double y , const double z, const double t, const double u, const double v, const double w, const double eps) 
{
	if(f.isDifferentiable(v_0) && f.d(v_0).isDifferentiable(v_1))
	{
		
		return eval(f.d(v_0).d(v_1), x, y, z, t, u, v, w) ;
	}
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
//			return ( deval(f, v_1, x-2.*eps, y, z, t, u, v, w)/12. -2./3. * deval(f, v_1, x-eps, y, z, t, u, v, w)- deval(f, v_1, x+2.*eps, y, z, t, u, v, w)/12. +2./3.* deval(f, v_1, x+eps, y, z, t, u, v, w) ) / eps ;
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
double tmp = -(-eval(f, x+2*eps, y, z, t+eps*2, u,v,w) + 8.*eval(f, x+2*eps, y, z, t+eps, u,v,w)-8.*eval(f, x+2*eps, y, z, t-eps, u,v,w)+ eval(f, x+2*eps, y, z, t-eps*2, u,v,w))/(12*eps) ;
tmp += 8* (-eval(f, x+eps, y, z, t+eps*2, u,v,w) + 8.*eval(f, x+eps, y, z, t+eps, u,v,w)-8.*eval(f, x+eps, y, z, t-eps, u,v,w)+ eval(f, x+eps, y, z, t-eps*2, u,v,w))/(12*eps) ;
tmp -= 8* (-eval(f, x-eps, y, z, t+eps*2, u,v,w) + 8.*eval(f, x-eps, y, z, t+eps, u,v,w)-8.*eval(f, x-eps, y, z, t-eps, u,v,w)+ eval(f, x-eps, y, z, t-eps*2, u,v,w))/(12*eps) ;
tmp += (-eval(f, x-2*eps, y, z, t+eps*2, u,v,w) + 8.*eval(f, x-2*eps, y, z, t+eps, u,v,w)-8.*eval(f, x-2*eps, y, z, t-eps, u,v,w)+ eval(f, x-2*eps, y, z, t-eps*2, u,v,w))/(12*eps) ;
					return tmp/(12*eps) ;
					return 0 ;( eval(f, x+eps, y, z, t+eps, u, v, w) - eval(f, x-eps, y, z, t+eps, u, v, w)- eval(f, x+eps, y, z, t-eps, u, v, w) + eval(f, x-eps, y, z, t-eps, u, v, w))/(4.*eps*eps) ;
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
//			return ( deval(f, v_1, x, y-2.*eps, z, t, u, v, w)/12. -2./3. * deval(f, v_1, x, y-eps, z, t, u, v, w)- deval(f, v_1, x, y+2.*eps, z, t, u, v, w)/12. +2./3.* deval(f, v_1, x, y+eps, z, t, u, v, w) ) / eps ;

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
double tmp = -(-eval(f, x, y+2*eps, z, t+eps*2, u,v,w) + 8.*eval(f, x, y+2*eps, z, t+eps, u,v,w)-8.*eval(f, x, y+2*eps, z, t-eps, u,v,w)+ eval(f, x, y+2*eps, z, t-eps*2, u,v,w))/(12*eps) ;
tmp += 8* (-eval(f, x, y+eps, z, t+eps*2, u,v,w) + 8.*eval(f, x, y+eps, z, t+eps, u,v,w)-8.*eval(f, x, y+eps, z, t-eps, u,v,w)+ eval(f, x, y+eps, z, t-eps*2, u,v,w))/(12*eps) ;
tmp -= 8* (-eval(f, x, y-eps, z, t+eps*2, u,v,w) + 8.*eval(f, x, y-eps, z, t+eps, u,v,w)-8.*eval(f, x, y-eps, z, t-eps, u,v,w)+ eval(f, x, y-eps, z, t-eps*2, u,v,w))/(12*eps) ;
tmp += (-eval(f, x, y-eps*2, z, t+eps*2, u,v,w) + 8.*eval(f, x, y-eps*2, z, t+eps, u,v,w)-8.*eval(f, x, y-eps*2, z, t-eps, u,v,w)+ eval(f, x, y-eps*2, z, t-eps*2, u,v,w))/(12*eps) ;
					return tmp/(12*eps) ;
					return ( eval(f, x, y+eps, z, t+eps, u, v, w) - eval(f, x, y-eps, z, t+eps, u, v, w)- eval(f, x, y+eps, z, t-eps, u, v, w) + eval(f, x, y-eps, z, t-eps, u, v, w))/(4.*eps*eps) ;
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
//			return ( deval(f, v_1, x, y, z-2.*eps, t, u, v, w)/12. -2./3. * deval(f, v_1, x, y, z-eps, t, u, v, w)- deval(f, v_1, x, y, z+2.*eps, t, u, v, w)/12. +2./3.* deval(f, v_1, x, y, z+eps, t, u, v, w) ) / eps ;
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
//			return ( deval(f, v_1, x, y, z, t-2.*eps, u, v, w)/12. -2./3. * deval(f, v_1, x, y, z, t-eps, u, v, w)- deval(f, v_1, x, y, z, t+2.*eps, u, v, w)/12. +2./3.* deval(f, v_1, x, y, z, t+eps, u, v, w) ) / eps ;
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
					return( eval(f, x, y+eps, z, t+eps, u, v, w) - eval(f, x, y-eps, z, t+eps, u, v, w)- eval(f, x, y+eps, z, t-eps, u, v, w) + eval(f, x, y-eps, z, t-eps, u, v, w))/(4.*eps*eps) ;
				}
			case ZETA:
				{
					return ( eval(f, x, y, z+eps, t+eps, u, v, w) - eval(f, x, y, z+eps, t-eps, u, v, w)- eval(f, x, y, z-eps, t+eps, u, v, w) + eval(f, x, y, z-eps, t-eps, u, v, w))/(4.*eps*eps) ;
				}
			case TIME_VARIABLE : 
				{
					std::cout << std::setprecision(16) << ( eval(f, x, y, z, t+eps, u, v, w) - 2.*eval(f, x, y, z, t, u, v, w) + eval(f, x, y, z, t-eps, u, v, w))/(4.*eps*eps) << std::endl ;
					return ( eval(f, x, y, z, t+eps, u, v, w) - 2.*eval(f, x, y, z, t, u, v, w) + eval(f, x, y, z, t-eps, u, v, w))/(4.*eps*eps) ;

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
			return ( deval(f, v_1, x, y, z, t, u-2.*eps, v, w)/12. -2./3. * deval(f, v_1, x, y, z, t, u-eps, v, w)- deval(f, v_1, x, y, z, t, u+2.*eps, v, w)/12. +2./3.* deval(f, v_1, x, y, z, t, u+eps, v, w) ) / eps ;

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
			return ( deval(f, v_1, x, y, z, t, u, v-2.*eps, w)/12. -2./3. * deval(f, v_1, x, y, z, t, u, v-eps, w)- deval(f, v_1, x, y, z, t, u, v+2.*eps, w)/12. +2./3.* deval(f, v_1, x, y, z, t, u, v+eps, w) ) / eps ;
			
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
			return ( deval(f, v_1, x, y, z, t, u, v, w-2.*eps)/12. -2./3. * deval(f, v_1, x, y, z, t, u, v, w-eps)- deval(f, v_1, x, y, z, t, u, v, w+2.*eps)/12. +2./3.* deval(f, v_1, x, y, z, t, u, v, w+eps) ) / eps ;
			
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

Vector VirtualMachine::dddeval(const Function&f, const Variable v_0, const Variable v_1, const Variable v_2, const GaussPointArray & gp, const double eps)
{
	Vector ret(gp.gaussPoints.size());
	for(size_t i = 0 ; i < ret.size() ; i++)
	{
		ret[i] = dddeval(f, v_0, v_1, v_2, gp.gaussPoints[i].first.x, gp.gaussPoints[i].first.y, gp.gaussPoints[i].first.z, gp.gaussPoints[i].first.t, 0., 0., 0., eps ) ;
	}	
	return ret ;
}

double VirtualMachine::dddeval(const Function &f, const Variable v_0, const Variable v_1,  const Variable v_2, const double x, const double y , const double z, const double t, const double u, const double v, const double w, const double eps) 
{
	#warning dddeval not fully implemented
	
	if(f.isDifferentiable(v_0) && f.d(v_0).isDifferentiable(v_1) && f.d(v_0).d(v_1).isDifferentiable(v_2) )
	{
		return eval( f.d(v_0).d(v_1).d(v_2), x,y,z,t,u,v,w) ;
	}
	
	double h = 1e-3 ;
	if(v_1 == v_2 && v_1 == TIME_VARIABLE)
	{
		double ret = 0. ;
		switch(v_0)
		{
			case XI:
				ret = eval(f, x+h, y, z, t+h+h, u, v, w) 
				    - eval(f, x-h, y, z, t+h+h, u, v, w) 
				    - 2.*eval(f, x+h, y, z, t, u, v, w) 
				    + 2.*eval(f, x-h, y, z, t, u, v, w) 
				    + eval(f, x+h, y, z, t-h-h, u, v, w)
				    - eval(f, x-h, y, z, t-h-h, u, v, w) ;
				return  ret / (8.*h*h*h) ;
			case ETA:
				ret = eval(f, x, y+h, z, t+h+h, u, v, w) 
				    - eval(f, x, y-h, z, t+h+h, u, v, w) 
				    - 2.*eval(f, x, y+h, z, t, u, v, w) 
				    + 2.*eval(f, x, y-h, z, t, u, v, w) 
				    + eval(f, x, y+h, z, t-h-h, u, v, w) 
				    - eval(f, x, y-h, z, t-h-h, u, v, w) ;
				return ret / (8.*h*h*h) ;		
			case ZETA:
				ret  = eval(f, x, y, z+h, t+h+h, u, v, w) ;
				ret -= eval(f, x, y, z-h, t+h+h, u, v, w) ;
				ret -= 2.*eval(f, x, y, z+h, t, u, v, w) ;
				ret += 2.*eval(f, x, y, z-h, t, u, v, w) ;
				ret += eval(f, x, y, z+h, t-h-h, u, v, w) ;
				ret -= eval(f, x, y, z-h, t-h-h, u, v, w) ;
				return ret / (8.*h*h*h) ;	
			case TIME_VARIABLE:
				return 0. ;
// 				ret = 0. ;
// 				ret += eval(f, x, y, z, t+h+h+h, u, v, w) ;
// 				ret -= 3.*eval(f, x, y, z, t+h, u, v, w) ;
// 				ret += 3.*eval(f, x, y, z, t-h, u, v, w) ;
// 				ret -= eval(f, x, y, z, t-h-h-h, u, v, w) ;
// 				return ret / (8.*h*h*h) ;			
		}
		return ret ;
	}
	return 0 ;
}

double VirtualMachine::deval(const Function &f, const Variable v_,  const double x, const double y , const double z,  const double t, const double u, const double v, const double w, const double eps) 
{
	if(f.isDifferentiable(v_))
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
				double h = eps ;
				volatile double temp = x+h ;
				h = temp - x ;
				return ( eval(f, x-2.*h, y, z, t, u, v, w)/12. -2./3. * eval(f, x-h, y, z, t, u, v, w)- eval(f, x+2.*h, y, z, t, u, v, w)/12. +2./3.* eval(f, x+h, y, z, t, u, v, w) ) / h ;
				return .5*( eval(f, x+h, y, z, t, u, v, w) - eval(f, x-h, y, z, t, u, v, w))/h ;
			}
		case ETA:
			{
				double h = eps ;
				volatile double temp = y+h ;
				h = temp - y ;
				return ( eval(f, x, y-2.*h, z, t, u, v, w)/12. -2./3. * eval(f, x, y-h, z, t, u, v, w)- eval(f, x, y+2.*h, z, t, u, v, w)/12. +2./3.* eval(f, x, y+h, z, t, u, v, w) ) / h ;
				return .5*( eval(f, x, y+h, z, t, u, v, w) - eval(f, x, y-h, z, t, u, v, w))/h ;
			}
		case ZETA:
			{
				double h = eps ;
				volatile double temp = z+h ;
				h = temp - z ;
				return ( eval(f, x, y, z-2.*h, t, u, v, w)/12. -2./3. * eval(f, x, y, z-h, t, u, v, w)- eval(f, x, y, z+2.*h, t, u, v, w)/12. +2./3.* eval(f, x, y, z+h, t, u, v, w) ) / h ;
				return .5*( eval(f, x, y, z+h, t, u, v, w) - eval(f, x, y, z-h, t, u, v, w))/h ;
			}
		case TIME_VARIABLE : 
			{
				double h = eps ;
				volatile double temp = t+h ;
				h = temp - t ;
				return ( eval(f, x, y, z, t-2.*h, u, v, w)/12. -2./3. * eval(f, x, y, z, t-h, u, v, w)- eval(f, x, y, z, t+2.*h, u, v, w)/12. +2./3.* eval(f, x, y, z, t+h, u, v, w) ) / h ;
				return .5*(eval(f, x, y, z, t+h, u, v, w) - eval(f, x, y, z, t-h, u, v, w))/(h) ;
			}
		case U_VARIABLE:
			{
				double h = eps ;
				volatile double temp = u+h ;
				h = temp -u ;
				return .5*( eval(f, x, y, z, t, u+h, v, w) - eval(f, x, y, z, t, u-h, v, w))/(h) ;
			}
		case V_VARIABLE:
			{
				double h = eps ;
				volatile double temp = v+h ;
				h = temp -v ;
				return .5*( eval(f, x, y, z, t, u, v+h, w) - eval(f, x, y, z, t, u, v-h, w))/(h) ;
			}
		case W_VARIABLE:
			{
				double h = eps ;
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

Matrix VirtualMachine::ieval(const FtF & f, const GaussPointArray &gp, const std::valarray<Matrix> & Jinv) 
{
	Matrix ret(Jinv[0].numRows(),Jinv[0].numRows()) ;
	for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
	{
	  	double d = eval(f.first*f.second, gp.gaussPoints[i].first)*gp.gaussPoints[i].second ;
		Matrix m = Jinv[i].transpose() ;
		ret += (Matrix) (m*Jinv[i])*d ;
	}
	
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
		std::valarray< std::pair<Point, double> > gaussPoints = gamma[i].first->getGaussPoints() ;
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

Vector VirtualMachine::ieval(const std::vector<Vector> &f, IntegrableEntity *e)
{
	Vector ret ; ret.resize(f[0].size()) ; ret = 0 ;
	GaussPointArray gp = e->getGaussPoints() ;
	for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
		ret += f[i]*gp.gaussPoints[i].second ;
	
	if(e->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
		ret /= e->area() ;
	else
		ret /= e->volume() ;
	
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

Matrix VirtualMachine::gdeval(const Function &f, IntegrableEntity *e, const std::vector<Variable> & vars, const double x, const double y, const double z, const double t, bool transpose)
{
	Matrix Jinv ; e->getInverseJacobianMatrix(Point(x,y,z,t), Jinv) ;
	return gdeval(f, Jinv, vars, x, y, z,t, transpose) ;
}

Matrix VirtualMachine::geval(const Gradient &f, IntegrableEntity *e, const std::vector<Variable> & vars, const double x, const double  y, const double z, const double t)
{
	return geval(f.f, e,vars, x, y, z,t, f.transpose) ;
}

Matrix VirtualMachine::gdeval(const GradientDot &f, IntegrableEntity *e, const std::vector<Variable> & vars, const double x, const double  y, const double z, const double t)
{
	return gdeval(f.f, e,vars, x, y, z,t, f.transpose) ;
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
				Vector dtau = deval(f, var[2], gp) ;
				
				for(size_t i = 0 ; i < ret.size() ; i++)
				{
					ret[i][0][0] = dxi[i]*m[i][0][0] + deta[i]*m[i][0][1] ;//+ dtau[i]*m[i][0][2] ;
					ret[i][1][1] = dxi[i]*m[i][1][0] + deta[i]*m[i][1][1] ;//+ dtau[i]*m[i][1][2]  ;
					ret[i][2][0] = ret[i][1][1] ;
					ret[i][2][1] = ret[i][0][0] ;
// 					ret[i] *= m[i][2][2] ;
				}
				
				return ret ;
			}
			else
			{
				std::valarray<Matrix> ret(Matrix(2,3), gp.gaussPoints.size()) ;

				Vector dxi = deval(f, var[0], gp) ;
				Vector deta = deval(f, var[1], gp) ;
				Vector dtau = deval(f, var[2], gp) ;
				
				for(size_t i = 0 ; i < ret.size() ; i++)
				{
					ret[i][0][0] = dxi[i]*m[i][0][0] + deta[i]*m[i][0][1] ;//+ dtau[i]*m[i][0][2] ;
					ret[i][1][1] = dxi[i]*m[i][1][0] + deta[i]*m[i][1][1] ;//+ dtau[i]*m[i][1][2] ;
					ret[i][0][2] = ret[i][1][1] ;
					ret[i][1][2] = ret[i][0][0] ;
// 					ret[i] *= m[i][2][2] ;
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
			Vector dteta = deval(f, var[3], gp) ;
			
			for(size_t i = 0 ; i < ret.size() ; i++)
			{
				ret[i][0][0] = dxi[i]*m[i][0][0] + deta[i]*m[i][0][1] + dzeta[i]*m[i][0][2] ;//+ dteta[i]*m[i][0][3];
				ret[i][1][1] = dxi[i]*m[i][1][0] + deta[i]*m[i][1][1] + dzeta[i]*m[i][1][2] ;//+ dteta[i]*m[i][1][3];
				ret[i][2][2] = dxi[i]*m[i][2][0] + deta[i]*m[i][2][1] + dzeta[i]*m[i][2][2] ;//+ dteta[i]*m[i][2][3];
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
			Vector dteta = deval(f, var[3], gp) ;
			
			for(size_t i = 0 ; i < ret.size() ; i++)
			{
				ret[i][0][0] = dxi[i]*m[i][0][0] + deta[i]*m[i][0][1] + dzeta[i]*m[i][0][2];//+ dteta[i]*m[i][0][3];
				ret[i][1][1] = dxi[i]*m[i][1][0] + deta[i]*m[i][1][1] + dzeta[i]*m[i][1][2];//+ dteta[i]*m[i][0][3];
				ret[i][2][2] = dxi[i]*m[i][2][0] + deta[i]*m[i][2][1] + dzeta[i]*m[i][2][2];//+ dteta[i]*m[i][0][3];
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
				Vector dtau = deval(f, var[2], gp) ;
				
				for(size_t i = 0 ; i < ret.size() ; i++)
				{
					ret[i][0][0] = dxi[i]*m[i][0][0] + deta[i]*m[i][0][1] ;//+ dtau[i]*m[i][0][2] ;
					ret[i][1][1] = dxi[i]*m[i][1][0] + deta[i]*m[i][1][1] ;//+ dtau[i]*m[i][1][2]  ;
					ret[i][2][0] = ret[i][1][1] ;
					ret[i][2][1] = ret[i][0][0] ;
				}
				
			}
			else
			{
				
				Vector dxi = deval(f, var[0], gp) ;
				Vector deta = deval(f, var[1], gp) ;
				Vector dtau = deval(f, var[2], gp) ;
				
				for(size_t i = 0 ; i < ret.size() ; i++)
				{
					ret[i][0][0] = dxi[i]*m[i][0][0] + deta[i]*m[i][0][1] ;//+ dtau[i]*m[i][0][2]  ;
					ret[i][1][1] = dxi[i]*m[i][1][0] + deta[i]*m[i][1][1] ;//+ dtau[i]*m[i][1][2]  ;
					ret[i][0][2] = ret[i][1][1] ;
					ret[i][1][2] = ret[i][0][0] ;
				}
//				ret[0].print() ;
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
			Vector dteta = deval(f, var[3], gp) ;
			
			for(size_t i = 0 ; i < ret.size() ; i++)
			{
				ret[i][0][0] = dxi[i]*m[i][0][0] + deta[i]*m[i][0][1] + dzeta[i]*m[i][0][2] ;//+ dteta[i]*m[i][0][3];
				ret[i][1][1] = dxi[i]*m[i][1][0] + deta[i]*m[i][1][1] + dzeta[i]*m[i][1][2] ;//+ dteta[i]*m[i][1][3];
				ret[i][2][2] = dxi[i]*m[i][2][0] + deta[i]*m[i][2][1] + dzeta[i]*m[i][2][2] ;//+ dteta[i]*m[i][2][3];
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
			Vector dteta = deval(f, var[3], gp) ;
			
			for(size_t i = 0 ; i < ret.size() ; i++)
			{
				ret[i][0][0] = dxi[i]*m[i][0][0] + deta[i]*m[i][0][1] + dzeta[i]*m[i][0][2];//+ dteta[i]*m[i][0][3];
				ret[i][1][1] = dxi[i]*m[i][1][0] + deta[i]*m[i][1][1] + dzeta[i]*m[i][1][2];//+ dteta[i]*m[i][0][3];
				ret[i][2][2] = dxi[i]*m[i][2][0] + deta[i]*m[i][2][1] + dzeta[i]*m[i][2][2];//+ dteta[i]*m[i][0][3];
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

std::valarray<Matrix> VirtualMachine::gddeval(const Function &f, const std::valarray<Matrix> & m, const std::vector<Variable> & var, const GaussPointArray & gp, bool transpose )
{
	if(var.size() == 3 )
	{
		if(var[2] == TIME_VARIABLE)
		{
			if(transpose)
			{
				std::valarray<Matrix> ret(Matrix(3,2), gp.gaussPoints.size()) ;
				
/*				Vector dxi = deval( f, var[0], gp) ;
				Vector deta = deval(f, var[1], gp) ;
				Vector dtau = deval(f, var[2], gp) ;
				
				Vector dxidxi = ddeval( f, XI, XI, gp, default_derivation_delta*100.) ;
				Vector detadeta = ddeval( f, ETA, ETA, gp, default_derivation_delta*100.) ;
				Vector dtaudtau = ddeval( f, TIME_VARIABLE, TIME_VARIABLE, gp, default_derivation_delta*100.) ;
				Vector dxideta = ddeval( f, XI, ETA, gp , default_derivation_delta*100.) ;
				Vector detadtau = ddeval( f, ETA, TIME_VARIABLE, gp , default_derivation_delta*100.) ;
				Vector dtaudxi = ddeval( f, TIME_VARIABLE, XI, gp , default_derivation_delta*100.) ;
				
				Vector dxidxidxi(0., gp.gaussPoints.size()) ;
				Vector detadetadeta(0., gp.gaussPoints.size()) ;
				Vector dtaudtaudtau(0., gp.gaussPoints.size()) ;
				Vector dxidxideta(0., gp.gaussPoints.size()) ;
				Vector detadetadxi(0., gp.gaussPoints.size()) ;
				Vector detadetadtau(0., gp.gaussPoints.size()) ;*/
				Vector deta = dddeval(f, ETA, TIME_VARIABLE, TIME_VARIABLE, gp, default_derivation_delta*100.) ;
				Vector dxi = dddeval(f, XI, TIME_VARIABLE, TIME_VARIABLE, gp, default_derivation_delta*100.) ;
				Vector dtau = dddeval(f, TIME_VARIABLE, TIME_VARIABLE, TIME_VARIABLE, gp, default_derivation_delta*100.) ;
/*				Vector dxidxidtau(0., gp.gaussPoints.size()) ;
				Vector dxidetadtau(0., gp.gaussPoints.size()) ;
				
				Vector dx(0.,3) ;
				Vector dxdx(0.,6) ;
				Vector dxdxdx(0.,10) ;
				
				Vector dXdXdX(0., 10) ;*/
				

				for(size_t i = 0 ; i < ret.size() ; i++)
				{
/*					dx[0] = dxi[i] ;
					dx[1] = deta[i] ;
					dx[2] = dtau[i] ;
				  
					dxdx[0] = dxidxi[i] ;
					dxdx[1] = detadeta[i] ;
					dxdx[2] = dtaudtau[i] ;
					dxdx[3] = dxideta[i] ;
					dxdx[4] = detadtau[i] ;
					dxdx[5] = dtaudxi[i] ;
					
					dxdxdx[0] = dxidxidxi[i] ;
					dxdxdx[1] = detadetadeta[i] ;
					dxdxdx[2] = dtaudtaudtau[i] ;
					dxdxdx[3] = dxidxideta[i] ;
					dxdxdx[4] = detadetadxi[i] ;
					dxdxdx[5] = detadetadtau[i] ;
					dxdxdx[6] = dtaudtaudeta[i] ;
					dxdxdx[7] = dtaudtaudxi[i] ;
					dxdxdx[8] = dxidxidtau[i] ;
					dxdxdx[9] = dxidetadtau[i] ;
										
					std::cout << dxi[i] << "\t" << deta[i] << std::endl ;*/
					ret[i][0][0] = dxi[i]*m[i][0][0] + deta[i]*m[i][0][1] + dtau[i]*m[i][0][2] ;
					ret[i][1][1] = dxi[i]*m[i][1][0] + deta[i]*m[i][1][1] + dtau[i]*m[i][1][2] ;
					ret[i][2][0] = ret[i][1][1] ;
					ret[i][2][1] = ret[i][0][0] ;
					ret[i] *= m[i][2][2] ;
					
				}
				
				return ret ;
			}
			else
			{
				std::valarray<Matrix> ret(Matrix(2,3), gp.gaussPoints.size()) ;
				size_t g = gp.gaussPoints.size() ;
/*				
				Vector dxi = deval( f, var[0], gp) ;
				Vector deta = deval(f, var[1], gp) ;
				Vector dtau = deval(f, var[2], gp) ;
				
				Vector dxidxi = ddeval( f, XI, XI, gp, default_derivation_delta*100) ;
				Vector detadeta = ddeval( f, ETA, ETA, gp, default_derivation_delta*100) ;
				Vector dtaudtau = ddeval( f, TIME_VARIABLE, TIME_VARIABLE, gp, default_derivation_delta*100) ;
				Vector dxideta = ddeval( f, XI, ETA, gp , default_derivation_delta*100) ;
				Vector detadtau = ddeval( f, ETA, TIME_VARIABLE, gp , default_derivation_delta*100) ;
				Vector dtaudxi = ddeval( f, XI, TIME_VARIABLE, gp , default_derivation_delta*100) ;
				
				Vector dxidxidxi(0., gp.gaussPoints.size()) ;
				Vector detadetadeta(0., gp.gaussPoints.size()) ;*/
				Vector dtaudtaudtau = dddeval(f, TIME_VARIABLE, TIME_VARIABLE, TIME_VARIABLE, gp, default_derivation_delta*100.) ;
/*				Vector dxidxideta(0., gp.gaussPoints.size()) ;
				Vector detadetadxi(0., gp.gaussPoints.size()) ;
				Vector detadetadtau(0., gp.gaussPoints.size()) ;*/
				Vector dtaudtaudeta = dddeval(f, ETA, TIME_VARIABLE, TIME_VARIABLE, gp, default_derivation_delta*100.) ;
				Vector dtaudtaudxi = dddeval(f, XI, TIME_VARIABLE, TIME_VARIABLE, gp, default_derivation_delta*100.) ;
/*				Vector dxidxidtau(0., gp.gaussPoints.size()) ;
				Vector dxidetadtau(0., gp.gaussPoints.size()) ;
				
				Vector dx(0.,3) ;
				Vector dxdx(0.,6) ;
				Vector dxdxdx(0.,10) ;
				
				Vector dXdXdX(0., 10) ;*/

				for(size_t i = 0 ; i < ret.size() ; i++)
				{
/*					dx[0] = dxi[i] ;
					dx[1] = deta[i] ;
					dx[2] = dtau[i] ;
				  
					dxdx[0] = dxidxi[i] ;
					dxdx[1] = detadeta[i] ;
					dxdx[2] = dtaudtau[i] ;
					dxdx[3] = dxideta[i] ;
					dxdx[4] = detadtau[i] ;
					dxdx[5] = dtaudxi[i] ;
					
					dxdxdx[0] = dxidxidxi[i] ;
					dxdxdx[1] = detadetadeta[i] ;
					dxdxdx[2] = dtaudtaudtau[i] ;
					dxdxdx[3] = dxidxideta[i] ;
					dxdxdx[4] = detadetadxi[i] ;
					dxdxdx[5] = detadetadtau[i] ;
					dxdxdx[6] = dtaudtaudeta[i] ;
					dxdxdx[7] = dtaudtaudxi[i] ;
					dxdxdx[8] = dxidxidtau[i] ;
					dxdxdx[9] = dxidetadtau[i] ;

					dXdXdX = ((Vector) (m[g*5+i]*dxdxdx)) ;//* m[i][2][2]* m[i][2][2] ;
 					dXdXdX += ((Vector) (m[g*4+i]*dxdx)) ;
 					dXdXdX += ((Vector) (m[g*3+i]*dx)) ;//* m[i][2][2] ;

					if(std::abs(dXdXdX[7]) < POINT_TOLERANCE_3D)
						dXdXdX[7] = 0 ;
					if(std::abs(dXdXdX[6]) < POINT_TOLERANCE_3D)
						dXdXdX[6] = 0 ;*/
				  
					ret[i][0][0] = dtaudtaudxi[i]*m[i][0][0] + dtaudtaudeta[i]*m[i][0][1] + dtaudtaudtau[i]*m[i][0][2] ;
					ret[i][1][1] = dtaudtaudxi[i]*m[i][1][0] + dtaudtaudeta[i]*m[i][1][1] + dtaudtaudtau[i]*m[i][1][2]  ;
					ret[i][0][2] = ret[i][1][1] ;
					ret[i][1][2] = ret[i][0][0] ;
 					ret[i] *= m[i][2][2]*m[i][2][2] ;
				}
				
				return ret ;
			}
		}
		else
		{
			if(transpose)
			{
				std::valarray<Matrix> ret(Matrix(6,3), gp.gaussPoints.size()) ;
				
				Vector dxi = dddeval(f, var[0],TIME_VARIABLE,TIME_VARIABLE, gp) ;
				Vector deta = dddeval(f, var[1],TIME_VARIABLE,TIME_VARIABLE, gp) ;
				Vector dzeta = dddeval(f, var[2],TIME_VARIABLE,TIME_VARIABLE,gp) ;
				
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
				
				Vector dxi = dddeval(f, var[0],TIME_VARIABLE,TIME_VARIABLE,gp) ;
				Vector deta = dddeval(f, var[1],TIME_VARIABLE,TIME_VARIABLE,gp) ;
				Vector dzeta = dddeval(f, var[2],TIME_VARIABLE,TIME_VARIABLE, gp) ;
				
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
			
			Vector dxi = dddeval(f, var[0],TIME_VARIABLE,TIME_VARIABLE, gp) ;
			Vector deta = dddeval(f, var[1],TIME_VARIABLE,TIME_VARIABLE, gp) ;
			
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
			
			Vector dxi = dddeval(f, var[0],TIME_VARIABLE,TIME_VARIABLE, gp) ;
			Vector deta = dddeval(f, var[1],TIME_VARIABLE,TIME_VARIABLE, gp) ;
			
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
			
			Vector dxi = dddeval(f, var[0],TIME_VARIABLE,TIME_VARIABLE, gp) ;
			Vector deta = dddeval(f, var[1],TIME_VARIABLE,TIME_VARIABLE, gp) ;
			Vector dzeta = dddeval(f, var[2],TIME_VARIABLE,TIME_VARIABLE, gp) ;
			Vector dteta = dddeval(f, var[3],TIME_VARIABLE,TIME_VARIABLE, gp) ;
			
			for(size_t i = 0 ; i < ret.size() ; i++)
			{
				ret[i][0][0] = dxi[i]*m[i][0][0] + deta[i]*m[i][0][1] + dzeta[i]*m[i][0][2] ;//+ dteta[i]*m[i][0][3];
				ret[i][1][1] = dxi[i]*m[i][1][0] + deta[i]*m[i][1][1] + dzeta[i]*m[i][1][2] ;//+ dteta[i]*m[i][1][3];
				ret[i][2][2] = dxi[i]*m[i][2][0] + deta[i]*m[i][2][1] + dzeta[i]*m[i][2][2] ;//+ dteta[i]*m[i][2][3];
				ret[i][3][1] = ret[i][2][2] ;
				ret[i][3][2] = ret[i][1][1] ;
				ret[i][4][0] = ret[i][2][2] ;
				ret[i][4][2] = ret[i][0][0] ;
				ret[i][5][0] = ret[i][1][1] ;
				ret[i][5][1] = ret[i][0][0] ;
				
				ret[i] *= m[i][3][3] ;
			}
			
			return ret ;
		}
		else
		{
			std::valarray<Matrix> ret(Matrix(3,6), gp.gaussPoints.size()) ;
			
			Vector dxi = dddeval(f, var[0],TIME_VARIABLE,TIME_VARIABLE, gp) ;
			Vector deta = dddeval(f, var[1],TIME_VARIABLE,TIME_VARIABLE, gp) ;
			Vector dzeta = dddeval(f, var[2],TIME_VARIABLE,TIME_VARIABLE, gp) ;
			Vector dteta = dddeval(f, var[3],TIME_VARIABLE,TIME_VARIABLE, gp) ;

			for(size_t i = 0 ; i < ret.size() ; i++)
			{
				ret[i][0][0] = dxi[i]*m[i][0][0] + deta[i]*m[i][0][1] + dzeta[i]*m[i][0][2];//+ dteta[i]*m[i][0][3];
				ret[i][1][1] = dxi[i]*m[i][1][0] + deta[i]*m[i][1][1] + dzeta[i]*m[i][1][2];//+ dteta[i]*m[i][0][3];
				ret[i][2][2] = dxi[i]*m[i][2][0] + deta[i]*m[i][2][1] + dzeta[i]*m[i][2][2];//+ dteta[i]*m[i][0][3];
				ret[i][1][3] = ret[i][2][2] ;
				ret[i][2][3] = ret[i][1][1] ;
				ret[i][0][4] = ret[i][2][2] ;
				ret[i][2][4] = ret[i][0][0] ;
				ret[i][0][5] = ret[i][1][1] ;
				ret[i][1][5] = ret[i][0][0] ;

				ret[i] *= m[i][3][3] ;
			}
			
			return ret ;
		}
	}
	
	return std::valarray<Matrix>(Matrix(0,0), gp.gaussPoints.size()) ; ;
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
				
//				size_t g = gp.gaussPoints.size() ;

//				Vector dxidxi = ddeval(f, var[0], var[0], gp,default_derivation_delta*100./**0.01*/) ;
//				Vector detadeta = ddeval(f, var[1], var[1], gp,default_derivation_delta*100./**0.01*/) ;
				Vector dtau = ddeval(f, var[2], var[2], gp,default_derivation_delta*100./**0.01*/) ;
				Vector dxi = ddeval(f, var[0], var[2], gp,default_derivation_delta*100./**0.01*/) ;
//				Vector detadtau = ddeval(f, var[1], var[2], gp,default_derivation_delta*100./**0.01*/) ;
				Vector deta = ddeval(f, var[1], var[2], gp,default_derivation_delta*100./**0.01*/) ;

//				Vector dxi = deval(f, var[0], gp) ;
//				Vector deta = deval(f, var[1], gp) ;
//				Vector dtau = deval(f, var[2], gp) ;
				
				for(size_t i = 0 ; i < ret.size() ; i++)
				{
				ret[i][0][0] = dxi[i]*m[i][0][0] + deta[i]*m[i][0][1] ;//+ dtau[i]*m[i][0][2] ;
				ret[i][1][1] = dxi[i]*m[i][1][0] + deta[i]*m[i][1][1] ;//+ dtau[i]*m[i][1][2]  ;
				ret[i][2][0] = ret[i][1][1] ;
				ret[i][2][1] = ret[i][0][0] ;
				ret[i] *= m[i][2][2] ;
/*					ret[i][0][0] = dxi[i]*m[g+i][5][0] + deta[i]*m[g+i][5][1] + dtau[i]*m[g+i][5][2] 
						+ dxidxi[i]*m[g*2+i][5][0] + detadeta[i]*m[g*2+i][5][1] + dtaudtau[i]*m[g*2+i][5][2] 
						+ dxideta[i]*m[g*2+i][5][3] + detadtau[i]*m[g*2+i][5][4] + dxidtau[i]*m[g*2+i][5][5] ;
					ret[i][1][1] = dxi[i]*m[g+i][4][0] + deta[i]*m[g+i][4][1] + dtau[i]*m[g+i][4][2] 
						+ dxidxi[i]*m[g*2+i][4][0] + detadeta[i]*m[g*2+i][4][1] + dtaudtau[i]*m[g*2+i][4][2] 
						+ dxideta[i]*m[g*2+i][4][3] + detadtau[i]*m[g*2+i][4][4] + dxidtau[i]*m[g*2+i][4][5] ;
					ret[i][2][0] = ret[i][1][1] ;
					ret[i][2][1] = ret[i][0][0] ;*/
				}
				
				return ret ;
			}
			else
			{
				std::valarray<Matrix> ret(Matrix(2,3), gp.gaussPoints.size()) ;
				
//				size_t g = gp.gaussPoints.size() ;

//				Vector dxidxi = ddeval(f, var[0], var[0], gp,default_derivation_delta*100./**0.01*/) ;
//				Vector detadeta = ddeval(f, var[1], var[1], gp,default_derivation_delta*100./**0.01*/) ;
				Vector dtau = ddeval(f, var[2], var[2], gp,default_derivation_delta*100./**0.01*/) ;
				Vector deta = ddeval(f, var[1], var[2], gp,default_derivation_delta*100./**0.01*/) ;
//				Vector detadtau = ddeval(f, var[1], var[2], gp,default_derivation_delta*100./**0.01*/) ;
				Vector dxi = ddeval(f, var[0], var[2], gp,default_derivation_delta*100./**0.01*/) ;

//				Vector dxi = deval(f, var[0], gp) ;
//				Vector deta = deval(f, var[1], gp) ;
//				Vector dtau = deval(f, var[2], gp) ;
				
				for(size_t i = 0 ; i < ret.size() ; i++)
				{
				ret[i][0][0] = dxi[i]*m[i][0][0] + deta[i]*m[i][0][1] ;//+ dtau[i]*m[i][0][2] ;
				ret[i][1][1] = dxi[i]*m[i][1][0] + deta[i]*m[i][1][1] ;//+ dtau[i]*m[i][1][2]  ;
				ret[i][0][2] = ret[i][1][1] ;
				ret[i][1][2] = ret[i][0][0] ;
				ret[i] *= m[i][2][2] ;
/*					ret[i][0][0] = dxi[i]*m[g+i][5][0] + deta[i]*m[g+i][5][1] + dtau[i]*m[g+i][5][2] 
						+ dxidxi[i]*m[g*2+i][5][0] + detadeta[i]*m[g*2+i][5][1] + dtaudtau[i]*m[g*2+i][5][2] 
						+ dxideta[i]*m[g*2+i][5][3] + detadtau[i]*m[g*2+i][5][4] + dxidtau[i]*m[g*2+i][5][5] ;
					ret[i][1][1] = dxi[i]*m[g+i][4][0] + deta[i]*m[g+i][4][1] + dtau[i]*m[g+i][4][2] 
						+ dxidxi[i]*m[g*2+i][4][0] + detadeta[i]*m[g*2+i][4][1] + dtaudtau[i]*m[g*2+i][4][2] 
						+ dxideta[i]*m[g*2+i][4][3] + detadtau[i]*m[g*2+i][4][4] + dxidtau[i]*m[g*2+i][4][5] ;
					ret[i][2][0] = ret[i][1][1] ;
					ret[i][2][1] = ret[i][0][0] ;*/
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
				ret[i][0][0] = dxi[i]*m[i][0][0] + deta[i]*m[i][0][1] + dzeta[i]*m[i][0][2] ;// + dteta[i]*m[i][0][3];
				ret[i][1][1] = dxi[i]*m[i][1][0] + deta[i]*m[i][1][1] + dzeta[i]*m[i][1][2] ;// + dteta[i]*m[i][1][3];
				ret[i][2][2] = dxi[i]*m[i][2][0] + deta[i]*m[i][2][1] + dzeta[i]*m[i][2][2] ;// + dteta[i]*m[i][2][3];
				ret[i][3][1] = ret[i][2][2] ;
				ret[i][3][2] = ret[i][1][1] ;
				ret[i][4][0] = ret[i][2][2] ;
				ret[i][4][2] = ret[i][0][0] ;
				ret[i][5][0] = ret[i][1][1] ;
				ret[i][5][1] = ret[i][0][0] ;
				
				ret[i] *= m[i][3][3] ;				
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
//			std::cout << dxi[0] << std::endl ;
			
			for(size_t i = 0 ; i < ret.size() ; i++)
			{
				ret[i][0][0] = dxi[i]*m[i][0][0] + deta[i]*m[i][0][1] + dzeta[i]*m[i][0][2];//+ dteta[i]*m[i][0][3];
				ret[i][1][1] = dxi[i]*m[i][1][0] + deta[i]*m[i][1][1] + dzeta[i]*m[i][1][2];//+ dteta[i]*m[i][0][3];
				ret[i][2][2] = dxi[i]*m[i][2][0] + deta[i]*m[i][2][1] + dzeta[i]*m[i][2][2];//+ dteta[i]*m[i][0][3];
				ret[i][1][3] = ret[i][2][2] ;
				ret[i][2][3] = ret[i][1][1] ;
				ret[i][0][4] = ret[i][2][2] ;
				ret[i][2][4] = ret[i][0][0] ;
				ret[i][0][5] = ret[i][1][1] ;
				ret[i][1][5] = ret[i][0][0] ;

				ret[i] *= m[i][3][3] ;
			  
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
				double dtau = deval(f, var[2], x,y,z,t) ;

				Matrix ret(3,2) ;
				ret[0][0] = dxi*m[0][0] + deta*m[0][1] ;//+ dtau*m[0][2] ;
				ret[1][1] = dxi*m[1][0] + deta*m[1][1] ;//+ dtau*m[1][2]  ;
				ret[2][0] = ret[1][1] ;
				ret[2][1] = ret[0][0] ;
				
				return ret ;
			}
			else
			{
				Matrix ret(2,3) ;
				double dxi = deval(f, var[0], x,y,z,t) ;
				double deta = deval(f, var[1], x,y,z,t) ;
				double dtau = deval(f, var[2], x,y,z,t) ;

				ret[0][0] = dxi*m[0][0] + deta*m[0][1] ;//+ dtau*m[0][2]  ;
				ret[1][1] = dxi*m[1][0] + deta*m[1][1] ;//+ dtau*m[1][2]  ;
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
			double dxi = deval(f, var[0], x,y,z,t) ;
			double deta = deval(f, var[1], x,y,z,t) ;
			ret[0][0] = dxi*m[0][0] + deta*m[0][1] ;
			ret[1][1] = dxi*m[1][0] + deta*m[1][1] ;
			ret[2][0] = ret[1][1] ;
			ret[2][1] = ret[0][0] ;
			
			return ret ;
		}
		else
		{
			Matrix ret(2,3) ;
			double dxi = deval(f, var[0], x,y,z,t) ;
			double deta = deval(f, var[1], x,y,z,t) ;
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
//			double dteta = deval(f, var[3], x,y,z,t) ;
			Matrix ret(6,3) ;
			ret[0][0] = dxi*m[0][0] + deta*m[0][1] + dzeta*m[0][2] ;//+ dteta*m[0][3];
			ret[1][1] = dxi*m[1][0] + deta*m[1][1] + dzeta*m[1][2] ;//+ dteta*m[1][3];
			ret[2][2] = dxi*m[2][0] + deta*m[2][1] + dzeta*m[2][2] ;//+ dteta*m[2][3];
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
//			double dteta = deval(f, var[3], x,y,z,t) ;
			ret[0][0] = dxi*m[0][0] + deta*m[0][1] + dzeta*m[0][2];//+ dteta*m[0][3];
			ret[1][1] = dxi*m[1][0] + deta*m[1][1] + dzeta*m[1][2];//+ dteta*m[0][3];
			ret[2][2] = dxi*m[2][0] + deta*m[2][1] + dzeta*m[2][2];//+ dteta*m[0][3];
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
				double dtau = deval(f, var[2], x,y,z,t) ;
				if(ret.isNull() || ret.numRows() !=3 ||ret.numCols() != 2)
					ret.resize(3,2) ;
				ret[0][0] = dxi*m[0][0] + deta*m[0][1] ;//+ dtau *m[0][2];
 				ret[1][1] = dxi*m[1][0] + deta*m[1][1] ;//+ dtau *m[1][2] ;
				ret[2][0] = ret[1][1] ;
				ret[2][1] = ret[0][0] ;
				
			}
			else
			{
				double dxi = deval(f, var[0], x,y,z,t) ;
				double deta = deval(f, var[1], x,y,z,t) ;
				double dtau = deval(f, var[2], x,y,z,t) ;
				if(ret.isNull() || ret.numRows() != 2 ||ret.numCols() != 3)
					ret.resize(2,3) ;
				ret[0][0] = dxi*m[0][0] + deta*m[0][1] ;//+ dtau *m[0][2] ;
				ret[1][1] = dxi*m[1][0] + deta*m[1][1] ;//+ dtau *m[1][2] ;
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
			double dxi = deval(f, var[0], x,y,z,t) ;
			double deta = deval(f, var[1], x,y,z,t) ;
			ret[0][0] = dxi*m[0][0] + deta*m[0][1] ;
			ret[1][1] = dxi*m[1][0] + deta*m[1][1] ;
			ret[2][0] = ret[1][1] ;
			ret[2][1] = ret[0][0] ;
		
			
		}
		else
		{
			if(ret.isNull()|| ret.numRows() != 2 ||ret.numCols() != 3)
				ret.resize(2,3) ;
			double dxi = deval(f, var[0], x,y,z,t) ;
			double deta = deval(f, var[1], x,y,z,t) ;
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
			double dteta = deval(f, var[3], x,y,z,t) ;
			ret[0][0] = dxi*m[0][0] + deta*m[0][1] + dzeta*m[0][2] ;//+ dteta*m[0][3];
			ret[1][1] = dxi*m[1][0] + deta*m[1][1] + dzeta*m[1][2] ;//+ dteta*m[1][3];
			ret[2][2] = dxi*m[2][0] + deta*m[2][1] + dzeta*m[2][2] ;//+ dteta*m[2][3];
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
			double dxi = deval(f, var[0], x,y,z,t) ;
			double deta = deval(f, var[1], x,y,z,t) ;
			double dzeta = deval(f, var[2], x,y,z,t) ;
			double dteta = deval(f, var[3], x,y,z,t) ;
			ret[0][0] = dxi*m[0][0] + deta*m[0][1] + dzeta*m[0][2];//+ dteta*m[0][3];
			ret[1][1] = dxi*m[1][0] + deta*m[1][1] + dzeta*m[1][2];//+ dteta*m[0][3];
			ret[2][2] = dxi*m[2][0] + deta*m[2][1] + dzeta*m[2][2];//+ dteta*m[0][3];
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
				double dxi = ddeval(f, var[0],TIME_VARIABLE, x,y,z,t,0,0,0,100.*default_derivation_delta) ;
				double deta = ddeval(f, var[1],TIME_VARIABLE, x,y,z,t,0,0,0,100.*default_derivation_delta) ;
				double dtau = ddeval(f, var[2], TIME_VARIABLE, x,y,z,t,0,0,0,100.*default_derivation_delta) ;
// 				std::cout << dxi*m[1][0] << "\t" << deta*m[1][1] << "\t" << dtau*m[1][2] << std::endl ;
//  				std::cout << dxi << "\t" << deta << "\t" << dtau << std::endl ;

				if(std::abs(dxi) < POINT_TOLERANCE_3D)
					dxi = 0 ;
				if(std::abs(deta) < POINT_TOLERANCE_3D)
					deta = 0 ;
				if(std::abs(dtau) < POINT_TOLERANCE_3D)
					dtau = 0 ;
				
				Matrix ret(3,2) ;
				ret[0][0] = dxi*m[0][0] + deta*m[0][1] ;//+ dtau*m[0][2];
				ret[1][1] = dxi*m[1][0] + deta*m[1][1] ;//+ dtau*m[1][2] ;
				ret[2][0] = ret[1][1] ;
				ret[2][1] = ret[0][0] ;
				ret *= m[2][2] ;
				
				return ret ;
			}
			else
			{
				Matrix ret(2,3) ;
				double dxi = ddeval(f, var[0],TIME_VARIABLE, x,y,z,t,0,0,0,100.*default_derivation_delta) ;
				double deta = ddeval(f, var[1],TIME_VARIABLE, x,y,z,t,0,0,0,100.*default_derivation_delta) ;
				double dtau = ddeval(f, var[2], TIME_VARIABLE, x,y,z,t,0,0,0,100.*default_derivation_delta) ;
//  				std::cout << dxi << "\t" << deta << "\t" << dtau << std::endl ;
				
				ret[0][0] = dxi*m[0][0] + deta*m[0][1] ;//+ dtau*m[0][2];
				ret[1][1] = dxi*m[1][0] + deta*m[1][1] ;//+ dtau*m[1][2];
				ret[0][2] = ret[1][1] ;
				ret[1][2] = ret[0][0] ;
				ret *= m[2][2] ;
				
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
			ret *= m[3][3] ;
			
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
			ret *= m[3][3] ;
			
			return ret ;
		}
	}
	return Matrix(0,0) ;
}

void VirtualMachine::gdeval(const Function &f, const Matrix & m, const std::vector<Variable> & var, const double x, const double y, const double  z, const double  t, bool transpose, Matrix & ret )
{
	if(var.size() == 3 )
	{
		if(var[2] == TIME_VARIABLE)
		{
			if(transpose)
			{
				double dxi = ddeval(f, var[0],TIME_VARIABLE, x,y,z,t, default_derivation_delta*100.) ;
				double deta = ddeval(f, var[1],TIME_VARIABLE, x,y,z,t, default_derivation_delta*100.) ;
				double dtau = ddeval(f, var[2], TIME_VARIABLE, x,y,z,t, default_derivation_delta*100.) ;
				
				if(ret.isNull() || ret.numRows() !=3 ||ret.numCols() != 2)
					ret.resize(3,2) ;
				
				ret[0][0] = dxi*m[0][0] + deta*m[0][1] ;// + dtau*m[0][2] ;
				ret[1][1] = dxi*m[1][0] + deta*m[1][1] ;//+ dtau*m[1][2]  ;
				ret[2][0] = ret[1][1] ;
				ret[2][1] = ret[0][0] ;
				ret *= m[2][2] ;
				
				return ;
			}
			else
			{
				double dxi = ddeval(f, var[0],TIME_VARIABLE, x,y,z, default_derivation_delta*100.) ;
				double deta = ddeval(f, var[1],TIME_VARIABLE, x,y,z, default_derivation_delta*100.) ;
				double dtau = ddeval(f, var[2], TIME_VARIABLE, x,y,z,t, default_derivation_delta*100.) ;
				if(ret.isNull() || ret.numRows() !=2 ||ret.numCols() != 3)
					ret.resize(2,3) ;
				
				ret[0][0] = dxi*m[0][0] + deta*m[0][1] ;// + dtau*m[0][2] ;;
				ret[1][1] = dxi*m[1][0] + deta*m[1][1] ;// + dtau*m[1][2] ;;
				ret[0][2] = ret[1][1] ;
				ret[1][2] = ret[0][0] ;
				ret *= m[2][2] ;
				
				return ;
			}
		}
		else
		{
			if(transpose)
			{
				double dxi = ddeval(f, var[0],TIME_VARIABLE, x,y,z,t) ;
				double deta = ddeval(f, var[1],TIME_VARIABLE, x,y,z,t) ;
				double dzeta = ddeval(f, var[2],TIME_VARIABLE, x,y,z,t) ;
				if(ret.isNull() || ret.numRows() !=6 ||ret.numCols() != 3)
					ret.resize(3,6) ;
				
				ret[0][0] = dxi*m[0][0] + deta*m[0][1] + dzeta*m[0][2] ;
				ret[1][1] = dxi*m[1][0] + deta*m[1][1] + dzeta*m[1][2] ;
				ret[2][2] = dxi*m[2][0] + deta*m[2][1] + dzeta*m[2][2] ;
				ret[3][1] = ret[2][2] ;
				ret[3][2] = ret[1][1] ;
				ret[4][0] = ret[2][2] ;
				ret[4][2] = ret[0][0] ;
				ret[5][0] = ret[1][1] ;
				ret[5][1] = ret[0][0] ;
				
				return ;
			}
			else
			{
				double dxi = ddeval(f, var[0],TIME_VARIABLE, x,y,z,t) ;
				double deta = ddeval(f, var[1],TIME_VARIABLE, x,y,z,t) ;
				double dzeta = ddeval(f, var[2],TIME_VARIABLE, x,y,z,t) ;
				if(ret.isNull() || ret.numRows() !=6 ||ret.numCols() != 3)
					ret.resize(6,3) ;
				ret[0][0] = dxi*m[0][0] + deta*m[0][1] + dzeta*m[0][2];
				ret[1][1] = dxi*m[1][0] + deta*m[1][1] + dzeta*m[1][2];
				ret[2][2] = dxi*m[2][0] + deta*m[2][1] + dzeta*m[2][2];
				ret[1][3] = ret[2][2] ;
				ret[2][3] = ret[1][1] ;
				ret[0][4] = ret[2][2] ;
				ret[2][4] = ret[0][0] ;
				ret[0][5] = ret[1][1] ;
				ret[1][5] = ret[0][0] ;
				
				return ;
			}
		}
	}
	else if (var.size() == 2)
	{
		if(transpose)
		{
			
			double dxi = ddeval(f, var[0],TIME_VARIABLE, x,y,z) ;
			double deta = ddeval(f, var[1],TIME_VARIABLE, x,y,z) ;
				if(ret.isNull() || ret.numRows() !=2 ||ret.numCols() != 2)
					ret.resize(3,2) ;
			ret[0][0] = dxi*m[0][0] + deta*m[0][1] ;
			ret[1][1] = dxi*m[1][0] + deta*m[1][1] ;
			ret[2][0] = ret[1][1] ;
			ret[2][1] = ret[0][0] ;
			
			return ;
		}
		else
		{
			double dxi = ddeval(f, var[0],TIME_VARIABLE, x,y,z) ;
			double deta = ddeval(f, var[1],TIME_VARIABLE, x,y,z) ;
				if(ret.isNull() || ret.numRows() !=2 ||ret.numCols() != 2)
					ret.resize(2,3) ;
			ret[0][0] = dxi*m[0][0] + deta*m[0][1] ;
			ret[1][1] = dxi*m[1][0] + deta*m[1][1] ;
			ret[0][2] = ret[1][1] ;
			ret[1][2] = ret[0][0] ;
			
			
			return ;
		}
	}
	else if(var.size() == 4 )
	{
		if(transpose)
		{
			double dxi = ddeval(f, var[0],TIME_VARIABLE, x,y,z,t) ;
			double deta = ddeval(f, var[1],TIME_VARIABLE, x,y,z,t) ;
			double dzeta = ddeval(f, var[2],TIME_VARIABLE, x,y,z,t) ;
			double dteta = ddeval(f, var[3],TIME_VARIABLE, x,y,z,t) ;
				if(ret.isNull() || ret.numRows() !=6 ||ret.numCols() != 3)
					ret.resize(6,3) ;
			ret[0][0] = dxi*m[0][0] + deta*m[0][1] + dzeta*m[0][2] ;//+ dteta*m[0][3];
			ret[1][1] = dxi*m[1][0] + deta*m[1][1] + dzeta*m[1][2] ;//+ dteta*m[1][3];
			ret[2][2] = dxi*m[2][0] + deta*m[2][1] + dzeta*m[2][2] ;//+ dteta*m[2][3];
			ret[3][1] = ret[2][2] ;
			ret[3][2] = ret[1][1] ;
			ret[4][0] = ret[2][2] ;
			ret[4][2] = ret[0][0] ;
			ret[5][0] = ret[1][1] ;
			ret[5][1] = ret[0][0] ;
			ret *= m[3][3] ;
			
			return ;
		}
		else
		{
			double dxi = ddeval(f, var[0],TIME_VARIABLE, x,y,z,t) ;
			double deta = ddeval(f, var[1],TIME_VARIABLE, x,y,z,t) ;
			double dzeta = ddeval(f, var[2],TIME_VARIABLE, x,y,z,t) ;
			double dteta = ddeval(f, var[3],TIME_VARIABLE, x,y,z,t) ;
				if(ret.isNull() || ret.numRows() !=3 ||ret.numCols() != 6)
					ret.resize(3,6) ;
			ret[0][0] = dxi*m[0][0] + deta*m[0][1] + dzeta*m[0][2];//+ dteta*m[0][3];
			ret[1][1] = dxi*m[1][0] + deta*m[1][1] + dzeta*m[1][2];//+ dteta*m[0][3];
			ret[2][2] = dxi*m[2][0] + deta*m[2][1] + dzeta*m[2][2];//+ dteta*m[0][3];
			ret[1][3] = ret[2][2] ;
			ret[2][3] = ret[1][1] ;
			ret[0][4] = ret[2][2] ;
			ret[2][4] = ret[0][0] ;
			ret[0][5] = ret[1][1] ;
			ret[1][5] = ret[0][0] ;
			ret *= m[3][3] ;
			
		}
	}
	return ;
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

Matrix VirtualMachine::ieval(const GtML &f, IntegrableEntity *e, const std::vector<Variable> & var)
{

	GaussPointArray gp = e->getGaussPoints() ;
	Matrix B (geval(f.first.f, e,var, gp.gaussPoints[0].first.x, gp.gaussPoints[0].first.y, gp.gaussPoints[0].first.z, gp.gaussPoints[0].first.t, f.first.transpose)*gp.gaussPoints[0].second) ;
	Matrix ret = B*f.second[0] ;
	for(size_t i = 1 ; i < gp.gaussPoints.size() ; i++)
	{
		B = geval(f.first.f, e,var, gp.gaussPoints[i].first.x, gp.gaussPoints[i].first.y, gp.gaussPoints[i].first.z, gp.gaussPoints[i].first.t, f.first.transpose)*gp.gaussPoints[i].second ;
		ret += B*f.second[i] ;
	}
	
	return ret ;
}

Vector VirtualMachine::ieval(const GtVL &f, IntegrableEntity *e, const std::vector<Variable> & var)
{

	GaussPointArray gp = e->getGaussPoints() ;
	Matrix B = geval(f.first.f, e,var, gp.gaussPoints[0].first.x, gp.gaussPoints[0].first.y, gp.gaussPoints[0].first.z, gp.gaussPoints[0].first.t, f.first.transpose) ;
	B *= gp.gaussPoints[0].second ;
//	B.print() ;
	Vector ret = B*f.second[0] ;
//	std::cout << gp.gaussPoints.size() << "\t" << f.second.size() << std::endl ;
	for(size_t i = 1 ; i < gp.gaussPoints.size() ; i++)
	{
		B = geval(f.first.f, e,var, gp.gaussPoints[i].first.x, gp.gaussPoints[i].first.y, gp.gaussPoints[i].first.z, gp.gaussPoints[i].first.t, f.first.transpose)*gp.gaussPoints[i].second ;
		ret += B*f.second[i] ;
	}
	return ret ;
}

Vector VirtualMachine::ieval(const GtV &f, IntegrableEntity *e, const std::vector<Variable> & var)
{

	GaussPointArray gp = e->getGaussPoints() ;
	Matrix B = geval(f.first.f, e,var, gp.gaussPoints[0].first.x, gp.gaussPoints[0].first.y, gp.gaussPoints[0].first.z, gp.gaussPoints[0].first.t, f.first.transpose) ; 
	B *= gp.gaussPoints[0].second ;
	Vector ret = B*f.second ;
	for(size_t i = 1 ; i < gp.gaussPoints.size() ; i++)
	{
		B = geval(f.first.f, e,var, gp.gaussPoints[i].first.x, gp.gaussPoints[i].first.y, gp.gaussPoints[i].first.z, gp.gaussPoints[i].first.t, f.first.transpose) ;
		B *= gp.gaussPoints[i].second ;
		ret += B*f.second ;
	}
	return ret ;
}

Vector VirtualMachine::ieval(const GDtV &f, IntegrableEntity *e, const std::vector<Variable> & var)
{

	GaussPointArray gp = e->getGaussPoints() ;
	Matrix B = gdeval(f.first.f, e,var, gp.gaussPoints[0].first.x, gp.gaussPoints[0].first.y, gp.gaussPoints[0].first.z, gp.gaussPoints[0].first.t, f.first.transpose) ;
	B *= gp.gaussPoints[0].second ;
	Vector ret = B*f.second ;
	for(size_t i = 1 ; i < gp.gaussPoints.size() ; i++)
	{
		B = gdeval(f.first.f, e,var, gp.gaussPoints[i].first.x, gp.gaussPoints[i].first.y, gp.gaussPoints[i].first.z, gp.gaussPoints[i].first.t, f.first.transpose)*gp.gaussPoints[i].second ;
		ret += B*f.second ;
	}
	return ret ;
}

Vector VirtualMachine::ieval(const GDtVL &f, IntegrableEntity *e, const std::vector<Variable> & var)
{

	GaussPointArray gp = e->getGaussPoints() ;
	Matrix B = gdeval(f.first.f, e,var, gp.gaussPoints[0].first.x, gp.gaussPoints[0].first.y, gp.gaussPoints[0].first.z, gp.gaussPoints[0].first.t, f.first.transpose) ;
	B *= gp.gaussPoints[0].second ;
	Vector ret = B*f.second[0] ;
	for(size_t i = 1 ; i < gp.gaussPoints.size() ; i++)
	{
		B = gdeval(f.first.f, e,var, gp.gaussPoints[i].first.x, gp.gaussPoints[i].first.y, gp.gaussPoints[i].first.z, gp.gaussPoints[i].first.t, f.first.transpose)*gp.gaussPoints[i].second ;
		ret += B*f.second[i] ;
	}
	return ret ;
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

std::valarray<Matrix> VirtualMachine::ievalDecomposed(const GtMtG &f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & var)
{
	std::valarray<Matrix> ret(Jinv.size()) ;
    
	geval(f.third.f, Jinv[0], var, gp.gaussPoints[0].first, f.third.transpose,B_) ;
	geval(f.first.f, Jinv[0], var, gp.gaussPoints[0].first, f.first.transpose,B) ;
	
	Matrix ret_(B*f.second*B_) ;
	ret[0] = ret_*gp.gaussPoints[0].second ;
	
	for(size_t i = 1 ; i  < gp.gaussPoints.size() ; i++)
	{
		geval(f.third.f, Jinv[i], var, gp.gaussPoints[i].first, f.third.transpose, B_) ;
		geval(f.first.f, Jinv[i], var, gp.gaussPoints[i].first, f.first.transpose, B) ;
		ret_ = B*f.second*B_ ;
		ret[i] = ret_*gp.gaussPoints[i].second ;
//		matrix_matrix_matrix_multiply_and_add(B, f.second, B_, gp.gaussPoints[i].second, ret) ;
	}

	return ret ;
}

std::valarray<Matrix> VirtualMachine::ievalDecomposed(const GtMtGD &f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & var)
{
	
  
  
  
	std::valarray<Matrix> ret(Jinv.size()) ;
    
	Matrix B_ = gdeval(f.third.f, Jinv[0], var, gp.gaussPoints[0].first.x, gp.gaussPoints[0].first.y, gp.gaussPoints[0].first.z, gp.gaussPoints[0].first.t, f.third.transpose) ;
	Matrix B = geval(f.first.f, Jinv[0], var, gp.gaussPoints[0].first.x, gp.gaussPoints[0].first.y, gp.gaussPoints[0].first.z, gp.gaussPoints[0].first.t, f.first.transpose) ;
	
	Matrix ret_(B*f.second*B_) ;
	ret[0] = ret_*gp.gaussPoints[0].second ;
	
	for(size_t i = 1 ; i  < gp.gaussPoints.size() ; i++)
	{
		B_ = gdeval(f.third.f, Jinv[i], var, gp.gaussPoints[i].first.x, gp.gaussPoints[i].first.y, gp.gaussPoints[i].first.z, gp.gaussPoints[i].first.t, f.third.transpose) ;
		B = geval(f.first.f, Jinv[i], var, gp.gaussPoints[i].first.x, gp.gaussPoints[i].first.y, gp.gaussPoints[i].first.z, gp.gaussPoints[i].first.t, f.first.transpose) ;
		ret_ = B*f.second*B_ ;
		ret[i] = ret_*gp.gaussPoints[i].second ;
//		matrix_matrix_matrix_multiply_and_add(B, f.second, B_, gp.gaussPoints[i].second, ret) ;
	}
	
	


	return ret ;
}

std::valarray<Matrix> VirtualMachine::ievalDecomposedDebug(const GtMtGD &f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & var)
{
	
  
  
  
	std::valarray<Matrix> ret(Jinv.size()) ;
    
	Matrix B_ = gdeval(f.third.f, Jinv[0], var, gp.gaussPoints[0].first.x, gp.gaussPoints[0].first.y, gp.gaussPoints[0].first.z, gp.gaussPoints[0].first.t, f.third.transpose) ;
	Matrix B = geval(f.first.f, Jinv[0], var, gp.gaussPoints[0].first.x, gp.gaussPoints[0].first.y, gp.gaussPoints[0].first.z, gp.gaussPoints[0].first.t, f.first.transpose) ;
	
	Matrix ret_(B*f.second*B_) ;
	ret[0] = ret_*gp.gaussPoints[0].second ;
	
	if(std::abs(ret[0][0][0]) > 1)
	{
		B.print() ;
		f.second.print() ;
		B_.print() ;
		if(f.third.transpose)
			std::cout << "transpose-third-" ;
		std::cout << "miaou " << 0 << std::endl ;
	}
	
	for(size_t i = 1 ; i  < gp.gaussPoints.size() ; i++)
	{
		B_ = gdeval(f.third.f, Jinv[i], var, gp.gaussPoints[i].first.x, gp.gaussPoints[i].first.y, gp.gaussPoints[i].first.z, gp.gaussPoints[i].first.t, f.third.transpose) ;
		B = geval(f.first.f, Jinv[i], var, gp.gaussPoints[i].first.x, gp.gaussPoints[i].first.y, gp.gaussPoints[i].first.z, gp.gaussPoints[i].first.t, f.first.transpose) ;
		ret_ = B*f.second*B_ ;
		ret[i] = ret_*gp.gaussPoints[i].second ;
		if(std::abs(ret[i][0][0]) > 1)
	{
		B.print() ;
		f.second.print() ;
		B_.print() ;
		if(f.third.transpose)
			std::cout << "transpose-third-" ;
		std::cout << "miaou " << i << std::endl ;
	}
//		matrix_matrix_matrix_multiply_and_add(B, f.second, B_, gp.gaussPoints[i].second, ret) ;
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

Matrix VirtualMachine::ieval(const GtMLtG &f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & var)
{

	geval(f.third.f, Jinv[0], var, gp.gaussPoints[0].first, f.third.transpose,B_) ;
	geval(f.first.f, Jinv[0], var, gp.gaussPoints[0].first, f.first.transpose,B) ;
	
	Matrix ret(B*f.second[0]*B_) ;
	ret *= gp.gaussPoints[0].second ;
	
	for(size_t i = 1 ; i  < gp.gaussPoints.size() ; i++)
	{
		geval(f.third.f, Jinv[i], var, gp.gaussPoints[i].first, f.third.transpose, B_) ;
		geval(f.first.f, Jinv[i], var, gp.gaussPoints[i].first, f.first.transpose, B) ;
		matrix_matrix_matrix_multiply_and_add(B, f.second[i], B_, gp.gaussPoints[i].second, ret) ;
	}

	return ret ;
}

void VirtualMachine::ieval(const GtMLtG &f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & var, Matrix & ret)
{
	geval(f.third.f, Jinv[0], var, gp.gaussPoints[0].first, f.third.transpose,B_) ;
	geval(f.first.f, Jinv[0], var, gp.gaussPoints[0].first, f.first.transpose,B) ;
	
	ret = B*f.second[0]*B_ ;
	ret *= gp.gaussPoints[0].second ;
	
	for(size_t i = 1 ; i  < gp.gaussPoints.size() ; i++)
	{
		geval(f.third.f, Jinv[i], var, gp.gaussPoints[i].first, f.third.transpose, B_) ;
		geval(f.first.f, Jinv[i], var, gp.gaussPoints[i].first, f.first.transpose, B) ;
		matrix_matrix_matrix_multiply_and_add(B, f.second[i], B_, gp.gaussPoints[i].second, ret) ;
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
//			x_a.print() ;
			
			for(size_t i = 0 ; i < gp_b.gaussPoints.size() ; i++)
			{
				gp_b.gaussPoints[i].first.t -= default_derivation_delta ;
			}
			
			Matrix x_b (ieval(d.second, gp_b, Jinv, vars)) ;
//			x_b.print() ;

			return (x_a-x_b)/(2.*default_derivation_delta) ;
		}
	default:
		{
			std::cerr << "operator not implemented" << std::endl ;
			return Matrix() ;
		}
	}
}

void VirtualMachine::ieval(const DdGtMtG & d, const GaussPointArray &gp_, const std::valarray<Matrix> &Jinv, const IntegrableEntity * e, const std::vector<Variable> & vars, Matrix & ret)
{
	GaussPointArray gp_a(gp_);
	gp_a.id = -1 ;
	GaussPointArray gp_b(gp_);
	gp_b.id = -1 ;
	switch(d.first)
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
			
			ret = (x_a-x_b)/(2.*default_derivation_delta) ;
			break ;
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
			
			ret = (x_a-x_b)/(2.*default_derivation_delta) ;
			break ;
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
			
			ret = (x_a-x_b)/(2.*default_derivation_delta) ;
			break ;
		}
	case TIME_VARIABLE :
		{
			for(size_t i = 0 ; i < gp_a.gaussPoints.size() ; i++)
				gp_a.gaussPoints[i].first.t += default_derivation_delta*10 ;			
			std::valarray<Matrix> t_a(ievalDecomposed(d.second, gp_a, Jinv, vars)) ;
			for(size_t i = 0 ; i < gp_a.gaussPoints.size() ; i++)
				gp_a.gaussPoints[i].first.t -= default_derivation_delta*10 ;

			
			for(size_t i = 0 ; i < gp_b.gaussPoints.size() ; i++)
				gp_b.gaussPoints[i].first.t -= default_derivation_delta*10 ;			
			std::valarray<Matrix> t_b (ievalDecomposed(d.second, gp_b, Jinv, vars)) ;
			for(size_t i = 0 ; i < gp_b.gaussPoints.size() ; i++)
				gp_b.gaussPoints[i].first.t += default_derivation_delta*10 ;			
			
			size_t dim = vars.size()-1 ;
			
			
			ret = ((Matrix) ((Matrix) (t_a[0]-t_b[0])*Jinv[0][dim][dim]))/(20.*default_derivation_delta) ;
			for(size_t i = 1 ; i < gp_a.gaussPoints.size() ; i++)
			{
				ret += ((Matrix) ((Matrix) (t_a[i]-t_b[i])*Jinv[i][dim][dim]))/(20.*default_derivation_delta) ;
			}

			

			
			
			for(size_t i = 0 ; i < gp_a.gaussPoints.size() ; i++)
				gp_a.gaussPoints[i].first.x += default_derivation_delta/10 ;			
			std::valarray<Matrix> x_a(ievalDecomposed(d.second, gp_a, Jinv, vars)) ;
			for(size_t i = 0 ; i < gp_a.gaussPoints.size() ; i++)
				gp_a.gaussPoints[i].first.x -= default_derivation_delta/10 ;

			
			for(size_t i = 0 ; i < gp_b.gaussPoints.size() ; i++)
				gp_b.gaussPoints[i].first.x -= default_derivation_delta/10 ;			
			std::valarray<Matrix> x_b (ievalDecomposed(d.second, gp_b, Jinv, vars)) ;
			for(size_t i = 0 ; i < gp_b.gaussPoints.size() ; i++)
				gp_b.gaussPoints[i].first.x += default_derivation_delta/10 ;			
			
			for(size_t i = 0 ; i < gp_a.gaussPoints.size() ; i++)
			{
				ret += ((Matrix) ((Matrix) (x_a[i]-x_b[i])*Jinv[i][dim][0]))/(20.*default_derivation_delta) ;
			}
			
			
			

			
			for(size_t i = 0 ; i < gp_a.gaussPoints.size() ; i++)
				gp_a.gaussPoints[i].first.y += default_derivation_delta*10 ;			
			std::valarray<Matrix> y_a(ievalDecomposed(d.second, gp_a, Jinv, vars)) ;
			for(size_t i = 0 ; i < gp_a.gaussPoints.size() ; i++)
				gp_a.gaussPoints[i].first.y -= default_derivation_delta*10 ;

			
			for(size_t i = 0 ; i < gp_b.gaussPoints.size() ; i++)
				gp_b.gaussPoints[i].first.y -= default_derivation_delta*10 ;			
			std::valarray<Matrix> y_b (ievalDecomposed(d.second, gp_b, Jinv, vars)) ;
			for(size_t i = 0 ; i < gp_b.gaussPoints.size() ; i++)
				gp_b.gaussPoints[i].first.y += default_derivation_delta*10 ;			
			
			for(size_t i = 0 ; i < gp_a.gaussPoints.size() ; i++)
			{
				ret += ((Matrix) ((Matrix) (y_a[i]-y_b[i])*Jinv[i][dim][1]))/(20.*default_derivation_delta) ;
			}
// 			ret.print() ;
			
			
			if(std::abs(Jinv[0][0][dim]) > POINT_TOLERANCE_2D || 
			  std::abs(Jinv[0][1][dim]) > POINT_TOLERANCE_2D || 
			  std::abs(Jinv[0][dim][0]) > POINT_TOLERANCE_2D || 
			  std::abs(Jinv[0][dim][1]) > POINT_TOLERANCE_2D)
			{
			
				std::valarray<Matrix> domega (ievalDecomposed(d.second, gp_, Jinv, vars)) ;
				for(size_t i = 0 ; i < gp_a.gaussPoints.size() ; i++)
					gp_a.gaussPoints[i].first.t += default_derivation_delta*10 ;			
				for(size_t i = 0 ; i < gp_b.gaussPoints.size() ; i++)
					gp_b.gaussPoints[i].first.t -= default_derivation_delta*10 ;			
				for(size_t i = 0 ; i < gp_.gaussPoints.size() ; i++)
					ret += (domega[i]*Jinv[i][dim][dim] /(20.*default_derivation_delta)) ;
			}
			break ;
		}
	default:
		{
			std::cerr << "operator not implemented" << std::endl ;
//			ret = Matrix() ;
		}
	}
}

void VirtualMachine::ieval(const DdGtMtGD & d, const GaussPointArray &gp_, const std::valarray<Matrix> &Jinv, const IntegrableEntity * e, const std::vector<Variable> & vars, Matrix & ret)
{
	GaussPointArray gp_a(gp_);
	gp_a.id = -1 ;
	GaussPointArray gp_b(gp_);
	gp_b.id = -1 ;
	switch(d.first)
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
			
			ret = (x_a-x_b)/(2.*default_derivation_delta) ;
			break ;
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
			
			ret = (x_a-x_b)/(2.*default_derivation_delta) ;
			break ;
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
			
			ret = (x_a-x_b)/(2.*default_derivation_delta) ;
			break ;
		}
	case TIME_VARIABLE :
		{
			for(size_t i = 0 ; i < gp_a.gaussPoints.size() ; i++)
				gp_a.gaussPoints[i].first.t += default_derivation_delta*100 ;			
			std::valarray<Matrix> t_a(ievalDecomposed(d.second, gp_a, Jinv, vars)) ;
			for(size_t i = 0 ; i < gp_a.gaussPoints.size() ; i++)
				gp_a.gaussPoints[i].first.t -= default_derivation_delta*100 ;

			
			for(size_t i = 0 ; i < gp_b.gaussPoints.size() ; i++)
				gp_b.gaussPoints[i].first.t -= default_derivation_delta*100 ;			
			std::valarray<Matrix> t_b (ievalDecomposed(d.second, gp_b, Jinv, vars)) ;
			for(size_t i = 0 ; i < gp_b.gaussPoints.size() ; i++)
				gp_b.gaussPoints[i].first.t += default_derivation_delta*100 ;			
			
			size_t dim = vars.size()-1 ;
			
			
			ret = ((Matrix) ((Matrix) (t_a[0]-t_b[0])*Jinv[0][dim][dim]))/(200.*default_derivation_delta) ;
			for(size_t i = 1 ; i < gp_a.gaussPoints.size() ; i++)
			{
				ret += ((Matrix) ((Matrix) (t_a[i]-t_b[i])*Jinv[i][dim][dim]))/(200.*default_derivation_delta) ;
			}

			

			
			
			for(size_t i = 0 ; i < gp_a.gaussPoints.size() ; i++)
				gp_a.gaussPoints[i].first.x += default_derivation_delta*100 ;			
			std::valarray<Matrix> x_a(ievalDecomposed(d.second, gp_a, Jinv, vars)) ;
			for(size_t i = 0 ; i < gp_a.gaussPoints.size() ; i++)
				gp_a.gaussPoints[i].first.x -= default_derivation_delta*100 ;

			
			for(size_t i = 0 ; i < gp_b.gaussPoints.size() ; i++)
				gp_b.gaussPoints[i].first.x -= default_derivation_delta*100 ;			
			std::valarray<Matrix> x_b (ievalDecomposed(d.second, gp_b, Jinv, vars)) ;
			for(size_t i = 0 ; i < gp_b.gaussPoints.size() ; i++)
				gp_b.gaussPoints[i].first.x += default_derivation_delta*100 ;			
			
			for(size_t i = 0 ; i < gp_a.gaussPoints.size() ; i++)
			{
 				ret += ((Matrix) ((Matrix) (x_a[i]-x_b[i])*Jinv[i][dim][0]))/(200.*default_derivation_delta) ;
			}
			
			
			

			
			for(size_t i = 0 ; i < gp_a.gaussPoints.size() ; i++)
				gp_a.gaussPoints[i].first.y += default_derivation_delta*100 ;			
			std::valarray<Matrix> y_a(ievalDecomposed(d.second, gp_a, Jinv, vars)) ;
			for(size_t i = 0 ; i < gp_a.gaussPoints.size() ; i++)
				gp_a.gaussPoints[i].first.y -= default_derivation_delta*100 ;

			
			for(size_t i = 0 ; i < gp_b.gaussPoints.size() ; i++)
				gp_b.gaussPoints[i].first.y -= default_derivation_delta*100 ;			
			std::valarray<Matrix> y_b (ievalDecomposed(d.second, gp_b, Jinv, vars)) ;
			for(size_t i = 0 ; i < gp_b.gaussPoints.size() ; i++)
				gp_b.gaussPoints[i].first.y += default_derivation_delta*100 ;			
			
			for(size_t i = 0 ; i < gp_a.gaussPoints.size() ; i++)
			{
 				ret += ((Matrix) ((Matrix) (y_a[i]-y_b[i])*Jinv[i][dim][1]))/(200.*default_derivation_delta) ;
			}
// 			std::cout << "+++++" << std::endl  ;
// 			ret.print() ;
// 			std::cout << "+++++" << std::endl  ;
			
			
			std::valarray<Matrix> domega (ievalDecomposed(d.second, gp_, Jinv, vars)) ;
			for(size_t i = 0 ; i < gp_a.gaussPoints.size() ; i++)
				gp_a.gaussPoints[i].first.t += default_derivation_delta*100 ;			
			for(size_t i = 0 ; i < gp_b.gaussPoints.size() ; i++)
				gp_b.gaussPoints[i].first.t -= default_derivation_delta*100 ;			
			
			
			if(std::abs(Jinv[0][0][dim]) > POINT_TOLERANCE_2D || 
			  std::abs(Jinv[0][1][dim]) > POINT_TOLERANCE_2D || 
			  std::abs(Jinv[0][dim][0]) > POINT_TOLERANCE_2D || 
			  std::abs(Jinv[0][dim][1]) > POINT_TOLERANCE_2D)
			{
			
				std::valarray<Matrix> domega (ievalDecomposed(d.second, gp_, Jinv, vars)) ;
				for(size_t i = 0 ; i < gp_a.gaussPoints.size() ; i++)
					gp_a.gaussPoints[i].first.t += default_derivation_delta*10 ;			
				for(size_t i = 0 ; i < gp_b.gaussPoints.size() ; i++)
					gp_b.gaussPoints[i].first.t -= default_derivation_delta*10 ;			
				for(size_t i = 0 ; i < gp_.gaussPoints.size() ; i++)
					ret += (domega[i]*Jinv[i][dim][dim] /(20.*default_derivation_delta)) ;
			}

			break ;
		}
	default:
		{
			std::cerr << "operator not implemented" << std::endl ;
//			ret = Matrix() ;
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

// 	B[0].print() ;
// 	B_[0].print() ;
		
	ret = (B[0]*f.second*B_[0]);
	ret *= gp.gaussPoints[0].second ;
	
	for(size_t i = 1 ; i  < gp.gaussPoints.size() ; i++)
	{
		ret += (Matrix)(B[i]*f.second*B_[i])*gp.gaussPoints[i].second ;
	}
	
}

void VirtualMachine::ieval(const GDtMLtGD & f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars, Matrix & ret)
{
	std::valarray<Matrix> B = gdeval(f.first.f, Jinv, vars, gp, f.first.transpose) ;
	std::valarray<Matrix> B_ = gdeval(f.third.f, Jinv, vars, gp, f.third.transpose) ;

// 	B[0].print() ;
// 	B_[0].print() ;
		
	ret = (B[0]*f.second[0]*B_[0]);
	ret *= gp.gaussPoints[0].second ;
	
	for(size_t i = 1 ; i  < gp.gaussPoints.size() ; i++)
	{
		ret += (Matrix)(B[i]*f.second[i]*B_[i])*gp.gaussPoints[i].second ;
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
		ret += (Matrix)(B[i]*f.second*B_[i])*gp.gaussPoints[i].second ;
	}
}

void VirtualMachine::ieval(const GDtMLtG & f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars, Matrix & ret)
{
	std::valarray<Matrix> B = gdeval(f.first.f, Jinv, vars, gp, f.first.transpose) ;
	std::valarray<Matrix> B_ = geval(f.third.f, Jinv, vars, gp, f.third.transpose) ;

	ret = (B[0]*f.second[0]*B_[0]);
	ret *= gp.gaussPoints[0].second ;
	
	for(size_t i = 1 ; i  < gp.gaussPoints.size() ; i++)
	{
		ret += (Matrix)(B[i]*f.second[i]*B_[i])*gp.gaussPoints[i].second ;
	}
}

void VirtualMachine::ieval(const GDDtMtG & f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars, Matrix & ret)
{
	std::valarray<Matrix> B = gddeval(f.first.f, Jinv, vars, gp, f.first.transpose) ;
	std::valarray<Matrix> B_ = geval(f.third.f, Jinv, vars, gp, f.third.transpose) ;
	
	ret = (B[0]*f.second*B_[0]);
	ret *= gp.gaussPoints[0].second ;
	
	for(size_t i = 1 ; i  < gp.gaussPoints.size() ; i++)
	{
		ret += (Matrix)(B[i]*f.second*B_[i])*gp.gaussPoints[i].second ;
	}
}

void VirtualMachine::ieval(const GDDtMLtG & f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars, Matrix & ret)
{
	std::valarray<Matrix> B = gddeval(f.first.f, Jinv, vars, gp, f.first.transpose) ;
// 	B[0].print() ;
	std::valarray<Matrix> B_ = geval(f.third.f, Jinv, vars, gp, f.third.transpose) ;

	ret = (B[0]*f.second[0]*B_[0]);
	ret *= gp.gaussPoints[0].second ;
	
	for(size_t i = 1 ; i  < gp.gaussPoints.size() ; i++)
	{
//		B[i].print() ;
		ret += (Matrix)(B[i]*f.second[i]*B_[i])*gp.gaussPoints[i].second ;
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

// 	f.second.print() ;
// 	std::cout << B.size() << "\t" << B_.size() << std::endl ;
// 	B[0].print() ;
// 	B_[0].print() ;
// 
	ret = (B[0]*f.second*B_[0]);
	ret *= gp.gaussPoints[0].second ;


	for(size_t i = 1 ; i  < gp.gaussPoints.size() ; i++)
	{
		ret += (Matrix)(B[i]*f.second*B_[i])*gp.gaussPoints[i].second ;
	}
}

void VirtualMachine::ieval(const GtMLtGD & f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars, Matrix & ret)
{
	std::valarray<Matrix> B = geval(f.first.f, Jinv, vars, gp, f.first.transpose) ;
	std::valarray<Matrix> B_ = gdeval(f.third.f, Jinv, vars, gp, f.third.transpose) ;

	ret = (B[0]*f.second[0]*B_[0]);
	ret *= gp.gaussPoints[0].second ;

	for(size_t i = 1 ; i  < gp.gaussPoints.size() ; i++)
	{
		ret += (Matrix)(B[i]*f.second[i]*B_[i])*gp.gaussPoints[i].second ;
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
		r_ *= gp.gaussPoints[i].second ;
		ret += r_ ;
	}
	
	return ret ;
}

Vector VirtualMachine::ieval(const GtVL &f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & var)
{
	if(!Jinv.size())
	{
		int size = var.size() - (var[var.size()-1] == TIME_VARIABLE) ;
		return Vector(double(0), size) ;
	}
	Matrix M ;
	geval(f.first.f, Jinv[0],var, gp.gaussPoints[0].first.x, gp.gaussPoints[0].first.y, gp.gaussPoints[0].first.z, gp.gaussPoints[0].first.t, f.first.transpose, M);
	Vector temp(double(0),f.second[0].size()) ;
	for(size_t i = 0 ; i < temp.size() ; i++)
		temp[i] = f.second[0][i] ;
	
	Vector r_  = M*temp ;
	Vector ret = r_ * gp.gaussPoints[0].second;
//	std::cout << gp.gaussPoints.size() << "\t" << f.second.size() << std::endl ;

	for(size_t i = 1 ; i < gp.gaussPoints.size() ; i++)
	{
		 M = geval(f.first.f, Jinv[i],var, gp.gaussPoints[i].first.x, gp.gaussPoints[i].first.y, gp.gaussPoints[i].first.z, gp.gaussPoints[i].first.t, f.first.transpose) ;
		
		for(size_t j = 0 ; j < temp.size() ; j++)
		{
			temp[j] = f.second[i][j] ;
		}
		r_ = 0. ;
		r_  = M*temp ;
		r_ *= gp.gaussPoints[i].second ;
		ret += r_ ;
	}
	
	return ret ;
}

Vector VirtualMachine::ieval(const GDtV &f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & var)
{
	if(!Jinv.size())
	{
		int size = var.size() - (var[var.size()-1] == TIME_VARIABLE) ;
		return Vector(double(0), size) ;
	}
	Matrix M ;
	gdeval(f.first.f, Jinv[0],var, gp.gaussPoints[0].first.x, gp.gaussPoints[0].first.y, gp.gaussPoints[0].first.z, gp.gaussPoints[0].first.t, f.first.transpose, M);
	Vector temp(double(0),f.second.size()) ;
	for(size_t i = 0 ; i < temp.size() ; i++)
		temp[i] = f.second[i] ;
	
	Vector r_  = M*temp ;
	Vector ret = r_ * gp.gaussPoints[0].second;
	for(size_t i = 1 ; i < gp.gaussPoints.size() ; i++)
	{
		  gdeval(f.first.f, Jinv[i],var, gp.gaussPoints[i].first.x, gp.gaussPoints[i].first.y, gp.gaussPoints[i].first.z, gp.gaussPoints[i].first.t, f.first.transpose, M) ;
		
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
/*		for(size_t j = 0 ; j < temp.size() ; j++)
			temp[j] = f.second[i*temp.size()+j] ;*/
		Vector r_  = B*temp ;
		ret += r_[0] * gp.gaussPoints[i].second ;
	}
	
	return ret ;
}

std::vector<Point> VirtualMachine::allHints(const Function &f0, const Function &f1,IntegrableEntity *e )
{
	std::vector<Point> hints(0) ;
	
	if ((e!= nullptr && f0.hasIntegrationHint()) || f1.hasIntegrationHint())
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
	

	if (e!= nullptr && f0.hasIntegrationHint() )
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
