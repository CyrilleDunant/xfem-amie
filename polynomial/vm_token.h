//
// C++ Interface: vm_token
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef VM_TOKEN_H
#define VM_TOKEN_H

#include <iostream>
#include <sstream>
#include <vector>

#include "variable.h"
#include "../geometry/geometry_base.h"

double sign(const double t) ;
double positivity(const double t) ;
double negativity(const double t) ;

namespace Mu
{

typedef double (*unaryFunctionPointer)(const double) ;
typedef double (*binaryFunctionPointer)(const double, double) ;
typedef double (*trinaryFunctionPointer)(const double, double, double) ;

// typedef std::vector<double> Memory ;

typedef enum 
{
	TOKEN,
	TOKEN_CONSTANT,
	TOKEN_NULL,
	TOKEN_POSITION,
	TOKEN_POSITION_OPERATOR,
	TOKEN_DOMAIN,
	TOKEN_WRITE_VARIABLE,
	TOKEN_READ_VARIABLE,
	TOKEN_PROJECTION,
	TOKEN_SIGNED_DISTANCE,
	TOKEN_POINT_DISTANCE_OPERATOR,
	TOKEN_POINT_SQUARE_DISTANCE_OPERATOR,
	TOKEN_X,
	TOKEN_Y,
	TOKEN_Z,
	TOKEN_T,
	TOKEN_U,
	TOKEN_V,
	TOKEN_W,
	TOKEN_VARIABLE,
	TOKEN_UNARY_FUNCTION,
	TOKEN_BINARY_FUNCTION,
	TOKEN_PLUS,
	TOKEN_MINUS,
	TOKEN_TIMES,
	TOKEN_DIVIDES,
	TOKEN_POWER,
	TOKEN_ABS,
	TOKEN_COS,
	TOKEN_SIN,
	TOKEN_TAN,
	TOKEN_COSH,
	TOKEN_SINH,
	TOKEN_TANH,
	TOKEN_EXP,
	TOKEN_SIGN,
	TOKEN_POSITIVITY,
	TOKEN_NEGATIVITY,
	TOKEN_LOG,
	TOKEN_SQRT,
	TOKEN_ATAN2,
	TOKEN_X_POWER_AND_MULTIPLY,
	TOKEN_Y_POWER_AND_MULTIPLY,
	TOKEN_Z_POWER_AND_MULTIPLY,
	TOKEN_T_POWER_AND_MULTIPLY,
	TOKEN_X_POWER_CONST,
	TOKEN_Y_POWER_CONST,
	TOKEN_Z_POWER_CONST,
	TOKEN_T_POWER_CONST,
	TOKEN_ADD_X_AND_CONST,
	TOKEN_ADD_Y_AND_CONST,
	TOKEN_ADD_Z_AND_CONST,
	TOKEN_ADD_T_AND_CONST,
	TOKEN_MULTIPLY_X_AND_CONST,
	TOKEN_MULTIPLY_Y_AND_CONST,
	TOKEN_MULTIPLY_Z_AND_CONST,
	TOKEN_MULTIPLY_T_AND_CONST,
	TOKEN_MULTIPLY_X_Y_AND_CONST,
	TOKEN_MULTIPLY_X_Z_AND_CONST,
	TOKEN_MULTIPLY_Y_Z_AND_CONST,
	TOKEN_SUBSTRACT_X_WITH_CONST,
	TOKEN_SUBSTRACT_Y_WITH_CONST,
	TOKEN_SUBSTRACT_Z_WITH_CONST,
	TOKEN_SUBSTRACT_T_WITH_CONST,
	TOKEN_SUBSTRACT_CONST_WITH_X,
	TOKEN_SUBSTRACT_CONST_WITH_Y,
	TOKEN_SUBSTRACT_CONST_WITH_Z,
	TOKEN_SUBSTRACT_CONST_WITH_T,
	TOKEN_DIVIDE_X_WITH_CONST,
	TOKEN_DIVIDE_Y_WITH_CONST,
	TOKEN_DIVIDE_Z_WITH_CONST,
	TOKEN_DIVIDE_T_WITH_CONST,
	TOKEN_DIVIDE_CONST_WITH_X,
	TOKEN_DIVIDE_CONST_WITH_Y,
	TOKEN_DIVIDE_CONST_WITH_Z,
	TOKEN_DIVIDE_CONST_WITH_T,
	TOKEN_READ_ADD_CONST,
	TOKEN_READ_ADD_READ,
	TOKEN_READ_SUBSTRACT_CONST,
	TOKEN_CONST_SUBSTRACT_READ,
	TOKEN_READ_MULTIPLY_CONST,
	TOKEN_READ_DIVIDE_CONST,
	TOKEN_CONST_DIVIDE_READ,
	TOKEN_READ_POWER_CONST,
	
} TokenTypeId ;

typedef std::pair<std::pair<TokenTypeId, short unsigned int>, double> TokenType ;

struct Memory
{
	Vector stack;
	Vector heap ;
	int top_pos ;
	Memory() : stack(0., 64), heap(0., 256), top_pos(-1)
	{
	} ;
	
	size_t size() const
	{
		return top_pos+1 ;
	}
	
	double & operator[](const size_t i)
	{
// 		if(i > 31)
// 			std::cout << "stack overflow !" << std::endl ;
		return stack[i] ;
	}
	
	const double & operator[](const size_t i) const
	{
// 		if(i > 31)
// 			std::cout << "stack overflow !" << std::endl ;
		return stack[i] ;
	}
	
	void push_back(const double & a)
	{
// 		if(top_pos > 30)
// 			std::cout << "stack overflow !" << std::endl ;
		stack[++top_pos] = a ;
	}
	void pop_back()
	{
		top_pos-- ;
	}
	
	void reset()
	{
		top_pos = -1 ;
	}
} ;

struct Context
{
	Memory &memory; 
	const double & x ; 
	const double & y ; 
	const double & z ;
	const double & t ; 
	const double & u ; 
	const double & v ; 
	const double & w ;
	Context(Memory &stack_, const double & x_, const double & y_ = 0, const double & z_ = 0, const double & t_ = 0, const double & u_ = 0, const double & v_ = 0, const double & w_ = 0) : memory(stack_), x(x_), y(y_), z(z_), t(t_), u(u_), v(v_), w(w_) { } ;
} ;

class Token
{
public:
	
	
	Token(bool null = false, TokenType t =  std::make_pair(std::make_pair(TOKEN,0),(double)(0))) : isNull(null), type(t) { ; };
	
	virtual ~Token() { } ;
	
	virtual void eval(Context & context) const = 0;
	
	virtual std::string print() const  = 0 ;
	
	const bool isNull ;
	const TokenType type ;
} ;

class ConstantToken : public Token
{
	

public:
	ConstantToken(double v, bool nul = false );
	
	virtual void eval(Context & context)  const 
	{
		context.memory.push_back(type.second) ;
	}

	virtual ~ConstantToken();
	
	virtual std::string print() const
	{
		std::ostringstream s ;
		s << type.second ;
		return s.str() ;
	}
	
} ;

class NullToken : public Token
{
public:
	NullToken() : Token(true, std::make_pair(std::make_pair(TOKEN_NULL,0),(double)(0))) { };
	virtual void eval(Context & context)  const 
	{
		context.memory.push_back(0) ;
	}
	virtual ~NullToken();
	virtual std::string print() const
	{
		std::ostringstream s ;
		s << "0" ;
		return s.str() ;
	}
} ;

class PositionToken : public Token
{
	std::vector<Segment> s ;
	Point w ;
public:
	PositionToken(Segment s_ ) : Token(false, std::make_pair(std::make_pair(TOKEN_POSITION, 0), 0))
	{
		Point vector(-s_.vector().y, s_.vector().x) ;
		w=( s_.midPoint()+vector*10.) ;
		s.push_back(s_)  ;
	}
	
	PositionToken(std::vector<Segment> s_ ) : Token(false, std::make_pair(std::make_pair(TOKEN_POSITION, 0), 0))
	{
		for(size_t i = 0 ; i < s_.size() ; i++)
		{

			s.push_back(s_[i])  ;
		}			
		Point vector(-s_[0].vector().y, s_[0].vector().x) ;
		w=( s_[0].midPoint()+vector*10.) ;
	}
	virtual void eval(Context & context) const
	{
		if(s.empty())
		{
			context.memory.push_back(0) ;
			return ;
		}

		Point test(context.x,context.y) ;
		size_t intersections = 0 ;
		for(size_t i = 0 ; i < s.size() ; i++)
		{
			if(s[i].intersects(test, w))
			{
				intersections++ ;
			}
		}
		
		context.memory.stack[context.memory.top_pos] =  (intersections & 1) * 2 - 1;
// 		if(intersections%2 == 1)
// 		    context.memory.stack[context.memory.top_pos] = -1 ;
// 		else
// 		    context.memory.stack[context.memory.top_pos] = 1 ;
	}
	
	virtual ~PositionToken();
	virtual std::string print() const
	{
		return std::string("position") ;
	}
} ;

class LineDistanceOperatorToken : public Token
{
	Line l ;
public:
		LineDistanceOperatorToken(const Line & l_ ) : Token(false, std::make_pair(std::make_pair(TOKEN_POSITION_OPERATOR, 0), (double)(0))), l(l_)
	{
	}
	
	virtual void eval(Context & context) const
	{
		
		Point test(context.memory.stack[context.memory.top_pos-1], context.memory.stack[context.memory.top_pos]) ;
	    context.memory.pop_back() ;
	
		
		context.memory.stack[context.memory.top_pos] =  dist(test,l.projection(test));
	}
	
	virtual ~LineDistanceOperatorToken() { };
	virtual std::string print() const
	{
		return std::string("distToLine") ;
	}
} ;

class PositionOperatorToken : public Token
{
	std::vector<Segment> s ;
	Point w ;
public:
	PositionOperatorToken(Segment s_ ) : Token(false, std::make_pair(std::make_pair(TOKEN_POSITION_OPERATOR, 0), (double)(0)))
	{
		Point vector(-s_.vector().y, s_.vector().x) ;
		w= s_.midPoint()+vector*1000. ;
		s.push_back(s_)  ;
	}
	
	PositionOperatorToken(std::vector<Segment> s_ ) : Token(false)
	{
		
		for(size_t i = 0 ; i < s_.size() ; i++)
		{

			s.push_back(s_[i])  ;
		}
		
		Point vector(-s[0].vector().y, s[0].vector().x) ;
		w= s[0].midPoint()+vector*100. ;
	}
	
	virtual void eval(Context & context) const
	{
		if(s.empty())
		{
			context.memory.pop_back() ;
			context.memory.stack[context.memory.top_pos] = 0 ;
			return ;
		}
		

		Point test(context.memory.stack[context.memory.top_pos-1], context.memory.stack[context.memory.top_pos]) ;
	    context.memory.pop_back() ;
		
		int intersections = 0 ;
		for(size_t i = 0 ; i < s.size() ; i++)
		{
			if(s[i].on(test))
			{
				context.memory.stack[context.memory.top_pos] = 0 ;
				return ;
			}
			if(s[i].intersects(test, w))
			{
				intersections++ ;
			}
			
		}
		
		if((intersections % 2)  != 0)
			context.memory.stack[context.memory.top_pos] = 1;
		else
			context.memory.stack[context.memory.top_pos] =  -1;

	}
	
	virtual ~PositionOperatorToken() { };
	virtual std::string print() const
	{
		return std::string("position") ;
	}
} ;

class DomainToken : public Token
{
	const Geometry* geo ;
public:
	DomainToken(const Geometry * g ) : Token(false, std::make_pair(std::make_pair(TOKEN_DOMAIN, 0), (double)(0)))
	{
		geo = g ;
	}
	
	virtual void eval(Context & context) const
	{
		Point p(context.x,context.y,context.z) ;
		if(geo->in(p))
		    context.memory.push_back(1) ;
		else
		    context.memory.push_back(-1) ;
	}
	virtual ~DomainToken();
	virtual std::string print() const
	{
		return std::string("domain") ;
	}
} ;

class DomainBinaryOperatorToken : public Token
{
	const Geometry* geo ;
public:
	DomainBinaryOperatorToken(const Geometry * g ) : Token(false, std::make_pair(std::make_pair(TOKEN_DOMAIN, 0), (double)(0)))
	{
		geo = g ;
	}
	
	virtual void eval(Context & context) const
	{
		
		Point p(context.memory.stack[context.memory.top_pos], context.memory.stack[context.memory.top_pos-1]) ;
		context.memory.pop_back() ;
		if(geo->in(p))
			context.memory.stack[context.memory.top_pos] = 1 ;
		else
			context.memory.stack[context.memory.top_pos] = -1 ;
	}
	virtual ~DomainBinaryOperatorToken() { };
	virtual std::string print() const
	{
		return std::string("domainBinOp") ;
	}
} ;

class PointDistanceBinaryOperatorToken : public Token
{
	Point base ;
public:
	PointDistanceBinaryOperatorToken(const Point & p ) : Token(false, std::make_pair(std::make_pair(TOKEN_POINT_DISTANCE_OPERATOR, 0), (double)(0)))
	{
		base = p ;
	}
	
	virtual void eval(Context & context) const
	{
		
		Point p(context.memory.stack[context.memory.top_pos], context.memory.stack[context.memory.top_pos-1]) ;
		context.memory.pop_back() ;

		context.memory.stack[context.memory.top_pos] = dist(p, base) ;

	}
	virtual ~PointDistanceBinaryOperatorToken() { };
	virtual std::string print() const
	{
		return std::string("pointDistBinOp") ;
	}
} ;

class PointSquareDistanceBinaryOperatorToken : public Token
{
	Point base ;
public:
	PointSquareDistanceBinaryOperatorToken(const Point & p ) : Token(false, std::make_pair(std::make_pair(TOKEN_POINT_SQUARE_DISTANCE_OPERATOR, 0), (double)(0)))
	{
		base = p ;
	}
	
	virtual void eval(Context & context) const
	{
		
		Point p(context.memory.stack[context.memory.top_pos], context.memory.stack[context.memory.top_pos-1]) ;
		context.memory.pop_back() ;

		context.memory.stack[context.memory.top_pos] = squareDist2D(p, base) ;

	}
	virtual ~PointSquareDistanceBinaryOperatorToken() { };
	virtual std::string print() const
	{
		return std::string("pointSqDistBinOp") ;
	}
} ;

class ReadHeapVariableToken : public Token
{

public:
	ReadHeapVariableToken(const size_t a) : Token(false, std::make_pair(std::make_pair(TOKEN_READ_VARIABLE, a), (double)(0))) { } ;
	
	virtual void eval(Context & context) const
	{
	    context.memory.push_back(context.memory.heap[type.first.second]) ;
	}
	
	virtual std::string print() const
	{
		std::ostringstream s ;
		s << "read@" << type.first.second;
		return s.str() ;
	}
} ;

class SetHeapVariableToken : public Token
{

public:
	SetHeapVariableToken(const size_t a) : Token(false, std::make_pair(std::make_pair(TOKEN_WRITE_VARIABLE, a),(double)(0))) { } ;

	
	virtual void eval(Context & context) const
	{
	    context.memory.heap[type.first.second] = context.memory.stack[context.memory.top_pos] ;
	    context.memory.pop_back() ;
	}
	
	virtual std::string print() const
	{
		std::ostringstream s ;
		s << "write@" << type.first.second;
		return s.str() ;
	}
} ;

class ProjectionToken : public Token
{
	Point normal ;
public:
	ProjectionToken(Segment s_ ) : Token(false, std::make_pair(std::make_pair(TOKEN_PROJECTION, 0),(double)(0)))
	
	{
		normal = Point(-s_.vector().y, s_.vector().x) ;
	}
	
	virtual void eval(Context & context) const
	{
	    context.memory.push_back(Point(context.x, context.y, context.z) * normal) ;
	}
	virtual ~ProjectionToken();
	virtual std::string print() const
	{
		return std::string("check position") ;
	}
} ;

class XToken : public Token
{
public:
	XToken() : Token(false, std::make_pair(std::make_pair(TOKEN_X, 0),(double)(0))) { };
	virtual void eval(Context & context) const {context.memory.push_back(context.x) ;}
	virtual ~XToken(){ };
	virtual std::string print() const { return std::string("x") ;}
} ;

class YToken : public Token
{
public:
	YToken() : Token(false, std::make_pair(std::make_pair(TOKEN_Y, 0),(double)(0))) { };
	virtual void eval(Context & context) const {context.memory.push_back(context.y) ;}
	virtual ~YToken(){ };
	virtual std::string print() const { return std::string("y") ;}
} ;

class ZToken : public Token
{
public:
	ZToken() : Token(false, std::make_pair(std::make_pair(TOKEN_Z, 0),(double)(0))) { };
	virtual void eval(Context & context) const {context.memory.push_back(context.z) ;}
	virtual ~ZToken(){ };
	virtual std::string print() const { return std::string("z") ;}
} ;

class TToken : public Token
{
public:
	TToken() : Token(false, std::make_pair(std::make_pair(TOKEN_T, 0),(double)(0))) { };
	virtual void eval(Context & context) const {context.memory.push_back(context.t) ;}
	virtual ~TToken(){ };
	virtual std::string print() const { return std::string("t") ;}
} ;

class UToken : public Token
{
public:
	UToken() : Token(false, std::make_pair(std::make_pair(TOKEN_U, 0),(double)(0))) { };
	virtual void eval(Context & context) const {context.memory.push_back(context.u) ;}
	virtual ~UToken(){ };
	virtual std::string print() const { return std::string("u") ;}
} ;

class VToken : public Token
{
public:
	VToken() : Token(false, std::make_pair(std::make_pair(TOKEN_V, 0),(double)(0))) { };
	virtual void eval(Context & context) const {context.memory.push_back(context.v) ;}
	virtual ~VToken(){ };
	virtual std::string print() const { return std::string("v") ;}
} ;

class WToken : public Token
{
public:
	WToken() : Token(false, std::make_pair(std::make_pair(TOKEN_W, 0),(double)(0))) { };
	virtual void eval(Context & context) const {context.memory.push_back(context.w) ;}
	virtual ~WToken(){ };
	virtual std::string print() const { return std::string("w") ;}
} ;

class VariableToken : public Token
{
protected:
	const Variable var ;
public:
	VariableToken(const Variable var , bool nul = false);
	
	virtual void eval(Context & context) const {
	double v ;
	switch (var)
	{
	case ONE:
		{
			v = 1 ;
			break ;
		}
	case XI:
		{
			v = context.x ;
			break ;
		}
	case ETA:
		{
			v = context.y ;
			break ;
		}
	case ZETA:
		{
			v = context.z ;
			break ;
		}
	case TIME_VARIABLE:
		{
			v = 1 ;
			break ;
		}
	case U_VARIABLE:
		{
			v = context.u ;
			break ;
		}
	case V_VARIABLE:
		{
			v = context.v ;
			break ;
		}
	case W_VARIABLE:
		{
			v = context.w ;
			break ;
		}
	default:
		{
			v = 1 ;
			break ;
		}
	}
    context.memory.push_back(v) ;
}
	
	virtual ~VariableToken();
	
	virtual std::string print() const
	{
		switch(var)
		{
		case ONE:
			{
				return (std::string("1")) ;
			}
		case XI:
			{
				return (std::string("x")) ;
			}
		case ETA:
			{
				return (std::string("y")) ;
			}
		case ZETA:
			{
				return (std::string("z")) ;
			}
		case TIME_VARIABLE:
			{
				return (std::string("t")) ;
			}
		case U_VARIABLE:
			{
				return (std::string("u")) ;
			}
		case V_VARIABLE:
			{
				return (std::string("v")) ;
			}
		case W_VARIABLE:
			{
				return (std::string("w")) ;
			}
		}
		
		return std::string("") ;
	}
} ;

class UnaryFunctionToken : public Token
{
protected:
	const unaryFunctionPointer fctPtr ;
	
public:
	UnaryFunctionToken(unaryFunctionPointer f, bool nul = false) ;
	
	virtual void eval(Context & context) const
	{
	    context.memory.stack[context.memory.top_pos] = fctPtr(context.memory.stack[context.memory.top_pos]) ;
	}
	
	virtual ~UnaryFunctionToken() ;
	
	virtual std::string print() const
	{
		return std::string("1-aryFct")  ;
	}
} ;

class BinaryFunctionToken : public Token
{
protected:
	const binaryFunctionPointer fctPtr ;
	
public:
	BinaryFunctionToken(binaryFunctionPointer f, bool nul = false);
	
	virtual void eval(Context & context) const
	{
		double new_val_0 = context.memory.stack[context.memory.top_pos-1] ;
		double new_val_1 = context.memory.stack[context.memory.top_pos] ;
	    context.memory.pop_back() ;
	    context.memory.stack[context.memory.top_pos] = fctPtr(new_val_0, new_val_1) ;
	}
	
	virtual ~BinaryFunctionToken() ;
	
	virtual std::string print() const
	{
		return std::string("2-aryFct")  ;
	}
} ;

class CosToken : public Token
{
public:
	CosToken() : Token(false, std::make_pair(std::make_pair(TOKEN_COS, 0),(double)(0))) { };
	
	virtual void eval(Context & context) const
	{
	    context.memory.stack[context.memory.top_pos] = cos(context.memory.stack[context.memory.top_pos]) ;
	}
	
	virtual ~CosToken() {} ;
	
	virtual std::string print() const
	{
		return std::string("cos")  ;
	}
} ;

class AbsToken : public Token
{
public:
	AbsToken() : Token(false, std::make_pair(std::make_pair(TOKEN_ABS, 0),(double)(0))) { };
	
	virtual void eval(Context & context) const
	{
		context.memory.stack[context.memory.top_pos] = std::abs(context.memory.stack[context.memory.top_pos]) ;
	}
	
	virtual ~AbsToken() {} ;
	
	virtual std::string print() const
	{
		return std::string("abs")  ;
	}
} ;

class TanToken : public Token
{
public:
	TanToken() : Token(false, std::make_pair(std::make_pair(TOKEN_TAN, 0),(double)(0))) { };
	
	virtual void eval(Context & context) const
	{
	    context.memory.stack[context.memory.top_pos] = tan(context.memory.stack[context.memory.top_pos]) ;
	}
	
	virtual ~TanToken() {} ;
	
	virtual std::string print() const
	{
		return std::string("tan")  ;
	}
} ;

class SinToken : public Token
{
public:
	SinToken() : Token(false, std::make_pair(std::make_pair(TOKEN_SIN, 0),(double)(0)) ) { };
	
	virtual void eval(Context & context) const
	{
	    context.memory.stack[context.memory.top_pos] = sin(context.memory.stack[context.memory.top_pos]) ;
	}
	
	virtual ~SinToken() {} ;
	
	virtual std::string print() const
	{
		return std::string("sin")  ;
	}
} ;

class ExpToken : public Token
{
public:
	ExpToken() : Token(false, std::make_pair(std::make_pair(TOKEN_EXP, 0),(double)(0)) ) { };
	
	virtual void eval(Context & context) const
	{
	    context.memory.stack[context.memory.top_pos] = exp(context.memory.stack[context.memory.top_pos]) ;
	}
	
	virtual ~ExpToken() {} ;
	
	virtual std::string print() const
	{
		return std::string("exp")  ;
	}
} ;

class SignFunctionToken : public Token
{
public:
	SignFunctionToken() : Token(false, std::make_pair(std::make_pair(TOKEN_SIGN, 0),(double)(0)) ) { };
	
	virtual void eval(Context & context) const
	{
	    context.memory.stack[context.memory.top_pos] = sign(context.memory.stack[context.memory.top_pos]) ;
	}
	
	virtual ~SignFunctionToken() {} ;
	
	virtual std::string print() const
	{
		return std::string("sign")  ;
	}
} ;

class PositivityFunctionToken : public Token
{
public:
	PositivityFunctionToken() : Token(false, std::make_pair(std::make_pair(TOKEN_SIGN, 0),(double)(0)) ) { };
	
	virtual void eval(Context & context) const
	{
		double s = sign(context.memory.stack[context.memory.top_pos]) ;
		if(s > 0)
			context.memory.stack[context.memory.top_pos] = 1 ;
		else
			context.memory.stack[context.memory.top_pos] = 0 ;
	}
	
	virtual ~PositivityFunctionToken() {} ;
	
	virtual std::string print() const
	{
		return std::string(">0")  ;
	}
} ;

class NegativityFunctionToken : public Token
{
public:
	NegativityFunctionToken() : Token(false, std::make_pair(std::make_pair(TOKEN_SIGN, 0),(double)(0)) ) { };
	
	virtual void eval(Context & context) const
	{
		double s = sign(context.memory.stack[context.memory.top_pos]) ;
		if(s < 0)
	    	context.memory.stack[context.memory.top_pos] = 1 ;
		else
			context.memory.stack[context.memory.top_pos] = 0 ;
	}
	
	virtual ~NegativityFunctionToken() {} ;
	
	virtual std::string print() const
	{
		return std::string("<0")  ;
	}
} ;

class LogToken : public Token
{
public:
	LogToken() : Token(false, std::make_pair(std::make_pair(TOKEN_LOG, 0),(double)(0)))  { };
	
	virtual void eval(Context & context) const
	{
	    context.memory.stack[context.memory.top_pos] = log(context.memory.stack[context.memory.top_pos]) ;
	}
	
	virtual ~LogToken() {} ;
	
	virtual std::string print() const
	{
		return std::string("log")  ;
	}
} ;

class CoshToken : public Token
{
public:
	CoshToken() : Token(false, std::make_pair(std::make_pair(TOKEN_COSH, 0),(double)(0)) ) { };
	
	virtual void eval(Context & context) const
	{
	    context.memory.stack[context.memory.top_pos] = cosh(context.memory.stack[context.memory.top_pos]) ;
	}
	
	virtual ~CoshToken() {} ;
	
	virtual std::string print() const
	{
		return std::string("cosh")  ;
	}
} ;

class SinhToken : public Token
{
public:
	SinhToken() : Token(false, std::make_pair(std::make_pair(TOKEN_SINH, 0),(double)(0)) ) { };
	
	virtual void eval(Context & context) const
	{
	    context.memory.stack[context.memory.top_pos] = sinh(context.memory.stack[context.memory.top_pos]) ;
	}
	
	virtual ~SinhToken() {} ;
	
	virtual std::string print() const
	{
		return std::string("sinh")  ;
	}
} ;

class TanhToken : public Token
{
public:
	TanhToken() : Token(false, std::make_pair(std::make_pair(TOKEN_TANH, 0),(double)(0)) ) { };
	
	virtual void eval(Context & context) const
	{
	    context.memory.stack[context.memory.top_pos] = tanhl(context.memory.stack[context.memory.top_pos]) ;
	}
	
	virtual ~TanhToken() {} ;
	
	virtual std::string print() const
	{
		return std::string("tanh")  ;
	}
} ;

class SqrtToken : public Token
{
public:
	SqrtToken() : Token(false, std::make_pair(std::make_pair(TOKEN_SQRT, 0),(double)(0)))  { };
	
	virtual void eval(Context & context) const
	{
	    context.memory.stack[context.memory.top_pos] = sqrt(context.memory.stack[context.memory.top_pos]) ;
	}
	
	virtual ~SqrtToken() {} ;
	
	virtual std::string print() const
	{
		return std::string("sqrt")  ;
	}
} ;

class Atan2Token : public Token
{
public:
	Atan2Token() : Token(false, std::make_pair(std::make_pair(TOKEN_ATAN2, 0),(double)(0)))  { };
	
	virtual void eval(Context & context) const
	{
		double new_val = atan2(context.memory.stack[context.memory.top_pos-1], context.memory.stack[context.memory.top_pos]) ;
	    context.memory.pop_back() ;
	    context.memory.stack[context.memory.top_pos] = new_val ;
	}
	
	virtual ~Atan2Token() {} ;
	
	virtual std::string print() const
	{
		return std::string("atan2")  ;
	}
} ;

class PlusOperatorToken : public Token
{
public:
	PlusOperatorToken( bool nul = false) ;
	
	virtual void eval(Context & context) const
	{
		double new_val = context.memory.stack[context.memory.top_pos-1]+context.memory.stack[context.memory.top_pos] ;
	    context.memory.pop_back() ;
	    context.memory.stack[context.memory.top_pos] = new_val ;
	}
	
	virtual ~PlusOperatorToken() ;
	
	virtual std::string print() const
	{
		return std::string("+")  ;
	}
} ;

class MinusOperatorToken : public Token
{
public:
	MinusOperatorToken( bool nul = false) ;
	
	virtual void eval(Context & context) const
	{
		double new_val = context.memory.stack[context.memory.top_pos-1]-context.memory.stack[context.memory.top_pos] ;
	    context.memory.pop_back() ;
	    context.memory.stack[context.memory.top_pos] = new_val ;
	}
	
	virtual ~MinusOperatorToken() ;
	
	virtual std::string print() const
	{
		return std::string("-")  ;
	}
} ;

class TimesOperatorToken : public Token
{
public:
	TimesOperatorToken(bool nul = false);
	
	virtual void eval(Context & context) const
	{
		double new_val = context.memory.stack[context.memory.top_pos-1] * context.memory.stack[context.memory.top_pos] ;
	    context.memory.pop_back() ;
	    context.memory.stack[context.memory.top_pos] = new_val ;
	}
	
	virtual ~TimesOperatorToken() ;
	
	virtual std::string print() const
	{
		return std::string("*")  ;
	}
} ;

class DivideOperatorToken : public Token
{
public:
	DivideOperatorToken(bool nul = false) ;
	
	virtual void eval(Context & context) const
	{
		double new_val = context.memory.stack[context.memory.top_pos-1] / context.memory.stack[context.memory.top_pos];
	    context.memory.pop_back() ;
	    context.memory.stack[context.memory.top_pos] = new_val ;
	}
	
	virtual ~DivideOperatorToken() ;
	
	virtual std::string print() const
	{
		return std::string("/")  ;
	}
} ;

class PowerOperatorToken : public Token
{
public:
	PowerOperatorToken(bool nul = false) ;
	
	virtual void eval(Context & context) const
	{
		double val = context.memory.stack[context.memory.top_pos-1] ;
		size_t pow = static_cast<size_t>(context.memory.stack[context.memory.top_pos]) -1;
	    context.memory.pop_back() ;
	    context.memory.stack[context.memory.top_pos] = val ;
		for(size_t i = 0 ; i < pow ; ++i)
		{
		    context.memory.stack[context.memory.top_pos] *= val;
		}
	}
	
	virtual ~PowerOperatorToken();
	
	virtual std::string print() const
	{
		return std::string("^")  ;
	}
} ;


class XPowerConstToken: public Token
{
public:
	XPowerConstToken(const double& pow) : Token(false, std::make_pair(std::make_pair(TOKEN_X_POWER_CONST,0), pow)) { } ;
	virtual ~XPowerConstToken() { } ;
	
	virtual void eval(Context & context) const
	{
		double val = context.x ;
		int pow = (int)type.second-1 ;
		
		for(int i = 0 ; i < pow ; ++i)
		{
			val *= context.x;
		}
	    context.memory.push_back(val) ;
	}
	virtual std::string print() const
	{
		std::ostringstream s ;
		s << "< x^" << type.second << " >";
		return s.str() ;
	}
} ;

class YPowerConstToken: public Token
{
public:
	YPowerConstToken(const double& pow) : Token(false,  std::make_pair(std::make_pair(TOKEN_Y_POWER_CONST,0), pow)) { } ;
	virtual ~YPowerConstToken() { } ;
	
	virtual void eval(Context & context) const
	{
		double val = context.y ;
		int pow = (int)type.second-1 ;
		
		for(int i = 0 ; i < pow ; ++i)
		{
			val *= context.y;
		}
	    context.memory.push_back(val) ;
	}
	
	virtual std::string print() const
	{
		std::ostringstream s ;
		s << "< y^" << type.second << " >";
		return s.str() ;
	}
} ;

class ZPowerConstToken: public Token
{
public:
	ZPowerConstToken(const double& pow) : Token(false,  std::make_pair(std::make_pair(TOKEN_Z_POWER_CONST,0), pow)) { } ;
	virtual ~ZPowerConstToken() { } ;
	
	virtual void eval(Context & context) const
	{
		double val = context.z ;
		int pow = (int)type.second-1 ;
		
		for(int i = 0 ; i < pow ; ++i)
		{
			val *= context.z;
		}
	    context.memory.push_back(val) ;
	}
	
	virtual std::string print() const
	{
		std::ostringstream s ;
		s << "< z^" << type.second << " >";
		return s.str() ;
	}
} ;

class TPowerConstToken: public Token
{
public:
	TPowerConstToken(const double& pow) : Token(false,  std::make_pair(std::make_pair(TOKEN_T_POWER_CONST,0), pow)) { } ;
	virtual ~TPowerConstToken() { } ;
	
	virtual void eval(Context & context) const
	{
		double val = context.t ;
		int pow = (int)type.second-1 ;
		
		for(int i = 0 ; i < pow ; ++i)
		{
			val *= context.t;
		}
	    context.memory.push_back(val) ;
	}
	
	virtual std::string print() const
	{
		std::ostringstream s ;
		s << "< t^" << type.second << " >";
		return s.str() ;
	}
} ;

class ReadPowerConstToken: public Token
{
public:
	ReadPowerConstToken(const size_t adress, const double& pow) : Token(false,  std::make_pair(std::make_pair(TOKEN_READ_POWER_CONST,adress), pow)) { } ;
	virtual ~ReadPowerConstToken() { } ;
	
	virtual void eval(Context & context) const
	{
		double val = context.memory.heap[type.first.second] ;
		double v = val ;
		int pow = (int)type.second-1 ;
		
		for(int i = 0 ; i < pow ; ++i)
		{
			val *= v;
		}
	    context.memory.push_back(val) ;
	}
	
	virtual std::string print() const
	{
		std::ostringstream s ;
		s << "< readd@" << type.first.second <<" ^" << type.second << " >";
		return s.str() ;
	}
} ;

class AddReadAndConstToken: public Token
{
public:
	AddReadAndConstToken(const size_t adress, const double& v) : Token(false,  std::make_pair(std::make_pair(TOKEN_READ_ADD_CONST,adress), v)) { } ;
	virtual ~AddReadAndConstToken() { } ;
	
	virtual void eval(Context & context) const
	{
		context.memory.push_back(type.second+context.memory.heap[type.first.second]) ;
	}
	
	virtual std::string print() const
	{
		std::ostringstream s ;
		s << "< read@" << type.first.second <<"+" << type.second << " >";
		return s.str() ;
	}
} ;

class AddReadAndReadToken: public Token
{
protected:
	size_t add ;
public:
	AddReadAndReadToken(const size_t adress0, const size_t adress1) : Token(false,  std::make_pair(std::make_pair(TOKEN_READ_ADD_READ,adress0), 0)), add(adress1) { } ;
	virtual ~AddReadAndReadToken() { } ;
	
	virtual void eval(Context & context) const
	{
		context.memory.push_back(context.memory.heap[add]+context.memory.heap[type.first.second]) ;
	}
	
	virtual std::string print() const
	{
		std::ostringstream s ;
		s << "< read@" << type.first.second <<"+read@" << add << " >";
		return s.str() ;
	}
} ;

class MultiplyReadAndReadToken: public Token
{
protected:
	size_t add ;
public:
	MultiplyReadAndReadToken(const size_t adress0, const size_t adress1) : Token(false,  std::make_pair(std::make_pair(TOKEN_READ_ADD_READ,adress0), 0)), add(adress1) { } ;
	virtual ~MultiplyReadAndReadToken() { } ;
	
	virtual void eval(Context & context) const
	{
		context.memory.push_back(context.memory.heap[add]*context.memory.heap[type.first.second]) ;
	}
	
	virtual std::string print() const
	{
		std::ostringstream s ;
		s << "< read@" << type.first.second <<"+read@" << add << " >";
		return s.str() ;
	}
} ;

class AddXAndConstToken: public Token
{
public:
	AddXAndConstToken(const double& v) : Token(false,  std::make_pair(std::make_pair(TOKEN_ADD_X_AND_CONST,0), v)) { } ;
	virtual ~AddXAndConstToken() { } ;
	
	virtual void eval(Context & context) const
	{
		context.memory.push_back(type.second+context.x) ;
	}
	
	virtual std::string print() const
	{
		std::ostringstream s ;
		s << "< x+" << type.second << " >";
		return s.str() ;
	}
} ;


class AddYAndConstToken: public Token
{
public:
	AddYAndConstToken(const double& v) : Token(false,  std::make_pair(std::make_pair(TOKEN_ADD_Y_AND_CONST,0), v)) { } ;
	virtual ~AddYAndConstToken() { } ;
	
	virtual void eval(Context & context) const
	{
		context.memory.push_back(type.second+context.y) ;
	}
	
	virtual std::string print() const
	{
		std::ostringstream s ;
		s << "< y+" << type.second << " >";
		return s.str() ;
	}
} ;

class AddZAndConstToken: public Token
{
public:
	AddZAndConstToken(const double& v) : Token(false,  std::make_pair(std::make_pair(TOKEN_ADD_Z_AND_CONST,0), v)) { } ;
	virtual ~AddZAndConstToken() { } ;
	
	virtual void eval(Context & context) const
	{
		context.memory.push_back(type.second+context.z) ;
	}
	
	virtual std::string print() const
	{
		std::ostringstream s ;
		s << "< z+" << type.second << " >";
		return s.str() ;
	}
} ;

class AddTAndConstToken: public Token
{
public:
	AddTAndConstToken(const double& v) : Token(false,  std::make_pair(std::make_pair(TOKEN_ADD_T_AND_CONST,0), v)) { } ;
	virtual ~AddTAndConstToken() { } ;
	
	virtual void eval(Context & context) const
	{
		context.memory.push_back(type.second+context.t) ;
	}
	
	virtual std::string print() const
	{
		std::ostringstream s ;
		s << "< t+" << type.second << " >";
		return s.str() ;
	}
} ;
	

class MultiplyReadAndConstToken: public Token
{
public:
	MultiplyReadAndConstToken(const size_t adress, const double& v) : Token(false,  std::make_pair(std::make_pair(TOKEN_READ_MULTIPLY_CONST,adress), v)) { } ;
	virtual ~MultiplyReadAndConstToken() { } ;
	
	virtual void eval(Context & context) const
	{
		context.memory.push_back(type.second*context.memory.heap[type.first.second]) ;
	}
	
	virtual std::string print() const
	{
		std::ostringstream s ;
		s << "< read@" << type.first.second <<"*" << type.second << " >";
		return s.str() ;
	}
} ;

class MultiplyXAndConstToken: public Token
{
public:
	MultiplyXAndConstToken(const double& v) : Token(false,  std::make_pair(std::make_pair(TOKEN_MULTIPLY_X_AND_CONST,0), v)) { } ;
	virtual ~MultiplyXAndConstToken() { } ;
	
	virtual void eval(Context & context) const
	{
		context.memory.push_back(type.second*context.x) ;
	}
	
	virtual std::string print() const
	{
		std::ostringstream s ;
		s << "< x*" << type.second << " >";
		return s.str() ;
	}
} ;


class MultiplyYAndConstToken: public Token
{
public:
	MultiplyYAndConstToken(const double& v) : Token(false,  std::make_pair(std::make_pair(TOKEN_MULTIPLY_Y_AND_CONST,0), v)) { } ;
	virtual ~MultiplyYAndConstToken() { } ;
	
	virtual void eval(Context & context) const
	{
		context.memory.push_back(type.second*context.y) ;
	}
	
	virtual std::string print() const
	{
		std::ostringstream s ;
		s << "< y*" << type.second << " >";
		return s.str() ;
	}
} ;

class MultiplyZAndConstToken: public Token
{
public:
	MultiplyZAndConstToken(const double& v) : Token(false,  std::make_pair(std::make_pair(TOKEN_MULTIPLY_Z_AND_CONST,0), v)) { } ;
	virtual ~MultiplyZAndConstToken() { } ;
	
	virtual void eval(Context & context) const
	{
		context.memory.push_back(type.second*context.z) ;
	}
	
	virtual std::string print() const
	{
		std::ostringstream s ;
		s << "< z*" << type.second << " >";
		return s.str() ;
	}
} ;

class MultiplyTAndConstToken: public Token
{
public:
	MultiplyTAndConstToken(const double& v) : Token(false,  std::make_pair(std::make_pair(TOKEN_MULTIPLY_T_AND_CONST,0), v)) { } ;
	virtual ~MultiplyTAndConstToken() { } ;
	
	virtual void eval(Context & context) const
	{
		context.memory.push_back(type.second*context.t) ;
	}
	
	virtual std::string print() const
	{
		std::ostringstream s ;
		s << "< t*" << type.second << " >";
		return s.str() ;
	}
} ;

class SubstractXWithConstToken: public Token
{
public:
	SubstractXWithConstToken(const double& v) : Token(false,  std::make_pair(std::make_pair(TOKEN_SUBSTRACT_X_WITH_CONST,0), v)) { } ;
	virtual ~SubstractXWithConstToken() { } ;
	
	virtual void eval(Context & context) const
	{
		context.memory.push_back(context.x-type.second) ;
	}
	
	virtual std::string print() const
	{
		std::ostringstream s ;
		s << "< x-" << type.second << " >";
		return s.str() ;
	}
} ;


class SubstractYWithConstToken: public Token
{
public:
	SubstractYWithConstToken(const double& v) : Token(false,  std::make_pair(std::make_pair(TOKEN_SUBSTRACT_Y_WITH_CONST,0), v)) { } ;
	virtual ~SubstractYWithConstToken() { } ;
	
	virtual void eval(Context & context) const
	{
		context.memory.push_back(context.y-type.second) ;
	}
	
	virtual std::string print() const
	{
		std::ostringstream s ;
		s << "< y-" << type.second << " >";
		return s.str() ;
	}
} ;

class SubstractZWithConstToken: public Token
{
public:
	SubstractZWithConstToken(const double& v) : Token(false,  std::make_pair(std::make_pair(TOKEN_SUBSTRACT_Z_WITH_CONST,0), v)) { } ;
	virtual ~SubstractZWithConstToken() { } ;
	
	virtual void eval(Context & context) const
	{
		context.memory.push_back(context.z-type.second) ;
	}
	
	virtual std::string print() const
	{
		std::ostringstream s ;
		s << "< z-" << type.second << " >";
		return s.str() ;
	}
} ;

class SubstractTWithConstToken: public Token
{
public:
	SubstractTWithConstToken(const double& v) : Token(false,  std::make_pair(std::make_pair(TOKEN_SUBSTRACT_T_WITH_CONST,0), v)) { } ;
	virtual ~SubstractTWithConstToken() { } ;
	
	virtual void eval(Context & context) const
	{
		context.memory.push_back(context.t-type.second) ;
	}
	
	virtual std::string print() const
	{
		std::ostringstream s ;
		s << "< t-" << type.second << " >";
		return s.str() ;
	}
} ;

class SubstractConstWithXToken: public Token
{
public:
	SubstractConstWithXToken(const double& v) : Token(false,  std::make_pair(std::make_pair(TOKEN_SUBSTRACT_CONST_WITH_X,0), v)) { } ;
	virtual ~SubstractConstWithXToken() { } ;
	
	virtual void eval(Context & context) const
	{
		context.memory.push_back(type.second-context.x) ;
	}
	
	virtual std::string print() const
	{
		std::ostringstream s ;
		s << "< " << type.second << "-x >";
		return s.str() ;
	}
} ;


class SubstractConstWithYToken: public Token
{
public:
	SubstractConstWithYToken(const double& v) : Token(false,  std::make_pair(std::make_pair(TOKEN_SUBSTRACT_CONST_WITH_Y,0), v)) { } ;
	virtual ~SubstractConstWithYToken() { } ;
	
	virtual void eval(Context & context) const
	{
		context.memory.push_back(type.second-context.y) ;
	}
	
	virtual std::string print() const
	{
		std::ostringstream s ;
		s << "< " << type.second << "-y >";
		return s.str() ;
	}
} ;

class SubstractConstWithTToken: public Token
{
public:
	SubstractConstWithTToken(const double& v) : Token(false,  std::make_pair(std::make_pair(TOKEN_SUBSTRACT_CONST_WITH_T,0), v)) { } ;
	virtual ~SubstractConstWithTToken() { } ;
	
	virtual void eval(Context & context) const
	{
		context.memory.push_back(type.second-context.t) ;
	}
	
	virtual std::string print() const
	{
		std::ostringstream s ;
		s << "< " << type.second << "-t >";
		return s.str() ;
	}
} ;

class SubstractConstWithReadToken: public Token
{
public:
	SubstractConstWithReadToken(const size_t adress, const double& v) : Token(false,  std::make_pair(std::make_pair(TOKEN_CONST_SUBSTRACT_READ,adress), v)) { } ;
	virtual ~SubstractConstWithReadToken() { } ;
	
	virtual void eval(Context & context) const
	{
		context.memory.push_back(type.second-context.memory.heap[type.first.second]) ;
	}
	
	virtual std::string print() const
	{
		std::ostringstream s ;
		s << "< " << type.second << "-"<<"read@"<< type.first.second << " >";
		return s.str() ;
	}
} ;

class SubstractReadWithConstToken: public Token
{
public:
	SubstractReadWithConstToken(const size_t adress, const double& v) : Token(false,  std::make_pair(std::make_pair(TOKEN_READ_SUBSTRACT_CONST,adress), v)) { } ;
	virtual ~SubstractReadWithConstToken() { } ;
	
	virtual void eval(Context & context) const
	{
	    context.memory.push_back(context.memory.heap[type.first.second]-type.second) ;
	}
	
	virtual std::string print() const
	{
		std::ostringstream s ;
		s << "< read@" << type.first.second << "-"<< type.second<< " >";
		return s.str() ;
	}
} ;

class DivideConstWithReadToken: public Token
{
public:
	DivideConstWithReadToken(const size_t adress, const double& v) : Token(false,  std::make_pair(std::make_pair(TOKEN_CONST_DIVIDE_READ,adress), v)) { } ;
	virtual ~DivideConstWithReadToken() { } ;
	
	virtual void eval(Context & context) const
	{
		context.memory.push_back(type.second/context.memory.heap[type.first.second]) ;
	}
	
	virtual std::string print() const
	{
		std::ostringstream s ;
		s << "< " << type.second << "/"<<"read@"<< type.first.second << " >";
		return s.str() ;
	}
} ;

class DivideReadWithConstToken: public Token
{
public:
	DivideReadWithConstToken(const size_t adress, const double& v) : Token(false,  std::make_pair(std::make_pair(TOKEN_READ_DIVIDE_CONST,adress), v)) { } ;
	virtual ~DivideReadWithConstToken() { } ;
	
	virtual void eval(Context & context) const
	{
	    context.memory.push_back(context.memory.heap[type.first.second]/type.second) ;
	}
	
	virtual std::string print() const
	{
		std::ostringstream s ;
		s << "< read@" << type.first.second << "/"<< type.second<< " >";
		return s.str() ;
	}
} ;

class SubstractConstWithZToken: public Token
{
public:
	SubstractConstWithZToken(const double& v) : Token(false,  std::make_pair(std::make_pair(TOKEN_SUBSTRACT_CONST_WITH_Z,0), v)) { } ;
	virtual ~SubstractConstWithZToken() { } ;
	
	virtual void eval(Context & context) const
	{
		context.memory.push_back(type.second-context.z) ;
	}
	
	virtual std::string print() const
	{
		std::ostringstream s ;
		s << "< " << type.second << "-z >";
		return s.str() ;
	}
} ;


class DivideXWithConstToken: public Token
{
public:
	DivideXWithConstToken(const double& v) : Token(false,  std::make_pair(std::make_pair(TOKEN_DIVIDE_X_WITH_CONST,0), v)) { } ;
	virtual ~DivideXWithConstToken() { } ;
	
	virtual void eval(Context & context) const
	{
		context.memory.push_back(context.x/type.second) ;
	}
	virtual std::string print() const
	{
		std::ostringstream s ;
		s << "< x/" << type.second << " >";
		return s.str() ;
	}
} ;


class DivideYWithConstToken: public Token
{
public:
	DivideYWithConstToken(const double& v) : Token(false,  std::make_pair(std::make_pair(TOKEN_DIVIDE_Y_WITH_CONST,0), v)) { } ;
	virtual ~DivideYWithConstToken() { } ;
	
	virtual void eval(Context & context) const
	{
		context.memory.push_back(context.y/type.second) ;
	}
	virtual std::string print() const
	{
		std::ostringstream s ;
		s << "< y/" << type.second << " >";
		return s.str() ;
	}
} ;

class DivideZWithConstToken: public Token
{
public:
	DivideZWithConstToken(const double& v) : Token(false,  std::make_pair(std::make_pair(TOKEN_DIVIDE_Z_WITH_CONST,0), v)) { } ;
	virtual ~DivideZWithConstToken() { } ;
	
	virtual void eval(Context & context) const
	{
		context.memory.push_back(context.z/type.second) ;
	}
	virtual std::string print() const
	{
		std::ostringstream s ;
		s << "< z/" << type.second << " >";
		return s.str() ;
	}
} ;

class DivideTWithConstToken: public Token
{
public:
	DivideTWithConstToken(const double& v) : Token(false,  std::make_pair(std::make_pair(TOKEN_DIVIDE_T_WITH_CONST,0), v)) { } ;
	virtual ~DivideTWithConstToken() { } ;
	
	virtual void eval(Context & context) const
	{
		context.memory.push_back(context.t/type.second) ;
	}
	virtual std::string print() const
	{
		std::ostringstream s ;
		s << "< t/" << type.second << " >";
		return s.str() ;
	}
} ;

class DivideConstWithXToken: public Token
{
public:
	DivideConstWithXToken(const double& v) : Token(false,  std::make_pair(std::make_pair(TOKEN_DIVIDE_CONST_WITH_X,0), v)) { } ;
	virtual ~DivideConstWithXToken() { } ;
	
	virtual void eval(Context & context) const
	{
		context.memory.push_back(type.second/context.x) ;
	}
	virtual std::string print() const
	{
		std::ostringstream s ;
		s << "< " << type.second << "/x >";
		return s.str() ;
	}
} ;


class DivideConstWithYToken: public Token
{
public:
	DivideConstWithYToken(const double& v) : Token(false,  std::make_pair(std::make_pair(TOKEN_DIVIDE_CONST_WITH_Y,0), v)) { } ;
	virtual ~DivideConstWithYToken() { } ;
	
	virtual void eval(Context & context) const
	{
		context.memory.push_back(type.second/context.y) ;
	}
	
	virtual std::string print() const
	{
		std::ostringstream s ;
		s << "< " << type.second << "/y >";
		return s.str() ;
	}
} ;

class DivideConstWithZToken: public Token
{
public:
	DivideConstWithZToken(const double& v) : Token(false,  std::make_pair(std::make_pair(TOKEN_DIVIDE_CONST_WITH_Z,0), v)) { } ;
	virtual ~DivideConstWithZToken() { } ;
	
	virtual void eval(Context & context) const
	{
		context.memory.push_back(type.second/context.z) ;
	}
	
	virtual std::string print() const
	{
		std::ostringstream s ;
		s << "< " << type.second << "/z >";
		return s.str() ;
	}
} ;

class DivideConstWithTToken: public Token
{
public:
	DivideConstWithTToken(const double& v) : Token(false,  std::make_pair(std::make_pair(TOKEN_DIVIDE_CONST_WITH_T,0), v)) { } ;
	virtual ~DivideConstWithTToken() { } ;
	
	virtual void eval(Context & context) const
	{
		context.memory.push_back(type.second/context.t) ;
	}
	
	virtual std::string print() const
	{
		std::ostringstream s ;
		s << "< " << type.second << "/t >";
		return s.str() ;
	}
} ;


class XPowerAndMultiplyToken: public Token
{
public:
	XPowerAndMultiplyToken(const double& v, const size_t pow) : Token(false,  std::make_pair(std::make_pair(TOKEN_X_POWER_AND_MULTIPLY,pow), v)) { } ;
	virtual ~XPowerAndMultiplyToken() { } ;
	
	virtual void eval(Context & context) const
	{
		
		double val = context.x ;
		int pow = (int)type.first.second-1 ;
		
		for(int i = 0 ; i < pow ; ++i)
		{
			val *= context.x;
		}
	    context.memory.push_back(val*type.second) ;
	}
	
	virtual std::string print() const
	{
		std::ostringstream s ;
		s << "< x^" <<  (int)type.first.second << "*" << type.second << " >";
		return s.str() ;
	}
} ;

class YPowerAndMultiplyToken: public Token
{
public:
	YPowerAndMultiplyToken(const double& v, const size_t pow) : Token(false,  std::make_pair(std::make_pair(TOKEN_Y_POWER_AND_MULTIPLY,pow), v)) { } ;
	virtual ~YPowerAndMultiplyToken() { } ;
	
	virtual void eval(Context & context) const
	{
		
		double val = context.y ;
		int pow = (int)type.first.second-1 ;
		
		for(int i = 0 ; i < pow ; ++i)
		{
			val *= context.y;
		}
	    context.memory.push_back(val*type.second) ;
	}
	
	virtual std::string print() const
	{
		std::ostringstream s ;
		s << "< y^" <<  (int)type.first.second << "*" << type.second << " >";
		return s.str() ;
	}
} ;

class ZPowerAndMultiplyToken: public Token
{
public:
	ZPowerAndMultiplyToken(const double& v, const size_t pow) : Token(false,  std::make_pair(std::make_pair(TOKEN_Z_POWER_AND_MULTIPLY,pow), v)) { } ;
	virtual ~ZPowerAndMultiplyToken() { } ;
	
	virtual void eval(Context & context) const
	{
		
		double val = context.z ;
		int pow = (int)type.first.second-1 ;
		
		for(int i = 0 ; i < pow ; ++i)
		{
			val *= context.z;
		}
	    context.memory.push_back(val*type.second) ;
	}
	
	virtual std::string print() const
	{
		std::ostringstream s ;
		s << "< z^" <<  (int)type.first.second << "*" << type.second << " >";
		return s.str() ;
	}
} ;

class TPowerAndMultiplyToken: public Token
{
public:
	TPowerAndMultiplyToken(const double& v, const size_t pow) : Token(false,  std::make_pair(std::make_pair(TOKEN_T_POWER_AND_MULTIPLY,pow), v)) { } ;
	virtual ~TPowerAndMultiplyToken() { } ;
	
	virtual void eval(Context & context) const
	{
		
		double val = context.t ;
		int pow = (int)type.first.second-1 ;
		
		for(int i = 0 ; i < pow ; ++i)
		{
			val *= context.t;
		}
	    context.memory.push_back(val*type.second) ;
	}
	
	virtual std::string print() const
	{
		std::ostringstream s ;
		s << "< t^" <<  (int)type.first.second << "*" << type.second << " >";
		return s.str() ;
	}
} ;


class XYConstToken: public Token
{
public:
	XYConstToken(const double& v) : Token(false,  std::make_pair(std::make_pair(TOKEN_MULTIPLY_X_Y_AND_CONST,0), v)) { } ;
	virtual ~XYConstToken() { } ;
	
	virtual void eval(Context & context) const
	{
		context.memory.push_back(type.second*context.x*context.y) ;
	}
	
	virtual std::string print() const
	{
		std::ostringstream s ;
		s << "< x*y*" << type.second << " >";
		return s.str() ;
	}
} ;

class XZConstToken: public Token
{
public:
	XZConstToken(const double& v) : Token(false,  std::make_pair(std::make_pair(TOKEN_MULTIPLY_X_Z_AND_CONST,0), v)) { } ;
	virtual ~XZConstToken() { } ;
	
	virtual void eval(Context & context) const
	{
		context.memory.push_back(type.second*context.x*context.z) ;
	}
	
	virtual std::string print() const
	{
		std::ostringstream s ;
		s << "< x*z*" << type.second << " >";
		return s.str() ;
	}
} ;

class YZConstToken: public Token
{
public:
	YZConstToken(const double& v) : Token(false,  std::make_pair(std::make_pair(TOKEN_MULTIPLY_Y_Z_AND_CONST,0), v)) { } ;
	virtual ~YZConstToken() { } ;
	
	virtual void eval(Context & context) const
	{
		context.memory.push_back(type.second*context.y*context.z) ;
	}
	
	virtual std::string print() const
	{
		std::ostringstream s ;
		s << "< y*z*" << type.second << " >";
		return s.str() ;
	}
} ;

} ;









#endif
