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
#ifdef HAVE_TR1
#include <tr1/cmath>
#endif
#include "variable.h"
#include "../geometry/geometry_base.h"
#include "../geometry/geometry_2D.h"


double sign(const double t) ;
double positivity(const double t) ;
double negativity(const double t) ;
double interpolate(const double a, const double b) ;

namespace Mu
{

template<class ETYPE, class EABSTRACTTYPE>
struct Mesh ;

struct ElementarySurface ;
struct ElementaryVolume ;
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
	TOKEN_PROJECTION_OPERATOR,
	TOKEN_SIGNED_DISTANCE,
	TOKEN_2D_TRANSFORM,
	TOKEN_3D_TRANSFORM,
	TOKEN_POINT_DISTANCE_OPERATOR,
	TOKEN_POINT_DISTANCE_TRI_OPERATOR,
	TOKEN_ROTATION_OPERATOR,
	TOKEN_ANGLE_OPERATOR,
	TOKEN_POINT_SQUARE_DISTANCE_OPERATOR,
	TOKEN_LINE_OF_SIGHT_OPERATOR,
	TOKEN_CURVILINEAR_X_OPERATOR,
	TOKEN_CURVILINEAR_Y_OPERATOR,
	TOKEN_X,
	TOKEN_Y,
	TOKEN_Z,
	TOKEN_T,
	TOKEN_U,
	TOKEN_V,
	TOKEN_W,
	TOKEN_M_X,
	TOKEN_M_Y,
	TOKEN_M_Z,
	TOKEN_M_T,
	TOKEN_M_U,
	TOKEN_M_V,
	TOKEN_M_W,
	TOKEN_NAMED,
	TOKEN_UNARY_FUNCTION,
	TOKEN_BINARY_FUNCTION,
	TOKEN_PLUS,
	TOKEN_MINUS,
	TOKEN_TIMES,
	TOKEN_DIVIDES,
	TOKEN_POWER,
	TOKEN_EQUALS,
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
	TOKEN_BESSEL,
	TOKEN_ATAN2,
	TOKEN_INTERPOLATE,
	TOKEN_MULTIPLE_INTERPOLATE_FROM_TOP_2D,
	TOKEN_MULTIPLE_INTERPOLATE_FROM_BOTTOM_2D,
	TOKEN_MULTIPLE_INTERPOLATE_FROM_TOP_3D,
	TOKEN_MULTIPLE_INTERPOLATE_FROM_BOTTOM_3D,
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
	TOKEN_MESH_DISPLACEMENT_FIELD_X,
	TOKEN_MESH_DISPLACEMENT_FIELD_Y,
	TOKEN_MESH_DISPLACEMENT_FIELD_Z
} TokenTypeId ;

typedef std::pair<std::pair<TokenTypeId, short unsigned int>, double> TokenType ;

/** \brief Memory structure for the VirtualMachine. It provides a stack and a heap.
*/
struct Memory
{
	double stack[64];
	double heap[256] ;
	std::map<std::string, double *> variables ;
	std::vector<double *> variable_register ;
	std::string current_variable ;
	double * top_pos ;
	double * prev_top_pos ;

	/** \brief Constructor. Initialises stack and heap to 0.
	 * 
	 */
	Memory() : top_pos(&stack[0]), prev_top_pos(&stack[0]-1)
	{
	} ;
	
	~Memory() 
	{ 
		for(size_t i = 0 ; i < variable_register.size() ; i++)
			delete variable_register[i] ;
	} ;

	/** \brief Number of used positions in the stack
	 * 
	 * @return  Number of used positions in the stack
	 */
	size_t size() const
	{
		return top_pos-&stack[0]+1 ;
	}
	
	/** \brief Return the ith value stored in the stack
	 * 
	 * @param i index
	 * @return referenced value
	 */
	double & operator[](const size_t i)
	{
		return stack[i] ;
	}
		
	/** \brief Return the ith value stored in the stack
	 * 
	 * @param i index
	 * @return referenced value
	 */
	const double & operator[](const size_t i) const
	{
		return stack[i] ;
	}
	
	/** \brief Push a value on the stack
	 * 
	 * @param a values to push
	 */
	void push_back(double a)
	{
// 		if(top_pos > 30)
// 			std::cout << "stack overflow !" << std::endl ;
		++top_pos ;
		++prev_top_pos ;
		(*top_pos)= a ;
		
	}

	/** \brief pop the last value of the stack
	 * 
	 */
	void pop_back()
	{
		--top_pos ;
		--prev_top_pos ;
	}
	
	/** \brief reset stack to its original values and positions
	 * 
	 */
	void reset()
	{
		top_pos = &stack[0] ;
		prev_top_pos =  top_pos-1 ;
		for(std::map<std::string, double *>::iterator i = variables.begin(); i!= variables.end() ; i++)
			*i->second = 0 ;
	}
	
	void setNamedVariable(std::string name, double val)
	{
		if(variables.find(name) != variables.end() )
		{
			*variables[name] = val ;
		}
		else
		{
			variable_register.push_back(new double(val));
			variables[name] = variable_register.back() ;
		}
	}
} ;

/** \brief Context for the evaluation of a token. Contains a Memory and the values for the function arguments */
struct Context
{
	Memory memory; 
	double x ; 
	double y ; 
	double z ;
	double t ; 
	double u ; 
	double v ; 
	double w ;

	/** \brief Constructor, initialises the argument values
	 * 
	 * @param x_ 
	 * @param y_ 
	 * @param z_ 
	 * @param t_ 
	 * @param u_ 
	 * @param v_ 
	 * @param w_ 
	 */
	Context(const double & x_, const double & y_ , const double & z_ , const double & t_ , const double & u_ , const double & v_ , const double & w_ ) : x(x_), y(y_), z(z_), t(t_), u(u_), v(v_), w(w_) 
	{ 
		memory.variables["x"] = &x ;
		memory.variables["y"] = &y ;
		memory.variables["z"] = &z ;
		memory.variables["t"] = &t ;
		memory.variables["u"] = &u ;
		memory.variables["v"] = &v ;
		memory.variables["w"] = &w ;
	} ;

	/** \brief Constructor, initialises the argument values. arguments not set are 0
	 * 
	 * @param x_ 
	 * @param y_ 
	 * @param z_ 
	 * @param t_ 
	 * @param u_ 
	 * @param v_ 
	 */
	Context(const double & x_, const double & y_, const double & z_, const double & t_, const double & u_, const double & v_ ) : x(x_), y(y_), z(z_), t(t_), u(u_), v(v_), w(0) 
	{ 
		memory.variables["x"] = &x ;
		memory.variables["y"] = &y ;
		memory.variables["z"] = &z ;
		memory.variables["t"] = &t ;
		memory.variables["u"] = &u ;
		memory.variables["v"] = &v ;
		memory.variables["w"] = &w ;
	} ;

	/** \brief Constructor, initialises the argument values. arguments not set are 0
	 * 
	 * @param x_ 
	 * @param y_ 
	 * @param z_ 
	 * @param t_ 
	 * @param u_ 
	 */
	Context(const double & x_, const double & y_, const double & z_, const double & t_, const double & u_) : x(x_), y(y_), z(z_), t(t_), u(u_), v(0), w(0) 
	{ 
		memory.variables["x"] = &x ;
		memory.variables["y"] = &y ;
		memory.variables["z"] = &z ;
		memory.variables["t"] = &t ;
		memory.variables["u"] = &u ;
		memory.variables["v"] = &v ;
		memory.variables["w"] = &w ;
	} ;

	/** \brief Constructor, initialises the argument values. arguments not set are 0
	 * 
	 * @param x_ 
	 * @param y_ 
	 * @param z_ 
	 * @param t_ 
	 */
	Context(const double & x_, const double & y_, const double & z_, const double & t_) : x(x_), y(y_), z(z_), t(t_), u(0), v(0), w(0) 
	{ 
		memory.variables["x"] = &x ;
		memory.variables["y"] = &y ;
		memory.variables["z"] = &z ;
		memory.variables["t"] = &t ;
		memory.variables["u"] = &u ;
		memory.variables["v"] = &v ;
		memory.variables["w"] = &w ;
	} ;

	/** \brief Constructor, initialises the argument values. arguments not set are 0
	 * 
	 * @param x_ 
	 * @param y_ 
	 * @param z_ 
	 */
	Context(const double & x_, const double & y_, const double & z_) : x(x_), y(y_), z(z_), t(0), u(0), v(0), w(0) 
	{ 
		memory.variables["x"] = &x ;
		memory.variables["y"] = &y ;
		memory.variables["z"] = &z ;
		memory.variables["t"] = &t ;
		memory.variables["u"] = &u ;
		memory.variables["v"] = &v ;
		memory.variables["w"] = &w ;
	} ;

	/** \brief Constructor, initialises the argument values. arguments not set are 0
	 * 
	 * @param x_ 
	 * @param y_ 
	 */
	Context(const double & x_, const double & y_) : x(x_), y(y_), z(0), t(0), u(0), v(0), w(0) 
	{ 
		memory.variables["x"] = &x ;
		memory.variables["y"] = &y ;
		memory.variables["z"] = &z ;
		memory.variables["t"] = &t ;
		memory.variables["u"] = &u ;
		memory.variables["v"] = &v ;
		memory.variables["w"] = &w ;
	} ;

	/** \brief Constructor, initialises the argument values. arguments not set are 0
	 * 
	 * @param x_ 
	 */
	Context(const double & x_) : x(x_), y(0), z(0), t(0), u(0), v(0), w(0) 
	{ 
		memory.variables["x"] = &x ;
		memory.variables["y"] = &y ;
		memory.variables["z"] = &z ;
		memory.variables["t"] = &t ;
		memory.variables["u"] = &u ;
		memory.variables["v"] = &v ;
		memory.variables["w"] = &w ;
	} ;

	/** \brief Constructor, initialises the arguments to 0
	 * 
	 */
	Context() : x(0), y(0), z(0), t(0), u(0), v(0), w(0) 
	{ 
		memory.variables["x"] = &x ;
		memory.variables["y"] = &y ;
		memory.variables["z"] = &z ;
		memory.variables["t"] = &t ;
		memory.variables["u"] = &u ;
		memory.variables["v"] = &v ;
		memory.variables["w"] = &w ;
	} ;

	/** \brief Sets the argument values.
	 * 
	 * @param x_ 
	 * @param y_ 
	 * @param z_ 
	 * @param t_ 
	 * @param u_ 
	 * @param v_ 
	 * @param w_ 
	 */
	void set(const double & x_, const double & y_ , const double & z_ , const double & t_ , const double & u_ , const double & v_ , const double & w_ )
	{
		x = x_; 
		y = y_;
		z = z_;
		t = t_; 
		u = u_; 
		v = v_; 
		w = w_;
	}

	/** \brief Sets the argument values. Arguments not set are kept.
	 * 
	 * @param x_ 
	 * @param y_ 
	 * @param z_ 
	 * @param t_ 
	 * @param u_ 
	 * @param v_ 
	 * @param w_ 
	 */
	void set(const double & x_, const double & y_, const double & z_ )
	{
		x = x_; 
		y = y_;
		z = z_;
	}

	/** \brief Sets the argument values. Arguments not set are kept.
	 * 
	 * @param x_ 
	 * @param y_ 
	 */
	void set(const double & x_, const double & y_)
	{
		x = x_; 
		y = y_;
	}

	/** \brief Sets the argument values. Arguments not set are kept.
	 * 
	 * @param x_ 
	 */
	void set(const double & x_)
	{
		x = x_; 
	}
} ;

/** \brief Bytecode Token for the constitutions of Functions. A token operates on a stack and heap given arguments*/
class Token
{
public:
	
	
	/** \brief Constructor. The Type of the Token contains three values: the Enum defining the type, and two values which can be used for the evaluation.
	 * 
	 * @param null true is token represents putting 0 on the stack
	 * @param t type.
	 */
	Token(bool null = false, TokenType t =  std::make_pair(std::make_pair(TOKEN,0),(double)(0))) : isNull(null), type(t) { ; };
	
	virtual ~Token() { } ;
	
	virtual void eval(Context & context) const = 0;
	
	virtual std::string print() const  = 0 ;
	
	const bool isNull ;
	const TokenType type ;
} ;

/** \brief Put a constant on the stack */
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

/** \brief Put a 0 on the stack */
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

/** \brief Put -1 or 1 on the stack depending on whether we are within or without a domain delimited by a set of Segment s */
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
		
		*context.memory.top_pos =  (intersections & 1) * 2 - 1;
// 		if(intersections%2 == 1)
// 		    *context.memory.top_pos = -1 ;
// 		else
// 		    *context.memory.top_pos = 1 ;
	}
	
	virtual ~PositionToken();
	virtual std::string print() const
	{
		return std::string("position") ;
	}
} ;

/** \brief Put the distance of a point to an origin along a segmented line*/
class CurvilinearXOperatorToken : public Token
{
	const SegmentedLine * line ;
	bool origin ;
public:
	CurvilinearXOperatorToken(const SegmentedLine * l, bool fromHead) ;
	
	virtual void eval(Context & context) const;
	virtual std::string print() const;
	
} ;

/** \brief Put the height of a point to an origin along a segmented line*/
class CurvilinearYOperatorToken: public Token
{
	const SegmentedLine * line ;
	bool origin ;
public:
	CurvilinearYOperatorToken(const SegmentedLine * l, bool fromHead) ;
	
	virtual void eval(Context & context) const;
	virtual std::string print() const;
	
} ;

/** \brief Put the distance between a point defined by the two last positions on the stack and its projection on a Line on the stack. */
class LineDistanceOperatorToken : public Token
{
	Line l ;
public:
		LineDistanceOperatorToken(const Line & l_ ) : Token(false, std::make_pair(std::make_pair(TOKEN_POSITION_OPERATOR, 0), (double)(0))), l(l_)
	{
	}
	
	virtual void eval(Context & context) const
	{
		
		Point test(*context.memory.top_pos, *context.memory.prev_top_pos) ;
	    	context.memory.pop_back() ;
	
		
		*context.memory.top_pos =  sqrt(squareDist2D(test,l.projection(test)));
	}
	
	virtual ~LineDistanceOperatorToken() { };
	virtual std::string print() const
	{
		return std::string("distToLine") ;
	}
} ;

/** \brief Put the distance between a point defined by the arguments and its projection on a Line on the stack. */
class PositionOperatorToken : public Token
{
	std::vector<Segment> s ;
	Point w ;
public:
	PositionOperatorToken(Segment s_ ) ;
	
	PositionOperatorToken(std::vector<Segment> s_ ) ;
	
	virtual void eval(Context & context) const;
	
	virtual ~PositionOperatorToken();
	virtual std::string print() const;
} ;

/** \brief Put 1 on the stack if Point taken from context is in the given Geometry, -1 otherwise */
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

/** \brief Put -1 point defined by the two last positions on the stack lies outside a given Geometry, 1 otherwise. */
class DomainBinaryOperatorToken : public Token
{
	const Geometry* geo ;
public:
	DomainBinaryOperatorToken(const Geometry * g ) : Token(false, std::make_pair(std::make_pair(TOKEN_DOMAIN, 0), (double)(0))), geo(g)
	{
	}
	
	virtual void eval(Context & context) const
	{
		
		Point p(*context.memory.top_pos, *context.memory.prev_top_pos) ;
		context.memory.pop_back() ;
		if(geo->in(p))
			*context.memory.top_pos = 1 ;
		else
			*context.memory.top_pos = -1 ;
	}
	virtual ~DomainBinaryOperatorToken() { };
	virtual std::string print() const
	{
		return std::string("domainBinOp") ;
	}
} ;

/** \brief Put on the stack two values corresponding to the transformation defined by a surface element, using the arguments as the original position. */
class Transform2DToken : public Token
{

	const ElementarySurface * e ;
public:
	Transform2DToken(const ElementarySurface * g ) ;
	
	virtual void eval(Context & context) const;
	virtual ~Transform2DToken() ;
	virtual std::string print() const;
} ;

/** \brief Put on the stack three values corresponding to the transformation defined by a volume element, using the arguments as the original position. */
class Transform3DToken : public Token
{
	const ElementaryVolume * e ;
public:
	Transform3DToken(const ElementaryVolume * g ) ;
	
	virtual void eval(Context & context) const;
	virtual ~Transform3DToken() ;
	virtual std::string print() const;
} ;

/** \brief Put on the stack the distance between a point defined by the last positions on the stack and a stored position. */
class PointDistanceBinaryOperatorToken : public Token
{
	double x0 ;
	double y0 ;
public:
	PointDistanceBinaryOperatorToken(const Point & p ) : Token(false, std::make_pair(std::make_pair(TOKEN_POINT_DISTANCE_OPERATOR, 0), (double)(0))), x0(p.x), y0(p.y)
	{
	}
	
	virtual void eval(Context & context) const
	{
		double x = *context.memory.top_pos-x0 ;
		double y = *context.memory.prev_top_pos-y0 ;
		context.memory.pop_back() ;
		*context.memory.top_pos = sqrt(x*x+y*y) ;

	}
	virtual ~PointDistanceBinaryOperatorToken() { };
	virtual std::string print() const
	{
		return std::string("pointDistBinOp") ;
	}
} ;

/** \brief Put on the stack the distance between a point defined by the last positions on the stack and a stored position. */
class PointDistanceTrinaryOperatorToken : public Token
{
	double x0 ;
	double y0 ;
	double z0 ;
public:
	PointDistanceTrinaryOperatorToken(const Point & p ) : Token(false, std::make_pair(std::make_pair(TOKEN_POINT_DISTANCE_TRI_OPERATOR, 0), (double)(0))), x0(p.x), y0(p.y), z0(p.z)
	{
	}
	
	virtual void eval(Context & context) const
	{
		double x = *context.memory.top_pos-x0 ;
		double y = *context.memory.prev_top_pos-y0 ;
		double z = *(context.memory.prev_top_pos-1)-z0 ;
		context.memory.pop_back() ;
		context.memory.pop_back() ;
		*context.memory.top_pos = sqrt(x*x+y*y+z*z) ;

	}
	virtual ~PointDistanceTrinaryOperatorToken() { };
	virtual std::string print() const
	{
		return std::string("pointDistTrinOp") ;
	}
} ;

/** \brief return the x displacement if (x, y, z) lies in the mesh, else return 0 */
template<class ETYPE, class EABSTRACTTYPE>
class MeshXDisplacementToken : public Token
{
	Mesh<ETYPE, EABSTRACTTYPE> * mesh ;
public:
	MeshXDisplacementToken(Mesh<ETYPE, EABSTRACTTYPE> * mesh ) : Token(false, 
																		std::make_pair
																		(
																			std::make_pair
																			(
																				TOKEN_MESH_DISPLACEMENT_FIELD_X, 
																				0
																			), 
																			0.
																		)
																		), mesh(mesh)
	{ }
	
	virtual void eval(Context & context) const
	{
		Point where(context.x, context.y, context.z) ;
		std::vector<ETYPE *> elements  = mesh->getConflictingElements(&where) ;
		ETYPE * target = NULL;
		for(size_t i = 0 ; i < elements.size() ;  i++)
		{
			if(elements[i]->in(where))
			{
				target = elements[i] ;
				break ;
			}
			
		}
		double x = target->getState().getDisplacements(where)[0] ;
		
		context.memory.push_back(x) ;

	}
	virtual ~MeshXDisplacementToken() { };
	virtual std::string print() const
	{
		return std::string("mesh x displacement field") ;
	}
} ;

/** \brief return the y displacement if (x, y, z) lies in the mesh, else return 0 */
template<class ETYPE, class EABSTRACTTYPE>
class MeshYDisplacementToken : public Token
{
	Mesh<ETYPE, EABSTRACTTYPE> * mesh ;
public:
	MeshYDisplacementToken(Mesh<ETYPE, EABSTRACTTYPE> * mesh ) : Token(false, std::make_pair(std::make_pair(TOKEN_MESH_DISPLACEMENT_FIELD_Y, 0), 0.)), mesh(mesh)
	{ }
	
	virtual void eval(Context & context) const
	{
		Point where(context.x, context.y, context.z) ;
		std::vector<ETYPE *> elements  = mesh->getConflictingElements(&where) ;
		ETYPE * target = NULL;
		for(size_t i = 0 ; i < elements.size() ;  i++)
		{
			if(elements[i]->in(where))
			{
				target = elements[i] ;
				break ;
			}
			
		}
		double x = target->getState().getDisplacements(where)[1] ;
		
		context.memory.push_back(x) ;

	}
	virtual ~MeshYDisplacementToken() { };
	virtual std::string print() const
	{
		return std::string("mesh y displacement field") ;
	}
} ;

/** \brief return the y displacement if (x, y, z) lies in the mesh, else return 0 */
template<class ETYPE, class EABSTRACTTYPE>
class MeshZDisplacementToken : public Token
{
	Mesh<ETYPE, EABSTRACTTYPE> * mesh ;
public:
	MeshZDisplacementToken(Mesh<ETYPE, EABSTRACTTYPE> * mesh ) : Token(false, std::make_pair(std::make_pair(TOKEN_MESH_DISPLACEMENT_FIELD_Z, 0), 0.)), mesh(mesh)
	{ }
	
	virtual void eval(Context & context) const
	{
		Point where(context.x, context.y, context.z) ;
		std::vector<ETYPE *> elements  = mesh->getConflictingElements(&where) ;
		ETYPE * target = NULL;
		for(size_t i = 0 ; i < elements.size() ;  i++)
		{
			if(elements[i]->in(where))
			{
				target = elements[i] ;
				break ;
			}
			
		}
		double x = target->getState().getDisplacements(where)[2] ;
		
		context.memory.push_back(x) ;

	}
	virtual ~MeshZDisplacementToken() { };
	virtual std::string print() const
	{
		return std::string("mesh z displacement field") ;
	}
} ;

/** \brief Put on the stack the two coordinates of a transformed point defined by the last positions on the stack given a rotation. */
class RotationBinaryOperatorToken : public Token
{
	double cangle ;
	double sangle ;
public:
	RotationBinaryOperatorToken(double a ) : Token(false, std::make_pair(std::make_pair(TOKEN_ROTATION_OPERATOR, 0), (double)(0))), cangle(cos(a)), sangle(sin(a))
	{
	}
	
	virtual void eval(Context & context) const
	{
		
		double x = *context.memory.top_pos ;
		double y =  *context.memory.prev_top_pos ;
		*context.memory.top_pos = x*cangle + y*sangle ;
		*context.memory.prev_top_pos = -x*sangle + y*cangle ;

	}
	virtual ~RotationBinaryOperatorToken() { };
	virtual std::string print() const
	{
		return std::string("rotate") ;
	}
} ;

/** \brief Put on the stack the angle between a point defined by the last positions on the stack and a stored position. */
class AngleBinaryOperatorToken : public Token
{
	double cangle ;
	double sangle ;
	Point pivot ;
public:
	AngleBinaryOperatorToken(double a, const Point & p ) : Token(false, std::make_pair(std::make_pair(TOKEN_ANGLE_OPERATOR, 0), (double)(0))), cangle(cos(a)), sangle(sin(a)), pivot(p.x*cos(a)+p.y*sin(a), -p.x*sin(a)+p.y*cos(a))
	{
	}
	
	virtual void eval(Context & context) const
	{
		
		double x = *context.memory.top_pos ;
		double y =  *context.memory.prev_top_pos ;
		double x_t = x*cangle + y*sangle ;
		double y_t = -x*sangle + y*cangle ;
		context.memory.pop_back() ;
		*context.memory.top_pos = atan2(y_t-pivot.y, x_t-pivot.x) ;

	}
	virtual ~AngleBinaryOperatorToken() { };
	virtual std::string print() const
	{
		return std::string("angle") ;
	}
} ;

/** \brief Put on the stack squared distance between a point defined by the last two positions on the stack and a tored position */
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
		
		Point p(*context.memory.top_pos, *context.memory.prev_top_pos) ;
		context.memory.pop_back() ;

		*context.memory.top_pos = squareDist2D(p, base) ;

	}
	virtual ~PointSquareDistanceBinaryOperatorToken() { };
	virtual std::string print() const
	{
		return std::string("pointSqDistBinOp") ;
	}
} ;

/** \brief Put 1 one the stack if a point defined by the two last positions on the stack is not visible from a stored position given an obstructing Geometry, 0 otherwise. */
class LineOfSightOperatorToken : public Token
{
	Point base ;
	const Geometry * obstruction ;
public:
	LineOfSightOperatorToken(const Point & p,  const Geometry * o) : Token(false, std::make_pair(std::make_pair(TOKEN_LINE_OF_SIGHT_OPERATOR, 0), (double)(0))), base(p), obstruction(o)
	{
	}
	
	virtual void eval(Context & context) const
	{
		
		Point p(*context.memory.top_pos, *context.memory.prev_top_pos) ;
		

		*context.memory.top_pos = Segment(base, p).intersects(obstruction) ;

	}
	virtual ~LineOfSightOperatorToken() { };
	virtual std::string print() const
	{
		return std::string("lineOfSightOp") ;
	}
} ;

/** \brief put on the stack a value stored in the heap. */
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

/** \brief set a heap value */
class SetHeapVariableToken : public Token
{

public:
	SetHeapVariableToken(const size_t a) : Token(false, std::make_pair(std::make_pair(TOKEN_WRITE_VARIABLE, a),(double)(0))) { } ;

	
	virtual void eval(Context & context) const
	{
	    context.memory.heap[type.first.second] = *context.memory.top_pos ;
	    context.memory.pop_back() ;
	}
	
	virtual std::string print() const
	{
		std::ostringstream s ;
		s << "write@" << type.first.second;
		return s.str() ;
	}
} ;

/** \brief put on the stack the distance between a point defined by the arguments and its projection on a stored Segment */
class ProjectionToken : public Token
{
	Segment s ;
public:
	ProjectionToken(Segment s_ ) : Token(false, std::make_pair(std::make_pair(TOKEN_PROJECTION, 0),(double)(0))), s(s_)
	{
	}
	
	virtual void eval(Context & context) const
	{
		
	    context.memory.push_back(dist(Point(context.x, context.y, context.z), s.project(Point(context.x, context.y, context.z)))) ;
	}
	virtual ~ProjectionToken();
	virtual std::string print() const
	{
		return std::string("check position") ;
	}
} ;

/** \brief put on the stack the distance between a point defined by the two last positions on the Stack and its projection on a stored Segment */
class ProjectionBinaryOperatorToken : public Token
{
	const Geometry * g ;
public:
	ProjectionBinaryOperatorToken(const Geometry * s_ ) : Token(false, std::make_pair(std::make_pair(TOKEN_PROJECTION_OPERATOR, 0),(double)(0))), g(s_)
	{
	}
	
	virtual void eval(Context & context) const
	{
		Point p(*context.memory.top_pos, *context.memory.prev_top_pos) ;
		Point p_(p) ;
		g->project(&p_) ;
		context.memory.pop_back() ;

		*context.memory.top_pos = sqrt(squareDist2D(p, p_)) ;

	}
	virtual ~ProjectionBinaryOperatorToken() { };
	virtual std::string print() const
	{
		return std::string("projection") ;
	}
} ;

/** \brief put on the stack the value of the x argument */
class XToken : public Token
{
public:
	XToken() : Token(false, std::make_pair(std::make_pair(TOKEN_X, 0),(double)(0))) { };
	virtual void eval(Context & context) const {context.memory.push_back(context.x) ; context.memory.current_variable = std::string("x") ;}
	virtual ~XToken(){ };
	virtual std::string print() const { return std::string("x") ;}
} ;

/** \brief put on the stack minus the value of the x argument */

class XMToken : public Token
{
public:
	XMToken() : Token(false, std::make_pair(std::make_pair(TOKEN_M_X, 0),(double)(0))) { };
	virtual void eval(Context & context) const {context.memory.push_back(-context.x) ;context.memory.current_variable = std::string("x") ;}
	virtual ~XMToken(){ };
	virtual std::string print() const { return std::string("-x") ;}
} ;

/** \brief put on the stack the value of the y argument */
class YToken : public Token
{
public:
	YToken() : Token(false, std::make_pair(std::make_pair(TOKEN_Y, 0),(double)(0))) { };
	virtual void eval(Context & context) const {context.memory.push_back(context.y) ;context.memory.current_variable = std::string("y") ;}
	virtual ~YToken(){ };
	virtual std::string print() const { return std::string("y") ;}
} ;

/** \brief put on the stack minus the value of the y argument */
class YMToken : public Token
{
public:
	YMToken() : Token(false, std::make_pair(std::make_pair(TOKEN_M_Y, 0),(double)(0))) { };
	virtual void eval(Context & context) const {context.memory.push_back(-context.y) ;context.memory.current_variable = std::string("y") ;}
	virtual ~YMToken(){ };
	virtual std::string print() const { return std::string("-y") ;}
} ;

/** \brief put on the stack the value of the z argument */
class ZToken : public Token
{
public:
	ZToken() : Token(false, std::make_pair(std::make_pair(TOKEN_Z, 0),(double)(0))) { };
	virtual void eval(Context & context) const {context.memory.push_back(context.z) ;context.memory.current_variable = std::string("z") ;}
	virtual ~ZToken(){ };
	virtual std::string print() const { return std::string("z") ;}
} ;

/** \brief put on the stack minus the value of the z argument */
class ZMToken : public Token
{
public:
	ZMToken() : Token(false, std::make_pair(std::make_pair(TOKEN_M_Z, 0),(double)(0))) { };
	virtual void eval(Context & context) const {context.memory.push_back(-context.z) ;context.memory.current_variable = std::string("z") ;}
	virtual ~ZMToken(){ };
	virtual std::string print() const { return std::string("-z") ;}
} ;

/** \brief put on the stack the value of the t argument */
class TToken : public Token
{
public:
	TToken() : Token(false, std::make_pair(std::make_pair(TOKEN_T, 0),(double)(0))) { };
	virtual void eval(Context & context) const {context.memory.push_back(context.t) ;context.memory.current_variable = std::string("t") ;}
	virtual ~TToken(){ };
	virtual std::string print() const { return std::string("t") ;}
} ;

/** \brief put on the stack minus the value of the t argument */
class TMToken : public Token
{
public:
	TMToken() : Token(false, std::make_pair(std::make_pair(TOKEN_M_T, 0),(double)(0))) { };
	virtual void eval(Context & context) const {context.memory.push_back(-context.t) ;context.memory.current_variable = std::string("t") ;}
	virtual ~TMToken(){ };
	virtual std::string print() const { return std::string("-t") ;}
} ;

/** \brief put on the stack the value of the u argument */
class UToken : public Token
{
public:
	UToken() : Token(false, std::make_pair(std::make_pair(TOKEN_U, 0),(double)(0))) { };
	virtual void eval(Context & context) const {context.memory.push_back(context.u) ;context.memory.current_variable = std::string("u") ;}
	virtual ~UToken(){ };
	virtual std::string print() const { return std::string("u") ;}
} ;


/** \brief put on the stack minus the value of the u argument */
class UMToken : public Token
{
public:
	UMToken() : Token(false, std::make_pair(std::make_pair(TOKEN_M_U, 0),(double)(0))) { };
	virtual void eval(Context & context) const {context.memory.push_back(-context.u) ;context.memory.current_variable = std::string("u") ;}
	virtual ~UMToken(){ };
	virtual std::string print() const { return std::string("-u") ;}
} ;

/** \brief put on the stack the value of the v argument */
class VToken : public Token
{
public:
	VToken() : Token(false, std::make_pair(std::make_pair(TOKEN_V, 0),(double)(0))) { };
	virtual void eval(Context & context) const {context.memory.push_back(context.v) ;context.memory.current_variable = std::string("v") ;}
	virtual ~VToken(){ };
	virtual std::string print() const { return std::string("v") ;}
} ;

/** \brief put on the stack minus the value of the v argument */
class VMToken : public Token
{
public:
	VMToken() : Token(false, std::make_pair(std::make_pair(TOKEN_M_V, 0),(double)(0))) { };
	virtual void eval(Context & context) const {context.memory.push_back(-context.v) ;context.memory.current_variable = std::string("v") ;}
	virtual ~VMToken(){ };
	virtual std::string print() const { return std::string("-v") ;}
} ;

/** \brief put on the stack the value of the w argument */
class WToken : public Token
{
public:
	WToken() : Token(false, std::make_pair(std::make_pair(TOKEN_W, 0),(double)(0))) { };
	virtual void eval(Context & context) const {context.memory.push_back(context.w) ;context.memory.current_variable = std::string("w") ;}
	virtual ~WToken(){ };
	virtual std::string print() const { return std::string("w") ;}
} ;

/** \brief put on the stack minus the value of the e argument */
class WMToken : public Token
{
public:
	WMToken() : Token(false, std::make_pair(std::make_pair(TOKEN_M_W, 0),(double)(0))) { };
	virtual void eval(Context & context) const {context.memory.push_back(-context.w) ;context.memory.current_variable = std::string("w") ;}
	virtual ~WMToken(){ };
	virtual std::string print() const { return std::string("-w") ;}
} ;

class NamedToken :public Token
{
public:
	std::string name ;
public:
	NamedToken(const std::string & name ) : Token(false, std::make_pair(std::make_pair(TOKEN_NAMED, 0),(double)(0))), name(name) { };
	virtual void eval(Context & context) const 
	{

		if(context.memory.variables.find(name) == context.memory.variables.end())
		{
			context.memory.variable_register.push_back(new double(0)) ;
			context.memory.variables[name] = context.memory.variable_register.back() ;
		}
		
		context.memory.push_back(*context.memory.variables[name]) ;
		context.memory.current_variable = name ;
		
	}
	virtual ~NamedToken(){ };
	virtual std::string print() const { return name ;}
} ;

/** \brief Apply unary function on the last value of the stack */
class UnaryFunctionToken : public Token
{
protected:
	const unaryFunctionPointer fctPtr ;
	
public:
	UnaryFunctionToken(unaryFunctionPointer f, bool nul = false) ;
	
	virtual void eval(Context & context) const
	{
	    *context.memory.top_pos = fctPtr(*context.memory.top_pos) ;
	}
	
	virtual ~UnaryFunctionToken() ;
	
	virtual std::string print() const
	{
		return std::string("1-aryFct")  ;
	}
} ;

/** \brief Apply binary function on the two last values of the stack */
class BinaryFunctionToken : public Token
{
protected:
	const binaryFunctionPointer fctPtr ;
	
public:
	BinaryFunctionToken(binaryFunctionPointer f, bool nul = false);
	
	virtual void eval(Context & context) const
	{
		double new_val_0 = *context.memory.prev_top_pos ;
		double new_val_1 = *context.memory.top_pos ;
	    context.memory.pop_back() ;
	    *context.memory.top_pos = fctPtr(new_val_0, new_val_1) ;
	}
	
	virtual ~BinaryFunctionToken() ;
	
	virtual std::string print() const
	{
		return std::string("2-aryFct")  ;
	}
} ;

/** \brief Take the cos of the last value of the stack */
class CosToken : public Token
{
public:
	CosToken() ;
	
	virtual void eval(Context & context) const ;
	
	virtual ~CosToken() ;
	
	virtual std::string print() const ;
} ;

/** \brief Take the abs of the last value of the stack */
class AbsToken : public Token
{
public:
	AbsToken() : Token(false, std::make_pair(std::make_pair(TOKEN_ABS, 0),(double)(0))) { };
	
	virtual void eval(Context & context) const
	{
		*context.memory.top_pos = std::abs(*context.memory.top_pos) ;
	}
	
	virtual ~AbsToken() {} ;
	
	virtual std::string print() const
	{
		return std::string("abs")  ;
	}
} ;

/** \brief Take the tan of the last value of the stack */
class TanToken : public Token
{
public:
	TanToken() : Token(false, std::make_pair(std::make_pair(TOKEN_TAN, 0),(double)(0))) { };
	
	virtual void eval(Context & context) const
	{
	    *context.memory.top_pos = tan(*context.memory.top_pos) ;
	}
	
	virtual ~TanToken() {} ;
	
	virtual std::string print() const
	{
		return std::string("tan")  ;
	}
} ;

/** \brief Take the sin of the last value of the stack */
class SinToken : public Token
{
public:
	SinToken() ;
	
	virtual void eval(Context & context) const ;
	
	virtual ~SinToken()  ;
	
	virtual std::string print() const ;
} ;

/** \brief Take the exp of the last value of the stack */
class ExpToken : public Token
{
public:
	ExpToken() : Token(false, std::make_pair(std::make_pair(TOKEN_EXP, 0),(double)(0)) ) { };
	
	virtual void eval(Context & context) const
	{
	    *context.memory.top_pos = exp(*context.memory.top_pos) ;
	}
	
	virtual ~ExpToken() {} ;
	
	virtual std::string print() const
	{
		return std::string("exp")  ;
	}
} ;

/** \brief Take the sign of the last value of the stack */
class SignFunctionToken : public Token
{
public:
	SignFunctionToken() : Token(false, std::make_pair(std::make_pair(TOKEN_SIGN, 0),(double)(0)) ) { };
	
	virtual void eval(Context & context) const
	{
	    *context.memory.top_pos = sign(*context.memory.top_pos) ;
	}
	
	virtual ~SignFunctionToken() {} ;
	
	virtual std::string print() const
	{
		return std::string("sign")  ;
	}
} ;

/** \brief Put one if the last value of the stack is positive, 0 otherwise*/
class PositivityFunctionToken : public Token
{
public:
	PositivityFunctionToken() : Token(false, std::make_pair(std::make_pair(TOKEN_SIGN, 0),(double)(0)) ) { };
	
	virtual void eval(Context & context) const
	{
		double s = sign(*context.memory.top_pos) ;
		if(s > 0)
			*context.memory.top_pos = 1 ;
		else
			*context.memory.top_pos = 0 ;
	}
	
	virtual ~PositivityFunctionToken() {} ;
	
	virtual std::string print() const
	{
		return std::string(">0")  ;
	}
} ;

/** \brief Put one if the last value of the stack is negative, 0 otherwise*/
class NegativityFunctionToken : public Token
{
public:
	NegativityFunctionToken() : Token(false, std::make_pair(std::make_pair(TOKEN_SIGN, 0),(double)(0)) ) { };
	
	virtual void eval(Context & context) const
	{
		double s = sign(*context.memory.top_pos) ;
		if(s < 0)
	    	*context.memory.top_pos = 1 ;
		else
			*context.memory.top_pos = 0 ;
	}
	
	virtual ~NegativityFunctionToken() {} ;
	
	virtual std::string print() const
	{
		return std::string("<0")  ;
	}
} ;

/** \brief Take the log of the last value of the stack */
class LogToken : public Token
{
public:
	LogToken() : Token(false, std::make_pair(std::make_pair(TOKEN_LOG, 0),(double)(0)))  { };
	
	virtual void eval(Context & context) const
	{
	    *context.memory.top_pos = log(*context.memory.top_pos) ;
	}
	
	virtual ~LogToken() {} ;
	
	virtual std::string print() const
	{
		return std::string("log")  ;
	}
} ;

/** \brief Take the cosh of the last value of the stack */
class CoshToken : public Token
{
public:
	CoshToken() : Token(false, std::make_pair(std::make_pair(TOKEN_COSH, 0),(double)(0)) ) { };
	
	virtual void eval(Context & context) const
	{
	    *context.memory.top_pos = cosh(*context.memory.top_pos) ;
	}
	
	virtual ~CoshToken() {} ;
	
	virtual std::string print() const
	{
		return std::string("cosh")  ;
	}
} ;

/** \brief Take the sinh of the last value of the stack */
class SinhToken : public Token
{
public:
	SinhToken() : Token(false, std::make_pair(std::make_pair(TOKEN_SINH, 0),(double)(0)) ) { };
	
	virtual void eval(Context & context) const
	{
	    *context.memory.top_pos = sinh(*context.memory.top_pos) ;
	}
	
	virtual ~SinhToken() {} ;
	
	virtual std::string print() const
	{
		return std::string("sinh")  ;
	}
} ;

/** \brief Take the tanh of the last value of the stack */
class TanhToken : public Token
{
public:
	TanhToken() : Token(false, std::make_pair(std::make_pair(TOKEN_TANH, 0),(double)(0)) ) { };
	
	virtual void eval(Context & context) const
	{
	    *context.memory.top_pos = tanhl(*context.memory.top_pos) ;
	}
	
	virtual ~TanhToken() {} ;
	
	virtual std::string print() const
	{
		return std::string("tanh")  ;
	}
} ;

/** \brief Take the sqrt of the last value of the stack */
class SqrtToken : public Token
{
public:
	SqrtToken() : Token(false, std::make_pair(std::make_pair(TOKEN_SQRT, 0),(double)(0)))  { };
	
	virtual void eval(Context & context) const
	{
	    *context.memory.top_pos = sqrt(*context.memory.top_pos) ;
	}
	
	virtual ~SqrtToken() {} ;
	
	virtual std::string print() const
	{
		return std::string("sqrt")  ;
	}
} ;

/** \brief Take the Bessel function of the last value of the stack */
class BesselToken : public Token
{
public:
	BesselToken(int i) : Token(false, std::make_pair(std::make_pair(TOKEN_BESSEL, i),(double)(0))) { };
	
	virtual void eval(Context & context) const
	{
		if(std::abs(*context.memory.top_pos) > 1e-8)
#ifdef HAVE_TR1
			*context.memory.top_pos = std::tr1::cyl_bessel_j(type.first.second, *context.memory.top_pos) ;
#else
			*context.memory.top_pos = 0 ;
#endif
		else
			*context.memory.top_pos = (type.first.second == 0) ;
	}
	
	virtual ~BesselToken() {} ;
	
	virtual std::string print() const
	{
		return std::string("bessel")  ;
	}
} ;

/** \brief Take the atan2 function of the last two values of the stack */
class Atan2Token : public Token
{
public:
	Atan2Token() ;
	
	virtual void eval(Context & context) const ;
	
	virtual ~Atan2Token() ;
	
	virtual std::string print() const ;
} ;

/** \brief Interpolate between the last two values on the stack given the x argument */
class InterpolationToken : public Token
{
public:
	InterpolationToken() : Token(false, std::make_pair(std::make_pair(TOKEN_INTERPOLATE, 0),(double)(0)))  { };
	
	virtual void eval(Context & context) const
	{
		double new_val = interpolate(*context.memory.top_pos, *context.memory.prev_top_pos) ;
	    context.memory.pop_back() ;
	    *context.memory.top_pos = new_val ;
	}
	
	virtual ~InterpolationToken() {} ;
	
	virtual std::string print() const
	{
		return std::string("interpolate")  ;
	}
} ;

/** \brief Return the normalised distance between n geometries, starting from outside. Use the last two values on the stack.*/
class MultipleInterpolationFromTopToken2D : public Token
{
	std::vector<Geometry *> geos ;
public:
	MultipleInterpolationFromTopToken2D() : Token(false, std::make_pair(std::make_pair(TOKEN_MULTIPLE_INTERPOLATE_FROM_TOP_2D, 0),(double)(0)))  { };
	
	virtual void eval(Context & context) const
	{
		Point test(*context.memory.top_pos, *context.memory.prev_top_pos) ;
		context.memory.pop_back() ;
		Geometry * lastIn = NULL ;
		Geometry * previousIn = NULL ;
		Geometry * lastOut = NULL ;
		Geometry * previousOut = NULL ;
		for(std::vector<Geometry *>::const_iterator i = geos.begin() ; i != geos.end() ; ++i)
		{
			if((*i)->in(test))
			{
				if(lastIn)
					previousIn = lastIn ;
				lastIn = *i ;
			}
			else
			{
				if(lastOut)
					previousOut = lastOut ;
				lastOut = *i ;
			}
		}
		
		if(lastOut == NULL) // I am in all geometries.
		{
			Point p0(test) ;
			lastIn->project(&p0) ;
			*context.memory.top_pos = sqrt(squareDist2D(test, p0)) ;
			
		}
		else if (lastIn != NULL) // I am in a geometry, but not in the next
		{
			Point p0(test) ;
			lastIn->project(&p0) ; double din = sqrt(squareDist2D(test, p0)) ;
			Point p1(test) ;
			lastOut->project(&p1) ; double dout = sqrt(squareDist2D(test, p1)) ;
			*context.memory.top_pos = interpolate(din, dout) ;
			
		}
		else if(previousOut != NULL) // I am between 2 geometries
		{
			Point p0(test) ;
			previousOut->project(&p0) ; double din = sqrt(squareDist2D(test, p0)) ;
			Point p1(test) ;
			lastOut->project(&p1) ; double dout = sqrt(squareDist2D(test, p1)) ;
			*context.memory.top_pos = interpolate(din, dout) ;
		}
		else if(!geos.empty()) // I am outside the only geometry
		{
			Point p0(test) ;
			lastOut->project(&p0) ;
			*context.memory.top_pos = sqrt(squareDist2D(test, p0)) ;
		}
		else // there is nothing
		{
			*context.memory.top_pos = 0 ;
		}
	}
	
	virtual ~MultipleInterpolationFromTopToken2D() {} ;
	
	virtual std::string print() const
	{
		return std::string("multiple interpolate from top")  ;
	}
} ;

/** \brief Return the normalised distance between n geometries, starting from inside, use the last two values on the stack.*/
class MultipleInterpolationFromBottomToken2D : public Token
{
	std::vector<Geometry *> geos ;
public:
	MultipleInterpolationFromBottomToken2D() : Token(false, std::make_pair(std::make_pair(TOKEN_MULTIPLE_INTERPOLATE_FROM_BOTTOM_2D, 0),(double)(0)))  { };
	
	virtual void eval(Context & context) const
	{
		Point test(*context.memory.top_pos, *context.memory.prev_top_pos);
		context.memory.pop_back() ;
		Geometry * lastIn = NULL ;
		Geometry * previousIn = NULL ;
		Geometry * lastOut = NULL ;
		Geometry * previousOut = NULL ;
		for(std::vector<Geometry *>::const_reverse_iterator i = geos.rbegin() ; i != geos.rend() ; ++i)
		{
			if((*i)->in(test))
			{
				if(lastIn)
					previousIn = lastIn ;
				lastIn = *i ;
			}
			else
			{
				if(lastOut)
					previousOut = lastOut ;
				lastOut = *i ;
			}
		}
		
		if(lastOut == NULL) // I am in all geometries.
		{
			Point p0(test) ;
			lastIn->project(&p0) ;
			
			*context.memory.top_pos = sqrt(squareDist2D(test, p0)) ;
			
		}
		else if (lastIn != NULL) // I am in a geometry, but not in the next
		{
			Point p0(test) ;
			lastIn->project(&p0) ; double din = sqrt(squareDist2D(test, p0)) ;
			Point p1(test) ;
			lastOut->project(&p1) ; double dout = sqrt(squareDist2D(test, p1)) ;
			*context.memory.top_pos = interpolate(din, dout) ;
			
		}
		else if(previousOut != NULL) // I am between 2 geometries
		{
			Point p0(test) ;
			previousOut->project(&p0) ; double din = sqrt(squareDist2D(test, p0)) ;
			Point p1(test) ;
			lastOut->project(&p1) ; double dout = sqrt(squareDist2D(test, p1)) ;
			*context.memory.top_pos = interpolate(din, dout) ;
		}
		else if(!geos.empty()) // I am outside the only geometry
		{
			Point p0(test) ;
			lastOut->project(&p0) ;
			*context.memory.top_pos = sqrt(squareDist2D(test, p0)) ;
		}
		else // there is nothing
		{
			*context.memory.top_pos = 0 ;
		}
	}
	
	virtual ~MultipleInterpolationFromBottomToken2D() {} ;
	
	virtual std::string print() const
	{
		return std::string("multiple interpolate from bottom")  ;
	}
} ;

/** \brief Return the normalised distance between n geometries, starting from outside. Use the last three values on the stack.*/
class MultipleInterpolationFromTopToken3D : public Token
{
	std::vector<Geometry *> geos ;
public:
	MultipleInterpolationFromTopToken3D() : Token(false, std::make_pair(std::make_pair(TOKEN_MULTIPLE_INTERPOLATE_FROM_TOP_3D, 0),(double)(0)))  { };
	
	virtual void eval(Context & context) const
	{
		Point test(*context.memory.top_pos, *context.memory.prev_top_pos) ;
		
		context.memory.pop_back() ;
		context.memory.pop_back() ;
		Geometry * lastIn = NULL ;
		Geometry * previousIn = NULL ;
		Geometry * lastOut = NULL ;
		Geometry * previousOut = NULL ;
		for(std::vector<Geometry *>::const_iterator i = geos.begin() ; i != geos.end() ; ++i)
		{
			if((*i)->in(test))
			{
				if(lastIn)
					previousIn = lastIn ;
				lastIn = *i ;
			}
			else
			{
				if(lastOut)
					previousOut = lastOut ;
				lastOut = *i ;
			}
		}
		
		if(lastOut == NULL) // I am in all geometries.
		{
			Point p0(test) ;
			lastIn->project(&p0) ;
			*context.memory.top_pos = sqrt(squareDist2D(test, p0)) ;
			
		}
		else if (lastIn != NULL) // I am in a geometry, but not in the next
		{
			Point p0(test) ;
			lastIn->project(&p0) ; double din = sqrt(squareDist2D(test, p0)) ;
			Point p1(test) ;
			lastOut->project(&p1) ; double dout = sqrt(squareDist2D(test, p1)) ;
			*context.memory.top_pos = interpolate(din, dout) ;
			
		}
		else if(previousOut != NULL) // I am between 2 geometries
		{
			Point p0(test) ;
			previousOut->project(&p0) ; double din = sqrt(squareDist2D(test, p0)) ;
			Point p1(test) ;
			lastOut->project(&p1) ; double dout = sqrt(squareDist2D(test, p1)) ;
			*context.memory.top_pos = interpolate(din, dout) ;
		}
		else if(!geos.empty()) // I am outside the only geometry
		{
			Point p0(test) ;
			lastOut->project(&p0) ;
			*context.memory.top_pos = sqrt(squareDist2D(test, p0)) ;
		}
		else // there is nothing
		{
			*context.memory.top_pos = 0 ;
		}
	}
	
	virtual ~MultipleInterpolationFromTopToken3D() {} ;
	
	virtual std::string print() const
	{
		return std::string("multiple interpolate from top")  ;
	}
} ;

/** \brief Return the normalised distance between n geometries, starting from inside. Use the last three values on the stack.*/
class MultipleInterpolationFromBottomToken3D : public Token
{
	std::vector<Geometry *> geos ;
public:
	MultipleInterpolationFromBottomToken3D() : Token(false, std::make_pair(std::make_pair(TOKEN_MULTIPLE_INTERPOLATE_FROM_BOTTOM_3D, 0),(double)(0)))  { };
	
	virtual void eval(Context & context) const
	{
		Point test(*context.memory.top_pos, *context.memory.prev_top_pos, *(context.memory.top_pos-2)) ;
		context.memory.pop_back() ;
		context.memory.pop_back() ;
		
		Geometry * lastIn = NULL ;
		Geometry * previousIn = NULL ;
		Geometry * lastOut = NULL ;
		Geometry * previousOut = NULL ;
		for(std::vector<Geometry *>::const_reverse_iterator i = geos.rbegin() ; i != geos.rend() ; ++i)
		{
			if((*i)->in(test))
			{
				if(lastIn)
					previousIn = lastIn ;
				lastIn = *i ;
			}
			else
			{
				if(lastOut)
					previousOut = lastOut ;
				lastOut = *i ;
			}
		}
		
		if(lastOut == NULL) // I am in all geometries.
		{
			Point p0(test) ;
			lastIn->project(&p0) ;
			*context.memory.top_pos = sqrt(squareDist2D(test, p0)) ;
			
		}
		else if (lastIn != NULL) // I am in a geometry, but not in the next
		{
			Point p0(test) ;
			lastIn->project(&p0) ; double din = sqrt(squareDist2D(test, p0)) ;
			Point p1(test) ;
			lastOut->project(&p1) ; double dout = sqrt(squareDist2D(test, p1)) ;
			*context.memory.top_pos = interpolate(din, dout) ;
			
		}
		else if(previousOut != NULL) // I am between 2 geometries
		{
			Point p0(test) ;
			previousOut->project(&p0) ; double din = sqrt(squareDist2D(test, p0)) ;
			Point p1(test) ;
			lastOut->project(&p1) ; double dout = sqrt(squareDist2D(test, p1)) ;
			*context.memory.top_pos = interpolate(din, dout) ;
		}
		else if(!geos.empty()) // I am outside the only geometry
		{
			Point p0(test) ;
			lastOut->project(&p0) ;
			*context.memory.top_pos = sqrt(squareDist2D(test, p0)) ;
		}
		else // there is nothing
		{
			*context.memory.top_pos = 0 ;
		}
	}
	
	virtual ~MultipleInterpolationFromBottomToken3D() {} ;
	
	virtual std::string print() const
	{
		return std::string("multiple interpolate from bottom")  ;
	}
} ;

/** \brief Add the last two values on the stack*/
class PlusOperatorToken : public Token
{
public:
	PlusOperatorToken( bool nul = false) ;
	
	virtual void eval(Context & context) const
	{
	    *context.memory.prev_top_pos += *context.memory.top_pos ;
		 context.memory.pop_back() ;
	}
	
	virtual ~PlusOperatorToken() ;
	
	virtual std::string print() const
	{
		return std::string("+")  ;
	}
} ;

class EqualsToken : public Token
{
public:
	EqualsToken(): Token(false, std::make_pair(std::make_pair(TOKEN_EQUALS,0), 0)) { } ;
	
	virtual void eval(Context & context) const
	{
		
	    *context.memory.variables[context.memory.current_variable] = *context.memory.top_pos ;
		 context.memory.pop_back() ;
		 context.memory.pop_back() ;
	}
	
	virtual ~EqualsToken() { };
	
	virtual std::string print() const
	{
		return std::string("=")  ;
	}
} ;

/** \brief Substract the last two values on the stack*/
class MinusOperatorToken : public Token
{
public:
	MinusOperatorToken( bool nul = false) ;
	
	virtual void eval(Context & context) const
	{
		*context.memory.prev_top_pos -= *context.memory.top_pos ; 
		context.memory.pop_back() ;
	}
	
	virtual ~MinusOperatorToken() ;
	
	virtual std::string print() const
	{
		return std::string("-")  ;
	}
} ;

/** \brief Multiply the last two values on the stack*/
class TimesOperatorToken : public Token
{
public:
	TimesOperatorToken(bool nul = false);
	
	virtual void eval(Context & context) const
	{
		*context.memory.prev_top_pos *= *context.memory.top_pos ;
		context.memory.pop_back() ;
	}
	
	virtual ~TimesOperatorToken() ;
	
	virtual std::string print() const
	{
		return std::string("*")  ;
	}
} ;

/** \brief Divide the last two values on the stack*/
class DivideOperatorToken : public Token
{
public:
	DivideOperatorToken(bool nul = false) ;
	
	virtual void eval(Context & context) const
	{
		*context.memory.prev_top_pos = 
			*context.memory.prev_top_pos 
			/ *context.memory.top_pos; ;
		context.memory.pop_back() ;
	}
	
	virtual ~DivideOperatorToken() ;
	
	virtual std::string print() const
	{
		return std::string("/")  ;
	}
} ;

/** \brief take the last value of the stack power of the second-to-last value of the stack.*/
class PowerOperatorToken : public Token
{
public:
	PowerOperatorToken(bool nul = false) ;
	
	virtual void eval(Context & context) const
	{
		double val = *context.memory.prev_top_pos ;
		size_t pow = static_cast<size_t>(*context.memory.top_pos) -1;
	    context.memory.pop_back() ;
	    *context.memory.top_pos = val ;
		for(size_t i = 0 ; i < pow ; ++i)
		{
		    *context.memory.top_pos *= val;
		}
	}
	
	virtual ~PowerOperatorToken();
	
	virtual std::string print() const
	{
		return std::string("^")  ;
	}
} ;

/** \brief put the x argument to the nth power on the stack*/
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

/** \brief put the y argument to the nth power on the stack*/
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

/** \brief put the z argument to the nth power on the stack*/
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

/** \brief put the t argument to the nth power on the stack*/
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

/** \brief put the nth power of a value read on the heap on the stack*/
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

/** \brief Add a double to a value read on the heap*/
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

/** \brief Add two values read on the heap*/
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

/** \brief Multiply two values read on the heap*/
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

/** \brief Add a value read on the heap to the x argument*/
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

/** \brief Add a value read on the heap to the y argument*/
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

/** \brief Add a value read on the heap to the z argument*/
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

/** \brief Add a value read on the heap to the t argument*/
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
	
/** \brief Multiply a value read on the heap with a constant*/
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

/** \brief Multiply a value read on the heap with the x argument*/
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

/** \brief Multiply a value read on the heap with the y argument*/
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

/** \brief Multiply a value read on the heap with the z argument*/
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

/** \brief Multiply a value read on the heap with the t argument*/
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

/** \brief Substract to the x argument with a constant*/
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

/** \brief Substract to the y argument with a constant*/
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

/** \brief Substract to the z argument with a constant*/
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

/** \brief Substract to the t argument with a constant*/
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

/** \brief Substract to a constant with the x argument*/
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

/** \brief Substract to a constant with the y argument*/
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

/** \brief Substract to a constant with the t argument*/
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

/** \brief Substract to a constant with a value read on the heap*/
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

/** \brief Substract to a value read on the heap with a constant*/
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

/** \brief Substract a constant with a value read on the heap*/
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

/** \brief Substract a value read on the heap with a constant*/
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

/** \brief Substract a value read on the heap with the z argument*/
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

/** \brief Divide the x argument with a constant*/
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

/** \brief Divide the y argument with a constant*/
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

/** \brief Divide the z argument with a constant*/
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

/** \brief Divide the t argument with a constant*/
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

/** \brief Divide a constant with the x argument*/
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

/** \brief Divide a constant with the y argument*/
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

/** \brief Divide a constant with the z argument*/
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

/** \brief Divide a constant with the t argument*/
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

/** \brief take the nth power of the x argument and multiply with a constant*/
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

/** \brief take the nth power of the y argument and multiply with a constant*/
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

/** \brief take the nth power of the z argument and multiply with a constant*/
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

/** \brief take the nth power of the t argument and multiply with a constant*/
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


/** \brief multiply the x and y arguments*/
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

/** \brief multiply the x and z arguments*/
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

/** \brief multiply the y and z arguments*/
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
