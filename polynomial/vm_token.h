//
// C++ Interface: vm_token
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
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
const size_t HEAP_SIZE = 8192 ;
const size_t FUNCTION_LENGTH = HEAP_SIZE ;

/** \brief Memory structure for the VirtualMachine. It provides a stack and a heap.
*/
struct Memory
{
	double stack[HEAP_SIZE];
	double heap[HEAP_SIZE] ;

	/** \brief Constructor. Initialises stack and heap to 0.
	 * 
	 */
	Memory()
	{
	} ;
	
	~Memory() 
	{ 
	} ;

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

	} ;

	/** \brief Constructor, initialises the argument values. arguments not set are 0
	 * 
	 * @param x_ 
	 * @param y_ 
	 * @param z_ 
	 */
	Context(const double & x_, const double & y_, const double & z_) : x(x_), y(y_), z(z_), t(0), u(0), v(0), w(0) 
	{ 

	} ;

	/** \brief Constructor, initialises the argument values. arguments not set are 0
	 * 
	 * @param x_ 
	 * @param y_ 
	 */
	Context(const double & x_, const double & y_) : x(x_), y(y_), z(0), t(0), u(0), v(0), w(0) 
	{ 

	} ;

	/** \brief Constructor, initialises the argument values. arguments not set are 0
	 * 
	 * @param x_ 
	 */
	Context(const double & x_) : x(x_), y(0), z(0), t(0), u(0), v(0), w(0) 
	{ 

	} ;

	/** \brief Constructor, initialises the arguments to 0
	 * 
	 */
	Context() : x(0), y(0), z(0), t(0), u(0), v(0), w(0) 
	{ 

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


typedef enum 
{
	TOKEN_OPERATION_CONSTANT,
	TOKEN_OPERATION_X,
	TOKEN_OPERATION_Y,
	TOKEN_OPERATION_Z,
	TOKEN_OPERATION_T,
	TOKEN_OPERATION_U,
	TOKEN_OPERATION_V,
	TOKEN_OPERATION_W,
	TOKEN_OPERATION_PLUS,
	TOKEN_OPERATION_INPLACE_PLUS,
	TOKEN_OPERATION_MINUS,
	TOKEN_OPERATION_INPLACE_MINUS,
	TOKEN_OPERATION_TIMES,
	TOKEN_OPERATION_INPLACE_TIMES,
	TOKEN_OPERATION_DIVIDES,
	TOKEN_OPERATION_INPLACE_DIVIDES,
	TOKEN_OPERATION_POWER,
	TOKEN_OPERATION_INPLACE_POWER,
	TOKEN_OPERATION_ABS,
	TOKEN_OPERATION_INPLACE_ABS,
	TOKEN_OPERATION_COS,
	TOKEN_OPERATION_INPLACE_COS,
	TOKEN_OPERATION_SIN,
	TOKEN_OPERATION_INPLACE_SIN,
	TOKEN_OPERATION_TAN,
	TOKEN_OPERATION_INPLACE_TAN,
	TOKEN_OPERATION_COSH,
	TOKEN_OPERATION_INPLACE_COSH,
	TOKEN_OPERATION_SINH,
	TOKEN_OPERATION_INPLACE_SINH,
	TOKEN_OPERATION_TANH,
	TOKEN_OPERATION_INPLACE_TANH,
	TOKEN_OPERATION_EXP,
	TOKEN_OPERATION_INPLACE_EXP,
	TOKEN_OPERATION_SIGN,
	TOKEN_OPERATION_INPLACE_SIGN,
	TOKEN_OPERATION_POSITIVITY,
	TOKEN_OPERATION_INPLACE_POSITIVITY,
	TOKEN_OPERATION_NEGATIVITY,
	TOKEN_OPERATION_INPLACE_NEGATIVITY,
	TOKEN_OPERATION_LOG,
	TOKEN_OPERATION_INPLACE_LOG,
	TOKEN_OPERATION_SQRT,
	TOKEN_OPERATION_INPLACE_SQRT,
	TOKEN_OPERATION_BESSEL,
	TOKEN_OPERATION_INPLACE_BESSEL,
	TOKEN_OPERATION_ATAN2,
	TOKEN_OPERATION_INPLACE_ATAN2,
	TOKEN_OPERATION_INTERPOLATE,
	TOKEN_OPERATION_INPLACE_INTERPOLATE,
	TOKEN_OPERATION_GEO_OPERATION
} TokenOperationType ;


class GeometryOperation
{
public:
	GeometryOperation() { } ;
	virtual void eval(double * a, double * b, double * c) const = 0;
	virtual GeometryOperation * getCopy() const = 0;
	virtual int adressOffset() const = 0 ;
} ;

/** \brief Put -1 or 1 on the stack depending on whether we are within or without a domain delimited by a set of Segment s */
class PositionOperation : public GeometryOperation
{
	std::vector<Segment> s ;
	Point w ;
public:
	PositionOperation(Segment s_ ) 
	{
		Point vector(-s_.vector().y, s_.vector().x) ;
		w=( s_.midPoint()+vector*10.) ;
		s.push_back(s_)  ;
	}
	
	PositionOperation(std::vector<Segment> s_ ) 
	{
		for(size_t i = 0 ; i < s_.size() ; i++)
		{

			s.push_back(s_[i])  ;
		}			
		Point vector(-s_[0].vector().y, s_[0].vector().x) ;
		w=( s_[0].midPoint()+vector*10.) ;
	}
	virtual void eval(double * a, double * b, double * c) const
	{
		Point test(*a,*b) ;
		size_t intersections = 0 ;
		for(size_t i = 0 ; i < s.size() ; i++)
		{
			if(s[i].intersects(test, w))
			{
				intersections++ ;
			}
		}
		
		*c =  (intersections & 1) * 2 - 1;
// 		if(intersections%2 == 1)
// 		    *context.memory.top_pos = -1 ;
// 		else
// 		    *context.memory.top_pos = 1 ;
	}
  virtual GeometryOperation * getCopy() const 
  {
		return new PositionOperation(s) ;
	};
	
	virtual int adressOffset() const
	{
		return -1 ;
	}
} ;

/** \brief Put the distance between a point defined by the two last positions on the stack and its projection on a Line on the stack. */
class LineDistanceOperation : public GeometryOperation
{
	Line l ;
public:
		LineDistanceOperation(const Line & l_ ) : l(l_)
	{
	}
	
	virtual void eval(double * a, double * b, double * c) const
	{
		
		Point test(*a, *b) ;
		*c=  sqrt(squareDist2D(test,l.projection(test)));
	}
	
	virtual GeometryOperation * getCopy() const 
  {
		return new LineDistanceOperation(l) ;
	};
	
	virtual int adressOffset() const
	{
		return -1 ;
	}
	
} ;

class InHomogeneousProjectionOperation : public GeometryOperation
{
	Geometry * inGeo ;
	std::vector<Segment> inProjector ;
	std::vector<Segment> outProjector ;
public:
	//by convention, the first point of each segment is assumed to be on the geometry
	InHomogeneousProjectionOperation(Geometry * inGeo, const std::vector<Segment> & inProjector, const std::vector<Segment> &outProjector) : inGeo(inGeo), inProjector(inProjector), outProjector(outProjector)
	{
	}

	virtual void eval(double * a, double * b, double * c) const ;

	virtual GeometryOperation * getCopy() const 
  {
		return new InHomogeneousProjectionOperation(inGeo, inProjector, outProjector) ;
	};
	
	virtual int adressOffset() const
	{
		return -1 ;
	}
} ;


/** \brief Put 1 on the stack if Point taken from context is in the given Geometry, -1 otherwise */
class DomainOperation : public GeometryOperation
{
	const Geometry* geo ;
public:
	DomainOperation(const Geometry * g ) 
	{
		geo = g ;
	}
	
	virtual void eval(double * a, double * b, double * c) const
	{
		Point p(*a,*b,*c) ;
		if(geo->in(p))
		    *c = 1 ;
		else
		    *c = -1 ;
	}
	
	virtual GeometryOperation * getCopy() const 
  {
		return new DomainOperation(geo) ;
	};
	
	virtual int adressOffset() const
	{
		return -2 ;
	}

} ;

/** \brief Put -1 point defined by the two last positions on the stack lies outside a given Geometry, 1 otherwise. */
class DomainBinaryOperation : public GeometryOperation
{
	const Geometry* geo ;
public:
	DomainBinaryOperation(const Geometry * g ) : geo(g)
	{
	}
	
	virtual void eval(double * a, double * b, double * c) const
	{
		
		Point p(*a, *b) ;
		if(geo->in(p))
			*c = 1 ;
		else
			*c = -1 ;
	}
	
	virtual GeometryOperation * getCopy() const 
  {
		return new DomainBinaryOperation(geo) ;
	};
	
	virtual int adressOffset() const
	{
		return -1 ;
	}

} ;


/** \brief Put on the stack the distance between a point defined by the last positions on the stack and a stored position. */
class PointDistanceBinaryOperation : public GeometryOperation
{
	double x0 ;
	double y0 ;
public:
	PointDistanceBinaryOperation(const Point & p ) :  x0(p.x), y0(p.y)
	{
	}
	
	virtual void eval(double * a, double * b, double * c) const
	{
		double x = *a-x0 ;
		double y = *b-y0 ;
		*c = sqrt(x*x+y*y) ;

	}
	
	virtual GeometryOperation * getCopy() const 
  {
		return new PointDistanceBinaryOperation(Point(x0, y0)) ;
	};
	
	virtual int adressOffset() const
	{
		return -1 ;
	}
} ;

/** \brief Put on the stack the distance between a point defined by the last positions on the stack and a stored position. */
class PointDistanceTrinaryOperation : public GeometryOperation
{
	double x0 ;
	double y0 ;
	double z0 ;
public:
	PointDistanceTrinaryOperation(const Point & p ) : x0(p.x), y0(p.y), z0(p.z)
	{
	}
	
	virtual void eval(double * a, double * b, double * c) const
	{
		double x = *a-x0 ;
		double y = *b-y0 ;
		double z = *c-z0 ;
		*c = sqrt(x*x+y*y+z*z) ;

	}
	
	virtual GeometryOperation * getCopy() const 
  {
		return new PointDistanceTrinaryOperation(Point(x0, y0, z0)) ;
	};
	
	virtual int adressOffset() const
	{
		return -2 ;
	}

} ;

/** \brief Put on the stack the two coordinates of a transformed point defined by the last positions on the stack given a rotation. */
class RotationBinaryOperation : public GeometryOperation
{
	double cangle ;
	double sangle ;
public:
	RotationBinaryOperation(double a ) : cangle(cos(a)), sangle(sin(a))
	{
	}
	
	virtual void eval(double * a, double * b, double * c) const
	{
		
		double x = *a ;
		double y =  *b ;
		*a = x*cangle + y*sangle ;
		*b = -x*sangle + y*cangle ;

	}
	virtual GeometryOperation * getCopy() const 
	{
		return new RotationBinaryOperation(acos(cangle)) ;
	}
	
	virtual int adressOffset() const
	{
		return 0 ;
	}
	
} ;

/** \brief Put on the stack the angle between a point defined by the last positions on the stack and a stored position. */
class AngleBinaryOperation : public GeometryOperation
{
	double cangle ;
	double sangle ;
	Point pivot ;
public:
	AngleBinaryOperation(double a, const Point & p ) :cangle(cos(a)), sangle(sin(a)), pivot(p.x*cos(a)+p.y*sin(a), -p.x*sin(a)+p.y*cos(a))
	{
	}
	
	virtual void eval(double * a, double * b, double * c) const
	{
		
		double x = *a ;
		double y =  *b ;
		double x_t = x*cangle + y*sangle ;
		double y_t = -x*sangle + y*cangle ;
		*c = atan2(y_t-pivot.y, x_t-pivot.x) ;

	}
	
	virtual GeometryOperation * getCopy() const 
	{
		return new AngleBinaryOperation(acos(cangle), pivot) ;
	}
	
	virtual int adressOffset() const
	{
		return -1 ;
	}

} ;

/** \brief Put on the stack squared distance between a point defined by the last two positions on the stack and a tored position */
class PointSquareDistanceBinaryOperation : public GeometryOperation
{
	Point base ;
public:
	PointSquareDistanceBinaryOperation(const Point & p ) 
	{
		base = p ;
	}
	
	virtual void eval(double * a, double * b, double * c) const
	{
		
		Point p(*a, *b) ;

		*c = squareDist2D(p, base) ;

	}
	
	virtual GeometryOperation * getCopy() const 
	{
		return new PointSquareDistanceBinaryOperation(base) ;
	}
	
	virtual int adressOffset() const
	{
		return -1 ;
	}
	
} ;

/** \brief Put 1 one the stack if a point defined by the two last positions on the stack is not visible from a stored position given an obstructing Geometry, 0 otherwise. */
class LineOfSightOperation : public GeometryOperation
{
	Point base ;
	const Geometry * obstruction ;
public:
	LineOfSightOperation(const Point & p,  const Geometry * o) :  base(p), obstruction(o)
	{
	}
	
	virtual void eval(double * a, double * b, double * c) const
	{
		
		Point p(*a, *b) ;
		

		*c = Segment(base, p).intersects(obstruction) ;

	}
	
	virtual GeometryOperation * getCopy() const 
	{
		return new LineOfSightOperation(base, obstruction) ;
	}
	
	virtual int adressOffset() const
	{
		return -1 ;
	}

} ;

/** \brief put on the stack the distance between a point defined by the arguments and its projection on a stored Segment */
class ProjectionOperation3D : public GeometryOperation
{
	Segment s ;
public:
	ProjectionOperation3D(Segment s_ ) : s(s_)
	{
	}
	
	virtual void eval(double * a, double * b, double * c) const
	{
		
	   *c = dist(Point(*a, *b, *c), s.project(Point(*a, *b, *c))) ;
	}
	
	virtual GeometryOperation * getCopy() const 
	{
		return new ProjectionOperation3D(s) ;
	}
	
	virtual int adressOffset() const
	{
		return -2 ;
	}
} ;

/** \brief put on the stack the distance between a point defined by the arguments and its projection on a stored Segment */
class ProjectionOperation2D : public GeometryOperation
{
	Segment s ;
public:
	ProjectionOperation2D(Segment s_ ) : s(s_)
	{
	}
	
	virtual void eval(double * a, double * b, double * c) const
	{
		
	    *c = dist(Point(*a, *b), s.project(Point(*a, *b))) ;
	}
	
	virtual GeometryOperation * getCopy() const 
	{
		return new ProjectionOperation2D(s) ;
	}
	
	virtual int adressOffset() const
	{
		return -1 ;
	}
} ;

/** \brief put on the stack the distance between a point defined by the two last positions on the Stack and its projection on a stored Segment */
class ProjectionBinaryOperation : public GeometryOperation
{
	const Geometry * g ;
public:
	ProjectionBinaryOperation(const Geometry * s_ ) :  g(s_)
	{
	}
	
	virtual void eval(double * a, double * b, double * c) const
	{
		Point p(*a, *b) ;
		Point p_(p) ;
		g->project(&p_) ;

		*c = sqrt(squareDist2D(p, p_)) ;

	}
	
	virtual GeometryOperation * getCopy() const 
	{
		return new ProjectionBinaryOperation(g) ;
	}
	
	virtual int adressOffset() const
	{
		return -1 ;
	}
} ;














} ;









#endif
