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




namespace Amie
{
    
double sign(const double t) ;
double positivity(const double t) ;
double negativity(const double t) ;
double interpolate(const double a, const double b) ;

typedef double (*unaryFunctionPointer)(const double) ;
typedef double (*binaryFunctionPointer)(const double, double) ;
typedef double (*trinaryFunctionPointer)(const double, double, double) ;

// typedef std::vector<double> Memory ;
const size_t HEAP_SIZE = 32768 ;

/** \brief Memory structure for the VirtualMachine. It provides a stack and a heap.
*/
struct Memory
{
	double stack[HEAP_SIZE];
	double heap[HEAP_SIZE] ;

	/** \brief Constructor. Initialises stack and heap to 0.
	 * 
	 */
	Memory();
	
	~Memory() ;

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
	Context(const double & x_, const double & y_ , const double & z_ , const double & t_ , const double & u_ , const double & v_ , const double & w_ ) ;
	/** \brief Constructor, initialises the argument values. arguments not set are 0
	 * 
	 * @param x_ 
	 * @param y_ 
	 * @param z_ 
	 * @param t_ 
	 * @param u_ 
	 * @param v_ 
	 */
	Context(const double & x_, const double & y_, const double & z_, const double & t_, const double & u_, const double & v_ ) ;

	/** \brief Constructor, initialises the argument values. arguments not set are 0
	 * 
	 * @param x_ 
	 * @param y_ 
	 * @param z_ 
	 * @param t_ 
	 * @param u_ 
	 */
	Context(const double & x_, const double & y_, const double & z_, const double & t_, const double & u_) ;

	/** \brief Constructor, initialises the argument values. arguments not set are 0
	 * 
	 * @param x_ 
	 * @param y_ 
	 * @param z_ 
	 * @param t_ 
	 */
	Context(const double & x_, const double & y_, const double & z_, const double & t_) ;

	/** \brief Constructor, initialises the argument values. arguments not set are 0
	 * 
	 * @param x_ 
	 * @param y_ 
	 * @param z_ 
	 */
	Context(const double & x_, const double & y_, const double & z_) ;

	/** \brief Constructor, initialises the argument values. arguments not set are 0
	 * 
	 * @param x_ 
	 * @param y_ 
	 */
	Context(const double & x_, const double & y_) ;

	/** \brief Constructor, initialises the argument values. arguments not set are 0
	 * 
	 * @param x_ 
	 */
	Context(const double & x_) ;

	/** \brief Constructor, initialises the arguments to 0
	 * 
	 */
	Context() ;

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
	void set(const double & x_, const double & y_ , const double & z_ , const double & t_ , const double & u_ , const double & v_ , const double & w_ );

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
	void set(const double & x_, const double & y_, const double & z_ );

	/** \brief Sets the argument values. Arguments not set are kept.
	 * 
	 * @param x_ 
	 * @param y_ 
	 */
	void set(const double & x_, const double & y_);

	/** \brief Sets the argument values. Arguments not set are kept.
	 * 
	 * @param x_ 
	 */
	void set(const double & x_);
} ;


typedef enum : unsigned char
{
	TOKEN_OPERATION_CONSTANT = 0,
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
	TOKEN_OPERATION_MIN,
	TOKEN_OPERATION_INPLACE_MIN,
	TOKEN_OPERATION_MAX,
	TOKEN_OPERATION_INPLACE_MAX,
	TOKEN_OPERATION_ATAN2,
	TOKEN_OPERATION_INPLACE_ATAN2,
	TOKEN_OPERATION_INTERPOLATE,
	TOKEN_OPERATION_INPLACE_INTERPOLATE,
	TOKEN_OPERATION_GEO_OPERATION
} TokenOperationType ;


class GeometryOperation
{
public:
	GeometryOperation() ;
    virtual ~GeometryOperation() {} ;
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
	PositionOperation(const Segment & s_ ) ;
	PositionOperation(const std::vector<Segment> & s_ ) ;
	virtual void eval(double * a, double * b, double * c) const;
  virtual GeometryOperation * getCopy() const ;
	virtual int adressOffset() const;
} ;

/** \brief Put the distance between a point defined by the two last positions on the stack and its projection on a Line on the stack. */
class LineDistanceOperation : public GeometryOperation
{
	Line l ;
public:
	LineDistanceOperation(const Line & l_ ) ;
	
	virtual void eval(double * a, double * b, double * c) const;
	
	virtual GeometryOperation * getCopy() const ;
	
	virtual int adressOffset() const;
	
} ;

class InHomogeneousProjectionOperation : public GeometryOperation
{
	Geometry * inGeo ;
	std::vector<Segment> inProjector ;
	std::vector<Segment> outProjector ;
public:
	//by convention, the first point of each segment is assumed to be on the geometry
	InHomogeneousProjectionOperation(Geometry * inGeo, const std::vector<Segment> & inProjector, const std::vector<Segment> &outProjector) ;
	virtual void eval(double * a, double * b, double * c) const ;

	virtual GeometryOperation * getCopy() const ;
	
	virtual int adressOffset() const;
} ;


/** \brief Put 1 on the stack if Point taken from context is in the given Geometry, -1 otherwise */
class DomainOperation : public GeometryOperation
{
	const Geometry* geo ;
public:
	DomainOperation(const Geometry * g ) ;
	virtual void eval(double * a, double * b, double * c) const;
	virtual GeometryOperation * getCopy() const ;
	
	virtual int adressOffset() const;

} ;

/** \brief Put -1 point defined by the two last positions on the stack lies outside a given Geometry, 1 otherwise. */
class DomainBinaryOperation : public GeometryOperation
{
	const Geometry* geo ;
public:
	DomainBinaryOperation(const Geometry * g ) ;
	
	virtual void eval(double * a, double * b, double * c) const;
	
	virtual GeometryOperation * getCopy() const ;
	
	virtual int adressOffset() const;
} ;


/** \brief Put on the stack the distance between a point defined by the last positions on the stack and a stored position. */
class PointDistanceBinaryOperation : public GeometryOperation
{
	double x0 ;
	double y0 ;
public:
	PointDistanceBinaryOperation(const Point & p ) ;
	virtual void eval(double * a, double * b, double * c) const;
	
	virtual GeometryOperation * getCopy() const ;
	
	virtual int adressOffset() const;
} ;

/** \brief Put on the stack the distance between a point defined by the last positions on the stack and a stored position. */
class PointDistanceTrinaryOperation : public GeometryOperation
{
	double x0 ;
	double y0 ;
	double z0 ;
public:
	PointDistanceTrinaryOperation(const Point & p ) ;
	
	virtual void eval(double * a, double * b, double * c) const;
	
	virtual GeometryOperation * getCopy() const ;
	
	virtual int adressOffset() const;

} ;

/** \brief Put on the stack the two coordinates of a transformed point defined by the last positions on the stack given a rotation. */
class RotationBinaryOperation : public GeometryOperation
{
	double cangle ;
	double sangle ;
public:
	RotationBinaryOperation(double a ) ;
	
	virtual void eval(double * a, double * b, double * c) const;
	virtual GeometryOperation * getCopy() const ;
	virtual int adressOffset() const;
} ;

/** \brief Put on the stack the angle between a point defined by the last positions on the stack and a stored position. */
class AngleBinaryOperation : public GeometryOperation
{
	double cangle ;
	double sangle ;
	Point pivot ;
public:
	AngleBinaryOperation(double a, const Point & p ) ;
	virtual void eval(double * a, double * b, double * c) const;
	
	virtual GeometryOperation * getCopy() const ;
	
	virtual int adressOffset() const;

} ;

/** \brief Put on the stack squared distance between a point defined by the last two positions on the stack and a tored position */
class PointSquareDistanceBinaryOperation : public GeometryOperation
{
	Point base ;
public:
	PointSquareDistanceBinaryOperation(const Point & p ) ;
	virtual void eval(double * a, double * b, double * c) const;
	
	virtual GeometryOperation * getCopy() const ;
	
	virtual int adressOffset() const;
} ;

/** \brief Put 1 one the stack if a point defined by the two last positions on the stack is not visible from a stored position given an obstructing Geometry, 0 otherwise. */
class LineOfSightOperation : public GeometryOperation
{
	Point base ;
	const Geometry * obstruction ;
public:
	LineOfSightOperation(const Point & p,  const Geometry * o) ;
	
	virtual void eval(double * a, double * b, double * c) const;
	
	virtual GeometryOperation * getCopy() const ;
	
	virtual int adressOffset() const;
} ;

/** \brief put on the stack the distance between a point defined by the arguments and its projection on a stored Segment */
class ProjectionOperation3D : public GeometryOperation
{
	Segment s ;
public:
	ProjectionOperation3D(Segment s_ ) ;
	
	virtual void eval(double * a, double * b, double * c) const;
	virtual GeometryOperation * getCopy() const ;
	virtual int adressOffset() const;
} ;

/** \brief put on the stack the distance between a point defined by the arguments and its projection on a stored Segment */
class ProjectionOperation2D : public GeometryOperation
{
	Segment s ;
public:
	ProjectionOperation2D(Segment s_ ) ;
	
	virtual void eval(double * a, double * b, double * c) const;
	virtual GeometryOperation * getCopy() const ;
	
	virtual int adressOffset() const;
} ;

/** \brief put on the stack the distance between a point defined by the two last positions on the Stack and its projection on a stored Segment */
class ProjectionBinaryOperation : public GeometryOperation
{
	const Geometry * g ;
public:
	ProjectionBinaryOperation(const Geometry * s_ ) ;
	
	virtual void eval(double * a, double * b, double * c) const;
	
	virtual GeometryOperation * getCopy() const ;
	
	virtual int adressOffset() const;
} ;


class HatEnrichment : public GeometryOperation
{
    const Geometry * g ;
    Point  p ;
    Segment s ;
protected:
    double compute(double a, double b) const;
public:
    HatEnrichment(const Geometry * g , const Point & p, const Segment & s) ;
    
    virtual void eval(double * a, double * b, double * c) const;
    
    virtual GeometryOperation * getCopy() const ;
    
    virtual int adressOffset() const;
} ;

class HatEnrichmentAlt : public GeometryOperation
{
    const Geometry * g ;
    Point  head ;
    Point p0 ;
    Point p1 ;
public:
    HatEnrichmentAlt(const Geometry * g , const Point & head,  const Point & p0, const Point & p1) ;
    
    virtual void eval(double * a, double * b, double * c) const;
    
    virtual GeometryOperation * getCopy() const ;
    
    virtual int adressOffset() const;
} ;

class HatEnrichmentDerivative : public GeometryOperation
{
    const Geometry * g ;
    Point  p ;
    Segment s ;
    Variable v ;
public:
    HatEnrichmentDerivative(const Geometry * g , const Point & p, const Segment & s, Variable v) ;
    
    virtual void eval(double * a, double * b, double * c) const;
    
    virtual GeometryOperation * getCopy() const ;
    
    virtual int adressOffset() const;
} ;

class HatEnrichment3D : public GeometryOperation
{
    const Geometry * g ;
    Point  p ;
    TriPoint s ;
public:
    HatEnrichment3D(const Geometry * g , const Point & p, const TriPoint & s) ;
    
    virtual void eval(double * a, double * b, double * c) const;
    
    virtual GeometryOperation * getCopy() const ;
    
    virtual int adressOffset() const;
} ;

}



#endif

