//
// C++ Interface: vm_function_base
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef VM_FUNCTION_BASE_H
#define VM_FUNCTION_BASE_H

#include <string>
#include <vector>
#include <valarray>
#include <iostream>
#include <cmath>

#include "vm_refcount_token.h"
#include "../elements/integrable_entity.h"
#include "../utilities/matrixops.h"
#include "../geometry/geometry_base.h"

namespace Mu
{

struct GaussPointArray ;

typedef enum
{
	POSITION_TOKEN,
	PROJECTION_TOKEN
} PositionTokenType ;

/** \brief Function object used for runtime generation of functions, or symbolic manipulation of mathematical expressions.
 * Function can be defined using helper functions an a variety of constructors, or by writing them in RPN as a string:
 * \code 
Function f() ;          // null function
Function g("x 1 -") ;   // x -> x -1
Function s = f_sin(g) ; // x -> sin(x-1)
\endcode
 * Functions can be composed using all operators, including themselves :
\code
Function f("x") ;     // x -> x
Function m("2") ;     // x -> 2
f = f*m ;             // x -> 2x
Function g("x sin") ; // x -> sin(x) 
g = g(f) ;            // x-> sin(2x)
\endcode
 * Functions can have the variables x, y, z, t, u, v, w, which will be passed at evaluation.
 * Functions depend on VirtualMachine for their evaluation.
 */
class Function
{
	
	std::valarray<Function> derivative ; 
	std::vector< Point > iPoint ;
	
protected:
	const RefCountedToken toToken(const std::string str) const ;

	bool isOperator(const char c) const ;
	
	bool isNumeral(const char c) const ;
	
	bool isSeparator(const char c) const ;
	
	std::pair<size_t, RefCountedToken > getNext(size_t init, const char * form) ; 
	std::pair<size_t, RefCountedToken > getNext(size_t init, const std::string & form) ;
	
	Point * ptID ;
	int dofID ;
	
protected:
	ByteCode byteCode ;
	
	bool e_diff ;
	
	bool isBinaryOperator(const Token * t) const ;

	bool isUnaryOperator(const Token * t) const ;

	std::map<int, Vector *> precalc ;
	std::map<int, std::map<Variable, Vector *> > dprecalc ;
// 	void vectorizeOne(std::vector<RefCountedToken>   &bytecode,  size_t &lastAddress , int & precalculatedEnd ) const ;
// 	void vectorizeTwo(std::vector<RefCountedToken>   &bytecode,  size_t &lastAddress , int & precalculatedEnd ) const ;
// 	void factorize(std::vector<RefCountedToken>   &bytecode,  size_t &lastAddress , int & precalculatedEnd ) const ;
	
	bool compiled ;
public:
	
	/** \brief Default constructor. Creates a null function.
	 * 
	 */
	Function();

	/** \brief Copy-constructor. Function is replaced by source
	 * 
	 * @param  s Function source.
	 */
	Function(const Function & s) ;

	/** \brief Create function from RPN expression
	 * 
	 * @param f string to parse to instanciate function.
	 */
	Function(const char *f) ;

	/** \brief Create testing for the side of a line, given a space transform (x, y) -> x(x, y), y(x, y)
	 * 
	 * @param l Line defining the boundary
	 * @param x x transform
	 * @param y y transform
	 */
	Function(const Line & l, Function x, Function y) ;

	/** \brief Create testing for the side of a line, given a space transform (x, y) -> x(x, y), y(x, y) deduced from an ElementarySurface
	 * 
	 * @param l Line defining the boundary
	 * @param s ElementarySurface defining the transform
	 */
	Function(const Line & l, ElementarySurface * s) ;

	/** \brief Function returning the distance to a point, given a space transform (x, y) -> x(x, y), y(x, y).
	 * 
	 * @param p Point the distance to which to compute.
	 * @param x x transform
	 * @param y y transform
	 */
	Function(const Point & p, Function x, Function y) ;

	/** \brief Function returning the distance to a point, given a space transform (x, y, z) -> x(x, y, z), y(x, y, z), z(x, y, z).
	 * 
	 * @param p Point the distance to which to compute.
	 * @param x x transform
	 * @param y y transform
	 * @param z z transform
	 */
	Function(const Point & p, const Function &x, const Function & y,const  Function &z) ;

	/** \brief Function returning the distance to a point, given a space transform (x, y) -> x(x, y), y(x, y) deduced from an ElementarySurface
	 * 
	 * @param p Point the distance to which to compute.
	 * @param s ElementarySurface defining the transform
	 */
	Function(const Point & p, ElementarySurface * s) ;

	/** \brief Function computing the rotation of a point(x, y) round 0, 0, given a space transform (x, y) -> x(x, y), y(x, y) deduced from an ElementarySurface. 
   * This Function is not useful as itself as when computed, it will return a double. However it can be composed with other operations.
	 * @param a angle
	 * @param s ElementarySurface defining the transform
	 */
	Function(double a, ElementarySurface * s) ;

	/** \brief Function computing the rotation of a point(x, y) round 0, 0, given a space transform (x, y) -> x(x, y), y(x, y)
   * This Function is not useful as itself as when computed, it will return a double. However it can be composed with other operations.
	 * @param a angle
	 * @param x x transform
	 * @param y y transform
	 */
	Function(double a, Function x, Function y) ;


	/** \brief Function computing the rotation of a point(x, y) round a Point given a space transform (x, y) -> x(x, y), y(x, y)
 * This Function is not useful as itself as when computed, it will return a double. However it can be composed with other operations.
	 * @param a angle
	 * @param p pivot of the rotation
	 * @param x x transform
	 * @param y y transform
	 */
	Function(double a,const Point & p,  Function x, Function y) ;

	/** \brief Function computing the rotation of a point(x, y) round a Point given a space transform (x, y) -> x(x, y), y(x, y) deduced from an ElementarySurface. 
 * This Function is not useful as itself as when computed, it will return a double. However it can be composed with other operations.
	 * @param a angle
	 * @param p pivot of the rotation
	 * @param s ElementarySurface defining the transform
	 */
	Function(double a,const Point & p,  ElementarySurface * s) ;

	/** \brief Create function from RPN expression
	 * 
	 * @param f string to parse to instanciate function.
	 */
	Function(const std::string &f) ;

	/** \brief Polynomial in x, y, z, t given a 4D tensor of coefficents
	 * 
	 * @param coeffs coefficients of the polynomial
	 * @param diff compute the differentials
	 */
	Function(const std::valarray< std::valarray<Matrix> > & coeffs, bool diff = true) ;

	/** \brief Polynomial in x, y, z given a 3D tensor of coefficents
	 * 
	 * @param coeffs coefficients of the polynomial
	 * @param diff compute the differentials
	 */
	Function(const std::valarray<Matrix> & coeffs, bool diff = true) ;

	/** \brief Polynomial in x, y given a Matrix of coefficents
	 * 
	 * @param coeffs coefficients of the polynomial
	 * @param diff compute the differentials
	 */
	Function(const Matrix & coeffs, bool diff = true) ;

	/** \brief Polynomial in x given a array of coefficents
	 * 
	 * @param coeffs coefficients of the polynomial
	 * @param diff compute the differentials
	 */
	Function(const std::valarray<double> & coeffs, bool diff = true) ;

	/** \brief Function returning 0 or 1 depending on the position of a point with repect to a segment, given a space transform (x, y) -> x(x, y), y(x, y), or the distance to the projection to the segment
	 * 
	 * @param s Segment defining a boundary
	 * @param x x transform
	 * @param y y transform
	 * @param t type of function, Position or Projection distance
	 */
	Function(const Segment s, const Function & x, const Function & y, PositionTokenType t = POSITION_TOKEN) ;

	/** \brief Function returning 0 or 1 depending on the position of a point with repect to a segment, given a space transform (x, y) -> x(x, y), y(x, y) deduced from an ElementarySurface.  or the distance to the projection to the segment
	 * 
	 * @param s_ Segment defining a boundary
	 * @param s ElementarySurface defining the transform
	 * @param t type of function, Position or Projection distance
	 */
	Function(const Segment s_,  ElementarySurface * s, PositionTokenType = POSITION_TOKEN) ;

	/** \brief Function returning 0 or 1 depending on the position of a point with repect to a set of segments, given a space transform (x, y) -> x(x, y), y(x, y), or the distance to the projection to the segments
	 * 
	 * @param s Segments defining a boundary
	 * @param x x transform
	 * @param y y transform
	 * @param t type of function, Position or Projection distance
	 */
	Function(const std::vector<Segment> s, const Function & x, const Function & y,PositionTokenType = POSITION_TOKEN) ;

	/** \brief Function returning 0 or 1 depending on the position of a point with repect to a set of segments, given a space transform (x, y) -> x(x, y), y(x, y) deduced from an ElementarySurface.  or the distance to the projection to the segment
	 * 
	 * @param s_ Segments defining a boundary
	 * @param s ElementarySurface defining the transform
	 * @param t type of function, Position or Projection distance
	 */
	Function(const std::vector<Segment> s_, ElementarySurface * s,PositionTokenType = POSITION_TOKEN) ;


	/** \brief Create a function from the concatenation of two ByteCode s.
	 * No differential is created by this process
	 * @param b_0 first ByteCode
	 * @param b_1 second ByteCode
	 */
	Function(const ByteCode & b_0, const ByteCode & b_1) ;

	/** \brief Check whether a point p is hidden by a geometry g, given a space transform (x, y) -> x(x, y), y(x, y)
	 * 
	 * @param p Point whose visibility is checked
	 * @param g Hiding Geometry
	 * @param x x transform
	 * @param y y transform
	 */
	Function(const Point & p , const Geometry * g, const Function & x, const Function & y) ;

	/** \brief Function returning 1 if x,y is in the Geometry g and -1 otherwise
	 * 
	 * @param g domain defining Geometry
	 */
	Function(const Geometry * g) ;

	/** \brief Function returning 1 if x,y is in the Geometry g and -1 otherwise, given a space transform (x, y) -> x(x, y), y(x, y)
	 * 
	 * @param g domain defining Geometry
	 * @param x x transform
	 * @param y y transform
	 */
	Function(const Geometry * g, Function x, Function y) ;

	/** \brief Function returning 1 if x,y is in the Geometry g and -1 otherwise, given a space transform (x, y) -> x(x, y), y(x, y)  deduced from an ElementarySurface
	 * 
	 * @param g domain defining Geometry
	 * @param s ElementarySurface defining the transform
	 */
	Function(const Geometry * g, const ElementarySurface * s) ;

public:

	/** \brief Construct a function from two bytecodes and a token. 
	 * This is typically useful to create a function from an operation between two functions
	 * 
	 * @param b_0 first ByteCode
	 * @param b_1 second ByteCode
	 * @param op operator linking the two functions
	 * @param diff compute the differential 
	 */
	Function(const ByteCode &b_0, const ByteCode &b_1, RefCountedToken op, const bool diff = false); 

	/** \brief Construct a function from a bytecodes and a double  and a token. 
	 * This is typically useful to create a function from an operation between a function and a double
	 * 
	 * @param b_0 first ByteCode
	 * @param a scalar
	 * @param op operator linking the two functions
	 * @param diff compute the differential 
	 */
	Function(const ByteCode &b_0, const double a, RefCountedToken op, const bool diff = false); 


	/** \brief Construct a function from a bytecodes and a double  and a token. 
	 * This is typically useful to create a function from an operation between a function and a double
	 * @param a scalar
	 * @param b_0 first ByteCode
	 * @param op operator linking the two functions
	 * @param diff compute the differential 
	 */
	Function(const double a, const ByteCode &b_0,  RefCountedToken op, const bool diff= false) ;
	
	template<class ETYPE, class EABSTRACTTYPE>
	Function(Mesh<ETYPE, EABSTRACTTYPE>  * mesh, Variable v):  derivative(2)  , byteCode(1), e_diff(true), compiled(false)
	{

			this->dofID =-1 ;
			this->ptID = NULL ;
			if(v == XI)
				byteCode[0] = RefCountedToken(new MeshXDisplacementToken<ETYPE, EABSTRACTTYPE>(mesh)) ;
			if(v == ETA)
				byteCode[0] = RefCountedToken(new MeshYDisplacementToken<ETYPE, EABSTRACTTYPE>(mesh)) ;
			if(v == ZETA)
				byteCode[0] = RefCountedToken(new MeshZDisplacementToken<ETYPE, EABSTRACTTYPE>(mesh)) ;
			
			derivative[XI] = Function() ;
			derivative[ETA] = Function() ;
	
	}
	
	virtual ~Function()  ;
	
	/** \brief Check wether the function is nil
	 * @return true if function is identically 0
	 */
	bool isNull() const ;
	
	/** \brief return a reference to this function's ByteCode
	 * 
	 * @return a reference to this function's ByteCode
	 */
	const ByteCode & getByteCode() const ;

	/** \brief return a reference to this function's ByteCode
	 * 
	 * @return a reference to this function's ByteCode
	 */
	ByteCode & getByteCode() ;

	/** \brief return the ith token of this functions ByteCode
	 * 
	 * @return a reference to the ith Token of this function's ByteCode
	 */
	const RefCountedToken& getToken(const size_t i) const ;
	
	/** \brief return the size of the ByteCode
	 * 
	 * @return the size of the ByteCode
	 */
	size_t size() const { return byteCode.size() ; }
	
	/** \brief return the symbolic differential of this function if it has been computed, or a nil function otherwise
	 * 
	 * @param v Variable with which to differentiate
	 * @return the symbolic differential of this function if it has been computed, or a nil function otherwise
	 */
	const Function & d(const Variable v) const ;

	/** \brief return the symbolic differential of this function if it has been computed, or a nil function otherwise
	 * 
	 * @param v Variable with which to differentiate
	 * @return the symbolic differential of this function if it has been computed, or a nil function otherwise
	 */
	Function &d(const Variable v) ;

	/** \brief return the symbolic differentials of this function if they have been computed, or a nil function otherwise
	 * 
	 * @param v Variable with which to differentiate
	 * @return the symbolic differential of this function if it has been computed, or a nil function otherwise
	 */
	std::valarray<Function> getDerivatives() const;

	/** \brief return the symbolic differentials of this function if they have been computed, or a nil function otherwise
	 * 
	 * @param v Variable with which to differentiate
	 * @return the symbolic differential of this function if it has been computed, or a nil function otherwise
	 */
	std::valarray<Function> & getDerivatives();


	/** \brief return true if the differentiable has been computed
	 * 
	 * @return true if the differentiable has been computed
	 */
	bool isDifferentiable() const ;
	
	/** \brief Return the set of points to use for sub-tesselation for integration
	 * 
	 * @return the set of points to use for sub-tesselation for integration
	 */
	const std::vector< Point > & getIntegrationHint() const ;

	/** \brief Return the ith point to use for sub-tesselation for integration
	 * 
	 * @param  i index of the Point to return
	 * @return the ith point to use for sub-tesselation for integration
	 */
	Point  getIntegrationHint(size_t i) const ;


	/** \brief Set the integration hints: the set of points to use for sub-tesselation for integration
	 * 
	 * @param  pts set of points to use for sub-tesselation for integration
	 */
	void setIntegrationHint(const std::vector< Point > pts) ;

	/** \brief Add an the integration hint: a points to use for sub-tesselation for integration
	 * 
	 * @param  p points to use for sub-tesselation for integration
	 */
	void addIntegrationHint(const Point p) ;

	/** \brief Return true if the vector of integration hints is not empty
	 * 
	 * @return true if the vector of integration hints is not empty
	 */
	bool hasIntegrationHint() const ;
	
	/** \brief If this function is attached to a DOF, return its ID, -1 otherwise
	 * 
	 * @return If this function is attached to a DOF, return its ID, -1 otherwise
	 */
	int getDofID() const;

	/** \brief Set the ID of the DOF to which this function is attached
	 */
	void setDofID(size_t) ;
	
	/** \brief Return the point to which the function is attached
	 * 
	 * @return the point to which the function is attached
	 */
	Point * getPoint() const;

	/** \brief Attach a Point to the Function
	 * 
	 * @param  p Point to attach
	 */
	void setPoint(Point * p) ;
	
	/** \brief Copy-constructor. This operator conserves the id of the lvalue function. 
	 * Therefore, it can have strange side effects for example when copying std::vectors of functions, such as not 
	 * overwriting the ids as would be expected.
	 * 
	 * @param f source
	 * @return self
	 */
	Function & operator=(const Function &f) ;
	
	/** \brief Multiply with a Function
	 * 
	 * @param f Function to multiply with
	 * @return a new function
	 */
	Function operator*(const Function &f) const ;
	
	/** \brief Modify the function such that the result is multiplied by -1 if the point defined by the arguments is not in the Geometry
	 * 
	 * @param f Geometry defining the domain
	 * @return a new Function
	 */
	Function operator*(const Geometry *f) const ;
	
	/** \brief Divide with a Function
	 * 
	 * @param f Function to divide with
	 * @return a new Function
	 */
	Function operator/(const Function &f) const ;
	
	/** \brief Add a function
	 * 
	 * @param f Function to add
	 * @return a new Function
	 */
	Function operator+(const Function &f) const ;
	
	/** \brief Substract a Function
	 * 
	 * @param f Function to substract
	 * @return a new Function
	 */
	Function operator-(const Function &f) const ;
	
	/** \brief Multiply by a double
	 * 
	 * @param a factor
	 * @return a new Function
	 */
	Function operator*(const double a) const ;
	
	/** \brief divide by a double
	 * 
	 * @param a divisor
	 * @return a new Function
	 */
	Function operator/(const double a) const ;
	
	/** \brief add a double
	 * 
	 * @param a double to add
	 * @return a new Function
	 */
	Function operator+(const double a) const ;
	
	/** \brief substract a double
	 * 
	 * @param a double to substract
	 * @return a new Function
	 */
	Function operator-(const double a) const ;
	
	/** \brief take the integer power of a function
	 * 
	 * @param a power
	 * @return a new Function
	 */
	Function operator^(const int a) const ;
	
	/** \brief Multiply with a Function and assign 
	 * 
	 * @param f factor
	 */
	void operator*=(const Function &f)  ;
	
	/** \brief Multiply with a Function returning -1 if the point defined by the arguments is outside the Geometry.
	 * 
	 * @param f Geometry defining the domain
	 */
	void operator*=(const Geometry *f)  ;
	
	/** \brief Divide by a Function and assign
	 * 
	 * @param f factor
	 */
	void operator/=(const Function &f)  ;
	
	/** \brief Add a Function and assign
	 * 
	 * @param f Function to add
	 */
	void operator+=(const Function &f)  ;
	
	/** \brief Substract by a Function and assign
	 * 
	 * @param f Function to substract
	 */
	void operator-=(const Function &f)  ;
	
	/** \brief Multiply by a double and assign
	 * 
	 * @param a factor
	 */
	void operator*=(const double a)  ;
	
	/** \brief divide by a double and assign
	 * 
	 * @param a divisor
	 */
	void operator/=(const double a)  ;
	
	/** \brief Add a double and assign
	 * 
	 * @param a double to add
	 */
	void operator+=(const double a)  ;
	
	/** \brief Substract a double and assign
	 * 
	 * @param a double to substract
	 */
	void operator-=(const double a)  ;
	
	/** \brief Take the composite of two functions
	 * 
	 * @param f Function with which to compose
	 * @return a new Function
	 */
	Function operator()(const Function & f) const ;
	
	/** \brief Modify the Context (stack and arguments) by evaluating the ith Token of the ByteCode
	 * 
	 * @param i index of the Token
	 * @param context context to modify
	 */
	void tokenEval(size_t i, Context & context) const
	{
		byteCode[i].eval(context ) ;
	}
	
	/** \brief Compile the function: operations are vectorised (a+x)*b becomes a single token to evaluate, for example, sub-expressions which are multiply occuring are calculated only once and their result re-read, etc. 
	 * 
	 */
	void compile() ;


	/** \brief Precalculate the function derivatives at a set of locations
	 * 
	 * @param gp GaussPointArray defining the locations
	 * @param var Variables with respect to which differentiation should be calculated
	 * @param eps delata for the numerical evaluation of derivatives
	 */
	void preCalculate(const GaussPointArray & gp , std::vector<Variable> & var, const double eps = default_derivation_delta) ;

	/** \brief Precalculate the function  at a set of locations
	 * 
	 * @param gp GaussPointArray defining the locations
	 */
	void preCalculate(const GaussPointArray & gp) ;	

	/** \brief Check whether the function has been precalculated at a set of positions
	 * 
	 * @param gp GaussPointArray defining the locations
	 * @return true if the function has been precalculated
	 */
	bool precalculated(const GaussPointArray & gp) const ;

	/** \brief Check whether the function derivative has been precalculated at a set of positions
	 * 
	 * @param gp GaussPointArray defining the locations
	 * @param v Variable with respect to which differentiation should be calculated
	 * @return 
	 */
	bool precalculated(const GaussPointArray & gp, Variable v) const ;

	/** \brief return the precalculated values at a set of points
	 * 
	 * @param gp GaussPointArray  defining the locations
	 * @return Vector of values at those points
	 */
	const Vector & getPrecalculatedValue(const GaussPointArray & gp) const ;

	/** \brief return the precalculated derivative values at a set of points
	 * 
	 * @param gp GaussPointArray  defining the locations
	 * @param v Variable with respect to which differentiation should be calculated
	 * @return Vector of values at those points
	 */
	const Vector & getPrecalculatedValue(const GaussPointArray & gp, Variable v) const ;
	
} ;

struct DtF ;
struct GtM ;
struct GtV ;
struct GtMtG ;
struct VGtM ;
struct VGtV ;
struct DtD ;
struct VGtMtVG ;
struct DtGtMtG ;
struct GDtM ;
struct GDtMtGD ;
struct GDtMtG ;
struct GtMtGD ;

/** \brief Structure used for the lazy computation of the differential of a function
 * This is equivalent to a derivative, but can be used on all functions, even if their 
 * derivatives have not been specified. In that case, a numerical derivative will be computed
*/
struct Differential
{
	const Function & f ;
	const Variable & v ;

	/** \brief Constructor, initalise the references
	 * 
	 * @param u Function
	 * @param m Variable
	 */
	Differential(const Function &u, const Variable & m) : f(u), v(m) { } ;

	/** \brief Create a structure for the lazy evaluation of the Differential * Function
	 * 
	 * @param f factor
	 * @return DtF
	 */
	DtF operator *(const Function & f) const ;


	/** \brief Create a structure for the lazy evaluation of the Differential * Differential
	 * 
	 * @param f factor
	 * @return DtD
	 */
	DtD operator *(const Differential & f) const ;


	/** \brief Create a structure for the lazy evaluation of a Gradient * Matrix * Gradient * Differential
	 * 
	 * @param g GtMtG
	 * @return DtGtMtG
	 */
	DtGtMtG operator *(const GtMtG & g) const ;
} ;

/** \brief Create a structure for the lazy evaluation of a Gradient 
 * The gradient of a function is evaluated as a Matrix, which dimensions are given by the 
 * list of specifiend space variables given at the time of evaluation
*/
struct Gradient
{
	const Function & f ;
	const bool transpose ;
	/** \brief Constructor, initalise the references
	 * 
	 * @param u Function
	 * @param t true if the result should be transposed
	 */
	Gradient(const Function &u, bool t = false) : f(u), transpose(t) { };

	/** \brief Create a structure for the lazy evaluation of a Gradient * Matrix
	 * 
	 * @param f Matrix
	 * @return GtM
	 */
	GtM operator *(const Matrix & f) const ;

	/** \brief Create a structure for the lazy evaluation of a Gradient * Vector
	 * 
	 * @param f Vector
	 * @return GtV
	 */
	GtV operator *(const Vector & f) const ;
} ;

/** \brief Create a structure for the lazy evaluation of a GradientDot 
 * The gradient of a function is evaluated as a Matrix, which dimensions are given by the 
 * list of specifiend space variables given at the time of evaluation. The result 
 * is further differentiated with respect to time
*/
struct GradientDot
{
	const Function & f ;
	const bool transpose ;

	/**  \brief Constructor, initalise the references
	 * 
	 * @param u Function
	 * @param t transpose result if true
	 */
	GradientDot(const Function &u, bool t = false) : f(u), transpose(t) { };

	/** \brief Create a structure for the lazy evaluation of a GradientDot * Matrix
	 * 
	 * @param f Matrix
	 * @return GDtM
	 */
	GDtM operator *(const Matrix & f) const ;
} ;

/** \brief Create a structure for the lazy evaluation of a VectorGradient 
 * The gradient of a function is evaluated as a single line Matrix, which dimension is given by the 
 * list of specifiend space variables given at the time of evaluation.
*/
struct VectorGradient
{
	const Function & f ;
	const bool transpose ;

	/**  \brief Constructor, initalise the references
	 * 
	 * @param u Function
	 * @param t transpose result if true
	 */
	VectorGradient(const Function &u, bool t = false) : f(u), transpose(t) { };

	/** \brief Create a structure for the lazy evaluation of a VectorGradient * Matrix
	 * 
	 * @param f Matrix
	 * @return VGtM
	 */
	VGtM operator *(const Matrix & f) const ;

	/** \brief Create a structure for the lazy evaluation of a VectorGradient * Vector
	 * 
	 * @param f Vector
	 * @return VGtV
	 */
	VGtV operator *(const Vector & f) const ;
} ;

/** \brief Structure for the lazy evaluation of a Gradient * Matrix structure*/
struct GtM
{
	const Gradient & first ;
	const Matrix & second ;
	
	/** \brief Constructor, initalise the references
	 * 
	 * @param g Gradient
	 * @param f Matrix
	 */
	GtM(const Gradient & g, const Matrix & f) : first(g), second(f) { };

	/** \brief Create a structure for the lazy evaluation of a Gradient * Matrix * Gradient
	 * 
	 * @param f Gradient
	 * @return GtMtG
	 */
	GtMtG operator*(const Mu::Gradient & f) const ;

	/** \brief Create a structure for the lazy evaluation of a Gradient * Matrix * GradientDot
	 * 
	 * @param f Gradient
	 * @return GtMtGD
	 */
	GtMtGD operator*(const Mu::GradientDot & f) const ;
} ;

/** \brief Structure for the lazy evaluation of a GradientDot * Matrix */
struct GDtM
{
	const GradientDot & first ;
	const Matrix & second ;

	/** \brief Constructor, initalise the references
	 * 
	 * @param g GradientDot
	 * @param f Matrix
	 */
	GDtM(const GradientDot & g, const Matrix & f) : first(g), second(f) { };

	/** \brief Create a structure for the lazy evaluation of a GradientDot * Matrix * GradientDot
	 * 
	 * @param f GradientDot
	 * @return GDtMtGD
	 */
	GDtMtGD operator*(const Mu::GradientDot & f) const ;

	/** \brief Create a structure for the lazy evaluation of a GradientDot * Matrix * Gradient
	 * 
	 * @param f Gradient
	 * @return GDtMtGD
	 */
	GDtMtG operator*(const Mu::Gradient & f) const ;
} ;

/** \brief Structure for the lazy evaluation of a GradientDot * Matrix * GradientDot */
struct GDtMtGD
{
	const GradientDot & first ;
	const Matrix & second ;
	const GradientDot & third ;
	/** \brief Constructor, initalise the references
	 * 
	 * @param g GDtM
	 * @param f GradientDot
	 */
	GDtMtGD(const GDtM & g, const GradientDot & p) : first(g.first), second(g.second), third(p) { };

} ;

/** \brief Structure for the lazy evaluation of a GradientDot * Matrix * Gradient */
struct GDtMtG
{
	const GradientDot & first ;
	const Matrix & second ;
	const Gradient & third ;
	/** \brief Constructor, initalise the references
	 * 
	 * @param g GDtM
	 * @param f GradientDot
	 */
	GDtMtG(const GDtM & g, const Gradient & p) : first(g.first), second(g.second), third(p) { };

} ;

/** \brief Structure for the lazy evaluation of a GradientDot * Matrix * Gradient */
struct GtMtGD
{
	const Gradient & first ;
	const Matrix & second ;
	const GradientDot & third ;
	/** \brief Constructor, initalise the references
	 * 
	 * @param g GtM
	 * @param f GradientDot
	 */
	GtMtGD(const GtM & g, const GradientDot & p) : first(g.first), second(g.second), third(p) { };

} ;

/** \brief Structure for the lazy evaluation of a Differential * Function */
struct DtF
{
	const Differential & d ;
	const Function & f ;

	/** \brief Constructor, initalise the references
	 * 
	 * @param d_ Differential
	 * @param f_ Function
	 */
	DtF(const Differential & d_, const Function & f_) : d(d_), f(f_) { } ;
} ;

/** \brief Structure for the lazy evaluation of a Differential * Differential */
struct DtD
{
	const Differential & d ;
	const Differential & f ;

	/** \brief Constructor, initalise the references
	 * 
	 * @param d_ Differential
	 * @param f_ Differential
	 */
	DtD(const Differential & d_, const Differential & f_) : d(d_), f(f_) { } ;
} ;

/** \brief Structure for the lazy evaluation of a VectorGradient * Matrix */
struct VGtM
{
	const VectorGradient & first ;
	const Matrix & second ;
	
	/** \brief Constructor, initalise the references
	 * 
	 * @param g VectorGradient
	 * @param f Matrix
	 */
	VGtM(const VectorGradient & g, const Matrix & f) : first(g), second(f) { };

	/** \brief Create a structure for the lazy evaluation of a VectorGradient * Matrix * VectorGradient
	 * 
	 * @param f VectorGradient
	 * @return VGtMtVG
	 */
	VGtMtVG operator*(const Mu::VectorGradient & f) const ;
} ;

/** \brief Structure for the lazy evaluation of a Gradient * Vector */
struct GtV
{
	const Gradient & first ;
	const Vector & second ;
	
	/**  \brief Constructor, initalise the references
	 * 
	 * @param g Gradient
	 * @param f Vector
	 */
	GtV(const Gradient & g, const Vector & f) : first(g), second(f) { };
} ;

/** \brief Structure for the lazy evaluation of a VectorGradient * Vector */
struct VGtV
{
	const VectorGradient & first ;
	const Vector & second ;
	
	/**  \brief Constructor, initalise the references
	 * 
	 * @param g VectorGradient
	 * @param f Vector
	 */
	VGtV(const VectorGradient & g, const Vector & f) : first(g), second(f) { };
} ;

/** \brief Structure for the lazy evaluation of a Gradient * Matrix * Gradient*/
struct GtMtG
{
	const Gradient & first ;
	const Matrix & second ;
	const Gradient & third ;
	
	/**  \brief Constructor, initalise the references
	 * 
	 * @param g Gradient
	 * @param f Matrix
	 * @param g_ Gradient
	 */
	GtMtG(const Gradient & g, const Matrix & f,const Gradient & g_) : first(g), second(f), third(g_) { };
	
} ;

/** \brief Structure for the lazy evaluation of a Differential * Gradient * Matrix * Gradient*/
struct DtGtMtG
{
	const Differential & first ;
	const GtMtG & second ;

	/** \brief Constructor, initalise the references
	 * 
	 * @param d Differential
	 * @param g GtMtG
	 */
	DtGtMtG(const Differential & d, const GtMtG & g) : first(d), second(g) {  };
} ;

/** \brief Structure for the lazy evaluation of a VectorGradient * Matrix * VectorGradient*/
struct VGtMtVG
{
	const VectorGradient & first ;
	const Matrix & second ;
	const VectorGradient & third ;
	
	/** \brief Constructor, initalise the references
	 * 
	 * @param g VectorGradient
	 * @param f Matrix
	 * @param g_ VectorGradient
	 */
	VGtMtVG(const VectorGradient & g, const Matrix & f,const VectorGradient & g_) : first(g), second(f), third(g_) { };
	
} ;

} ;

/** \brief Create a Function from the substraction of a Function to a double
 * 
 * @param a double
 * @param f Function
 * @return a new Function
 */
Mu::Function operator-(const double & a, const Mu::Function &f) ;

/** \brief Create a Function from the multiplication of a double to a function
 * 
 * @param a double
 * @param f Function
 * @return a new Function
 */
Mu::Function operator*(const double & a, const Mu::Function &f) ;

/** \brief Create a Function from the addition of a Function to a double
 * 
 * @param a double
 * @param f Function
 * @return a new Function
 */
Mu::Function operator+(const double & a, const Mu::Function &f) ;

/** \brief Create a Function from the division by a Function of a double
 * 
 * @param a double
 * @param f Function
 * @return a new Function
 */
Mu::Function operator/(const double & a, const Mu::Function &f) ;

/** \brief Helper function to create a Function which is the square root of the argument
 * 
 * @param f Function
 * @return a new Function
 */
Mu::Function f_sqrt(const Mu::Function &f) ;

/** \brief Helper function to create a Function which is the exponential of the argument
 * 
 * @param f Function
 * @return a new Function
 */
Mu::Function f_exp(const Mu::Function &f) ;

/** \brief Helper function to create a Function which is the absolute value of the argument
 * 
 * @param f Function
 * @return a new Function
 */
Mu::Function f_abs(const Mu::Function &f) ;

/**  \brief Helper function to create a Function which is the log of the argument
 * 
 * @param f Function
 * @return a new Function
 */
Mu::Function f_log(const Mu::Function &f) ;

/** \brief Helper function to create a Function which is arc tangeant of the arguments
 * 
 * @param f0 Function
 * @param f1 Function
 * @return a new Function
 */
Mu::Function f_atan2(const Mu::Function &f0, const Mu::Function &f1) ;

/** \brief Helper function to create a Function which is the sine of the argument
 * 
 * @param f Function
 * @return a new Function
 */
Mu::Function f_sin(const Mu::Function &f) ;

/** \brief Helper function to create a Function which is the cosine of the argument
 * 
 * @param f Function
 * @return a new Function
 */
Mu::Function f_cos(const Mu::Function &f) ;

/** \brief Helper function to create a Function which is the sign of the argument
 * 
 * @param f Function
 * @return a new Function
 */
Mu::Function f_sign(const Mu::Function &f) ;

/** \brief Helper function to create a Function which is the ith Cylindrical Bessel Function of the argument
 * 
 * @param i index of the function
 * @param f Function
 * @return a new Function
 */
Mu::Function f_cyl_bessel_j(int i, const Mu::Function &f) ;

/** \brief Helper function to create a Function which is 1 if the argument is positive, 0 otherwise
 * 
 * @param f Function
 * @return a new Function
 */
Mu::Function f_positivity(const Mu::Function &f) ;

/** \brief Helper function to create a Function which is 1 if the argument is negative, 0 otherwise
 * 
 * @param f Function
 * @return a new Function
 */
Mu::Function f_negativity(const Mu::Function &f) ;

/** \brief  Helper function to create a Function which is the interpolated value between two functions. 0 as an argument will yield the value of the first function and 1 the second, intermediate values will yield the weighted average at that point.
 * 
 * @param f0 Function
 * @param f1 Function
 * @return a new Function
 */
Mu::Function f_interpolate(const Mu::Function &f0, const Mu::Function &f1) ;

/** \brief Helper function to create a Function which computes the distance between a point, given by the arguments and its pojection on a given Geometry. The coordinates are transformed using the provided functions
 * 
 * @param g Geometry
 * @param x x transform Function
 * @param y y transform Function
 * @return a new Function
 */
Mu::Function f_project(const Mu::Geometry *g, const Mu::Function &x, const Mu::Function &y) ;

/** \brief Helper function to compute the curvilinear absciss on a SegmentedLine
 *
 * @param s SegmentedLine serving as source for the computation
 * @param fromHead compute starting from the Head
 * @param x x transform Function
 * @param y y transform Function
 * @return a new Function
*/
Mu::Function f_curvilinear_x(const Mu::SegmentedLine * s, bool fromHead,   const Mu::Function &x, const Mu::Function &y) ;

/** \brief Helper function to compute the curvilinear ordinate on a SegmentedLine
 *
 * @param s SegmentedLine serving as source for the computation
 * @param fromHead compute starting from the Head
 * @param x x transform Function
 * @param y y transform Function
 * @return a new Function
*/
Mu::Function f_curvilinear_y(const Mu::SegmentedLine * s, bool fromHead,   const Mu::Function &x, const Mu::Function &y) ;
#endif

