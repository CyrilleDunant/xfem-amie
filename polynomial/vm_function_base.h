//
// C++ Interface: vm_function_base
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
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
// #include <boost/tuple/tuple.hpp>

#include "../elements/integrable_entity.h"
#include "../utilities/matrixops.h"
#include "vm_token.h"
#include "../geometry/geometry_base.h"

namespace Amie
{

struct GaussPointArray ;
struct ElementarySurface ;
struct ElementaryVolume ;

const size_t HEAP_VARIABLE_TRANSFORM_OFFSET = 4096 ;

struct functionParseElement {
    TokenOperationType first ;
    double second ;
    std::string third ;
} ;

typedef enum {
    POSITION_TOKEN,
    PROJECTION_TOKEN
} PositionTokenType ;

typedef enum : unsigned char {
    NO_TEMPORARY,
    SET_TEMPORARY,
    SET_GET_TEMPORARY_A,
    SET_GET_TEMPORARY_B,
    GET_TEMPORARY_A,
    GET_TEMPORARY_B,
} TemporayUsageType ;

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
class Function final
{
protected :
    std::vector< Point > iPoint ;
    std::valarray<Function *> * derivative ;
    std::vector<Function *> * transforms ;
    bool derivativeTransformed = false ;

    /*void defaultInitialise()
    {
    	derivative(nullptr),
    	e_diff(false) ;
    	byteCodeSize(0) ;
    	constNumber(0) ;
    	byteCode.resize(TOKEN_OPERATION_CONSTANT,FUNCTION_LENGTH, );
    	geo_op.resize((GeometryOperation *)nullptr,FUNCTION_LENGTH);
    	use_temp.resize(NO_TEMPORARY,FUNCTION_LENGTH);
    	values.resize(0.,FUNCTION_LENGTH) ;
    	adress_a.resize(0,FUNCTION_LENGTH*4) ;
    	dofID =-1 ;
    	ptID = nullptr ;
    	hasGeoOp = false ;
    }*/

protected:
    functionParseElement toToken ( const std::string & str, int iter, std::vector<double> & v ) const ;

    bool isOperator ( const char c ) const ;

    bool isNumeral ( const char c ) const ;

    bool isSeparator ( const char c ) const ;

    std::pair<size_t, functionParseElement> getNext ( size_t init, const char * form, int iter, std::vector<double> & val ) ;
    std::pair<size_t, functionParseElement> getNext ( size_t init, const std::string & form, int iter, std::vector<double> & val ) ;

    Point * ptID;
    int dofID;

protected:


    bool e_diff = false;
    bool isBinaryOperator ( TokenOperationType t ) const ;
    bool isUnaryOperator ( TokenOperationType t ) const ;
    bool isBinaryVectorisableOperator ( TokenOperationType t ) const ;
    bool isTrinaryVectorisableOperator ( TokenOperationType t ) const ;
    bool isTrinaryOperator ( TokenOperationType t ) const ;
    std::map<int, Vector *> precalc ;
    std::map<int, std::map<Variable, Vector *> > dprecalc ;

public:
    std::vector<Variable> transformed ;
    bool hasGeoOp ;
    void initialiseAdresses ( size_t offset = 0 ) ;
    std::vector<TokenOperationType> byteCode ;
    std::valarray<GeometryOperation *> geo_op ;
    std::vector<double> values ;
    std::vector<unsigned short int> adress_a ;
    std::vector<unsigned short int> adress_t ;
// 	Function * xtransform ;
// 	Function * ytransform ;
// 	Function * ztransform ;
// 	Function * ttransform ;
    /** \brief Default constructor. Creates a null function.
     *
     */
    Function();

    /** \brief Copy-constructor. Function is replaced by source
     *
     * @param  s Function source.
     */
    Function ( const Function & s ) ;

    /** \brief Create function from RPN expression
     *
     * @param f string to parse to instanciate function.
     */
    Function ( const char *f ) ;

    /** \brief Create function from RPN expression
     *
     * @param f string to parse to instanciate function.
     * @param val double value for "constant" token
     */
    Function ( const char *f, double val ) ;

    /** \brief Create function from RPN expression
     *
     * @param f string to parse to instanciate function.
     * @param val double values for "constant" tokens
     */
    Function ( const char *f, std::vector<double> & val ) ;

    /** \brief Create testing for the side of a line, given a space transform (x, y) -> x(x, y), y(x, y) deduced from an ElementarySurface
     *
     * @param l Line defining the boundary
     * @param s ElementarySurface defining the transform
     */
    Function ( const Line & l, ElementarySurface * s ) ;

    /** \brief Function returning the distance to a point, given a space transform (x, y) -> x(x, y), y(x, y) deduced from an ElementarySurface
     *
     * @param p Point the distance to which to compute.
     * @param s ElementarySurface defining the transform
     */
    Function ( const Point & p, ElementarySurface * s ) ;

    /** \brief Function returning the distance to a point, given a space transform (x, y) -> x(x, y), y(x, y) deduced from an ElementarySurface
     *
     * @param p Point the distance to which to compute.
     * @param s ElementarySurface defining the transform
     */
    Function ( const Point & p, ElementaryVolume * s ) ;

    /** \brief Function computing the rotation of a point(x, y) round 0, 0, given a space transform (x, y) -> x(x, y), y(x, y) deduced from an ElementarySurface.
    * This Function is not useful as itself as when computed, it will return a double. However it can be composed with other operations.
     * @param a angle
     * @param s ElementarySurface defining the transform
     */
    Function ( double a, ElementarySurface * s ) ;


    /** \brief Function computing the rotation of a point(x, y) round a Point given a space transform (x, y) -> x(x, y), y(x, y) deduced from an ElementarySurface.
    * This Function is not useful as itself as when computed, it will return a double. However it can be composed with other operations.
     * @param a angle
     * @param p pivot of the rotation
     * @param s ElementarySurface defining the transform
     */
    Function ( double a,const Point & p,  ElementarySurface * s ) ;

    /** \brief Create function from RPN expression
     *
     * @param f string to parse to instanciate function.
     */
    Function ( const std::string &f ) ;

    Function ( const std::string &f, double v ) ;

    Function ( const std::string &f, std::vector<double> & v ) ;

    /** \brief Function returning 0 or 1 depending on the position of a point with repect to a segment, given a space transform (x, y) -> x(x, y), y(x, y) deduced from an ElementarySurface.  or the distance to the projection to the segment
     *
     * @param s_ Segment defining a boundary
     * @param s ElementarySurface defining the transform
     * @param t type of function, Position or Projection distance
     */
    Function ( const Segment s_,  ElementarySurface * s, PositionTokenType = POSITION_TOKEN ) ;


    /** \brief Function returning 0 or 1 depending on the position of a point with repect to a set of segments, given a space transform (x, y) -> x(x, y), y(x, y) deduced from an ElementarySurface.  or the distance to the projection to the segment
     *
     * @param s_ Segments defining a boundary
     * @param s ElementarySurface defining the transform
     * @param t type of function, Position or Projection distance
     */
    Function ( const std::vector<Segment> s_, ElementarySurface * s,PositionTokenType = POSITION_TOKEN ) ;

    /** \brief Function returning 1 if x,y is in the Geometry g and -1 otherwise, given a space transform (x, y) -> x(x, y), y(x, y)  deduced from an ElementarySurface
     *
     * @param g domain defining Geometry
     * @param s ElementarySurface defining the transform
     */
    Function ( const Geometry * g, const ElementarySurface * s ) ;

public:


    /** \brief Construct a function doing barycentric interpolation from a set of segments and a boundary
     */
    Function ( Geometry * inGeo, const std::vector<Segment> & inProjector, const std::vector<Segment> &outProjector, const Function &x, const Function &y ) ;

    virtual ~Function() final ;

// 	void setTransform(const ElementarySurface * s) ;
// 	void setTransform(const ElementaryVolume * s) ;
// 	void setTransform(const Function & x, const Function & y) ;
// 	void setTransform(const Function & x, const Function & y, const Function & z) ;

    /** \brief return the size of the ByteCode
     *
     * @return the size of the ByteCode
     */
    size_t size() const {
        return byteCode.size() ;
    }

    /** \brief return the symbolic differential of this function if it has been computed, or a nil function otherwise
     *
     * @param v Variable with which to differentiate
     * @return the symbolic differential of this function if it has been computed, or a nil function otherwise
     */
    const Function & d ( const Variable v ) const ;

    /** \brief return the symbolic differential of this function if it has been computed, or a nil function otherwise
     *
     * @param v Variable with which to differentiate
     * @return the symbolic differential of this function if it has been computed, or a nil function otherwise
     */
    Function  & d ( const Variable v ) ;

    /** \brief return the symbolic differentials of this function if they have been computed, or a nil function otherwise
     *
     * @param v Variable with which to differentiate
     * @return the symbolic differential of this function if it has been computed, or a nil function otherwise
     */
    std::valarray<Function *>  & getDerivatives() const;

    /** \brief return the symbolic differentials of this function if they have been computed, or a nil function otherwise
     *
     * @param v Variable with which to differentiate
     * @return the symbolic differential of this function if it has been computed, or a nil function otherwise
     */
    std::valarray<Function *> & getDerivatives();

    void setNumberOfDerivatives ( int n ) ;
    int getNumberOfDerivatives() const ;
    void setDerivative ( const Variable v, Function & f ) ;

    size_t derivationDepth() const ;

    /** \brief return true if the differentiable has been computed
     *
     * @return true if the differentiable has been computed
     */
    bool isDifferentiable() const ;

    bool isDifferentiable ( const Variable v ) const ;
    bool isDifferentiable ( size_t i ) const ;

    bool hasVariableTransform() const {
        return transforms ;
    }

    void setVariableTransform ( const Variable v, Function & f, bool replacetoken = true ) ;

    void makeVariableTransformDerivative() ;
    bool hasVariableTransformedDerivatives() const {
        return derivativeTransformed ;
    }
    Function & transform ( size_t i ) ;
    const Function & transform ( size_t i ) const ;

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
    const Point & getIntegrationHint ( size_t i ) const ;


    /** \brief Set the integration hints: the set of points to use for sub-tesselation for integration
     *
     * @param  pts set of points to use for sub-tesselation for integration
     */
    void setIntegrationHint ( const std::vector< Point > pts ) ;

    /** \brief Add an the integration hint: a points to use for sub-tesselation for integration
     *
     * @param  p points to use for sub-tesselation for integration
     */
    void addIntegrationHint ( const Point p ) ;

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
    void setDofID ( size_t ) ;

    /** \brief Return the point to which the function is attached
     *
     * @return the point to which the function is attached
     */
    Point * getPoint() const;

    /** \brief Attach a Point to the Function
     *
     * @param  p Point to attach
     */
    void setPoint ( Point * p ) ;

    /** \brief Copy-constructor. This operator conserves the id of the lvalue function.
     * Therefore, it can have strange side effects for example when copying std::vectors of functions, such as not
     * overwriting the ids as would be expected.
     *
     * @param f source
     * @return self
     */
    Function & operator= ( const Function &f ) ;

    /** \brief Multiply with a Function
     *
     * @param f Function to multiply with
     * @return a new function
     */
    Function operator* ( const Function &f ) const ;


    /** \brief Divide with a Function
     *
     * @param f Function to divide with
     * @return a new Function
     */
    Function operator/ ( const Function &f ) const ;

    /** \brief Add a function
     *
     * @param f Function to add
     * @return a new Function
     */
    Function operator+ ( const Function &f ) const ;

    /** \brief Substract a Function
     *
     * @param f Function to substract
     * @return a new Function
     */
    Function operator- ( const Function &f ) const ;

    /** \brief Multiply by a double
     *
     * @param a factor
     * @return a new Function
     */
    Function operator* ( const double a ) const ;

    /** \brief divide by a double
     *
     * @param a divisor
     * @return a new Function
     */
    Function operator/ ( const double a ) const ;

    /** \brief add a double
     *
     * @param a double to add
     * @return a new Function
     */
    Function operator+ ( const double a ) const ;

    /** \brief substract a double
     *
     * @param a double to substract
     * @return a new Function
     */
    Function operator- ( const double a ) const ;

    Function operator-() const  ;
    /** \brief take the integer power of a function
     *
     * @param a power
     * @return a new Function
     */
    Function operator^ ( const int a ) const ;

    /** \brief Multiply with a Function and assign
     *
     * @param f factor
     */
    void operator*= ( const Function &f )  ;

    /** \brief Divide by a Function and assign
     *
     * @param f factor
     */
    void operator/= ( const Function &f )  ;

    /** \brief Add a Function and assign
     *
     * @param f Function to add
     */
    void operator+= ( const Function &f )  ;

    /** \brief Substract by a Function and assign
     *
     * @param f Function to substract
     */
    void operator-= ( const Function &f )  ;

    /** \brief Multiply by a double and assign
     *
     * @param a factor
     */
    void operator*= ( const double a )  ;

    /** \brief divide by a double and assign
     *
     * @param a divisor
     */
    void operator/= ( const double a )  ;

    /** \brief Add a double and assign
     *
     * @param a double to add
     */
    void operator+= ( const double a )  ;

    /** \brief Substract a double and assign
     *
     * @param a double to substract
     */
    void operator-= ( const double a )  ;

    /** \brief Create a Function from the substraction of a Function to a double
    *
    * @param a double
    * @param f Function
    * @return a new Function
    */
    friend Function operator- ( double a, Amie::Function f ) {
        Amie::Function ret ( f ) ;
        ret.values.clear();
        ret.values.push_back ( a ) ;
        ret.values.insert ( ret.values.end(), f.values.begin(), f.values.end() ) ;
        if ( f.hasGeoOp ) {
            ret.geo_op.resize ( HEAP_SIZE, ( GeometryOperation* ) nullptr );
        }

        for ( size_t i = 0 ; i < f.byteCode.size() ; i++ ) {
            if ( f.hasGeoOp ) {
                if ( f.geo_op[i] ) {
                    ret.geo_op[i] = f.geo_op[i]->getCopy() ;
                }
            }

            if ( f.adress_a[i*4] >= HEAP_SIZE-f.values.size()-2 ) {
                ret.adress_a[ ( i ) *4] = f.adress_a[i*4]-1 ;
            }

            if ( f.adress_a[i*4+1] >= HEAP_SIZE-f.values.size()-2 ) {
                ret.adress_a[ ( i ) *4+1] = f.adress_a[i*4+1]-1 ;
            }

            if ( f.adress_a[i*4+2] >= HEAP_SIZE-f.values.size()-2 ) {
                ret.adress_a[ ( i ) *4+2] = f.adress_a[i*4+2]-1 ;
            }

            if ( f.adress_a[i*4] >= 8 && f.adress_a[i*4] < HEAP_SIZE-1-f.values.size() ) {
                ret.adress_a[ ( i ) *4] = f.adress_a[i*4] ;
            }

            if ( f.adress_a[i*4+1] >= 8 && f.adress_a[i*4+1] < HEAP_SIZE-1-f.values.size() ) {
                ret.adress_a[ ( i ) *4+1] = f.adress_a[i*4+1] ;
            }

            if ( f.adress_a[i*4+2] >= 8 && f.adress_a[i*4+2] < HEAP_SIZE-1-f.values.size() ) {
                ret.adress_a[ ( i ) *4+2] = f.adress_a[i*4+2] ;
            }

            if ( f.adress_a[i*4] < 8 ) {
                ret.adress_a[ ( i ) *4] = f.adress_a[i*4] ;
            }

            if ( f.adress_a[i*4+1] < 8 ) {
                ret.adress_a[ ( i ) *4+1] = f.adress_a[i*4+1] ;
            }

            if ( f.adress_a[i*4+2] < 8 ) {
                ret.adress_a[ ( i ) *4+2] = f.adress_a[i*4+2] ;
            }

        }
        ret.byteCode.push_back ( TOKEN_OPERATION_MINUS );
        ret.adress_a.push_back ( 0 );
        ret.adress_a.push_back ( 0 );
        ret.adress_a.push_back ( 0 );
        ret.adress_a.push_back ( 0 );
        ret.adress_a[ ( ret.byteCode.size()-1 ) *4+2] = 8 ;
        ret.adress_a[ ( ret.byteCode.size()-1 ) *4+1] = 8 ;
        ret.adress_a[ ( ret.byteCode.size()-1 ) *4] = HEAP_SIZE-1 ;


        size_t n = f.getNumberOfDerivatives() ;
        if ( n > 0 ) {
            ret.setNumberOfDerivatives ( n ) ;
            for ( size_t i = 0 ; i < n ; i++ ) {
                if ( f.isDifferentiable ( i ) ) {
                    Amie::Function diff =  f.d ( ( const Variable ) i ) * -1;
                    ret.setDerivative ( ( const Variable ) i, diff ) ;
                }
            }
        }
        ret.hasGeoOp = f.hasGeoOp ;
        return ret ;

    }


    /** \brief Create a Function from the multiplication of a double to a function
    *
    * @param a double
    * @param f Function
    * @return a new Function
    */
    friend Function operator* ( const double & a, const Amie::Function &f ) {
        Amie::Function ret ( f ) ;
        ret.values.clear();
        ret.values.push_back ( a ) ;
        ret.values.insert ( ret.values.end(), f.values.begin(), f.values.end() ) ;
        if ( f.hasGeoOp ) {
            ret.geo_op.resize ( HEAP_SIZE, ( GeometryOperation* ) nullptr );
        }

        for ( size_t i = 0 ; i < f.byteCode.size() ; i++ ) {
            if ( f.hasGeoOp )
                if ( f.geo_op[i] ) {
                    ret.geo_op[i] = f.geo_op[i]->getCopy() ;
                }

            if ( f.adress_a[i*4] >= HEAP_SIZE-2-f.values.size() ) {
                ret.adress_a[ ( i ) *4] = f.adress_a[i*4]-1 ;
            }

            if ( f.adress_a[i*4+1] >= HEAP_SIZE-2-f.values.size() ) {
                ret.adress_a[ ( i ) *4+1] = f.adress_a[i*4+1]-1 ;
            }

            if ( f.adress_a[i*4+2] >= HEAP_SIZE-2-f.values.size() ) {
                ret.adress_a[ ( i ) *4+2] = f.adress_a[i*4+2]-1 ;
            }

            if ( f.adress_a[i*4] >= 8 && f.adress_a[i*4] < HEAP_SIZE-f.values.size()-1 ) {
                ret.adress_a[ ( i ) *4] = f.adress_a[i*4]+1 ;
            }
            if ( f.adress_a[i*4+1] >= 8 && f.adress_a[i*4+1] < HEAP_SIZE-f.values.size()-1 ) {
                ret.adress_a[ ( i ) *4+1] = f.adress_a[i*4+1]+1 ;
            }
            if ( f.adress_a[i*4+2] >= 8 && f.adress_a[i*4+2] < HEAP_SIZE-f.values.size()-1 ) {
                ret.adress_a[ ( i ) *4+2] = f.adress_a[i*4+2]+1 ;
            }
            if ( f.adress_a[i*4] < 8 ) {
                ret.adress_a[ ( i ) *4] = f.adress_a[i*4] ;
            }
            if ( f.adress_a[i*4+1] < 8 ) {
                ret.adress_a[ ( i ) *4+1] = f.adress_a[i*4+1] ;
            }
            if ( f.adress_a[i*4+2] < 8 ) {
                ret.adress_a[ ( i ) *4+2] = f.adress_a[i*4+2] ;
            }

        }
        ret.byteCode.push_back ( TOKEN_OPERATION_TIMES );
        ret.adress_a.push_back ( 0 );
        ret.adress_a.push_back ( 0 );
        ret.adress_a.push_back ( 0 );
        ret.adress_a.push_back ( 0 );
        ret.adress_a[ ( ret.byteCode.size()-1 ) *4+2] = 8 ;
        ret.adress_a[ ( ret.byteCode.size()-1 ) *4+1] = 8 ;
        ret.adress_a[ ( ret.byteCode.size()-1 ) *4] = HEAP_SIZE-1 ;


        size_t n = f.getNumberOfDerivatives() ;
        if ( n > 0 ) {
            ret.setNumberOfDerivatives ( n ) ;
            for ( size_t i = 0 ; i < n ; i++ ) {
                if ( f.isDifferentiable ( i ) ) {
                    Amie::Function diff =  f.d ( ( const Variable ) i ) * a;
                    ret.setDerivative ( ( const Variable ) i, diff ) ;
                }
            }
        }

        return ret ;
    }

    /** \brief Create a Function from the addition of a Function to a double
    *
    * @param a double
    * @param f Function
    * @return a new Function
    */
    friend Function operator+ ( const double & a, const Amie::Function &f ) {
        Amie::Function ret ( f ) ;

        ret.values.clear();
        ret.values.push_back ( a ) ;
        ret.values.insert ( ret.values.end(), f.values.begin(), f.values.end() ) ;
        if ( f.hasGeoOp ) {
            ret.geo_op.resize ( HEAP_SIZE, ( GeometryOperation* ) nullptr );
        }

        for ( size_t i = 0 ; i < f.byteCode.size() ; i++ ) {

            if ( f.hasGeoOp )
                if ( f.geo_op[i] ) {
                    ret.geo_op[i] = f.geo_op[i]->getCopy() ;
                }

            if ( f.adress_a[i*4] >= HEAP_SIZE-2-f.values.size() ) {
                ret.adress_a[ ( i ) *4] = f.adress_a[i*4]-1 ;
            }

            if ( f.adress_a[i*4+1] >= HEAP_SIZE-2-f.values.size() ) {
                ret.adress_a[ ( i ) *4+1] = f.adress_a[i*4+1]-1 ;
            }

            if ( f.adress_a[i*4+2] >= HEAP_SIZE-2-f.values.size() ) {
                ret.adress_a[ ( i ) *4+2] = f.adress_a[i*4+2]-1 ;
            }

            if ( f.adress_a[i*4] >= 8 && f.adress_a[i*4] < HEAP_SIZE-f.values.size()-1 ) {
                ret.adress_a[ ( i ) *4] = f.adress_a[i*4]+1 ;
            }
            if ( f.adress_a[i*4+1] >= 8 && f.adress_a[i*4+1] < HEAP_SIZE-f.values.size()-1 ) {
                ret.adress_a[ ( i ) *4+1] = f.adress_a[i*4+1]+1 ;
            }
            if ( f.adress_a[i*4+2] >= 8 && f.adress_a[i*4+2] < HEAP_SIZE-f.values.size()-1 ) {
                ret.adress_a[ ( i ) *4+2] = f.adress_a[i*4+2]+1 ;
            }
            if ( f.adress_a[i*4] < 8 ) {
                ret.adress_a[ ( i ) *4] = f.adress_a[i*4] ;
            }
            if ( f.adress_a[i*4+1] < 8 ) {
                ret.adress_a[ ( i ) *4+1] = f.adress_a[i*4+1] ;
            }
            if ( f.adress_a[i*4+2] < 8 ) {
                ret.adress_a[ ( i ) *4+2] = f.adress_a[i*4+2] ;
            }

        }
        ret.byteCode.push_back ( TOKEN_OPERATION_PLUS );
        ret.adress_a.push_back ( 0 );
        ret.adress_a.push_back ( 0 );
        ret.adress_a.push_back ( 0 );
        ret.adress_a.push_back ( 0 );
        ret.adress_a[ ( ret.byteCode.size()-1 ) *4+2] = 8 ;
        ret.adress_a[ ( ret.byteCode.size()-1 ) *4+1] = 8 ;
        ret.adress_a[ ( ret.byteCode.size()-1 ) *4] = HEAP_SIZE-1 ;


        size_t n = f.getNumberOfDerivatives() ;
        if ( n > 0 ) {
            ret.setNumberOfDerivatives ( n ) ;
            for ( size_t i = 0 ; i < n ; i++ ) {
                if ( f.isDifferentiable ( i ) ) {
                    Amie::Function diff =  f.d ( ( const Variable ) i );
                    ret.setDerivative ( ( const Variable ) i, diff ) ;
                }
            }
        }

        return ret ;
    }

    /** \brief Create a Function from the division by a Function of a double
    *
    * @param a double
    * @param f Function
    * @return a new Function
    */
    friend Function operator/ ( const double & a, const Amie::Function &f ) {
        Amie::Function ret ( f ) ;
        ret.values.clear();
        ret.values.push_back ( a ) ;
        ret.values.insert ( ret.values.end(), f.values.begin(), f.values.end() ) ;
        if ( f.hasGeoOp ) {
            ret.geo_op.resize ( HEAP_SIZE, ( GeometryOperation* ) nullptr );
        }

        for ( size_t i = 0 ; i < f.byteCode.size() ; i++ ) {

            if ( f.hasGeoOp )
                if ( f.geo_op[i] ) {
                    ret.geo_op[i] = f.geo_op[i]->getCopy() ;
                }

            if ( f.adress_a[i*4] >= HEAP_SIZE-2-f.values.size() ) {
                ret.adress_a[ ( i ) *4] = f.adress_a[i*4]-1 ;
            }

            if ( f.adress_a[i*4+1] >= HEAP_SIZE-2-f.values.size() ) {
                ret.adress_a[ ( i ) *4+1] = f.adress_a[i*4+1]-1 ;
            }

            if ( f.adress_a[i*4+2] >= HEAP_SIZE-2-f.values.size() ) {
                ret.adress_a[ ( i ) *4+2] = f.adress_a[i*4+2]-1 ;
            }

            if ( f.adress_a[i*4] >= 8 && f.adress_a[i*4] < HEAP_SIZE-f.values.size()-1 ) {
                ret.adress_a[ ( i ) *4] = f.adress_a[i*4]+1 ;
            }
            if ( f.adress_a[i*4+1] >= 8 && f.adress_a[i*4+1] < HEAP_SIZE-f.values.size()-1 ) {
                ret.adress_a[ ( i ) *4+1] = f.adress_a[i*4+1]+1 ;
            }
            if ( f.adress_a[i*4+2] >= 8 && f.adress_a[i*4+2] < HEAP_SIZE-f.values.size()-1 ) {
                ret.adress_a[ ( i ) *4+2] = f.adress_a[i*4+2]+1 ;
            }
            if ( f.adress_a[i*4] < 8 ) {
                ret.adress_a[ ( i ) *4] = f.adress_a[i*4] ;
            }
            if ( f.adress_a[i*4+1] < 8 ) {
                ret.adress_a[ ( i ) *4+1] = f.adress_a[i*4+1] ;
            }
            if ( f.adress_a[i*4+2] < 8 ) {
                ret.adress_a[ ( i ) *4+2] = f.adress_a[i*4+2] ;
            }

        }
        ret.byteCode.push_back ( TOKEN_OPERATION_DIVIDES );
        ret.adress_a.push_back ( 0 );
        ret.adress_a.push_back ( 0 );
        ret.adress_a.push_back ( 0 );
        ret.adress_a.push_back ( 0 );
        ret.adress_a[ ( ret.byteCode.size()-1 ) *4+2] = 8 ;
        ret.adress_a[ ( ret.byteCode.size()-1 ) *4+1] = 8 ;
        ret.adress_a[ ( ret.byteCode.size()-1 ) *4] = HEAP_SIZE-1 ;


        size_t n = f.getNumberOfDerivatives() ;
        if ( n > 0 ) {
            ret.setNumberOfDerivatives ( n ) ;
            for ( size_t i = 0 ; i < n ; i++ ) {
                if ( f.isDifferentiable ( i ) ) {
                    Amie::Function diff = ( -a ) * f.d ( ( const Variable ) i ) / ( f*f ) ;
                    ret.setDerivative ( ( const Variable ) i, diff ) ;
                }
            }
        }

        return ret ;

    }


    /** \brief Take the composite of two functions
     *
     * @param f Function with which to compose
     * @return a new Function
     *//*
Function operator()(const Function & f) const ;

Function operator()(const Function &f0, const Function &f1) const;*/


    /** \brief Precalculate the function derivatives at a set of locations
    *
    * @param gp GaussPointArray defining the locations
    * @param var Variables with respect to which differentiation should be calculated
    * @param eps delata for the numerical evaluation of derivatives
    */
    void preCalculate ( const GaussPointArray & gp , std::vector<Variable> & var, const double eps = default_derivation_delta ) ;

    /** \brief Precalculate the function  at a set of locations
     *
     * @param gp GaussPointArray defining the locations
     */
    void preCalculate ( const GaussPointArray & gp ) ;

    /** \brief Check whether the function has been precalculated at a set of positions
     *
     * @param gp GaussPointArray defining the locations
     * @return true if the function has been precalculated
     */
    bool precalculated ( const GaussPointArray & gp ) const ;

    /** \brief Check whether the function derivative has been precalculated at a set of positions
     *
     * @param gp GaussPointArray defining the locations
     * @param v Variable with respect to which differentiation should be calculated
     * @return
     */
    bool precalculated ( const GaussPointArray & gp, Variable v ) const ;

    /** \brief return the precalculated values at a set of points
     *
     * @param gp GaussPointArray  defining the locations
     * @return Vector of values at those points
     */
    const Vector & getPrecalculatedValue ( const GaussPointArray & gp ) const ;

    /** \brief return the precalculated derivative values at a set of points
     *
     * @param gp GaussPointArray  defining the locations
     * @param v Variable with respect to which differentiation should be calculated
     * @return Vector of values at those points
     */
    const Vector & getPrecalculatedValue ( const GaussPointArray & gp, Variable v ) const ;
} ;

struct FtF ;
struct DtF ;
struct DDtF ;
struct GtM ;
struct GtML ;
struct GtV ;
struct GtVL ;
struct GtMtG ;
struct GtMLtG ;
struct VGtM ;
struct VGtV ;
struct VGDtM ;
struct VGDtMtVG ;
struct DtD ;
struct DtV ;
struct DtVL ;
struct VGtMtVG ;
struct VGtMtVGD ;
struct DtGtMtG ;
struct GDtM ;
struct GDtML ;
struct GDtV ;
struct GDtVL ;
struct GDtMtGD ;
struct GDtMtG ;
struct GDDtMtG ;
struct GDDtM ;
struct GDDtML ;
struct GtMtGD ;
struct GDtMLtGD ;
struct GDDtMLtG ;
struct GtMLtGD ;
struct GDtMLtG ;
struct DdGtMtG ;
struct DdGtMtGD ;

/** \brief Structure used for the lazy computation of the differential of a function
 * This is equivalent to a derivative, but can be used on all functions, even if their
 * derivatives have not been specified. In that case, a numerical derivative will be computed
*/
struct Differential {
    const Function & f ;
    const Variable & v ;

    /** \brief Constructor, initalise the references
     *
     * @param u Function
     * @param m Variable
     */
    Differential ( const Function &u, const Variable & m ) : f ( u ), v ( m ) { } ;

    /** \brief Create a structure for the lazy evaluation of the Differential * Function
     *
     * @param f factor
     * @return DtF
     */
    DtF operator * ( const Function & f ) const ;


    /** \brief Create a structure for the lazy evaluation of the Differential * Differential
     *
     * @param f factor
     * @return DtD
     */
    DtD operator * ( const Differential & f ) const ;

    DtV operator * ( const Vector & f ) const ;
    DtVL operator * ( const std::vector<Vector> & f ) const ;

    /** \brief Create a structure for the lazy evaluation of a Gradient * Matrix * Gradient * Differential
     *
     * @param g GtMtG
     * @return DtGtMtG
     */
    DtGtMtG operator * ( const GtMtG & g ) const ;
} ;

struct DoubleDifferential {
    const Function & f ;
    const Variable & v1 ;
    const Variable & v2 ;

    DoubleDifferential ( const Function & u, const Variable & v1_, const Variable & v2_ ) : f ( u ), v1 ( v1_ ), v2 ( v2_ ) { } ;

    DDtF operator* ( const Function & f ) const ;
} ;

/** \brief Create a structure for the lazy evaluation of a Gradient
 * The gradient of a function is evaluated as a Matrix, which dimensions are given by the
 * list of specifiend space variables given at the time of evaluation
*/
struct Gradient {
    const Function & f ;
    const bool transpose ;
    /** \brief Constructor, initalise the references
     *
     * @param u Function
     * @param t true if the result should be transposed
     */
    Gradient ( const Function &u, bool t = false ) : f ( u ), transpose ( t ) { };

    /** \brief Create a structure for the lazy evaluation of a Gradient * Matrix
     *
     * @param f Matrix
     * @return GtM
     */
    GtM operator * ( const Matrix & f ) const ;

    /** \brief Create a structure for the lazy evaluation of a Gradient * Matrix
     *
     * @param f Matrix
     * @return GtM
     */
    GtML operator * ( const std::vector<Matrix> & f ) const ;

    /** \brief Create a structure for the lazy evaluation of a Gradient * Vector
     *
     * @param f Vector
     * @return GtV
     */
    GtV operator * ( const Vector & f ) const ;

    /** \brief Create a structure for the lazy evaluation of a Gradient * Vector
     *
     * @param f Vector
     * @return GtVL
     */
    GtVL operator * ( const std::vector<Vector> & f ) const ;

} ;

/** \brief Create a structure for the lazy evaluation of a GradientDot
 * The gradient of a function is evaluated as a Matrix, which dimensions are given by the
 * list of specifiend space variables given at the time of evaluation. The result
 * is further differentiated with respect to time
*/
struct GradientDot {
    const Function & f ;
    const bool transpose ;

    /**  \brief Constructor, initalise the references
     *
     * @param u Function
     * @param t transpose result if true
     */
    GradientDot ( const Function &u, bool t = false ) : f ( u ), transpose ( t ) { };

    /** \brief Create a structure for the lazy evaluation of a GradientDot * Matrix
     *
     * @param f Matrix
     * @return GDtM
     */
    GDtM operator * ( const Matrix & f ) const ;
    GDtML operator * ( const std::vector<Matrix> & f ) const ;
    GDtV operator * ( const Vector & v ) const ;
    GDtVL operator * ( const std::vector<Vector> & f ) const ;


} ;

struct GradientDotDot {
    const Function & f ;
    const bool transpose ;

    /**  \brief Constructor, initalise the references
     *
     * @param u Function
     * @param t transpose result if true
     */
    GradientDotDot ( const Function &u, bool t = false ) : f ( u ), transpose ( t ) { };

    /** \brief Create a structure for the lazy evaluation of a GradientDot * Matrix
     *
     * @param f Matrix
     * @return GDtM
     */
    GDDtM operator * ( const Matrix & f ) const ;
    GDDtML operator * ( const std::vector<Matrix> & f ) const ;


} ;

/** \brief Create a structure for the lazy evaluation of a VectorGradient
 * The gradient of a function is evaluated as a single line Matrix, which dimension is given by the
 * list of specifiend space variables given at the time of evaluation.
*/
struct VectorGradient {
    const Function & f ;
    const bool transpose ;

    /**  \brief Constructor, initalise the references
     *
     * @param u Function
     * @param t transpose result if true
     */
    VectorGradient ( const Function &u, bool t = false ) : f ( u ), transpose ( t ) { };

    /** \brief Create a structure for the lazy evaluation of a VectorGradient * Matrix
     *
     * @param f Matrix
     * @return VGtM
     */
    VGtM operator * ( const Matrix & f ) const ;

    /** \brief Create a structure for the lazy evaluation of a VectorGradient * Vector
     *
     * @param f Vector
     * @return VGtV
     */
    VGtV operator * ( const Vector & f ) const ;
} ;

struct VectorGradientDot {
    const Function & f ;
    const bool transpose ;

    VectorGradientDot ( const Function &u, bool t = false ) : f ( u ), transpose ( t ) { } ;

    VGDtM operator * ( const Matrix & f ) const ;
} ;

struct FtF {
    const Function & first ;
    const Function & second ;

    FtF ( const Function &f, const Function &s ) : first ( f ), second ( s ) { } ;
} ;

/** \brief Structure for the lazy evaluation of a Gradient * Matrix structure*/
struct GtM {
    const Gradient & first ;
    const Matrix & second ;

    /** \brief Constructor, initalise the references
     *
     * @param g Gradient
     * @param f Matrix
     */
    GtM ( const Gradient & g, const Matrix & f ) : first ( g ), second ( f ) { };

    /** \brief Create a structure for the lazy evaluation of a Gradient * Matrix * Gradient
     *
     * @param f Gradient
     * @return GtMtG
     */
    GtMtG operator* ( const Amie::Gradient & f ) const ;

    /** \brief Create a structure for the lazy evaluation of a Gradient * Matrix * GradientDot
     *
     * @param f Gradient
     * @return GtMtGD
     */
    GtMtGD operator* ( const Amie::GradientDot & f ) const ;
} ;

/** \brief Structure for the lazy evaluation of a Gradient * Matrix at list of Gauss Points structure*/
struct GtML {
    const Gradient & first ;
    const std::vector<Matrix> & second ;

    /** \brief Constructor, initalise the references
     *
     * @param g Gradient
     * @param f std::vector<Matrix>
     */
    GtML ( const Gradient & g, const std::vector<Matrix> & f ) : first ( g ), second ( f ) { };

    /** \brief Create a structure for the lazy evaluation of a Gradient * Matrix * Gradient
     *
     * @param f Gradient
     * @return GtMtG
     */
    GtMLtG operator* ( const Amie::Gradient & f ) const ;

    GtMLtGD operator* ( const Amie::GradientDot & f ) const ;

} ;

/** \brief Structure for the lazy evaluation of a GradientDot * Matrix */
struct GDtM {
    const GradientDot & first ;
    const Matrix & second ;

    /** \brief Constructor, initalise the references
     *
     * @param g GradientDot
     * @param f Matrix
     */
    GDtM ( const GradientDot & g, const Matrix & f ) : first ( g ), second ( f ) { };

    /** \brief Create a structure for the lazy evaluation of a GradientDot * Matrix * GradientDot
     *
     * @param f GradientDot
     * @return GDtMtGD
     */
    GDtMtGD operator* ( const Amie::GradientDot & f ) const ;

    /** \brief Create a structure for the lazy evaluation of a GradientDot * Matrix * Gradient
     *
     * @param f Gradient
     * @return GDtMtGD
     */
    GDtMtG operator* ( const Amie::Gradient & f ) const ;
} ;

/** \brief Structure for the lazy evaluation of a GradientDot * Matrix */
struct GDtML {
    const GradientDot & first ;
    const std::vector<Matrix> & second ;

    /** \brief Constructor, initalise the references
     *
     * @param g GradientDot
     * @param f Matrix
     */
    GDtML ( const GradientDot & g, const std::vector<Matrix> & f ) : first ( g ), second ( f ) { };

    /** \brief Create a structure for the lazy evaluation of a GradientDot * Matrix * GradientDot
     *
     * @param f GradientDot
     * @return GDtMtGD
     */
    GDtMLtGD operator* ( const Amie::GradientDot & f ) const ;

    /** \brief Create a structure for the lazy evaluation of a GradientDot * Matrix * Gradient
     *
     * @param f Gradient
     * @return GDtMtGD
     */
    GDtMLtG operator* ( const Amie::Gradient & f ) const ;
} ;

/** \brief Structure for the lazy evaluation of a GradientDotDot * Matrix */
struct GDDtM {
    const GradientDotDot & first ;
    const Matrix & second ;

    /** \brief Constructor, initalise the references
     *
     * @param g GradientDot
     * @param f Matrix
     */
    GDDtM ( const GradientDotDot & g, const Matrix & f ) : first ( g ), second ( f ) { };

    /** \brief Create a structure for the lazy evaluation of a GradientDotDot * Matrix
     *
     * @param f GradientDot
     * @return GDtMtGD
     */
    GDDtMtG operator* ( const Amie::Gradient & f ) const ;
} ;

/** \brief Structure for the lazy evaluation of a GradientDotDot * Matrix */
struct GDDtML {
    const GradientDotDot & first ;
    const std::vector<Matrix> & second ;

    /** \brief Constructor, initalise the references
     *
     * @param g GradientDot
     * @param f Matrix
     */
    GDDtML ( const GradientDotDot & g, const std::vector<Matrix> & f ) : first ( g ), second ( f ) { };

    /** \brief Create a structure for the lazy evaluation of a GradientDotDot * Matrix
     *
     * @param f GradientDot
     * @return GDtMtGD
     */
    GDDtMLtG operator* ( const Amie::Gradient & f ) const ;
} ;

/** \brief Structure for the lazy evaluation of a GradientDot * Matrix * GradientDot */
struct GDtMtGD {
    const GradientDot & first ;
    const Matrix & second ;
    const GradientDot & third ;
    /** \brief Constructor, initalise the references
     *
     * @param g GDtM
     * @param f GradientDot
     */
    GDtMtGD ( const GDtM & g, const GradientDot & p ) : first ( g.first ), second ( g.second ), third ( p ) { };

} ;

/** \brief Structure for the lazy evaluation of a GradientDot * Matrix * GradientDot */
struct GDtMLtGD {
    const GradientDot & first ;
    const std::vector<Matrix> & second ;
    const GradientDot & third ;
    /** \brief Constructor, initalise the references
     *
     * @param g GDtM
     * @param f GradientDot
     */
    GDtMLtGD ( const GDtML & g, const GradientDot & p ) : first ( g.first ), second ( g.second ), third ( p ) { };

} ;

/** \brief Structure for the lazy evaluation of a GradientDot * Matrix * Gradient */
struct GDtMtG {
    const GradientDot & first ;
    const Matrix & second ;
    const Gradient & third ;
    /** \brief Constructor, initalise the references
     *
     * @param g GDtM
     * @param f GradientDot
     */
    GDtMtG ( const GDtM & g, const Gradient & p ) : first ( g.first ), second ( g.second ), third ( p ) { };

} ;

/** \brief Structure for the lazy evaluation of a GradientDot * Matrix * Gradient */
struct GDtMLtG {
    const GradientDot & first ;
    const std::vector<Matrix> & second ;
    const Gradient & third ;
    /** \brief Constructor, initalise the references
     *
     * @param g GDtM
     * @param f GradientDot
     */
    GDtMLtG ( const GDtML & g, const Gradient & p ) : first ( g.first ), second ( g.second ), third ( p ) { };

} ;

/** \brief Structure for the lazy evaluation of a GradientDotDot * Matrix * Gradient */
struct GDDtMtG {
    const GradientDotDot & first ;
    const Matrix & second ;
    const Gradient & third ;
    /** \brief Constructor, initalise the references
     *
     * @param g GDtM
     * @param f GradientDot
     */
    GDDtMtG ( const GDDtM & g, const Gradient & p ) : first ( g.first ), second ( g.second ), third ( p ) { };

} ;

/** \brief Structure for the lazy evaluation of a GradientDotDot * Matrix * Gradient */
struct GDDtMLtG {
    const GradientDotDot & first ;
    const std::vector<Matrix> & second ;
    const Gradient & third ;
    /** \brief Constructor, initalise the references
     *
     * @param g GDtM
     * @param f GradientDot
     */
    GDDtMLtG ( const GDDtML & g, const Gradient & p ) : first ( g.first ), second ( g.second ), third ( p ) { };

} ;

/** \brief Structure for the lazy evaluation of a GradientDot * Matrix * Gradient */
struct GtMLtGD {
    const Gradient & first ;
    const std::vector<Matrix> & second ;
    const GradientDot & third ;
    /** \brief Constructor, initalise the references
     *
     * @param g GtM
     * @param f GradientDot
     */
    GtMLtGD ( const GtML & g, const GradientDot & p ) : first ( g.first ), second ( g.second ), third ( p ) { };

} ;

/** \brief Structure for the lazy evaluation of a Differential * Function */
struct DtF {
    const Differential & d ;
    const Function & f ;

    /** \brief Constructor, initalise the references
     *
     * @param d_ Differential
     * @param f_ Function
     */
    DtF ( const Differential & d_, const Function & f_ ) : d ( d_ ), f ( f_ ) { } ;
} ;

struct DDtF {
    const DoubleDifferential & d ;
    const Function & f ;

    DDtF ( const DoubleDifferential & dd, const Function & u ) : d ( dd ), f ( u ) { } ;
} ;

/** \brief Structure for the lazy evaluation of a Differential * Differential */
struct DtD {
    const Differential & d ;
    const Differential & f ;

    /** \brief Constructor, initalise the references
     *
     * @param d_ Differential
     * @param f_ Differential
     */
    DtD ( const Differential & d_, const Differential & f_ ) : d ( d_ ), f ( f_ ) { } ;
} ;

struct DtV {
    const Differential & d ;
    const Vector & f ;

    /** \brief Constructor, initalise the references
     *
     * @param d_ Differential
     * @param f_ Differential
     */
    DtV ( const Differential & d_, const Vector & f_ ) : d ( d_ ), f ( f_ ) { } ;
} ;

struct DtVL {
    const Differential & d ;
    const std::vector<Vector> & f ;

    /** \brief Constructor, initalise the references
     *
     * @param d_ Differential
     * @param f_ Differential
     */
    DtVL ( const Differential & d_, const std::vector<Vector> & f_ ) : d ( d_ ), f ( f_ ) { } ;
} ;

/** \brief Structure for the lazy evaluation of a VectorGradient * Matrix */
struct VGtM {
    const VectorGradient & first ;
    const Matrix & second ;

    /** \brief Constructor, initalise the references
     *
     * @param g VectorGradient
     * @param f Matrix
     */
    VGtM ( const VectorGradient & g, const Matrix & f ) : first ( g ), second ( f ) { };

    /** \brief Create a structure for the lazy evaluation of a VectorGradient * Matrix * VectorGradient
     *
     * @param f VectorGradient
     * @return VGtMtVG
     */
    VGtMtVG operator* ( const Amie::VectorGradient & f ) const ;

    VGtMtVGD operator* ( const Amie::VectorGradientDot & f ) const ;
} ;

struct VGDtM {
    const VectorGradientDot & first ;
    const Matrix & second ;

    VGDtM ( const VectorGradientDot & g, const Matrix & f ) : first ( g ), second ( f ) { } ;

    VGDtMtVG operator* ( const Amie::VectorGradient & d ) const ;
} ;

/** \brief Structure for the lazy evaluation of a Gradient * Vector */
struct GtV {
    const Gradient & first ;
    const Vector & second ;

    /**  \brief Constructor, initalise the references
     *
     * @param g Gradient
     * @param f Vector
     */
    GtV ( const Gradient & g, const Vector & f ) : first ( g ), second ( f ) { };
} ;

/** \brief Structure for the lazy evaluation of a Gradient * Vector at series of points*/
struct GtVL {
    const Gradient & first ;
    const std::vector<Vector> & second ;

    /**  \brief Constructor, initalise the references
     *
     * @param g Gradient
     * @param f Vector
     */
    GtVL ( const Gradient & g, const std::vector<Vector> & f ) : first ( g ), second ( f ) { };
} ;

struct GDtV {
    const GradientDot & first ;
    const Vector & second ;

    /**  \brief Constructor, initalise the references
     *
     * @param g Gradient
     * @param f Vector
     */
    GDtV ( const GradientDot & g, const Vector & f ) : first ( g ), second ( f ) { };
} ;

struct GDtVL {
    const GradientDot & first ;
    const std::vector<Vector> & second ;

    /**  \brief Constructor, initalise the references
     *
     * @param g Gradient
     * @param f vector of Vector
     */
    GDtVL ( const GradientDot & g, const std::vector<Vector> & f ) : first ( g ), second ( f ) { };
} ;

/** \brief Structure for the lazy evaluation of a VectorGradient * Vector */
struct VGtV {
    const VectorGradient & first ;
    const Vector & second ;

    /**  \brief Constructor, initalise the references
     *
     * @param g VectorGradient
     * @param f Vector
     */
    VGtV ( const VectorGradient & g, const Vector & f ) : first ( g ), second ( f ) { };
} ;

struct DdGtMtG {
    const Variable first ;
    const GtMtG & second ;

    DdGtMtG ( Variable v, const GtMtG & g ) : first ( v ), second ( g ) { } ;
} ;

struct DdGtMtGD {
    const Variable first ;
    const GtMtGD & second ;

    DdGtMtGD ( Variable v, const GtMtGD & g ) : first ( v ), second ( g ) { } ;
} ;

/** \brief Structure for the lazy evaluation of a Gradient * Matrix * Gradient*/
struct GtMtG {
    const Gradient & first ;
    const Matrix & second ;
    const Gradient & third ;

    /**  \brief Constructor, initalise the references
     *
     * @param g Gradient
     * @param f Matrix
     * @param g_ Gradient
     */
    GtMtG ( const Gradient & g, const Matrix & f,const Gradient & g_ ) : first ( g ), second ( f ), third ( g_ ) { };

    DdGtMtG operator() ( Variable v ) {
        return DdGtMtG ( v, *this ) ;
    }

} ;

/** \brief Structure for the lazy evaluation of a GradientDot * Matrix * Gradient */
struct GtMtGD {
    const Gradient & first ;
    const Matrix & second ;
    const GradientDot & third ;
    /** \brief Constructor, initalise the references
     *
     * @param g GtM
     * @param f GradientDot
     */
    GtMtGD ( const GtM & g, const GradientDot & p ) : first ( g.first ), second ( g.second ), third ( p ) { };

    DdGtMtGD operator() ( Variable v ) {
        return DdGtMtGD ( v, *this ) ;
    }

} ;


/** \brief Structure for the lazy evaluation of a Gradient * Matrix * Gradient at a list of Gauss Points*/
struct GtMLtG {
    const Gradient & first ;
    const std::vector<Matrix> & second ;
    const Gradient & third ;

    /**  \brief Constructor, initalise the references
     *
     * @param g Gradient
     * @param f std::vector<Matrix>
     * @param g_ Gradient
     */
    GtMLtG ( const Gradient & g, const std::vector<Matrix> & f,const Gradient & g_ ) : first ( g ), second ( f ), third ( g_ ) { };

} ;

/** \brief Structure for the lazy evaluation of a Differential * Gradient * Matrix * Gradient*/
struct DtGtMtG {
    const Differential & first ;
    const GtMtG & second ;

    /** \brief Constructor, initalise the references
     *
     * @param d Differential
     * @param g GtMtG
     */
    DtGtMtG ( const Differential & d, const GtMtG & g ) : first ( d ), second ( g ) {  };
} ;

/** \brief Structure for the lazy evaluation of a VectorGradient * Matrix * VectorGradient*/
struct VGtMtVG {
    const VectorGradient & first ;
    const Matrix & second ;
    const VectorGradient & third ;

    /** \brief Constructor, initalise the references
     *
     * @param g VectorGradient
     * @param f Matrix
     * @param g_ VectorGradient
     */
    VGtMtVG ( const VectorGradient & g, const Matrix & f,const VectorGradient & g_ ) : first ( g ), second ( f ), third ( g_ ) { };

} ;

struct VGDtMtVG {
    const VectorGradientDot & first ;
    const Matrix & second ;
    const VectorGradient & third ;

    VGDtMtVG ( const VectorGradientDot & gd, const Matrix & f, const VectorGradient & g ) : first ( gd ), second ( f ), third ( g ) { } ;
} ;

struct VGtMtVGD {
    const VectorGradient & first ;
    const Matrix & second ;
    const VectorGradientDot & third ;

    VGtMtVGD ( const VectorGradient & g, const Matrix & f, const VectorGradientDot & gd ) : first ( g ), second ( f ), third ( gd ) { } ;
} ;

void concatenateFunctions(const Function & src0, const Function & src1, Function & dst) ;

void concatenateFunctions(const Function & src0_, const Function & src1_, const Function & src2_, Function & dst) ;


}

#endif

