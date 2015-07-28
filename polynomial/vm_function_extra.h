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

#ifndef VM_FUNCTION_EXTRA_H
#define VM_FUNCTION_EXTRA_H


#include "vm_token.h"
#include "vm_base.h"

namespace Amie {
/** \brief Helper function to create a Function which is the square root of the argument
 * 
 * @param f Function
 * @return a new Function
 */
Function f_sqrt(const Function &f, bool differentiate = true ) ;

/** \brief Helper function to create a Function which is the exponential of the argument
 * 
 * @param f Function
 * @return a new Function
 */
Function f_exp(const Function &f) ;

/** \brief Helper function to create a Function which is the absolute value of the argument
 * 
 * @param f Function
 * @return a new Function
 */
Function f_abs(const Function &f, bool differentiate = true) ;

/**  \brief Helper function to create a Function which is the log of the argument
 * 
 * @param f Function
 * @return a new Function
 */
Function f_log(const Function &f) ;

/** \brief Helper function to create a Function which is arc tangeant of the arguments
 * 
 * @param f0 Function
 * @param f1 Function
 * @return a new Function
 */
Function f_atan2(const Function &f0, const Function &f1) ;

/** \brief Helper function to create a Function which is the sine of the argument
 * 
 * @param f Function
 * @return a new Function
 */
Function f_sin(const Function &f) ;

/** \brief Helper function to create a Function which is the cosine of the argument
 * 
 * @param f Function
 * @return a new Function
 */
Function f_cos(const Function &f) ;

/** \brief Helper function to create a Function which is the hyperbolic sine of the argument
 * 
 * @param f Function
 * @return a new Function
 */
Function f_sinh(const Function &f) ;

/** \brief Helper function to create a Function which is the hyperbolic cosine of the argument
 * 
 * @param f Function
 * @return a new Function
 */
Function f_cosh(const Function &f) ;

/** \brief Helper function to create a Function which is the tangent of the argument
 * 
 * @param f Function
 * @return a new Function
 */
Function f_tan(const Function &f) ;

/** \brief Helper function to create a Function which is the hyperbolic tangeant of the argument
 * 
 * @param f Function
 * @return a new Function
 */
Function f_tanh(const Function &f) ;

/** \brief Helper function to create a Function which is the sign of the argument
 * 
 * @param f Function
 * @return a new Function
 */
Function f_sign(const Function &f) ;

/** \brief Helper function to create a Function which is 1 if the argument is positive, 0 otherwise
 * 
 * @param f Function
 * @return a new Function
 */
Function f_positivity(const Function& f, bool differentiate = true) ;

/** \brief Helper function to create a Function which is 1 if the argument is negative, 0 otherwise
 * 
 * @param f Function
 * @return a new Function
 */
Function f_negativity(const Function& f, bool differentiate = true) ;


Function f_range(const Function &f, double min, double max, bool differentiate = false);


/** \brief  Helper function to create a Function which is the interpolated value between two functions. 0 as an argument will yield the value of the first function and 1 the second, intermediate values will yield the weighted average at that point.
 * 
 * @param f0 Function
 * @param f1 Function
 * @return a new Function
 */
Function f_interpolate(const Function &f0, const Function &f1) ;

/** \brief Helper function to create a Function which computes the distance between a point, given by the arguments and its pojection on a given Geometry. The coordinates are transformed using the provided functions
 * 
 * @param g Geometry
 * @param x x transform Function
 * @param y y transform Function
 * @return a new Function
 */
Function f_project(const Geometry *g, const Function &x, const Function &y) ;

/** \brief Helper function to compute the curvilinear absciss on a SegmentedLine
 *
 * @param s SegmentedLine serving as source for the computation
 * @param fromHead compute starting from the Head
 * @param x x transform Function
 * @param y y transform Function
 * @return a new Function
*/
Function f_curvilinear_x(const SegmentedLine * s, bool fromHead,   const Function &x, const Function &y) ;

/** \brief Helper function to compute the curvilinear ordinate on a SegmentedLine
 *
 * @param s SegmentedLine serving as source for the computation
 * @param fromHead compute starting from the Head
 * @param x x transform Function
 * @param y y transform Function
 * @return a new Function
*/
Function f_curvilinear_y(const SegmentedLine * s, bool fromHead,   const Function &x, const Function &y) ;

}

#endif

