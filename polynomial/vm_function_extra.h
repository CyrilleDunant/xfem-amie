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

/** \brief Helper function to create a Function which is the square root of the argument
 * 
 * @param f Function
 * @return a new Function
 */
Amie::Function f_sqrt(const Amie::Function &f, bool differentiate = true ) ;

/** \brief Helper function to create a Function which is the exponential of the argument
 * 
 * @param f Function
 * @return a new Function
 */
Amie::Function f_exp(const Amie::Function &f) ;

/** \brief Helper function to create a Function which is the absolute value of the argument
 * 
 * @param f Function
 * @return a new Function
 */
Amie::Function f_abs(const Amie::Function &f, bool differentiate = true) ;

/**  \brief Helper function to create a Function which is the log of the argument
 * 
 * @param f Function
 * @return a new Function
 */
Amie::Function f_log(const Amie::Function &f) ;

/** \brief Helper function to create a Function which is arc tangeant of the arguments
 * 
 * @param f0 Function
 * @param f1 Function
 * @return a new Function
 */
Amie::Function f_atan2(const Amie::Function &f0, const Amie::Function &f1) ;

/** \brief Helper function to create a Function which is the sine of the argument
 * 
 * @param f Function
 * @return a new Function
 */
Amie::Function f_sin(const Amie::Function &f) ;

/** \brief Helper function to create a Function which is the cosine of the argument
 * 
 * @param f Function
 * @return a new Function
 */
Amie::Function f_cos(const Amie::Function &f) ;

/** \brief Helper function to create a Function which is the sign of the argument
 * 
 * @param f Function
 * @return a new Function
 */
Amie::Function f_sign(const Amie::Function &f) ;

/** \brief Helper function to create a Function which is 1 if the argument is positive, 0 otherwise
 * 
 * @param f Function
 * @return a new Function
 */
Amie::Function f_positivity(const Amie::Function& f, bool differentiate = true) ;

/** \brief Helper function to create a Function which is 1 if the argument is negative, 0 otherwise
 * 
 * @param f Function
 * @return a new Function
 */
Amie::Function f_negativity(const Amie::Function& f, bool differentiate = true) ;

/** \brief  Helper function to create a Function which is the interpolated value between two functions. 0 as an argument will yield the value of the first function and 1 the second, intermediate values will yield the weighted average at that point.
 * 
 * @param f0 Function
 * @param f1 Function
 * @return a new Function
 */
Amie::Function f_interpolate(const Amie::Function &f0, const Amie::Function &f1) ;

/** \brief Helper function to create a Function which computes the distance between a point, given by the arguments and its pojection on a given Geometry. The coordinates are transformed using the provided functions
 * 
 * @param g Geometry
 * @param x x transform Function
 * @param y y transform Function
 * @return a new Function
 */
Amie::Function f_project(const Amie::Geometry *g, const Amie::Function &x, const Amie::Function &y) ;

/** \brief Helper function to compute the curvilinear absciss on a SegmentedLine
 *
 * @param s SegmentedLine serving as source for the computation
 * @param fromHead compute starting from the Head
 * @param x x transform Function
 * @param y y transform Function
 * @return a new Function
*/
Amie::Function f_curvilinear_x(const Amie::SegmentedLine * s, bool fromHead,   const Amie::Function &x, const Amie::Function &y) ;

/** \brief Helper function to compute the curvilinear ordinate on a SegmentedLine
 *
 * @param s SegmentedLine serving as source for the computation
 * @param fromHead compute starting from the Head
 * @param x x transform Function
 * @param y y transform Function
 * @return a new Function
*/
Amie::Function f_curvilinear_y(const Amie::SegmentedLine * s, bool fromHead,   const Amie::Function &x, const Amie::Function &y) ;


#endif
