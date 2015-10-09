// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef VM_BASE_H
#define VM_BASE_H

#include "../elements/integrable_entity.h"
#include "vm_function_base.h"
#include "vm_function_matrix.h"

namespace Amie
{

struct FtF ;
struct GtFMtG ;
struct VGtMtVG ;
struct VGDtMtVG ;
struct VGtMtVGD ;
struct GtM ;
struct GtML ;
struct GtMtG ;
struct GtMLtG ;
struct GtV ;
struct GtVL ;
struct VGtM ;
struct VGtV ;
struct FMtMtFM ;
struct DtGtMtG ;
struct DdGtMtG ;
struct DdGtMtGD ;
struct DdGtMLtG ;
struct DdGtMLtGD ;
struct DdGtMLtG ;
struct DdGtMLtGD ;
struct GDtMtGD ;
struct GDDtMtG ;
struct GtMtGDD ;
struct GtMtGD ;
struct GDtMtG ;
struct GDtMLtGD ;
struct GDDtMLtG ;
struct GtMLtGDD ;
struct GtMLtGD ;
struct GDtMLtG ;
struct GDtV ;
struct GDtVL ;
struct IntegrableEntity ;
struct GaussPointArray ;
class Function ;
struct DeprecatedFunction ;
class FunctionMatrix ;
struct Differential ;
struct DtF ;
struct DtD ;
struct DDtF ;
struct DtV ;
struct DtVL ;
struct Gradient ;
struct GradientDot ;
struct VectorGradient ;
struct VectorGradientDot ;

/** \brief A virtual machine class for computation of symbolic formulas.
 *
 * The VirtualMachine is a fully-descending stack-based machine optimised for mathematical Function computation in the context
 * of finite element modeling. Integro-differential operators are provided. Caching of temporaties used across a typical simulation
 * run is provided to minimise memory iniitalisation and release.
 */
class VirtualMachine
{
    Context stack ;
    Matrix B ;
    Matrix B_ ;
public:

    /** \brief Constructor.
    The constructor created and initialises the stack
    */
    VirtualMachine() ;

    /** \brief Evaluate Function f at (x, y, z, t, u, v, w).
    @param f Function to evaluate.
    @param x coordinate at which to evaluate.
    @param y coordinate at which to evaluate.
    @param z coordinate at which to evaluate.
    @param t coordinate at which to evaluate.
    @param u coordinate at which to evaluate.
    @param v coordinate at which to evaluate.
    @param w coordinate at which to evaluate.
    */
    double eval(const Function &f, const double x = 0, const double y = 0, const double z = 0 ,const double t = 0, const double u = 0, const double v = 0, const double w = 0) ;

    /** \brief Evaluate FunctionMatrix f at (x, 0...).
    @param f FunctionMatrix to evaluate.
    @param x coordinate at which to evaluate.
    @param y coordinate at which to evaluate.
    @param z coordinate at which to evaluate.
    @param t coordinate at which to evaluate.
    @param u coordinate at which to evaluate.
    @param v coordinate at which to evaluate.
    @param w coordinate at which to evaluate.
    */
    Matrix eval(const FunctionMatrix &f, const double x = 0 , const double y = 0, const double z = 0,const double t = 0, const double u = 0, const double v = 0, const double w = 0)  ;

    /** \brief Evaluate Function f at (p.getX(), p.getY(), p.getZ(), p.getT(), p_.getX(), p_.getY(), p_.getZ(),).
    @param f Function to evaluate.
    @param p Point giving the first four coordinates at which to evaluate.
    @param p_ Point givint the last three coordinates at which to evaluate.
    */
    double eval(const Function &f, const Point & p, const Point &p_ = Point())  ;

    /** \brief Evaluate Function f at (p->getX(), p->getY(), p->getZ(), p->getT(), p_->getX(), p_->getY(), p_->getZ(),).
    @param f Function to evaluate.
    @param p Point giving the first four coordinates at which to evaluate.
    @param p_ Point givint the last three coordinates at which to evaluate.
    */
    double eval(const Function &f, const Point * p, const Point *p_ = nullptr)  ;

    /** \brief Evaluate FunctionMatrix f at (p.getX(), p.getY(), p.getZ(), p.getT(), p_.getX(), p_.getY(), p_.getZ(),).
    @param f FunctionMatrix to evaluate.
    @param p Point giving the first four coordinates at which to evaluate.
    @param p_ Point givint the last three coordinates at which to evaluate.
    */
    Matrix eval(const FunctionMatrix &f, const Point & p, const Point &p_ = Point())  ;

    /** \brief Evaluate FunctionMatrix f at (p->getX(), p->getY(), p->getZ(), p->getT(), p_->getX(), p_->getY(), p_->getZ(),).
    @param f FunctionMatrix to evaluate.
    @param p Point giving the first four coordinates at which to evaluate.
    @param p_ Point givint the last three coordinates at which to evaluate.
    */
    Matrix eval(const FunctionMatrix &f, const Point * p, const Point *p_ = nullptr)  ;

    /** \brief Evaluate Function f at all points of the given GaussPointArray gp
    @param f Function to evaluate.
    @param gp GaussPointArray giving the points at which to evaluate the Function f
    */
    Vector eval(const Function &f, const GaussPointArray &gp) ;
    Vector eval(const DeprecatedFunction &f, const GaussPointArray &gp) ;

    /** \brief Evaluate Function f at all points of the given GaussPointArray gp, each offset by "offset".
    This method is used internally to compute the numerical differential of a set of points.
    @param f Function to evaluate.
    @param gp GaussPointArray giving the points at which to evaluate the Function f
    @param offset Point used to compute the offset
    */
    Vector oeval(const Function &f, const GaussPointArray &gp, const Point & offset) ;

    /** \brief Evaluate df/dv_ at (x, y, z, t, u, v, w) where f is a Function.
    @param f Function to evaluate.
    @param v_ Variable with respect to which the differential should be calculated
    @param x coordinate at which to evaluate.
    @param y coordinate at which to evaluate.
    @param z coordinate at which to evaluate.
    @param t coordinate at which to evaluate.
    @param u coordinate at which to evaluate.
    @param v coordinate at which to evaluate.
    @param w coordinate at which to evaluate.
    @param eps is the epsilon used to perform the numerical differential.
    */
    double deval(const Function &f, const Variable v_,  const double x, const double y = 0, const double z = 0,const double t=0, const double u = 0, const double v = 0, const double w = 0, const double eps = default_derivation_delta)  ;

    /** \brief Evaluate df/dv_ at (x, y, z, t, u, v, w), where f is a FunctionMatrix
    @param f FunctionMatrix to evaluate.
    @param v_ Variable with respect to which the differential should be calculated
    @param x coordinate at which to evaluate.
    @param y coordinate at which to evaluate.
    @param z coordinate at which to evaluate.
    @param t coordinate at which to evaluate.
    @param u coordinate at which to evaluate.
    @param v coordinate at which to evaluate.
    @param w coordinate at which to evaluate.
     * @param eps is the epsilon used to perform the numerical differential.
     */
    Matrix deval(const FunctionMatrix &f, const Variable v_,  const double x, const double y = 0, const double z = 0,const double t=0, const double u = 0, const double v = 0, const double w = 0, const double eps= default_derivation_delta)  ;

    /** \brief Evaluate df/dv at (p.getX(), p.getY(), p.getZ(), p.getT(), p_.getX(), p_.getY(), p_.getZ(),).
    @param f Function to evaluate.
    @param v Variable with respect to which the differential should be calculated
    @param p Point giving the first four coordinates at which to evaluate.
    @param p_ Point givint the last three coordinates at which to evaluate.
    @param eps is the epsilon used to perform the numerical differential.
    */
    double deval(const Function &f, const Variable v, const Point & p, const Point & p_ = Point(), const double eps= default_derivation_delta)  ;

    /** \brief Evaluate df/dv at (p.getX(), p.getY(), p.getZ(), p.getT(), p_.getX(), p_.getY(), p_.getZ(),).
    @param f FunctionMatrix to evaluate.
    @param v Variable with respect to which the differential should be calculated
    @param p Point giving the first four coordinates at which to evaluate.
    @param p_ Point givint the last three coordinates at which to evaluate.
    @param eps is the epsilon used to perform the numerical differential.
    */
    Matrix deval(const FunctionMatrix &f, const Variable v, const Point & p, const Point & p_ = Point(), const double eps= default_derivation_delta)  ;

    /** \brief Evaluate \f$ (df/dx, df/dy, df/dz) \cdot (p.getX(), p.getY(), p.getZ())\f$ at (x, y, z, t).
    This is the directional derivative.
    @param f Function to evaluate.
    @param p Direction on which the derivative vector should be projected
    @param x coordinate at which to evaluate.
    @param y coordinate at which to evaluate.
    @param z coordinate at which to evaluate.
    @param t coordinate at which to evaluate.
    @param eps is the epsilon used to perform the numerical differential.
    @param normed true if the direction vector should be normalised
    */
    double deval(const Function &f, const Point&p,  const double x, const double y = 0, const double z = 0, const double t = 0,const double eps = default_derivation_delta, bool normed = false)  ;

    /** \brief Evaluate \f$ (df/dx, df/dy, df/dz) \cdot (p.getX(), p.getY(), p.getZ()) \f$ at (x, y, z, t) for each element of the FunctionMatrix.
    This is the directional derivative.
    @param f FunctionMatrix to evaluate.
    @param p Direction on which the derivative vector should be projected
    @param x coordinate at which to evaluate.
    @param y coordinate at which to evaluate.
    @param z coordinate at which to evaluate.
    @param t coordinate at which to evaluate.
    @param eps is the epsilon used to perform the numerical differential.
    @param normed true if the direction vector should be normalised
    */
    Matrix deval(const FunctionMatrix &f,const Point&p,  const double x, const double y = 0, const double z = 0, const double t = 0,const double eps= default_derivation_delta, bool normed = false)  ;

    /** \brief Evaluate \f$ (df/dx, df/dy, df/dz) \cdot (p_.getX(), p_.getY(), p_.getZ()) \f$ at (p.getX(), p.getY(), p.getZ(), p.getT()).
    This is the directional derivative.
    @param f Function to evaluate.
    @param p_ Direction on which the derivative vector should be projected
    @param p Point at which to evaluate
    * @param eps is the epsilon used to perform the numerical differential.
    * @param normed normalise the direction vector if set to true.
    */
    double deval(const Function &f,  const Point&p_, const Point & p, const double eps= default_derivation_delta, bool normed = false)  ;

    /** \brief Evaluate \f$ (df/dx, df/dy, df/dz) \cdot (p_.getX(), p_.getY(), p_.getZ()) \f$ at (p.getX(), p.getY(), p.getZ(), p.getT()).
    This is the directional derivative.
    @param f FunctionMatrix to evaluate.
    @param p_ Direction on which the derivative vector should be projected
    @param p Point at which to evaluate
    * @param eps is the epsilon used to perform the numerical differential.
    * @param normed normalise the direction vector if set to true.
    */
    Matrix deval(const FunctionMatrix &f,  const Point&p_, const Point & p, const double eps = default_derivation_delta, bool normed = false)  ;

    /** \brief Evaluate d²f/d(v_0v_1) at (x, y, z, t, u, v, w).
    @param f Function to evaluate.
    @param v_ first Variable with respect to which the differential should be calculated
    @param v_ second Variable with respect to which the differential should be calculated
    @param x coordinate at which to evaluate.
    @param y coordinate at which to evaluate.
    @param z coordinate at which to evaluate.
    @param t coordinate at which to evaluate.
    @param u coordinate at which to evaluate.
    @param v coordinate at which to evaluate.
    @param w coordinate at which to evaluate.
    @param eps is the epsilon used to perform the numerical differential.
     */
    double ddeval(const Function &f, const Variable v_0, const Variable v_1,  const double x=0, const double y=0 , const double z=0, const double t=0, const double u=0, const double v=0, const double w=0, const double eps = default_derivation_delta) ;

    double ddeval(const Function &f, const Variable v_0, const Variable v_1,  const Point & p, const double eps = default_derivation_delta) ;

    /** \brief Evaluate d²f/d(v_0v_1) at (x, y, z, t, u, v, w) for each point for the GaussPointArray.
    @param f Function to evaluate.
    @param v_0 first Variable with respect to which the differential should be calculated
    @param v_1 second Variable with respect to which the differential should be calculated
    @param gp GaussPointArray defining the coordinates at which to evaluate.
     * @param eps is the epsilon used to perform the numerical differential.
     */
    Vector ddeval(const Function&f, const Variable v_0, const Variable v_1, const GaussPointArray & gp, const double eps= default_derivation_delta) ;

    double dddeval(const Function &f, const Variable v_0, const Variable v_1,  const Variable v_2, const double x=0, const double y=0 , const double z=0, const double t=0, const double u=0, const double v=0, const double w=0, const double eps = default_derivation_delta) ;

    Vector dddeval(const Function&f, const Variable v_0, const Variable v_1, const Variable v_2, const GaussPointArray & gp, const double eps= default_derivation_delta) ;

    /** \brief Evaluate df/dv_ at (x, y, z, t, u, v, w) for each point for the GaussPointArray.
    @param f Function to evaluate.
    @param v_ Variable with respect to which the differential should be calculated
    @param gp GaussPointArray defining the coordinates at which to evaluate.
    * @param eps is the epsilon used to perform the numerical differential.
    */
    Vector deval(const Function&f, const Variable v_, const GaussPointArray & gp, const double eps= default_derivation_delta) ;

    /** \brief Compute the numerical integration of function f over the IntegrableEntity e.
    The function is assumed to be expressed in the local coordinates of e.
    @param f Function to integrate.
    @param e IntegrableEntity on which to integrate
    */
    double ieval(const Function &f, IntegrableEntity *e)  ;

    Matrix ieval(const FtF & f, const GaussPointArray &gp, const std::valarray<Matrix> & Jinv) ;

    /** \brief Compute the numerical integration of function f using the GaussPointArray to perform the quadrature.
    @param f Function to integrate.
    @param gp GaussPointArray defining the quadrature.
    */
    double ieval(const Function &f, const GaussPointArray &gp)  ;

    /** \brief Compute the sum of the values of f, weighted by the quadrature given by the IntegrableEntity e.
    The function is assumed to be expressed in the local coordinates of e.
    @param f Vector of values to sum.
    @param e IntegrableEntity on which to integrate
    */
    double ieval(Vector &f, IntegrableEntity *e) ;

    /** \brief Compute the sum of the values of f, weighted by the quadrature given by the IntegrableEntity e.
    The function is assumed to be expressed in the local coordinates of e.
    @param f Vector of values to sum.
    @param e IntegrableEntity on which to integrate
    */
    Vector ieval(const std::vector<Vector> &f, IntegrableEntity *e) ;


    /** \brief Compute the numerical integration of FunctionMatrix f over the IntegrableEntity e.
    The function is assumed to be expressed in the local coordinates of e.
    @param f FunctionMatrix to integrate.
    @param e IntegrableEntity on which to integrate
    */
    Matrix ieval(const FunctionMatrix &f, IntegrableEntity *e)  ;

    /** \brief Overloaded function to compute the integral of a FunctionMatrix times a Matrix times a second FunctionMatrix over the IntegrableEntity e.
    The function is assumed to be expressed in the local coordinates of e.
    @param f FMtMtFM to integrate.
    @param e IntegrableEntity on which to integrate
    */
    Matrix ieval(const FMtMtFM &f, IntegrableEntity *e) ;

    /** \brief Overloaded function to compute the integral of a Gradient times a Matrix  over the IntegrableEntity e, with variables defined by vars.
    @param f GtM to integrate.
    @param e IntegrableEntity on which to integrate
    @param vars std::vector of space Variable s
    The function is assumed to be expressed in the local coordinates of e.
    */
    Matrix ieval(const GtM &f, IntegrableEntity *e, const std::vector<Variable> & vars) ;

    Matrix ieval(const GtML &f, IntegrableEntity *e, const std::vector<Variable> & vars) ;

    Vector ieval(const GtVL &f, IntegrableEntity *e, const std::vector<Variable> & vars) ;

    Vector ieval(const GtV &f, IntegrableEntity *e, const std::vector<Variable> & vars) ;

    Vector ieval(const GDtV &f, IntegrableEntity *e, const std::vector<Variable> & vars) ;

    Vector ieval(const GDtVL &f, IntegrableEntity *e, const std::vector<Variable> & vars) ;

    /** \brief Overloaded function to compute the integral of a Gradient times a Matrix times a Gradient over the IntegrableEntity e, with variables defined by vars.
    The function is assumed to be expressed in the local coordinates of e.
    Gradient is the usual \f$ \nabla\otimes \f$ operator.
    @param f GtMtG to integrate.
    @param e IntegrableEntity on which to integrate
    @param vars std::vector of space Variable s
    */
    Matrix ieval(const GtMtG &f, IntegrableEntity *e, const std::vector<Variable> & vars) ;

    /** \brief Overloaded function to compute the integral of a Gradient times a Matrix times a Gradient using the inverse Jacobian matrices given by Jinv and the Gauss points in gp, with variables defined by vars.
    Gradient is the usual \f$ \nabla\otimes \f$ operator.
    @param f GtMtG to integrate.
    @param gp GaussPointArray defining the quadrature
    @param Jinv Inverse Jacobian Matrices to compute the gradients
    @param vars std::vector of space Variable s
    */
    Matrix ieval(const GtMtG &f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars) ;

    std::valarray<Matrix> ievalDecomposed(const GtMtG &f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars) ;

    std::valarray<Matrix> ievalDecomposed(const GtMtGD &f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars) ;
    std::valarray<Matrix> ievalDecomposedDebug(const GtMtGD &f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars) ;

    /** \brief Overloaded function to compute the integral of a Gradient times a Matrix times a Gradient using the inverse Jacobian matrices given by Jinv and the Gauss points in gp, with variables defined by vars.
    The result is stored in ret. The version of the function can be used to minimise initialisation of memory.
    Gradient is the usual \f$ \nabla\otimes \f$ operator.
    @param f GtMtG to integrate.
    @param gp GaussPointArray defining the quadrature
    @param Jinv Inverse Jacobian Matrices to compute the gradients
    @param vars std::vector of space Variable s
    */
    void ieval(const GtMtG &f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars, Matrix & ret) ;

    /** \brief Overloaded function to compute the integral of a Gradient times a Matrix times a Gradient using the inverse Jacobian matrices given by Jinv and the Gauss points in gp, with variables defined by vars.
    Gradient is the usual \f$ \nabla\otimes \f$ operator.
    @param f GtMtG to integrate.
    @param gp GaussPointArray defining the quadrature
    @param Jinv Inverse Jacobian Matrices to compute the gradients
    @param vars std::vector of space Variable s
    */
    Matrix ieval(const GtMLtG &f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars) ;

    /** \brief Overloaded function to compute the integral of a Gradient times a Matrix times a Gradient using the inverse Jacobian matrices given by Jinv and the Gauss points in gp, with variables defined by vars.
    The result is stored in ret. The version of the function can be used to minimise initialisation of memory.
    Gradient is the usual \f$ \nabla\otimes \f$ operator.
    @param f GtMtG to integrate.
    @param gp GaussPointArray defining the quadrature
    @param Jinv Inverse Jacobian Matrices to compute the gradients
    @param vars std::vector of space Variable s
    */
    void ieval(const GtMLtG &f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars, Matrix & ret) ;

    /** \brief Overloaded function to compute the integral of a Gradient times a Vector using the inverse Jacobian matrices given by Jinv and the Gauss points in gp, with variables defined by vars.
    Gradient is the usual \f$ \nabla\otimes \f$ operator.
    @param f GtV to integrate.
    @param gp GaussPointArray defining the quadrature
    @param Jinv Inverse Jacobian Matrices to compute the gradients
    @param vars std::vector of space Variable s
    */
    Vector ieval(const GtV &f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars) ;

    Vector ieval(const GtVL &f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars) ;
    Vector ieval(const GDtVL &f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars) ;


    Vector ieval(const GDtV &f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars) ;

    /** \brief Overloaded function to compute the integral of a Gradient times a FunctionMatrix times a Gradient over the IntegrableEntity e, with variables defined by vars.
    Gradient is the usual \f$ \nabla\otimes \f$ operator.
    @param f GtFMtG to integrate.
    @param e IntegrableEntity on which to integrate
    @param vars std::vector of space Variable s
    */
    Matrix ieval(const GtFMtG &f, IntegrableEntity *e, const std::vector<Variable> & vars) ;

    /** \brief Overloaded function to compute the integral of a Gradient times a FunctionMatrix times a Gradient using the inverse Jacobian matrices given by Jinv and the Gauss points in gp, with variables defined by vars.
    Gradient is the usual \f$ \nabla\otimes \f$ operator.
    @param f GtV to integrate.
    @param gp GaussPointArray defining the quadrature
    @param Jinv Inverse Jacobian Matrices to compute the gradients
    @param vars std::vector of space Variable s
    */
    Matrix ieval(const GtFMtG &f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars) ;

    /** \brief Overloaded function to compute the integral of a Gradient times a FunctionMatrix times a Gradient using the inverse Jacobian matrices given by Jinv and the Gauss points in gp, with variables defined by vars.
    The result is stored in ret. The version of the function can be used to minimise initialisation of memory.
    Gradient is the usual \f$ \nabla \otimes \f$ operator.
    @param f GtFMtG to integrate.
    @param gp GaussPointArray defining the quadrature
    @param Jinv Inverse Jacobian Matrices to compute the gradients
    @param vars std::vector of space Variable s
    @param ret Matrix in which to store the results
    */
    void ieval(const GtFMtG &f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars, Matrix & ret) ;

    /** \brief Overloaded function to compute the integral of a GradientDot times a Matrix times a Gradient over the IntegrableEntity e, with variables defined by vars.
    Gradient is the usual \f$ \nabla\otimes \f$ operator.
    @param f GDtMtG to integrate.
    @param e IntegrableEntity on which to integrate
    @param vars std::vector of space Variable s
    */
    Matrix ieval(const GDtMtG & f, IntegrableEntity * e , const std::vector<Variable> & vars) ;


    /** \brief Overloaded function to compute the integral of a Gradient times a Matrix times a GradientDot over the IntegrableEntity e, with variables defined by vars.
    Gradient is the usual \f$ \nabla\otimes \f$ operator.
    @param f GtMtGD to integrate.
    @param e IntegrableEntity on which to integrate
    @param vars std::vector of space Variable s
    */
    Matrix ieval(const GtMtGD & f, IntegrableEntity * e , const std::vector<Variable> & vars) ;

    /** \brief Overloaded function to compute the integral of a GradientDot times a FunctionMatrix times a GradientDot using the inverse Jacobian matrices given by Jinv and the Gauss points in gp, with variables defined by vars.
     * The result is stored in ret. The version of the function can be used to minimise initialisation of memory.
     * GradientDot is the operator \f$ \dot{\nabla}\otimes(\cdot) \f$ operator.
    @param f GDtMtGD to integrate.
    @param gp GaussPointArray defining the quadrature
    @param Jinv Inverse Jacobian Matrices to compute the gradients
    @param vars std::vector of space Variable s
    @param ret Matrix in which to store the results
    */
    void ieval(const GDtMtGD &f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars, Matrix & ret) ;

    void ieval(const GDtMLtGD &f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars, Matrix & ret) ;

    /** \brief Overloaded function to compute the integral of a GradientDot times a FunctionMatrix times a GradientDot using the inverse Jacobian matrices given by Jinv and the Gauss points in gp, with variables defined by vars.
     * GradientDot is the operator \f$ \dot{\nabla}\otimes \f$ operator.
    @param f GDtMtGD to integrate.
    @param gp GaussPointArray defining the quadrature
    @param Jinv Inverse Jacobian Matrices to compute the gradients
    @param vars std::vector of space Variable s
    */
    Matrix ieval(const GDtMtGD &f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars) ;

    /** \brief Overloaded function to compute the integral of a GradientDot times a Matrix times a Gradient using the inverse Jacobian matrices given by Jinv and the Gauss points in gp, with variables defined by vars.
     * The result is stored in ret. The version of the function can be used to minimise initialisation of memory.
     * GradientDot is the operator \f$ \dot{\nabla}\otimes(\cdot) \f$ operator.
    @param f GDtMtG to integrate.
    @param gp GaussPointArray defining the quadrature
    @param Jinv Inverse Jacobian Matrices to compute the gradients
    @param vars std::vector of space Variable s
    @param ret Matrix in which to store the results
    */
    void ieval(const GDtMtG &f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars, Matrix & ret) ;

    void ieval(const GDDtMtG &f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars, Matrix & ret) ;

    void ieval(const GtMtGDD &f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars, Matrix & ret) ;

    void ieval(const GDtMLtG &f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars, Matrix & ret) ;

    void ieval(const GDDtMLtG &f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars, Matrix & ret) ;

    void ieval(const GtMLtGDD &f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars, Matrix & ret) ;

    /** \brief Overloaded function to compute the integral of a GradientDot times a Matrix times a Gradient using the inverse Jacobian matrices given by Jinv and the Gauss points in gp, with variables defined by vars.
     * GradientDot is the operator \f$ \dot{\nabla}\otimes \f$ operator.
    @param f GDtMtG to integrate.
    @param gp GaussPointArray defining the quadrature
    @param Jinv Inverse Jacobian Matrices to compute the gradients
    @param vars std::vector of space Variable s
    */
    Matrix ieval(const GDtMtG &f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars) ;

    Matrix ieval(const GDtMLtG &f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars) ;

    /** \brief Overloaded function to compute the integral of a Gradient times a Matrix times a GradientDot using the inverse Jacobian matrices given by Jinv and the Gauss points in gp, with variables defined by vars.
     * The result is stored in ret. The version of the function can be used to minimise initialisation of memory.
     * GradientDot is the operator \f$ \dot{\nabla}\otimes(\cdot) \f$ operator.
    @param f GtMtGD to integrate.
    @param gp GaussPointArray defining the quadrature
    @param Jinv Inverse Jacobian Matrices to compute the gradients
    @param vars std::vector of space Variable s
    @param ret Matrix in which to store the results
    */
    void ieval(const GtMtGD &f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars, Matrix & ret) ;

    void ieval(const GtMLtGD &f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars, Matrix & ret) ;

    /** \brief Overloaded function to compute the integral of a Gradient times a Matrix times a GradientDot using the inverse Jacobian matrices given by Jinv and the Gauss points in gp, with variables defined by vars.
     * GradientDot is the operator \f$ \dot{\nabla}\otimes \f$ operator.
    @param f GtMtGD to integrate.
    @param gp GaussPointArray defining the quadrature
    @param Jinv Inverse Jacobian Matrices to compute the gradients
    @param vars std::vector of space Variable s
    */
    Matrix ieval(const GtMtGD &f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars) ;

    /** \brief Overloaded function to compute the integral of a VectorGradient times a FunctionMatrix using the inverse Jacobian matrices given by Jinv and the Gauss points in gp, with variables defined by vars.
     The function is assumed to be expressed in the local coordinates of e.
    	VectorGradient is the usual \f$ \nabla \cdot \f$ operator.
    @param f VGtM to integrate.
    @param e IntegrableEntity on which to integrate
    @param vars std::vector of space Variable s
    */
    Matrix ieval(const VGtM &f, IntegrableEntity *e, const std::vector<Variable> & vars) ;

    /** \brief Overloaded function to compute the integral of a Differential times a Gradient times a Matrix times a Gradient over the IntegrableEntity e, with variables defined by vars.
    The function is assumed to be expressed in the local coordinates of e.
    @param f DtGtMtG to integrate.
    @param e IntegrableEntity on which to integrate
    @param var std::vector of space Variable s
    */
    Matrix ieval(const DtGtMtG & d, IntegrableEntity *e, const std::vector<Variable> & var) ;

    /** \brief Overloaded function to compute the integral of a Differential times a Gradient times a Matrix times a Gradient using the inverse Jacobian matrices given by Jinv and the Gauss points in gp, with variables defined by vars.
    Gradient is the usual \f$ \nabla\otimes \f$ operator, and Differential is the derivative of a function.
    @param f DtGtMtG to integrate.
    @param gp GaussPointArray defining the quadrature
    @param Jinv Inverse Jacobian Matrices to compute the gradients
    @param vars std::vector of space Variable s
    */
    Matrix ieval(const DtGtMtG & d, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars) ;

    void ieval(const DdGtMtG & d, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const IntegrableEntity * e, const std::vector<Variable> & vars, Matrix & ret) ;

    void ieval(const DdGtMtGD & d, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const IntegrableEntity * e, const std::vector<Variable> & vars, Matrix & ret) ;

    void ieval(const DdGtMtG & d, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars, Matrix & ret) ;

    void ieval(const DdGtMtGD & d, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars, Matrix & ret) ;

    void ieval(const DdGtMLtG & d, const std::vector<Matrix> & dmat, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars, Matrix & ret) ;

    void ieval(const DdGtMLtGD & d, const std::vector<Matrix> & dmat, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars, Matrix & ret) ;

    /** \brief Overloaded function to compute the integral of a VectorGradient times a Matrix times a VectorGradient over the IntegrableEntity e, with variables defined by vars.
    The function is assumed to be expressed in the local coordinates of e.
    VectorGradient is the usual \f$ \nabla\ \cdot \f$ operator.
    @param f VGtMtVG to integrate.
    @param e IntegrableEntity on which to integrate
    @param vars std::vector of space Variable s
    */
    double ieval(const VGtMtVG &f, IntegrableEntity *e, const std::vector<Variable> & vars) ;

    /** \brief Overloaded function to compute the line integral of a function, on the seg of Segment supported by the corresponding IntegrableEntity defined by Gamma.
    @param f Function to integrate.
    @param gamma std::vector of pairs of Segment *, IntegrableEntity * defining a support on which to integrate.
    */
    double ieval(const Function &f, const std::vector<std::pair<Segment *, IntegrableEntity *> > & gamma) ;

    /** \brief Overloaded function to compute the integral of a VectorGradient times a Matrix times a VectorGradient using the inverse Jacobian matrices given by Jinv and the Gauss points in gp, with variables defined by vars.
    VectorGradient is the usual \f$ \nabla \cdot \f$ operator.
    @param f VGtMtVG to integrate.
    @param gp GaussPointArray defining the quadrature
    @param Jinv Inverse Jacobian Matrices to compute the gradients
    @param vars std::vector of space Variable s
    */
    double ieval(const VGtMtVG &f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars) ;

    double ieval(const VGDtMtVG &f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars) ;

    double ieval(const VGtMtVGD &f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars) ;

    /** \brief Overloaded function to compute the integral of a VectorGradient times a Matrix times a Vector using the inverse Jacobian matrices given by Jinv and the Gauss points in gp, with variables defined by vars.
    VectorGradient is the usual \f$ \nabla \cdot \f$ operator.
    @param f VGtV to integrate.
    @param gp GaussPointArray defining the quadrature
    @param Jinv Inverse Jacobian Matrices to compute the gradients
    @param vars std::vector of space Variable s
    */
    double ieval(const VGtV &f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars) ;

    /** \brief Overloaded function to compute the integral of a Differential using the inverse Jacobian matrices given by Jinv and the Gauss points in gp, with variables defined by vars.
    Differential is the usual \f$ d \cdot \f$ operator.
    @param d Differential to integrate.
    @param gp GaussPointArray defining the quadrature
    @param Jinv Inverse Jacobian Matrices to compute the gradients
    @param vars std::vector of space Variable s
    */
    double ieval(const Differential & d, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars) ;

    /** \brief Overloaded function to compute the integral of a Differential times a Function using the inverse Jacobian matrices given by Jinv and the Gauss points in gp, with variables defined by vars.
    Differential is the usual \f$ d \cdot \f$ operator.
    @param d DtF to integrate.
    @param gp GaussPointArray defining the quadrature
    @param Jinv Inverse Jacobian Matrices to compute the gradients
    @param vars std::vector of space Variable s
    */
    double ieval(const DtF & d, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars) ;

    /** \brief Overloaded function to compute the integral of a Differential times a Differential using the inverse Jacobian matrices given by Jinv and the Gauss points in gp, with variables defined by vars.
    Differential is the usual \f$ d \cdot \f$ operator.
    @param d DtD to integrate.
    @param gp GaussPointArray defining the quadrature
    @param Jinv Inverse Jacobian Matrices to compute the gradients
    @param vars std::vector of space Variable s
    */
    double ieval(const DtD & d, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars) ;

    double ieval(const DDtF & d, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars) ;

    /** \brief Overloaded function to compute the integral of a Differential over the IntegrableEntity e, with variables defined by vars.
    Differential is the usual \f$ d \cdot \f$ operator.
    @param d Differential to integrate.
    @param e IntegrableEntity on which to integrate
    @param var std::vector of space Variable s
    */
    double ieval(const Differential & d, IntegrableEntity *e, const std::vector<Variable> & var) ;

    /** \brief Overloaded function to compute the integral of a Differential time a Function over the IntegrableEntity e, with variables defined by vars.
    The function is assumed to be expressed in the local coordinates of e.
    Differential is the usual \f$ d \cdot \f$ operator.
    @param d DtF to integrate.
    @param e IntegrableEntity on which to integrate
    @param var std::vector of space Variable s
    */
    double ieval(const DtF & d, IntegrableEntity *e, const std::vector<Variable> & var) ;

    /** \brief Overloaded function to compute the integral of a Differential time a Differential over the IntegrableEntity e, with variables defined by vars.
    The function is assumed to be expressed in the local coordinates of e.
    Differential is the usual \f$ d \cdot \f$ operator.
    @param d DtD to integrate.
    @param e IntegrableEntity on which to integrate
    @param var std::vector of space Variable s
    */
    double ieval(const DtD & d, IntegrableEntity *e, const std::vector<Variable> & var) ;

    /** \brief Overloaded function to compute the value of the gradient of the Function f in the IntegrableEntity e, with variables defined by vars, at point x, y, z, t.
    The function is assumed to be expressed in the local coordinates of e.
    @param f Function the gradient of which should be computed.
    @param e IntegrableEntity on which to perform the computation
    @param vars std::vector of space Variable s
    @param x coordinate at which to evaluate.
    @param y coordinate at which to evaluate.
    @param z coordinate at which to evaluate.
    @param t coordinate at which to evaluate.
    @param transpose transpose the result if true.
    */
    Matrix geval(const Function &f, IntegrableEntity *e, const std::vector<Variable> & vars, const double x, const double y = 0, const double z = 0, const double t = 0, bool transpose = false) ;

    /** \brief Overloaded function to compute the value of the Gradient f in the IntegrableEntity e, with variables defined by vars, at point x, y, z, t.
    The function is assumed to be expressed in the local coordinates of e.
    Gradient is the usual \f$ \nabla\otimes \f$ operator.
    @param f Gradient to evaluate.
    @param e IntegrableEntity on which to perform the computation
    @param vars std::vector of space Variable s
    @param x coordinate at which to evaluate.
    @param y coordinate at which to evaluate.
    @param z coordinate at which to evaluate.
    @param t coordinate at which to evaluate.
    @param transpose transpose the result if true.
    */
    Matrix geval(const Gradient &f, IntegrableEntity *e, const std::vector<Variable> & vars, const double x, const double y = 0, const double z = 0, const double t = 0) ;

    /** \brief Overloaded function to compute the value of the gradient of the Function f using the Matrix m as an inverse Jacobian, with variables defined by vars, at point x, y, z, t.
    Gradient is the usual \f$ \nabla\otimes \f$ operator.
    @param f Function the gradient of which should be computed.
    @param m Inverse Jacobian to use for the computation
    @param vars std::vector of space Variable s
    @param x coordinate at which to evaluate.
    @param y coordinate at which to evaluate.
    @param z coordinate at which to evaluate.
    @param t coordinate at which to evaluate.
    @param transpose transpose the result if true.
    */
    Matrix geval(const Function &f, const Matrix & m, const std::vector<Variable> & vars, const double x, const double y = 0, const double z = 0, const double t = 0, bool transpose = false) ;

    /** \brief Overloaded function to compute the value of the gradient of the Function f using the Matrix m as an inverse Jacobian, with variables defined by vars, at point x, y, z, t.
    The result is stored in ret. The version of the function can be used to minimise initialisation of memory.
    Gradient is the usual \f$ \nabla\otimes \f$ operator.
    @param f Function the gradient of which should be computed.
    @param m Inverse Jacobian to use for the computation
    @param vars std::vector of space Variable s
    @param x coordinate at which to evaluate.
    @param y coordinate at which to evaluate.
    @param z coordinate at which to evaluate.
    @param t coordinate at which to evaluate.
    @param transpose transpose the result if true.
    @param ret Matrix in which to store the result
    */
    void geval(const Function &f, const Matrix & m, const std::vector<Variable> & vars, const double x, const double y , const double z , const double t , bool transpose, Matrix & ret) ;

    /** \brief Overloaded function to compute the value of the gradient of the Function f using the Matrix m as an inverse Jacobian, with variables defined by vars, at point p.getX(), p.getY(), p.getZ(), p.getT().
    The result is stored in ret. The version of the function can be used to minimise initialisation of memory.
    Gradient is the usual \f$ \nabla\otimes \f$ operator.
    @param f Function the gradient of which should be computed.
    @param m Inverse Jacobian to use for the computation
    @param vars std::vector of space Variable s
    @param p Point at which to evaluate
    @param transpose transpose the result if true.
    @param ret Matrix in which to store the result
    */
    void geval(const Function &f, const Matrix & m, const std::vector<Variable> & vars, const Point &p , bool transpose, Matrix & ret) ;

    /** \brief Overloaded function to compute the value of the gradient of the Function f using the Matrix m as an inverse Jacobian, with variables defined by vars, at point p.getX(), p.getY(), p.getZ(), p.getT().
    Gradient is the usual \f$ \nabla\otimes \f$ operator.
    @param f Function the gradient of which should be computed.
    @param m Inverse Jacobian to use for the computation
    @param vars std::vector of space Variable s
    @param p Point at which to evaluate
    @param transpose transpose the result if true.
    */
    Matrix geval(const Function &f, const Matrix & m, const std::vector<Variable> & vars, const Point& p, bool transpose ) ;

    /** \brief Overloaded function to compute the value of the Gradient f using the Matrix m as an inverse Jacobian, with variables defined by vars, at point x, y, z, t.
    Gradient is the usual \f$ \nabla\otimes \f$ operator.
    @param f Gradient to compute.
    @param m Inverse Jacobian to use for the computation
    @param vars std::vector of space Variable s
    @param p Point at which to evaluate
    */
    Matrix geval(const Gradient &f, const Matrix & m, const std::vector<Variable> & vars, const double x, const double y = 0, const double z = 0, const double t = 0) ;

    /** \brief Overloaded function to compute the value of the vector gradient of the function f in the IntegrableEntity e, with variables defined by vars, at point x, y, z, t.
    The function is assumed to be expressed in the local coordinates of e.
    VectorGradient is the usual \f$ \nabla \cdot \f$ operator.
    @param f Function the vector gradient of which should be computed.
    @param e IntegrableEntity on which to perform the computation
    @param vars std::vector of space Variable s
    @param x coordinate at which to evaluate.
    @param y coordinate at which to evaluate.
    @param z coordinate at which to evaluate.
    @param t coordinate at which to evaluate.
    @param transpose transpose the result if true.
    */
    Matrix gveval(const Function &f, IntegrableEntity *e, const std::vector<Variable> & vars, const double x, const double y = 0, const double z = 0, const double t = 0, bool transpose = false) ;

    /** \brief Overloaded function to compute the value of the VectorGradient f in the IntegrableEntity e, with variables defined by vars, at point x, y, z, t.
    The function is assumed to be expressed in the local coordinates of e.
    VectorGradient is the usual \f$ \nabla \cdot \f$ operator.
    @param f VectorGradient which should be computed.
    @param e IntegrableEntity on which to perform the computation
    @param vars std::vector of space Variable s
    @param x coordinate at which to evaluate.
    @param y coordinate at which to evaluate.
    @param z coordinate at which to evaluate.
    @param t coordinate at which to evaluate.
    */
    Matrix gveval(const VectorGradient &f, IntegrableEntity *e, const std::vector<Variable> & vars, const double x, const double y = 0, const double z = 0, const double t = 0) ;

    /** \brief Overloaded function to compute the value of the vector gradient of Function f using the Matrix m as the inverse Jacobian, with variables defined by vars, at point x, y, z, t.
    VectorGradient is the usual \f$ \nabla \cdot \f$ operator.
    @param f Function the vector gradient of which should be computed.
    @param m Inverse Jacobian to use for the computation
    @param vars std::vector of space Variable s
    @param x coordinate at which to evaluate.
    @param y coordinate at which to evaluate.
    @param z coordinate at which to evaluate.
    @param t coordinate at which to evaluate.
    @param transpose transpose the result if true.
    */
    Matrix gveval(const Function &f, const Matrix & m, const std::vector<Variable> & vars, const double x, const double y = 0, const double z = 0, const double t = 0, bool transpose = false) ;

    Matrix gvdeval(const Function &f, const Matrix & m, const std::vector<Variable> & vars, const double x, const double y = 0, const double z = 0, const double t = 0, bool transpose = false) ;

    /** \brief Overloaded function to compute the value of the VectorGradient f using the Matrix m as the inverse Jacobian, with variables defined by vars, at point x, y, z, t.
    VectorGradient is the usual \f$ \nabla \cdot \f$ operator.
    @param f VectorGradient which should be computed.
    @param m Inverse Jacobian to use for the computation
    @param vars std::vector of space Variable s
    @param x coordinate at which to evaluate.
    @param y coordinate at which to evaluate.
    @param z coordinate at which to evaluate.
    @param t coordinate at which to evaluate.
    */
    Matrix gveval(const VectorGradient &f, const Matrix & m, const std::vector<Variable> & vars, const double x, const double y = 0, const double z = 0, const double t = 0) ;

    Matrix gvdeval(const VectorGradientDot &f, const Matrix & m, const std::vector<Variable> & vars, const double x, const double y = 0, const double z = 0, const double t = 0) ;

    /** \brief Overloaded function to compute the value of the GradientDot f using the Matrix m as the inverse Jacobian, with variables defined by vars, at point x, y, z, t.
    GradientDot is the usual \f$ \dot{\nabla}\otimes(\cdot) \f$ operator.
    @param f the function of which the gradient-dot which should be computed.
    @param m Inverse Jacobian to use for the computation
    @param var std::vector of space Variable s
    @param x coordinate at which to evaluate.
    @param y coordinate at which to evaluate.
    @param z coordinate at which to evaluate.
    @param t coordinate at which to evaluate.
    @param transpose transpose the result if true.
    */
    Matrix gdeval(const Function &f, const Matrix & m, const std::vector<Variable> & var, const double x, const double y, const double  z, const double  t, bool transpose ) ;

    Matrix gdeval(const Function &f, IntegrableEntity *e, const std::vector<Variable> & vars, const double x, const double y = 0, const double z = 0, const double t = 0, bool transpose = false) ;

    void gdeval(const Function &f, const Matrix & m, const std::vector<Variable> & var, const double x, const double y, const double  z, const double  t, bool transpose, Matrix & ret ) ;

    Matrix gdeval(const GradientDot &f, IntegrableEntity *e, const std::vector<Variable> & vars, const double x, const double y = 0, const double z = 0, const double t = 0) ;

    /** \brief Overloaded function to compute the values of the GradientDot f using the array of Matrix m as the inverse Jacobian, with variables defined by vars, at points defined by the GaussPointArray gp.
    GradientDot is the usual \f$ \dot{\nabla}\otimes(\cdot) \f$ operator.
    @param f the function of which the gradient-dot which should be computed.
    @param m array of inverse Jacobians to use for the computation
    @param var std::vector of space Variable s
    @param gp GaussPointArray defining the coordinates at which to compute the operator
    @param transpose transpose the result if true.
    */
    std::valarray<Matrix> gdeval(const Function &f, const std::valarray<Matrix> & m, const std::vector<Variable> & var, const GaussPointArray & gp, bool transpose );

    std::valarray<Matrix> gddeval(const Function &f, const std::valarray<Matrix> & m, const std::vector<Variable> & var, const GaussPointArray & gp, bool transpose );

    /** \brief Overloaded function to compute the values of the gradients of f using array of Matrix m as the inverse Jacobians, with variables defined by vars, at points defined by the GaussPointArray gp.
    Gradient is the usual \f$ \nabla\otimes(\cdot) \f$ operator.
    @param f the function of which the gradient which should be computed.
    @param m array of inverse Jacobians to use for the computation
    @param var std::vector of space Variable s
    @param gp GaussPointArray defining the coordinates at which to compute the operator
    @param transpose transpose the result if true.
    */
    std::valarray<Matrix> geval(const Function &f, const std::valarray<Matrix> & m, const std::vector<Variable> & var, const GaussPointArray & gp, bool transpose ) ;

    /** \brief Overloaded function to compute the values of the gradients of f using array of Matrix m as the inverse Jacobians, with variables defined by vars, at points defined by the GaussPointArray gp.
    The result is stored in ret. The version of the function can be used to minimise initialisation of memory.
    Gradient is the usual \f$ \dot{\nabla}\otimes(\cdot) \f$ operator.
    @param f the function of which the gradient which should be computed.
    @param m array of inverse Jacobians to use for the computation
    @param var std::vector of space Variable s
    @param gp GaussPointArray defining the coordinates at which to compute the operator
    @param transpose transpose the result if true.
    @param ret Matrix in which to store the result
    */
    void geval(const Function &f, const std::valarray<Matrix> & m, const std::vector<Variable> & var, const GaussPointArray & gp, bool transpose,  std::valarray<Matrix> & ret) ;

    /** \brief Return all points to be used for sub-tesselation of the IntegrableEntity e, given by two functions
    @param f0 First source of points
    @param f1 Second source of points
    @param e defining surface
    */
    std::vector<Point> allHints(const Function &f0, const Function &f1,IntegrableEntity *e ) ;

    /** \brief Return all points to be used for sub-tesselation of the IntegrableEntity e, given by a functions
    @param f source of points
    @param e defining surface
    */
    std::vector<Point> allHints(const Function &f,IntegrableEntity *e ) ;

    /** \brief Print the ByteCode of the Function f.
    This is sometimes useful for debug, if functions have been generated on-the-fly.
    @param f function to print
    */
    void print(const Function &f) const ;

    /** \brief Print the ByteCode f.
    @param f bytecode to print
    */
} ;


}


#endif
