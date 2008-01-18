// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef VM_BASE_H
#define VM_BASE_H

#include "../elements/integrable_entity.h"
#include "vm_function_base.h"
#include "vm_function_matrix.h"


namespace Mu
{

struct GtFMtG ;
struct VGtMtVG ;
struct GtM ;
struct GtMtG ;
struct GtV ;
struct VGtM ;
struct VGtV ;
struct FMtMtFM ;
struct DtGtMtG ;
struct GDtMtGD ;
struct IntegrableEntity ;
struct GaussPointArray ;
struct Function ;
struct FunctionMatrix ;
struct Differential ;
struct DtF ;
struct Gradient ;
struct VectorGradient ;

class VirtualMachine
{
	Memory stack ;
public:
	VirtualMachine() ;
	
	double eval(const Function &f, const double x, const double y = 0, const double z = 0,const double t=0, const double u = 0, const double v = 0, const double w = 0)   ;
	double eval(const std::vector<RefCountedToken> &f, const double x, const double y = 0, const double z = 0,const double t=0, const double u = 0, const double v = 0, const double w = 0)   ;
	Matrix eval(const FunctionMatrix &f, const double x, const double y = 0, const double z = 0,const double t=0, const double u = 0, const double v = 0, const double w = 0)  ;
	
	double eval(const Function &f, const Point & p, const Point &p_ = Point())  ;
	double eval(const Function &f, const Point * p, const Point *p_ = NULL)  ;
	Matrix eval(const FunctionMatrix &f, const Point & p, const Point &p_ = Point())  ;
	Matrix eval(const FunctionMatrix &f, const Point * p, const Point *p_ = NULL)  ;
	Vector eval(const Function &f, const GaussPointArray &gp) ;
	
	Vector oeval(const Function &f, const GaussPointArray &gp, const Point & offset) ;
	
	double deval(const Function &f, const Variable v,  const double x, const double y = 0, const double z = 0,const double t=0, const double u = 0, const double v = 0, const double w = 0, const double eps = default_derivation_delta)  ;
	Matrix deval(const FunctionMatrix &f, const Variable v,  const double x, const double y = 0, const double z = 0,const double t=0, const double u = 0, const double v = 0, const double w = 0, const double eps= default_derivation_delta)  ;
	double deval(const Function &f, const Variable v, const Point p, const Point p_ = Point(), const double eps= default_derivation_delta)  ;
	Matrix deval(const FunctionMatrix &f, const Variable v, const Point p, const Point p_ = Point(), const double eps= default_derivation_delta)  ;
	
	double deval(const Function &f, const Point&p,  const double x, const double y = 0, const double z = 0, const double t = 0,const double eps = default_derivation_delta, bool normed = false)  ;
	Matrix deval(const FunctionMatrix &f,const Point&p,  const double x, const double y = 0, const double z = 0, const double t = 0,const double eps= default_derivation_delta, bool normed = false)  ;
	double deval(const Function &f,  const Point&p, const Point p, const double eps= default_derivation_delta, bool normed = false)  ;
	Matrix deval(const FunctionMatrix &f,  const Point&p, const Point p, const double eps = default_derivation_delta, bool normed = false)  ;

	double ddeval(const Function &f, const Variable v_0, const Variable v_1,  const double x=0, const double y=0 , const double z=0, const double t=0, const double u=0, const double v=0, const double w=0, const double eps = default_derivation_delta) ;
	Vector ddeval(const Function&f, const Variable v_0, const Variable v_1, const GaussPointArray & gp, const double eps= default_derivation_delta) ;

	Vector deval(const Function&f, const Variable v_, const GaussPointArray & gp, const double eps= default_derivation_delta) ;
	double ieval(const Function &f, const IntegrableEntity *e)  ;
	double ieval(const Function &f, const GaussPointArray &gp)  ;
	double ieval(Vector &f, const IntegrableEntity *e) ;
// 	
	Matrix ieval(const FunctionMatrix &f, const IntegrableEntity *e)  ;
	Matrix ieval(const FMtMtFM &f, const IntegrableEntity *e) ;
	Matrix ieval(const GtM &f, const IntegrableEntity *e, const std::vector<Variable> & vars) ;
	Matrix ieval(const GtMtG &f, const IntegrableEntity *e, const std::vector<Variable> & vars) ;
	Matrix ieval(const GtMtG &f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars) ;
	Vector ieval(const GtV &f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars) ;
	Matrix ieval(const GtFMtG &f, const IntegrableEntity *e, const std::vector<Variable> & vars) ;
	Matrix ieval(const GtFMtG &f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars) ;
	Matrix ieval(const GDtMtGD &f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars) ;
	Matrix ieval(const VGtM &f, const IntegrableEntity *e, const std::vector<Variable> & vars) ;
	Matrix ieval(const DtGtMtG & d, const IntegrableEntity *e, const std::vector<Variable> & var) ;
	Matrix ieval(const DtGtMtG & d, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars) ;
	double ieval(const VGtMtVG &f, const IntegrableEntity *e, const std::vector<Variable> & vars) ;
	double ieval(const Function &f, const std::vector<std::pair<Segment *, IntegrableEntity *> > & gamma) ;
	double ieval(const VGtMtVG &f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars) ;
	double ieval(const VGtV &f, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars) ;
	double ieval(const Differential & d, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars) ;
	double ieval(const DtF & d, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars) ;
	double ieval(const Differential & d, const IntegrableEntity *e, const std::vector<Variable> & var) ;
	double ieval(const DtF & d, const IntegrableEntity *e, const std::vector<Variable> & var) ;

	Matrix geval(const Function &f, const IntegrableEntity *e, const std::vector<Variable> & vars, const double x, const double y = 0, const double z = 0, const double t = 0, bool transpose = false) ;	
	Matrix geval(const Gradient &f, const IntegrableEntity *e, const std::vector<Variable> & vars, const double x, const double y = 0, const double z = 0, const double t = 0) ;
	Matrix geval(const Function &f, const Matrix & m, const std::vector<Variable> & vars, const double x, const double y = 0, const double z = 0, const double t = 0, bool transpose = false) ;	
	Matrix geval(const Function &f, const Matrix & m, const std::vector<Variable> & vars, const Point& p, bool transpose ) ;
	Matrix geval(const Gradient &f, const Matrix & m, const std::vector<Variable> & vars, const double x, const double y = 0, const double z = 0, const double t = 0) ;
	Matrix gveval(const Function &f, const IntegrableEntity *e, const std::vector<Variable> & vars, const double x, const double y = 0, const double z = 0, const double t = 0, bool transpose = false) ;	
	Matrix gveval(const VectorGradient &f, const IntegrableEntity *e, const std::vector<Variable> & vars, const double x, const double y = 0, const double z = 0, const double t = 0) ;
	Matrix gveval(const Function &f, const Matrix & m, const std::vector<Variable> & vars, const double x, const double y = 0, const double z = 0, const double t = 0, bool transpose = false) ;	
	Matrix gveval(const VectorGradient &f, const Matrix & m, const std::vector<Variable> & vars, const double x, const double y = 0, const double z = 0, const double t = 0) ;

	Matrix gdeval(const Function &f, const Matrix & m, const std::vector<Variable> & var, const double x, const double y, const double  z, const double  t, bool transpose ) ;
	std::valarray<Matrix> gdeval(const Function &f, const std::valarray<Matrix> & m, const std::vector<Variable> & var, const GaussPointArray & gp, bool transpose );

	std::valarray<Matrix> geval(const Function &f, const std::valarray<Matrix> & m, const std::vector<Variable> & var, const GaussPointArray & gp, bool transpose ) ;

	std::vector<Point> allHints(const Function &, const Function &,const IntegrableEntity *e ) ;
	std::vector<Point> allHints(const Function &,const IntegrableEntity *e ) ;
	
	void print(const Function &f) const ;
	void print(const std::vector<RefCountedToken> &f) const ;
} ;


} ;


#endif
