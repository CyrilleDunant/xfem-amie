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

const double default_derivation_delta= 2e-8 ;

struct GtFMtG ;
struct IntegrableEntity ;

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
	
	double deval(const Function &f, const Variable v,  const double x, const double y = 0, const double z = 0,const double t=0, const double u = 0, const double v = 0, const double w = 0, const double eps = default_derivation_delta)  ;
	Matrix deval(const FunctionMatrix &f, const Variable v,  const double x, const double y = 0, const double z = 0,const double t=0, const double u = 0, const double v = 0, const double w = 0, const double eps= default_derivation_delta)  ;
	double deval(const Function &f, const Variable v, const Point p, const Point p_ = Point(), const double eps= default_derivation_delta)  ;
	Matrix deval(const FunctionMatrix &f, const Variable v, const Point p, const Point p_ = Point(), const double eps= default_derivation_delta)  ;
	
	double deval(const Function &f, const Point&p,  const double x, const double y = 0, const double z = 0, const double t = 0,const double eps = default_derivation_delta, bool normed = false)  ;
	Matrix deval(const FunctionMatrix &f,const Point&p,  const double x, const double y = 0, const double z = 0, const double t = 0,const double eps= default_derivation_delta, bool normed = false)  ;
	double deval(const Function &f,  const Point&p, const Point p, const double eps= default_derivation_delta, bool normed = false)  ;
	Matrix deval(const FunctionMatrix &f,  const Point&p, const Point p, const double eps= default_derivation_delta, bool normed = false)  ;
	
	double ieval(const Function &f, const IntegrableEntity *e)  ;
	double ieval(const Function &f, const std::valarray< std::pair<Point, double> > &gp)  ;
	double ieval(Vector &f, const IntegrableEntity *e) ;
// 	
	Matrix ieval(const FunctionMatrix &f, const IntegrableEntity *e)  ;
	Matrix ieval(const FMtMtFM &f, const IntegrableEntity *e) ;
	Matrix ieval(const GtM &f, const IntegrableEntity *e, const std::vector<Variable> & vars) ;
	Matrix ieval(const GtMtG &f, const IntegrableEntity *e, const std::vector<Variable> & vars) ;
	Matrix ieval(const GtMtG &f, const std::valarray< std::pair<Point, double> > &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars) ;
	Vector ieval(const GtV &f, const std::valarray< std::pair<Point, double> > &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars) ;
	Matrix ieval(const GtFMtG &f, const IntegrableEntity *e, const std::vector<Variable> & vars) ;
	Matrix ieval(const GtFMtG &f, const std::valarray< std::pair<Point, double> > &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars) ;
	Matrix ieval(const VGtM &f, const IntegrableEntity *e, const std::vector<Variable> & vars) ;
	double ieval(const VGtMtVG &f, const IntegrableEntity *e, const std::vector<Variable> & vars) ;
	double ieval(const Function &f, const std::vector<std::pair<Segment *, IntegrableEntity *> > & gamma) ;
	double ieval(const VGtMtVG &f, const std::valarray< std::pair<Point, double> > &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars) ;
	double ieval(const VGtV &f, const std::valarray< std::pair<Point, double> > &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars) ;
	double ieval(const Differential & d, const std::valarray< std::pair<Point, double> > &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars) ;
	double ieval(const DtF & d, const std::valarray< std::pair<Point, double> > &gp, const std::valarray<Matrix> &Jinv, const std::vector<Variable> & vars) ;
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

	std::vector<Point> allHints(const Function &, const Function &,const IntegrableEntity *e ) ;
	std::vector<Point> allHints(const Function &,const IntegrableEntity *e ) ;
	
	void print(const Function &f) const ;
	void print(const std::vector<RefCountedToken> &f) const ;
} ;

struct GtFM
{
	const Gradient first ;
	const FunctionMatrix second ;
	
	GtFM(const Gradient & g, const FunctionMatrix & f) : first(g), second(f) { };
	GtFMtG operator*(const Mu::Gradient & f) const ;
} ;

struct GtFMtG
{
	const Gradient first ;
	const FunctionMatrix second ;
	const Gradient third ; ;
	
	GtFMtG(const Gradient & g, const FunctionMatrix & f,const Gradient & g_) : first(g), second(f), third(g_) { };
	
} ;

} ;

Mu::GtFM operator *(const Mu::Gradient g, const Mu::FunctionMatrix m) ;

#endif
