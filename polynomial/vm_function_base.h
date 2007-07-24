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
#include "../matrixops.h"
#include "../geometry/geometry_base.h"

namespace Mu
{
typedef enum
{
	POSITION_TOKEN,
	PROJECTION_TOKEN
} PositionTokenType ;

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
	
	int ptID ;
	int dofID ;
	
protected:
	ByteCode byteCode ;
	
	bool e_diff ;
	
	bool isBinaryOperator(const Token * t) const ;

	bool isUnaryOperator(const Token * t) const ;

// 	void vectorizeOne(std::vector<RefCountedToken>   &bytecode,  size_t &lastAddress , int & precalculatedEnd ) const ;
// 	void vectorizeTwo(std::vector<RefCountedToken>   &bytecode,  size_t &lastAddress , int & precalculatedEnd ) const ;
// 	void factorize(std::vector<RefCountedToken>   &bytecode,  size_t &lastAddress , int & precalculatedEnd ) const ;
	
	
public:
	Function();
	Function(const char *f) ;
	Function(const Line & l, Function x, Function y) ;
	Function(const Point & p, Function x, Function y) ;
	Function(const std::string &f) ;
	Function(const std::valarray< std::valarray<Matrix> > & coeffs, bool diff = true) ;
	Function(const std::valarray<Matrix> & coeffs, bool diff = true) ;
	Function(const Matrix & coeffs, bool diff = true) ;
	Function(const std::valarray<double> & coeffs, bool diff = true) ;
	Function(const Segment s, const Function & x, const Function & y, PositionTokenType = POSITION_TOKEN) ;
	Function(const std::vector<Segment> s, const Function & x, const Function & y,PositionTokenType = POSITION_TOKEN) ;
	Function(const ByteCode & b_0, const ByteCode & b_1) ;
	Function(Geometry * g) ;
	Function(Geometry * g, Function x, Function y) ;

public:
	Function(const ByteCode &b_0, const ByteCode &b_1, RefCountedToken op, const bool diff = false); 
	Function(const ByteCode &b_0, const double a, RefCountedToken op, const bool diff = false); 
	Function(const double a, const ByteCode &b_0,  RefCountedToken op, const bool diff= false) ;
	
	virtual ~Function()  ;
	
	bool isNull() const ;
	
	const ByteCode & getByteCode() const ;
	ByteCode & getByteCode() ;
	const RefCountedToken& getToken(const size_t i) const ;
	
	size_t size() const ;
	
	const Function & d(const Variable v) const ;
	Function &d(const Variable v) ;
	std::valarray<Function> getDerivatives() const;
	std::valarray<Function> & getDerivatives();
	bool isDifferentiable() const ;
	
	const std::vector< Point > & getIntegrationHint() const ;
	Point  getIntegrationHint(size_t) const ;
	void setIntegrationHint(const std::vector< Point >) ;
	void addIntegrationHint(const Point ) ;
	bool hasIntegrationHint() const ;
	
	int getDofID() const;
	void setDofID(size_t) ;
	
	int getPointID() const;
	void setPointID(size_t) ;
	
	Function & operator=(const Function &f) ;
	
	Function operator*(const Function &f) const ;
	
	Function operator*(const Geometry *f) const ;
	
	Function operator/(const Function &f) const ;
	
	Function operator+(const Function &f) const ;
	
	Function operator-(const Function &f) const ;
	
	Function operator*(const double a) const ;
	
	Function operator/(const double a) const ;
	
	Function operator+(const double a) const ;
	
	Function operator-(const double a) const ;
	
	Function operator^(const int a) const ;
	
	void operator*=(const Function &f)  ;
	
	void operator*=(const Geometry *f)  ;
	
	void operator/=(const Function &f)  ;
	
	void operator+=(const Function &f)  ;
	
	void operator-=(const Function &f)  ;
	
	void operator*=(const double a)  ;
	
	void operator/=(const double a)  ;
	
	void operator+=(const double a)  ;
	
	void operator-=(const double a)  ;
	
	Function operator()(const Function & f) const ;
	
	void tokenEval(const size_t i, Context & context) const
	{
		byteCode[i].eval(context ) ;
	}
	
	void compile() ;
	
	
	
} ;

struct DtF ;
struct GtM ;
struct GtV ;
struct GtMtG ;
struct VGtM ;
struct VGtV ;
struct VGtMtVG ;

struct Differential
{
	const Function & f ;
	const Variable & v ;
	Differential(const Function &u, const Variable & m) : f(u), v(m) { } ;
	DtF operator *(const Function & f) const ;
} ;

struct Gradient
{
	const Function & f ;
	const bool transpose ;
	Gradient(const Function &u, bool t = false) : f(u), transpose(t) { };
	GtM operator *(const Matrix & f) const ;
	GtV operator *(const Vector & f) const ;
} ;

struct VectorGradient
{
	const Function & f ;
	const bool transpose ;
	VectorGradient(const Function &u, bool t = false) : f(u), transpose(t) { };
	VGtM operator *(const Matrix & f) const ;
	VGtV operator *(const Vector & f) const ;
} ;

struct GtM
{
	const Gradient & first ;
	const Matrix & second ;
	
	GtM(const Gradient & g, const Matrix & f) : first(g), second(f) { };
	GtMtG operator*(const Mu::Gradient & f) const ;
} ;

struct DtF
{
	const Differential & d ;
	const Function & f ;
	DtF(const Differential & d_, const Function & f_) : d(d_), f(f_) { } ;
} ;

struct VGtM
{
	const VectorGradient & first ;
	const Matrix & second ;
	
	VGtM(const VectorGradient & g, const Matrix & f) : first(g), second(f) { };
	VGtMtVG operator*(const Mu::VectorGradient & f) const ;
} ;

struct GtV
{
	const Gradient & first ;
	const Vector & second ;
	
	GtV(const Gradient & g, const Vector & f) : first(g), second(f) { };
} ;

struct VGtV
{
	const VectorGradient & first ;
	const Vector & second ;
	
	VGtV(const VectorGradient & g, const Vector & f) : first(g), second(f) { };
} ;

struct GtMtG
{
	const Gradient & first ;
	const Matrix & second ;
	const Gradient & third ;
	
	GtMtG(const Gradient & g, const Matrix & f,const Gradient & g_) : first(g), second(f), third(g_) { };
	
} ;

struct VGtMtVG
{
	const VectorGradient & first ;
	const Matrix & second ;
	const VectorGradient & third ;
	
	VGtMtVG(const VectorGradient & g, const Matrix & f,const VectorGradient & g_) : first(g), second(f), third(g_) { };
	
} ;


} ;


Mu::Function operator-(const double & a, const Mu::Function &f) ;
Mu::Function operator*(const double & a, const Mu::Function &f) ;
Mu::Function operator+(const double & a, const Mu::Function &f) ;
Mu::Function operator/(const double & a, const Mu::Function &f) ;

Mu::Function f_sqrt(const Mu::Function &f) ;
Mu::Function f_exp(const Mu::Function &f) ;
Mu::Function f_abs(const Mu::Function &f) ;
Mu::Function f_log(const Mu::Function &f) ;
Mu::Function f_atan2(const Mu::Function &f0, const Mu::Function &f1) ;
Mu::Function f_sin(const Mu::Function &f) ;
Mu::Function f_cos(const Mu::Function &f) ;
Mu::Function f_sign(const Mu::Function &f) ;
Mu::Function f_positivity(const Mu::Function &f) ;
Mu::Function f_negativity(const Mu::Function &f) ;


#endif

