//
// C++ Implementation: vm_function_base
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "vm_function_base.h"
#include "../elements/elements.h"
#include <string.h>

using namespace Mu ;


void concatenateFunctions(const Function & src0_, const Function & src1_, Function & dst)
{
	Function src0(src0_) ;
	Function src1(src1_) ;
	dst.byteCodeSize = src0.byteCodeSize+src1.byteCodeSize ;
	dst.adress_a = src0.adress_a ;
	dst.use_temp = src0.use_temp ;
	dst.byteCode = src0.byteCode ;
	
	for(size_t i = 0 ; i < src0.constNumber ; i++)
		dst.values[i] = src0.values[i] ;
	for(size_t i = 0 ; i < src1.constNumber ; i++)
		dst.values[i+src0.constNumber] = src1.values[i] ;
	
	dst.constNumber = src0.constNumber ;
	
	for(size_t i = 0 ; i < src0.byteCodeSize ; i++)
	{
		if(src0.geo_op[i])
			dst.geo_op[i] = src0.geo_op[i]->getCopy() ;
	}
	
	for(size_t i = 0 ; i < src1.byteCodeSize ; i++)
	{
		dst.use_temp[i+src0.byteCodeSize] = src1.use_temp[i] ;
		dst.byteCode[i+src0.byteCodeSize] = src1.byteCode[i] ;
		
		if(src1.geo_op[i])
			dst.geo_op[i+src0.byteCodeSize] = src1.geo_op[i]->getCopy() ;
			
		if(src1.adress_a[i*4] >= FUNCTION_LENGTH-1-src1.constNumber)
		{
			dst.adress_a[(i+src0.byteCodeSize)*4] = src1.adress_a[i*4]-src0.constNumber ;
			dst.constNumber++ ;
		}
		
		if(src1.adress_a[i*4+1] >= FUNCTION_LENGTH-1-src1.constNumber)
		{
			dst.adress_a[(i+src0.byteCodeSize)*4+1] = src1.adress_a[i*4+1]-src0.constNumber ;
			dst.constNumber++ ;
		}
		
		if(src1.adress_a[i*4+2] >= FUNCTION_LENGTH-1-src1.constNumber)
		{
			dst.adress_a[(i+src0.byteCodeSize)*4+2] = src1.adress_a[i*4+2]-src0.constNumber ;
			dst.constNumber++ ;
		}
		
		if(src1.adress_a[i*4] >= 8 && src1.adress_a[i*4] < FUNCTION_LENGTH-1-src1.constNumber)
			dst.adress_a[(i+src0.byteCodeSize)*4] = src1.adress_a[i*4]+1 ;
		if(src1.adress_a[i*4+1] >= 8 && src1.adress_a[i*4+1] < FUNCTION_LENGTH-1-src1.constNumber)
			dst.adress_a[(i+src0.byteCodeSize)*4+1] = src1.adress_a[i*4+1]+1 ;
		if(src1.adress_a[i*4+2] >= 8 && src1.adress_a[i*4+2] < FUNCTION_LENGTH-1-src1.constNumber)
			dst.adress_a[(i+src0.byteCodeSize)*4+2] = src1.adress_a[i*4+2]+1 ;
		
		if(src1.adress_a[i*4] < 8 )
			dst.adress_a[(i+src0.byteCodeSize)*4] = src1.adress_a[i*4] ;
		if(src1.adress_a[i*4+1] < 8 )
			dst.adress_a[(i+src0.byteCodeSize)*4+1] = src1.adress_a[i*4+1] ;
		if(src1.adress_a[i*4+2] < 8 )
			dst.adress_a[(i+src0.byteCodeSize)*4+2] = src1.adress_a[i*4+2] ;
	}
}

void concatenateFunctions(const Function & src0_, const Function & src1_, const Function & src2_, Function & dst)
{
	Function src0(src0_) ;
	Function src1(src1_) ;
	Function src2(src2_) ;
	dst.byteCodeSize = src0.byteCodeSize+src1.byteCodeSize+src2.byteCodeSize ;
	dst.adress_a = src0.adress_a ;
	dst.use_temp = src0.use_temp ;
	dst.byteCode = src0.byteCode ;
	
	for(size_t i = 0 ; i < src0.constNumber ; i++)
		dst.values[i] = src0.values[i] ;
	for(size_t i = 0 ; i < src1.constNumber ; i++)
		dst.values[i+src0.constNumber] = src1.values[i] ;
	for(size_t i = 0 ; i < src2.constNumber ; i++)
		dst.values[i+src0.constNumber+src1.constNumber] = src2.values[i] ;
	
	dst.constNumber = src0.constNumber ;
	
	for(size_t i = 0 ; i < src0.byteCodeSize ; i++)
	{
		if(src0.geo_op[i])
			dst.geo_op[i] = src0.geo_op[i]->getCopy() ;
	}
	
	for(size_t i = 0 ; i < src1.byteCodeSize ; i++)
	{
		dst.use_temp[i+src0.byteCodeSize] = src1.use_temp[i] ;
		dst.byteCode[i+src0.byteCodeSize] = src1.byteCode[i] ;
		
		if(src1.geo_op[i])
			dst.geo_op[i+src0.byteCodeSize] = src1.geo_op[i]->getCopy() ;
			
		if(src1.adress_a[i*4] >= FUNCTION_LENGTH-1-src1.constNumber)
		{
			dst.adress_a[(i+src0.byteCodeSize)*4] = src1.adress_a[i*4]-src0.constNumber ;
			dst.constNumber++ ;
		}
		
		if(src1.adress_a[i*4+1] >= FUNCTION_LENGTH-1-src1.constNumber)
		{
			dst.adress_a[(i+src0.byteCodeSize)*4+1] = src1.adress_a[i*4+1]-src0.constNumber ;
			dst.constNumber++ ;
		}
		
		if(src1.adress_a[i*4+2] >= FUNCTION_LENGTH-1-src1.constNumber)
		{
			dst.adress_a[(i+src0.byteCodeSize)*4+2] = src1.adress_a[i*4+2]-src0.constNumber ;
			dst.constNumber++ ;
		}
		
		if(src1.adress_a[i*4] >= 8 && src1.adress_a[i*4] < FUNCTION_LENGTH-1-src1.constNumber)
			dst.adress_a[(i+src0.byteCodeSize)*4] = src1.adress_a[i*4]+1 ;
		
		if(src1.adress_a[i*4+1] >= 8 && src1.adress_a[i*4+1] < FUNCTION_LENGTH-1-src1.constNumber)
			dst.adress_a[(i+src0.byteCodeSize)*4+1] = src1.adress_a[i*4+1]+1 ;
		
		if(src1.adress_a[i*4+2] >= 8 && src1.adress_a[i*4+2] < FUNCTION_LENGTH-1-src1.constNumber)
			dst.adress_a[(i+src0.byteCodeSize)*4+2] = src1.adress_a[i*4+2]+1 ;
		
		if(src1.adress_a[i*4] < 8 )
			dst.adress_a[(i+src0.byteCodeSize)*4] = src1.adress_a[i*4] ;
		
		if(src1.adress_a[i*4+1] < 8 )
			dst.adress_a[(i+src0.byteCodeSize)*4+1] = src1.adress_a[i*4+1] ;
		
		if(src1.adress_a[i*4+2] < 8 )
			dst.adress_a[(i+src0.byteCodeSize)*4+2] = src1.adress_a[i*4+2] ;
	}
	
	for(size_t i = 0 ; i < src2.byteCodeSize ; i++)
	{
		dst.use_temp[i+src0.byteCodeSize+src1.byteCodeSize] = src2.use_temp[i] ;
		dst.byteCode[i+src0.byteCodeSize+src1.byteCodeSize] = src2.byteCode[i] ;
		
		if(src2.geo_op[i])
			dst.geo_op[i+src0.byteCodeSize+src1.byteCodeSize] = src2.geo_op[i]->getCopy() ;
			
		if(src2.adress_a[i*4] >= FUNCTION_LENGTH-1-src2.constNumber)
		{
			dst.adress_a[(i+src0.byteCodeSize+src1.byteCodeSize)*4] = src2.adress_a[i*4]-src0.constNumber-src1.constNumber ;
			dst.constNumber++ ;
		}
		
		if(src2.adress_a[i*4+1] >= FUNCTION_LENGTH-1-src2.constNumber)
		{
			dst.adress_a[(i+src0.byteCodeSize+src1.byteCodeSize)*4+1] = src2.adress_a[i*4+1]-src0.constNumber-src1.constNumber ;
			dst.constNumber++ ;
		}
		
		if(src2.adress_a[i*4+2] >= FUNCTION_LENGTH-1-src2.constNumber)
		{
			dst.adress_a[(i+src0.byteCodeSize+src1.byteCodeSize)*4+2] = src2.adress_a[i*4+2]-src0.constNumber-src1.constNumber ;
			dst.constNumber++ ;
		}
		
		if(src2.adress_a[i*4] >= 8 && src2.adress_a[i*4] < FUNCTION_LENGTH-1-src2.constNumber)
			dst.adress_a[(i+src0.byteCodeSize+src1.byteCodeSize)*4] = src2.adress_a[i*4]+2 ;
		
		if(src2.adress_a[i*4+1] >= 8 && src2.adress_a[i*4+1] < FUNCTION_LENGTH-1-src2.constNumber)
			dst.adress_a[(i+src0.byteCodeSize+src1.byteCodeSize)*4+1] = src2.adress_a[i*4+1]+2 ;
		
		if(src2.adress_a[i*4+2] >= 8 && src2.adress_a[i*4+2] < FUNCTION_LENGTH-1-src2.constNumber)
			dst.adress_a[(i+src0.byteCodeSize+src1.byteCodeSize)*4+2] = src2.adress_a[i*4+2]+2 ;
		
		if(src2.adress_a[i*4] < 8 )
			dst.adress_a[(i+src0.byteCodeSize+src1.byteCodeSize)*4] = src2.adress_a[i*4] ;
		
		if(src2.adress_a[i*4+1] < 8 )
			dst.adress_a[(i+src0.byteCodeSize+src1.byteCodeSize)*4+1] = src2.adress_a[i*4+1] ;
		
		if(src2.adress_a[i*4+2] < 8 )
			dst.adress_a[(i+src0.byteCodeSize+src1.byteCodeSize)*4+2] = src2.adress_a[i*4+2] ;
	}
}


void  Function::preCalculate(const GaussPointArray & gp , std::vector<Variable> & var, const double eps)
{
	//this needs to be in to stages, otherwise memory gets accessed too early.
	Vector * newVal = new Vector(VirtualMachine().eval(*this, gp)) ;
	if(precalc.find(gp.id) != precalc.end())
		delete precalc[gp.id] ;
	precalc[gp.id] = newVal ;
	std::map<Variable, Vector *> val ;
	for(size_t i = 0 ; i < var.size() ; i++)
	{
		if(val.find(var[i]) == val.end())
			val[var[i]] = new Vector(VirtualMachine().deval(*this,var[i] ,gp, eps)) ;
	}
	if(dprecalc.find(gp.id) != dprecalc.end())
	{
		for(auto i = dprecalc[gp.id].begin() ;  i != dprecalc[gp.id].end() ; ++i)
			delete i->second ;
	}

	dprecalc[gp.id] = val ;
}

void  Function::preCalculate(const GaussPointArray & gp )
{
	if(precalc.find(gp.id) != precalc.end())
		delete precalc[gp.id] ;
	precalc[gp.id] = new Vector(VirtualMachine().eval(*this, gp)) ;
}

const Vector & Function::getPrecalculatedValue(const GaussPointArray &gp ) const
{
	return *(precalc.find(gp.id)->second) ;
}

const Vector & Function::getPrecalculatedValue(const GaussPointArray &gp, Variable v) const
{
	return *dprecalc.find(gp.id)->second.find(v)->second ;
}

bool Function::precalculated(const GaussPointArray & gp) const
{
	return (precalc.find(gp.id) != precalc.end()) ;
}

bool Function::precalculated(const GaussPointArray & gp, Variable v) const
{
	auto precalculatedDerivative = dprecalc.find(gp.id) ;
	return precalculatedDerivative != dprecalc.end() && precalculatedDerivative->second.find(v) != precalculatedDerivative->second.end() ;
}


GtM Gradient::operator*( const Matrix & f) const
{
	return GtM(*this, f) ;
}

GDDtM GradientDotDot::operator*( const Matrix & f) const
{
	return GDDtM(*this, f) ;
}

GDDtML GradientDotDot::operator*( const std::vector<Matrix> & f) const
{
	return GDDtML(*this, f) ;
}

GDDtMtG GDDtM::operator*( const Gradient & f) const
{
	return GDDtMtG(*this, f) ;
}

GDDtMLtG GDDtML::operator*( const Gradient & f) const
{
	return GDDtMLtG(*this, f) ;
}

GtML Gradient::operator*( const std::vector<Matrix> & f) const
{
	return GtML(*this, f) ;
}

GDtML GradientDot::operator*( const std::vector<Matrix> & f) const
{
	return GDtML(*this, f) ;
}

GDtMLtG GDtML::operator*( const Gradient & f) const
{
	return GDtMLtG(*this, f) ;
}

GDtMLtGD GDtML::operator*( const GradientDot & f) const
{
	return GDtMLtGD(*this, f) ;
}


GtV Gradient::operator*( const Vector & f) const
{
	return GtV(*this, f) ;
}

GtVL Gradient::operator*( const std::vector<Vector> & f) const
{
	return GtVL(*this, f) ;
}

GDtVL GradientDot::operator*( const std::vector<Vector> & f) const
{
	return GDtVL(*this, f) ;
}

GtMtG GtM::operator*(const Gradient & f) const
{
	return GtMtG(this->first, this->second, f) ;
}

GtMLtG GtML::operator*(const Gradient & f) const
{
	return GtMLtG(this->first, this->second, f) ;
}

GtMLtGD GtML::operator*(const GradientDot & f) const
{
	return GtMLtGD(*this, f) ;
}

GtMtGD GtM::operator*(const GradientDot & f) const
{
	return GtMtGD(*this, f) ;
}

DtF Differential::operator *(const Function & f) const
{
	return DtF(*this, f) ;
}

DtD Differential::operator *(const Differential & f) const
{
	return DtD(*this, f) ;
}

DtV Differential::operator *(const Vector & f) const
{
	return DtV(*this, f) ;
}

DtVL Differential::operator *(const std::vector<Vector> & f) const
{
	return DtVL(*this, f) ;
}

DtGtMtG Differential::operator *(const GtMtG & g) const
{
	return DtGtMtG(*this, g) ;
}

GDtMtGD GDtM::operator*(const Mu::GradientDot & f) const
{
	return GDtMtGD(*this, f) ;
}

GDtMtG GDtM::operator*(const Mu::Gradient & f) const
{
	return GDtMtG(*this, f) ;
}

GDtM GradientDot::operator *(const Matrix & f) const
{
	return GDtM(*this, f) ;
}

GDtV GradientDot::operator *(const Vector & v) const
{
	return GDtV(*this, v) ;
}

VGtM VectorGradient::operator *(const Matrix & f) const 
{ 
	return VGtM(*this, f) ; 
}

VGtV VectorGradient::operator *(const Vector & f) const 
{ 
	return VGtV(*this, f) ; 
}


VGtMtVG VGtM::operator*(const Mu::VectorGradient & f) const 
{
	return VGtMtVG(first,second, f) ;
}


Function::Function() 
{
	defaultInitialise();
}

	
Function Function::operator*(const Geometry *f) const
{

	Function f_(*this) ;
	f_.byteCode[byteCodeSize] = TOKEN_OPERATION_GEO_OPERATION ;
	f_.geo_op[byteCodeSize] = new DomainOperation(f) ;
	f_.byteCode[byteCodeSize+1] = TOKEN_OPERATION_TIMES ;
	f_.byteCodeSize += 2;
	return f_ ;
}

void Function::operator*=(const Geometry *f) 
{
	byteCode[byteCodeSize] = TOKEN_OPERATION_GEO_OPERATION ;
	geo_op[byteCodeSize] = new DomainOperation(f) ;
	byteCode[byteCodeSize+1] = TOKEN_OPERATION_TIMES ;
	byteCodeSize += 2;
}

Function & Function::operator=(const Function &f)
{
	if(f.derivative && false)
	{
		derivative = new std::valarray<Function *>((Function *)nullptr, f.derivative->size()) ;
		e_diff = true ;
		for(size_t i = 0 ; i < f.derivative->size() ; i++)
			(*derivative)[i] = new Function(*(*f.derivative)[i]) ;
	}
	else
	{
		derivative = nullptr ;
		e_diff = false ;
	}
	
	ptID = f.ptID ;
	dofID = f.dofID ;
	iPoint = f.iPoint ;
	byteCode = f.byteCode ;
	values = f.values ;
	byteCodeSize = f.byteCodeSize ;
	use_temp = f.use_temp ;
	adress_a = f.adress_a ;

	constNumber = f.constNumber ;

	for(size_t i = 0 ; i < f.byteCodeSize ; i++)
	{
		delete geo_op[i] ;
		if(f.geo_op[i])
			geo_op[i] = f.geo_op[i]->getCopy() ;
	}
// 	initialiseAdresses();

	for(auto i = f.precalc.begin() ; i != f.precalc.end() ; ++i)
		precalc[i->first] = new Vector(*i->second) ;
	for(auto i = f.dprecalc.begin() ; i != f.dprecalc.end() ; ++i)
	{
		dprecalc[i->first] =  std::map<Variable, Vector *>() ;
		for(auto j = i->second.begin() ; j != i->second.end() ; ++j)
		{
			dprecalc[i->first][j->first] = new Vector(*j->second) ;
		}
	}
	return *this ;
}

int Function::getDofID() const
{
	return this->dofID ;
}
void Function::setDofID(size_t id)
{
	this->dofID = id ;
}

Point * Function::getPoint() const
{
	return this->ptID ;
}
void Function::setPoint(Point * id)
{
	this->ptID = id ;
}

boost::tuple<TokenOperationType, double, std::string> Function::toToken(const std::string & str) const
{
	if(isNumeral(str[0]) || (str[0] == '-' && isNumeral(str[1])))
	{
		return boost::make_tuple(TOKEN_OPERATION_CONSTANT, atof(str.c_str()), "") ;
	}
	else if(str == std::string("sin"))
	{
		return boost::make_tuple(TOKEN_OPERATION_SIN, 0., "") ;
	}
	else if(str == std::string("cos"))
	{
		return boost::make_tuple(TOKEN_OPERATION_COS, 0., "") ;
	}
	else if(str == std::string("tan"))
	{
		return boost::make_tuple(TOKEN_OPERATION_TAN, 0., "") ;
	}
	else if(str == std::string("sinh"))
	{
		return boost::make_tuple(TOKEN_OPERATION_SINH, 0., "") ;
	}
	else if(str == std::string("cosh"))
	{
		return boost::make_tuple(TOKEN_OPERATION_COSH, 0., "") ;
	}
	else if(str == std::string("tanh"))
	{
		return boost::make_tuple(TOKEN_OPERATION_TANH, 0., "") ;
	}
	else if(str == std::string("exp"))
	{
		return boost::make_tuple(TOKEN_OPERATION_EXP, 0., "") ;
	}
	else if(str == std::string("abs"))
	{
		return boost::make_tuple(TOKEN_OPERATION_ABS, 0., "") ;
	}
	else if(str == std::string("log"))
	{
		return boost::make_tuple(TOKEN_OPERATION_LOG, 0., "") ;
	}
	else if(str == std::string("sqrt"))
	{
		return boost::make_tuple(TOKEN_OPERATION_SQRT, 0., "") ;
	}
	else if(str == std::string("sign"))
	{
		return boost::make_tuple(TOKEN_OPERATION_SIGN, 0., "") ;
	}
	else if(str == std::string("atan2"))
	{
		return boost::make_tuple(TOKEN_OPERATION_ATAN2, 0., "") ;
	}
	else if(str == std::string("sign"))
	{
		return boost::make_tuple(TOKEN_OPERATION_SIGN, 0., "") ;
	}
	else if(str == std::string("positive"))
	{
		return boost::make_tuple(TOKEN_OPERATION_POSITIVITY, 0., "") ;
	}
	else if(str == std::string("negative"))
	{
		return boost::make_tuple(TOKEN_OPERATION_NEGATIVITY, 0., "") ;
	}
	else if(str == std::string("+"))
	{
		return boost::make_tuple(TOKEN_OPERATION_PLUS, 0., "") ;
	}
	else if(str == std::string("-"))
	{
		return boost::make_tuple(TOKEN_OPERATION_MINUS, 0., "") ;
	}
	else if(str == std::string("*"))
	{
		return boost::make_tuple(TOKEN_OPERATION_TIMES, 0., "") ;
	}
	else if(str == std::string("/"))
	{
		return boost::make_tuple(TOKEN_OPERATION_DIVIDES, 0., "") ;
	}
	else if(str == std::string("^"))
	{
		return boost::make_tuple(TOKEN_OPERATION_POWER, 0., "") ;
	}
	else if(str == std::string("x"))
	{
		return boost::make_tuple(TOKEN_OPERATION_X, 0., "") ;
	}
	else if(str == std::string("y"))
	{
		return boost::make_tuple(TOKEN_OPERATION_Y, 0., "") ;
	}
	else if(str == std::string("z"))
	{
		return boost::make_tuple(TOKEN_OPERATION_Z, 0., "") ;
	}
	else if(str == std::string("t"))
	{
		return boost::make_tuple(TOKEN_OPERATION_T, 0., "") ;
	}
	else if(str == std::string("u"))
	{
		return boost::make_tuple(TOKEN_OPERATION_U, 0., "") ;
	}
	else if(str == std::string("v"))
	{
		return boost::make_tuple(TOKEN_OPERATION_V, 0., "") ;
	}
	else 
	{
		return boost::make_tuple(TOKEN_OPERATION_W, 0., "") ;
	}
}
	
bool Function::isOperator(const char c) const
{
	if(
		c == '+' ||
		c == '-' ||
		c == '/' ||
		c == '*' ||
		c == '^'
		)
		return true ;
	return false ;
} ;
	
bool Function::isNumeral(const char c) const
{
	if(
		c == '0' ||
		c == '1' ||
		c == '2' ||
		c == '3' ||
		c == '4' ||
		c == '5' ||
		c == '6' ||
		c == '7' ||
		c == '8' ||
		c == '9' 
		)
		return true ;
	return false ;
} ;
	
bool Function::isSeparator(const char c) const
{
	if(
		c == ' ' ||
		isOperator(c)
		)
		return true ;
	return false ;
} ;

std::pair<size_t, boost::tuple< TokenOperationType, double, std::string>> Function::getNext(size_t init, const std::string & form)
{
	return getNext(init, form.c_str()) ;
}

std::pair<size_t, boost::tuple<TokenOperationType, double, std::string>> Function::getNext(size_t init, const char * form)
{
	std::string cur_token("") ;
	char cur_char  = form[init] ;
	while(isSeparator(cur_char) && !isOperator(cur_char))
	{
		init++ ;
		cur_char  = form[init] ;
	}
	
	if(isOperator(cur_char))
	{
		if(cur_char != '-')
		{
			cur_token+=cur_char ;
			init++ ;
			return std::make_pair(init, toToken(cur_token)) ;
		}
		else if (isNumeral(form[init+1]))
		{
			cur_token+=cur_char ;
			init++ ;
			cur_char  = form[init] ;
			while(!isSeparator(cur_char))
			{
				if(init < strlen(form))
				{
					cur_token+=cur_char ;
					init++ ;
					cur_char  = form[init] ;
				}
				else
					break ;
			}
			return std::make_pair(init, toToken(cur_token)) ;
		}
		else
		{
			cur_token+=cur_char ;
			init++ ;
			return std::make_pair(init, toToken(cur_token)) ;
		}
	}
	while(!isSeparator(cur_char))
	{
		if(init < strlen(form))
		{
			cur_token+=cur_char ;
			init++ ;
			cur_char  = form[init] ;
		}
		else
			break ;
	}
	return std::make_pair(init, toToken(cur_token)) ;
}

void Function::initialiseAdresses(size_t offset)
{
	unsigned short int da = 8 + (offset>0);
	int counter = offset ;
	for(size_t i = offset ; i < byteCodeSize ; i++)
	{
		switch(byteCode[i])
		{
			case TOKEN_OPERATION_CONSTANT:
			{
				adress_a[4*counter+2] = da ;
				adress_a[4*counter++] = da++ ;
				break ;
			}
			case TOKEN_OPERATION_X:
			{
				adress_a[4*counter+2] = da ;
				adress_a[4*counter++] = da++ ;
				break ;
			}
			case TOKEN_OPERATION_Y:
			{
				adress_a[4*counter+2] = da ;
				adress_a[4*counter++] = da++ ;
				break ;
			}
			case TOKEN_OPERATION_Z:
			{
				adress_a[4*counter+2] = da ;
				adress_a[4*counter++] = da++ ;
				break ;
			}
			case TOKEN_OPERATION_T:
			{
				adress_a[4*counter+2] = da ;
				adress_a[4*counter++] = da++ ;
				break ;
			}
			case TOKEN_OPERATION_U:
			{
				adress_a[4*counter+2] = da ;
				adress_a[4*counter++] = da++ ;
				break ;
			}
			case TOKEN_OPERATION_V:
			{
				adress_a[4*counter+2] = da ;
				adress_a[4*counter++] = da++ ;
				break ;
			}
			case TOKEN_OPERATION_W:
			{
				adress_a[4*counter+2] = da ;
				adress_a[4*counter++] = da++ ;
				break ;
			}
			case TOKEN_OPERATION_COS:
			{
				adress_a[4*counter+2] = da-1 ;
				adress_a[4*counter++] = da-1 ;
				break ;
			}
			case TOKEN_OPERATION_ABS:
			{
				adress_a[4*counter+2] = da-1 ;
				adress_a[4*counter++] = da-1 ;
				break ;
			}
			case TOKEN_OPERATION_TAN:
			{
				adress_a[4*counter+2] = da-1 ;
				adress_a[4*counter++] = da-1 ;
				break ;
			}
			case TOKEN_OPERATION_SIN:
			{
				adress_a[4*counter+2] = da-1 ;
				adress_a[4*counter++] = da-1 ;
				break ;
			}
			case TOKEN_OPERATION_EXP:
			{
				adress_a[4*counter+2] = da-1 ;
				adress_a[4*counter++] = da-1 ;
				break ;
			}
			case TOKEN_OPERATION_SIGN:
			{
				adress_a[4*counter+2] = da-1 ;
				adress_a[4*counter++] = da-1 ;
				break ;
			}
			case TOKEN_OPERATION_POSITIVITY:
			{
				adress_a[4*counter+2] = da-1 ;
				adress_a[4*counter++] = da-1 ;
				break ;
			}
			case TOKEN_OPERATION_NEGATIVITY:
			{
				adress_a[4*counter+2] = da-1 ;
				adress_a[4*counter++] = da-1 ;
				break ;
			}
			case TOKEN_OPERATION_LOG:
			{
				adress_a[4*counter+2] = da-1 ;
				adress_a[4*counter++] = da-1 ;
				break ;
			}
			case TOKEN_OPERATION_COSH:
			{
				adress_a[4*counter+2] = da-1 ;
				adress_a[4*counter++] = da-1 ;
				break ;
			}
			case TOKEN_OPERATION_SINH:
			{
				adress_a[4*counter+2] = da-1 ;
				adress_a[4*counter++] = da -1;
				break ;
			}
			case TOKEN_OPERATION_TANH:
			{
				adress_a[4*counter+2] = da-1 ;
				adress_a[4*counter++] = da-1 ;
				break ;
			}
			case TOKEN_OPERATION_SQRT:
			{
				adress_a[4*counter+2] = da-1 ;
				adress_a[4*counter++] = da-1 ;
				break ;
			}
			case TOKEN_OPERATION_BESSEL:
			{
				adress_a[counter+2] = da-1 ;
				adress_a[counter++] = da-1 ;
				break ;
			}
			case TOKEN_OPERATION_ATAN2:
			{
				adress_a[4*counter] = da-2 ;
				adress_a[4*counter+1] = da-1 ;
				adress_a[4*counter++] = da-2 ;
				--da ;
				break ;
			}
			case TOKEN_OPERATION_INTERPOLATE:
			{
				adress_a[4*counter+2] = da-2 ;
				adress_a[4*counter+1] = da-1 ;
				adress_a[4*counter++] = da-2 ;
				--da ;
				break ;
			}
			case TOKEN_OPERATION_PLUS:
			{
				adress_a[4*counter+2] = da-2 ;
				adress_a[4*counter+1] = da-1 ;
				adress_a[4*counter++] = da-2 ;
				--da ;
				break ;
			}
			case TOKEN_OPERATION_MINUS:
			{
				adress_a[4*counter+2] = da-2 ;
				adress_a[4*counter+1] = da-1 ;
				adress_a[4*counter++] = da-2 ;
				--da ;
				break ;
			}
			case TOKEN_OPERATION_TIMES:
			{
				adress_a[4*counter+2] = da-2 ;
				adress_a[4*counter+1] = da-1 ;
				adress_a[4*counter++] = da-2 ;
				--da ;
				break ;
			}
			case TOKEN_OPERATION_DIVIDES:
			{
				adress_a[4*counter+2] = da-2 ;
				adress_a[4*counter+1] = da-1 ;
				adress_a[4*counter++] = da-2 ;
				--da ;
				break ;
			}
			case TOKEN_OPERATION_POWER:
			{
				adress_a[4*counter+2] = da-2 ;
				adress_a[4*counter+1] = da-1 ;
				adress_a[4*counter++] = da-2 ;
				--da ;
				break ;
			}
			case TOKEN_OPERATION_GEO_OPERATION:
			{
				int offset = geo_op[i]->adressOffset() ;
				if(offset == -1) 
				{
					adress_a[4*counter+2] = da-2 ;
					adress_a[4*counter+1] = da-1 ;
					adress_a[4*counter++] = da-2 ;
					--da ;
				}
				
				if(offset == -2)
				{
					adress_a[4*counter+2] = da-3 ;
					adress_a[4*counter+1] = da-2 ;
					adress_a[4*counter++] = da-1 ;
					da -= 2 ;
				}
				if(offset == 0) 
				{
					adress_a[4*counter+2] = da-3 ;
					adress_a[4*counter+1] = da-2 ;
					adress_a[4*counter++] = da-1 ;
				}
			}
			default:
			break ;
		}
	}
	
	std::valarray<TokenOperationType> newbyteCode ;
	std::valarray<GeometryOperation *> newgeo_op ;
	std::valarray<double> newvalues ;
	std::valarray<short unsigned int> newadress_a ;
	
	newbyteCode.resize(FUNCTION_LENGTH, TOKEN_OPERATION_CONSTANT); 
	newgeo_op.resize(FUNCTION_LENGTH, (GeometryOperation *)nullptr); 
	newvalues.resize(FUNCTION_LENGTH, 0.) ;
	newadress_a.resize(FUNCTION_LENGTH*4, 0) ;
	newbyteCode = byteCode ; 
	newgeo_op= geo_op; 
	newvalues= values ;
	newadress_a = adress_a;
	size_t constcounter = constNumber ;
	for(size_t i = offset ; i < byteCodeSize ; i++)
	{
		switch(byteCode[i])
		{
			case TOKEN_OPERATION_X:
			{
				bool found = false ;
				for(size_t j = i+1 ; j < byteCodeSize ; j++)
				{
					if(adress_a[i*4] == adress_a[j*4])
					{
						adress_a[j*4] = 1 ;
						found = true ;
						break ;
					}
					if(adress_a[i*4] == adress_a[j*4+1])
					{
						adress_a[j*4+1] = 1 ;
						found = true ;
						break ;
					}
					if(adress_a[i*4] == adress_a[j*4+2])
					{
						adress_a[j*4+2] = 1 ;
						found = true ;
						break ;
					}
				}
				if(!found)
				{
					adress_a[i*4] = 1 ;
				}
				break ;
			}
			case TOKEN_OPERATION_Y:
			{
				bool found = false ;
				for(size_t j = i+1 ; j < byteCodeSize ; j++)
				{
					if(adress_a[i*4] == adress_a[j*4])
					{
						adress_a[j*4] = 2 ;
						found = true ;
						break ;
					}
					if(adress_a[i*4] == adress_a[j*4+1])
					{
						adress_a[j*4+1] = 2 ;
						found = true ;
						break ;
					}
					if(adress_a[i*4] == adress_a[j*4+2])
					{
						adress_a[j*4+2] = 2 ;
						found = true ;
						break ;
					}
				}
				if(!found)
				{
					adress_a[i*4] = 2 ;
				}
				break ;
			}
			case TOKEN_OPERATION_Z:
			{
				bool found = false ;
				for(size_t j = i+1 ; j < byteCodeSize ; j++)
				{
					if(adress_a[i*4] == adress_a[j*4])
					{
						adress_a[j*4] = 3 ;
						found = true ;
						break ;
					}
					if(adress_a[i*4] == adress_a[j*4+1])
					{
						adress_a[j*4+1] = 3 ;
						found = true ;
						break ;
					}
					if(adress_a[i*4] == adress_a[j*4+2])
					{
						adress_a[j*4+2] = 3 ;
						found = true ;
						break ;
					}
				}
				if(!found)
				{
					adress_a[i*4] = 3 ;
				}
				break ;
			}
			case TOKEN_OPERATION_T:
			{
				bool found = false ;
				for(size_t j = i+1 ; j < byteCodeSize ; j++)
				{
					if(adress_a[i*4] == adress_a[j*4])
					{
						adress_a[j*4] = 4 ;
						found = true ;
						break ;
					}
					if(adress_a[i*4] == adress_a[j*4+1])
					{
						adress_a[j*4+1] = 4 ;
						found = true ;
						break ;
					}
					if(adress_a[i*4] == adress_a[j*4+2])
					{
						adress_a[j*4+2] = 4 ;
						found = true ;
						break ;
					}
				}
				if(!found)
				{
					adress_a[i*4] = 4 ;
				}
				break ;
			}
			case TOKEN_OPERATION_U:
			{
				bool found = false ;
				for(size_t j = i+1 ; j < byteCodeSize ; j++)
				{
					if(adress_a[i*4] == adress_a[j*4])
					{
						adress_a[j*4] = 5 ;
						found = true ;
						break ;
					}
					if(adress_a[i*4] == adress_a[j*4+1])
					{
						adress_a[j*4+1] = 5 ;
						found = true ;
						break ;
					}
					if(adress_a[i*4] == adress_a[j*4+2])
					{
						adress_a[j*4+2] = 5 ;
						found = true ;
						break ;
					}
				}
				if(!found)
				{
					adress_a[i*4] = 5 ;
				}
				break ;
			}
			case TOKEN_OPERATION_V:
			{
				bool found = false ;
				for(size_t j = i+1 ; j < byteCodeSize ; j++)
				{
				if(adress_a[i*4] == adress_a[j*4])
					{
						adress_a[j*4] = 6 ;
						found = true ;
						break ;
					}
					if(adress_a[i*4] == adress_a[j*4+1])
					{
						adress_a[j*4+1] = 6 ;
						found = true ;
						break ;
					}
					if(adress_a[i*4] == adress_a[j*4+2])
					{
						adress_a[j*4+2] = 6 ;
						found = true ;
						break ;
					}
				}
				if(!found)
				{
					adress_a[i*4] = 6 ;
				}
				break ;
			}
			case TOKEN_OPERATION_W:
			{
				bool found = false ;
				for(size_t j = i+1 ; j < byteCodeSize ; j++)
				{
					if(adress_a[i*4] == adress_a[j*4])
					{
						adress_a[j*4] = 7 ;
						found = true ;
						break ;
					}
					if(adress_a[i*4] == adress_a[j*4+1])
					{
						adress_a[j*4+1] = 7 ;
						found = true ;
						break ;
					}
					if(adress_a[i*4] == adress_a[j*4+2])
					{
						adress_a[j*4+2] = 7 ;
						found = true ;
						break ;
					}
				}
				if(!found)
				{
					adress_a[i*4] = 7 ;
				}
				
				break ;
			}
			case TOKEN_OPERATION_CONSTANT:
			{
				bool found = false ;
				for(size_t j = i+1 ; j < byteCodeSize ; j++)
				{
					if(adress_a[i*4] == adress_a[j*4])
					{
						newvalues[constcounter] = values[i] ;
						adress_a[j*4] = FUNCTION_LENGTH-1-constcounter++ ;
						found = true ;
						break ;
					}
					else if(adress_a[i*4] == adress_a[j*4+1])
					{
						newvalues[constcounter] = values[i] ;
						adress_a[j*4+1] = FUNCTION_LENGTH-1-constcounter++ ;
						found = true ;
						break ;
					}
					else if(adress_a[i*4] == adress_a[j*4+2])
					{
						newvalues[constcounter] = values[i] ;
						adress_a[j*4+2] = FUNCTION_LENGTH-1-constcounter++ ;
						found = true ;
						break ;
					}
				}
				if(!found)
				{
					newvalues[constcounter] = values[i] ;
					adress_a[i*4] = FUNCTION_LENGTH-1-constcounter++ ;
				}
				
				break ;
			}
		}
	}
	constNumber = constcounter ;
	
	size_t newByteCodeSize = byteCodeSize ;
	counter = offset ;
	if(byteCodeSize > 1)
	{
		for(size_t i = offset ; i < byteCodeSize ; i++)
		{
			if(byteCode[i] == TOKEN_OPERATION_X ||
				byteCode[i] == TOKEN_OPERATION_Y ||
				byteCode[i] == TOKEN_OPERATION_Z ||
				byteCode[i] == TOKEN_OPERATION_T ||
				byteCode[i] == TOKEN_OPERATION_W ||
				byteCode[i] == TOKEN_OPERATION_U ||
				byteCode[i] == TOKEN_OPERATION_V ||
				byteCode[i] == TOKEN_OPERATION_CONSTANT 
			)
			{
				newByteCodeSize-- ;
				continue ;
			}
			if(adress_a[i*4] == adress_a[i*4+2] && i < byteCodeSize-1 && i < newByteCodeSize)
			{
				if(byteCode[i] == TOKEN_OPERATION_PLUS)
					byteCode[i] = TOKEN_OPERATION_INPLACE_PLUS ;
				if(byteCode[i] == TOKEN_OPERATION_MINUS)
					byteCode[i] = TOKEN_OPERATION_INPLACE_MINUS ;
				if(byteCode[i] == TOKEN_OPERATION_TIMES)
					byteCode[i] = TOKEN_OPERATION_INPLACE_TIMES ;
				if(byteCode[i] == TOKEN_OPERATION_DIVIDES)
					byteCode[i] = TOKEN_OPERATION_INPLACE_DIVIDES ;
				if(byteCode[i] == TOKEN_OPERATION_POWER)
					byteCode[i] = TOKEN_OPERATION_INPLACE_POWER ;
			}

			newbyteCode[counter] = byteCode[i]; 
			newgeo_op[counter]=geo_op[i]; 
			newadress_a[counter*4]= adress_a[i*4];
			newadress_a[counter*4+1]= adress_a[i*4+1];
			newadress_a[counter*4+2]= adress_a[i*4+2];
			counter++ ;

		}
		
		byteCodeSize = newByteCodeSize ;
		byteCode = newbyteCode; 
		geo_op= newgeo_op; 
		values= newvalues;
		adress_a= newadress_a;
		return ;
		for(size_t i = offset ; i < byteCodeSize ; i++)
		{
			if(byteCode[i] < TOKEN_OPERATION_PLUS || byteCode[i] > TOKEN_OPERATION_INPLACE_POWER || byteCode[i+1] < TOKEN_OPERATION_PLUS || byteCode[i+1] > TOKEN_OPERATION_INPLACE_POWER || i >= byteCodeSize-1)
				continue ;
			if((adress_a[i*4+2] == adress_a[(i+1)*4] || 
				adress_a[i*4+2] == adress_a[(i+1)*4+1] )&&
				use_temp[i] == NO_TEMPORARY
			)
			{
				use_temp[i] = SET_TEMPORARY ;
				if(adress_a[i*4+2] == adress_a[(i+1)*4])
					use_temp[i+1] = GET_TEMPORARY_A ;
				else
					use_temp[i+1] = GET_TEMPORARY_B ;
			}
			
			if((adress_a[i*4+2] == adress_a[(i+1)*4] || 
				adress_a[i*4+2] == adress_a[(i+1)*4+1] )&&
				use_temp[i] == GET_TEMPORARY_A
			)
			{
				use_temp[i] = SET_GET_TEMPORARY_A ;
				if(adress_a[i*4+2] == adress_a[(i+1)*4])
					use_temp[i+1] = GET_TEMPORARY_A ;
				else
					use_temp[i+1] = GET_TEMPORARY_B ;
			}
			
			if((adress_a[i*4+2] == adress_a[(i+1)*4] || 
				adress_a[i*4+2] == adress_a[(i+1)*4+1] )&&
				use_temp[i] == GET_TEMPORARY_B
			)
			{
				use_temp[i] = SET_GET_TEMPORARY_B ;
				if(adress_a[i*4+2] == adress_a[(i+1)*4])
					use_temp[i+1] = GET_TEMPORARY_A ;
				else
					use_temp[i+1] = GET_TEMPORARY_B ;
			}
		}
	}
}
Function::Function(const char *f)
{	
	defaultInitialise() ;

	size_t init = 0 ;

	while(init < strlen(f))
	{
		std::pair<size_t, boost::tuple<TokenOperationType, double, std::string> > temp = getNext(init, f) ;
		byteCode[byteCodeSize] = boost::get<0>(temp.second) ;
		values[byteCodeSize++]   = boost::get<1>(temp.second) ;
		init = temp.first ;
	}
	initialiseAdresses();
}

Function::Function(const std::string &f) 
{
	defaultInitialise() ;
	size_t init = 0 ;
		
	while(init < f.length())
	{
		std::pair<size_t, boost::tuple<TokenOperationType, double, std::string> > temp = getNext(init, f) ;
		byteCode[byteCodeSize] = boost::get<0>(temp.second) ;
		values[byteCodeSize++] = boost::get<1>(temp.second) ;
		init = temp.first ;
	}
	initialiseAdresses();
}

Function::Function(const std::valarray<Matrix> & coeffs, bool diff)
{
	defaultInitialise() ;
	e_diff = diff ;
	
	
	bool first = true ;
	
	for(size_t i =  0;  i < coeffs.size() ; i++)
	{
		for(size_t k =  0;  k < coeffs[i].numRows() ; k++)
		{
			for(size_t j = 0 ; j < coeffs[i].numCols() ;  j++)
			{

					if(coeffs[i][j][k] != 0)
					{
						if(std::abs(coeffs[i][j][k]-1) > POINT_TOLERANCE_2D || i == 0 && j == 0 && k == 0)
						{
							byteCode[byteCodeSize] = TOKEN_OPERATION_CONSTANT ;
							values[byteCodeSize] = coeffs[i][j][k] ;
							byteCodeSize++;
						}
						
						if(i > 0)
						{
							byteCode[byteCodeSize++] = TOKEN_OPERATION_X ;
							if(i > 1)
							{
								byteCode[byteCodeSize] = TOKEN_OPERATION_CONSTANT ;
								values[byteCodeSize++] = i ;
								byteCode[byteCodeSize++] = TOKEN_OPERATION_POWER ;
							}
							if(std::abs(coeffs[i][j][k]-1) > POINT_TOLERANCE_2D)
							{
								byteCode[byteCodeSize] = TOKEN_OPERATION_TIMES;
								byteCodeSize++ ;
							}
						}
						
						if(j > 0)
						{
							byteCode[byteCodeSize++] = TOKEN_OPERATION_Y ;
							if(j > 1)
							{
								byteCode[byteCodeSize] = TOKEN_OPERATION_CONSTANT ;
								values[byteCodeSize++] = j ;
								byteCode[byteCodeSize++] = TOKEN_OPERATION_POWER ;
							}
							if(std::abs(coeffs[i][j][k]-1) > POINT_TOLERANCE_2D || i)
							{
								byteCode[byteCodeSize] = TOKEN_OPERATION_TIMES;
								byteCodeSize++ ;
							}
						}
						
						if(k > 0)
						{
							byteCode[byteCodeSize++] = TOKEN_OPERATION_Z ;
							if(k > 1)
							{
								byteCode[byteCodeSize] = TOKEN_OPERATION_CONSTANT ;
								values[byteCodeSize++] = k ;
								byteCode[byteCodeSize++] = TOKEN_OPERATION_POWER ;
							}
							if(std::abs(coeffs[i][j][k]-1) > POINT_TOLERANCE_2D || i || j)
							{
								byteCode[byteCodeSize] = TOKEN_OPERATION_TIMES;
								byteCodeSize++ ;
							}
						}
						
						if(!first)
							byteCode[byteCodeSize++] = TOKEN_OPERATION_PLUS ;
						else
							first = false ;
						
				}
			}
		}
	}
	

	if(first)
	{
		byteCode[0] = TOKEN_OPERATION_CONSTANT ;
		values[0] = 0 ;
		byteCodeSize = 1 ;
	}
	
	if(e_diff)
	{
		derivative = new std::valarray<Function *>((Function *)nullptr, 3) ;
		std::valarray<Matrix> dx(Matrix (coeffs[0].numRows()-1,coeffs[0].numCols() ),coeffs.size()-1) ;
		for(size_t i = 0 ; i < dx.size() ; i++)
		{
			dx[i] = coeffs[i+1]*(i+1) ;
		}
		
		(*derivative)[XI] = new Function(dx, false) ;
		std::valarray<Matrix> dy(Matrix (coeffs[0].numRows()-1,coeffs[0].numCols() ), coeffs.size() ) ;
		for(size_t i = 0 ; i < dy.size() ; i++)
		{
			Matrix dy_(coeffs[i].numRows()-1,coeffs[i].numCols() ) ;
			for(size_t j = 0 ; j < coeffs[i].numRows()-1 ;  j++)
			{
				for(size_t k = 0 ; k < coeffs[i].numCols() ;  k++)
				{
					dy_[j][k] = coeffs[i][j+1][k]*(j+1) ;
				}
			}
			dy[i] = dy_ ;
		}
		(*derivative)[ETA] = new Function(dy, false) ;
		
		std::valarray<Matrix> dz(Matrix(coeffs[0].numRows(),coeffs[0].numCols()-1 ), coeffs.size()) ;
		for(size_t i = 0 ; i < dz.size() ; i++)
		{
			Matrix dz_(coeffs[i].numRows(),coeffs[i].numCols()-1 ) ;
			for(size_t j = 0 ; j < coeffs[i].numRows() ;  j++)
			{
				for(size_t k = 0 ; k < coeffs[i].numCols()-1 ;  k++)
				{
					dz_[j][k] = coeffs[i][j][k+1]*(k+1) ;
				}
			}
			dz[i] = dz_ ;
		}
		(*derivative)[ZETA] = new Function(dz, false) ;
	}
	initialiseAdresses();
}


Function::Function(const std::valarray<std::valarray<Matrix> > & coeffs, bool diff) 
{
	defaultInitialise() ;
	e_diff = diff ;
	
	bool first = true ;
	
	for(size_t i =  0;  i < coeffs.size() ; i++)
	{
		for(size_t j =  0;  j < coeffs[i].size() ; j++)
		{
			for(size_t k = 0 ; k < coeffs[i][j].numRows() ;  k++)
			{
				for(size_t l = 0 ; l < coeffs[i][j].numCols() ;  l++)
				{
					if(coeffs[i][j][k][l] != 0)
					{
						byteCode[byteCodeSize] = TOKEN_OPERATION_CONSTANT ;
						values[byteCodeSize++] = coeffs[i][j][k][l] ;
						if(i > 0)
						{
							byteCode[byteCodeSize++] = TOKEN_OPERATION_X ;
							if(i > 1)
							{
								byteCode[byteCodeSize] = TOKEN_OPERATION_CONSTANT ;
								values[byteCodeSize++] = i ;
								byteCode[byteCodeSize++] = TOKEN_OPERATION_POWER ;
							}
							byteCode[byteCodeSize++] = TOKEN_OPERATION_TIMES;
						}
						
						if(j > 0)
						{
							byteCode[byteCodeSize++] = TOKEN_OPERATION_Y ;
							if(j > 1)
							{
								byteCode[byteCodeSize] = TOKEN_OPERATION_CONSTANT ;
								values[byteCodeSize++] = j ;
								byteCode[byteCodeSize++] = TOKEN_OPERATION_POWER ;
							}
							byteCode[byteCodeSize++] = TOKEN_OPERATION_TIMES;
						}
						
						if(k > 0)
						{
							byteCode[byteCodeSize++] = TOKEN_OPERATION_Z ;
							if(k > 1)
							{
								byteCode[byteCodeSize] = TOKEN_OPERATION_CONSTANT ;
								values[byteCodeSize++] = k ;
								byteCode[byteCodeSize++] = TOKEN_OPERATION_POWER ;
							}
							byteCode[byteCodeSize++] = TOKEN_OPERATION_TIMES;
						}
						
						if(l > 0)
						{
							byteCode[byteCodeSize++] = TOKEN_OPERATION_T ;
							if(l > 1)
							{
								byteCode[byteCodeSize] = TOKEN_OPERATION_CONSTANT ;
								values[byteCodeSize++] = l ;
								byteCode[byteCodeSize++] = TOKEN_OPERATION_POWER ;
							}
							byteCode[byteCodeSize++] = TOKEN_OPERATION_TIMES;
						}
						
						if(!first)
							byteCode[byteCodeSize++] = TOKEN_OPERATION_PLUS;
						else
							first = false ;
	
					}
				}
			}
		}
	}
	
	if(first)
	{
		byteCode[0] = TOKEN_OPERATION_CONSTANT ;
		values[0] = 0 ;
		byteCodeSize = 1 ;
	}
	
	if(e_diff)
	{
		derivative = new std::valarray<Function *>((Function *)nullptr, 4*diff);
		std::valarray< std::valarray<Matrix> > dx(coeffs.size()-1) ;
		for(size_t i = 0 ; i < dx.size() ; i++)
		{
			std::valarray<Matrix> dx_(Matrix(coeffs[i][0].numRows(),coeffs[i][0].numCols() ), coeffs.size()) ;
			for(size_t j = 0 ; j < dx_.size() ; j++)
			{
				Matrix dx__(coeffs[i][j].numRows(),coeffs[i][j].numCols() ) ;
				for(size_t k = 0 ; k < coeffs[i][j].numRows() ;  k++)
				{
					for(size_t l = 0 ; l < coeffs[i][j].numCols() ;  l++)
					{
						dx__[k][l] = coeffs[i+1][j][k][l]*(i+1) ;
					}
				}
				dx_[j] = dx__ ;
			}
			dx[i] = dx_ ;
		}
		(*derivative)[XI] = new Function(dx, false) ;
		
		std::valarray< std::valarray<Matrix> > dy(coeffs.size()) ;
		for(size_t i = 0 ; i < dy.size() ; i++)
		{
			std::valarray<Matrix> dy_(Matrix(coeffs[i][0].numRows(),coeffs[i][0].numCols() ), coeffs[i].size()-1) ;
			for(size_t j = 0 ; j < dy_.size() ; j++)
			{
				Matrix dy__(coeffs[i][j].numRows(),coeffs[i][j].numCols() ) ;
				for(size_t k = 0 ; k < coeffs[i][j].numRows()-1 ;  k++)
				{
					for(size_t l = 0 ; l < coeffs[i][j].numCols() ;  l++)
					{
						dy__[k][l] = coeffs[i][j+1][k][l]*(j+1) ;
					}
				}
				dy_[j] = dy__ ;
			}
			dy[i] = dy_ ;
		}
		(*derivative)[ETA] = new Function(dy, false) ;
		
		std::valarray< std::valarray<Matrix> > dz(coeffs.size()) ;
		for(size_t i = 0 ; i < dz.size() ; i++)
		{
			std::valarray<Matrix> dz_(Matrix(coeffs[i][0].numRows()-1,coeffs[i][0].numCols() ), coeffs[i].size()) ;
			for(size_t j = 0 ; j < dz_.size() ; j++)
			{
				Matrix dz__(coeffs[i][j].numRows()-1,coeffs[i][j].numCols() ) ;
				for(size_t k = 0 ; k < coeffs[i][j].numRows()-1 ;  k++)
				{
					for(size_t l = 0 ; l < coeffs[i][j].numCols() ;  l++)
					{
						dz__[k][l] = coeffs[i][j][k+1][l]*(k+1) ;
					}
				}
				dz[i] = dz_ ;
			}
		}
		(*derivative)[ZETA] = new Function(dz, false) ;
		
		std::valarray< std::valarray<Matrix> > dt(coeffs.size()) ;
		for(size_t i = 0 ; i < dt.size() ; i++)
		{
			std::valarray<Matrix> dt_(Matrix(coeffs[i][0].numRows(),coeffs[i][0].numCols()-1 ),coeffs[i].size()) ;
			
			for(size_t j = 0 ;  j <  dt_.size() ; j++)
			{
				
				Matrix dt__(coeffs[i][j].numRows(),coeffs[i][j].numCols()-1 ) ;
				for(size_t k = 0 ; k < coeffs[i][j].numRows() ;  k++)
				{
					for(size_t l = 0 ; l < coeffs[i][j].numCols()-1 ;  l++)
					{
						dt__[k][l] = coeffs[i][j][k][l+1]*(l+1) ;
					}
				}
				dt_[j] = dt__ ;
			}
			
			dt[i] = dt_ ;
		}
		(*derivative)[TIME_VARIABLE] = new Function(dt, false) ;
	}
	initialiseAdresses();
}

Function::Function(const Matrix & coeffs, bool diff)
{
	defaultInitialise() ;
	e_diff = diff ;

	
	bool first = true ;
	
	for(size_t j = 0 ; j < coeffs.numRows() ;  j++)
	{
		for(size_t k = 0 ; k < coeffs.numCols() ;  k++)
		{
			if(std::abs(coeffs[j][k]) > POINT_TOLERANCE_2D)
			{
				if(std::abs(coeffs[j][k]-1) > POINT_TOLERANCE_2D || j == 0 && k == 0)
				{
					byteCode[byteCodeSize] = TOKEN_OPERATION_CONSTANT ;
					values[byteCodeSize] = coeffs[j][k] ;
					byteCodeSize++;
				}
				
				if(j > 0)
				{
					byteCode[byteCodeSize] = TOKEN_OPERATION_X ;
					byteCodeSize++;
					if(j > 1)
					{
						byteCode[byteCodeSize] = TOKEN_OPERATION_CONSTANT ;
						values[byteCodeSize] = j ;
						byteCodeSize++ ;
						byteCode[byteCodeSize] = TOKEN_OPERATION_POWER ;
						byteCodeSize++ ;
					}
					if(std::abs(coeffs[j][k]-1) > POINT_TOLERANCE_2D)
					{
						byteCode[byteCodeSize] = TOKEN_OPERATION_TIMES;
						byteCodeSize++ ;
					}
				}
				
				if(k > 0)
				{
					byteCode[byteCodeSize] = TOKEN_OPERATION_Y ;
					byteCodeSize++ ; 
					if(k > 1)
					{
						byteCode[byteCodeSize] = TOKEN_OPERATION_CONSTANT ;
						values[byteCodeSize] = k ;
						byteCodeSize++ ;
						byteCode[byteCodeSize] = TOKEN_OPERATION_POWER ;
						byteCodeSize++ ;
					}
					if(std::abs(coeffs[j][k]-1) > POINT_TOLERANCE_2D || j)
					{
						byteCode[byteCodeSize] = TOKEN_OPERATION_TIMES;
						byteCodeSize++ ;
					}
				}
				
				if(!first)
				{
					byteCode[byteCodeSize] = TOKEN_OPERATION_PLUS;
					byteCodeSize++ ;
				}
				else
					first = false ;
			}
		}
	}
/*	
		std::cout << "initial bytecode with no adresses-------" << std::endl ;
		VirtualMachine().print(*this);
		std::cout << "-------initial bytecode with no adresses"<< std::endl ;*/
		
	if(e_diff)
	{
		derivative = new std::valarray<Function *>((Function *)nullptr,2*diff);
		Matrix dx(coeffs.numRows()-1, coeffs.numCols()) ;
		size_t imax = dx.numRows() ;
		size_t jmax = dx.numCols() ;
		for(size_t i = 0 ; i < imax ; i++)
		{
			for(size_t j = 0 ; j < jmax ; j++)
			{
				dx[i][j] = coeffs[i+1][j]*(i+1) ;
			}
			
		}
		
		(*derivative)[XI] = new Function(dx, false) ;
		
		Matrix dy(coeffs.numRows(), coeffs.numCols()-1) ;
		
		imax = dy.numRows() ;
		jmax = dy.numCols() ;
		
		for(size_t i = 0 ; i < imax ; i++)
		{
			for(size_t j = 0 ; j < jmax ; j++)
			{
				dy[i][j] = coeffs[i][j+1]*(j+1) ;
			}
		}
		
		(*derivative)[ETA] = new Function(dy, false) ;
	}
	initialiseAdresses();
}

bool Function::isNull() const
{
	return (byteCodeSize == 0) ;
}

Function::Function(const std::valarray<double> & coeffs, bool diff)
{

	defaultInitialise() ;
	e_diff = diff ;
	
	
	for(size_t j = 0 ; j < coeffs.size() ;  j++)
	{

			if(coeffs[j] != 0)
			{
				byteCode[byteCodeSize] = TOKEN_OPERATION_CONSTANT ;
				values[byteCodeSize++] = coeffs[j] ;
				
				if(j > 0)
				{
					byteCode[byteCodeSize++] = TOKEN_OPERATION_X ;
					if(j > 1)
					{
						byteCode[byteCodeSize] = TOKEN_OPERATION_CONSTANT ;
						values[byteCodeSize++] = j ;
						byteCode[byteCodeSize++] = TOKEN_OPERATION_POWER ;
					}
					byteCode[byteCodeSize++] = TOKEN_OPERATION_TIMES;
				}

				byteCode[byteCodeSize++] = TOKEN_OPERATION_PLUS;

			}

	}

	
	if(e_diff)
	{
		derivative = new std::valarray<Function *>((Function *)nullptr, 1) ;
		size_t ds = 0 ;
		if((int)coeffs.size()-1 > 0)
			ds = coeffs.size()-1 ;
		std::valarray<double> dx(ds) ;
// 		std::cout << "coeffs.size()-1 = "<< coeffs.size()-1 << std::endl ;
		for(size_t i = 1 ; i < dx.size()+1 ; i++)
		{
			dx[i] = coeffs[i]*(i) ;
		}
		
		derivative[XI] = new Function(dx, false) ;
	}
	
	initialiseAdresses();
}

// Function Function::operator()(const Function & f) const
// {
// 	Function ret(*this) ;
// 	return ret ;
// 	Function * currentTransform = ret.xtransform ;
// 	if(currentTransform)
// 	{
// 		while(currentTransform)
// 		{
// 			if(!currentTransform->xtransform)
// 			{
// 				currentTransform->xtransform = new Function(f) ;
// 				currentTransform->transformed = true ;
// 			}
// 		}
// 	}
// 	else
// 	{
// 		ret.xtransform = new Function(f) ;
// 	}
// 	ret.transformed = true ;
// 	
// }

// Function Function::operator()(const Function & f0, const Function & f1) const
// {
// 	Function ret(*this) ;
// 	return ret ;
// 	Function * currentTransform = ret.xtransform ;
// 	if(currentTransform)
// 	{
// 		while(currentTransform)
// 		{
// 			if(!currentTransform->xtransform)
// 			{
// 				currentTransform->xtransform = new Function(f0) ;
// 				currentTransform->transformed = true ;
// 			}
// 		}
// 	}
// 	else
// 	{
// 		ret.xtransform = new Function(f0) ;
// 	}
// 	currentTransform = ret.ytransform ;
// 	if(currentTransform)
// 	{
// 		while(currentTransform)
// 		{
// 			if(!currentTransform->ytransform)
// 			{
// 				currentTransform->ytransform = new Function(f1) ;
// 				currentTransform->transformed = true ;
// 			}
// 		}
// 	}
// 	else
// 	{
// 		ret.ytransform = new Function(f1) ;
// 	}
// 
// 	ret.transformed = true ;
// 	
// }

Function::Function(const Line & l, ElementarySurface * s) 
{
	defaultInitialise() ;
	Function g = s->getXTransform() ;
	Function f = s->getYTransform() ;
	concatenateFunctions(g, f, *this);
	byteCodeSize++ ;
	adress_a[(byteCodeSize-1)*4+2] = 8 ;
	adress_a[(byteCodeSize-1)*4+1] = 9 ;
	adress_a[(byteCodeSize-1)*4] = 8 ;
	geo_op[byteCodeSize-1] = new LineDistanceOperation(l) ;
	initialiseAdresses();
}

Function::Function(const Point & l,  ElementarySurface * s) 
{
	defaultInitialise() ;
	Function g = s->getXTransform() ;
	Function f = s->getYTransform() ;
	concatenateFunctions(g, f, *this);
	byteCodeSize++ ;
	adress_a[(byteCodeSize-1)*4+2] = 8 ;
	adress_a[(byteCodeSize-1)*4+1] = 9 ;
	adress_a[(byteCodeSize-1)*4] = 8 ;
	byteCode[byteCodeSize-1] = TOKEN_OPERATION_GEO_OPERATION ;
	geo_op[byteCodeSize-1] = new PointDistanceBinaryOperation(l) ;
}

Function::Function(const Point & l,  ElementaryVolume * s) 
{

	defaultInitialise() ;
	Function g = s->getXTransform() ;
	Function f = s->getYTransform() ;
	Function h = s->getZTransform() ;
	concatenateFunctions(g, f, h,  *this);
	byteCodeSize++ ;
	adress_a[(byteCodeSize-1)*4+2] = 8 ;
	adress_a[(byteCodeSize-1)*4+1] = 9 ;
	adress_a[(byteCodeSize-1)*4] = 10 ;
	byteCode[byteCodeSize-1] = TOKEN_OPERATION_GEO_OPERATION;
	geo_op[byteCodeSize-1] = new PointDistanceTrinaryOperation(l) ;
}

Function::Function(double a,  ElementarySurface * s) 
{
	defaultInitialise() ;
	byteCodeSize = 3 ;
	byteCode[byteCodeSize-3] = TOKEN_OPERATION_X;
	byteCode[byteCodeSize-2] = TOKEN_OPERATION_Y;
	byteCode[byteCodeSize-1] = TOKEN_OPERATION_GEO_OPERATION;
	geo_op[byteCodeSize-1] = new RotationBinaryOperation(a) ;
}


Function::Function(double a,const Point & p,   ElementarySurface * s)
{
	defaultInitialise() ;
	Function g = s->getXTransform() ;
	Function f = s->getYTransform() ;
	concatenateFunctions(g, f, *this);
	byteCodeSize++ ;
	adress_a[(byteCodeSize-1)*4+2] = 8 ;
	adress_a[(byteCodeSize-1)*4+1] = 9 ;
	adress_a[(byteCodeSize-1)*4] = 8 ;
	byteCode[byteCodeSize-1] = TOKEN_OPERATION_GEO_OPERATION ;
	geo_op[byteCodeSize-1] = new AngleBinaryOperation(a,p) ;
}

Function::Function( const Geometry * geo, const ElementarySurface * s) 
{
	defaultInitialise() ;
	Function g = s->getXTransform() ;
	Function f = s->getYTransform() ;
	concatenateFunctions(g, f, *this);
	byteCodeSize++ ;
	adress_a[(byteCodeSize-1)*4+2] = 8 ;
	adress_a[(byteCodeSize-1)*4+1] = 9 ;
	adress_a[(byteCodeSize-1)*4] = 8 ;
	byteCode[byteCodeSize-1] = TOKEN_OPERATION_GEO_OPERATION ;
	geo_op[byteCodeSize-1] = new DomainBinaryOperation(geo) ;
}

Function f_project(const Geometry *g, const Function &x, const Function &y)
{
	Function ret ;
	ret.byteCodeSize = y.byteCodeSize +x.byteCodeSize + 1 ;
		for(size_t i = 0 ; i < y.byteCodeSize ; i++)
	{
		ret.byteCode[i] = y.byteCode[i] ;
		ret.values[i] = y.values[i] ;
		if(y.geo_op[i])
			ret.geo_op[i] = y.geo_op[i]->getCopy() ;
		
	}
	for(size_t i = 0 ; i < x.byteCodeSize ; i++)
	{
		ret.byteCode[i+y.byteCodeSize] = x.byteCode[i] ;
		ret.values[i+y.byteCodeSize] = x.values[i] ;
		if(x.geo_op[i])
			ret.geo_op[i+y.byteCodeSize] = x.geo_op[i]->getCopy() ;
	}
	ret.byteCode[ret.byteCodeSize-1] = TOKEN_OPERATION_GEO_OPERATION ;
	ret.geo_op[ret.byteCode.size()-1] = new ProjectionBinaryOperation(g) ;
	ret.initialiseAdresses();
	return ret ;
}

Function::Function(const std::vector<Segment> s , ElementarySurface * u, PositionTokenType t) 
{
	defaultInitialise();

	switch(t)
	{
	case POSITION_TOKEN :
		{
			byteCodeSize = 3 ;
			byteCode[byteCodeSize-3] = TOKEN_OPERATION_X;
			byteCode[byteCodeSize-2] = TOKEN_OPERATION_Y;
			byteCode[byteCodeSize-1] = TOKEN_OPERATION_GEO_OPERATION;
			geo_op[byteCode.size()-1] = new PositionOperation(s) ;

			
			break ;
		}
	case PROJECTION_TOKEN :
		{

			byteCodeSize = s.size()+s.size()-1 ;
			for(size_t i = 0 ; i < s.size() ; i++)
			{
				byteCode[i] = TOKEN_OPERATION_GEO_OPERATION ;
				geo_op[i] = new ProjectionOperation2D(s[i]) ;
			}
			for(size_t i = s.size() ; i < byteCode.size() ; i++)
				byteCode[i] = TOKEN_OPERATION_PLUS ;
			
			break ;
		}
	}
	initialiseAdresses();
}

Function::Function(const Function &f)
{
	if(f.derivative && false)
	{
		derivative = new std::valarray<Function *>((Function *)nullptr, f.derivative->size()) ;
		e_diff = true ;
		for(size_t i = 0 ; i < f.derivative->size() ; i++)
			(*derivative)[i] = new Function(*(*f.derivative)[i]) ;
	}
	else
	{
		derivative = nullptr ;
		e_diff = false ;
	}

	byteCode.resize(FUNCTION_LENGTH, TOKEN_OPERATION_CONSTANT); 
	geo_op.resize(FUNCTION_LENGTH, (GeometryOperation *)nullptr); 
	use_temp.resize(FUNCTION_LENGTH, NO_TEMPORARY);
	values.resize(FUNCTION_LENGTH, 0.) ;
	adress_a.resize(FUNCTION_LENGTH*4, 0) ;
	ptID = f.ptID ;
	dofID = f.dofID ;
	
	byteCode = f.byteCode ;
	values = f.values ;
	byteCodeSize = f.byteCodeSize ;
	use_temp = f.use_temp ;
	adress_a = f.adress_a ;
	constNumber = f.constNumber ;
	iPoint = f.iPoint ;

	for(size_t i = 0 ; i < f.byteCodeSize ; i++)
	{
		delete geo_op[i] ;
		if(f.geo_op[i])
			geo_op[i] = f.geo_op[i]->getCopy() ;
	}
	
	for(auto i = f.precalc.begin() ; i != f.precalc.end() ; ++i)
		precalc[i->first] = new Vector(*i->second) ;
	for(auto i = f.dprecalc.begin() ; i != f.dprecalc.end() ; ++i)
	{
		dprecalc[i->first] =  std::map<Variable, Vector *>() ;
		for(auto j = i->second.begin() ; j != i->second.end() ; ++j)
		{
			dprecalc[i->first][j->first] = new Vector(*j->second) ;
		}
	}
	
}


Function::~Function()
{
	
	for(size_t i = 0 ; i < geo_op.size() ; i++)
	{
		delete geo_op[i] ;
	}
	if(derivative)
	{
		for(size_t i = 0 ; i < derivative->size() ; i++)
		{
			delete (*derivative)[i] ;
		}
		delete derivative ;
	}

	for(auto i = precalc.begin() ; i != precalc.end() ; ++i)
		delete i->second ;
	for(auto i = dprecalc.begin() ; i != dprecalc.end() ; ++i)
	{
		for(auto j = i->second.begin() ; j != i->second.end() ; ++j)
		{
			delete j->second ;
		}
	}
}


bool Function::isDifferentiable() const 
{
	return e_diff ;
}


Function Function::operator*(const Function &f) const
{
	if(f.isNull())
	{
		return Function() ;
	}
	Function ret ;
	concatenateFunctions(*this, f, ret);
	ret.byteCodeSize++ ;
	ret.adress_a[(ret.byteCodeSize-1)*4+2] = 8 ;
	ret.adress_a[(ret.byteCodeSize-1)*4+1] = 9 ;
	ret.adress_a[(ret.byteCodeSize-1)*4] = 8 ;
	ret.byteCode[ret.byteCodeSize-1] = TOKEN_OPERATION_TIMES ;
	return ret ;
}
	
Function Function::operator/(const Function &f) const 
{
	if(f.isNull())
	{
		std::cout  << "Divide By Zero " << std::endl ;
		exit(0) ;
		return *this ;
	}
		
	Function ret ;
	concatenateFunctions(*this, f, ret);
	ret.byteCodeSize++ ;
	ret.adress_a[(ret.byteCodeSize-1)*4+2] = 8 ;
	ret.adress_a[(ret.byteCodeSize-1)*4+1] = 9 ;
	ret.adress_a[(ret.byteCodeSize-1)*4] = 8 ;
	ret.byteCode[ret.byteCodeSize-1] = TOKEN_OPERATION_DIVIDES ;
	return ret ;
}
	
Function Function::operator+(const Function &f) const
{
	if(f.isNull())
		return *this ;
	Function ret ;
	concatenateFunctions(*this, f, ret);
	ret.byteCodeSize++ ;
	ret.adress_a[(ret.byteCodeSize-1)*4+2] = 8 ;
	ret.adress_a[(ret.byteCodeSize-1)*4+1] = 9 ;
	ret.adress_a[(ret.byteCodeSize-1)*4] = 8 ;
	ret.byteCode[ret.byteCodeSize-1] = TOKEN_OPERATION_PLUS ;
	return ret ;
}
	
Function operator-(const double & a, const Function &f)
{
	Function ret ;
	ret.byteCodeSize = f.byteCodeSize+1 ;
	ret.constNumber = 1+f.constNumber ;
	ret.values[0] = a ;
	ret.use_temp = f.use_temp ;
	ret.byteCode = f.byteCode ;
	for(size_t i = 0 ; i < f.constNumber ; i++)
		ret.values[i+1] = f.values[i] ;
	
	for(size_t i = 0 ; i < f.byteCodeSize ; i++)
	{
		if(f.geo_op[i])
			ret.geo_op[i] = f.geo_op[i]->getCopy() ;
			
		if(f.adress_a[i*4] >= FUNCTION_LENGTH-1-f.constNumber)
			ret.adress_a[(i)*4] = f.adress_a[i*4]-1 ;
		
		if(f.adress_a[i*4+1] >= FUNCTION_LENGTH-1-f.constNumber)
			ret.adress_a[(i)*4+1] = f.adress_a[i*4+1]-1 ;
		
		if(f.adress_a[i*4+2] >= FUNCTION_LENGTH-1-f.constNumber)
			ret.adress_a[(i)*4+2] = f.adress_a[i*4+2]-1 ;
		
		if(f.adress_a[i*4] >= 8 && f.adress_a[i*4] < FUNCTION_LENGTH-1-f.constNumber)
			ret.adress_a[(i)*4] = f.adress_a[i*4] ;
		
		if(f.adress_a[i*4+1] >= 8 && f.adress_a[i*4+1] < FUNCTION_LENGTH-1-f.constNumber)
			ret.adress_a[(i)*4+1] = f.adress_a[i*4+1] ;
		
		if(f.adress_a[i*4+2] >= 8 && f.adress_a[i*4+2] < FUNCTION_LENGTH-1-f.constNumber)
			ret.adress_a[(i)*4+2] = f.adress_a[i*4+2] ;
		
		if(f.adress_a[i*4] < 8 )
			ret.adress_a[(i)*4] = f.adress_a[i*4] ;
		
		if(f.adress_a[i*4+1] < 8 )
			ret.adress_a[(i)*4+1] = f.adress_a[i*4+1] ;
		
		if(f.adress_a[i*4+2] < 8 )
			ret.adress_a[(i)*4+2] = f.adress_a[i*4+2] ;

	}
	ret.adress_a[(ret.byteCodeSize-1)*4+2] = 8 ;
	ret.adress_a[(ret.byteCodeSize-1)*4+1] = 8 ;
	ret.adress_a[(ret.byteCodeSize-1)*4] = FUNCTION_LENGTH-1 ;
	ret.byteCode[ret.byteCodeSize-1] = TOKEN_OPERATION_MINUS ;
	return ret ;
	
}

Function operator*(const double & a, const Function &f)
{
	Function ret ;
	ret.byteCodeSize = f.byteCodeSize+1 ;
	ret.constNumber = 1+f.constNumber ;
	ret.values[0] = a ;
	ret.use_temp = f.use_temp ;
	ret.byteCode = f.byteCode ;
	for(size_t i = 0 ; i < f.constNumber ; i++)
		ret.values[i+1] = f.values[i] ;
	
	for(size_t i = 0 ; i < f.byteCodeSize ; i++)
	{

		if(f.geo_op[i])
			ret.geo_op[i] = f.geo_op[i]->getCopy() ;
			
		if(f.adress_a[i*4] >= FUNCTION_LENGTH-1-f.constNumber)
			ret.adress_a[(i)*4] = f.adress_a[i*4]-1 ;
		
		if(f.adress_a[i*4+1] >= FUNCTION_LENGTH-1-f.constNumber)
			ret.adress_a[(i)*4+1] = f.adress_a[i*4+1]-1 ;
		
		if(f.adress_a[i*4+2] >= FUNCTION_LENGTH-1-f.constNumber)
			ret.adress_a[(i)*4+2] = f.adress_a[i*4+2]-1 ;
		
		if(f.adress_a[i*4] >= 8 && f.adress_a[i*4] < FUNCTION_LENGTH-1-f.constNumber)
			ret.adress_a[(i)*4] = f.adress_a[i*4]+1 ;
		if(f.adress_a[i*4+1] >= 8 && f.adress_a[i*4+1] < FUNCTION_LENGTH-1-f.constNumber)
			ret.adress_a[(i)*4+1] = f.adress_a[i*4+1]+1 ;
		if(f.adress_a[i*4+2] >= 8 && f.adress_a[i*4+2] < FUNCTION_LENGTH-1-f.constNumber)
			ret.adress_a[(i)*4+2] = f.adress_a[i*4+2]+1 ;
		if(f.adress_a[i*4] < 8 )
			ret.adress_a[(i)*4] = f.adress_a[i*4] ;
		if(f.adress_a[i*4+1] < 8 )
			ret.adress_a[(i)*4+1] = f.adress_a[i*4+1] ;
		if(f.adress_a[i*4+2] < 8 )
			ret.adress_a[(i)*4+2] = f.adress_a[i*4+2] ;

	}
	ret.adress_a[(ret.byteCodeSize-1)*4+2] = 8 ;
	ret.adress_a[(ret.byteCodeSize-1)*4+1] = 8 ;
	ret.adress_a[(ret.byteCodeSize-1)*4] = FUNCTION_LENGTH-1 ;
	ret.byteCode[ret.byteCodeSize-1] = TOKEN_OPERATION_TIMES ;
	return ret ;
}

Function operator+(const double & a, const Function &f)
{
	Function ret ;
	ret.byteCodeSize = f.byteCodeSize+1 ;
	ret.constNumber = 1+f.constNumber ;
	ret.values[0] = a ;
	ret.use_temp = f.use_temp ;
	ret.byteCode = f.byteCode ;
	for(size_t i = 0 ; i < f.constNumber ; i++)
		ret.values[i+1] = f.values[i] ;
	
	for(size_t i = 0 ; i < f.byteCodeSize ; i++)
	{

		if(f.geo_op[i])
			ret.geo_op[i] = f.geo_op[i]->getCopy() ;
			
		if(f.adress_a[i*4] >= FUNCTION_LENGTH-1-f.constNumber)
			ret.adress_a[(i)*4] = f.adress_a[i*4]-1 ;
		
		if(f.adress_a[i*4+1] >= FUNCTION_LENGTH-1-f.constNumber)
			ret.adress_a[(i)*4+1] = f.adress_a[i*4+1]-1 ;
		
		if(f.adress_a[i*4+2] >= FUNCTION_LENGTH-1-f.constNumber)
			ret.adress_a[(i)*4+2] = f.adress_a[i*4+2]-1 ;
		
		if(f.adress_a[i*4] >= 8 && f.adress_a[i*4] < FUNCTION_LENGTH-1-f.constNumber)
			ret.adress_a[(i)*4] = f.adress_a[i*4]+1 ;
		if(f.adress_a[i*4+1] >= 8 && f.adress_a[i*4+1] < FUNCTION_LENGTH-1-f.constNumber)
			ret.adress_a[(i)*4+1] = f.adress_a[i*4+1]+1 ;
		if(f.adress_a[i*4+2] >= 8 && f.adress_a[i*4+2] < FUNCTION_LENGTH-1-f.constNumber)
			ret.adress_a[(i)*4+2] = f.adress_a[i*4+2]+1 ;
		if(f.adress_a[i*4] < 8 )
			ret.adress_a[(i)*4] = f.adress_a[i*4] ;
		if(f.adress_a[i*4+1] < 8 )
			ret.adress_a[(i)*4+1] = f.adress_a[i*4+1] ;
		if(f.adress_a[i*4+2] < 8 )
			ret.adress_a[(i)*4+2] = f.adress_a[i*4+2] ;

	}
	ret.adress_a[(ret.byteCodeSize-1)*4+2] = 8 ;
	ret.adress_a[(ret.byteCodeSize-1)*4+1] = 8 ;
	ret.adress_a[(ret.byteCodeSize-1)*4] = FUNCTION_LENGTH-1 ;
	ret.byteCode[ret.byteCodeSize-1] = TOKEN_OPERATION_PLUS ;
	return ret ;
}

Function operator/(const double & a, const Function &f)
{
	Function ret ;
	ret.byteCodeSize = f.byteCodeSize+1 ;
	ret.constNumber = 1+f.constNumber ;
	ret.values[0] = a ;
	ret.use_temp = f.use_temp ;
	ret.byteCode = f.byteCode ;
	for(size_t i = 0 ; i < f.constNumber ; i++)
		ret.values[i+1] = f.values[i] ;
	
	for(size_t i = 0 ; i < f.byteCodeSize ; i++)
	{

		if(f.geo_op[i])
			ret.geo_op[i] = f.geo_op[i]->getCopy() ;
			
		if(f.adress_a[i*4] >= FUNCTION_LENGTH-1-f.constNumber)
			ret.adress_a[(i)*4] = f.adress_a[i*4]-1 ;
		
		if(f.adress_a[i*4+1] >= FUNCTION_LENGTH-1-f.constNumber)
			ret.adress_a[(i)*4+1] = f.adress_a[i*4+1]-1 ;
		
		if(f.adress_a[i*4+2] >= FUNCTION_LENGTH-1-f.constNumber)
			ret.adress_a[(i)*4+2] = f.adress_a[i*4+2]-1 ;
		
		if(f.adress_a[i*4] >= 8 && f.adress_a[i*4] < FUNCTION_LENGTH-1-f.constNumber)
			ret.adress_a[(i)*4] = f.adress_a[i*4]+1 ;
		if(f.adress_a[i*4+1] >= 8 && f.adress_a[i*4+1] < FUNCTION_LENGTH-1-f.constNumber)
			ret.adress_a[(i)*4+1] = f.adress_a[i*4+1]+1 ;
		if(f.adress_a[i*4+2] >= 8 && f.adress_a[i*4+2] < FUNCTION_LENGTH-1-f.constNumber)
			ret.adress_a[(i)*4+2] = f.adress_a[i*4+2]+1 ;
		if(f.adress_a[i*4] < 8 )
			ret.adress_a[(i)*4] = f.adress_a[i*4] ;
		if(f.adress_a[i*4+1] < 8 )
			ret.adress_a[(i)*4+1] = f.adress_a[i*4+1] ;
		if(f.adress_a[i*4+2] < 8 )
			ret.adress_a[(i)*4+2] = f.adress_a[i*4+2] ;

	}
	ret.adress_a[(ret.byteCodeSize-1)*4+2] = 8 ;
	ret.adress_a[(ret.byteCodeSize-1)*4+1] = 8 ;
	ret.adress_a[(ret.byteCodeSize-1)*4] = FUNCTION_LENGTH-1 ;
	ret.byteCode[ret.byteCodeSize-1] = TOKEN_OPERATION_DIVIDES ;
	return ret ;
	
}

Function Function::operator-(const Function &f) const 
{
	if(f.isNull())
		return *this ;
	Function ret ;
	concatenateFunctions(*this, f, ret);
	ret.byteCodeSize++ ;
	ret.adress_a[(ret.byteCodeSize-1)*4+2] = 8 ;
	ret.adress_a[(ret.byteCodeSize-1)*4+1] = 9 ;
	ret.adress_a[(ret.byteCodeSize-1)*4] = 8 ;
	ret.byteCode[ret.byteCodeSize-1] = TOKEN_OPERATION_MINUS ;
	return ret ;

}

Function Function::operator*(const double a) const
{
	if(std::abs(a) < POINT_TOLERANCE_2D)
	{
		return Function();
	}
	if(std::abs(a-1) < POINT_TOLERANCE_2D)
	{
		return *this;
	}
	Function ret ;
	ret.byteCodeSize = byteCodeSize+1 ;
	ret.adress_a = adress_a ;
	ret.use_temp = use_temp ;
	ret.byteCode = byteCode ;
	for(size_t i = 0 ; i < constNumber ; i++)
		ret.values[i] = values[i] ;
	
	ret.values[constNumber] = a ;

	for(size_t i = 0 ; i < byteCodeSize ; i++)
	{
		if(geo_op[i])
			ret.geo_op[i] = geo_op[i]->getCopy() ;
	}
	ret.constNumber = constNumber+1 ;
	ret.adress_a[(ret.byteCodeSize-1)*4+2] = 8 ;
	ret.adress_a[(ret.byteCodeSize-1)*4+1] = FUNCTION_LENGTH-1-constNumber ;
	ret.adress_a[(ret.byteCodeSize-1)*4] = 8 ;
	ret.byteCode[ret.byteCodeSize-1] = TOKEN_OPERATION_TIMES ;
	return ret ;
}

Function Function::operator/(const double a) const 
{
	if(std::abs(a) < POINT_TOLERANCE_2D)
	{
		std::cout << "Divide By Zero!" << std::endl ;
		exit(0) ;
		return Function() ;
	}
	if(std::abs(a-1) < POINT_TOLERANCE_2D)
	{
		return *this;
	}
	
	Function ret ;
	ret.byteCodeSize = byteCodeSize+1 ;
	ret.adress_a = adress_a ;
	ret.use_temp = use_temp ;
	ret.byteCode = byteCode ;
	for(size_t i = 0 ; i < constNumber ; i++)
		ret.values[i] = values[i] ;
	
	ret.values[constNumber] = a ;

	for(size_t i = 0 ; i < byteCodeSize ; i++)
	{

		if(geo_op[i])
			ret.geo_op[i] = geo_op[i]->getCopy() ;
	}
	ret.constNumber = constNumber+1 ;
	ret.adress_a[(ret.byteCodeSize-1)*4+2] = 8 ;
	ret.adress_a[(ret.byteCodeSize-1)*4+1] = FUNCTION_LENGTH-1-constNumber ;
	ret.adress_a[(ret.byteCodeSize-1)*4] = 8 ;
	ret.byteCode[ret.byteCodeSize-1] = TOKEN_OPERATION_DIVIDES ;
	return ret ;
}

Function Function::operator+(const double a) const
{
	if(std::abs(a) < POINT_TOLERANCE_2D)
	{
		return *this ;
	}
	Function ret ;
	ret.byteCodeSize = byteCodeSize+1 ;
	ret.adress_a = adress_a ;
	ret.use_temp = use_temp ;
	ret.byteCode = byteCode ;
	for(size_t i = 0 ; i < constNumber ; i++)
		ret.values[i] = values[i] ;
	
	ret.values[constNumber] = a ;

	for(size_t i = 0 ; i < byteCodeSize ; i++)
	{
		if(geo_op[i])
			ret.geo_op[i] = geo_op[i]->getCopy() ;
	}
	ret.constNumber = constNumber+1 ;
	ret.adress_a[(ret.byteCodeSize-1)*4+2] = 8 ;
	ret.adress_a[(ret.byteCodeSize-1)*4+1] = FUNCTION_LENGTH-1-constNumber ;
	ret.adress_a[(ret.byteCodeSize-1)*4] = 8 ;
	ret.byteCode[ret.byteCodeSize-1] = TOKEN_OPERATION_PLUS ;
	return ret ;
}

Function Function::operator-(const double a) const 
{
	if(std::abs(a) < POINT_TOLERANCE_2D)
	{
		return *this;
	}
	
	Function ret ;
	ret.byteCodeSize = byteCodeSize+1 ;
	ret.adress_a = adress_a ;
	ret.use_temp = use_temp ;
	ret.byteCode = byteCode ;
	ret.values = values ;
	
	ret.values[constNumber] = a ;

	for(size_t i = 0 ; i < byteCodeSize ; i++)
	{
		if(geo_op[i])
			ret.geo_op[i] = geo_op[i]->getCopy() ;
	}
	ret.constNumber = constNumber+1 ;
	ret.adress_a[(ret.byteCodeSize-1)*4+2] = 8 ;
	ret.adress_a[(ret.byteCodeSize-1)*4+1] = FUNCTION_LENGTH-1-constNumber ;
	ret.adress_a[(ret.byteCodeSize-1)*4] = 8 ;
	ret.byteCode[ret.byteCodeSize-1] = TOKEN_OPERATION_MINUS ;
	return ret ;
}

Function  Function::operator^(const int a) const
{
	Function ret ;
	ret.byteCodeSize = byteCodeSize+1 ;
	ret.adress_a = adress_a ;
	ret.use_temp = use_temp ;
	ret.byteCode = byteCode ;
	for(size_t i = 0 ; i < constNumber ; i++)
		ret.values[i] = values[i] ;
	
	ret.values[constNumber] = a ;

	for(size_t i = 0 ; i < byteCodeSize ; i++)
	{
		if(geo_op[i])
			ret.geo_op[i] = geo_op[i]->getCopy() ;
	}
	ret.constNumber = constNumber+1 ;
	ret.adress_a[(ret.byteCodeSize-1)*4+2] = 8 ;
	ret.adress_a[(ret.byteCodeSize-1)*4+1] = FUNCTION_LENGTH-1-constNumber ;
	ret.adress_a[(ret.byteCodeSize-1)*4] = 8 ;
	ret.byteCode[ret.byteCodeSize-1] = TOKEN_OPERATION_POWER ;
	return ret ;
	return ret ;
}

void Function::operator*=(const Function &f) 
{
	if(isNull())
	{
		return ;
	}
	if(f.isNull())
	{
		defaultInitialise();
		return ;
	}
		
	concatenateFunctions(*this, f, *this);
	byteCodeSize++ ;
	adress_a[(byteCodeSize-1)*4+2] = 8 ;
	adress_a[(byteCodeSize-1)*4+1] = 9 ;
	adress_a[(byteCodeSize-1)*4] = 8 ;
	byteCode[byteCodeSize-1] = TOKEN_OPERATION_TIMES ;
}


void Function::operator/=(const Function &f)  
{
	if(isNull())
	{
		std::cout << "Divide By Zero" << std::endl ;
		exit(0) ;
		return ;
	}
	concatenateFunctions(*this, f, *this);
	byteCodeSize++ ;
	adress_a[(byteCodeSize-1)*4+2] = 8 ;
	adress_a[(byteCodeSize-1)*4+1] = 9 ;
	adress_a[(byteCodeSize-1)*4] = 8 ;
	byteCode[byteCodeSize-1] = TOKEN_OPERATION_DIVIDES ;
}

void Function::operator+=(const Function &f) 
{
	if(isNull())
	{
		*this = f ;
		return ;
	}
	if(f.isNull())
		return ;

	concatenateFunctions(*this, f, *this);
	byteCodeSize++ ;
	adress_a[(byteCodeSize-1)*4+2] = 8 ;
	adress_a[(byteCodeSize-1)*4+1] = 9 ;
	adress_a[(byteCodeSize-1)*4] = 8 ;
	byteCode[byteCodeSize-1] = TOKEN_OPERATION_PLUS ;
}

void Function::operator-=(const Function &f)  
{
	if(isNull())
	{
		*this = f*-1 ;
		return ;
	}
	concatenateFunctions(*this, f, *this);
	byteCodeSize++ ;
	adress_a[(byteCodeSize-1)*4+2] = 8 ;
	adress_a[(byteCodeSize-1)*4+1] = 9 ;
	adress_a[(byteCodeSize-1)*4] = 8 ;
	byteCode[byteCodeSize-1] = TOKEN_OPERATION_MINUS ;
}

void Function::operator*=(const double a) 
{
	if(std::abs(a) < POINT_TOLERANCE_2D)
	{
		defaultInitialise();
		return ;
	}
	byteCodeSize++ ;
	values[constNumber] = a ;
	constNumber++;
	adress_a[(byteCodeSize-1)*4+2] = 8 ;
	adress_a[(byteCodeSize-1)*4+1] = FUNCTION_LENGTH-1-constNumber+1 ;
	adress_a[(byteCodeSize-1)*4] = 8 ;
	byteCode[byteCodeSize-1] = TOKEN_OPERATION_TIMES ;
}

void Function::operator/=(const double a)  
{
	if(std::abs(a-1) < POINT_TOLERANCE_2D)
		return ;
	if(std::abs(a) < POINT_TOLERANCE_2D)
	{
		std::cout << "DivideBy Zero" << std::endl ;
		exit(0) ;
	}
	byteCodeSize++ ;
	values[constNumber] = a ;
	constNumber++;
	adress_a[(byteCodeSize-1)*4+2] = 8 ;
	adress_a[(byteCodeSize-1)*4+1] = FUNCTION_LENGTH-1-constNumber+1 ;
	adress_a[(byteCodeSize-1)*4] = 8 ;
	byteCode[byteCodeSize-1] = TOKEN_OPERATION_DIVIDES ;
}

void Function::operator+=(const double a) 
{
	if(std::abs(a) < POINT_TOLERANCE_2D)
		return ;
	byteCodeSize++ ;
	values[constNumber] = a ;
	constNumber++;
	adress_a[(byteCodeSize-1)*4+2] = 8 ;
	adress_a[(byteCodeSize-1)*4+1] = FUNCTION_LENGTH-1-constNumber+1 ;
	adress_a[(byteCodeSize-1)*4] = 8 ;
	byteCode[byteCodeSize-1] = TOKEN_OPERATION_PLUS ;
}

void Function::operator-=(const double a)  
{
	if(std::abs(a) < POINT_TOLERANCE_2D)
		return ;
	byteCodeSize++ ;
	values[constNumber] = a ;
	constNumber++;
	adress_a[(byteCodeSize-1)*4+2] = 8 ;
	adress_a[(byteCodeSize-1)*4+1] = FUNCTION_LENGTH-1-constNumber+1 ;
	adress_a[(byteCodeSize-1)*4] = 8 ;
	byteCode[byteCodeSize-1] = TOKEN_OPERATION_MINUS ;
}

const Function & Function::d(const Variable v) const
{
	if(derivative &&  derivative->size() >= v)
	{
		return *(*derivative)[v] ;
	}
}

Function & Function::d(const Variable v) 
{
	if(derivative &&  derivative->size() >= v)
	{
		return *(*derivative)[v] ;
	}
	
}

std::valarray<Function *> & Function::getDerivatives() const
{
	return *derivative ;
}

std::valarray<Function *> & Function::getDerivatives()
{
	return *derivative ;
}


Function f_exp(const Function &f)
{
	Function ret = f ;
	ret.byteCodeSize++ ;
	ret.adress_a[(ret.byteCodeSize-1)*4+2] = 8 ;
	ret.adress_a[(ret.byteCodeSize-1)*4+1] = 9 ;
	ret.adress_a[(ret.byteCodeSize-1)*4] = 8 ;
	ret.byteCode[ret.byteCodeSize-1] = TOKEN_OPERATION_EXP ;
	return ret ;
}

Function f_abs(const Function &f)
{
	Function ret = f ;
	ret.byteCodeSize++ ;
	ret.adress_a[(ret.byteCodeSize-1)*4+2] = 8 ;
	ret.adress_a[(ret.byteCodeSize-1)*4+1] = 9 ;
	ret.adress_a[(ret.byteCodeSize-1)*4] = 8 ;
	ret.byteCode[ret.byteCodeSize-1] = TOKEN_OPERATION_ABS ;
	return ret ;
}

Function f_log(const Function &f)
{
	Function ret = f ;
	ret.byteCodeSize++ ;
	ret.adress_a[(ret.byteCodeSize-1)*4+2] = 8 ;
	ret.adress_a[(ret.byteCodeSize-1)*4+1] = 9 ;
	ret.adress_a[(ret.byteCodeSize-1)*4] = 8 ;
	ret.byteCode[ret.byteCodeSize-1] = TOKEN_OPERATION_LOG ;
	return ret ;
}

Function f_sqrt(const Function &f)
{
	Function ret = f ;
	ret.byteCodeSize++ ;
	ret.adress_a[(ret.byteCodeSize-1)*4+2] = 8 ;
	ret.adress_a[(ret.byteCodeSize-1)*4+1] = 9 ;
	ret.adress_a[(ret.byteCodeSize-1)*4] = 8 ;
	ret.byteCode[ret.byteCodeSize-1] = TOKEN_OPERATION_SQRT ;
	return ret ;
}

Function f_atan2(const Function &f0, const Function &f1)
{
	Function ret ;
	concatenateFunctions(f0, f1, ret);
	ret.byteCodeSize++ ;
	ret.adress_a[(ret.byteCodeSize-1)*4+2] = 8 ;
	ret.adress_a[(ret.byteCodeSize-1)*4+1] = 9 ;
	ret.adress_a[(ret.byteCodeSize-1)*4] = 8 ;
	ret.byteCode[ret.byteCodeSize-1] = TOKEN_OPERATION_ATAN2 ;
	return ret ;
}

Function f_sin(const Function &f)
{
	Function ret = f ;
	ret.byteCodeSize++ ;
	ret.adress_a[(ret.byteCodeSize-1)*4+2] = 8 ;
	ret.adress_a[(ret.byteCodeSize-1)*4+1] = 9 ;
	ret.adress_a[(ret.byteCodeSize-1)*4] = 8 ;
	ret.byteCode[ret.byteCodeSize-1] = TOKEN_OPERATION_SIN ;
	return ret ;
}


Mu::Function f_sign(const Mu::Function &f)
{
	Function ret = f ;
	ret.byteCodeSize++ ;
	ret.adress_a[(ret.byteCodeSize-1)*4+2] = 8 ;
	ret.adress_a[(ret.byteCodeSize-1)*4+1] = 9 ;
	ret.adress_a[(ret.byteCodeSize-1)*4] = 8 ;
	ret.byteCode[ret.byteCodeSize-1] = TOKEN_OPERATION_SIGN ;
	return ret ;
}



Mu::Function f_positivity(const Mu::Function &f)
{
	Function ret = f ;
	ret.byteCodeSize++ ;
	ret.adress_a[(ret.byteCodeSize-1)*4+2] = 8 ;
	ret.adress_a[(ret.byteCodeSize-1)*4+1] = 9 ;
	ret.adress_a[(ret.byteCodeSize-1)*4] = 8 ;
	ret.byteCode[ret.byteCodeSize-1] = TOKEN_OPERATION_POSITIVITY ;
	return ret ;
}


Mu::Function f_negativity(const Mu::Function &f)
{
	Function ret = f ;
	ret.byteCodeSize++ ;
	ret.adress_a[(ret.byteCodeSize-1)*4+2] = 8 ;
	ret.adress_a[(ret.byteCodeSize-1)*4+1] = 9 ;
	ret.adress_a[(ret.byteCodeSize-1)*4] = 8 ;
	ret.byteCode[ret.byteCodeSize-1] = TOKEN_OPERATION_NEGATIVITY ;
	return ret ;
}


Function f_cos(const Function &f)
{
	Function ret = f ;
	ret.byteCodeSize++ ;
	ret.adress_a[(ret.byteCodeSize-1)*4+2] = 8 ;
	ret.adress_a[(ret.byteCodeSize-1)*4+1] = 9 ;
	ret.adress_a[(ret.byteCodeSize-1)*4] = 8 ;
	ret.byteCode[ret.byteCodeSize-1] = TOKEN_OPERATION_COS ;
	return ret ;
} 

const std::vector< Point > & Function::getIntegrationHint() const
{
	return iPoint ;
}

const Point & Function::getIntegrationHint(size_t i) const
{
	return iPoint[i] ;
}

void Function::setIntegrationHint(const std::vector< Point > v)
{
	iPoint.clear() ;
	iPoint = v ;
}

void Function::addIntegrationHint(const Point  p)
{
	iPoint.push_back(p) ;
}

bool Function::hasIntegrationHint() const
{
	return iPoint.size() > 0 ;
}





bool Function::isBinaryOperator(TokenOperationType t) const
{
	return ( t== TOKEN_OPERATION_PLUS || 
	         t== TOKEN_OPERATION_MINUS || 
	         t== TOKEN_OPERATION_TIMES || 
	         t== TOKEN_OPERATION_DIVIDES || 
	         t== TOKEN_OPERATION_POWER ||
	         t== TOKEN_OPERATION_INTERPOLATE ||
	         t== TOKEN_OPERATION_ATAN2) ;
}

bool Function::isBinaryVectorisableOperator(TokenOperationType t) const
{
	return ( t == TOKEN_OPERATION_PLUS || 
	         t == TOKEN_OPERATION_MINUS || 
	         t == TOKEN_OPERATION_DIVIDES || 
	         t == TOKEN_OPERATION_POWER ) ;
}

bool Function::isUnaryOperator(TokenOperationType t) const
{
	return ( t == TOKEN_OPERATION_SIN || 
	         t == TOKEN_OPERATION_COS || 
	         t == TOKEN_OPERATION_SQRT || 
	         t == TOKEN_OPERATION_COSH || 
	         t == TOKEN_OPERATION_SINH ||
	         t == TOKEN_OPERATION_TANH || 
	         t == TOKEN_OPERATION_EXP || 
	         t == TOKEN_OPERATION_LOG || 
	         t == TOKEN_OPERATION_TAN  
	       ) ;
}
