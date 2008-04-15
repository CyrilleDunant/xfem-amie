//
// C++ Implementation: vm_function_base
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "vm_function_base.h"

using namespace Mu ;

Function::Function(const ByteCode & b_0, const ByteCode & b_1, RefCountedToken op, const bool diff) : iPoint(0), byteCode(b_0.size() + b_1.size()+1) , e_diff(diff)
{
	this->dofID =-1 ;
	this->ptID =-1 ;
	
	std::copy(&b_0[0], &b_0[b_0.size()], &byteCode[0]) ;
	std::copy(&b_1[0], &b_1[b_1.size()], &byteCode[b_0.size()]) ;
	byteCode[b_0.size() + b_1.size()] = op ;
}

Function::Function(const ByteCode & b_0, const ByteCode & b_1) : iPoint(0), byteCode(b_0.size() + b_1.size()) , e_diff(false)
{
	this->dofID =-1 ;
	this->ptID =-1 ;
	
	std::copy(&b_0[0], &b_0[b_0.size()], &byteCode[0]) ;
	std::copy(&b_1[0], &b_1[b_1.size()], &byteCode[b_0.size()]) ;
}

Function::Function(const ByteCode &b_0, const double a, RefCountedToken op, const bool diff) :  iPoint(0), byteCode(b_0.size() + 2) , e_diff(diff)
{
	this->dofID =-1 ;
	this->ptID =-1 ;
	
	std::copy(&b_0[0], &b_0[b_0.size()], &byteCode[0]) ;
	byteCode[b_0.size()] = RefCountedToken(new ConstantToken(a)) ;
	byteCode[b_0.size()+1] = op ;
}

Function::Function(const double a, const ByteCode &b_0,  RefCountedToken op, const bool diff) :  iPoint(0), byteCode(b_0.size() + 2) , e_diff(diff)
{
	this->dofID =-1 ;
	this->ptID =-1 ;
	
	byteCode[0] = RefCountedToken(new ConstantToken(a)) ;
	std::copy(&b_0[0], &b_0[b_0.size()], &byteCode[1]) ;
	byteCode[b_0.size()+1] = op ;
}
	
Function Function::operator*(const Geometry *f) const
{
	ByteCode bc(1) ;
	bc[0]= new DomainToken(f) ;
	return Function(byteCode, bc, new TimesOperatorToken()) ;
}

void Function::operator*=(const Geometry *f) 
{
	ByteCode bc(1) ;
	bc[0]= new DomainToken(f) ;
	*this =  Function(byteCode, bc, new TimesOperatorToken()) ;
}


Function & Function::operator=(const Function &f)
{
	this->byteCode.resize(f.getByteCode().size()) ;
	for(size_t i = 0 ; i < f.getByteCode().size() ; i++)
	{
		this->byteCode[i] = f.getToken(i) ;
	}
	
	if(f.e_diff)
	{
		this->derivative.resize(f.derivative.size()) ;
		this->derivative = f.derivative ;
		this->e_diff = true ;
	}
	else
	{
		this->e_diff = false ;
	}
	
	return *this ;
}

void  Function::preCalculate(const GaussPointArray & gp , std::vector<Variable> & var, const double eps)
{
	//this needs to be in to stages, otherwise memory gets accessed too early.
	Vector * newVal = new Vector(VirtualMachine().eval(*this, gp)) ;
	precalc[gp.id] = newVal ;
	std::map<Variable, Vector *> val ;
	for(size_t i = 0 ; i < var.size() ; i++)
	{
		val[var[i]] = new Vector(VirtualMachine().deval(*this,var[i] ,gp, eps)) ;
	}
	dprecalc[gp.id] = val ;
}

void  Function::preCalculate(const GaussPointArray & gp )
{
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
	return dprecalc.find(gp.id) != dprecalc.end() && dprecalc.find(gp.id)->second.find(v) != dprecalc.find(gp.id)->second.end() ;
}

int Function::getDofID() const
{
	return this->dofID ;
}
void Function::setDofID(size_t id)
{
	this->dofID = id ;
}

int Function::getPointID() const
{
	return this->ptID ;
}
void Function::setPointID(size_t id)
{
	this->ptID = id ;
}

const RefCountedToken Function::toToken(const std::string str) const
{
	if(isNumeral(str[0]) || (str[0] == '-' && isNumeral(str[1])))
	{
		return new ConstantToken(atof(str.c_str())) ;
	}
	else if(str == std::string("sin"))
	{
		return new SinToken() ;
	}
	else if(str.c_str() == std::string("cos"))
	{
		return new CosToken() ;
	}
	else if(str.c_str() == std::string("tan"))
	{
		return new TanToken() ;
	}
	else if(str.c_str() == std::string("sinh"))
	{
		return new SinhToken() ;
	}
	else if(str.c_str() == std::string("cosh"))
	{
		return new CoshToken() ;
	}
	else if(str.c_str() == std::string("tanh"))
	{
		return new TanhToken() ;
	}
	else if(str.c_str() == std::string("exp"))
	{
		return new ExpToken() ;
	}
	else if(str.c_str() == std::string("log"))
	{
		return new LogToken() ;
	}
	else if(str.c_str() == std::string("sqrt"))
	{
		return new SqrtToken() ;
	}
	else if(str.c_str() == std::string("sign"))
	{
		return new SignFunctionToken() ;
	}
	else if(str.c_str() == std::string("atan2"))
	{
		return new Atan2Token() ;
	}
	else if(str.c_str() == std::string("+"))
	{
		return new PlusOperatorToken() ;
	}
	else if(str.c_str() == std::string("-"))
	{
		return new MinusOperatorToken() ;
	}
	else if(str.c_str() == std::string("*"))
	{
		return new TimesOperatorToken() ;
	}
	else if(str.c_str() == std::string("/"))
	{
		return new DivideOperatorToken() ;
	}
	else if(str.c_str() == std::string("^"))
	{
		return new PowerOperatorToken() ;
	}
	else if(str.c_str() == std::string("x"))
	{
		return new XToken() ;
	}
	else if(str.c_str() == std::string("y"))
	{
		return new YToken();
	}
	else if(str.c_str() == std::string("z"))
	{
		return new ZToken() ;
	}
	else if(str.c_str() == std::string("t"))
	{
		return new TToken() ;
	}
	else if(str.c_str() == std::string("u"))
	{
		return new UToken();
	}
	else if(str.c_str() == std::string("v"))
	{
		return new VToken() ;
	}
	else if(str.c_str() == std::string("w"))
	{
		return new WToken() ;
	}
	else
	{
		return NULL ;
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

std::pair<size_t, RefCountedToken > Function::getNext(size_t init, const std::string & form)
{
	return getNext(init, form.c_str()) ;
}

std::pair<size_t, RefCountedToken > Function::getNext(size_t init, const char * form)
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
	

const RefCountedToken& Function::getToken(const size_t i) const
{
	return byteCode[i] ;
}

Function::Function() :  derivative(0), e_diff(false)
{
	this->dofID =-1 ;
	this->ptID =-1 ;
	this->byteCode.resize(1) ;
	this->byteCode[0] = RefCountedToken(new NullToken()) ;
	compiled = false ;
}

Function::Function(const char *f) :  derivative(0),iPoint(0), e_diff(false)
{	
	this->dofID =-1 ;
	this->ptID =-1 ;
	
	size_t init = 0 ;
	
	std::vector<RefCountedToken> tempb ;
	
	while(init < strlen(f))
	{
		std::pair<size_t, RefCountedToken> temp = getNext(init, f) ;
		if((const Token *)temp.second != NULL)
			tempb.push_back(temp.second) ;
		init = temp.first ;
	}
	
	byteCode.resize(tempb.size()) ;
	std::copy(tempb.begin(), tempb.end(), &byteCode[0]) ;
	compiled = false ;
}

Function::Function(const std::string &f) : derivative(0), iPoint(0), e_diff(false)
{	
	this->dofID =-1 ;
	this->ptID =-1 ;
	
	size_t init = 0 ;
	
	std::vector<RefCountedToken> tempb ;
	
	while(init < f.size())
	{
		std::pair<size_t, RefCountedToken> temp = getNext(init, f) ;
		if((const Token *)temp.second != NULL)
			tempb.push_back(temp.second) ;
		init = temp.first ;
	}
	
	byteCode.resize(tempb.size()) ;
	std::copy(tempb.begin(), tempb.end(), &byteCode[0]) ;
	compiled = false ;
}

Function::Function(const std::valarray<Matrix> & coeffs, bool diff) :iPoint(0) , e_diff(diff)
{
	this->dofID =-1 ;
	this->ptID =-1 ;
	
	std::vector<RefCountedToken> tempb ;
	
	bool first = true ;
	
	for(size_t i =  0;  i < coeffs.size() ; i++)
	{
		for(size_t j =  0;  j < coeffs[i].numRows() ; j++)
		{
			for(size_t k = 0 ; k < coeffs[i].numCols() ;  k++)
			{

					if(coeffs[i][j][k] != 0)
					{
						tempb.push_back(new ConstantToken(coeffs[i][j][k])) ;
						
						if(i > 0)
						{
							tempb.push_back(new XToken()) ;
							if(i > 1)
							{
								tempb.push_back(new ConstantToken(i)) ;
								tempb.push_back(new PowerOperatorToken()) ;
							}
							tempb.push_back(new TimesOperatorToken()) ;
						}
						
						if(j > 0)
						{
							tempb.push_back(new YToken()) ;
							if(j > 1)
							{
								tempb.push_back(new ConstantToken(j)) ;
								tempb.push_back(new PowerOperatorToken()) ;
							}
							tempb.push_back(new TimesOperatorToken()) ;
						}
						
						if(k > 0)
						{
							tempb.push_back(new ZToken()) ;
							if(k > 1)
							{
								tempb.push_back(new ConstantToken(k)) ;
								tempb.push_back(new PowerOperatorToken()) ;
							}
							tempb.push_back(new TimesOperatorToken()) ;
						}
						
						if(!first)
							tempb.push_back(new PlusOperatorToken()) ;
						else
							first = false ;
						
				}
			}
		}
	}
	
	byteCode.resize(tempb.size()) ;
	std::copy(tempb.begin(), tempb.end(), &byteCode[0]) ;
	if(first)
	{
		byteCode.resize(1) ;
		byteCode[0] = new NullToken() ;
	}
	
	if(e_diff)
	{
		derivative.resize(3) ;
		std::valarray<Matrix> dx(coeffs.size()-1) ;
		for(size_t i = 0 ; i < dx.size() ; i++)
		{
			dx[i] = coeffs[i+1]*(i+1) ;
		}
		
		derivative[XI] = Function(dx, false) ;
		std::valarray<Matrix> dy(coeffs.size()) ;
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
		derivative[ETA] = Function(dy, false) ;
		
		std::valarray<Matrix> dz(coeffs.size()) ;
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
		derivative[ZETA] = Function(dz, false) ;
	}
	compiled = false ;
}


Function::Function(const std::valarray<std::valarray<Matrix> > & coeffs, bool diff) :iPoint(0) , e_diff(diff)
{
	this->dofID =-1 ;
	this->ptID =-1 ;
	
	std::vector<RefCountedToken> tempb ;
	
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
						tempb.push_back(new ConstantToken(coeffs[i][j][k][l])) ;
						
						if(i > 0)
						{
							tempb.push_back(new XToken()) ;
							if(i > 1)
							{
								tempb.push_back(new ConstantToken(i)) ;
								tempb.push_back(new PowerOperatorToken()) ;
							}
							tempb.push_back(new TimesOperatorToken()) ;
						}
						
						if(j > 0)
						{
							tempb.push_back(new YToken()) ;
							if(j > 1)
							{
								tempb.push_back(new ConstantToken(j)) ;
								tempb.push_back(new PowerOperatorToken()) ;
							}
							tempb.push_back(new TimesOperatorToken()) ;
						}
						
						if(k > 0)
						{
							tempb.push_back(new ZToken()) ;
							if(k > 1)
							{
								tempb.push_back(new ConstantToken(k)) ;
								tempb.push_back(new PowerOperatorToken()) ;
							}
							tempb.push_back(new TimesOperatorToken()) ;
						}
						
						if(l > 0)
						{
							tempb.push_back(new TToken()) ;
							if(l > 1)
							{
								tempb.push_back(new ConstantToken(l)) ;
								tempb.push_back(new PowerOperatorToken()) ;
							}
							tempb.push_back(new TimesOperatorToken()) ;
						}
						
						if(!first)
							tempb.push_back(new PlusOperatorToken()) ;
						else
							first = false ;
	
					}
				}
			}
		}
	}
	
	byteCode.resize(tempb.size()) ;
	std::copy(tempb.begin(), tempb.end(), &byteCode[0]) ;
	if(first)
	{
		byteCode.resize(1) ;
		byteCode[0] = new NullToken() ;
	}
	
	if(e_diff)
	{
		derivative.resize(4) ;
		std::valarray< std::valarray<Matrix> > dx(coeffs.size()-1) ;
		for(size_t i = 0 ; i < dx.size() ; i++)
		{
			std::valarray<Matrix> dx_(coeffs.size()) ;
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
		derivative[XI] = Function(dx, false) ;
		
		std::valarray< std::valarray<Matrix> > dy(coeffs.size()) ;
		for(size_t i = 0 ; i < dy.size() ; i++)
		{
			std::valarray<Matrix> dy_(coeffs[i].size()-1) ;
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
		derivative[ETA] = Function(dy, false) ;
		
		std::valarray< std::valarray<Matrix> > dz(coeffs.size()) ;
		for(size_t i = 0 ; i < dz.size() ; i++)
		{
			std::valarray<Matrix> dz_(coeffs[i].size()) ;
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
		derivative[ZETA] = Function(dz, false) ;
		
		std::valarray< std::valarray<Matrix> > dt(coeffs.size()) ;
		for(size_t i = 0 ; i < dt.size() ; i++)
		{
			std::valarray<Matrix> dt_(coeffs[i].size()) ;
			
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
		derivative[TIME_VARIABLE] = Function(dt, false) ;
	}
	compiled = false ;
}

Function::Function(const Matrix & coeffs, bool diff) :iPoint(0) , e_diff(diff)
{
	this->dofID =-1 ;
	this->ptID =-1 ;
	compiled = false ;
	
	std::vector<RefCountedToken> tempb ;
	
	bool first = true ;
	
	for(size_t j = 0 ; j < coeffs.numRows() ;  j++)
	{
		for(size_t k = 0 ; k < coeffs.numCols() ;  k++)
		{
			if(coeffs[j][k] != 0)
			{
				tempb.push_back(new ConstantToken(coeffs[j][k])) ;
				
				
				if(j > 0)
				{
					tempb.push_back(new XToken()) ;
					if(j > 1)
					{
						tempb.push_back(new ConstantToken(j)) ;
						tempb.push_back(new PowerOperatorToken()) ;
					}
					tempb.push_back(new TimesOperatorToken()) ;
				}
				
				if(k > 0)
				{
					tempb.push_back(new YToken()) ;
					if(k > 1)
					{
						tempb.push_back(new ConstantToken(k)) ;
						tempb.push_back(new PowerOperatorToken()) ;
					}
					tempb.push_back(new TimesOperatorToken()) ;
				}
				
				if(!first)
					tempb.push_back(new PlusOperatorToken()) ;
				else
					first = false ;
			}
		}
	}
	
	byteCode.resize(tempb.size()) ;
	std::copy(tempb.begin(), tempb.end(), &byteCode[0]) ;
	if(first)
	{
		byteCode.resize(1) ;
		byteCode[0] = new NullToken() ;
	}
	
	if(e_diff)
	{
		derivative.resize(2) ;
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
		
		derivative[XI] = Function(dx, false) ;
		
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
		
		derivative[ETA] = Function(dy, false) ;
	}
}

bool Function::isNull() const
{
	return ((this->byteCode.size() == 1 )&& this->byteCode[0]->isNull) ;
}

Function::Function(const std::valarray<double> & coeffs, bool diff) :iPoint(0) , e_diff(diff)
{
	this->dofID =-1 ;
	this->ptID =-1 ;
	compiled = false ;
	
	std::vector<RefCountedToken> tempb ;
	
	for(size_t j = 0 ; j < coeffs.size() ;  j++)
	{

			if(coeffs[j] != 0)
			{
				tempb.push_back(new ConstantToken(coeffs[j])) ;
				
				
				if(j > 0)
				{
					tempb.push_back(new XToken()) ;
					if(j > 1)
					{
						tempb.push_back(new ConstantToken(j)) ;
						tempb.push_back(new PowerOperatorToken()) ;
					}
					tempb.push_back(new TimesOperatorToken()) ;
				}
				
				
				if(tempb.size() > 1)
					tempb.push_back(new PlusOperatorToken()) ;
			}

	}
	
	byteCode.resize(tempb.size()) ;
	std::copy(tempb.begin(), tempb.end(), &byteCode[0]) ;
	
	if(e_diff)
	{
		derivative.resize(1) ;
		size_t ds = 0 ;
		if((int)coeffs.size()-1 > 0)
			ds = coeffs.size()-1 ;
		std::valarray<double> dx(ds) ;
// 		std::cout << "coeffs.size()-1 = "<< coeffs.size()-1 << std::endl ;
		for(size_t i = 1 ; i < dx.size()+1 ; i++)
		{
			dx[i] = coeffs[i]*(i) ;
		}
		
		derivative[XI] = Function(dx, false) ;
	}
}

Function::Function(const Segment s, const Function & x, const Function & y, PositionTokenType t) :  derivative(2)  , byteCode(x.getByteCode().size()+y.getByteCode().size()+1), e_diff(false)
{
	compiled = false ;
	switch(t)
	{
	case POSITION_TOKEN :
		{
			this->dofID =-1 ;
			this->ptID =-1 ;
			
			for(size_t i = 0 ; i < x.getByteCode().size() ; i++)
				byteCode[i] = x.getByteCode()[i] ;
			for(size_t i = 0 ; i < y.getByteCode().size() ; i++)
				byteCode[i+x.getByteCode().size()] = y.getByteCode()[i] ;
			
			byteCode[byteCode.size()-1] = RefCountedToken(new PositionOperatorToken(s)) ;
			
			derivative[XI] = Function() ;
			derivative[ETA] = Function() ;
			
			break ;
		}
	case PROJECTION_TOKEN :
		{
			derivative.resize(0) ;
			this->dofID =-1 ;
			this->ptID =-1 ;
			
			byteCode[0] = RefCountedToken(new ProjectionToken(s)) ;
			
			break ;
		}
	}
}

Function::Function(Geometry * g) :  derivative(2)  , byteCode(1), e_diff(true)
{

			this->dofID =-1 ;
			this->ptID =-1 ;
			
			byteCode[0] = RefCountedToken(new DomainBinaryOperatorToken(g)) ;
			
			derivative[XI] = Function() ;
			derivative[ETA] = Function() ;
	
}

Function Function::operator()(const Function & f) const
{
	return Function(f.getByteCode(), getByteCode()) ;
}

Function::Function(const Line & l, Function x, Function y) : derivative(2), byteCode(x.getByteCode().size()+y.getByteCode().size()+1), e_diff(false)
{
	this->dofID =-1 ;
	this->ptID =-1 ;
	for(size_t i = 0 ; i < x.getByteCode().size() ; i++)
		byteCode[i] = x.getByteCode()[i] ;
	for(size_t i = 0 ; i < y.getByteCode().size() ; i++)
		byteCode[i+x.getByteCode().size()] = y.getByteCode()[i] ;
	byteCode[byteCode.size()-1] = RefCountedToken(new LineDistanceOperatorToken(l)) ;
}


Function::Function(const Point & l, Function x, Function y) : derivative(2), byteCode(x.getByteCode().size()+y.getByteCode().size()+1), e_diff(false)
{
	this->dofID =-1 ;
	this->ptID =-1 ;
	for(size_t i = 0 ; i < y.getByteCode().size() ; i++)
		byteCode[i] = y.getByteCode()[i] ;
	for(size_t i = 0 ; i < x.getByteCode().size() ; i++)
		byteCode[i+y.getByteCode().size()] = x.getByteCode()[i] ;
	byteCode[byteCode.size()-1] = RefCountedToken(new PointDistanceBinaryOperatorToken(l)) ;
}

Function::Function(double a, Function x, Function y) : derivative(2), byteCode(x.getByteCode().size()+y.getByteCode().size()+1), e_diff(false)
{
	this->dofID =-1 ;
	this->ptID =-1 ;
	for(size_t i = 0 ; i < y.getByteCode().size() ; i++)
		byteCode[i] = y.getByteCode()[i] ;
	for(size_t i = 0 ; i < x.getByteCode().size() ; i++)
		byteCode[i+y.getByteCode().size()] = x.getByteCode()[i] ;
	byteCode[byteCode.size()-1] = RefCountedToken(new RotationBinaryOperatorToken(a)) ;
}

Function::Function(double a,const Point & p,  Function x, Function y): derivative(2), byteCode(x.getByteCode().size()+y.getByteCode().size()+1), e_diff(false)
{
	this->dofID =-1 ;
	this->ptID =-1 ;
	for(size_t i = 0 ; i < y.getByteCode().size() ; i++)
		byteCode[i] = y.getByteCode()[i] ;
	for(size_t i = 0 ; i < x.getByteCode().size() ; i++)
		byteCode[i+y.getByteCode().size()] = x.getByteCode()[i] ;
	byteCode[byteCode.size()-1] = RefCountedToken(new AngleBinaryOperatorToken(a,p)) ;
}

Function::Function( Geometry * g, Function x, Function y) : derivative(2), byteCode(x.getByteCode().size()+y.getByteCode().size()+1), e_diff(false)
{
	this->dofID =-1 ;
	this->ptID =-1 ;
	for(size_t i = 0 ; i < y.getByteCode().size() ; i++)
		byteCode[i] = y.getByteCode()[i] ;
	for(size_t i = 0 ; i < x.getByteCode().size() ; i++)
		byteCode[i+y.getByteCode().size()] = x.getByteCode()[i] ;
	byteCode[byteCode.size()-1] = RefCountedToken(new DomainBinaryOperatorToken(g)) ;
}

Function::Function(const std::vector<Segment> s , const Function & x, const Function & y, PositionTokenType t) :  derivative(2) , byteCode(x.getByteCode().size()+y.getByteCode().size()+1), e_diff(false)
{
	switch(t)
	{
	case POSITION_TOKEN :
		{
			this->dofID =-1 ;
			this->ptID =-1 ;
			for(size_t i = 0 ; i < x.getByteCode().size() ; i++)
				byteCode[i] = x.getByteCode()[i] ;
			for(size_t i = 0 ; i < y.getByteCode().size() ; i++)
				byteCode[i+x.getByteCode().size()] = y.getByteCode()[i] ;
			
			byteCode[byteCode.size()-1] = RefCountedToken(new PositionOperatorToken(s)) ;
			
			derivative[XI] = Function() ;
			derivative[ETA] = Function() ;
			
			break ;
		}
	case PROJECTION_TOKEN :
		{
			derivative.resize(0) ;
			this->dofID =-1 ;
			this->ptID =-1 ;
			byteCode.resize(s.size()+s.size()-1) ;
			for(size_t i = 0 ; i < s.size() ; i++)
			{
				byteCode[i] = RefCountedToken(new ProjectionToken(s[i])) ;
			}
			for(size_t i = s.size() ; i < byteCode.size() ; i++)
				byteCode[i] = RefCountedToken(new PlusOperatorToken()) ;
			
			break ;
		}
	}
}

Function::Function(const Function &f): derivative(f.derivative), iPoint(f.iPoint), ptID(f.ptID), dofID(f.dofID), byteCode(f.byteCode), e_diff(f.e_diff)
{
	for(std::map<int, Vector *>::const_iterator i = f.precalc.begin() ; i != f.precalc.end() ; ++i)
		precalc[i->first] = new Vector(*i->second) ;

	for(std::map<int, std::map<Variable, Vector *> >::const_iterator i = dprecalc.begin() ; i != dprecalc.end() ; ++i)
	{
		std::map<Variable, Vector *> v ;
		for(std::map<Variable, Vector *>::const_iterator j = i->second.begin() ; j != i->second.end() ; ++j)
		{
			v[j->first] = new Vector(*j->second) ;
		}
		dprecalc[i->first] = v ;
	}
}


Function::~Function()
{
	for(std::map<int, Vector *>::iterator i = precalc.begin() ; i != precalc.end() ; ++i)
		delete i->second ;

	for(std::map<int, std::map<Variable, Vector *> >::iterator i = dprecalc.begin() ; i != dprecalc.end() ; ++i)
	{
		for(std::map<Variable, Vector *>::iterator j = i->second.begin() ; j != i->second.end() ; ++j)
		{
			delete j->second ;
		}
	}
}


bool Function::isDifferentiable() const 
{
	return e_diff ;
}

const ByteCode & Function::getByteCode() const
{
	return byteCode ;
}

size_t Function::size() const
{
	return byteCode.size() ;
}
	
ByteCode & Function::getByteCode()
{
	return byteCode ;
}

Function Function::operator*(const Function &f) const
{
// 	if(f.isNull())
// 		return Function() ;
		
	return Function(byteCode, f.getByteCode(),  RefCountedToken(new TimesOperatorToken())) ;
}
	
Function Function::operator/(const Function &f) const 
{
	return Function(byteCode, f.getByteCode(), RefCountedToken(new DivideOperatorToken())) ;
}
	
Function Function::operator+(const Function &f) const
{
// 	if(f.isNull())
// 		return Function(*this) ;
		
	return Function(byteCode, f.getByteCode(), RefCountedToken(new PlusOperatorToken())) ;
}
	
Function operator-(const double & a, const Function &f)
{
	return Function(a, f.getByteCode(),  RefCountedToken(new MinusOperatorToken())) ;
}

Function operator*(const double & a, const Function &f)
{
	return Function(a, f.getByteCode(),  RefCountedToken(new TimesOperatorToken())) ;
}

Function operator+(const double & a, const Function &f)
{
	return Function(a, f.getByteCode(),  RefCountedToken(new PlusOperatorToken())) ;
}

Function operator/(const double & a, const Function &f)
{
	return Function(a, f.getByteCode(),  RefCountedToken(new DivideOperatorToken())) ;
}

Function Function::operator-(const Function &f) const 
{
// 	if(f.isNull())
// 		return Function(*this) ;
	return Function(byteCode, f.getByteCode(), RefCountedToken(new MinusOperatorToken())) ;
}

Function Function::operator*(const double a) const
{
// 	if(a == 0)
// 		return Function() ;
	return Function(byteCode, a, RefCountedToken(new TimesOperatorToken())) ;
}

Function Function::operator/(const double a) const 
{
	return Function(byteCode, a, RefCountedToken(new DivideOperatorToken())) ;
}

Function Function::operator+(const double a) const
{
// 	if(a == 0)
// 		return Function(*this)  ;
	return Function(byteCode, a, RefCountedToken(new PlusOperatorToken())) ;
}

Function Function::operator-(const double a) const 
{
// 	if(a == 0)
// 		return Function(*this) ;
	return Function(byteCode, a, RefCountedToken(new MinusOperatorToken())) ;
}

Function  Function::operator^(const int a) const
{
	return Function(byteCode, a, RefCountedToken(new PowerOperatorToken())) ;
}

void Function::operator*=(const Function &f) 
{
// 	if(f.isNull())
// 	{
// 		*this =  Function() ;
// 		return ;
// 	}
	*this =  Function(byteCode, f.getByteCode(), RefCountedToken(new TimesOperatorToken())) ;
	
}

void Function::operator/=(const Function &f)  
{
	*this =  Function(byteCode, f.getByteCode(), RefCountedToken(new DivideOperatorToken())) ;
}

void Function::operator+=(const Function &f) 
{
// 	if(f.isNull())
// 		return ;
	*this =  Function(byteCode, f.getByteCode(), RefCountedToken(new PlusOperatorToken())) ;
}

void Function::operator-=(const Function &f)  
{
// 	if(f.isNull())
// 		return ;
	*this =  Function(byteCode, f.getByteCode(), RefCountedToken(new MinusOperatorToken())) ;
}

void Function::operator*=(const double a) 
{
// 	if(a == 0)
// 	{
// 		*this = Function() ;
// 		return ;
// 	}
	*this =  Function(byteCode, a,RefCountedToken( new TimesOperatorToken())) ;
}

void Function::operator/=(const double a)  
{
	*this =  Function(byteCode, a, RefCountedToken(new DivideOperatorToken())) ;
}

void Function::operator+=(const double a) 
{
	if(a == 0)
	{
		return ;
	}
	*this =  Function(byteCode, a, RefCountedToken(new PlusOperatorToken())) ;
}

void Function::operator-=(const double a)  
{
	if(a == 0)
	{
		return ;
	}
	*this =  Function(byteCode, a, RefCountedToken(new MinusOperatorToken())) ;
}

const Function & Function::d(const Variable v) const
{
	return derivative[v] ;
}

Function & Function::d(const Variable v) 
{
	if(derivative.size() == 0)
	{
		derivative.resize(3) ;
		derivative = Function() ;
		e_diff = true ;
	}
	
	return derivative[v] ;
}

std::valarray<Function> Function::getDerivatives() const
{
	return derivative ;
}

std::valarray<Function> & Function::getDerivatives()
{
	return derivative ;
}


GtM Gradient::operator*( const Matrix & f) const
{
	return GtM(*this, f) ;
}

GtV Gradient::operator*( const Vector & f) const
{
	return GtV(*this, f) ;
}

GtMtG GtM::operator*(const Gradient & f) const
{
	return GtMtG(this->first, this->second, f) ;
}

Function f_exp(const Function &f)
{
	ByteCode b(f.getByteCode().size() + 1) ;
	for(size_t i = 0 ;  i < f.getByteCode().size() ; i++)
	{
		b[i] = f.getToken(i) ;
	}	
	b[f.getByteCode().size()] = new ExpToken() ;
	Function ret ;
	ret.getByteCode().resize(b.size()) ;
	ret.getByteCode() = b ;
	return ret ;
}

Function f_abs(const Function &f)
{
	ByteCode b(f.getByteCode().size() + 1) ;
	for(size_t i = 0 ;  i < f.getByteCode().size() ; i++)
	{
		b[i] = f.getToken(i) ;
	}	
	b[f.getByteCode().size()] = RefCountedToken(new AbsToken()) ;
	Function ret ;
	ret.getByteCode().resize(b.size()) ;
	ret.getByteCode() = b ;
	return ret ;
}

Function f_log(const Function &f)
{
	ByteCode b(f.getByteCode().size() + 1) ;
	for(size_t i = 0 ;  i < f.getByteCode().size() ; i++)
	{
		b[i] = f.getToken(i) ;
	}	
	b[f.getByteCode().size()] = RefCountedToken(new LogToken()) ;
	Function ret ;
	ret.getByteCode().resize(b.size()) ;
	ret.getByteCode() = b ;
	return ret ;
}

Function f_sqrt(const Function &f)
{
	ByteCode b(f.getByteCode().size() + 1) ;
	for(size_t i = 0 ;  i < f.getByteCode().size() ; i++)
	{
		b[i] = f.getToken(i) ;
	}	
	b[f.getByteCode().size()] = RefCountedToken(new SqrtToken()) ;
	Function ret ;
	ret.getByteCode().resize(b.size()) ;
	ret.getByteCode() = b ;
	return ret ;
}

Function f_atan2(const Function &f0, const Function &f1)
{
	ByteCode b(f0.getByteCode().size() + f0.getByteCode().size() + 1) ;
	for(size_t i = 0 ;  i < f0.getByteCode().size() ; i++)
	{
		b[i] = f0.getToken(i) ;
	}
	for(size_t i = 0 ;  i < f1.getByteCode().size() ; i++)
	{
		b[i+f0.getByteCode().size()] = f1.getToken(i) ;
	}
	b[f0.getByteCode().size()+ f1.getByteCode().size()] = RefCountedToken(new Atan2Token()) ;
	Function ret ;
	ret.getByteCode().resize(b.size()) ;
	ret.getByteCode() = b ;
	return ret ;
}

Function f_sin(const Function &f)
{
	ByteCode b(f.getByteCode().size() + 1) ;
	for(size_t i = 0 ;  i < f.getByteCode().size() ; i++)
	{
		b[i] = f.getToken(i) ;
	}	
	b[f.getByteCode().size()] = RefCountedToken(new SinToken()) ;
	Function ret ;
	ret.getByteCode().resize(b.size()) ;
	ret.getByteCode() = b ;
	return ret ;
}

Mu::Function f_sign(const Mu::Function &f)
{
	ByteCode b(f.getByteCode().size() + 1) ;
	for(size_t i = 0 ;  i < f.getByteCode().size() ; i++)
	{
		b[i] = f.getToken(i) ;
	}	
	b[f.getByteCode().size()] = RefCountedToken(new SignFunctionToken()) ;
	Function ret ;
	ret.getByteCode().resize(b.size()) ;
	ret.getByteCode() = b ;
	return ret ;
}

Mu::Function f_positivity(const Mu::Function &f)
{
	ByteCode b(f.getByteCode().size() + 1) ;
	for(size_t i = 0 ;  i < f.getByteCode().size() ; i++)
	{
		b[i] = f.getToken(i) ;
	}	
	b[f.getByteCode().size()] = RefCountedToken(new PositivityFunctionToken()) ;
	Function ret ;
	ret.getByteCode().resize(b.size()) ;
	ret.getByteCode() = b ;
	return ret ;
}

Mu::Function f_negativity(const Mu::Function &f)
{
	ByteCode b(f.getByteCode().size() + 1) ;
	for(size_t i = 0 ;  i < f.getByteCode().size() ; i++)
	{
		b[i] = f.getToken(i) ;
	}
	
	b[f.getByteCode().size()] = RefCountedToken(new NegativityFunctionToken()) ;
	Function ret ;
	ret.getByteCode().resize(b.size()) ;
	ret.getByteCode() = b ;
	return ret ;
}

Function f_cos(const Function &f)
{
	ByteCode b(f.getByteCode().size() + 1) ;
	for(size_t i = 0 ;  i < f.getByteCode().size() ; i++)
	{
		b[i] = f.getToken(i) ;
	}	
	b[f.getByteCode().size()] = new CosToken() ;
	Function ret ;
	ret.getByteCode().resize(b.size()) ;
	ret.getByteCode() = b ;
	return ret ;
}

const std::vector< Point > & Function::getIntegrationHint() const
{
	return iPoint ;
}

Point Function::getIntegrationHint(size_t i) const
{
	return iPoint[i] ;
}

void Function::setIntegrationHint(const std::vector< Point > v)
{
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





bool Function::isBinaryOperator(const Token * t) const
{
	return ( t->type.first.first == TOKEN_PLUS || 
	         t->type.first.first == TOKEN_MINUS || 
	         t->type.first.first == TOKEN_TIMES || 
	         t->type.first.first == TOKEN_DIVIDES || 
	         t->type.first.first == TOKEN_POWER || 
	         t->type.first.first == TOKEN_ATAN2) ;
}

bool Function::isUnaryOperator(const Token * t) const
{
	return ( t->type.first.first == TOKEN_SIN || 
	         t->type.first.first == TOKEN_COS || 
	         t->type.first.first == TOKEN_SQRT || 
	         t->type.first.first == TOKEN_COSH || 
	         t->type.first.first == TOKEN_SINH ||
	         t->type.first.first == TOKEN_TANH || 
	         t->type.first.first == TOKEN_EXP || 
	         t->type.first.first == TOKEN_LOG || 
	         t->type.first.first == TOKEN_TAN  
	       ) ;
}
void Function::compile()
{
// 	if(this->e_diff)
// 	{
		for(size_t i = 0 ; i < this->derivative.size() ; i++)
			if(this->derivative[i].getByteCode().size())
				derivative[i].compile() ;
// 	}
	
	size_t lastAddress = 0 ;
	std::vector<RefCountedToken> newByteCode ;
	std::vector<RefCountedToken> bytecode(byteCode.size()) ; 
	std::copy(&byteCode[0], &byteCode[byteCode.size()], bytecode.begin()) ;
	int precalculatedEnd = 0 ;
	
	std::vector<std::vector<RefCountedToken> > subexpressions ;
	std::vector<size_t> adresses ;
	
	//first, we vectorize what we can
	for(size_t i = 0 ; i < bytecode.size() ; i++)
	{
		bool foundBinOp = false ;
		int subexpressionEnd = i;
		int subexpressionStart = i ;
		int byteCodeSize = bytecode.size() ;
		while(!foundBinOp && subexpressionStart < byteCodeSize)
		{
			subexpressionEnd++ ;
			
			if(isBinaryOperator(bytecode[subexpressionStart]) && 
			   bytecode[subexpressionStart]->type.first.first != TOKEN_ATAN2)
				foundBinOp = true ;
			else
				subexpressionStart++ ;
			
		}
		
		if(foundBinOp && subexpressionEnd < byteCodeSize+1)
		{
			subexpressionStart = subexpressionEnd-3 ;
		}

		
		if(foundBinOp)
		{
			for(int l = 0 ; l < subexpressionStart ; l++)
			{
				newByteCode.push_back(bytecode[l]) ;
			}
			
			if(bytecode[subexpressionEnd-1]->type.first.first == TOKEN_PLUS)
			{
				if(bytecode[subexpressionEnd-3]->type.first.first == TOKEN_CONSTANT)
				{
					if(bytecode[subexpressionEnd-2]->type.first.first == TOKEN_X)
					{
						if(std::abs(bytecode[subexpressionEnd-3]->type.second) > 1e-18)
							newByteCode.push_back(new AddXAndConstToken(bytecode[subexpressionEnd-3]->type.second)) ;
						else
							newByteCode.push_back(new XToken()) ;
					}
					else if (bytecode[subexpressionEnd-2]->type.first.first == TOKEN_Y)
					{
						if(std::abs(bytecode[subexpressionEnd-3]->type.second) > 1e-18)
							newByteCode.push_back(new AddYAndConstToken(bytecode[subexpressionEnd-3]->type.second)) ;
						else
							newByteCode.push_back(new YToken()) ;
					}
					else if (bytecode[subexpressionEnd-2]->type.first.first == TOKEN_Z)
					{
						if(std::abs(bytecode[subexpressionEnd-3]->type.second) > 1e-18)
							newByteCode.push_back(new AddZAndConstToken(bytecode[subexpressionEnd-3]->type.second)) ;
						else
							newByteCode.push_back(new ZToken()) ;
					}
					else if (bytecode[subexpressionEnd-2]->type.first.first == TOKEN_T)
					{
						if(std::abs(bytecode[subexpressionEnd-3]->type.second) > 1e-18)
							newByteCode.push_back(new AddTAndConstToken(bytecode[subexpressionEnd-3]->type.second)) ;
						else
							newByteCode.push_back(new TToken()) ;
					}
					else
					{
						subexpressionEnd = subexpressionStart ;
					}
				}
				else if(bytecode[subexpressionEnd-2]->type.first.first == TOKEN_CONSTANT)
				{
					if(bytecode[subexpressionEnd-3]->type.first.first == TOKEN_X)
					{
						if(std::abs(bytecode[subexpressionEnd-2]->type.second) > 1e-18)
							newByteCode.push_back(new AddXAndConstToken(bytecode[subexpressionEnd-2]->type.second)) ;
						else
							newByteCode.push_back(new XToken()) ;
					}
					else if (bytecode[subexpressionEnd-3]->type.first.first == TOKEN_Y)
					{
						if(std::abs(bytecode[subexpressionEnd-2]->type.second) > 1e-18)
							newByteCode.push_back(new AddYAndConstToken(bytecode[subexpressionEnd-2]->type.second)) ;
						else
							newByteCode.push_back(new YToken()) ;
					}
					else if (bytecode[subexpressionEnd-3]->type.first.first == TOKEN_Z)
					{
						if(std::abs(bytecode[subexpressionEnd-2]->type.second) > 1e-18)
							newByteCode.push_back(new AddZAndConstToken(bytecode[subexpressionEnd-2]->type.second)) ;
						else
							newByteCode.push_back(new ZToken()) ;
					}
					else if (bytecode[subexpressionEnd-3]->type.first.first == TOKEN_T)
					{
						if(std::abs(bytecode[subexpressionEnd-2]->type.second) > 1e-18)
							newByteCode.push_back(new AddTAndConstToken(bytecode[subexpressionEnd-2]->type.second)) ;
						else
							newByteCode.push_back(new TToken()) ;
					}
					else
					{
						subexpressionEnd = subexpressionStart ;
					}
				}
				else
				{
					subexpressionEnd = subexpressionStart ;
				}
			}
			else if(bytecode[subexpressionEnd-1]->type.first.first == TOKEN_MINUS)
			{
				if(bytecode[subexpressionEnd-3]->type.first.first == TOKEN_CONSTANT)
				{
					if(bytecode[subexpressionEnd-2]->type.first.first == TOKEN_X)
					{
						newByteCode.push_back(new SubstractConstWithXToken(bytecode[subexpressionEnd-3]->type.second)) ;
					}
					else if (bytecode[subexpressionEnd-2]->type.first.first == TOKEN_Y)
					{
						newByteCode.push_back(new SubstractConstWithYToken(bytecode[subexpressionEnd-3]->type.second)) ;
					}
					else if (bytecode[subexpressionEnd-2]->type.first.first == TOKEN_Z)
					{
						newByteCode.push_back(new SubstractConstWithZToken(bytecode[subexpressionEnd-3]->type.second)) ;
					}
					else if (bytecode[subexpressionEnd-2]->type.first.first == TOKEN_T)
					{
						newByteCode.push_back(new SubstractConstWithTToken(bytecode[subexpressionEnd-3]->type.second)) ;
					}
					else
					{
						subexpressionEnd = subexpressionStart ;
					}
				}
				else if(bytecode[subexpressionEnd-2]->type.first.first == TOKEN_CONSTANT)
				{
					if(bytecode[subexpressionEnd-3]->type.first.first == TOKEN_X)
					{
						if(std::abs(bytecode[subexpressionEnd-2]->type.second) >1e-18)
							newByteCode.push_back(new SubstractXWithConstToken(bytecode[subexpressionEnd-2]->type.second)) ;
						else
							newByteCode.push_back(new XToken()) ;
					}
					else if (bytecode[subexpressionEnd-3]->type.first.first == TOKEN_Y)
					{
						if(std::abs(bytecode[subexpressionEnd-2]->type.second) >1e-18)
							newByteCode.push_back(new SubstractYWithConstToken(bytecode[subexpressionEnd-2]->type.second)) ;
						else
							newByteCode.push_back(new YToken()) ;
					}
					else if (bytecode[subexpressionEnd-3]->type.first.first == TOKEN_Z)
					{
						if(std::abs(bytecode[subexpressionEnd-2]->type.second) >1e-18)
							newByteCode.push_back(new SubstractZWithConstToken(bytecode[subexpressionEnd-2]->type.second)) ;
						else
							newByteCode.push_back(new ZToken()) ;
					}
					else if (bytecode[subexpressionEnd-3]->type.first.first == TOKEN_T)
					{
						if(std::abs(bytecode[subexpressionEnd-2]->type.second) >1e-18)
							newByteCode.push_back(new SubstractTWithConstToken(bytecode[subexpressionEnd-2]->type.second)) ;
						else
							newByteCode.push_back(new TToken()) ;
					}
					else
					{
						subexpressionEnd = subexpressionStart ;
					}
				}
				else
				{
					subexpressionEnd = subexpressionStart ;
				}
			}
			else if(bytecode[subexpressionEnd-1]->type.first.first == TOKEN_TIMES)
			{
				if(bytecode[subexpressionEnd-3]->type.first.first == TOKEN_CONSTANT)
				{
					if(bytecode[subexpressionEnd-2]->type.first.first == TOKEN_X)
					{
						if(std::abs(bytecode[subexpressionEnd-3]->type.second) >1e-18)
							newByteCode.push_back(new MultiplyXAndConstToken(bytecode[subexpressionEnd-3]->type.second)) ;
						else
							newByteCode.push_back(new ConstantToken(0)) ;
					}
					else if (bytecode[subexpressionEnd-2]->type.first.first == TOKEN_Y)
					{
						if(std::abs(bytecode[subexpressionEnd-3]->type.second) >1e-18)
							newByteCode.push_back(new MultiplyYAndConstToken(bytecode[subexpressionEnd-3]->type.second)) ;
						else
							newByteCode.push_back(new ConstantToken(0)) ;
					}
					else if (bytecode[subexpressionEnd-2]->type.first.first == TOKEN_Z)
					{
						if(std::abs(bytecode[subexpressionEnd-3]->type.second) >1e-18)
							newByteCode.push_back(new MultiplyZAndConstToken(bytecode[subexpressionEnd-3]->type.second)) ;
						else
							newByteCode.push_back(new ConstantToken(0)) ;
					}
					else if (bytecode[subexpressionEnd-2]->type.first.first == TOKEN_T)
					{
						if(std::abs(bytecode[subexpressionEnd-3]->type.second) >1e-18)
							newByteCode.push_back(new MultiplyTAndConstToken(bytecode[subexpressionEnd-3]->type.second)) ;
						else
							newByteCode.push_back(new ConstantToken(0)) ;
					}
					else
					{
						subexpressionEnd = subexpressionStart ;
					}
				}
				else if(bytecode[subexpressionEnd-2]->type.first.first == TOKEN_CONSTANT)
				{
					if(bytecode[subexpressionEnd-3]->type.first.first == TOKEN_X)
					{
						if(std::abs(bytecode[subexpressionEnd-2]->type.second) >1e-18 && std::abs(bytecode[subexpressionEnd-2]->type.second-1) >1e-18)
							newByteCode.push_back(new MultiplyXAndConstToken(bytecode[subexpressionEnd-2]->type.second)) ;
						else if(std::abs(bytecode[subexpressionEnd-2]->type.second-1) <= 1e-18)
							newByteCode.push_back(new XToken()) ;
						else
							newByteCode.push_back(new ConstantToken(0)) ;
					}
					else if (bytecode[subexpressionEnd-3]->type.first.first == TOKEN_Y)
					{
						if(std::abs(bytecode[subexpressionEnd-2]->type.second) >1e-18 && std::abs(bytecode[subexpressionEnd-2]->type.second-1) >1e-18)
							newByteCode.push_back(new MultiplyYAndConstToken(bytecode[subexpressionEnd-2]->type.second)) ;
						else if(std::abs(bytecode[subexpressionEnd-2]->type.second-1) <= 1e-18)
							newByteCode.push_back(new YToken()) ;
						else
							newByteCode.push_back(new ConstantToken(0)) ;
					}
					else if (bytecode[subexpressionEnd-3]->type.first.first == TOKEN_Z)
					{
						if(std::abs(bytecode[subexpressionEnd-2]->type.second) >1e-18 && std::abs(bytecode[subexpressionEnd-2]->type.second-1) >1e-18)
							newByteCode.push_back(new MultiplyZAndConstToken(bytecode[subexpressionEnd-2]->type.second)) ;
						else if(std::abs(bytecode[subexpressionEnd-2]->type.second-1) <= 1e-18)
							newByteCode.push_back(new ZToken()) ;
						else
							newByteCode.push_back(new ConstantToken(0)) ;
					}
					else if (bytecode[subexpressionEnd-3]->type.first.first == TOKEN_T)
					{
						if(std::abs(bytecode[subexpressionEnd-2]->type.second) >1e-18 && std::abs(bytecode[subexpressionEnd-2]->type.second-1) >1e-18)
							newByteCode.push_back(new MultiplyTAndConstToken(bytecode[subexpressionEnd-2]->type.second)) ;
						else if(std::abs(bytecode[subexpressionEnd-2]->type.second-1) <= 1e-18)
							newByteCode.push_back(new TToken()) ;
						else
							newByteCode.push_back(new ConstantToken(0)) ;
					}
					else
					{
						subexpressionEnd = subexpressionStart ;
					}
				}
				else
				{
					subexpressionEnd = subexpressionStart ;
				}
			}
			else if(bytecode[subexpressionEnd-1]->type.first.first == TOKEN_DIVIDES)
			{
				if(bytecode[subexpressionEnd-3]->type.first.first == TOKEN_CONSTANT)
				{
					if(bytecode[subexpressionEnd-2]->type.first.first == TOKEN_X)
					{
						if(std::abs(bytecode[subexpressionEnd-3]->type.second) >1e-18)
							newByteCode.push_back(new DivideConstWithXToken(bytecode[subexpressionEnd-3]->type.second)) ;
						else
							newByteCode.push_back(new ConstantToken(0)) ;
					}
					else if (bytecode[subexpressionEnd-2]->type.first.first == TOKEN_Y)
					{
						if(std::abs(bytecode[subexpressionEnd-3]->type.second) >1e-18)
							newByteCode.push_back(new DivideConstWithYToken(bytecode[subexpressionEnd-3]->type.second)) ;
						else
							newByteCode.push_back(new ConstantToken(0)) ;
					}
					else if (bytecode[subexpressionEnd-2]->type.first.first == TOKEN_Z)
					{
						if(std::abs(bytecode[subexpressionEnd-3]->type.second) >1e-18)
							newByteCode.push_back(new DivideConstWithZToken(bytecode[subexpressionEnd-3]->type.second)) ;
						else
							newByteCode.push_back(new ConstantToken(0)) ;
					}
					else if (bytecode[subexpressionEnd-2]->type.first.first == TOKEN_T)
					{
						if(std::abs(bytecode[subexpressionEnd-3]->type.second) >1e-18)
							newByteCode.push_back(new DivideConstWithTToken(bytecode[subexpressionEnd-3]->type.second)) ;
						else
							newByteCode.push_back(new ConstantToken(0)) ;
					}
					else
					{
						subexpressionEnd = subexpressionStart ;
					}
				}
				else if(bytecode[subexpressionEnd-2]->type.first.first == TOKEN_CONSTANT)
				{
					if(bytecode[subexpressionEnd-3]->type.first.first == TOKEN_X)
					{
						newByteCode.push_back(new DivideXWithConstToken(bytecode[subexpressionEnd-2]->type.second)) ;
					}
					else if (bytecode[subexpressionEnd-3]->type.first.first == TOKEN_Y)
					{
						newByteCode.push_back(new DivideYWithConstToken(bytecode[subexpressionEnd-2]->type.second)) ;
					}
					else if (bytecode[subexpressionEnd-3]->type.first.first == TOKEN_Z)
					{
						newByteCode.push_back(new DivideZWithConstToken(bytecode[subexpressionEnd-2]->type.second)) ;
					}
					else if (bytecode[subexpressionEnd-3]->type.first.first == TOKEN_T)
					{
						newByteCode.push_back(new DivideTWithConstToken(bytecode[subexpressionEnd-2]->type.second)) ;
					}
					else
					{
						subexpressionEnd = subexpressionStart ;
					}
				}
				else
				{
					subexpressionEnd = subexpressionStart ;
				}
			}
			else if(bytecode[subexpressionEnd-1]->type.first.first == TOKEN_POWER)
			{
				if(bytecode[subexpressionEnd-2]->type.first.first == TOKEN_CONSTANT)
				{
					if(bytecode[subexpressionEnd-3]->type.first.first == TOKEN_X)
					{
						newByteCode.push_back(new XPowerConstToken(bytecode[subexpressionEnd-2]->type.second)) ;
					}
					else if (bytecode[subexpressionEnd-3]->type.first.first == TOKEN_Y)
					{
						newByteCode.push_back(new YPowerConstToken(bytecode[subexpressionEnd-2]->type.second)) ;
					}
					else if (bytecode[subexpressionEnd-3]->type.first.first == TOKEN_Z)
					{
						newByteCode.push_back(new ZPowerConstToken(bytecode[subexpressionEnd-2]->type.second)) ;
					}
					else if (bytecode[subexpressionEnd-3]->type.first.first == TOKEN_T)
					{
						newByteCode.push_back(new TPowerConstToken(bytecode[subexpressionEnd-2]->type.second)) ;
					}
					else
					{
						subexpressionEnd = subexpressionStart ;
					}
				}
				else
				{
					subexpressionEnd = subexpressionStart ;
				}
			}
			
			for(int j = subexpressionEnd ; j < byteCodeSize ; j++)
			{
				newByteCode.push_back(bytecode[j]) ;
			}

			bytecode=newByteCode ;
			newByteCode.clear() ;
		}
		
	}
	
	//Then, we factorize
	while(true)
	{
		bool foundBinOp = false ;
		bool foundUnOp = false ;
		int subexpressionEnd = precalculatedEnd+1;
		int subexpressionStart = precalculatedEnd+1 ;
		int byteCodeSize = bytecode.size() ;
		while((!foundBinOp && !foundUnOp) && subexpressionStart < byteCodeSize-2)
		{
			subexpressionEnd++ ;
			
			if(isBinaryOperator(bytecode[subexpressionStart]))
				foundBinOp = true ;
			else if(isUnaryOperator(bytecode[subexpressionStart]))
				foundUnOp = true ;
			else
				subexpressionStart++ ;
		}
		
		if(foundBinOp && subexpressionEnd < byteCodeSize-1)
		{
			subexpressionStart = subexpressionEnd-3 ;
		}
		else if(foundUnOp&& subexpressionEnd < byteCodeSize-1)
		{
			subexpressionStart = subexpressionEnd-2 ;
		}
		else
			break ;
		
		std::vector<RefCountedToken> subexpression ; 
		
			//we keep the precalculated stuff ;
		for(int l = 0 ; l < precalculatedEnd ; l++)
		{
			newByteCode.push_back(bytecode[l]) ;
		}
		for(int l = subexpressionStart ; l < subexpressionEnd ; l++)
		{
			subexpression.push_back(bytecode[l]) ;
			newByteCode.push_back(bytecode[l]) ;
		}
		newByteCode.push_back(new SetHeapVariableToken(lastAddress)) ;
		size_t replacementStart = precalculatedEnd ;
		precalculatedEnd = newByteCode.size() ;
		
		size_t occurences = 0 ;
		for(int l = replacementStart ; l < byteCodeSize ; l++)
		{
			size_t idx = 0 ;
			
			if(subexpression[0]->type == bytecode[l]->type)
			{
				idx++ ;
				while(idx < subexpression.size() && subexpression[idx]->type == bytecode[idx+l]->type)
				{
					idx++ ;
				}
				
				if(subexpression.size()==idx)
				{
					l+=subexpression.size()-1 ;
					occurences++ ;
					newByteCode.push_back(new ReadHeapVariableToken(lastAddress)) ;
				}
				else 
					newByteCode.push_back(bytecode[l]) ;
			}
			else
				newByteCode.push_back(bytecode[l]) ;
		}
		
		if(occurences < 3) // no point reading and writing for single occurences
		{
			adresses.push_back(lastAddress) ;
			subexpressions.push_back(subexpression) ;
		}
		
		lastAddress++ ;
		bytecode=newByteCode ;
		newByteCode.clear() ;
	}
	
	//remove duplicates
	if(!adresses.empty())
	{
		for(int i = adresses.size()-1 ; i >=0 ; i--)
		{
			TokenType sin(std::make_pair(TOKEN_READ_VARIABLE,adresses[i]), (double)(0)) ;
			TokenType sout(std::make_pair(TOKEN_WRITE_VARIABLE,adresses[i]), (double)(0)) ;	
			size_t subsize = subexpressions[i].size() ;

			for(std::vector<RefCountedToken>::iterator j = bytecode.begin() ; j != bytecode.end() ; ++j)
			{
				if((*j)->type == sout)
				{
					bytecode.erase(j-subsize, j+1) ;
					j = j - subsize;
				}
				if((*j)->type == sin)
				{
					bytecode.erase(j) ;
					bytecode.insert(j,subexpressions[i].begin(), subexpressions[i].end() ) ;
					j = j + subsize ;
				}
			}
		}
	}

	//last, we vectorize again
	for(size_t i = 0 ; i < bytecode.size() ; i++)
	{
		bool foundBinOp = false ;
		int subexpressionEnd = i;
		int subexpressionStart = i ;
		int byteCodeSize = bytecode.size() ;
		while(!foundBinOp && subexpressionStart < byteCodeSize)
		{
			subexpressionEnd++ ;
			
			if(isBinaryOperator(bytecode[subexpressionStart]) && 
			   bytecode[subexpressionStart]->type.first.first != TOKEN_ATAN2)
				foundBinOp = true ;
			else
				subexpressionStart++ ;
			
		}
		
		if(foundBinOp && subexpressionEnd < byteCodeSize+1)
		{
			subexpressionStart = subexpressionEnd-3 ;
		}
		
		
		if(foundBinOp)
		{
			for(int l = 0 ; l < subexpressionStart ; l++)
			{
				newByteCode.push_back(bytecode[l]) ;
			}
			
			switch(bytecode[subexpressionEnd-1]->type.first.first)
			{
			case TOKEN_PLUS :
			{
				switch(bytecode[subexpressionEnd-3]->type.first.first)
				{
				case TOKEN_CONSTANT :
				{
					switch(bytecode[subexpressionEnd-2]->type.first.first)
					{
					case TOKEN_READ_VARIABLE :
						{
							if(std::abs(bytecode[subexpressionEnd-3]->type.second)> 1e-18)
								newByteCode.push_back(new AddReadAndConstToken(bytecode[subexpressionEnd-2]->type.first.second, bytecode[subexpressionEnd-3]->type.second)) ;
							else
								newByteCode.push_back(new ReadHeapVariableToken(bytecode[subexpressionEnd-2]->type.first.second)) ;
							break ;
						}
					case TOKEN_NULL :
						{
							newByteCode.push_back(bytecode[subexpressionEnd-3]) ;
							break ;
						}
					default :
						{
							subexpressionEnd = subexpressionStart ;
							break ;
						}
					}
					break ;
				}
				case TOKEN_READ_VARIABLE :
				{
					switch(bytecode[subexpressionEnd-2]->type.first.first)
					{
						case TOKEN_CONSTANT :
						{
							if(std::abs(bytecode[subexpressionEnd-2]->type.second) > 1e-18)
								newByteCode.push_back(new AddReadAndConstToken(bytecode[subexpressionEnd-3]->type.first.second, bytecode[subexpressionEnd-2]->type.second)) ;
							else
								newByteCode.push_back(new ReadHeapVariableToken(bytecode[subexpressionEnd-3]->type.first.second)) ;
							break ;
						}
						case TOKEN_READ_VARIABLE :
						{
							newByteCode.push_back(new AddReadAndReadToken(bytecode[subexpressionEnd-3]->type.first.second, bytecode[subexpressionEnd-2]->type.first.second)) ;
							break ;
						}
						case TOKEN_NULL :
						{
							newByteCode.push_back(bytecode[subexpressionEnd-3]) ;
							break ;
						}
						default :
						{
							subexpressionEnd = subexpressionStart ;
							break ;
						}
					}
					break ;
				}
				case TOKEN_NULL :
				{
					if(bytecode[subexpressionEnd-2]->type.first.first == TOKEN_CONSTANT ||
					   bytecode[subexpressionEnd-2]->type.first.first == TOKEN_READ_VARIABLE)
					{
						newByteCode.push_back(bytecode[subexpressionEnd-2]) ;
					}
					else
					{
						subexpressionEnd = subexpressionStart ;
					}
					break ;
				}
				default :
				{
					subexpressionEnd = subexpressionStart ;
					break ;
				}
				}
				break ;
			}
			case TOKEN_POWER :
			{
				if(bytecode[subexpressionEnd-2]->type.first.first == TOKEN_CONSTANT)
				{
					if(bytecode[subexpressionEnd-3]->type.first.first == TOKEN_READ_VARIABLE)
					{
						newByteCode.push_back(new ReadPowerConstToken(bytecode[subexpressionEnd-3]->type.first.second, bytecode[subexpressionEnd-2]->type.second)) ;
						
					}
					else
					{
						subexpressionEnd = subexpressionStart ;
					}
				}
				else
				{
					subexpressionEnd = subexpressionStart ;
				}
				break ;
			}
			case TOKEN_TIMES :
			{
				switch(bytecode[subexpressionEnd-3]->type.first.first)
				{
				case TOKEN_CONSTANT :
				{
					switch(bytecode[subexpressionEnd-2]->type.first.first)
					{
					case TOKEN_READ_VARIABLE :
						{
							if(std::abs(bytecode[subexpressionEnd-3]->type.second) > 1e-18 
							   && std::abs(bytecode[subexpressionEnd-3]->type.second - 1 ) > 1e-18)
								newByteCode.push_back(new MultiplyReadAndConstToken(bytecode[subexpressionEnd-2]->type.first.second, bytecode[subexpressionEnd-3]->type.second)) ;
							else if (std::abs(bytecode[subexpressionEnd-3]->type.second - 1 ) <= 1e-18)
								newByteCode.push_back(new ReadHeapVariableToken(bytecode[subexpressionEnd-2]->type.first.second)) ;
							else
								newByteCode.push_back(new ConstantToken(0)) ;
							break ;
						}
					case TOKEN_X_POWER_CONST :
						{
							newByteCode.push_back(new XPowerAndMultiplyToken(bytecode[subexpressionEnd-3]->type.second, (size_t)bytecode[subexpressionEnd-2]->type.second )) ;
							break ;
						}
					case TOKEN_Y_POWER_CONST :
						{
							newByteCode.push_back(new YPowerAndMultiplyToken(bytecode[subexpressionEnd-3]->type.second, (size_t)bytecode[subexpressionEnd-2]->type.second )) ;
							break ;
						}
					case TOKEN_Z_POWER_CONST :
						{
							newByteCode.push_back(new ZPowerAndMultiplyToken(bytecode[subexpressionEnd-3]->type.second, (size_t)bytecode[subexpressionEnd-2]->type.second )) ;
							break ;
						}
					case TOKEN_T_POWER_CONST :
						{
							newByteCode.push_back(new TPowerAndMultiplyToken(bytecode[subexpressionEnd-3]->type.second, (size_t)bytecode[subexpressionEnd-2]->type.second )) ;
							break ;
						}
					default :
						{
							subexpressionEnd = subexpressionStart ;
							break ;
						}
					}
					break ;
				}
				case TOKEN_MULTIPLY_X_AND_CONST :
				{
					switch(bytecode[subexpressionEnd-2]->type.first.first)
					{
						case TOKEN_Y :
						{
							newByteCode.push_back(new XYConstToken(bytecode[subexpressionEnd-3]->type.second)) ;
							break ;
						}
						case TOKEN_Z :
						{
							newByteCode.push_back(new XZConstToken(bytecode[subexpressionEnd-3]->type.second)) ;
							break ;
						}
						default :
						{
							subexpressionEnd = subexpressionStart ;
							break ;
						}
					}
					break ;
				}
				case TOKEN_MULTIPLY_Y_AND_CONST :
				{
					switch(bytecode[subexpressionEnd-2]->type.first.first)
					{
						case TOKEN_X :
						{
							newByteCode.push_back(new XYConstToken(bytecode[subexpressionEnd-3]->type.second)) ;
							break ;
						}
						case TOKEN_Z :
						{
							newByteCode.push_back(new YZConstToken(bytecode[subexpressionEnd-3]->type.second)) ;
							break ;
						}
						default :
						{
							subexpressionEnd = subexpressionStart ;
							break ;
						}
					}
					break ;
				}
				case TOKEN_MULTIPLY_Z_AND_CONST :
				{
					switch(bytecode[subexpressionEnd-2]->type.first.first)
					{
						case TOKEN_X :
						{
							newByteCode.push_back(new XYConstToken(bytecode[subexpressionEnd-3]->type.second)) ;
							break ;
						}
						case TOKEN_Y :
						{
							newByteCode.push_back(new YZConstToken(bytecode[subexpressionEnd-3]->type.second)) ;
							break ;
						}
						default :
						{
							subexpressionEnd = subexpressionStart ;
							break ;
						}
					}
					break ;
				}
				case TOKEN_READ_VARIABLE :
				{
					switch(bytecode[subexpressionEnd-2]->type.first.first)
					{
						case TOKEN_CONSTANT :
						{
							if(std::abs(bytecode[subexpressionEnd-2]->type.second) > 1e-18 && std::abs(bytecode[subexpressionEnd-2]->type.second -1 ) > 1e-18)
								newByteCode.push_back(new MultiplyReadAndConstToken(bytecode[subexpressionEnd-3]->type.first.second, bytecode[subexpressionEnd-2]->type.second)) ;
							else if (std::abs(bytecode[subexpressionEnd-2]->type.second -1 ) <= 1e-18 )
								newByteCode.push_back(new ReadHeapVariableToken(bytecode[subexpressionEnd-3]->type.first.second)) ;
							else
								newByteCode.push_back(new ConstantToken(0)) ;
							break ;
						}
						case TOKEN_READ_VARIABLE :
						{
							newByteCode.push_back(new MultiplyReadAndReadToken(bytecode[subexpressionEnd-3]->type.first.second, bytecode[subexpressionEnd-2]->type.first.second)) ;
							break ;
						}
						default :
						{
							subexpressionEnd = subexpressionStart ;
							break ;
						}
					}
					break ;
				}
				case TOKEN_X_POWER_CONST :
				{
					if(bytecode[subexpressionEnd-2]->type.first.first == TOKEN_CONSTANT)
					{
						newByteCode.push_back(new XPowerAndMultiplyToken(bytecode[subexpressionEnd-2]->type.second, (size_t)bytecode[subexpressionEnd-3]->type.second )) ;
					}
					else
					{
						subexpressionEnd = subexpressionStart ;
					}
					break ;
				}
				case TOKEN_Y_POWER_CONST :
				{
					if(bytecode[subexpressionEnd-2]->type.first.first == TOKEN_CONSTANT)
					{
						newByteCode.push_back(new YPowerAndMultiplyToken(bytecode[subexpressionEnd-2]->type.second, (size_t)bytecode[subexpressionEnd-3]->type.second )) ;
					}
					else
					{
						subexpressionEnd = subexpressionStart ;
					}
					break ;
				}
				case TOKEN_Z_POWER_CONST :
				{
					if(bytecode[subexpressionEnd-2]->type.first.first == TOKEN_CONSTANT)
					{
						newByteCode.push_back(new ZPowerAndMultiplyToken(bytecode[subexpressionEnd-2]->type.second, (size_t)bytecode[subexpressionEnd-3]->type.second )) ;
					}
					else
					{
						subexpressionEnd = subexpressionStart ;
					}
					break ;
				}
				case TOKEN_T_POWER_CONST :
					{
						if(bytecode[subexpressionEnd-2]->type.first.first == TOKEN_CONSTANT)
						{
							newByteCode.push_back(new TPowerAndMultiplyToken(bytecode[subexpressionEnd-2]->type.second, (size_t)bytecode[subexpressionEnd-3]->type.second )) ;
						}
						else
						{
							subexpressionEnd = subexpressionStart ;
						}
						break ;
					}
				default :
				{
					subexpressionEnd = subexpressionStart ;
					break ;
				}
				}
				break ;
			}
			case TOKEN_MINUS :
			{
				if(bytecode[subexpressionEnd-3]->type.first.first == TOKEN_CONSTANT)
				{
					if(bytecode[subexpressionEnd-2]->type.first.first == TOKEN_READ_VARIABLE)
					{
						newByteCode.push_back(new SubstractReadWithConstToken(bytecode[subexpressionEnd-2]->type.first.second, bytecode[subexpressionEnd-3]->type.second)) ;
					}
					else
					{
						subexpressionEnd = subexpressionStart ;
					}
				}
				else if(bytecode[subexpressionEnd-2]->type.first.first == TOKEN_CONSTANT)
				{
					if(bytecode[subexpressionEnd-3]->type.first.first == TOKEN_READ_VARIABLE)
					{
						newByteCode.push_back(new SubstractConstWithReadToken(bytecode[subexpressionEnd-3]->type.first.second, bytecode[subexpressionEnd-2]->type.second)) ;
					}
					else
					{
						subexpressionEnd = subexpressionStart ;
					}
				}
				else
				{
					subexpressionEnd = subexpressionStart ;
				}
				break ;
			}
			case TOKEN_DIVIDES :
			{
				if(bytecode[subexpressionEnd-3]->type.first.first == TOKEN_CONSTANT)
				{
					if(bytecode[subexpressionEnd-2]->type.first.first == TOKEN_READ_VARIABLE)
					{
						newByteCode.push_back(new DivideReadWithConstToken(bytecode[subexpressionEnd-2]->type.first.second, bytecode[subexpressionEnd-3]->type.second)) ;
					}
					else
					{
						subexpressionEnd = subexpressionStart ;
					}
				}
				else if(bytecode[subexpressionEnd-2]->type.first.first == TOKEN_CONSTANT)
				{
					if(bytecode[subexpressionEnd-3]->type.first.first == TOKEN_READ_VARIABLE)
					{
						newByteCode.push_back(new DivideConstWithReadToken(bytecode[subexpressionEnd-3]->type.first.second, bytecode[subexpressionEnd-2]->type.second)) ;
					}
					else
					{
						subexpressionEnd = subexpressionStart ;
					}
				}
				else
				{
					subexpressionEnd = subexpressionStart ;
				}
				break ;
			}
			default:
				break ;
			}

			for(int j = subexpressionEnd ; j < byteCodeSize ; j++)
			{
				newByteCode.push_back(bytecode[j]) ;
			}
			
			bytecode=newByteCode ;
			newByteCode.clear() ;

		}
	}
	
	byteCode.resize(bytecode.size()) ;
	std::copy(bytecode.begin(), bytecode.end(), &byteCode[0]) ;

}


DtF Differential::operator *(const Function & f) const
{
	return DtF(*this, f) ;
}

DtGtMtG Differential::operator *(const GtMtG & g) const
{
	return DtGtMtG(*this, g) ;
}

GDtMtGD GDtM::operator*(const Mu::GradientDot & f) const
{
	return GDtMtGD(*this, f) ;
}

GDtM GradientDot::operator *(const Matrix & f) const
{
	return GDtM(*this, f) ;
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



