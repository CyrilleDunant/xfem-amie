//
// C++ Implementation: vm_token
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "vm_token.h"

using namespace Mu ;



ConstantToken::ConstantToken(double v , bool nul) : Token(nul, std::make_pair(std::make_pair(TOKEN_CONSTANT, 0), (double)(v))){ };
	
ConstantToken::~ConstantToken() { } ;

NullToken::~NullToken() { }

VariableToken::VariableToken(const Variable v, bool nul) :Token(nul, std::make_pair(std::make_pair(TOKEN_VARIABLE, 0),(double)(0))), var(v)
{
};

VariableToken::~VariableToken() { } ;

PositionToken::~PositionToken() { } ;

DomainToken::~DomainToken() { } ;

ProjectionToken::~ProjectionToken() { } ;

UnaryFunctionToken::UnaryFunctionToken(unaryFunctionPointer f, bool nul) :Token(nul, std::make_pair(std::make_pair(TOKEN_UNARY_FUNCTION, 0),(double)(0))),fctPtr(f) { };

UnaryFunctionToken::~UnaryFunctionToken() { } ;

BinaryFunctionToken::BinaryFunctionToken(binaryFunctionPointer f, bool nul) :Token(nul, std::make_pair(std::make_pair(TOKEN_BINARY_FUNCTION, 0),(double)(0))), fctPtr(f) { };

BinaryFunctionToken::~BinaryFunctionToken() { } ;

PlusOperatorToken::PlusOperatorToken( bool nul) :Token(nul, std::make_pair(std::make_pair(TOKEN_PLUS, 0),(double)(0))) { };

PlusOperatorToken::~PlusOperatorToken() { } ;

MinusOperatorToken::MinusOperatorToken( bool nul) :Token(nul, std::make_pair(std::make_pair(TOKEN_MINUS, 0),(double)(0))) { };

MinusOperatorToken::~MinusOperatorToken() { } ;

TimesOperatorToken ::TimesOperatorToken( bool nul) :Token(nul, std::make_pair(std::make_pair(TOKEN_TIMES, 0),(double)(0))){ };

TimesOperatorToken::~TimesOperatorToken() { } ;

DivideOperatorToken::DivideOperatorToken( bool nul) :Token(nul, std::make_pair(std::make_pair(TOKEN_DIVIDES, 0),(double)(0))){ };

DivideOperatorToken::~DivideOperatorToken() { } ;

PowerOperatorToken::PowerOperatorToken( bool nul) :Token(nul, std::make_pair(std::make_pair(TOKEN_POWER, 0),(double)(0))) { };

PowerOperatorToken::~PowerOperatorToken() { } ;


double sign(const double t)
{
	if(t < 0)
		return -1 ;
	if(t > 0)
		return 1 ;
	return 0 ;
} 

double positivity(const double t)
{
	return t > 0 ;
} 

double negativity(const double t)
{
	return t < 0 ;
} 
