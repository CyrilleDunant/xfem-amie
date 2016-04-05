// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "function_parser.h"
#include "../font.h"
#include "../enumeration_translator.h"
#include "../../polynomial/vm_function_extra.h"
#include <sys/stat.h>

using namespace Amie ;

FunctionParser::FunctionParser(std::vector<std::string> f, std::map<std::string, char> coordinates) : isLinking(false)
{
	renewExpression(f, coordinates) ;
}

FunctionParser::FunctionParser(std::string f, std::map<std::string, char> coordinates) : isLinking(false)
{
	renewExpression( FunctionParser::breakString(f), coordinates) ;
}

std::vector<std::string> FunctionParser::breakString(std::string f)
{
	std::string current ;
	std::vector<std::string> str ;
	size_t i = 0 ;
	while(i < f.length())
	{
		bool cut = false ;
		char test = f[i] ;
		if(test == ' ' || test == '\t')
			cut = true ;
		else
			current += test ;
		if(cut)
		{
			if(current.length() > 0)
			{
				str.push_back( current ) ;
				current.clear() ;
			}
		}
		i++ ;
	}
	if(current.length() > 0)
		str.push_back( current ) ;
	return str ;
}

bool isOpenBracket( std::string s) 
{
	return s== "(" ;
}

bool isCloseBracket(std::string s) 
{
	return s==")" ;
}

std::string FunctionParser::translateCoordinate( std::string test, std::map<std::string, char> coordinates)
{
	if(coordinates.find(test) != coordinates.end())
		return std::string(1, coordinates[test]) ;
	return test ;
}

void FunctionParser::renewExpression(std::vector<std::string> str, std::map<std::string, char> coordinates) 
{
	size_t i = 0 ;
	std::vector<std::string> embeddedExpression ;
	size_t openCount = 0 ;
	
	while(i < str.size())
	{
		std::string translated = translateCoordinate( str[i], coordinates ) ;
		if(openCount == 0 && !isOpenBracket( translated )  )
		{
			if(translated.size() > 1 && translated[0] == '-' && !FunctionParserHelper::isNumeral(translated[1]))
			{
				std::vector<std::string> tmp ;
				tmp.push_back("-1") ;
				tmp.push_back("*") ;
				tmp.push_back(translateCoordinate( str[i].substr(1), coordinates )) ;
				tokens.push_back( new FunctionParser( tmp ) ) ;
			}
			else
				tokens.push_back( new FunctionParserToken( translated ) ) ;
		}
		else
		{
			if( isOpenBracket( translated ) )
			{
				if(openCount == 0)
					embeddedExpression.clear() ;
				else
					embeddedExpression.push_back( translated ) ;
				openCount++ ;
			}
			else if(isCloseBracket( str[i] ) )
			{
				openCount-- ;
				if(openCount == 0)
					tokens.push_back( new FunctionParser( embeddedExpression) ) ;
				else
					embeddedExpression.push_back( translated ) ;
			}
			else
				embeddedExpression.push_back( translated ) ;
		}
		i++ ;
	}

	this->link() ;
}

void FunctionParser::linkLeftAndRightToken( TokenOperationType op)
{
	for(size_t i = 0 ; i < tokens.size() ; i++)
	{
		if(tokens[i]->isFinalToken())
		{
			if( dynamic_cast<FunctionParserToken *>(tokens[i])->toFunctionToken() == op )
			{
				size_t left = i-1 ;
				size_t c = 0 ;
				while(roots[left] > -1 && c < tokens.size())
				{
					left = roots[left] ;
					c++;
				}
				size_t right = i+1 ;
				c = 0 ;
				while(roots[right] > -1 && c < tokens.size())
				{
					right = roots[right] ;
					c++;
				}
				if(roots[right] == -1 && right != i)
					roots[right] = i ;
				if(roots[left] == -1 && left != i)
					roots[left] = i ;
			}
		}
	}
}

void FunctionParser::link()
{
	if(isLinking || isFinalToken() )
		return ;

	isLinking=true ;

	roots.resize(tokens.size(), -1) ;

	for(size_t i = 0 ; i < tokens.size() ; i++)
	{
		if(tokens[i]->isFinalToken())
		{
			TokenOperationType elem = dynamic_cast<FunctionParserToken *>(tokens[i])->toFunctionToken() ;
			switch(elem)
			{
				case TOKEN_OPERATION_ABS:
				case TOKEN_OPERATION_COS:
				case TOKEN_OPERATION_SIN:
				case TOKEN_OPERATION_TAN:
				case TOKEN_OPERATION_COSH:
				case TOKEN_OPERATION_SINH:
				case TOKEN_OPERATION_TANH:
				case TOKEN_OPERATION_EXP:
				case TOKEN_OPERATION_SIGN:
				case TOKEN_OPERATION_POSITIVITY:
				case TOKEN_OPERATION_NEGATIVITY:
				case TOKEN_OPERATION_LOG:
				case TOKEN_OPERATION_SQRT:
					roots[i+1] = i ;
					break ;
				case TOKEN_OPERATION_MIN:
				case TOKEN_OPERATION_MAX:
					roots[i+1] = i ;
					roots[i+2] = i ;
					break ;
				default:
					break ;
			}
		}
//		else
//			tokens[i]->link() ;
	}

	linkLeftAndRightToken( TOKEN_OPERATION_POWER ) ;
	linkLeftAndRightToken( TOKEN_OPERATION_TIMES ) ;
	linkLeftAndRightToken( TOKEN_OPERATION_DIVIDES ) ;
	linkLeftAndRightToken( TOKEN_OPERATION_MINUS ) ;
	linkLeftAndRightToken( TOKEN_OPERATION_PLUS ) ;

//	printRoots() ;

//	for(size_t i = 0 ; i < tokens.size() ; i++)
//		std::cout << roots[i] << std::endl ;

	isLinking = false ;

}

int FunctionParser::getLeftTokenIndex( size_t i ) const 
{
	for(size_t j = 0 ; j < roots.size() ; j++)
	{
		if(roots[j] == (int) i)
			return j ;
	}
	return -1 ;
}

int FunctionParser::getRightTokenIndex( size_t i ) const 
{
	for(size_t j = 0 ; j < roots.size() ; j++)
	{
		if(roots[roots.size()-1-j] == (int) i)
			return roots.size()-1-j ;
	}
	return -1 ;
}

Function FunctionParser::getLeftFunction( size_t i ) const 
{
	int j = getLeftTokenIndex(i) ;
	if(j >= 0)
		return getFunction(j) ;
	return Function("0") ;
}

Function FunctionParser::getRightFunction( size_t i ) const 
{
	int j = getRightTokenIndex(i) ;
	if(j >= 0)
		return getFunction(j) ;
	return Function("0") ;
}

Function FunctionParser::getFunction() const
{
	for(size_t i = 0 ; i < roots.size() ; i++)
	{
		if(roots[i] == -1)
		{
			return getFunction(i) ;
		}
	}
	return Function("0") ;
}

Function FunctionParser::getFunction( size_t i ) const 
{
	if( tokens[i]->isFinalToken())
	{
		FunctionParserToken * local = dynamic_cast<FunctionParserToken *>( tokens[i] ) ; 
		TokenOperationType type = local->toFunctionToken() ;
		switch(type)
		{
			case TOKEN_OPERATION_CONSTANT:
				return Function( (local->data).c_str() ) ;
			case TOKEN_OPERATION_X:
				return Function("x") ;
			case TOKEN_OPERATION_Y:
				return Function("y") ;
			case TOKEN_OPERATION_Z:
				return Function("z") ;
			case TOKEN_OPERATION_T:
				return Function("t") ;
			case TOKEN_OPERATION_U:
				return Function("u") ;
			case TOKEN_OPERATION_V:
				return Function("v") ;
			case TOKEN_OPERATION_W:
				return Function("w") ;
			case TOKEN_OPERATION_PLUS:
			{
				if(getLeftTokenIndex(i) == getRightTokenIndex(i))
					return getRightFunction(i) ;
				return getLeftFunction(i)+getRightFunction(i) ;
			}
			case TOKEN_OPERATION_MINUS:
			{
				if(getLeftTokenIndex(i) == getRightTokenIndex(i))
					return -getRightFunction(i) ;
				return getLeftFunction(i)-getRightFunction(i) ;
			}
			case TOKEN_OPERATION_TIMES:
			{
				if(getLeftTokenIndex(i) == getRightTokenIndex(i))
					return getRightFunction(i) ;
				return getLeftFunction(i)*getRightFunction(i) ;
			}
			case TOKEN_OPERATION_DIVIDES:
			{
				if(getLeftTokenIndex(i) == getRightTokenIndex(i))
					return getRightFunction(i)^(-1) ;
				return getLeftFunction(i)/getRightFunction(i) ;
			}
			case TOKEN_OPERATION_POWER:
			{
				if(getLeftTokenIndex(i) == getRightTokenIndex(i))
					return getRightFunction(i) ;
				return getLeftFunction(i)^getRightFunction(i) ;
			}
			case TOKEN_OPERATION_ABS:
				return f_abs( getLeftFunction(i) ) ;
			case TOKEN_OPERATION_COS:
				return f_cos( getLeftFunction(i) ) ;
			case TOKEN_OPERATION_SIN:
				return f_sin( getLeftFunction(i) ) ;
			case TOKEN_OPERATION_TAN:
				return f_tan( getLeftFunction(i) ) ;
			case TOKEN_OPERATION_COSH:
				return f_cosh( getLeftFunction(i) ) ;
			case TOKEN_OPERATION_SINH:
				return f_sinh( getLeftFunction(i) ) ;
			case TOKEN_OPERATION_TANH:
				return f_tanh( getLeftFunction(i) ) ;
			case TOKEN_OPERATION_EXP:
				return f_exp( getLeftFunction(i) ) ;
			case TOKEN_OPERATION_SIGN:
				return f_sign( getLeftFunction(i) ) ;
			case TOKEN_OPERATION_POSITIVITY:
				return f_positivity( getLeftFunction(i) ) ;
			case TOKEN_OPERATION_NEGATIVITY:
				return f_negativity( getLeftFunction(i) ) ;
			case TOKEN_OPERATION_LOG:
				return f_log( getLeftFunction(i) ) ;
			case TOKEN_OPERATION_SQRT:
				return f_sqrt( getLeftFunction(i) ) ;
			case TOKEN_OPERATION_MIN:
				return f_min( getLeftFunction(i), getRightFunction(i) ) ;
			case TOKEN_OPERATION_MAX:
				return f_max( getLeftFunction(i), getRightFunction(i) ) ;
			default:
				return Function("0") ;
		}
	}
	return tokens[i]->getFunction() ;
		
}

Function FunctionParser::getFunction( std::string f, std::map<std::string, char> coordinates )
{
	FunctionParser expr(f, coordinates) ;
//	expr.print() ;
//	expr.printRoots() ;
	return expr.getFunction() ;
}

void FunctionParser::printRoots() const
{
	for(size_t i = 0 ; i < roots.size() ; i++)
		std::cout << roots[i] << " " ;
	std::cout << std::endl ;
}

void FunctionParser::print(bool end) const
{
	std::cout << "( " ;
	for(size_t i = 0 ; i < tokens.size() ; i++)
		tokens[i]->print(false) ;
	std::cout << " ) " ;
	if(end)
		std::cout << std::endl ;
}

void FunctionParserToken::print(bool end) const
{
	std::cout << " " << data << " " ;
	if(end)
		std::cout << std::endl ;
}

TokenOperationType FunctionParserToken::toFunctionToken() const
{
	std::vector<double> v ;
	return FunctionParserHelper::toToken( data, 0, v ).first ;
}
