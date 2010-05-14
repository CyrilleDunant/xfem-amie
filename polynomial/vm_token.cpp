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
#include "../elements/elements.h"

using namespace Mu ;



CurvilinearXOperatorToken::CurvilinearXOperatorToken(const SegmentedLine * l, bool fromHead) : Token(false, std::make_pair(std::make_pair(TOKEN_CURVILINEAR_X_OPERATOR, 0), (double)(0))), 
																																				line(l), 
																																				origin(fromHead)
	{
	} ;
	
void CurvilinearXOperatorToken::eval(Context & context) const
{
	Point p(*context.memory.top_pos, *context.memory.prev_top_pos) ;
	
	if(origin)
	{

		Point p_(p) ;
		line->project(&p_) ;
		
		if(p_ == *line->getHead())
		{
			Line l(p_, *line->getHead()-line->getBoundingPoint(1)) ;
			p_ = l.projection(p) ;
			context.memory.pop_back() ;
			*context.memory.top_pos =  sqrt(squareDist2D(&p_,line->getHead()));
			return ;
		}
		p_ = p ;
		line->project(&p_) ;
		if(p_ == *line->getTail())
		{
			Line l(p_, *line->getTail()-line->getBoundingPoint(line->getBoundingPoints().size()-2)) ;
			p_ = l.projection(p) ;
			context.memory.pop_back() ;
			*context.memory.top_pos = - sqrt(squareDist2D(&p_,line->getTail()));
			return ;
		}
		p_ = p ;
		line->project(&p_) ;
		double d = 0 ;
		for(size_t i = 0 ; i < line->getBoundingPoints().size()-1 ; i++)
		{
			Segment s(line->getBoundingPoint(i), line->getBoundingPoint(i+1)) ;
			if(s.on(p_))
			{
				context.memory.pop_back() ;
				*context.memory.top_pos =  -d - sqrt(squareDist2D(p_,line->getBoundingPoint(i)));
				return ;
			}
			else
				d += s.norm() ;
		}
	}
	else
	{
		Point p_(p) ;
		line->project(&p_) ;
		
		if(p_ == *line->getHead())
		{
			Line l(p_, *line->getHead()-line->getBoundingPoint(1)) ;
			p_ = l.projection(p) ;
			context.memory.pop_back() ;
			*context.memory.top_pos =  -sqrt(squareDist2D(&p_,line->getHead()));
			return ;
		}
		p_ = p ;
		line->project(&p_) ;

		if(p_ == *line->getTail())
		{

			Line l(p_, *line->getTail()-line->getBoundingPoint(line->getBoundingPoints().size()-2)) ;
			p_ = l.projection(p) ;
			context.memory.pop_back() ;
			*context.memory.top_pos =  sqrt(squareDist2D(&p_,line->getTail()));
			return ;
		}


		p_ = p ;
		line->project(&p_) ;
		double d = 0 ;
		for(size_t i = line->getBoundingPoints().size()-1 ; i > 0 ; i--)
		{
			Segment s(line->getBoundingPoint(i), line->getBoundingPoint(i-1)) ;
			if(s.on(p_))
			{
				context.memory.pop_back() ;
				*context.memory.top_pos =  - d - sqrt(squareDist2D(p_,line->getBoundingPoint(i)));
				return ;
			}
			else
				d += s.norm() ;
		}
	}
}
	
std::string CurvilinearXOperatorToken::print() const
{
	return std::string("s-x") ;
}


CurvilinearYOperatorToken::CurvilinearYOperatorToken(const SegmentedLine * l, bool fromHead) : Token(false, std::make_pair(std::make_pair(TOKEN_CURVILINEAR_Y_OPERATOR, 0), (double)(0))), 
																																			line(l), 
																																			origin(fromHead)
{
} ;
	
void CurvilinearYOperatorToken::eval(Context & context) const
{
	Point p( *context.memory.top_pos, *context.memory.prev_top_pos) ;

	Point p_(p) ;
	line->project(&p_) ;
	Point o(*line->getHead()) ;
	if(!origin)
		o = *line->getTail() ;

	double s = 1 ;
	if(p_ == *line->getHead())
	{
		Line l(*line->getHead(), *line->getHead()-line->getBoundingPoint(1)) ;
		p_ = l.projection(p) ;
		if(atan2(o.y-p.y, o.x-p.x) < atan2(o.y-p_.y, o.x-p_.x))
			s = -1 ;
		
		context.memory.pop_back() ;
		*context.memory.top_pos =  s * sqrt(squareDist2D(p_,p));
		return ;
	}
	else if(p_ == *line->getTail())
	{
// 		std::cout << "pluf" << std::endl ;
		Line l(*line->getTail(), *line->getTail()-line->getBoundingPoint(line->getBoundingPoints().size()-2)) ;
		p_ = l.projection(p) ;
		s = -1 ;
		if(atan2(p.y-o.y, p.x-o.x) > atan2(p_.y-o.y, p_.x-o.x))
			s = 1 ;
		context.memory.pop_back() ;
		*context.memory.top_pos = s * sqrt(squareDist2D(p_,p));
		return ;
	}


	if(atan2(o.y-p.y, o.x-p.x) > atan2(o.y-p_.y, o.x-p_.x))
		s = -1 ;
	context.memory.pop_back() ;
	*context.memory.top_pos = s * sqrt(squareDist2D(p_,p));

}
	
std::string CurvilinearYOperatorToken::print() const
{
	return std::string("s-y") ;
}


ConstantToken::ConstantToken(double v , bool nul) : Token(nul, std::make_pair(std::make_pair(TOKEN_CONSTANT, 0), (double)(v))){ };
	
ConstantToken::~ConstantToken() { } ;

NullToken::~NullToken() { }

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



Transform2DToken::Transform2DToken(const ElementarySurface * g ) : Token(false, std::make_pair(std::make_pair(TOKEN_2D_TRANSFORM, 0),
	(double)(0))), e(g)
{
}

void Transform2DToken::eval(Context & context) const
{
	Point p = coordinateTransform(Point(context.x,context.y), e->getBoundingPoints(), e->getShapeFunctions()) ;
	context.memory.push_back(p.y);
	context.memory.push_back(p.x) ;
}
Transform2DToken::~Transform2DToken() { };
std::string Transform2DToken::print() const
{
	return std::string("transform2D") ;
}



Transform3DToken::Transform3DToken(const ElementaryVolume * g ) : Token(false, std::make_pair(std::make_pair(TOKEN_3D_TRANSFORM, 0),
	(double)(0))), e(g)
{
}

void Transform3DToken::eval(Context & context) const
{
	Point p = coordinateTransform(Point(context.x,context.y,context.z), e->getBoundingPoints(), e->getShapeFunctions()) ;
	context.memory.push_back(p.z) ;
	context.memory.push_back(p.y) ;
	context.memory.push_back(p.x) ;
}

Transform3DToken::~Transform3DToken() { };
std::string Transform3DToken::print() const
{
	return std::string("transform3D") ;
}

double interpolate(const double a, const double b)
{	
	return a/(b + a) ;
}


PositionOperatorToken::PositionOperatorToken(Segment s_ ) : Token(false, std::make_pair(std::make_pair(TOKEN_POSITION_OPERATOR, 0), (double)(0)))
{
	Point vector(-s_.vector().y, s_.vector().x) ;
	w= s_.midPoint()+vector*0.00001 ;
	s.push_back(s_)  ;
}

PositionOperatorToken::PositionOperatorToken(std::vector<Segment> s_ ) : Token(false)
{
	
	for(size_t i = 0 ; i < s_.size() ; i++)
	{

		s.push_back(s_[i])  ;
	}
	
	Point vector(-s[0].vector().y, s[0].vector().x) ;
	w= s[0].midPoint()+vector*0.00001 ;
}

void PositionOperatorToken::eval(Context & context) const
{
	if(s.empty())
	{
		context.memory.pop_back() ;
		*context.memory.top_pos = 0 ;
		return ;
	}
	

	Point test(*context.memory.top_pos, *context.memory.prev_top_pos) ;
	context.memory.pop_back() ;
	int intersections = 0 ;
	for(size_t i = 0 ; i < s.size() ; i++)
	{
		if(s[i].on(test))
		{
			*context.memory.top_pos = 0 ;
			return ;
		}
		if(s[i].intersects(test, w))
		{
			intersections++ ;
		}
		
	}
	
	if((intersections % 2)  != 0)
		*context.memory.top_pos =  -1 ;
	else
		*context.memory.top_pos = 1 ;

}

PositionOperatorToken::~PositionOperatorToken() { };
std::string PositionOperatorToken::print() const
{
	return std::string("positionOp") ;
}



Atan2Token::Atan2Token() : Token(false, std::make_pair(std::make_pair(TOKEN_ATAN2, 0),(double)(0)))  { };

void Atan2Token::eval(Context & context) const
{

// 	double y = *context.memory.prev_top_pos ;
// 	double x = *context.memory.top_pos ;
// 	context.memory.pop_back() ;
// 	if (x>=0)
// 	{
// 		double r = (x - std::abs(y)) / (x + std::abs(y));
// 		*context.memory.top_pos = 0.1963 * r*r*r - 0.9817 * r + M_PI*.25;
// 	}
// 	else
// 	{
// 		double r = (x + std::abs(y)) / (std::abs(y) - x);
// 		*context.memory.top_pos = 0.1963 * r*r*r - 0.9817 * r + M_PI*.75 ;
// 	}
// 	if (y < 0)
// 		*context.memory.top_pos = -(*context.memory.top_pos);     // negate if in quad III or IV
// 	
	double new_val = atan2(*context.memory.prev_top_pos, *context.memory.top_pos) ;
	
	context.memory.pop_back() ;
	*context.memory.top_pos = new_val ;
}

Atan2Token::~Atan2Token() {} ;

std::string Atan2Token::print() const
{
	return std::string("atan2")  ;
}


CosToken::CosToken() : Token(false, std::make_pair(std::make_pair(TOKEN_COS, 0),(double)(0))) { };

void CosToken::eval(Context & context) const
{
// 	const double B = 4. * M_1_PI;
// 	const double C = -4. * M_1_PI * M_1_PI;
// 	const double theta = *context.memory.top_pos + M_PI_2 ;
// 	double y = B * theta + C * theta * abs(theta);
// 	
// //  const float Q = 0.775;
// 	const double P = 0.225;
// 	
// 	*context.memory.top_pos = P * (y * abs(y) - y) + y;   // Q * y + P * y * abs(y)
	*context.memory.top_pos = cos(*context.memory.top_pos) ;
}

CosToken::~CosToken() {} ;

std::string CosToken::print() const
{
	return std::string("cos")  ;
}


SinToken::SinToken() : Token(false, std::make_pair(std::make_pair(TOKEN_SIN, 0),(double)(0)) ) { };

void SinToken::eval(Context & context) const
{
	
// 	const double B = 4. * M_1_PI;
// 	const double C = -4. * M_1_PI * M_1_PI;
// 	
// 	double y = B * (*context.memory.top_pos) + C * (*context.memory.top_pos) * abs((*context.memory.top_pos));
// 
// //  const float Q = 0.775;
// 	const double P = 0.225;
// 		
// 	*context.memory.top_pos = P * (y * abs(y) - y) + y;   // Q * y + P * y * abs(y)
// 	
// 		
	*context.memory.top_pos = sin(*context.memory.top_pos) ;
}

SinToken::~SinToken() {} ;

std::string SinToken::print() const
{
	return std::string("sin")  ;
}


