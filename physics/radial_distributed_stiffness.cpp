//
// C++ Implementation: RadialDistributedStiffness
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
//

#include "radial_distributed_stiffness.h"
#include "../utilities/xml.h"
#include "physics_base.h"
#include "stiffness.h"


using namespace Mu ;

RadialDistributedStiffness::RadialDistributedStiffness(std::vector<std::pair<double, Matrix> > rig) : LinearForm(rig[0].second, false, false, rig[0].second.numRows()/3+1) 
{
	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() > 9)
		v.push_back(ZETA);

	stiff = rig ;
	angle = 0 ;
} ;

RadialDistributedStiffness::RadialDistributedStiffness(XMLTree * xml) : LinearForm(Matrix(3,3), false, false, 3) 
{
	param = Matrix(3,3) ;
	angle = 0 ;

	if(xml->match("radial distributed stiffness"))
	{
		for(size_t i = 0 ; i < xml->nChildren() ; i++)
		{
			stiff.push_back(std::make_pair(xml->getChild(i)->buildDouble().second,
						       xml->getChild(i)->getChild(0)->buildMatrix().second)) ;
		}
		param = stiff[0].second ;

	}

	v.push_back(XI);
	v.push_back(ETA);
	if(param.size() > 9)
		v.push_back(ZETA);

}

RadialDistributedStiffness::~RadialDistributedStiffness() { } ;

void RadialDistributedStiffness::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{
// std::cout << "a--" << std::endl ;
// Jinv[0].print() ;
// std::cout << "--b" << std::endl ;
	vm->ieval(Gradient(p_i) * param * Gradient(p_j, true), gp, Jinv,v, ret) ;
// 	ret.print() ;
}

bool RadialDistributedStiffness::fractured() const
{
	return false ;
}

XMLTree * RadialDistributedStiffness::toXML()
{
	XMLTree * rds = new XMLTree("radial distributed stiffness") ;
	std::vector<XMLTree *> s ;
	for(size_t i = 0 ; i < stiff.size() ; i++)
	{
		s.push_back(new XMLTree("angle",stiff[i].first)) ;
		s[i]->addChild(new XMLTree("stiffness",stiff[i].second)) ;
	}
	rds->addChild(s) ;
	return rds ;
}

void RadialDistributedStiffness::setAngle(double a)
{
	angle = a ;
	while(angle > 2 * M_PI)
		angle = angle - 2*M_PI ;
}

Form * RadialDistributedStiffness::getCopy() const 
{
	size_t i = 0 ;
	while(angle < stiff[i].first && i < stiff.size())
		i++ ;
	return new Stiffness(stiff[i].second) ;
}

RadialInclusion::RadialInclusion(Feature * f, double r, Point c) : Circle(r,c), Inclusion(f,r,c)
{
}

RadialInclusion::RadialInclusion(Feature * f, double r, double x, double y) : Circle(r,x,y), Inclusion(f,r,x,y)
{
}

RadialInclusion::RadialInclusion(double r, Point c) : Circle(r,c), Inclusion(r,c) 
{
}

RadialInclusion::RadialInclusion(double r, double x, double y) : Circle(r,x,y), Inclusion(r,x,y) 
{
}

RadialInclusion::RadialInclusion(XMLTree * xml) : Circle(xml->getChild(0)), Inclusion(Circle(xml->getChild(0)))
{
	this->setBehaviour(new RadialDistributedStiffness(xml->getChild(1))) ;
}

Form * RadialInclusion::getBehaviour( const Point & p )
{
	static_cast<RadialDistributedStiffness *>(behaviour)->setAngle(p.angle()) ;
	return behaviour->getCopy() ;
}		


