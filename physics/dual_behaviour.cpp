#include "dual_behaviour.h"

using namespace Mu ;

BimaterialInterface::BimaterialInterface(Geometry * in, Form * inbehaviour, Form * outbehaviour) : LinearForm(Matrix(), false, true, 2), inGeometry(in),inBehaviour(inbehaviour), outBehaviour(outbehaviour)  { }

BimaterialInterface::~BimaterialInterface() { } ;

void BimaterialInterface::transform(const Function & x, const Function & y)
{
	xtransform = x ;
	ytransform = y ;
}

Matrix BimaterialInterface::getTensor(const Point & p) const
{
	VirtualMachine vm ;
	Point test = Point(vm.eval(xtransform, p), vm.eval(ytransform, p)) ;
	
	if(inGeometry->in(test))
		return inBehaviour->getTensor(p) ;
	
	return outBehaviour->getTensor(p) ;
	
// 	FunctionMatrix C(3,3) ;
// 	
// 	Function domain(inGeometry, xtransform, ytransform) ;
// 	
// 	for(size_t i = 0 ; i < 3 ; i++)
// 	{
// 		for(size_t j = 0 ; j < 3 ; j++)
// 		{
// 			C[i][j] = f_positivity(domain)*inBehaviour->getTensor(p)[i][j] + f_negativity(domain)*outBehaviour->getTensor(p)[i][j];
// 		}
// 	}
// 	return vm.eval(C, p.x, p.y) ;
}

Matrix BimaterialInterface::apply(const Function & p_i, const Function & p_j, const IntegrableEntity *e) const
{
	std::valarray< std::pair<Point,double> > gp = e->getGaussPoints();
	bool allin = true ;
	bool allout = true ;
	for(size_t i = 0 ; i < gp.size() ; i++)
	{
		if(inGeometry->in(Point(VirtualMachine().eval(xtransform,gp[i].first), VirtualMachine().eval(ytransform,gp[i].first))))
			allout = false ;
		else
			allin = false ;
			
	}
	
	if(allin)
		return inBehaviour->apply(p_i, p_j, e) ;
	else if(allout)
		return outBehaviour->apply(p_i, p_j, e) ;
	else
	{
		std::valarray<Matrix> Jinv(gp.size()) ;
		for(size_t i = 0 ;i < gp.size() ; i++)
		{
			Jinv[i] = e->getInverseJacobianMatrix(gp[i].first) ;
		}
		std::valarray< std::pair<Point,double> > gpIn(gp) ;
		std::valarray< std::pair<Point,double> > gpOut(gp) ;
		for(size_t i = 0 ; i < gp.size() ; i++)
		{
			if(inGeometry->in(Point(VirtualMachine().eval(xtransform,gp[i].first), VirtualMachine().eval(ytransform,gp[i].first))))
				gpOut[i].second = 0 ;
			else
				gpIn[i].second = 0 ;
		}
		
		return inBehaviour->apply(p_i, p_j, gpIn, Jinv) + outBehaviour->apply(p_i, p_j, gpOut, Jinv) ;
		
	}
	
}

Matrix BimaterialInterface::apply(const Function & p_i, const Function & p_j, const std::valarray< std::pair<Point,double> > &gp, const std::valarray<Matrix> &Jinv) const
{
	bool allin = true ;
	bool allout = true ;
	for(size_t i = 0 ; i < gp.size() ; i++)
	{
		if(inGeometry->in(Point(VirtualMachine().eval(xtransform,gp[i].first), VirtualMachine().eval(ytransform,gp[i].first))))
			allout = false ;
		else
			allin = false ;
		
	}

	if(allin)
	{
		return inBehaviour->apply(p_i, p_j, gp, Jinv) ;
	}
	else if(allout)
	{
		return outBehaviour->apply(p_i, p_j, gp, Jinv) ;
	}
	

	std::valarray< std::pair<Point,double> > gpIn(gp) ;
	std::valarray< std::pair<Point,double> > gpOut(gp) ;
	for(size_t i = 0 ; i < gp.size() ; i++)
	{
		if(inGeometry->in(Point(VirtualMachine().eval(xtransform,gp[i].first), VirtualMachine().eval(ytransform,gp[i].first))))
			gpOut[i].second = 0 ;
		else
			gpIn[i].second = 0 ;
	}
	
	return inBehaviour->apply(p_i, p_j, gpIn, Jinv) + outBehaviour->apply(p_i, p_j, gpOut, Jinv) ;
		

}

bool BimaterialInterface::fractured() const
{
	return false;
}

Form * BimaterialInterface::getCopy() const 
{
	return new BimaterialInterface(*this) ;
}


Vector BimaterialInterface::getForces(const ElementState * s, const Function & p_i, const Function & p_j, const std::valarray< std::pair<Point, double> > &gp, const std::valarray<Matrix> &Jinv) const 
{
	std::valarray< std::pair<Point,double> > gpIn(gp) ;
	std::valarray< std::pair<Point,double> > gpOut(gp) ;
	for(size_t i = 0 ; i < gp.size() ; i++)
	{
		if(inGeometry->in(Point(VirtualMachine().eval(xtransform,gp[i].first), VirtualMachine().eval(ytransform,gp[i].first))))
			gpOut[i].second = 0 ;
		else
			gpIn[i].second = 0 ;
	}
	
	Vector inForces = inBehaviour->getForces(s, p_i,p_j, gpIn, Jinv) ;
	Vector outForces = outBehaviour->getForces(s, p_i,p_j, gpOut, Jinv) ;
	
	if(inForces.size() == outForces.size())
		return inForces+outForces ;
	else if(inForces.size() > outForces.size())
		return inForces ;
	
	return outForces ;
		
}



bool BimaterialInterface::hasInducedForces() const
{
	return inBehaviour->hasInducedForces() || outBehaviour->hasInducedForces() ;
}








