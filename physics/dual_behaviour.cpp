#include "dual_behaviour.h"

using namespace Mu ;

BimaterialInterface::BimaterialInterface(Geometry * in, Form * inbehaviour, Form * outbehaviour) : LinearForm(Matrix(), true, true, 2), inGeometry(in),inBehaviour(inbehaviour), outBehaviour(outbehaviour)  { }

BimaterialInterface::~BimaterialInterface() { } ;

void BimaterialInterface::transform(const Function & x, const Function & y)
{
	xtransform = x ;
	xtransform.compile() ;
	ytransform = y ;
	ytransform.compile() ;
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
	GaussPointArray gp = e->getGaussPoints();
	bool allin = true ;
	bool allout = true ;
	Vector x = VirtualMachine().eval(xtransform,gp) ;
	Vector y = VirtualMachine().eval(ytransform,gp) ;
	for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
	{
		if(inGeometry->in(Point(x[i],y[i])))
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
		std::valarray<Matrix> Jinv(gp.gaussPoints.size()) ;
		for(size_t i = 0 ;i < gp.gaussPoints.size() ; i++)
		{
			Jinv[i] = e->getInverseJacobianMatrix(gp.gaussPoints[i].first) ;
		}
		GaussPointArray gpIn(gp) ;
		gpIn.id = -1 ;
		GaussPointArray gpOut(gp) ;
		gpOut.id = -1 ;
		for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
		{
			if(inGeometry->in(Point(x[i], y[i])))
				gpOut.gaussPoints[i].second = 0 ;
			else
				gpIn.gaussPoints[i].second = 0 ;
		}
		
		return inBehaviour->apply(p_i, p_j, gpIn, Jinv) + outBehaviour->apply(p_i, p_j, gpOut, Jinv) ;
		
	}
	
}

Matrix BimaterialInterface::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const
{
	bool allin = true ;
	bool allout = true ;
	Vector x = VirtualMachine().eval(xtransform,gp) ;
	Vector y = VirtualMachine().eval(ytransform,gp) ;
	std::valarray<bool> inIn(false, x.size()) ;
	int inCount = 0;
	for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
	{
		if(inGeometry->in(Point(x[i],y[i])))
		{
			allout = false ;
			inIn[i] = true ;
			inCount++ ;
		}
		else
		{
			allin = false ;
		}
			
	}

	if(allin)
	{
		return inBehaviour->apply(p_i, p_j, gp, Jinv) ;
	}
	else if(allout)
	{
		return outBehaviour->apply(p_i, p_j, gp, Jinv) ;
	}
	
	std::valarray<std::pair<Point, double> > inArray(inCount) ;
	std::valarray<Matrix> inMatrixArray(inCount) ;
	std::valarray<std::pair<Point, double> > outArray(gp.gaussPoints.size()-inCount) ;
	std::valarray<Matrix> outMatrixArray(gp.gaussPoints.size()-inCount) ;
	GaussPointArray gpIn(inArray, -1) ;
	GaussPointArray gpOut(outArray, -1) ;

	int outIterator = 0;
	int inIterator = 0 ;
	for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
	{
		if(inIn[i])
		{
			inMatrixArray[inIterator] = Jinv[i] ;
			gpIn.gaussPoints[inIterator++] = gp.gaussPoints[i] ;
		}
		else
		{
			outMatrixArray[outIterator] = Jinv[i] ;
			gpOut.gaussPoints[outIterator++] = gp.gaussPoints[i] ;
		}
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


Vector BimaterialInterface::getForces(const ElementState & s, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const 
{

	bool allin = true ;
	bool allout = true ;
	Vector x = VirtualMachine().eval(xtransform,gp) ;
	Vector y = VirtualMachine().eval(ytransform,gp) ;
	std::valarray<bool> inIn(false, x.size()) ;
	int inCount = 0;
	for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
	{
		if(inGeometry->in(Point(x[i],y[i])))
		{
			allout = false ;
			inIn[i] = true ;
			inCount++ ;
		}
		else
		{
			allin = false ;
		}
			
	}

	if(allin)
	{
		return inBehaviour->getForces(s, p_i, gp, Jinv) ;
	}
	else if(allout)
	{
		return outBehaviour->getForces(s, p_i, gp, Jinv) ;
	}

	std::valarray<std::pair<Point, double> > inArray(inCount) ;
	std::valarray<Matrix> inMatrixArray(inCount) ;
	std::valarray<std::pair<Point, double> > outArray(gp.gaussPoints.size()-inCount) ;
	std::valarray<Matrix> outMatrixArray(gp.gaussPoints.size()-inCount) ;
	GaussPointArray gpIn(inArray, -1) ;
	GaussPointArray gpOut(outArray, -1) ;

	int outIterator = 0;
	int inIterator = 0 ;
	for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
	{
		if(inIn[i])
		{
			inMatrixArray[inIterator] = Jinv[i] ;
			gpIn.gaussPoints[inIterator++] = gp.gaussPoints[i] ;
		}
		else
		{
			outMatrixArray[outIterator] = Jinv[i] ;
			gpOut.gaussPoints[outIterator++] = gp.gaussPoints[i] ;
		}
	}
	
	Vector inForces = inBehaviour->getForces(s, p_i, gpIn, inMatrixArray) ;
	Vector outForces = outBehaviour->getForces(s, p_i, gpOut, outMatrixArray) ;
	
	if(inForces.size() == outForces.size())
		return inForces+outForces ;
	else if(inForces.size() > outForces.size())
		return inForces ;
	
	return outForces ;
		
}

void BimaterialInterface::step(double timestep, ElementState & currentState)
{
	inBehaviour->step(timestep, currentState) ;
	outBehaviour->step(timestep, currentState) ;
}

bool BimaterialInterface::hasInducedForces() const
{
	return inBehaviour->hasInducedForces() || outBehaviour->hasInducedForces() ;
}








