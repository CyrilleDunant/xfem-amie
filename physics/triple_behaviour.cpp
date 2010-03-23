#include "triple_behaviour.h"

using namespace Mu ;

TrimaterialInterface::TrimaterialInterface(Geometry * in,Geometry * out, Form * inbehaviour, Form * midbehaviour, Form * outbehaviour) : LinearForm(Matrix(), true, true, 2), inGeometry(in), outGeometry(out),inBehaviour(inbehaviour),midBehaviour(midbehaviour), outBehaviour(outbehaviour)  { }

TrimaterialInterface::~TrimaterialInterface() { } ;

void TrimaterialInterface::transform(const Function & x, const Function & y)
{
	xtransform = x ;
	xtransform.compile() ;
	ytransform = y ;
	ytransform.compile() ;
}

Matrix TrimaterialInterface::getTensor(const Point & p) const
{
	VirtualMachine vm ;
	Point test = Point(vm.eval(xtransform, p), vm.eval(ytransform, p)) ;
	
	if(inGeometry->in(test))
		return inBehaviour->getTensor(p) ;
	else if(outGeometry->in(test))
		return midBehaviour->getTensor(p) ;
	
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

Matrix TrimaterialInterface::apply(const Function & p_i, const Function & p_j, const IntegrableEntity *e) const
{
	VirtualMachine vm ;
	GaussPointArray gp = e->getGaussPoints();
	bool allin = true ;
	bool allout = true ;
	bool allmid = true ;
	Vector x = vm.eval(xtransform,gp) ;
	Vector y = vm.eval(ytransform,gp) ;
	std::valarray<bool> inIn(false, x.size()) ;
	std::valarray<bool> inMid(false, x.size()) ;
	int inCount = 0;
	int midCount = 0;
	for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
	{
		if(inGeometry->in(Point(x[i],y[i])))
		{
			allout = false ;
			allmid = false ;
			inIn[i] = true ;
			inCount++ ;
		}
		else if (outGeometry->in(Point(x[i],y[i])))
		{
			allout = false ;
			inMid[i] = true ;
			midCount++ ;
		}
		else
		{
			allmid = false ;
			allin = false ;
		}
			
	}
	
	if(allin)
		return inBehaviour->apply(p_i, p_j, e) ;
	else if(allout)
		return outBehaviour->apply(p_i, p_j, e) ;
	else if(allmid)
		return midBehaviour->apply(p_i, p_j, e) ;

	std::valarray<Matrix> Jinv(Matrix(), gp.gaussPoints.size()) ;
	for(size_t i = 0 ; i < gp.gaussPoints.size() ;  i++)
	{
		e->getInverseJacobianMatrix( gp.gaussPoints[i].first, Jinv[i] ) ;
	}

	std::valarray<std::pair<Point, double> > inArray(inCount) ;
	std::valarray<Matrix> inMatrixArray(inCount) ;
	std::valarray<std::pair<Point, double> > midArray(midCount) ;
	std::valarray<Matrix> midMatrixArray(midCount) ;
	std::valarray<std::pair<Point, double> > outArray(gp.gaussPoints.size()-inCount-midCount) ;
	std::valarray<Matrix> outMatrixArray(gp.gaussPoints.size()-inCount-midCount) ;
	GaussPointArray gpIn(inArray, -1) ;
	GaussPointArray gpMid(midArray, -1) ;
	GaussPointArray gpOut(outArray, -1) ;
	
	int outIterator = 0;
	int inIterator = 0 ;
	int midIterator = 0 ;
	for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
	{
		if(inIn[i])
		{
			inMatrixArray[inIterator] = Jinv[i] ;
			gpIn.gaussPoints[inIterator++] = gp.gaussPoints[i] ;
		}
		else if(inMid[i])
		{
			midMatrixArray[midIterator] = Jinv[i] ;
			gpMid.gaussPoints[midIterator++] = gp.gaussPoints[i] ;
		}
		else
		{
			outMatrixArray[outIterator] = Jinv[i] ;
			gpOut.gaussPoints[outIterator++] = gp.gaussPoints[i] ;
		}
	}
	
	Matrix retIn(inBehaviour->getNumberOfDegreesOfFreedom(), inBehaviour->getNumberOfDegreesOfFreedom()) ;
	Matrix retMid(midBehaviour->getNumberOfDegreesOfFreedom(), midBehaviour->getNumberOfDegreesOfFreedom()) ;
	Matrix retOut(outBehaviour->getNumberOfDegreesOfFreedom(), outBehaviour->getNumberOfDegreesOfFreedom()) ;
	inBehaviour->apply(p_i, p_j, gpIn, inMatrixArray, retIn, &vm) ;
	midBehaviour->apply(p_i, p_j, gpMid, midMatrixArray, retMid, &vm) ;
	outBehaviour->apply(p_i, p_j, gpOut, outMatrixArray, retOut, &vm) ;

	return retIn+retMid+retOut ;

}

Vector TrimaterialInterface::getImposedStress(const Point & p) const
{
	if(inGeometry->in(p))
		return inBehaviour->getImposedStress(p) ;
	return outBehaviour->getImposedStress(p) ;
}

void TrimaterialInterface::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix &ret, VirtualMachine * vm) const
{
	bool allin = true ;
	bool allout = true ;
	bool allmid = true ;
	Vector x = vm->eval(xtransform,gp) ;
	Vector y = vm->eval(ytransform,gp) ;
	std::valarray<bool> inIn(false, x.size()) ;
	std::valarray<bool> inMid(false, x.size()) ;
	int inCount = 0;
	int midCount = 0;
	for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
	{
		Point test(x[i],y[i]) ;
		if(inGeometry->in(test))
		{
			allout = false ;
			allmid = false ;
			inIn[i] = true ;
			inCount++ ;
		}
		else if (outGeometry->in(test))
		{
			allout = false ;
			inMid[i] = true ;
			midCount++ ;
		}
		else
		{
			allmid = false ;
			allin = false ;
		}
		
	}
	
	if(allin)
	{
		inBehaviour->apply(p_i, p_j, gp, Jinv,ret,vm) ;
		return ;
	}
	else if(allmid)
	{
		midBehaviour->apply(p_i, p_j, gp, Jinv,ret,vm) ;
		return ;
	}
	else if(allout)
	{
		outBehaviour->apply(p_i, p_j, gp, Jinv,ret,vm) ;
		return ;
	}
	
	std::valarray<std::pair<Point, double> > inArray(inCount) ;
	std::valarray<Matrix> inMatrixArray(inCount) ;
	std::valarray<std::pair<Point, double> > midArray(midCount) ;
	std::valarray<Matrix> midMatrixArray(midCount) ;
	std::valarray<std::pair<Point, double> > outArray(gp.gaussPoints.size()-inCount-midCount) ;
	std::valarray<Matrix> outMatrixArray(gp.gaussPoints.size()-inCount-midCount) ;
	GaussPointArray gpIn(inArray, -1) ;
	GaussPointArray gpMid(midArray, -1) ;
	GaussPointArray gpOut(outArray, -1) ;
	
	int outIterator = 0;
	int inIterator = 0 ;
	int midIterator = 0 ;
	for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
	{
		if(inIn[i])
		{
			inMatrixArray[inIterator] = Jinv[i] ;
			gpIn.gaussPoints[inIterator++] = gp.gaussPoints[i] ;
		}
		else if(inMid[i])
		{
			midMatrixArray[midIterator] = Jinv[i] ;
			gpMid.gaussPoints[midIterator++] = gp.gaussPoints[i] ;
		}
		else
		{
			outMatrixArray[outIterator] = Jinv[i] ;
			gpOut.gaussPoints[outIterator++] = gp.gaussPoints[i] ;
		}
	}
	
	Matrix temp(ret) ;
	if(inMatrixArray.size())
	{
		inBehaviour->apply(p_i, p_j, gpIn, inMatrixArray,temp,vm) ;
		ret+=temp ;
	}
	if(midMatrixArray.size())
	{
		midBehaviour->apply(p_i, p_j, gpMid, midMatrixArray,temp,vm) ;
		ret+=temp ;
	}
	if(outMatrixArray.size())
	{
		outBehaviour->apply(p_i, p_j, gpOut, outMatrixArray,temp,vm) ;
		ret +=temp ;
	}


}

bool TrimaterialInterface::fractured() const
{
	return false;
}

Form * TrimaterialInterface::getCopy() const 
{
	return new TrimaterialInterface(*this) ;
}


void TrimaterialInterface::getForces(const ElementState & s, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Vector & f) const 
{
	bool allin = true ;
	bool allout = true ;
	bool allmid = true ;
	Vector x = VirtualMachine().eval(xtransform,gp) ;
	Vector y = VirtualMachine().eval(ytransform,gp) ;
	std::valarray<bool> inIn(false, x.size()) ;
	std::valarray<bool> inMid(false, x.size()) ;
	int inCount = 0;
	int midCount = 0;
	for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
	{
		Point test(x[i],y[i]) ;
		if(inGeometry->in(test))
		{
			allout = false ;
			allmid = false ;
			inIn[i] = true ;
			inCount++ ;
		}
		else if (outGeometry->in(test))
		{
			allout = false ;
			inMid[i] = true ;
			midCount++ ;
		}
		else
		{
			allmid = false ;
			allin = false ;
		}
		
	}
	
	if(allin && inBehaviour->hasInducedForces())
	{
		inBehaviour->getForces(s, p_i, gp, Jinv, f) ; return ;
	}
	else if(allout && outBehaviour->hasInducedForces())
	{
		outBehaviour->getForces(s, p_i, gp, Jinv, f) ; return ;
	}
	else if(allmid&& midBehaviour->hasInducedForces())
	{
		midBehaviour->getForces(s, p_i, gp, Jinv, f) ; return ;
	}

	std::valarray<std::pair<Point, double> > inArray(inCount) ;
	std::valarray<Matrix> inMatrixArray(inCount) ;
	std::valarray<std::pair<Point, double> > midArray(midCount) ;
	std::valarray<Matrix> midMatrixArray(midCount) ;
	std::valarray<std::pair<Point, double> > outArray(gp.gaussPoints.size()-inCount-midCount) ;
	std::valarray<Matrix> outMatrixArray(gp.gaussPoints.size()-inCount-midCount) ;
	GaussPointArray gpIn(inArray, -1) ;
	GaussPointArray gpMid(midArray, -1) ;
	GaussPointArray gpOut(outArray, -1) ;
	
	int outIterator = 0;
	int inIterator = 0 ;
	int midIterator = 0 ;
	for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
	{
		if(inIn[i])
		{
			inMatrixArray[inIterator] = Jinv[i] ;
			gpIn.gaussPoints[inIterator++] = gp.gaussPoints[i] ;
		}
		else if(inMid[i])
		{
			midMatrixArray[midIterator] = Jinv[i] ;
			gpMid.gaussPoints[midIterator++] = gp.gaussPoints[i] ;
		}
		else
		{
			outMatrixArray[outIterator] = Jinv[i] ;
			gpOut.gaussPoints[outIterator++] = gp.gaussPoints[i] ;
		}
	}
	Vector temp(f) ;
	if(inBehaviour->hasInducedForces() && outBehaviour->hasInducedForces()&& midBehaviour->hasInducedForces())
	{
		inBehaviour->getForces(s, p_i, gpIn, inMatrixArray, f) ;
		midBehaviour->getForces(s, p_i, gpMid, midMatrixArray, temp) ;f += temp ;
		outBehaviour->getForces(s, p_i, gpOut, outMatrixArray, temp) ;f += temp ;
		return ;
	}
	else if(inBehaviour->hasInducedForces() && outBehaviour->hasInducedForces())
	{
		inBehaviour->getForces(s, p_i, gpIn, inMatrixArray, f) ;
		outBehaviour->getForces(s, p_i, gpOut, outMatrixArray, temp) ;f += temp ;
		return ;
	}
	else if(inBehaviour->hasInducedForces() && midBehaviour->hasInducedForces())
	{
		inBehaviour->getForces(s, p_i, gpIn, inMatrixArray, f) ;
		midBehaviour->getForces(s, p_i, gpMid, midMatrixArray, temp) ;f += temp ;
		return  ;
	}
	else if(outBehaviour->hasInducedForces() && midBehaviour->hasInducedForces())
	{
		midBehaviour->getForces(s, p_i, gpMid, midMatrixArray, f) ;
		outBehaviour->getForces(s, p_i, gpOut, outMatrixArray, temp) ;	f += temp ;
	return ;
	}
	else if(inBehaviour->hasInducedForces())
	{
		inBehaviour->getForces(s, p_i, gpIn, inMatrixArray, f) ;
		return ;
	}
	else if(midBehaviour->hasInducedForces())
	{
		midBehaviour->getForces(s, p_i, gpMid, midMatrixArray, f) ;
		return ;
	}
	else if(outBehaviour->hasInducedForces())
	{
		outBehaviour->getForces(s, p_i, gpOut, outMatrixArray, f) ;
		return ;
	}
}

void TrimaterialInterface::step(double timestep, ElementState & currentState)
{
	inBehaviour->step(timestep, currentState) ;
	midBehaviour->step(timestep, currentState) ;
	outBehaviour->step(timestep, currentState) ;
}

void TrimaterialInterface::artificialDamageStep(double d)
{
	inBehaviour->artificialDamageStep(d) ;
	midBehaviour->artificialDamageStep(d) ;
	outBehaviour->artificialDamageStep(d) ;
}


bool TrimaterialInterface::hasInducedForces() const
{
	return inBehaviour->hasInducedForces() || outBehaviour->hasInducedForces()|| midBehaviour->hasInducedForces() ;
}








