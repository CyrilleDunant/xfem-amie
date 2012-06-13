// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011

#include "triple_behaviour.h"
#include "fracturecriteria/fracturecriterion.h"
#include "damagemodels/damagemodel.h"

using namespace Mu ;

TrimaterialInterface::TrimaterialInterface(Geometry * in,Geometry * out, Form * inbehaviour, Form * midbehaviour, Form * outbehaviour) : LinearForm(Matrix(3,3), true, true, 2), inGeometry(in), outGeometry(out),inBehaviour(inbehaviour),midBehaviour(midbehaviour), outBehaviour(outbehaviour)  { }

TrimaterialInterface::~TrimaterialInterface() { } ;

void TrimaterialInterface::transform(const Function & x, const Function & y)
{
	xtransform = x ;
	xtransform.compile() ;
	ytransform = y ;
	ytransform.compile() ;
}

Matrix TrimaterialInterface::getTensor(const Point & p, IntegrableEntity * e) const
{
	VirtualMachine vm ;
	Point test = Point(vm.eval(xtransform, p), vm.eval(ytransform, p)) ;
	
	if(inGeometry->in(test))
		return inBehaviour->getTensor(p, e) ;
	else if(outGeometry->in(test))
		return midBehaviour->getTensor(p, e) ;
	
	return outBehaviour->getTensor(p, e) ;
	
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

Vector TrimaterialInterface::getImposedStress(const Point & p, IntegrableEntity * e) const
{
	if(inGeometry->in(p))
		return inBehaviour->getImposedStress(p, e) ;
	return outBehaviour->getImposedStress(p, e) ;
}

Vector TrimaterialInterface::getImposedStrain(const Point & p, IntegrableEntity * e) const
{
	if(inGeometry->in(p))
		return inBehaviour->getImposedStrain(p,e) ;
	return outBehaviour->getImposedStrain(p,e) ;
}

bool TrimaterialInterface::changed() const { return false ;}

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
		return inBehaviour->fractured() || outBehaviour->fractured() || midBehaviour->fractured();
}

Form * TrimaterialInterface::getCopy() const 
{
	return new TrimaterialInterface(*this) ;
}

std::vector<BoundaryCondition * > TrimaterialInterface::getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const
{
	std::vector<BoundaryCondition * > ret ;

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
			inIn[i] = true ;
			inCount++ ;
		}
		else if (outGeometry->in(test))
		{
			inMid[i] = true ;
			midCount++ ;
		}
		
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


	std::vector<BoundaryCondition * > temp = inBehaviour->getBoundaryConditions(s,id, p_i, gpIn, inMatrixArray) ;
	ret.insert(ret.end(), temp.begin(), temp.end()) ;
	temp = midBehaviour->getBoundaryConditions(s,id, p_i, gpMid, midMatrixArray) ;
	ret.insert(ret.end(), temp.begin(), temp.end()) ;
	temp = outBehaviour->getBoundaryConditions(s,id, p_i, gpOut, outMatrixArray) ;
	ret.insert(ret.end(), temp.begin(), temp.end()) ;

	return ret ;
}

void TrimaterialInterface::step(double timestep, ElementState & currentState)
{
	inBehaviour->step(timestep, currentState) ;
	midBehaviour->step(timestep, currentState) ;
	outBehaviour->step(timestep, currentState) ;
}

DamageModel * TrimaterialInterface::getDamageModel() const
{
	double max = -2 ;
	int ret = 0 ;

	double inScore = 0. ;
	DamageModel * inCriterion = inBehaviour->getDamageModel() ;
	DamageModel * midCriterion = midBehaviour->getDamageModel() ;
	DamageModel * outCriterion = outBehaviour->getDamageModel() ;
	if(inCriterion)
	{
		max = inCriterion->getState().max() ;
		ret = 1 ;
	}

	double midScore = 0. ;
	
	if(midCriterion)
	{
		midScore = midCriterion->getState().max() ;
		if(midScore > max || (inCriterion == NULL && outCriterion == NULL))
		{
			max = midScore ;
			ret = 2 ;
		}
	}
	
	double outScore = 0. ;
	
	if(outCriterion)
	{
		outScore = outCriterion->getState().max() ;
		if(outScore > max || (inCriterion == NULL && midCriterion == NULL))
			ret = 3 ;
	}
		
	switch(ret)
	{
	case 0:
		return NULL ;
	case 1:
		return inCriterion ;
	case 2:
		return midCriterion ;
	case 3:
		return outCriterion ;
	}
	return NULL ;
}

FractureCriterion * TrimaterialInterface::getFractureCriterion() const
{
	double max = -2 ;
	int ret = 0 ;

	double inScore = 0. ;
	FractureCriterion * inCriterion = inBehaviour->getFractureCriterion() ;
	FractureCriterion * midCriterion = midBehaviour->getFractureCriterion() ;
	FractureCriterion * outCriterion = outBehaviour->getFractureCriterion() ;
	if(inCriterion)
	{
		max = inCriterion->getScoreAtState() ;
		ret = 1 ;
	}

	double midScore = 0. ;
	
	if(midCriterion)
	{
		midScore = midCriterion->getScoreAtState() ;
		if(midScore > max || (inCriterion == NULL && outCriterion == NULL))
		{
			max = midScore ;
			ret = 2 ;
		}
	}
	
	double outScore = 0. ;
	
	if(outCriterion)
	{
		outScore = outCriterion->getScoreAtState() ;
		if(outScore > max || (inCriterion == NULL && midCriterion == NULL))
			ret = 3 ;
	}
		
	switch(ret)
	{
	case 0:
		return NULL ;
	case 1:
		return inCriterion ;
	case 2:
		return midCriterion ;
	case 3:
		return outCriterion ;
	}
	return NULL ;
}







