#include "dual_behaviour.h"
#include "fracturecriteria/fracturecriterion.h"

using namespace Mu ;

BimaterialInterface::BimaterialInterface(Geometry * in, Form * inbehaviour, Form * outbehaviour) : LinearForm(Matrix(outbehaviour->param), true, true, outbehaviour->getNumberOfDegreesOfFreedom()), inGeometry(in),inBehaviour(inbehaviour), outBehaviour(outbehaviour)  { }

BimaterialInterface::~BimaterialInterface() { } ;

void BimaterialInterface::transform(const Function & x, const Function & y)
{
	xtransform = x ;
	xtransform.compile() ;
	ytransform = y ;
	ytransform.compile() ;
}

void BimaterialInterface::transform(const Function & x, const Function & y, const Function & z)
{
	xtransform = x ;
	xtransform.compile() ;
	ytransform = y ;
	ytransform.compile() ;
	ztransform = z ;
	ztransform.compile() ;
}

Matrix BimaterialInterface::getTensor(const Point & p) const
{
	VirtualMachine vm ;
	Point test = Point(vm.eval(xtransform, p.x, p.y, p.z), vm.eval(ytransform,  p.x, p.y, p.z), vm.eval(ztransform,  p.x, p.y, p.z)) ;
	
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

Vector BimaterialInterface::getImposedStress(const Point & p) const
{
	if(inGeometry->in(p))
		return inBehaviour->getImposedStress(p) ;
	return outBehaviour->getImposedStress(p) ;
}

void BimaterialInterface::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const
{
	bool allin = true ;
	bool allout = true ;
	Vector x = vm->eval(xtransform,gp) ;
	Vector y = vm->eval(ytransform,gp) ;
	Vector z = vm->eval(ztransform,gp) ;
	std::valarray<bool> inIn(false, x.size()) ;
	int inCount = 0;
	for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
	{
		if(inGeometry->in(Point(x[i],y[i],z[i])))
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
		inBehaviour->apply(p_i, p_j, gp, Jinv,ret,vm) ;
		return ;
	}
	else if(allout)
	{
		outBehaviour->apply(p_i, p_j, gp, Jinv,ret,vm) ;
		return ;
	}
	
	std::valarray<std::pair<Point, double> > inArray(inCount) ;
	std::valarray<Matrix> inMatrixArray( Matrix(Jinv[0].numRows(), Jinv[0].numCols()), inCount) ;
	std::valarray<std::pair<Point, double> > outArray(gp.gaussPoints.size()-inCount) ;
	std::valarray<Matrix> outMatrixArray(Matrix(Jinv[0].numRows(), Jinv[0].numCols()) , gp.gaussPoints.size()-inCount) ;
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
	Matrix retIn(ret) ;
	inBehaviour->apply(p_i, p_j, gpIn, inMatrixArray, retIn,vm) ; 
	outBehaviour->apply(p_i, p_j, gpOut, outMatrixArray,ret,vm) ;
	ret += retIn ;
}

bool BimaterialInterface::fractured() const
{
	return inBehaviour->fractured() || outBehaviour->fractured();
}

Form * BimaterialInterface::getCopy() const 
{
	return new BimaterialInterface(*this) ;
}

std::vector<BoundaryCondition * > BimaterialInterface::getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const
{
	std::vector<BoundaryCondition * > ret ;

	Vector x = VirtualMachine().eval(xtransform,gp) ;
	Vector y = VirtualMachine().eval(ytransform,gp) ;
	Vector z = VirtualMachine().eval(ztransform,gp) ;
	std::valarray<bool> inIn(false, x.size()) ;
	int inCount = 0;
	for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
	{
		if(inGeometry->in(Point(x[i],y[i],z[i])))
		{
			inIn[i] = true ;
			inCount++ ;
		}
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
			gpIn.gaussPoints[inIterator] = gp.gaussPoints[i] ;
			inIterator++ ;
		}
		else
		{
			outMatrixArray[outIterator] = Jinv[i] ;
			gpOut.gaussPoints[outIterator] = gp.gaussPoints[i] ;
			outIterator++ ;
		}
	}


	std::vector<BoundaryCondition * > temp = inBehaviour->getBoundaryConditions(s,id, p_i, gpIn, inMatrixArray) ;
	ret.insert(ret.end(), temp.begin(), temp.end()) ;
	temp = outBehaviour->getBoundaryConditions(s,id, p_i, gpOut, outMatrixArray) ;
	ret.insert(ret.end(), temp.begin(), temp.end()) ;

	return ret ;
}


void BimaterialInterface::step(double timestep, ElementState & currentState)
{
	inBehaviour->step(timestep, currentState) ;
	if(inBehaviour->getFractureCriterion())
		inBehaviour->getFractureCriterion()->step(currentState);
	
	outBehaviour->step(timestep, currentState) ;
	if(outBehaviour->getFractureCriterion())
		outBehaviour->getFractureCriterion()->step(currentState);
}

void BimaterialInterface::artificialDamageStep(double d)
{
	inBehaviour->artificialDamageStep(d) ;
	outBehaviour->artificialDamageStep(d) ;
}

FractureCriterion * BimaterialInterface::getFractureCriterion() const
{
	double max = -1 ;
	int ret = 0 ;

	double inScore = 0. ;
	FractureCriterion * inCriterion = inBehaviour->getFractureCriterion() ;
	if(inCriterion)
	{
		max = inCriterion->getSteppedScore() ;
		ret = 1 ;
	}

	double outScore = 0. ;
	FractureCriterion * outCriterion = outBehaviour->getFractureCriterion() ;
	if(outCriterion)
	{
		outScore = outCriterion->getSteppedScore() ;
		if(outScore > max || inCriterion == NULL)
			ret = 2 ;
	}
		
	switch(ret)
	{
	case 0:
		return NULL ;
	case 1:
		return inCriterion ;
	case 2:
		return outCriterion ;
	}
	return NULL ;
}









