#include "dual_behaviour.h"

using namespace Mu ;

BimaterialInterface::BimaterialInterface(Geometry * in, Form * inbehaviour, Form * outbehaviour) : LinearForm(Matrix(outbehaviour->param), true, true, inbehaviour->getNumberOfDegreesOfFreedom()), inGeometry(in),inBehaviour(inbehaviour), outBehaviour(outbehaviour)  { }

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
	Point test = Point(vm.eval(xtransform, p), vm.eval(ytransform, p), vm.eval(ztransform, p)) ;
	
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
	Vector z = VirtualMachine().eval(ztransform,gp) ;
	std::valarray<bool> inIn(false, x.size()) ;
	size_t inCount = 0;
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
		return inBehaviour->apply(p_i, p_j, e) ;
	else if(allout)
		return outBehaviour->apply(p_i, p_j, e) ;
	else
	{
		 std::valarray<Matrix> Jinv(Matrix(), gp.gaussPoints.size()) ;
		for(size_t i = 0 ; i < gp.gaussPoints.size() ;  i++)
		{
			 e->getInverseJacobianMatrix( gp.gaussPoints[i].first, Jinv[i]) ;
		}

		std::valarray<std::pair<Point, double> > inArray(inCount) ;
		std::valarray<Matrix> inMatrixArray(Matrix(Jinv[0].numRows(), Jinv[0].numCols()), inCount ) ;
		std::valarray<std::pair<Point, double> > outArray(gp.gaussPoints.size()-inCount) ;
		std::valarray<Matrix> outMatrixArray( Matrix(Jinv[0].numRows(), Jinv[0].numCols()), gp.gaussPoints.size()-inCount) ;
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
		VirtualMachine vm ;
		Matrix retIn(inBehaviour->getNumberOfDegreesOfFreedom(), inBehaviour->getNumberOfDegreesOfFreedom()) ;
		Matrix retOut(inBehaviour->getNumberOfDegreesOfFreedom(), inBehaviour->getNumberOfDegreesOfFreedom()) ;
		inBehaviour->apply(p_i, p_j, gpIn, inMatrixArray, retIn, &vm ) ; 
		outBehaviour->apply(p_i, p_j, gpOut, outMatrixArray, retOut, &vm) ;
		return retIn+retOut ;
	}
	//shut up the compiler
	return Matrix() ;
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


void BimaterialInterface::getForces(const ElementState & s, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Vector & f) const 
{
	bool allin = true ;
	bool allout = true ;
	Vector x = VirtualMachine().eval(xtransform,gp) ;
	Vector y = VirtualMachine().eval(ytransform,gp) ;
	Vector z = VirtualMachine().eval(ztransform,gp) ;
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

	if(allin && inBehaviour->hasInducedForces())
	{
		inBehaviour->getForces(s, p_i, gp, Jinv, f) ; 
		return ;
	}
	else if(allout && outBehaviour->hasInducedForces())
	{
		outBehaviour->getForces(s, p_i, gp, Jinv, f) ;
		return ;
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

	if(inBehaviour->hasInducedForces() && outBehaviour->hasInducedForces())
	{
		Vector temp(f) ;
		inBehaviour->getForces(s, p_i, gpIn, inMatrixArray, f) ;
		outBehaviour->getForces(s, p_i, gpOut, outMatrixArray, temp) ;f += temp ;
		return ;
	}
	else if(inBehaviour->hasInducedForces())
	{
		inBehaviour->getForces(s, p_i, gpIn, inMatrixArray, f) ;
		return ;
	}
	else if(outBehaviour->hasInducedForces())
	{
		outBehaviour->getForces(s, p_i, gpOut, outMatrixArray, f) ;
		return  ;
	}
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








