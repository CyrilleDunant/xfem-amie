// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
// Author: Ines Jaouadi <ines.jaouadi@epfl.ch>, (C) 2005-2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "elements.h"

using namespace Mu ;

	
ElementState * ElementarySurface::getState() const
{
	return state ;
}


ElementState * ElementaryVolume::getState() const
{
	return state ;
}


ElementarySurface::~ElementarySurface()
{
	if(isFather)
		delete this->shapefunc ;
	
	delete this->state ;

}

ElementarySurface::ElementarySurface(bool f ) : isFather(f)
{
	this->state = new ElementState(this) ;
	this->nonlinbehaviour = NULL ;
	this->behaviour = NULL ;
}

void ElementarySurface::setOrder(Order o)
{
	order = o ;
}


Form * ElementarySurface::getBehaviour() const
{
	return this->behaviour ;
}

NonLinearForm * ElementarySurface::getNonLinearBehaviour() const
{
	return this->nonlinbehaviour ;
}

NonLinearForm * ElementaryVolume::getNonLinearBehaviour() const
{
	return this->nonlinbehaviour ;
}

void ElementarySurface::setBehaviour(Form * f)
{
	this->behaviour = f ;
}

void ElementarySurface::setNonLinearBehaviour(NonLinearForm * f)
{
	this->nonlinbehaviour = f ;
}

void ElementarySurface::step(double dt, Vector *displacements)
{
	getState()->step(dt, displacements) ;
	getBehaviour()->updateElementState(dt, getState()) ;
}

void ElementarySurface::nonLinearStep(double dt, Vector *displacements)
{
	getState()->step(dt, displacements) ;
	
// 	if(getNonLinearBehaviour() !=NULL)
		getNonLinearBehaviour()->step(dt,getState()) ;
}

void ElementaryVolume::step(double dt, Vector *displacements)
{
	getState()->step(dt, displacements) ;
	getBehaviour()->updateElementState(dt, getState()) ;
}

void ElementaryVolume::nonLinearStep(double dt, Vector *displacements)
{
	
	this->state->step(dt, displacements) ;
	
	if(this->nonlinbehaviour !=NULL)
		this->nonlinbehaviour->step(dt, this->state) ;
}

const std::vector<std::pair<size_t, Function> > ElementarySurface::getDofs() const
{
	std::vector<std::pair<size_t, Function> > ret ;
	for (size_t i = 0 ; i < getShapeFunctions().size() ; i++)
	{
		ret.push_back(std::make_pair(getBoundingPoint(i).id, getShapeFunction(i))) ;
	}
	for (size_t i = 0 ; i < getEnrichmentFunctions().size() ; i++)
	{
		ret.push_back(std::make_pair(getEnrichmentFunction(i).second.getDofID(),getEnrichmentFunction(i).second)) ;
	}
	
	return ret ;
}

const std::vector< size_t > ElementarySurface::getDofIds() const
{
	std::vector<size_t> ret ;
	for (size_t i = 0 ; i < getShapeFunctions().size() ; i++)
	{
		if(getBoundingPoint(i).id >= 0)
			ret.push_back(getBoundingPoint(i).id) ;
		else
			std::cout << "negative ID, check numbering !" << std::endl ;
	}
	for (size_t i = 0 ; i < getEnrichmentFunctions().size() ; i++)
	{
		if(getEnrichmentFunction(i).second.getDofID() >= 0)
			ret.push_back(getEnrichmentFunction(i).second.getDofID()) ;
		else
			std::cout << "negative ID, check numbering !" << std::endl ;
	}
	
	return ret ;
}



std::valarray< std::pair<Point, double> > TriElement::genGaussPoints() const
{
	size_t ordre = 0;
	std::valarray< std::pair<Point, double> > fin ;
	switch(order)
	{
	case CONSTANT:
	case LINEAR:
		{
			ordre = 1 ;
			fin.resize(ordre);
			fin[0] = std::pair<Point, double>(Point(0.333333333333333, 0.333333333333333), 0.5) ;
			break ;
		}
	case QUADRATIC:
	case CUBIC:
		{
			ordre = 4 ;
			fin.resize(ordre);
			fin[0] = std::pair<Point, double>(Point(0.2, 0.2), 0.260416666666667) ;
			fin[1] = std::pair<Point, double>(Point(0.6, 0.2), 0.260416666666667) ;
			fin[2] = std::pair<Point, double>(Point(0.2, 0.6), 0.260416666666667) ;
			fin[3] = std::pair<Point, double>(Point(1./3., 1./3.), -0.28125) ;
			break ;
		}
	case QUADRIC:
	case QUINTIC:
		{
			ordre = 7 ;
			fin.resize(ordre);
			fin[0] = std::pair<Point, double>(Point(0.101286507323456, 0.101286507323456), 0.062969590272413) ;
			fin[1] = std::pair<Point, double>(Point(0.797426985353087, 0.101286507323456), 0.062969590272413) ;
			fin[2] = std::pair<Point, double>(Point(0.101286507323456, 0.797426985353087), 0.062969590272413) ;
			fin[3] = std::pair<Point, double>(Point(0.470142064105115, 0.059715871789770), 0.066197076394253) ;
			fin[4] = std::pair<Point, double>(Point(0.470142064105115, 0.470142064105115), 0.066197076394253) ;
			fin[5] = std::pair<Point, double>(Point(0.059715871789770, 0.470142064105115), 0.066197076394253) ;
			fin[6] = std::pair<Point, double>(Point(0.333333333333333, 0.333333333333333), 0.1125) ;
			break;
		}
	case CONSTANT_TIME_LINEAR:
		{
			ordre = 1 ;
			fin.resize(ordre);
			fin[0] = std::pair<Point, double>(Point(0.333333333333333, 0.333333333333333,0,.5), 1) ;
			break;
		}
	case CONSTANT_TIME_QUADRATIC:
		{
			ordre = 2 ;
			fin.resize(ordre);
			fin[0] = std::pair<Point, double>(Point(0.333333333333333, 0.333333333333333,0,-0.577350269189626), 0.5) ;
			fin[1] = std::pair<Point, double>(Point(0.333333333333333, 0.333333333333333,0,-0.577350269189626), 0.5) ;
			break;
		}
	case LINEAR_TIME_LINEAR:
		{
			ordre = 1 ;
			fin.resize(ordre);
			fin[0] = std::pair<Point, double>(Point(0.333333333333333, 0.333333333333333,0,0), 1) ;
			break ;
		}
	case LINEAR_TIME_QUADRATIC:
		{
			ordre = 2 ;
			fin.resize(ordre);
			fin[0] = std::pair<Point, double>(Point(0.333333333333333, 0.333333333333333,0,-0.577350269189626), 0.5) ;
			fin[1] = std::pair<Point, double>(Point(0.333333333333333, 0.333333333333333,0,0.577350269189626), 0.5) ;
			break ;
		}
	case QUADRATIC_TIME_LINEAR:
		{
			ordre = 4 ;
			fin.resize(ordre);
			fin[0] = std::pair<Point, double>(Point(0.2, 0.2,0,0), 0.260416666666667*2) ;
			fin[1] = std::pair<Point, double>(Point(0.6, 0.2,0,0), 0.260416666666667*2) ;
			fin[2] = std::pair<Point, double>(Point(0.2, 0.6,0,0), 0.260416666666667*2) ;
			fin[3] = std::pair<Point, double>(Point(1./3., 1./3.,0,0), -0.28125*2) ;
			break;
		}
	case QUADRATIC_TIME_QUADRATIC:
		{
			ordre = 8 ;
			fin.resize(ordre);
			fin[0] = std::pair<Point, double>(Point(0.2, 0.2,0,-0.577350269189626), 0.260416666666667) ;
			fin[1] = std::pair<Point, double>(Point(0.6, 0.2,0,-0.577350269189626), 0.260416666666667) ;
			fin[2] = std::pair<Point, double>(Point(0.2, 0.6,0,-0.577350269189626), 0.260416666666667) ;
			fin[3] = std::pair<Point, double>(Point(1./3., 1./3.,0,-0.577350269189626), -0.28125) ;
			fin[4] = std::pair<Point, double>(Point(0.2, 0.2,0,0.577350269189626), 0.260416666666667) ;
			fin[5] = std::pair<Point, double>(Point(0.6, 0.2,0,0.577350269189626), 0.260416666666667) ;
			fin[6] = std::pair<Point, double>(Point(0.2, 0.6,0,0.577350269189626), 0.260416666666667) ;
			fin[7] = std::pair<Point, double>(Point(1./3., 1./3.,0,0.577350269189626), -0.28125) ;
			break;
		}
	case CUBIC_TIME_LINEAR:
		{
			ordre = 4 ;
			fin.resize(ordre);
			fin[0] = std::pair<Point, double>(Point(0.2, 0.2), 0.260416666666667*2) ;
			fin[1] = std::pair<Point, double>(Point(0.6, 0.2), 0.260416666666667*2) ;
			fin[2] = std::pair<Point, double>(Point(0.2, 0.6), 0.260416666666667*2) ;
			fin[3] = std::pair<Point, double>(Point(1./3., 1./3.), -0.28125*2) ;
			break ;
		}
	case CUBIC_TIME_QUADRATIC:
		{
			ordre = 8 ;
			fin.resize(ordre);
			fin[0] = std::pair<Point, double>(Point(0.2, 0.2,0,-0.577350269189626), 0.260416666666667) ;
			fin[1] = std::pair<Point, double>(Point(0.6, 0.2,0,-0.577350269189626), 0.260416666666667) ;
			fin[2] = std::pair<Point, double>(Point(0.2, 0.6,0,-0.577350269189626), 0.260416666666667) ;
			fin[3] = std::pair<Point, double>(Point(1./3., 1./3.,0,-0.577350269189626), -0.28125) ;
			fin[4] = std::pair<Point, double>(Point(0.2, 0.2,0,0.577350269189626), 0.260416666666667) ;
			fin[5] = std::pair<Point, double>(Point(0.6, 0.2,0,0.577350269189626), 0.260416666666667) ;
			fin[6] = std::pair<Point, double>(Point(0.2, 0.6,0,0.577350269189626), 0.260416666666667) ;
			fin[7] = std::pair<Point, double>(Point(1./3., 1./3.,0,0.577350269189626), -0.28125) ;
			break ;
		}
	case QUADRIC_TIME_LINEAR:
		{
			ordre = 7 ;
			fin.resize(ordre);
			fin[0] = std::pair<Point, double>(Point(0.101286507323456, 0.101286507323456), 0.062969590272413*2) ;
			fin[1] = std::pair<Point, double>(Point(0.797426985353087, 0.101286507323456), 0.062969590272413*2) ;
			fin[2] = std::pair<Point, double>(Point(0.101286507323456, 0.797426985353087), 0.062969590272413*2) ;
			fin[3] = std::pair<Point, double>(Point(0.470142064105115, 0.059715871789770), 0.066197076394253*2) ;
			fin[4] = std::pair<Point, double>(Point(0.470142064105115, 0.470142064105115), 0.066197076394253*2) ;
			fin[5] = std::pair<Point, double>(Point(0.059715871789770, 0.470142064105115), 0.066197076394253*2) ;
			fin[6] = std::pair<Point, double>(Point(0.333333333333333, 0.333333333333333), 0.1125*2) ;
			break ;
		}
	case QUADRIC_TIME_QUADRATIC:
		{
			ordre = 14 ;
			fin.resize(ordre);
			fin[0] = std::pair<Point, double>(Point(0.101286507323456, 0.101286507323456,0,-0.577350269189626), 0.062969590272413) ;
			fin[1] = std::pair<Point, double>(Point(0.797426985353087, 0.101286507323456,0,-0.577350269189626), 0.062969590272413) ;
			fin[2] = std::pair<Point, double>(Point(0.101286507323456, 0.797426985353087,0,-0.577350269189626), 0.062969590272413) ;
			fin[3] = std::pair<Point, double>(Point(0.470142064105115, 0.059715871789770,0,-0.577350269189626), 0.066197076394253) ;
			fin[4] = std::pair<Point, double>(Point(0.470142064105115, 0.470142064105115,0,-0.577350269189626), 0.066197076394253) ;
			fin[5] = std::pair<Point, double>(Point(0.059715871789770, 0.470142064105115,0,-0.577350269189626), 0.066197076394253) ;
			fin[6] = std::pair<Point, double>(Point(0.333333333333333, 0.333333333333333,0,-0.577350269189626), 0.1125) ;
			fin[7] = std::pair<Point, double>(Point(0.101286507323456, 0.101286507323456,0,0.577350269189626), 0.062969590272413) ;
			fin[8] = std::pair<Point, double>(Point(0.797426985353087, 0.101286507323456,0,0.577350269189626), 0.062969590272413) ;
			fin[9] = std::pair<Point, double>(Point(0.101286507323456, 0.797426985353087,0,0.577350269189626), 0.062969590272413) ;
			fin[10] = std::pair<Point, double>(Point(0.470142064105115, 0.059715871789770,0,0.577350269189626), 0.066197076394253) ;
			fin[11] = std::pair<Point, double>(Point(0.470142064105115, 0.470142064105115,0,0.577350269189626), 0.066197076394253) ;
			fin[12] = std::pair<Point, double>(Point(0.059715871789770, 0.470142064105115,0,0.577350269189626), 0.066197076394253) ;
			fin[13] = std::pair<Point, double>(Point(0.333333333333333, 0.333333333333333,0,0.577350269189626), 0.1125) ;
			break ;
		}
	case QUINTIC_TIME_LINEAR:
		{
			ordre = 7 ;
			fin.resize(ordre);
			fin[0] = std::pair<Point, double>(Point(0.101286507323456, 0.101286507323456), 0.062969590272413*2) ;
			fin[1] = std::pair<Point, double>(Point(0.797426985353087, 0.101286507323456), 0.062969590272413*2) ;
			fin[2] = std::pair<Point, double>(Point(0.101286507323456, 0.797426985353087), 0.062969590272413*2) ;
			fin[3] = std::pair<Point, double>(Point(0.470142064105115, 0.059715871789770), 0.066197076394253*2) ;
			fin[4] = std::pair<Point, double>(Point(0.470142064105115, 0.470142064105115), 0.066197076394253*2) ;
			fin[5] = std::pair<Point, double>(Point(0.059715871789770, 0.470142064105115), 0.066197076394253*2) ;
			fin[6] = std::pair<Point, double>(Point(0.333333333333333, 0.333333333333333), 0.1125*2) ;
			break;
		}
	case QUINTIC_TIME_QUADRATIC:
		{
			ordre = 14 ;
			fin.resize(ordre);
			fin[0] = std::pair<Point, double>(Point(0.101286507323456, 0.101286507323456,0,-0.577350269189626), 0.062969590272413) ;
			fin[1] = std::pair<Point, double>(Point(0.797426985353087, 0.101286507323456,0,-0.577350269189626), 0.062969590272413) ;
			fin[2] = std::pair<Point, double>(Point(0.101286507323456, 0.797426985353087,0,-0.577350269189626), 0.062969590272413) ;
			fin[3] = std::pair<Point, double>(Point(0.470142064105115, 0.059715871789770,0,-0.577350269189626), 0.066197076394253) ;
			fin[4] = std::pair<Point, double>(Point(0.470142064105115, 0.470142064105115,0,-0.577350269189626), 0.066197076394253) ;
			fin[5] = std::pair<Point, double>(Point(0.059715871789770, 0.470142064105115,0,-0.577350269189626), 0.066197076394253) ;
			fin[6] = std::pair<Point, double>(Point(0.333333333333333, 0.333333333333333,0,-0.577350269189626), 0.1125) ;
			fin[7] = std::pair<Point, double>(Point(0.101286507323456, 0.101286507323456,0,0.577350269189626), 0.062969590272413) ;
			fin[8] = std::pair<Point, double>(Point(0.797426985353087, 0.101286507323456,0,0.577350269189626), 0.062969590272413) ;
			fin[9] = std::pair<Point, double>(Point(0.101286507323456, 0.797426985353087,0,0.577350269189626), 0.062969590272413) ;
			fin[10] = std::pair<Point, double>(Point(0.470142064105115, 0.059715871789770,0,0.577350269189626), 0.066197076394253) ;
			fin[11] = std::pair<Point, double>(Point(0.470142064105115, 0.470142064105115,0,0.577350269189626), 0.066197076394253) ;
			fin[12] = std::pair<Point, double>(Point(0.059715871789770, 0.470142064105115,0,0.577350269189626), 0.066197076394253) ;
			fin[13] = std::pair<Point, double>(Point(0.333333333333333, 0.333333333333333,0,0.577350269189626), 0.1125) ;
			break;
		}
	}

// 	if(moved)
// 	{
		for(size_t i = 0 ; i < fin.size() ; i++)
		{
			fin[i].second*=jacobianAtPoint(fin[i].first);
		}
// 	}
// 	else
// 	{
// 		double j = jacobianAtPoint(Point(1./3., 1./3.)) ;
// 		for(size_t i = 0 ; i < fin.size() ; i++)
// 		{
// 			fin[i].second*=j;
// 		}
// 	}
	
	return fin ;
} ;

void TriElement::computeCenter()
{
	this->Triangle::computeCenter() ;
}


TriElement::TriElement( Point * p0,  Point * p1,  Point * p2, bool father ) : Triangle(p0, p1, p2), ElementarySurface(father), moved(false) { };
	
TriElement::TriElement(Order order , bool father ): ElementarySurface(father),moved(false) 
{
	this->order = order ;
	
	switch(order)
	{
	case CONSTANT :
		{
			shapefunc = new std::valarray<Function>(1) ;
			Matrix m(1, 1) ;
			m[0][0] = 1 ;
			(*shapefunc)[0] = Function(m) ;
			break ;
		}
	case LINEAR :
		{
			shapefunc = new std::valarray<Function>(3) ;
			Matrix xi(2,2) ; xi[1][0] = 1 ;
			Matrix eta(2,2) ; eta[0][1] = 1 ;
			Matrix one(2,2) ; one[0][0] = 1 ;
		//0
			(*shapefunc)[0] = Function(eta) ;
		//1
			(*shapefunc)[1] = Function(one-xi-eta) ;
		//2
			(*shapefunc)[2] = Function(xi) ;
			break ;
		}
	case QUADRATIC :
		{
			shapefunc = new std::valarray<Function>(6) ;
			
			
			Matrix xi(3,3) ; xi[1][0] = 1 ;
			Matrix eta(3,3) ; eta[0][1] = 1 ;
			Matrix one(3,3) ; one[0][0] = 1 ;
			Matrix xi_eta(3,3) ; xi_eta[1][1] = 1 ;
			Matrix xi_xi(3,3) ; xi_xi[2][0] = 1 ;
			Matrix eta_eta(3,3) ; eta_eta[0][2] = 1 ;
		//0
			(*shapefunc)[0] = Function(eta_eta*2 - eta) ;
		//1
			(*shapefunc)[1] = Function(eta*4 - xi_eta*4 - eta_eta*4) ;
		//2
			(*shapefunc)[2] = Function(one - xi*3 - eta*3 + xi_eta*4 + xi_xi*2 + eta_eta*2) ;
		//3
			(*shapefunc)[3] = Function(xi*4 - xi_xi*4 - xi_eta*4) ;
		//4
			(*shapefunc)[4] = Function(xi_xi*2 - xi) ;
		//5
			(*shapefunc)[5] = Function(xi_eta*4) ;
			break ;
		}
	case CONSTANT_TIME_LINEAR :
		{
			shapefunc = new std::valarray<Function>(2) ;
			Matrix m(1, 2) ;
			std::valarray<Matrix> v(m,1) ; v[0] = m ;
			std::valarray<std::valarray<Matrix> > vv(v,1) ;
			vv[0] = v ;
			
			vv[0][0][0][1] = 1 ;
			(*shapefunc)[0] = Function(m) ;
			vv[0][0][0][0] = 1 ;
			vv[0][0][0][1] = -1 ;
			(*shapefunc)[1] = Function(m) ;
			break ;
		}
	case CONSTANT_TIME_QUADRATIC :
		{
			shapefunc = new std::valarray<Function>(3) ;
			Matrix m(1, 3) ;
			std::valarray<Matrix> v(m,1) ; v[0] = m ;
			std::valarray<std::valarray<Matrix> > vv(v,1) ;
			vv[0] = v ;
			
			vv[0][0][0][2] = .5 ;
			vv[0][0][0][1] = -.5 ;
			(*shapefunc)[0] = Function(m) ;
			vv[0][0][0][0] = 1 ;
			vv[0][0][0][1] = 0 ;
			vv[0][0][0][2] = -1 ;
			(*shapefunc)[1] = Function(m) ;
			vv[0][0][0][0] = 0.5 ;
			vv[0][0][0][1] = 0 ;
			vv[0][0][0][2] = 0.5 ;
			(*shapefunc)[1] = Function(m) ;
			break ;
		}
	case LINEAR_TIME_LINEAR :
		{
			shapefunc = new std::valarray<Function>(6) ;
			Matrix m(2, 2) ;
			std::valarray<Matrix> v(m, 2) ; v[0] = m ; v[1] = m ;
			std::valarray<std::valarray<Matrix> >vv(v, 2) ; vv[0] = v ;  vv[1] = v ;

		//0
			vv[0][1][0][0] = .5  ;
			vv[0][1][0][1] = -.5  ;
			(*shapefunc)[0] = Function(vv) ;
			vv[0][1][0][0] = 0  ;
			vv[0][1][0][1] = 0  ;
		//1
			vv[0][0][0][0] = .5  ;
			vv[1][0][0][0] = -.5  ;
			vv[0][1][0][0] = -.5  ;
			vv[0][0][0][1] = -.5  ;
			vv[1][0][0][1] = .5  ;
			vv[0][1][0][1] = .5  ;
			(*shapefunc)[1] = Function(vv) ;
			vv[0][0][0][0] = 0  ;
			vv[1][0][0][0] = 0  ;
			vv[0][1][0][0] = 0  ;
			vv[0][0][0][1] = 0  ;
			vv[1][0][0][1] = 0  ;
			vv[0][1][0][1] = 0  ;
		//2
			vv[1][0][0][0] = .5  ;
			vv[1][0][0][1] = -.5  ;
			(*shapefunc)[2] = Function(vv) ;
			vv[1][0][0][0] = 0  ;
			vv[1][0][0][1] = 0  ;
		//3
			vv[0][1][0][0] = .5  ;
			vv[0][1][0][1] = .5  ;
			(*shapefunc)[3] = Function(vv) ;
			vv[0][1][0][0] = 0  ;
			vv[0][1][0][1] = 0  ;
		//4
			vv[0][0][0][0] = .5  ;
			vv[1][0][0][0] = -.5  ;
			vv[0][1][0][0] = -.5  ;
			vv[0][0][0][1] = .5  ;
			vv[1][0][0][1] = -.5  ;
			vv[0][1][0][1] = -.5  ;
			(*shapefunc)[4] = Function(vv) ;
			vv[0][0][0][0] = 0  ;
			vv[1][0][0][0] = 0  ;
			vv[0][1][0][0] = 0  ;
			vv[0][0][0][1] = 0  ;
			vv[1][0][0][1] = 0  ;
			vv[0][1][0][1] = 0  ;
		//5
			vv[1][0][0][0] = .5  ;
			vv[1][0][0][1] = .5  ;
			(*shapefunc)[5] = Function(vv) ;
			break ;
		}
	case LINEAR_TIME_QUADRATIC :
		{
			shapefunc = new std::valarray<Function>(9) ;
			
		//0
			(*shapefunc)[0] = Function("y 0.5 t t * * 0.5 t * - *") ;
		//1
			(*shapefunc)[1] = Function("1 x y - - 0.5 t t * * 0.5 t * - *") ;
		//2
			(*shapefunc)[2] = Function("x 0.5 t t * * 0.5 t * - *") ;
		//3
			(*shapefunc)[3] = Function("y 1 t t * - *") ;
		//4
			(*shapefunc)[4] = Function("1 x y - - 1 t t * - *") ;
		//5
			(*shapefunc)[5] = Function("x 1 t t * - *") ;
		//6
			(*shapefunc)[6] = Function("y 0.5 t t * * 0.5 t * + *") ;
		//7
			(*shapefunc)[7] = Function("1 x y - - 0.5 t t * * 0.5 t * + *") ;
		//8
			(*shapefunc)[8] = Function("x 0.5 t t * * 0.5 t * + *") ;
			break ;
		}
	case QUADRATIC_TIME_LINEAR:
		{
			shapefunc = new std::valarray<Function>(12) ;
			
		//0
			(*shapefunc)[0] = Function("y y 2 * * y - 0.5 0.5 t * - *") ;
		//1
			(*shapefunc)[1] = Function("y 4 * x y 4 * * - y y 4 * * - 0.5 0.5 t * - *") ;
		//2
			(*shapefunc)[2] = Function("1 x 3 * - y 3 * - x y 4 * * + x x 2 * * + y y 2 * * + 0.5 0.5 t * - *") ;
		//3
			(*shapefunc)[3] = Function("4 x *  x x 4 * * - x y 4 * * - 0.5 0.5 t * - *") ;
		//4
			(*shapefunc)[4] = Function("x x 2 * * x - 0.5 0.5 t * - *") ;
		//5
			(*shapefunc)[5] = Function("x y 4 * * 0.5 0.5 t * - *") ;
		//0
			(*shapefunc)[6] = Function("y y 2 * * y - 0.5 0.5 t * + *") ;
		//1
			(*shapefunc)[7] = Function("y 4 * x y 4 * * - y y 4 * * - 0.5 0.5 t * + *") ;
		//2
			(*shapefunc)[8] = Function("1 x 3 * - y 3 * - x y 4 * * + x x 2 * * + y y 2 * * + 0.5 0.5 t * + *") ;
		//3
			(*shapefunc)[9] = Function("4 x *  x x 4 * * - x y 4 * * - 0.5 0.5 t * + *") ;
		//4
			(*shapefunc)[10] = Function("x x 2 * * x - 0.5 0.5 t * + *") ;
		//5
			(*shapefunc)[11] = Function("x y 4 * * 0.5 0.5 t * + *") ;
			break ;
		}
	case QUADRATIC_TIME_QUADRATIC:
		{
			shapefunc = new std::valarray<Function>(18) ;
			
			(*shapefunc)[0] = Function("y y 2 * * y - 0.5 t t * * 0.5 t * - *") ;
			(*shapefunc)[1] = Function("y 4 * x y 4 * * - y y 4 * * - 0.5 t t * * 0.5 t * - *") ;
			(*shapefunc)[2] = Function("1 x 3 * - y 3 * - x y 4 * * + x x 2 * * + y y 2 * * + 0.5 t t * * 0.5 t * - *") ;
			(*shapefunc)[3] = Function("4 x *  x x 4 * * - x y 4 * * - 0.5 t t * * 0.5 t * - *") ;
			(*shapefunc)[4] = Function("x x 2 * * x - 0.5 t t * * 0.5 t * - *") ;
			(*shapefunc)[5] = Function("x y 4 * * 0.5 t t * * 0.5 t * - *") ;

			(*shapefunc)[6] = Function("y y 2 * * y - 1 t t * - *") ;
			(*shapefunc)[7] = Function("y 4 * x y 4 * * - y y 4 * * - 1 t t * - *") ;
			(*shapefunc)[8] = Function("1 x 3 * - y 3 * - x y 4 * * + x x 2 * * + y y 2 * * + 1 t t * - *") ;
			(*shapefunc)[9] = Function("4 x *  x x 4 * * - x y 4 * * - 1 t t * - *") ;
			(*shapefunc)[10] = Function("x x 2 * * x - 1 t t * - *") ;
			(*shapefunc)[11] = Function("x y 4 * * 1 t t * - *") ;

			(*shapefunc)[12] = Function("y y 2 * * y - 0.5 t t * * 0.5 t * + *") ;
			(*shapefunc)[13] = Function("y 4 * x y 4 * * - y y 4 * * - 0.5 t t * * 0.5 t * + *") ;
			(*shapefunc)[14] = Function("1 x 3 * - y 3 * - x y 4 * * + x x 2 * * + y y 2 * * + 0.5 t t * * 0.5 t * + *") ;
			(*shapefunc)[15] = Function("4 x *  x x 4 * * - x y 4 * * - 0.5 t t * * 0.5 t * + *") ;
			(*shapefunc)[16] = Function("x x 2 * * x - 0.5 t t * * 0.5 t * + *") ;
			(*shapefunc)[17] = Function("x y 4 * * 0.5 t t * * 0.5 t * + *") ;
			break ;
		}
	default:
		{
			assert(false) ;
			break ;
		}
	}
	
	
	for(size_t i = 0 ; i < this->Triangle::getBoundingPoints().size() ; i++)
		this->Triangle::getBoundingPoint(i).id = -1 ;
	
}
	
void TriElement::refresh(const TriElement * parent)
{
	
	this->order =  parent->getOrder() ;
	this->shapefunc = parent->shapefunc ;
	
}
	
std::vector<std::vector<Matrix> > TriElement::getElementaryMatrix() const 
{
	return std::vector<std::vector<Matrix> >() ;
}

std::vector<std::vector<Matrix> > TriElement::getNonLinearElementaryMatrix(Vector * state) const 
{
	return std::vector<std::vector<Matrix> >() ;
}
	
Function TriElement::jacobian() const 
{
	
	Function xdxi = this->getdXTransform(XI) ;
	Function ydxi = this->getdYTransform(XI) ;
	Function xdeta = this->getdXTransform(ETA) ;
	Function ydeta = this->getdYTransform(ETA) ;
	
	Function ret = ydeta*xdxi - ydxi*xdeta  ;
	
	return ret ;
	
}
	
double  TriElement::jacobianAtPoint(const Point p) const 
{
	if(order < CONSTANT_TIME_LINEAR)
	{
		double xdxi = this->getdXTransform(XI, p) ;
		double ydxi = this->getdYTransform(XI, p) ;
		double xdeta = this->getdXTransform(ETA, p) ;
		double ydeta = this->getdYTransform(ETA, p) ;
		
		return ydeta*xdxi - ydxi*xdeta ;
	}
	else
	{
		double xdxi = this->getdXTransform(XI, p) ;
		double ydxi = this->getdYTransform(XI, p) ;
		double zdxi = this->getdTTransform(XI, p) ;
		double xdeta = this->getdXTransform(ETA, p) ;
		double ydeta = this->getdYTransform(ETA, p) ;
		double zdeta = this->getdTTransform(ETA, p) ;
		double xdzeta = this->getdXTransform(TIME_VARIABLE,p) ;
		double ydzeta = this->getdYTransform(TIME_VARIABLE,p) ;
		double zdzeta = this->getdTTransform(TIME_VARIABLE,p) ;
		
		return xdxi*ydeta*zdzeta + zdeta*xdzeta*ydxi + ydzeta*zdxi*xdeta  -
			xdxi*ydzeta*zdeta - xdeta*ydxi*zdzeta - xdzeta*ydeta*zdxi ;
	}
	
}
	
Matrix TriElement::getInverseJacobianMatrix(const Point & p) const
{
	if(order < CONSTANT_TIME_LINEAR)
	{
		Matrix  J0(2,2) ;
		
		double xdxi = this->getdXTransform(XI,p) ;
		double ydxi = this->getdYTransform(XI,p) ;
		double xdeta = this->getdXTransform(ETA,p) ;
		double ydeta = this->getdYTransform(ETA,p) ;
		
		J0[0][0] = xdxi ; J0[0][1] = ydxi ; 
		J0[1][0] = xdeta ; J0[1][1] = ydeta ;
		invert2x2Matrix(J0) ;
		return J0 ;
	}
	else
	{
		Matrix  J0(3,3) ;
		double xdxi = this->getdXTransform(XI, p) ;
		double ydxi = this->getdYTransform(XI, p) ;
		double zdxi = this->getdTTransform(XI, p) ;
		double xdeta = this->getdXTransform(ETA, p) ;
		double ydeta = this->getdYTransform(ETA, p) ;
		double zdeta = this->getdTTransform(ETA, p) ;
		double xdzeta = this->getdXTransform(TIME_VARIABLE,p) ;
		double ydzeta = this->getdYTransform(TIME_VARIABLE,p) ;
		double zdzeta = this->getdTTransform(TIME_VARIABLE,p) ;
		
		J0[0][0] = xdxi ; J0[0][1] = ydxi ; J0[0][2] = zdxi ; 
		J0[1][0] = xdeta ; J0[1][1] = ydeta ; J0[1][2] = zdeta ;
		J0[2][0] = xdzeta ; J0[2][1] = ydzeta ; J0[2][2] = zdzeta ;
		invert3x3Matrix(J0) ;
		return J0 ;
	}
}
	
std::valarray< std::pair<Point, double> > TriElement::getGaussPoints() const
{
	return genGaussPoints() ;
}

bool TriElement::isMoved() const
{
	return moved ;
}

void TriElement::print() const
{
	std::cout << "coucou !" << std::endl ;
}

Point TriElement::inLocalCoordinates(const Point &p) const
{
		// in barycentric coordinates, we have the following :
		
	size_t factor = 1 ;
	
	if(order == QUADRATIC || order == QUADRATIC_TIME_LINEAR || order == QUADRATIC_TIME_QUADRATIC)
		factor = 2 ;
	
	Matrix S(3,3) ;
	S[0][0] = this->getBoundingPoint(0).x ; S[0][1] = this->getBoundingPoint(factor).x ;  S[0][2] = this->getBoundingPoint(factor*2).x ; 
	S[1][0] = this->getBoundingPoint(0).y ; S[1][1] = this->getBoundingPoint(factor).y ;  S[1][2] = this->getBoundingPoint(factor*2).y ; 
	S[2][0] = 1 ; S[2][1] = 1 ;  S[2][2] = 1 ; 
	
	Vector v(3) ; 
	v[0] = p.x ;
	v[1] = p.y ;
	v[2] = 1 ;
	
	Vector coeff = inverse3x3Matrix( S) * v ;
	
	return Point(0,1,0,p.t)*coeff[0] + Point(0,0,0,p.t)*coeff[1] + Point(1,0,0,p.t)*coeff[2] ; 
}


std::vector<std::vector<Matrix> > TriElement::getNonLinearElementaryMatrix() const 
{
	return std::vector<std::vector<Matrix> >() ;
}

Vector TriElement::getNonLinearForces() const 
{
	
	return Vector(0) ;
}


std::valarray< std::pair<Point, double> > TetrahedralElement::genGaussPoints() const
{
	size_t ordre=0 ;
	if(order == LINEAR || order == LINEAR_TIME_LINEAR)
		ordre = 1 ;
	if(order == LINEAR_TIME_QUADRATIC)
		ordre = 2 ;
	else if (order == CUBIC || order == QUADRATIC || order == CUBIC_TIME_LINEAR || order == QUADRATIC_TIME_LINEAR)
		ordre = 5 ;
	else
		ordre = 10 ;
	
	std::valarray< std::pair<Point, double> > fin(ordre);
	
	if(order == LINEAR )
	{
		fin[0] = std::pair<Point, double>(Point(0.25, 0.25, 0.25), 0.1666666666666667) ;
	}
	else if (order == CUBIC || order == QUADRATIC )
	{
		fin[0] = std::pair<Point, double>(Point(0.25, 0.25, 0.25), -0.133333333333333) ;
		fin[1] = std::pair<Point, double>(Point(0.166666666666667, 0.166666666666667, 0.166666666666667), 0.075) ;
		fin[2] = std::pair<Point, double>(Point(0.5, 0.166666666666667, 0.166666666666667), 0.075) ;
		fin[3] = std::pair<Point, double>(Point(0.166666666666667, 0.5, 0.166666666666667), 0.075) ;
		fin[4] = std::pair<Point, double>(Point(0.166666666666667, 0.166666666666667, 0.5), 0.075) ;
	}
	else if(order == LINEAR_TIME_LINEAR )
	{
		fin[0] = std::pair<Point, double>(Point(0.25, 0.25, 0.25), 0.1666666666666667*2.) ;
	}
	else if(order == LINEAR_TIME_QUADRATIC )
	{
		fin[0] = std::pair<Point, double>(Point(0.25, 0.25, 0.25,-0.577350269189626), 0.1666666666666667) ;
		fin[0] = std::pair<Point, double>(Point(0.25, 0.25, 0.25,0.577350269189626), 0.1666666666666667) ;
	}
	else if (order == CUBIC_TIME_LINEAR || order == QUADRATIC_TIME_LINEAR )
	{
		fin[0] = std::pair<Point, double>(Point(0.25, 0.25, 0.25), -0.133333333333333*2.) ;
		fin[1] = std::pair<Point, double>(Point(0.166666666666667, 0.166666666666667, 0.166666666666667), 0.075*2.) ;
		fin[2] = std::pair<Point, double>(Point(0.5, 0.166666666666667, 0.166666666666667), 0.075*2.) ;
		fin[3] = std::pair<Point, double>(Point(0.166666666666667, 0.5, 0.166666666666667), 0.075*2.) ;
		fin[4] = std::pair<Point, double>(Point(0.166666666666667, 0.166666666666667, 0.5), 0.075*2.) ;
	}
	else if (order == CUBIC_TIME_QUADRATIC || order == QUADRATIC_TIME_QUADRATIC )
	{
		fin[0] = std::pair<Point, double>(Point(0.25, 0.25, 0.25,-0.577350269189626), -0.133333333333333) ;
		fin[1] = std::pair<Point, double>(Point(0.166666666666667, 0.166666666666667, 0.166666666666667,-0.577350269189626), 0.075) ;
		fin[2] = std::pair<Point, double>(Point(0.5, 0.166666666666667, 0.166666666666667,-0.577350269189626), 0.075) ;
		fin[3] = std::pair<Point, double>(Point(0.166666666666667, 0.5, 0.166666666666667,-0.577350269189626), 0.075) ;
		fin[4] = std::pair<Point, double>(Point(0.166666666666667, 0.166666666666667, 0.5,-0.577350269189626), 0.075) ;
		fin[5] = std::pair<Point, double>(Point(0.25, 0.25, 0.25, 0.577350269189626), -0.133333333333333) ;
		fin[6] = std::pair<Point, double>(Point(0.166666666666667, 0.166666666666667, 0.166666666666667, 0.577350269189626), 0.075) ;
		fin[7] = std::pair<Point, double>(Point(0.5, 0.166666666666667, 0.166666666666667, 0.577350269189626), 0.075) ;
		fin[8] = std::pair<Point, double>(Point(0.166666666666667, 0.5, 0.166666666666667, 0.577350269189626), 0.075) ;
		fin[9] = std::pair<Point, double>(Point(0.166666666666667, 0.166666666666667, 0.5, 0.577350269189626), 0.075) ;
	}
	else
	{
		std::cout << "this set of Gauss points is not implemented" << std::endl ;
		assert(false) ;
	}
	
// 	if( moved)
// 	{
		for(size_t i = 0 ; i < fin.size() ; i++)
		{
			fin[i].second*=this->jacobianAtPoint(fin[i].first) ;
		}
// 	}
// 	else
// 	{
// 		double j = this->jacobianAtPoint(Point(.25, .25, .25)) ;
// 		for(size_t i = 0 ; i < fin.size() ; i++)
// 		{
// 			fin[i].second*=j ;
// 		}
// 	}
	
	return fin ;
}

void TetrahedralElement::computeCenter()
{
	this->Tetrahedron::computeCenter() ;
}

TetrahedralElement::TetrahedralElement( Point * p0,  Point * p1,  Point * p2, Point * p3, bool father ) : Tetrahedron(p0, p1, p2, p3), ElementaryVolume(father), moved(false) 
{	}

TetrahedralElement::TetrahedralElement(Order order , bool father): ElementaryVolume(father),moved(false)
{
	this->order = order ;
	
	if(order == LINEAR)
	{
		this->Tetrahedron::sampleSurface(4) ;
		shapefunc = new std::valarray<Function>(4) ;
		Matrix zero ;
		std::valarray<Matrix> xi(zero, 2) ;
		xi[1][0][0] = 1 ;
		std::valarray<Matrix> eta(zero,2) ;
		eta[0][1][0] = 1 ;
		std::valarray<Matrix> zeta(zero,2) ;
		zeta[0][0][1] = 1 ;
		std::valarray<Matrix> f(zero,2) ;
		f[0][0][0] = 1 ;
		f[1][0][0] = -1 ;
		f[0][1][0] = -1 ;
		f[0][0][1] = -1 ;
			
			//0
		(*shapefunc)[0] = Function(xi) ;
			//1
		(*shapefunc)[1] = Function(eta) ;
			//2
		(*shapefunc)[2] = Function(zeta) ;
			//3
		(*shapefunc)[3] = Function(f) ;
	}
	else if(order == QUADRATIC)
	{
		
		shapefunc = new std::valarray<Function>(10) ;
		
		Matrix zero(3,3) ;
		std::valarray<Matrix> f0(zero, 3) ;
		f0[0][0][1] = -1 ; // z
		f0[0][0][2] =  2 ; // z^2
		
		std::valarray<Matrix> f1(zero,3) ;
		f1[0][0][1] = 4 ; // z
		f1[0][0][2] = -4 ; // z^2
		f1[1][0][1] = -4 ; // zx
		f1[0][1][1] = -4 ; // zy
		
		std::valarray<Matrix> f2(zero,3) ;
		f2[0][0][0] = 1 ; // 1
		f2[0][0][1] = -3 ; // z
		f2[0][0][2] = 2 ; // z^2
		f2[1][0][1] = 4 ; // zx
		f2[0][1][1] = 4 ; // zy
		f2[1][0][0] = -3 ; // x
		f2[2][0][0] = 2 ; // x^2
		f2[1][1][0] = 4 ; // xy
		f2[0][1][0] = -3 ; // y
		f2[0][2][0] = 2 ; // y^2
		
		std::valarray<Matrix> f3(zero,3) ;
		f3[1][0][0] = 4 ; // x
		f3[2][0][0] = -4 ; // x^2
		f3[1][0][1] = -4 ; // zx
		f3[1][1][0] = -4 ; // xy
		
		std::valarray<Matrix> f4(zero, 3) ;
		f4[1][0][0] = -1 ; // x
		f4[2][0][0] = 2 ; // x^2
		
		std::valarray<Matrix> f5(zero,3) ;
		f5[1][1][0] = 4 ; // xy
		
		std::valarray<Matrix> f6(zero,3) ;
		f6[0][1][0] = -1 ; // y
		f6[0][2][0] = 2 ; // y^2
		
		std::valarray<Matrix> f7(zero,3) ;
		f7[0][1][1] = 4 ; // yz
		
		std::valarray<Matrix> f8(zero,3) ;
		f8[0][1][0] = 4 ; // y
		f8[0][2][0] = -4 ; // y^2
		f8[0][1][1] = -4 ; // zy
		f8[1][1][0] = -4 ; // xy
		
		
		std::valarray<Matrix> f9(zero,3) ;
		f9[1][0][1] = 4 ; // xz
			
		(*shapefunc)[0] = Function(f0) ;//z- z*2*(one-x-y-z) - x*z*2 - y*z*2 ; 
		(*shapefunc)[1] = Function(f1) ; //z*4*(one-x-y-z) ;
		(*shapefunc)[2] = Function(f2) ; //one-x-y-z-(one-x-y-z)*(x+y+z)*2 ;
		(*shapefunc)[3] = Function(f3) ; //x*4*(one-x-y-z) ;
		(*shapefunc)[4] = Function(f4) ; //x- x*2*(one-x-y-z) - x*z*2 - y*x*2 ; 
		(*shapefunc)[5] = Function(f5) ; //x*y*4 ; 
		(*shapefunc)[6] = Function(f6) ; //y- y*2*(one-x-y-z) - y*z*2 - y*x*2 ; 
		(*shapefunc)[7] = Function(f7) ; //y*z*4 ;
		(*shapefunc)[8] = Function(f8) ; //y*4*(one-x-y-z) ;
		(*shapefunc)[9] = Function(f9) ; //x*z*4 ;
	}
	else if(order == LINEAR_TIME_LINEAR)
	{
		shapefunc = new std::valarray<Function>(8) ;
			//0
		(*shapefunc)[0] = Function("z 0.5 0.5 t * - *") ;
			//1
		(*shapefunc)[1] = Function("1 x - y - z - 0.5 0.5 t * - *") ;
			//2
		(*shapefunc)[2] = Function("x 0.5 0.5 t * - *") ;
			//3
		(*shapefunc)[3] = Function("y 0.5 0.5 t * - *") ;
			//4
		(*shapefunc)[4] = Function("z 0.5 0.5 t * + *") ;
			//5
		(*shapefunc)[5] = Function("1 x - y - z - 0.5 0.5 t * + *") ;
			//6
		(*shapefunc)[6] = Function("x 0.5 0.5 t * + *") ;
			//7
		(*shapefunc)[7] = Function("y 0.5 0.5 t * + *") ;
	}
	else if(order == LINEAR_TIME_QUADRATIC)
	{
		shapefunc = new std::valarray<Function>(12) ;
			//0
		(*shapefunc)[0] = Function("x t 1 - t * 0.5 * *") ;
			//1
		(*shapefunc)[1] = Function("y t 1 - t * 0.5 * *") ;
			//2
		(*shapefunc)[2] = Function("z t 1 - t * 0.5 * *") ;
			//3
		(*shapefunc)[3] = Function("1 x - y - z - t 1 - t * 0.5 * *") ;
			//4
		(*shapefunc)[4] = Function("x 1 t t * - *") ;
			//5
		(*shapefunc)[5] = Function("y 1 t t * - *") ;
			//6
		(*shapefunc)[6] = Function("z 1 t t * - *") ;
			//7
		(*shapefunc)[7] = Function("1 x - y - z - 1 t t * - *") ;
		    //8
		(*shapefunc)[8] = Function("x t 1 + t * 0.5 * *") ;
			//9
		(*shapefunc)[9] = Function("y t 1 + t * 0.5 * *") ;
			//10
		(*shapefunc)[10] = Function("z t 1 + t * 0.5 * *") ;
			//11
		(*shapefunc)[11] = Function("1 x - y - z - t 1 + t * 0.5 * *") ;
	}
	else if(order == QUADRATIC_TIME_LINEAR)
	{
		shapefunc = new std::valarray<Function>(20) ;

			//0
		(*shapefunc)[0] = Function("2 z * 1 - z * 1 t - 0.5 * *") ;
			//1
		(*shapefunc)[1] = Function("1 x - y - z - z * 4 * 1 t - 0.5 * *") ;
			//2
		(*shapefunc)[2] = Function("1 z y x + + 3 * - z z* x x * y y * + + 2 * + z x *  z y * x y * + + 4 * + 1 t - 0.5 * *") ;
			//3
		(*shapefunc)[3] = Function("1 x - y - z - x * 4 * 1 t - 0.5 * *") ;
			//4
		(*shapefunc)[4] = Function("2 x * 1 - x * 1 t - 0.5 * *") ;
			//5
		(*shapefunc)[5] = Function("x y 4 * * 1 t - 0.5 * *") ;
			//6
		(*shapefunc)[6] = Function("2 y * 1 - y * 1 t - 0.5 * *") ;
			//7
		(*shapefunc)[7] = Function("y z 4 * * 1 t - 0.5 * *") ;
		    //8
		(*shapefunc)[8] = Function("1 x - y - z - y * 4 * 1 t - 0.5 * *") ;
			//9
		(*shapefunc)[9] = Function("x z 4 * * 1 t - 0.5 * *") ;
		    //10
		(*shapefunc)[10] = Function("2 z * 1 - z * 1 t + 0.5 * *") ;
			//11
		(*shapefunc)[11] = Function("1 x - y - z - z * 4 * 1 t + 0.5 * *") ;
			//12
		(*shapefunc)[12] = Function("1 z y x + + 3 * - z z* x x * y y * + + 2 * + z x *  z y * x y * + + 4 * + 1 t + 0.5 * *") ;
			//13
		(*shapefunc)[13] = Function("1 x - y - z - x * 4 * 1 t + 0.5 * *") ;
			//14
		(*shapefunc)[14] = Function("2 x * 1 - x * 1 t + 0.5 * *") ;
			//15
		(*shapefunc)[15] = Function("x y 4 * * 1 t + 0.5 * *") ;
			//16
		(*shapefunc)[16] = Function("2 y * 1 - y * 1 t + 0.5 * *") ;
			//17
		(*shapefunc)[17] = Function("y z 4 * * 1 t + 0.5 * *") ;
		    //18
		(*shapefunc)[18] = Function("1 x - y - z - y * 4 * 1 t + 0.5 * *") ;
			//19
		(*shapefunc)[19] = Function("x z 4 * * 1 t + 0.5 * *") ;
	}
	else if(order == QUADRATIC_TIME_QUADRATIC)
	{
		shapefunc = new std::valarray<Function>(30) ;
		
			//0
		(*shapefunc)[0] = Function("2 z * 1 - z * t 1 - t * 0.5 * *") ;
			//1
		(*shapefunc)[1] = Function("1 x - y - z - z * 4 * t 1 - t * 0.5 * *") ;
			//2
		(*shapefunc)[2] = Function("1 z y x + + 3 * - z z* x x * y y * + + 2 * + z x *  z y * x y * + + 4 * + t 1 - t * 0.5 * *") ;
			//3
		(*shapefunc)[3] = Function("1 x - y - z - x * 4 * t 1 - t * 0.5 * *") ;
			//4
		(*shapefunc)[4] = Function("2 x * 1 - x * t 1 - t * 0.5 * *") ;
			//5
		(*shapefunc)[5] = Function("x y * t 1 - t * 0.5 * *") ;
			//6
		(*shapefunc)[6] = Function("2 y * 1 - y * t 1 - t * 0.5 * *") ;
			//7
		(*shapefunc)[7] = Function("y z * t 1 - t * 0.5 * *") ;
		    //8
		(*shapefunc)[8] = Function("1 x - y - z - y * 4 * t 1 - t * 0.5 * *") ;
			//9
		(*shapefunc)[9] = Function("x z * t 1 - t * 0.5 * *") ;
		    //10
		(*shapefunc)[10] = Function("2 z * 1 - z * 1 t t * - *") ;
			//11
		(*shapefunc)[11] = Function("1 x - y - z - z * 4 * 1 t t * - *") ;
			//12
		(*shapefunc)[12] = Function("1 z y x + + 3 * - z z* x x * y y * + + 2 * + z x *  z y * x y * + + 4 * + 1 t t * - *") ;
			//13
		(*shapefunc)[13] = Function("1 x - y - z - x * 4 * 1 t t * - *") ;
			//14
		(*shapefunc)[14] = Function("2 x * 1 - x * 1 t t * - *") ;
			//15
		(*shapefunc)[15] = Function("x y * 1 t t * - *") ;
			//16
		(*shapefunc)[16] = Function("2 y * 1 - y * 1 t t * - *") ;
			//17
		(*shapefunc)[17] = Function("y z * 1 t t * - *") ;
		    //18
		(*shapefunc)[18] = Function("1 x - y - z - y * 4 * 1 t t * - *") ;
			//19
		(*shapefunc)[19] = Function("x z * 1 t t * - *") ;
			//20
		(*shapefunc)[20] = Function("2 z * 1 - z * t 1 + t * 0.5 * *") ;
			//21
		(*shapefunc)[21] = Function("1 x - y - z - z * 4 * t 1 + t * 0.5 * *") ;
			//22
		(*shapefunc)[22] = Function("1 z y x + + 3 * - z z* x x * y y * + + 2 * + z x *  z y * x y * + + 4 * + t 1 + t * 0.5 * *") ;
			//23
		(*shapefunc)[23] = Function("1 x - y - z - x * 4 * t 1 + t * 0.5 * *") ;
			//24
		(*shapefunc)[24] = Function("2 x * 1 - x * t 1 + t * 0.5 * *") ;
			//25
		(*shapefunc)[25] = Function("x y * t 1 + t * 0.5 * *") ;
			//26
		(*shapefunc)[26] = Function("2 y * 1 - y * t 1 + t * 0.5 * *") ;
			//27
		(*shapefunc)[27] = Function("y z * t 1 + t * 0.5 * *") ;
		    //28
		(*shapefunc)[28] = Function("1 x - y - z - y * 4 * t 1 + t * 0.5 * *") ;
			//29
		(*shapefunc)[29] = Function("x z * t 1 + t * 0.5 * *") ;
	}
	else
	{
		std::cout << "this time-space order combination is not implemented" << std::endl ;
		assert(false) ;
	}
	assert(this->jacobianAtPoint(*this->Tetrahedron::getCenter()) > 0) ;
}
	
TetrahedralElement::TetrahedralElement(TetrahedralElement * parent, Tetrahedron * t)
{
	this->order =  parent->getOrder() ;
	
	this->shapefunc = parent->shapefunc ;
// 		std::copy(&parent->getShapeFunctions()[0], &parent->getShapeFunctions()[parent->getShapeFunctions().size()], &this->shapefunc[0]) ;
	for(size_t i =  0  ; i < this->size() ; i++)
	{
		delete &this->getPoint(i) ;
	}
	
	this->getInPoints().resize(t->getInPoints().size()) ;
	this->getBoundingPoints().resize(t->getBoundingPoints().size()) ;
	
	
	
	std::copy(&t->getInPoints()[0], &t->getInPoints()[t->getInPoints().size()],&this->getInPoints()[0] ) ;
	std::copy(&t->getBoundingPoints()[0], &t->getBoundingPoints()[t->getBoundingPoints().size()],&this->getBoundingPoints()[0] ) ;
}
	
std::vector<std::vector<Matrix> > TetrahedralElement::getElementaryMatrix() const 
{
	return std::vector<std::vector<Matrix> >() ;
}

std::vector<std::vector<Matrix> > TetrahedralElement::getNonLinearElementaryMatrix() const 
{
	return std::vector<std::vector<Matrix> >() ;
}

Vector TetrahedralElement::getForces() const 
{
	return Vector(0) ;
}

Vector TetrahedralElement::getNonLinearForces() const 
{
	return Vector(0) ;
}
	
void TetrahedralElement::refresh(const TetrahedralElement * parent)
{
	
	this->order =  parent->getOrder() ;
	this->shapefunc = parent->shapefunc ;
}

void TetrahedralElement::print()  const
{
	std::cout << "I am a tetrahedron..." << std::endl ;
}

const Function ElementarySurface::getXTransform() const
{
	return XTransform( this->getBoundingPoints(), this->getShapeFunctions()) ;
}

const Function ElementarySurface::getYTransform() const
{
	return YTransform( this->getBoundingPoints(), this->getShapeFunctions()) ;
}

const Function ElementarySurface::getTTransform() const
{
	return TTransform( this->getBoundingPoints(), this->getShapeFunctions()) ;
}

const Function ElementarySurface::getdXTransform(Variable v) const
{
	return dXTransform( this->getBoundingPoints(), this->getShapeFunctions(),v) ;
}

const Function ElementarySurface::getdYTransform(Variable v) const
{
	return dYTransform( this->getBoundingPoints(), this->getShapeFunctions(),v) ;
}

const Function ElementarySurface::getdTTransform(Variable v) const
{
	return dTTransform( this->getBoundingPoints(), this->getShapeFunctions(),v) ;
}

const double ElementarySurface::getdXTransform(Variable v, const Point p) const
{
	return dXTransform( this->getBoundingPoints(), this->getShapeFunctions(),v, p) ;
}

const double ElementarySurface::getdYTransform(Variable v, const Point p) const
{
	return dYTransform( this->getBoundingPoints(), this->getShapeFunctions(),v, p) ;
}

const double ElementarySurface::getdTTransform(Variable v, const Point p) const
{
	return dTTransform( this->getBoundingPoints(), this->getShapeFunctions(),v, p) ;
}

void ElementarySurface::setEnrichment(std::pair<size_t, Function>  p)
{
	bool unique = true ;
	for(size_t i = 0 ;  i < enrichfunc.size() ; i++)
	{
		if (getEnrichmentFunction(i).second.getDofID() == p.second.getDofID())
		{
			unique = false ;
			break ;
		}
	}
	if(unique)
		enrichfunc.push_back(p) ;
}

const  std::pair<size_t, Function> & ElementarySurface::getEnrichmentFunction(size_t i) const
{
	return this->enrichfunc[i];
}

std::pair<size_t, Function> & ElementarySurface::getEnrichmentFunction(size_t i)
{
	return this->enrichfunc[i];
}





const std::valarray< Function >  & ElementarySurface::getShapeFunctions() const
{
	return *this->shapefunc   ;
}


const std::vector< std::pair<size_t, Function>  > & ElementarySurface::getEnrichmentFunctions() const
{
	return this->enrichfunc;
}

std::vector< std::pair<size_t, Function>  >  & ElementarySurface::getEnrichmentFunctions()
{
	return this->enrichfunc;
}

// const std::vector<Function> ElementarySurface::getEnrichment(size_t i) const
// {
// 	std::vector<Function> ret ;
// 	for(size_t j = 0 ; j < enrichfunc.size() ; j++)
// 	{
// 		if(enrichfunc[j].first == i)
// 			ret.push_back(enrichfunc[j].second) ;
// 	}
// 	
// 	return ret ;
// }
// 
// 
// std::vector<Function> ElementarySurface::getEnrichment(size_t i) 
// {
// 	std::vector<Function> ret ;
// 	for(size_t j = 0 ; j < enrichfunc.size() ; j++)
// 	{
// 		if(enrichfunc[j].first == i)
// 			ret.push_back(enrichfunc[j].second) ;
// 	}
// 	
// 	return ret ;
// }


Point TetrahedralElement::inLocalCoordinates(const Point & p) const
{
	
	size_t factor = 1 ;
	
	if(order == QUADRATIC || order == QUADRATIC_TIME_LINEAR || order == QUADRATIC_TIME_QUADRATIC)
		factor = 2 ;
	
	Matrix S(4,4) ;
	S[0][0] = this->getBoundingPoint(factor*2).x ; 
	S[0][1] = this->getBoundingPoint(factor*3).x ;  
	S[0][2] = this->getBoundingPoint(0      ).x;
	S[0][3]=  this->getBoundingPoint(factor  ).x; 
	
	S[1][0] = this->getBoundingPoint(factor*2).y ; 
	S[1][1] = this->getBoundingPoint(factor*3).y ;  
	S[1][2] = this->getBoundingPoint(0      ).y ;
	S[1][3]=  this->getBoundingPoint(factor  ).y;

	S[2][0] = this->getBoundingPoint(factor*2).z ; 
	S[2][1] = this->getBoundingPoint(factor*3).z ;  
	S[2][2] = this->getBoundingPoint(0      ).z ;
	S[2][3]=  this->getBoundingPoint(factor  ).z;
	
	S[3][0] = 1 ; S[3][1] = 1 ;  S[3][2] = 1 ; S[3][3]=1;
	
	Vector v(4) ; 
	v[0] = p.x ;
	v[1] = p.y ;
	v[2] = p.z ;
	v[3] = 1 ;
	Vector coeff = inverse4x4Matrix(S) * v ;
	return Point(1,0,0)*coeff[0] + Point(0,1,0)*coeff[1] + Point(0,0,1)*coeff[2] + Point(0,0,0,p.t); 
}

Point HexahedralElement::inLocalCoordinates(const Point& p) const
{
	
	
	Matrix S(4,4) ;
	S[0][0] = this->getBoundingPoint(0).x;
	S[0][1] = this->getBoundingPoint(2).x;
	S[0][2] = this->getBoundingPoint(4).x;
	S[0][3]=  this->getBoundingPoint(1).x;
	
	
	S[1][0] = this->getBoundingPoint(0).y;
	S[1][1] = this->getBoundingPoint(2).y;
	S[1][2] = this->getBoundingPoint(4).y;
	S[1][3]=  this->getBoundingPoint(1).y;
	
	
	S[2][0] = this->getBoundingPoint(0).z;
	S[2][1] = this->getBoundingPoint(2).z;
	S[2][2] = this->getBoundingPoint(4).z;
	S[2][3]=  this->getBoundingPoint(1).z;
	
	
	S[3][0] = 1 ; S[3][1] = 1 ;  S[3][2] = 1 ; S[3][3]=1;
	
	Vector v(4) ; 
	v[0] = p.x ;
	v[1] = p.y ;
	v[2] = p.z ;
	v[3] = 1 ;
	
	Vector coeff = inverse4x4Matrix(S) * v ;
	
	return Point(-1,-1,-1)*coeff[0] + Point(-1,-1,1)*coeff[1] + Point(-1,1,-1)*coeff[2]+  Point(-1,-1,1)*coeff[3];

}


const Function &  ElementarySurface::getShapeFunction(size_t i) const
{
	return (*shapefunc)[i] ;
}

Function & ElementarySurface::getShapeFunction(size_t i) 
{
	return (*shapefunc)[i] ;
}

Order ElementarySurface::getOrder() const
{
	return order ;
}

const Function ElementaryVolume::getdXTransform(Variable v) const
{
	return dXTransform( this->getBoundingPoints(), this->getShapeFunctions(),v) ;
}

const Function ElementaryVolume::getdYTransform(Variable v) const
{
	return dYTransform( this->getBoundingPoints(), this->getShapeFunctions(),v) ;
}

const Function ElementaryVolume::getdZTransform(Variable v) const
{
	return dZTransform( this->getBoundingPoints(), this->getShapeFunctions(),v) ;
}

const Function ElementaryVolume::getdTTransform(Variable v) const
{
	return dTTransform( this->getBoundingPoints(), this->getShapeFunctions(),v) ;
}


const double ElementaryVolume::getdXTransform(Variable v, const Point p) const
{
	return dXTransform( this->getBoundingPoints(), this->getShapeFunctions(),v, p) ;
}

const double ElementaryVolume::getdYTransform(Variable v, const Point p) const
{
	return dYTransform( this->getBoundingPoints(), this->getShapeFunctions(),v, p) ;
}

const double ElementaryVolume::getdZTransform(Variable v, const Point p) const
{
	return dZTransform( this->getBoundingPoints(), this->getShapeFunctions(),v, p) ;
}

const double ElementaryVolume::getdTTransform(Variable v, const Point p) const
{
	return dTTransform( this->getBoundingPoints(), this->getShapeFunctions(),v, p) ;
}


ElementaryVolume::ElementaryVolume(bool f )  : isFather(f)
{
	this->behaviour = NULL ;
	this->nonlinbehaviour = NULL ;
	this->state = new ElementState(this) ;
}
////

ElementaryVolume::~ElementaryVolume()
{
	if(isFather)
		delete shapefunc ;
	delete this->state ;
}

Function ElementaryVolume::jacobian() const 
{
	
	Function xdxi = this->getXTransform().d(XI) ;
	Function ydxi = this->getYTransform().d(XI) ;
	Function zdxi = this->getZTransform().d(XI) ;
	Function xdeta = this->getXTransform().d(ETA) ;
	Function ydeta = this->getYTransform().d(ETA) ;
	Function zdeta = this->getZTransform().d(ETA) ;
	Function xdzeta = this->getXTransform().d(ZETA) ;
	Function ydzeta = this->getYTransform().d(ZETA) ;
	Function zdzeta = this->getZTransform().d(ZETA) ;
	
	Function ret = xdxi*ydeta*zdzeta + zdeta*xdzeta*ydxi + ydzeta*zdxi*xdeta  -
						xdxi*ydzeta*zdeta - xdeta*ydxi*zdzeta - xdzeta*ydeta*zdxi ;
	
	return ret ;
}

double ElementaryVolume::jacobianAtPoint(const Point p) const 
{
	
	if(order < CONSTANT_TIME_LINEAR)
	{
		double xdxi = this->getdXTransform(XI, p) ;
		double ydxi = this->getdYTransform(XI, p) ;
		double zdxi = this->getdZTransform(XI, p) ;
		double xdeta = this->getdXTransform(ETA, p) ;
		double ydeta = this->getdYTransform(ETA, p) ;
		double zdeta = this->getdZTransform(ETA, p) ;
		double xdzeta = this->getdXTransform(ZETA,p) ;
		double ydzeta = this->getdYTransform(ZETA,p) ;
		double zdzeta = this->getdZTransform(ZETA,p) ;
		
		return xdxi*ydeta*zdzeta + zdeta*xdzeta*ydxi + ydzeta*zdxi*xdeta  -
			xdxi*ydzeta*zdeta - xdeta*ydxi*zdzeta - xdzeta*ydeta*zdxi ;
	}
	else
	{
		Matrix  J0(4,4) ;
		
		double xdxi = this->getdXTransform(XI,p) ;
		double ydxi = this->getdYTransform(XI,p) ;
		double zdxi = this->getdZTransform(XI,p) ;
		double tdxi = this->getdTTransform(XI,p) ;
		
		double xdeta = this->getdXTransform(ETA,p) ;
		double ydeta = this->getdYTransform(ETA,p) ;
		double zdeta = this->getdZTransform(ETA,p) ;
		double tdeta = this->getdTTransform(ETA,p) ;
		
		double xdzeta = this->getdXTransform(ZETA,p) ;
		double ydzeta = this->getdYTransform(ZETA,p) ;
		double zdzeta = this->getdZTransform(ZETA,p) ;
		double tdzeta = this->getdTTransform(ZETA,p) ;
		
		double xdtheta = this->getdXTransform(TIME_VARIABLE,p) ;
		double ydtheta = this->getdYTransform(TIME_VARIABLE,p) ;
		double zdtheta = this->getdZTransform(TIME_VARIABLE,p) ;
		double tdtheta = this->getdTTransform(TIME_VARIABLE,p) ;
		
		J0[0][0] = xdxi ; J0[0][1] = ydxi ; J0[0][2] = zdxi ; J0[0][3] = tdxi; 
		J0[1][0] = xdeta ; J0[1][1] = ydeta ; J0[1][2] = zdeta ; J0[1][3] = tdeta;
		J0[2][0] = xdzeta ; J0[2][1] = ydzeta ; J0[2][2] = zdzeta ; J0[2][3] = tdzeta;
		J0[3][0] = xdtheta ; J0[3][1] = ydtheta ; J0[3][2] = zdtheta ; J0[3][3] = tdtheta;
		
		return det(J0) ;
	}
	
	
	VirtualMachine vm ;
	return vm.eval(jacobian(), p) ;
}

const Function ElementaryVolume::getXTransform() const
{
	return XTransform( this->getBoundingPoints(), this->getShapeFunctions()) ;
}

const Function  ElementaryVolume::getYTransform() const
{
	return YTransform( this->getBoundingPoints(), this->getShapeFunctions()) ;
}

const Function  ElementaryVolume::getZTransform() const
{
	return ZTransform( this->getBoundingPoints(), this->getShapeFunctions()) ;
}

const Function  ElementaryVolume::getTTransform() const
{
	return TTransform( this->getBoundingPoints(), this->getShapeFunctions()) ;
}

void ElementaryVolume::setEnrichment(std::pair<size_t, Function>  p)
{
	enrichfunc.push_back(p) ;
}

const std::pair<size_t, Function> & ElementaryVolume::getEnrichmentFunction(size_t i)  const
{
	return this->enrichfunc[i] ;
}

std::pair<size_t, Function> & ElementaryVolume::getEnrichmentFunction(size_t i) 
{
	return this->enrichfunc[i] ;
}

const std::valarray<Function  >  & ElementaryVolume::getShapeFunctions() const
{
	return *this->shapefunc ;
}
//const std::valarray< Function >  ElementaryVolume::getShapeFunctions() const
//{
//	return *this->shapefunc ;
//}

const std::vector<std::pair<size_t, Function> > ElementaryVolume::getDofs() const
{
	std::vector<std::pair<size_t, Function> > ret ;
	for (size_t i = 0 ; i < getShapeFunctions().size() ; i++)
	{
		ret.push_back(std::make_pair(getBoundingPoint(i).id, getShapeFunction(i))) ;
	}
	for (size_t i = 0 ; i < getEnrichmentFunctions().size() ; i++)
	{
		ret.push_back(std::make_pair(getEnrichmentFunction(i).second.getDofID(),getEnrichmentFunction(i).second)) ;
	}
	
	return ret ;
}

bool ElementaryVolume::isMoved() const
{
	return false;
}

Form * ElementaryVolume::getBehaviour() const
{
	return this->behaviour ;
}

const std::vector< size_t > ElementaryVolume::getDofIds() const
{
	
	std::vector<size_t> ret ;
	for (size_t i = 0 ; i < getBoundingPoints().size() ; i++)
	{
		if(getBoundingPoint(i).id >= 0)
			ret.push_back(getBoundingPoint(i).id) ;
		else
			std::cout << "negative ID, check numbering !" << std::endl ;
	}
	for (size_t i = 0 ; i < getEnrichmentFunctions().size() ; i++)
	{
		if(getEnrichmentFunction(i).second.getDofID() >= 0)
			ret.push_back(getEnrichmentFunction(i).second.getDofID()) ;
		else
			std::cout << "negative ID, check numbering !" << std::endl ;
	}
	
	return ret ;
}

const std::vector< std::pair<size_t, Function>  > & ElementaryVolume::getEnrichmentFunctions() const
{
	return this->enrichfunc;
}

std::vector< std::pair<size_t, Function>  > & ElementaryVolume::getEnrichmentFunctions() 
{
	return this->enrichfunc;
}

std::vector<std::vector<Matrix> > HexahedralElement::getElementaryMatrix() const 
	{
	
		std::vector<std::vector<Matrix > > mother ;
		std::valarray<std::pair<Point, double> > gp  = this->getGaussPoints(); 
		std::valarray<Matrix> Jinv ;
		std::vector<std::pair<size_t, Function> > dofs = getDofs() ;
/*	VirtualMachine vm ;*/
		if(getEnrichmentFunctions().size() > 0 /*&& getEnrichmentFunction(0).second.getIntegrationHint().size() > 0*/ )
		{

			Matrix J = this->getInverseJacobianMatrix( gp[0].first ) ;
			Jinv.resize(gp.size()) ;
			
			for(size_t i = 0 ; i < gp.size() ;  i++)
			{
				Jinv[i] = J ;	
			}
		}
		else
		{
			Matrix J = this->getInverseJacobianMatrix( gp[0].first ) ;
			Jinv.resize(gp.size()) ;
				
			for(size_t i = 0 ; i < gp.size() ;  i++)
			{
				Jinv[i] = J ;
			}
		}
		
		for(size_t i = 0 ; i < dofs.size() ; i++)
		{
			std::vector< Matrix > v_j ;
			
			for(size_t j = 0 ; j < dofs.size() ; j++)
			{
				v_j.push_back(Matrix(3,3)) ;
				
			}
			
			mother.push_back(v_j) ;
		}
		
		
		for(size_t i = 0 ; i < dofs.size() ; i++)
		{
			mother[i][i] =behaviour->apply(dofs[i].second, dofs[i].second,gp, Jinv) ;
			
			for(size_t j = i+1 ; j < dofs.size() ; j++)
			{
				mother[i][j] = behaviour->apply(dofs[i].second, dofs[j].second,gp, Jinv) ;
				mother[j][i] = mother[i][j].transpose() ;
				
			}
		}

// 		}
// 	}
		return mother ;
// 		}
		

	}

//const std::vector< std::pair<size_t, Function> >  ElementaryVolume::getEnrichmentFunctions() const
//{
//	return this->enrichfunc;
//}

// std::valarray< Point *> * ElementaryVolume::getBoundingPoints() const 
// {
// 	return this->getBoundingPoints() ;
// }
// 
// std::valarray< Point *> * ElementaryVolume::getInPoints() const 
// {
// 	return this->getInPoints() ;
// }

 Vector HexahedralElement::getForces() const 
	{

		std::valarray<std::pair<Point, double> > gp  ;
		std::valarray<Matrix> Jinv ;
		std::vector<std::pair<size_t, Function> > dofs = getDofs() ;
		Vector forces (dofs.size()*3);
/*	VirtualMachine vm ;*/
		if(getEnrichmentFunctions().size() > 0)
		{

			std::valarray<std::pair<Point, double> > gp_alt(this->getGaussPoints()) ;
			Matrix J = this->getInverseJacobianMatrix( gp_alt[0].first ) ;
			gp.resize(gp_alt.size()) ;
			Jinv.resize(gp_alt.size()) ;
			std::copy(&gp_alt[0], &gp_alt[gp_alt.size()], &gp[0]);
			
			for(size_t i = 0 ; i < gp.size() ;  i++)
			{
				Jinv[i] = J ;	
			}
		}
		else
		{
			std::valarray<std::pair<Point, double> > gp_alt(this->getGaussPoints()) ;
			Matrix J = this->getInverseJacobianMatrix( gp_alt[0].first ) ;
			gp.resize(gp_alt.size()) ;
			Jinv.resize(gp_alt.size()) ;
			std::copy(&gp_alt[0], &gp_alt[gp_alt.size()], &gp[0]);
				
			for(size_t i = 0 ; i < gp.size() ;  i++)
			{
				Jinv[i] = J ;
			}
		}
			
		for(size_t i = 0 ; i < dofs.size() ; i++)
		{
			for(size_t j = 0 ; j < dofs.size() ; j++)
			{
				Vector f = behaviour->getForces(this->getState(), dofs[i].second ,dofs[j].second,gp, Jinv) ;
				
				forces[i*3]+=f[0];
				forces[i*3+1]+=f[1];
				forces[i*3+2]+=f[2];
			}
		}

		return forces ;

	}




const Function  & ElementaryVolume::getShapeFunction(size_t i) const
{
	return (*shapefunc)[i] ;
}

Function & ElementaryVolume::getShapeFunction(size_t i) 
{
	return (*shapefunc)[i] ;
}
//const std::pair<size_t, Function> getEnrichmentFunction(size_t i) const ;

Matrix ElementaryVolume::getInverseJacobianMatrix(const Point & p) const
{
	if(order < CONSTANT_TIME_LINEAR)
	{
		Matrix  J0(3,3) ;
		
		double xdxi = this->getdXTransform(XI,p) ;
		double ydxi = this->getdYTransform(XI,p) ;
		double zdxi = this->getdZTransform(XI,p) ;
		
		double xdeta = this->getdXTransform(ETA,p) ;
		double ydeta = this->getdYTransform(ETA,p) ;
		double zdeta = this->getdZTransform(ETA,p) ;
		
		double xdzeta = this->getdXTransform(ZETA,p) ;
		double ydzeta = this->getdYTransform(ZETA,p) ;
		double zdzeta = this->getdZTransform(ZETA,p) ;
		
		J0[0][0] = xdxi ; J0[0][1] = ydxi ; J0[0][2] = zdxi ; 
		J0[1][0] = xdeta ; J0[1][1] = ydeta ; J0[1][2] = zdeta ;
		J0[2][0] = xdzeta ; J0[2][1] = ydzeta ; J0[2][2] = zdzeta ;
	
		invert3x3Matrix(J0) ;
		return J0 ;
	}
	else
	{
		Matrix  J0(4,4) ;
		
		double xdxi = this->getdXTransform(XI,p) ;
		double ydxi = this->getdYTransform(XI,p) ;
		double zdxi = this->getdZTransform(XI,p) ;
		double tdxi = this->getdTTransform(XI,p) ;
		
		double xdeta = this->getdXTransform(ETA,p) ;
		double ydeta = this->getdYTransform(ETA,p) ;
		double zdeta = this->getdZTransform(ETA,p) ;
		double tdeta = this->getdTTransform(ETA,p) ;
		
		double xdzeta = this->getdXTransform(ZETA,p) ;
		double ydzeta = this->getdYTransform(ZETA,p) ;
		double zdzeta = this->getdZTransform(ZETA,p) ;
		double tdzeta = this->getdTTransform(ZETA,p) ;
		
		double xdtheta = this->getdXTransform(TIME_VARIABLE,p) ;
		double ydtheta = this->getdYTransform(TIME_VARIABLE,p) ;
		double zdtheta = this->getdZTransform(TIME_VARIABLE,p) ;
		double tdtheta = this->getdTTransform(TIME_VARIABLE,p) ;
		
		J0[0][0] = xdxi ; J0[0][1] = ydxi ; J0[0][2] = zdxi ; J0[0][3] = tdxi; 
		J0[1][0] = xdeta ; J0[1][1] = ydeta ; J0[1][2] = zdeta ; J0[1][3] = tdeta;
		J0[2][0] = xdzeta ; J0[2][1] = ydzeta ; J0[2][2] = zdzeta ; J0[2][3] = tdzeta;
		J0[3][0] = xdtheta ; J0[3][1] = ydtheta ; J0[3][2] = zdtheta ; J0[3][3] = tdtheta;
		
		J0 = inverse4x4Matrix(J0) ;
		return J0 ;
	}

}

 Vector HexahedralElement::getNonLinearForces() const 
{
	return Vector(0);
}

 std::vector<std::vector<Mu::Matrix> >  HexahedralElement::getNonLinearElementaryMatrix()  const 
{
	return std::vector<std::vector<Mu::Matrix> >(0);
}

void ElementaryVolume::setBehaviour(Form * f)
{
	behaviour =  f ;
}

Order ElementaryVolume::getOrder() const
{
	return order ;
}


Function XTransform(const std::valarray<Mu::Point*> & points, const std::valarray<Function > & basis)
{
	Function ret ;
	Function der_x ;
	Function der_y ;
	Function der_z ;
	
	assert(points->size() == basis.size()) ;
	for(size_t i = 0 ; i < points.size() ; i++)
	{
		ret += basis[i]*points[i]->x ;
		der_x += basis[i].d(XI)*points[i]->x ;
		if(basis[0].getDerivatives().size() > 1) 
		{
			der_y += basis[i].d(ETA)*points[i]->x ;
			if(basis[0].getDerivatives().size() > 2) 
			{
				der_z += basis[i].d(ZETA)*points[i]->x ;
			}
		}
	}
	ret.getDerivatives().resize(basis[0].getDerivatives().size()) ;
	
	ret.getDerivatives()[0] = der_x ;
	if(basis[0].getDerivatives().size() > 1) 
	{
		ret.getDerivatives()[1] = der_y ;
		if(basis[0].getDerivatives().size() > 2) 
		{
			ret.getDerivatives()[2] = der_z ;
		}
	}
	
	return ret ;
}

Function YTransform(const std::valarray<Mu::Point*> & points, const std::valarray<Function > & basis)
{
	Function ret ;
	Function der_x ;
	Function der_y ;
	Function der_z ;
	
	assert(points.size() == basis.size()) ;
	for(size_t i = 0 ; i < points.size() ; i++)
	{
		ret += basis[i]*points[i]->y ;
		der_x += basis[i].d(XI)*points[i]->y ;
		der_y += basis[i].d(ETA)*points[i]->y ;
		if(basis[0].getDerivatives().size() > 2) 
		{
			der_z += basis[i].d(ZETA)*points[i]->y ;
		}
	}
	ret.getDerivatives().resize(basis[0].getDerivatives().size()) ;
	
	ret.getDerivatives()[0] = der_x ;

	ret.getDerivatives()[1] = der_y ;
	if(basis[0].getDerivatives().size() > 2) 
	{
		ret.getDerivatives()[2] = der_z ;
	}
	
	return ret ;
}

Function ZTransform(const std::valarray<Mu::Point*> & points, const std::valarray<Function > & basis)
{
	Function ret ;
	Function der_x ;
	Function der_y ;
	Function der_z ;
	
	assert(points.size() == basis.size()) ;
	for(size_t i = 0 ; i < points.size() ; i++)
	{
		ret += basis[i]*points[i]->z ;
		der_x += basis[i].d(XI)*points[i]->z ;
		der_y += basis[i].d(ETA)*points[i]->z ;
		der_z += basis[i].d(ZETA)*points[i]->z ;

	}
	
	ret.getDerivatives().resize(basis[0].getDerivatives().size()) ;
	ret.getDerivatives()[0] = der_x ;
	ret.getDerivatives()[1] = der_y ;
	ret.getDerivatives()[2] = der_z ;

	
	return ret ;
}

Function TTransform(const std::valarray<Mu::Point*> & points, const std::valarray<Function > & basis)
{
	Function ret ;
	Function der_x ;
	Function der_y ;
	Function der_z ;
	Function der_t ;
	
	assert(points.size() == basis.size()) ;
	for(size_t i = 0 ; i < points.size() ; i++)
	{
		ret += basis[i]*points[i]->t ;
		der_x += basis[i].d(XI)*points[i]->t ;
		der_y += basis[i].d(ETA)*points[i]->t ;
		der_z += basis[i].d(ZETA)*points[i]->t ;
		der_t += basis[i].d(ZETA)*points[i]->t ;
		
	}
	
	ret.getDerivatives().resize(basis[0].getDerivatives().size()) ;
	ret.getDerivatives()[0] = der_x ;
	ret.getDerivatives()[1] = der_y ;
	ret.getDerivatives()[2] = der_z ;
	ret.getDerivatives()[3] = der_t ;
	
	
	return ret ;
}

Mu::Function dXTransform(const std::valarray<Mu::Point*> & points ,const std::valarray< Mu::Function> &basis, Mu::Variable v)
{
	switch(v)
	{
	case XI:
		{
			Function der_x ;
			
			assert(points->size() == basis.size()) ;
			for(size_t i = 0 ; i < points.size() ; i++)
			{
				der_x += basis[i].d(XI)*points[i]->x ;
			}

			return der_x ;
		}
	case ETA:
		{
			Function der_y ;
			
			assert(points.size() == basis.size()) ;
			for(size_t i = 0 ; i < points.size() ; i++)
			{
				der_y += basis[i].d(ETA)*points[i]->x ;
			}
			
			return der_y ;
		}
	case ZETA:
		{
			Function der_z ;
			
			assert(points.size() == basis.size()) ;
			for(size_t i = 0 ; i < points.size() ; i++)
			{
				der_z += basis[i].d(ZETA)*points[i]->x ;
			}
			
			return der_z ;
		}
	case TIME_VARIABLE:
		{
			Function der_t ;
			
			assert(points.size() == basis.size()) ;
			for(size_t i = 0 ; i < points.size() ; i++)
			{
				der_t += basis[i].d(TIME_VARIABLE)*points[i]->x ;
			}
			
			return der_t ;
		}
	default:
		{
			std::cout << "*** dXTransform error ***" << std::endl ;
			return Function() ;
		}
	}
}

Mu::Function dYTransform(const std::valarray<Mu::Point*> & points ,const std::valarray< Mu::Function> &basis, Mu::Variable v)
{
	switch(v)
	{
	case XI:
		{
			Function der_x ;
			
			assert(points.size() == basis.size()) ;
			for(size_t i = 0 ; i < points.size() ; i++)
			{
				der_x += basis[i].d(XI)*points[i]->y ;
			}
			
			return der_x ;
		}
	case ETA:
		{
			Function der_y ;
			
			assert(points->size() == basis.size()) ;
			for(size_t i = 0 ; i < points.size() ; i++)
			{
				der_y += basis[i].d(ETA)*points[i]->y ;
			}
			
			return der_y ;
		}
	case ZETA:
		{
			Function der_z ;
			
			assert(points.size() == basis.size()) ;
			for(size_t i = 0 ; i < points.size() ; i++)
			{
				der_z += basis[i].d(ZETA)*points[i]->y ;
			}
			
			return der_z ;
		}
	case TIME_VARIABLE:
		{
			Function der_t ;
			
			assert(points.size() == basis.size()) ;
			for(size_t i = 0 ; i < points.size() ; i++)
			{
				der_t += basis[i].d(TIME_VARIABLE)*points[i]->y ;
			}
			
			return der_t ;
		}
	default:
		{
			std::cout << "*** dYTransform error ***" << std::endl ;
			return Function() ;
		}
	}
}

Mu::Function dZTransform(const std::valarray<Mu::Point*> & points ,const std::valarray< Mu::Function> &basis, Mu::Variable v)
{
	switch(v)
	{
	case XI:
		{
			Function der_x ;
			
			assert(points.size() == basis.size()) ;
			for(size_t i = 0 ; i < points.size() ; i++)
			{
				der_x += basis[i].d(XI)*points[i]->z ;
			}
			
			return der_x ;
		}
	case ETA:
		{
			Function der_y ;
			
			assert(points.size() == basis.size()) ;
			for(size_t i = 0 ; i < points.size() ; i++)
			{
				der_y += basis[i].d(ETA)*points[i]->z ;
			}
			
			return der_y ;
		}
	case ZETA:
		{
			Function der_z ;
			
			assert(points.size() == basis.size()) ;
			for(size_t i = 0 ; i < points.size() ; i++)
			{
				der_z += basis[i].d(ZETA)*points[i]->z ;
			}
			
			return der_z ;
		}
	case TIME_VARIABLE:
		{
			Function der_t ;
			
			assert(points.size() == basis.size()) ;
			for(size_t i = 0 ; i < points.size() ; i++)
			{
				der_t += basis[i].d(TIME_VARIABLE)*points[i]->t ;
			}
			
			return der_t ;
		}
	default:
		{
			std::cout << "*** dZTransform error ***" << std::endl ;
			return Function() ;
		}
	}
}

Mu::Function dTTransform(const std::valarray<Mu::Point*> & points ,const std::valarray< Mu::Function> &basis, Mu::Variable v)
{
	switch(v)
	{
	case XI:
		{
			Function der_x ;
			
			assert(points.size() == basis.size()) ;
			for(size_t i = 0 ; i < points.size() ; i++)
			{
				der_x += basis[i].d(XI)*points[i]->t ;
			}
			
			return der_x ;
		}
	case ETA:
		{
			Function der_y ;
			
			assert(points.size() == basis.size()) ;
			for(size_t i = 0 ; i < points.size() ; i++)
			{
				der_y += basis[i].d(ETA)*points[i]->t ;
			}
			
			return der_y ;
		}
	case ZETA:
		{
			Function der_z ;
			
			assert(points.size() == basis.size()) ;
			for(size_t i = 0 ; i < points.size() ; i++)
			{
				der_z += basis[i].d(ZETA)*points[i]->t ;
			}
			
			return der_z ;
		}
	case TIME_VARIABLE:
		{
			Function der_t ;
			
			assert(points.size() == basis.size()) ;
			for(size_t i = 0 ; i < points.size() ; i++)
			{
				der_t += basis[i].d(ZETA)*points[i]->t ;
			}
			
			return der_t ;
		}
	default:
		{
			std::cout << "*** dTTransform error ***" << std::endl ;
			return Function() ;
		}
	}
}

double dXTransform(const std::valarray<Mu::Point*> & points ,const std::valarray< Mu::Function> &basis, Mu::Variable v, const Point & p)
{
	switch(v)
	{
	case XI:
		{
			VirtualMachine vm ;
			double der_x = 0;
			
			for(size_t i = 0 ; i < basis.size() ; i++)
			{
				der_x += vm.deval(basis[i],XI,p)*points[i]->x ;
			}
			return der_x ;
		}
	case ETA:
		{
			VirtualMachine vm ;
			double der_y = 0 ;

			for(size_t i = 0 ; i < basis.size() ; i++)
			{
				der_y += vm.deval(basis[i],ETA, p)*points[i]->x ;
			}
			return der_y ;
		}
	case ZETA:
		{
			VirtualMachine vm ;
			double der_z = 0 ;
			
			for(size_t i = 0 ; i < basis.size() ; i++)
			{
				der_z += vm.deval(basis[i],ZETA,p)*points[i]->x ;
			}
			
			return der_z ;
		}
	case TIME_VARIABLE:
		{
			VirtualMachine vm ;
			double der_t = 0 ;
			
			for(size_t i = 0 ; i < basis.size() ; i++)
			{
				der_t += vm.deval(basis[i],TIME_VARIABLE,p)*points[i]->x ;
			}
			
			return der_t ;
		}
	default:
		{
			std::cout << "*** dXTransform error ***" << std::endl ;
			return 0 ;
		}
	}
}

double dYTransform(const std::valarray<Mu::Point*> & points ,const std::valarray< Mu::Function> &basis, Mu::Variable v, const Point & p)
{
	switch(v)
	{
	case XI:
		{
			
			VirtualMachine vm ;
			double der_x = 0;
// 			std::cout << "-----" << std::endl ;
// 			p.print() ;
			for(size_t i = 0 ; i < basis.size() ; i++)
			{
// 				vm.print(basis[i]) ;
// 				std::cout << vm.deval(basis[i],XI, p) << " X "<< points[i]->y << std::endl ;
				der_x += vm.deval(basis[i],XI,p)*points[i]->y ;
			}
// 			std::cout << "-----" << std::endl ;
			return der_x ;
		}
	case ETA:
		{
			VirtualMachine vm ;
			double der_y = 0 ;
			
			for(size_t i = 0 ; i < basis.size() ; i++)
			{
				der_y += vm.deval(basis[i],ETA, p)*points[i]->y ;
			}
			
			return der_y ;
		}
	case ZETA:
		{
			VirtualMachine vm ;
			double der_z = 0 ;
			
			for(size_t i = 0 ; i < basis.size() ; i++)
			{
				der_z += vm.deval(basis[i],ZETA,p)*points[i]->y ;
			}
			
			return der_z ;
		}
	case TIME_VARIABLE:
		{
			VirtualMachine vm ;
			double der_t = 0 ;
			

			for(size_t i = 0 ; i < basis.size() ; i++)
			{
				der_t += vm.deval(basis[i],TIME_VARIABLE,p)*points[i]->y ;
			}
			
			return der_t ;
		}
	default:
		{
			std::cout << "*** dYTransform error ***" << std::endl ;
			return 0 ;
		}
	}
}

double dZTransform(const std::valarray<Mu::Point*> & points ,const std::valarray< Mu::Function> &basis, Mu::Variable v, const Point & p)
{
	switch(v)
	{
	case XI:
		{
			VirtualMachine vm ;
			double der_x = 0;
			
			for(size_t i = 0 ; i < points.size() ; i++)
			{
				der_x += vm.deval(basis[i],XI,p)*points[i]->z ;
			}
			
			return der_x ;
		}
	case ETA:
		{
			VirtualMachine vm ;
			double der_y = 0 ;
			
			for(size_t i = 0 ; i < basis.size() ; i++)
			{
				der_y += vm.deval(basis[i],ETA, p)*points[i]->z ;
			}
			
			return der_y ;
		}
	case ZETA:
		{
			VirtualMachine vm ;
			double der_z = 0 ;
			
			for(size_t i = 0 ; i < points.size() ; i++)
			{
				der_z += vm.deval(basis[i],ZETA,p)*points[i]->z ;
			}
			
			return der_z ;
		}
	case TIME_VARIABLE:
		{
			VirtualMachine vm ;
			double der_t = 0 ;
			
			for(size_t i = 0 ; i < points.size() ; i++)
			{
				der_t += vm.deval(basis[i],TIME_VARIABLE,p)*points[i]->z ;
			}
			
			return der_t ;
		}
	default:
		{
			std::cout << "*** dZTransform error ***" << std::endl ;
			return 0 ;
		}
	}
}

double dTTransform(const std::valarray<Mu::Point*> & points ,const std::valarray< Mu::Function> &basis, Mu::Variable v, const Point & p)
{
	switch(v)
	{
	case XI:
		{
			VirtualMachine vm ;
			double der_x = 0;
			for(size_t i = 0 ; i < basis.size() ; i++)
			{
				der_x += vm.deval(basis[i],XI,p)*points[i]->t ;
			}
			return der_x ;
		}
	case ETA:
		{
			VirtualMachine vm ;
			double der_y = 0 ;
			
			for(size_t i = 0 ; i < points.size() ; i++)
			{
				der_y += vm.deval(basis[i],ETA, p)*points[i]->t ;
			}
			
			return der_y ;
		}
	case ZETA:
		{
			VirtualMachine vm ;
			double der_z = 0 ;
			
			for(size_t i = 0 ; i < points.size() ; i++)
			{
				der_z += vm.deval(basis[i],ZETA,p)*points[i]->t ;
			}
			
			return der_z ;
		}
	case TIME_VARIABLE:
		{
			VirtualMachine vm ;
			double der_t = 0 ;
			
			for(size_t i = 0 ; i < points.size() ; i++)
			{
				der_t += vm.deval(basis[i],TIME_VARIABLE,p)*points[i]->t ;
			}
			
			return der_t ;
		}
	default:
		{
			std::cout << "*** dTTransform error ***" << std::endl ;
			return 0 ;
		}
	}
}

std::valarray< std::pair<Point, double> > HexahedralElement::genGaussPoints() const
{
	size_t ordre=0 ;
	if(order <= LINEAR)
		ordre = 1 ;
	else if (order <= CUBIC)
		ordre = 8 ;
	else if (order <= QUINTIC)
		ordre = 27 ;
	
	std::valarray< std::pair<Point, double> > fin(ordre);
	
	if(order <= LINEAR)
	{
		fin[0] = std::pair<Point, double>(Point(0, 0, 0), 2.0) ;
	}
	else if (order <= CUBIC)
	{
		fin[0] = std::pair<Point, double>(Point(0.577350269189626,0.577350269189626 ,0.577350269189626 ), 1.0) ;
		fin[1] = std::pair<Point, double>(Point(-0.577350269189626,0.577350269189626 ,0.577350269189626 ), 1.0) ;
		fin[2] = std::pair<Point, double>(Point(0.577350269189626,-0.577350269189626 , 0.577350269189626 ), 1.0) ;
		fin[3] = std::pair<Point, double>(Point(0.577350269189626,0.577350269189626 , -0.577350269189626 ), 1.0) ;
		fin[4] = std::pair<Point, double>(Point(-0.577350269189626,-0.577350269189626 ,0.577350269189626 ), 1.0) ;
		fin[5] = std::pair<Point, double>(Point(-0.577350269189626,0.577350269189626 ,-0.577350269189626 ), 1.0) ;
		fin[6] = std::pair<Point, double>(Point(0.577350269189626,-0.577350269189626 ,-0.577350269189626 ), 1.0) ;
		fin[7] = std::pair<Point, double>(Point(-0.577350269189626,-0.577350269189626 ,-0.577350269189626 ), 1.0) ;
	}
	else if (order <= QUINTIC)
	{
		assert(false) ;
	}
	else
	{
		assert(false) ;
	}
	Function  J = this->jacobian() ;
	
	VirtualMachine vm ;
	for(size_t i = 0 ; i < fin.size() ; i++)
	{
		fin[i].second*=vm.eval(J, fin[i].first) ;
	}
	
	
	return fin ;
}
	
void HexahedralElement::computeCenter()
{
	this->Hexahedron::computeCenter() ;
}

HexahedralElement::HexahedralElement(Order order, bool f ) : ElementaryVolume( f)
{
	this->order = order ;
	
	if(order == LINEAR)
	{
		
		this->Hexahedron::sampleSurface(8) ;
		shapefunc = new std::valarray<Function>(8) ;
		Matrix zero ;
		std::valarray<Matrix> f1(zero, 2) ;
		f1[1][0][0] = -0.125 ;
		f1[0][1][0] = -0.125 ;
		f1[0][0][1] = -0.125 ;
		f1[0][0][0] = 0.125 ;
		f1[1][1][0] = 0.125 ;
		f1[0][1][1] = 0.125 ;
		f1[1][0][1] = 0.125 ;
		f1[1][1][1] =-0.125 ;
		
		std::valarray<Matrix> f2(zero,2) ;
		f2[1][0][0] = 0.125 ;
		f2[0][1][0] = -0.125 ;
		f2[0][0][1] = -0.125 ;
		f2[0][0][0] = 0.125 ;
		f2[1][1][0] = -0.125 ;
		f2[0][1][1] = 0.125 ;
		f2[1][0][1] = -0.125 ;
		f2[1][1][1] = 0.125 ;
	
		std::valarray<Matrix> f3(zero,2) ;
		f3[1][0][0] = 0.125 ;
		f3[0][1][0] = 0.125 ;
		f3[0][0][1] = -0.125 ;
		f3[0][0][0] = 0.125 ;
		f3[1][1][0] = 0.125 ;
		f3[0][1][1] = -0.125 ;
		f3[1][0][1] = -0.125 ;
		f3[1][1][1] = -0.125 ;
		
		std::valarray<Matrix> f4(zero,2) ;
		f4[1][0][0] = -0.125 ;
		f4[0][1][0] = 0.125 ;
		f4[0][0][1] = -0.125 ;
		f4[0][0][0] = 0.125 ;
		f4[1][1][0] = -0.125 ;
		f4[0][1][1] = -0.125 ;
		f4[1][0][1] = 0.125 ;
		f4[1][1][1] = 0.125 ;
			
		std::valarray<Matrix> f5(zero,2) ;
		f5[1][0][0] = -0.125 ;
		f5[0][1][0] = -0.125 ;
		f5[0][0][1] = 0.125 ;
		f5[0][0][0] = 0.125 ;
		f5[1][1][0] = 0.125 ;
		f5[0][1][1] = -0.125 ;
		f5[1][0][1] = -0.125 ;
		f5[1][1][1] = 0.125 ;
		
		std::valarray<Matrix> f6(zero,2) ;
		f6[1][0][0] = 0.125 ;
		f6[0][1][0] = -0.125 ;
		f6[0][0][1] = 0.125 ;
		f6[0][0][0] = 0.125 ;
		f6[1][1][0] = -0.125 ;
		f6[0][1][1] = -0.125 ;
		f6[1][0][1] = 0.125 ;
		f6[1][1][1] = -0.125 ;
			
		std::valarray<Matrix> f7(zero,2) ;
		f7[1][0][0] = 0.125 ;
		f7[0][1][0] = 0.125 ;
		f7[0][0][1] = 0.125 ;
		f7[0][0][0] = 0.125 ;
		f7[1][1][0] = 0.125 ;
		f7[0][1][1] = 0.125 ;
		f7[1][0][1] = 0.125 ;
		f7[1][1][1] = 0.125 ;
		
		std::valarray<Matrix> f8(zero,2) ;
		f8[1][0][0] = -0.125 ;
		f8[0][1][0] = 0.125 ;
		f8[0][0][1] = 0.125 ;
		f8[0][0][0] = 0.125 ;
		f8[1][1][0] = -0.125 ;
		f8[0][1][1] = 0.125 ;
		f8[1][0][1] = -0.125 ;
		f8[1][1][1] = -0.125 ;
			
			//0
		(*shapefunc)[0] = Function(f1) ;
			//1
		(*shapefunc)[1] = Function(f2) ;
			//2
		(*shapefunc)[2] = Function(f3) ;
			//3
		(*shapefunc)[3] = Function(f4) ;
			
			//4
		(*shapefunc)[4] = Function(f5) ;
			//5
		(*shapefunc)[5] = Function(f6) ;
			//6
		(*shapefunc)[6] = Function(f7) ;
			//7
		(*shapefunc)[7] = Function(f8) ;
		
		
	}
	else if(order == QUADRATIC)
	{
		assert(false) ;
	}
	else
	{
		assert(false) ;
	}
	assert(this->jacobianAtPoint(*this->Hexahedron::getCenter()) > 0) ;
	
	this->size_x  = 2 ;
	this->size_y  = 2 ;
	this->size_z  = 2 ;
}
	
HexahedralElement::HexahedralElement(HexahedralElement * parent,Hexahedron * t)
{
	this->order =  parent->getOrder() ;
	
	this->shapefunc = parent->shapefunc ;

	for(size_t i =  0  ; i < this->size() ; i++)
	{
		delete &this->getPoint(i) ;
	}
	
	this->getInPoints().resize(t->getInPoints().size()) ;
	this->getBoundingPoints().resize(t->getBoundingPoints().size()) ;
	
	std::copy(&t->getInPoints()[0], &t->getInPoints()[t->getInPoints().size()],&this->getInPoints()[0] ) ;
	std::copy(&t->getBoundingPoints()[0], &t->getBoundingPoints()[t->getBoundingPoints().size()],&this->getBoundingPoints()[0] ) ;
	
	this->size_x  = t->getXSize() ;
	this->size_y  = t->getYSize() ;
	this->size_z  = t->getZSize() ;
}
	
void HexahedralElement::refresh(const HexahedralElement * parent)
{
	
	this->order =  parent->getOrder() ;
	this->shapefunc = parent->shapefunc ;
}

void HexahedralElement::print()  const
{
	std::cout << "I am an hexahedron..." << std::endl ;
}

std::vector<int> validNeighbours(int i, int j, int k, int length)
{
	std::vector<int> ret ;
	for(int ii = i-1 ; ii < i+2 ; ii++)
	{
		for(int jj = j-1 ; jj < j+2 ; jj++)
		{
			for(int kk = k-1 ; kk < k+2 ; kk++)
			{
				if((ii >= 0 && ii < length)
				   && (jj >= 0&& jj < length)
				   && (kk >= 0&& kk < length)
				   && !( ii == i && jj== j && kk==k)
				  )
					ret.push_back(kk + jj*length + ii*length*length) ;
			}
		}
	}
	return ret ;
}

void computeNeighbourhoodForStructuredHexahedralMesh(std::vector<Mu::HexahedralElement *> & vec)
{
	int length = (int)round(std::pow(vec.size(), 1./3.)) ;
	
	for(int i = 0 ; i < length ; i++)
	{
		for(int j = 0 ; j < length ; j++)
		{
			for(int k = 0 ; k < length ; j++)
			{
				std::vector<int> neighbours = validNeighbours(i,j,k,length) ;
				int current_index = k + j*length + i*length*length ;
				
				for(size_t l = 0 ; l < neighbours.size() ; l++)
				{
					vec[current_index]->neighbourhood.push_back(vec[neighbours[l]]) ;
				}
			}
		}
	}
}
