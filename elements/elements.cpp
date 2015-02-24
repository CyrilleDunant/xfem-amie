
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
// Author: Ines Jaouadi <ines.jaouadi@epfl.ch>, (C) 2005-2007
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2009-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "elements.h"
#include "../physics/fracturecriteria/fracturecriterion.h"
#include "../utilities/tensor.h"

namespace Amie
{


ElementarySurface::~ElementarySurface()
{
	delete behaviour ;
}

TriElement::~TriElement()
{
	if(isFather && shapefunc)
	{
		delete shapefunc ;
		shapefunc = nullptr ;
	}
}

TetrahedralElement::~TetrahedralElement()
{
	if(isFather && shapefunc)
		delete shapefunc ;
}

ElementarySurface::ElementarySurface()
{
	nonlinbehaviour = nullptr ;
	behaviour = nullptr ;
	enrichmentUpdated = true ;
	behaviourUpdated = true ;
	cachedGps = nullptr ;
}

void ElementarySurface::setOrder(Order o)
{

	order = o ;

	switch(order)
	{
	case CONSTANT:
		{
			break ;
		}
	case LINEAR :
		{
			break ;
		}
	case QUADRATIC :
		{
			break ;
		}
	case CUBIC :
		{
			break ;
		}
	case QUADRIC :
		{
			break ;
		}
	case QUINTIC :
		{
			break ;
		}
	case QUADTREE_REFINED :
		{
			break ;
		}
	case CONSTANT_TIME_LINEAR :
		{
			timePlanes() = 2 ;
			break ;
		}
	case CONSTANT_TIME_QUADRATIC :
		{
			timePlanes() = 3 ;
			break ;
		}
	case LINEAR_TIME_LINEAR :
		{
			timePlanes() = 2 ;
			break ;
		}
	case LINEAR_TIME_QUADRATIC :
		{
			timePlanes() = 3 ;
			break ;
		}
	case QUADRATIC_TIME_LINEAR :
		{
			timePlanes() = 2 ;
			break ;
		}
	case QUADRATIC_TIME_QUADRATIC :
		{
			timePlanes() = 3 ;
			break ;
		}
	case CUBIC_TIME_LINEAR :
		{
			timePlanes() = 3 ;
			break ;
		}
	case CUBIC_TIME_QUADRATIC	 :
		{
			timePlanes() = 3 ;
			break ;
		}
	case QUADRIC_TIME_LINEAR :
		{
			timePlanes() = 2 ;
			break ;
		}
	case QUADRIC_TIME_QUADRATIC :
		{
			timePlanes() = 3 ;
			break ;
		}
	case QUINTIC_TIME_LINEAR :
		{
			timePlanes() = 2 ;
			break ;
		}
	case QUINTIC_TIME_QUADRATIC :
		{
			timePlanes() = 3 ;
			break ;
		}
	default:
		{
			break ;
		}
	}
}


Form * ElementarySurface::getBehaviour() const
{
/*	if(this->behaviour->param[0][0] == 0)
		this->behaviour->param[2048][2048] = 0 ;*/
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

void ElementarySurface::setBehaviour( Mesh<DelaunayTriangle,DelaunayTreeItem> * msh, Form * f)
{	
	bool init = false ;
	if(state)
	{	
 		init = this->getState().getDisplacements().size() > 0 ;
		delete state ;
	}
	state = f->createElementState( this ) ;
//	Form * old = behaviour ;
	
	delete behaviour ;
	behaviour = f ;
	if(init)
	{
		state->initialize(msh) ;
	}
//	delete old ;
}

void ElementarySurface::setNonLinearBehaviour(NonLinearForm * f)
{
	this->nonlinbehaviour = f ;
}

void ElementarySurface::step(double dt, const Vector * displacements)
{
	getState().step(dt, displacements) ;
	if(getBehaviour())
		getBehaviour()->updateElementState(dt, getState()) ;
}


void ElementarySurface::nonLinearStep(double dt, const Vector *displacements)
{
	getState().step(dt, displacements) ;
	
	if(getNonLinearBehaviour())
		getNonLinearBehaviour()->step(dt,getState(), -1) ;
}

void ElementaryVolume::step(double dt, const Vector *displacements)
{
	getState().step(dt, displacements) ;
	if(getBehaviour())
		getBehaviour()->updateElementState(dt, getState()) ;
}


void ElementaryVolume::nonLinearStep(double dt, const Vector * displacements)
{
	
	this->getState().step(dt, displacements) ;
	
	if(getNonLinearBehaviour())
		getNonLinearBehaviour()->step(dt, getState(), -1) ;
	
}

// const std::vector<std::pair<size_t,const Function &> > ElementarySurface::getDofs() const
// {
// 	std::vector<std::pair<size_t,const Function &> > ret ;
// 	for (size_t i = 0 ; i < getShapeFunctions().size() ; i++)
// 	{
// 		ret.push_back(std::make_pair(getBoundingPoint(i).getId(), getShapeFunction(i))) ;
// 	}
// 	for (size_t i = 0 ; i < getEnrichmentFunctions().size() ; i++)
// 	{
// 		ret.push_back(std::make_pair(getEnrichmentFunction(i).getDofID(),getEnrichmentFunction(i))) ;
// 	}
// 	
// 	return ret ;
// }

const std::vector< size_t > ElementarySurface::getDofIds() const
{
	std::vector<size_t> ret ;

	for (size_t i = 0 ; i < getShapeFunctions().size() ; i++)
	{
		if(getBoundingPoint(i).getId() >= 0)
			ret.push_back(getBoundingPoint(i).getId()) ;
		else
			std::cout << "negative point ID, check numbering !" << std::endl ;
	}
	for (size_t i = 0 ; i < getEnrichmentFunctions().size() ; i++)
	{
		if(getEnrichmentFunction(i).getDofID() >= 0)
			ret.push_back(getEnrichmentFunction(i).getDofID()) ;
		else
			std::cout << "negative enrichment ID, check numbering !" << std::endl ;
	}
	
	return ret ;
}



GaussPointArray gaussPointSet(Order order, const TriElement * t)
{
	size_t ordre = 0;
	std::valarray< std::pair<Point, double> > fin ;
	switch(order)
	{
	case CONSTANT:
	    break;
	case LINEAR:
		{
// 			ordre = 4 ;
// 			fin.resize(ordre);
// 			fin[0] = std::pair<Point, double>(Point(0.2, 0.2), 0.260416666666667) ;
// 			fin[1] = std::pair<Point, double>(Point(0.6, 0.2), 0.260416666666667) ;
// 			fin[2] = std::pair<Point, double>(Point(0.2, 0.6), 0.260416666666667) ;
// 			fin[3] = std::pair<Point, double>(Point(1./3., 1./3.), -0.28125) ;
// 			break ;
			
			ordre = 1 ;
			fin.resize(ordre);
			fin[0] = std::pair<Point, double>(Point(0.333333333333333, 0.333333333333333), 0.5) ;
			break ;
		}
	case QUADRATIC:
	{

			ordre = 4 ;
			fin.resize(ordre);
			fin[0] = std::pair<Point, double>(Point(0.2, 0.2), 0.260416666666667) ;
			fin[1] = std::pair<Point, double>(Point(0.6, 0.2), 0.260416666666667) ;
			fin[2] = std::pair<Point, double>(Point(0.2, 0.6), 0.260416666666667) ;
			fin[3] = std::pair<Point, double>(Point(1./3., 1./3.), -0.28125) ;
			break ;
	}
	case QUADTREE_REFINED:
	{

			ordre = 4 ;
			fin.resize(ordre);
			fin[0] = std::pair<Point, double>(Point(1./6., 1./6.), 0.125) ;
			fin[1] = std::pair<Point, double>(Point(1./6., 2./3.), 0.125) ;
			fin[2] = std::pair<Point, double>(Point(2./3., 1./6.), 0.125) ;
			fin[3] = std::pair<Point, double>(Point(1./3., 1./3.), 0.125) ;
			break ;
	}
	case CUBIC:
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
// 			ordre = 4 ;
// 			fin.resize(ordre);
// 			fin[0] = std::pair<Point, double>(Point(0.2, 0.2), 0.260416666666667) ;
// 			fin[1] = std::pair<Point, double>(Point(0.6, 0.2), 0.260416666666667) ;
// 			fin[2] = std::pair<Point, double>(Point(0.2, 0.6), 0.260416666666667) ;
// 			fin[3] = std::pair<Point, double>(Point(1./3., 1./3.), -0.28125) ;
// 			break ;
		}
	case QUADRIC:
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
			fin[1] = std::pair<Point, double>(Point(0.333333333333333, 0.333333333333333,0,0.577350269189626), 0.5) ;
			break;
		}
	case LINEAR_TIME_LINEAR:
		{
// 			ordre = 8 ;
// 			fin.resize(ordre);
// 			fin[0] = std::pair<Point, double>(Point(0.2, 0.2,0,-0.577350269189626), 0.2604166666666666667) ;
// 			fin[1] = std::pair<Point, double>(Point(0.6, 0.2,0,-0.577350269189626), 0.2604166666666666667) ;
// 			fin[2] = std::pair<Point, double>(Point(0.2, 0.6,0,-0.577350269189626), 0.2604166666666666667) ;
// 			fin[3] = std::pair<Point, double>(Point(0.333333333333333333333, 0.333333333333333333333,0,-0.577350269189626), -0.28125) ;
// 			fin[4] = std::pair<Point, double>(Point(0.2, 0.2,0,0.577350269189626), 0.26041666666666666667) ;
// 			fin[5] = std::pair<Point, double>(Point(0.6, 0.2,0,0.577350269189626), 0.26041666666666666667) ;
// 			fin[6] = std::pair<Point, double>(Point(0.2, 0.6,0,0.577350269189626), 0.26041666666666666667) ;
// 			fin[7] = std::pair<Point, double>(Point(0.333333333333333333333, 0.333333333333333333333,0,0.577350269189626), -0.28125) ;
// 			ordre = 1 ;
// 			fin.resize(ordre);
// 			fin[0] = std::pair<Point, double>(Point(0.333333333333333, 0.333333333333333,0,0), 1) ;
/* 			ordre = 3 ;
 			fin.resize(ordre);
 			fin[0] = std::pair<Point, double>(Point(0.333333333333333, 0.333333333333333,0.,-0.774596669241483), 0.5555555555555556*0.5) ;
 			fin[1] = std::pair<Point, double>(Point(0.333333333333333, 0.333333333333333,0.,0.), 0.8888888888888889*0.5) ;
 			fin[2] = std::pair<Point, double>(Point(0.333333333333333, 0.333333333333333,0., 0.774596669241483), 0.5555555555555556*0.5) ;*/
 			ordre = 1 ;
 			fin.resize(ordre);
 			fin[0] = std::pair<Point, double>(Point(0.333333333333333, 0.333333333333333,0.,0.), 1.0) ;
/*			ordre = 5 ;
			fin.resize(ordre);
			double a = 1./3.*sqrt(5.-2*sqrt(10./7.)) ;
			double b = 1./3.*sqrt(5.+2*sqrt(10./7.)) ;
			fin[0] = std::pair<Point, double>(Point(1./3., 1./3.,0,-a), (322.+13.*sqrt(70.))/1800.) ;
			fin[1] = std::pair<Point, double>(Point(1./3., 1./3.,0,-b), (322.-13.*sqrt(70.))/1800.) ;
			fin[2] = std::pair<Point, double>(Point(1./3., 1./3.,0,0), 128./450.) ;
			fin[3] = std::pair<Point, double>(Point(1./3., 1./3.,0,b), (322.-13.*sqrt(70.))/1800.) ;
			fin[4] = std::pair<Point, double>(Point(1./3., 1./3.,0,a), (322.+13.*sqrt(70.))/1800.) ;*/
			
			break ;
		}
	case LINEAR_TIME_QUADRATIC:
		{
 /*			ordre = 3 ;
 			fin.resize(ordre);
 			fin[0] = std::pair<Point, double>(Point(0.333333333333333, 0.333333333333333,0.,-0.774596669241483), 0.5555555555555556*0.5) ;
 			fin[1] = std::pair<Point, double>(Point(0.333333333333333, 0.333333333333333,0.,0.), 0.8888888888888889*0.5) ;
 			fin[2] = std::pair<Point, double>(Point(0.333333333333333, 0.333333333333333,0.,0.774596669241483), 0.5555555555555556*0.5) ;*/
/* 			ordre = 6 ;
 			fin.resize(ordre);
 			fin[0] = std::pair<Point, double>(Point(0.666666666666666, 0.166666666666667,0.,-0.577350269189626), 0.333333333333333333) ;
 			fin[1] = std::pair<Point, double>(Point(0.166666666666667, 0.666666666666666,0.,-0.577350269189626), 0.333333333333333333) ;
 			fin[2] = std::pair<Point, double>(Point(0.166666666666667, 0.166666666666667,0.,-0.577350269189626), 0.333333333333333333) ;
 			fin[3] = std::pair<Point, double>(Point(0.666666666666666, 0.166666666666667,0., 0.577350269189626), 0.333333333333333333) ;
 			fin[4] = std::pair<Point, double>(Point(0.166666666666667, 0.666666666666666,0., 0.577350269189626), 0.333333333333333333) ;
 			fin[5] = std::pair<Point, double>(Point(0.166666666666667, 0.166666666666667,0., 0.577350269189626), 0.333333333333333333) ;*/

			ordre = 5 ;
			fin.resize(ordre);
			double a = 1./3.*sqrt(5.-2*sqrt(10./7.)) ;
			double b = 1./3.*sqrt(5.+2*sqrt(10./7.)) ;
			fin[0] = std::pair<Point, double>(Point(1./3., 1./3.,0,-a), (322.+13.*sqrt(70.))/1800.) ;
			fin[1] = std::pair<Point, double>(Point(1./3., 1./3.,0,-b), (322.-13.*sqrt(70.))/1800.) ;
			fin[2] = std::pair<Point, double>(Point(1./3., 1./3.,0,0), 128./450.) ;
			fin[3] = std::pair<Point, double>(Point(1./3., 1./3.,0,b), (322.-13.*sqrt(70.))/1800.) ;
			fin[4] = std::pair<Point, double>(Point(1./3., 1./3.,0,a), (322.+13.*sqrt(70.))/1800.) ;

			break ;
		}
	case QUADRATIC_TIME_LINEAR:
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
	case QUADRATIC_TIME_QUADRATIC:
		{
			ordre = 8 ;
			fin.resize(ordre);
			fin[0] = std::pair<Point, double>(Point(0.2, 0.2,0,-0.577350269189626), 0.2604166666666666667) ;
			fin[1] = std::pair<Point, double>(Point(0.6, 0.2,0,-0.577350269189626), 0.2604166666666666667) ;
			fin[2] = std::pair<Point, double>(Point(0.2, 0.6,0,-0.577350269189626), 0.2604166666666666667) ;
			fin[3] = std::pair<Point, double>(Point(0.333333333333333333333, 0.333333333333333333333,0,-0.577350269189626), -0.28125) ;
			fin[4] = std::pair<Point, double>(Point(0.2, 0.2,0,0.577350269189626), 0.26041666666666666667) ;
			fin[5] = std::pair<Point, double>(Point(0.6, 0.2,0,0.577350269189626), 0.26041666666666666667) ;
			fin[6] = std::pair<Point, double>(Point(0.2, 0.6,0,0.577350269189626), 0.26041666666666666667) ;
			fin[7] = std::pair<Point, double>(Point(0.333333333333333333333, 0.333333333333333333333,0,0.577350269189626), -0.28125) ;

			break ;
			
			
/*			ordre = 14 ;
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
			break;*/
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

	if(t->moved)
	{
		for(size_t i = 0 ; i < fin.size() ; i++)
		{
			double j = t->jacobianAtPoint(fin[i].first);
			fin[i].second *= j;
		}
	}
	else
	{
// 		Matrix J ;
// 		getInverseJacobianMatrix(fin[0].first, J);
		double j = t->area() ;//1./det(J) ;
                if(order < CONSTANT_TIME_LINEAR)
                    j *= 2. ;
		for(size_t i = 0 ; i < fin.size() ; i++)
		{
			fin[i].second *= j;
		}
	}
	return GaussPointArray(fin, order)  ;
}

GaussPointArray gaussPointSet(Order order, const TetrahedralElement * t)
{
	size_t ordre=0 ;
	if(order == LINEAR)
		ordre = 1 ;
	else if(order == LINEAR_TIME_QUADRATIC  || order == LINEAR_TIME_LINEAR)
		ordre = 3 ;
	else if (order == CUBIC || order == QUADRATIC || order == CUBIC_TIME_LINEAR || order == QUADRATIC_TIME_LINEAR)
		ordre = 5 ;
	else if (order == QUADRIC || order == QUINTIC )
		ordre = 17 ;
	else
		ordre = 10 ;
	
	std::valarray< std::pair<Point, double> > fin(ordre);
	
	if(order == LINEAR )
	{
		fin[0] = std::pair<Point, double>(Point(0.25, 0.25, 0.25), 0.1666666666666667) ;
	}
	else if (order == CUBIC || order == QUADRATIC )
	{
		fin[0] = std::pair<Point, double>(Point(0.25, 0.25, 0.25), 0.1666666666666667*-.8) ;
		fin[1] = std::pair<Point, double>(Point(0.166666666666667, 0.166666666666667, 0.166666666666667), 0.1666666666666667*0.45) ;
		fin[2] = std::pair<Point, double>(Point(0.5, 0.166666666666667, 0.166666666666667), 0.1666666666666667*0.45) ;
		fin[3] = std::pair<Point, double>(Point(0.166666666666667, 0.5, 0.166666666666667), 0.1666666666666667*0.45) ;
		fin[4] = std::pair<Point, double>(Point(0.166666666666667, 0.166666666666667, 0.5), 0.1666666666666667*0.45) ;
	}
	else if (order == QUADRIC || order == QUINTIC )
	{
		fin[0] = std::pair<Point, double>(Point(0.25, 0.25, 0.25), 0.031403092789502) ;
		fin[1] = std::pair<Point, double>(Point(0.731636907957618, 0.089454364014127, 0.089454364014127), 0.011173097287728) ;
		fin[2] = std::pair<Point, double>(Point(0.089454364014127, 0.089454364014127, 0.089454364014127), 0.011173097287728) ;
		fin[3] = std::pair<Point, double>(Point(0.089454364014127, 0.731636907957618, 0.089454364014127), 0.011173097287728) ;
		fin[4] = std::pair<Point, double>(Point(0.089454364014127, 0.089454364014127, 0.731636907957618), 0.011173097287728) ;		
		fin[5] = std::pair<Point, double>(Point(0.132581099938466, 0.024540037929030, 0.421439431066252), 0.007547598727224) ;
		fin[6] = std::pair<Point, double>(Point(0.132581099938466, 0.421439431066252, 0.024540037929030), 0.007547598727224) ;
		fin[7] = std::pair<Point, double>(Point(0.132581099938466, 0.421439431066252, 0.421439431066252), 0.007547598727224) ;
		fin[8] = std::pair<Point, double>(Point(0.024540037929030, 0.132581099938466, 0.421439431066252), 0.007547598727224) ;
		fin[9] = std::pair<Point, double>(Point(0.024540037929030, 0.421439431066252, 0.132581099938466), 0.007547598727224) ;
		fin[10] = std::pair<Point, double>(Point(0.024540037929030, 0.421439431066252, 0.421439431066252), 0.007547598727224) ;
		fin[11] = std::pair<Point, double>(Point(0.421439431066252, 0.132581099938466, 0.024540037929030), 0.007547598727224) ;
		fin[12] = std::pair<Point, double>(Point(0.421439431066252, 0.132581099938466, 0.421439431066252), 0.007547598727224) ;
		fin[13] = std::pair<Point, double>(Point(0.421439431066252, 0.024540037929030, 0.132581099938466), 0.007547598727224) ;
		fin[14] = std::pair<Point, double>(Point(0.421439431066252, 0.024540037929030, 0.421439431066252), 0.007547598727224) ;
		fin[15] = std::pair<Point, double>(Point(0.421439431066252, 0.421439431066252, 0.132581099938466), 0.007547598727224) ;
		fin[16] = std::pair<Point, double>(Point(0.421439431066252, 0.421439431066252, 0.024540037929030), 0.007547598727224) ;
	}
	else if(order == LINEAR_TIME_LINEAR )
	{
		fin[0] = std::pair<Point, double>(Point(0.25, 0.25, 0.25,-std::sqrt(0.6)), 0.1666666666666667*5./9) ;
		fin[1] = std::pair<Point, double>(Point(0.25, 0.25, 0.25,0.), 0.1666666666666667*8./9) ;
		fin[2] = std::pair<Point, double>(Point(0.25, 0.25, 0.25,std::sqrt(0.6)), 0.1666666666666667*5./9) ;
	}
	else if(order == LINEAR_TIME_QUADRATIC )
	{
		fin[0] = std::pair<Point, double>(Point(0.25, 0.25, 0.25,-std::sqrt(0.6)), 0.1666666666666667*5./9) ;
		fin[1] = std::pair<Point, double>(Point(0.25, 0.25, 0.25,0.), 0.1666666666666667*8./9) ;
		fin[2] = std::pair<Point, double>(Point(0.25, 0.25, 0.25,std::sqrt(0.6)), 0.1666666666666667*5./9) ;
	}
	else if (order == CUBIC_TIME_LINEAR || order == QUADRATIC_TIME_LINEAR )
	{
		fin[0] = std::pair<Point, double>(Point(0.25, 0.25, 0.25), -0.5333333333333333333) ;
		fin[1] = std::pair<Point, double>(Point(0.16666666666666666667, 0.166666666666667, 0.1666666666666666667), 0.3) ;
		fin[2] = std::pair<Point, double>(Point(0.5, 0.166666666666666666667, 0.166666666666666666667), 0.3) ;
		fin[3] = std::pair<Point, double>(Point(0.16666666666666666667, 0.5, 0.166666666666667), 0.3) ;
		fin[4] = std::pair<Point, double>(Point(0.16666666666666666667, 0.1666666666666666667, 0.5), 0.3) ;
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
	
	if( false && t->isMoved())
	{
		for(size_t i = 0 ; i < fin.size() ; i++)
		{
			fin[i].second *= t->jacobianAtPoint(fin[i].first) ;
		}
	}
	else
	{
		double j = t->volume()*6. ;
		for(size_t i = 0 ; i < fin.size() ; i++)
		{
			fin[i].second *= j;
		}
	}
	return GaussPointArray(fin, order)  ;
}

const GaussPointArray & TriElement::genGaussPoints() 
{
	if(getCachedGaussPoints() && !moved)
		return *getCachedGaussPoints() ;
	
	size_t ordre = 0;
	std::valarray< std::pair<Point, double> > fin ;
	switch(order)
	{
	case CONSTANT:
	    break;
	case LINEAR:
		{
// 			ordre = 4 ;
// 			fin.resize(ordre);
// 			fin[0] = std::pair<Point, double>(Point(0.2, 0.2), 0.260416666666667) ;
// 			fin[1] = std::pair<Point, double>(Point(0.6, 0.2), 0.260416666666667) ;
// 			fin[2] = std::pair<Point, double>(Point(0.2, 0.6), 0.260416666666667) ;
// 			fin[3] = std::pair<Point, double>(Point(1./3., 1./3.), -0.28125) ;
// 			break ;
			
			ordre = 1 ;
			fin.resize(ordre);
			fin[0] = std::pair<Point, double>(Point(0.333333333333333, 0.333333333333333), 0.5) ;
			break ;
		}
	case QUADRATIC:
	{

			ordre = 4 ;
			fin.resize(ordre);
			fin[0] = std::pair<Point, double>(Point(0.2, 0.2), 0.260416666666667) ;
			fin[1] = std::pair<Point, double>(Point(0.6, 0.2), 0.260416666666667) ;
			fin[2] = std::pair<Point, double>(Point(0.2, 0.6), 0.260416666666667) ;
			fin[3] = std::pair<Point, double>(Point(1./3., 1./3.), -0.28125) ;
			break ;
	}
	case QUADTREE_REFINED:
	{

			ordre = 4 ;
			fin.resize(ordre);
			fin[0] = std::pair<Point, double>(Point(1./6., 1./6.), 0.125) ;
			fin[1] = std::pair<Point, double>(Point(1./6., 2./3.), 0.125) ;
			fin[2] = std::pair<Point, double>(Point(2./3., 1./6.), 0.125) ;
			fin[3] = std::pair<Point, double>(Point(1./3., 1./3.), 0.125) ;
			break ;
	}
	case CUBIC:
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
// 			ordre = 4 ;
// 			fin.resize(ordre);
// 			fin[0] = std::pair<Point, double>(Point(0.2, 0.2), 0.260416666666667) ;
// 			fin[1] = std::pair<Point, double>(Point(0.6, 0.2), 0.260416666666667) ;
// 			fin[2] = std::pair<Point, double>(Point(0.2, 0.6), 0.260416666666667) ;
// 			fin[3] = std::pair<Point, double>(Point(1./3., 1./3.), -0.28125) ;
// 			break ;
		}
	case QUADRIC:
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
			fin[1] = std::pair<Point, double>(Point(0.333333333333333, 0.333333333333333,0,0.577350269189626), 0.5) ;
			break;
		}
	case LINEAR_TIME_LINEAR:
		{
 			ordre = 3 ;
 			fin.resize(ordre);
 			fin[0] = std::pair<Point, double>(Point(0.333333333333333, 0.333333333333333,0.,-std::sqrt(0.6)), 5./18) ;
 			fin[1] = std::pair<Point, double>(Point(0.333333333333333, 0.333333333333333,0.,0.), 8./18) ;
 			fin[2] = std::pair<Point, double>(Point(0.333333333333333, 0.333333333333333,0.,std::sqrt(0.6)), 5./18) ;
			break ;
		}
	case LINEAR_TIME_QUADRATIC:
		{
 			ordre = 3 ;
 			fin.resize(ordre);
 			fin[0] = std::pair<Point, double>(Point(0.333333333333333, 0.333333333333333,0.,-std::sqrt(0.6)), 5./18) ;
 			fin[1] = std::pair<Point, double>(Point(0.333333333333333, 0.333333333333333,0.,0.), 8./18) ;
 			fin[2] = std::pair<Point, double>(Point(0.333333333333333, 0.333333333333333,0.,std::sqrt(0.6)), 5./18) ;
			break ;
		}
	case QUADRATIC_TIME_LINEAR:
		{
			ordre = 12 ;
			fin.resize(ordre);
			fin[0] = std::pair<Point, double>(Point(0.2, 0.2,0,-std::sqrt(0.6)), 0.260416666666667*10/18) ;
			fin[1] = std::pair<Point, double>(Point(0.6, 0.2,0,-std::sqrt(0.6)), 0.260416666666667*10/18) ;
			fin[2] = std::pair<Point, double>(Point(0.2, 0.6,0,-std::sqrt(0.6)), 0.260416666666667*10/18) ;
			fin[3] = std::pair<Point, double>(Point(1./3., 1./3.,-std::sqrt(0.6)), -0.28125*10/18) ;
			fin[4] = std::pair<Point, double>(Point(0.2, 0.2,0,0.), 0.260416666666667*16/18) ;
			fin[5] = std::pair<Point, double>(Point(0.6, 0.2,0,0.), 0.260416666666667*16/18) ;
			fin[6] = std::pair<Point, double>(Point(0.2, 0.6,0,0.), 0.260416666666667*16/18) ;
			fin[7] = std::pair<Point, double>(Point(1./3., 1./3.,0,0.), -0.28125*16/18) ;
			fin[8] = std::pair<Point, double>(Point(0.2, 0.2,0,std::sqrt(0.6)), 0.260416666666667*10/18) ;
			fin[9] = std::pair<Point, double>(Point(0.6, 0.2,0,std::sqrt(0.6)), 0.260416666666667*10/18) ;
			fin[10] = std::pair<Point, double>(Point(0.2, 0.6,0,std::sqrt(0.6)), 0.260416666666667*10/18) ;
			fin[11] = std::pair<Point, double>(Point(1./3., 1./3.,0,std::sqrt(0.6)), -0.28125*10/18) ;
			break ;
		}
	case QUADRATIC_TIME_QUADRATIC:
		{
			ordre = 12 ;
			fin.resize(ordre);
			fin[0] = std::pair<Point, double>(Point(0.2, 0.2,0,-std::sqrt(0.6)), 0.260416666666667*10/18) ;
			fin[1] = std::pair<Point, double>(Point(0.6, 0.2,0,-std::sqrt(0.6)), 0.260416666666667*10/18) ;
			fin[2] = std::pair<Point, double>(Point(0.2, 0.6,0,-std::sqrt(0.6)), 0.260416666666667*10/18) ;
			fin[3] = std::pair<Point, double>(Point(1./3., 1./3.,-std::sqrt(0.6)), -0.28125*10/18) ;
			fin[4] = std::pair<Point, double>(Point(0.2, 0.2,0,0.), 0.260416666666667*16/18) ;
			fin[5] = std::pair<Point, double>(Point(0.6, 0.2,0,0.), 0.260416666666667*16/18) ;
			fin[6] = std::pair<Point, double>(Point(0.2, 0.6,0,0.), 0.260416666666667*16/18) ;
			fin[7] = std::pair<Point, double>(Point(1./3., 1./3.,0,0.), -0.28125*16/18) ;
			fin[8] = std::pair<Point, double>(Point(0.2, 0.2,0,std::sqrt(0.6)), 0.260416666666667*10/18) ;
			fin[9] = std::pair<Point, double>(Point(0.6, 0.2,0,std::sqrt(0.6)), 0.260416666666667*10/18) ;
			fin[10] = std::pair<Point, double>(Point(0.2, 0.6,0,std::sqrt(0.6)), 0.260416666666667*10/18) ;
			fin[11] = std::pair<Point, double>(Point(1./3., 1./3.,0,std::sqrt(0.6)), -0.28125*10/18) ;
			break ;
			
			
/*			ordre = 14 ;
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
			break;*/
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

	if(moved)
	{
		for(size_t i = 0 ; i < fin.size() ; i++)
		{
			double j = jacobianAtPoint(fin[i].first);
			fin[i].second *= j;
		}
	}
	else
	{
// 		Matrix J ;
// 		getInverseJacobianMatrix(fin[0].first, J);
		double j = area() ;//1./jacobianAtPoint(fin[0].first) ;
                if(order < CONSTANT_TIME_LINEAR)
                    j *= 2. ;
		for(size_t i = 0 ; i < fin.size() ; i++)
		{
			fin[i].second *= j;
		}
	}
	setCachedGaussPoints(new GaussPointArray(fin, order))  ;
	return *getCachedGaussPoints() ;
} ;


TriElement::TriElement( Point * p0,  Point * p1,  Point * p2) : Triangle(p0, p1, p2), moved(false) 
{ 
	isFather = false ;
	setOrder(LINEAR) ; 
	shapefunc = nullptr ; //new std::valarray<Function>(Function(),0) ;
/*	Matrix xi(2,2) ; xi[1][0] = 1 ;
	Matrix eta(2,2) ; eta[0][1] = 1 ;
	Matrix one(2,2) ; one[0][0] = 1 ;
//0
	(*shapefunc)[0] = Function(eta) ;
//1
	(*shapefunc)[1] = Function(one-xi-eta) ;
//2
	(*shapefunc)[2] = Function(xi) ;*/
};
	
TriElement::TriElement(Order order_ ): moved(false) 
{
	isFather = true ;
	setOrder( order_ );
	
	switch(order_)
	{
	case CONSTANT :
		{
			std::cout << "element order not implemented" << std::endl ;
			exit(0) ;
// 			shapefunc = new std::valarray<Function>(Function(), 1) ;
// 			Matrix m(1, 1) ;
// 			m[0][0] = 1 ;
// 			(*shapefunc)[0] = Function(m) ;
			break ;
		}
	case LINEAR :
		{
			shapefunc = new std::valarray<Function>(Function(),3) ;
// 			Matrix xi(2,2) ; xi[1][0] = 1 ;
// 			Matrix eta(2,2) ; eta[0][1] = 1 ;
// 			Matrix one(2,2) ; one[0][0] = 1 ;
			Function zero("0") ;
			Function one("1") ;
			Function mone("-1") ;
		//0
			(*shapefunc)[0] = Function("y") ;
			(*shapefunc)[0].setNumberOfDerivatives(2) ;
			(*shapefunc)[0].setDerivative( XI, zero) ;
			(*shapefunc)[0].setDerivative( ETA, one) ;
		//1
			(*shapefunc)[1] = Function("1 x - y -") ;
			(*shapefunc)[1].setNumberOfDerivatives(2) ;
			(*shapefunc)[1].setDerivative( XI, mone) ;
			(*shapefunc)[1].setDerivative( ETA, mone) ;
		//2
			(*shapefunc)[2] = Function("x") ;
			(*shapefunc)[2].setNumberOfDerivatives(2) ;
			(*shapefunc)[2].setDerivative( ETA, zero) ;
			(*shapefunc)[2].setDerivative( XI, one) ;
			break ;
		}
	case QUADRATIC :
		{
			shapefunc = new std::valarray<Function>(Function(),6) ;

			(*shapefunc)[0] = Function("y 2 ^ 2 * y -") ;
			(*shapefunc)[0].setNumberOfDerivatives(2) ;
			Function d("4 y * 1 -") ;
      (*shapefunc)[0].setDerivative( ETA, d) ;
			d = Function("0") ;
       (*shapefunc)[0].setDerivative( XI, d) ;


		//1
			(*shapefunc)[1] = Function("y 4 * y x 4 * * - y y 4 * * -") ;
			(*shapefunc)[1].setNumberOfDerivatives(2) ;
			d = Function("4 x 4 * - y 8 * -") ;
			(*shapefunc)[1].setDerivative( ETA, d) ;
			d = Function("y -4 *") ;
      (*shapefunc)[1].setDerivative( XI, d) ;

		//2
			(*shapefunc)[2] = Function("1 x 3 * - y 3 * - x y 4 * * + x x 2 * * + y y 2 * * +") ;
			(*shapefunc)[2].setNumberOfDerivatives(2) ;
			d = Function("x 4 * y 4 * + 3 -") ;
			(*shapefunc)[2].setDerivative( ETA, d ) ;
			d = Function("y 4 * x 4 * + 3 -") ;
			(*shapefunc)[2].setDerivative( XI, d) ;

		//3
			(*shapefunc)[3] = Function("x 4 * x x 4 * * - x y 4 * * -") ;
			(*shapefunc)[3].setNumberOfDerivatives(2) ;
			d = Function("x -4 * ") ;
			(*shapefunc)[3].setDerivative( ETA, d) ;
			d = Function("4 8 x * - y 4 * -") ;
			(*shapefunc)[3].setDerivative( XI, d) ;

		//4
			(*shapefunc)[4] = Function("x x 2 * * x -") ;
			(*shapefunc)[4].setNumberOfDerivatives(2) ;
			d = Function("0") ;
			(*shapefunc)[4].setDerivative( ETA, d) ;
			d = Function("4 x * 1 -") ;
			(*shapefunc)[4].setDerivative( XI, d) ;

		//5
			(*shapefunc)[5] = Function("x y 4 * *") ;
			(*shapefunc)[5].setNumberOfDerivatives(2) ;
			d = Function("4 x *") ;
			(*shapefunc)[5].setDerivative( ETA, d) ;
			d = Function("4 y *") ;
			(*shapefunc)[5].setDerivative( XI, d) ;

			break ;
		}
	case CONSTANT_TIME_LINEAR :
		{
			std::cout << "element order not implemented" << std::endl ;
			exit(0) ;
// 			shapefunc = new std::valarray<Function>(Function(),2) ;
// 			Matrix m(1, 2) ;
// 			std::valarray<Matrix> v(m,1) ; v[0] = m ;
// 			std::valarray<std::valarray<Matrix> > vv(v,1) ;
// 			vv[0] = v ;
// 			
// 			vv[0][0][0][1] = 1 ;
// 			(*shapefunc)[0] = Function(m) ;
// 			vv[0][0][0][0] = 1 ;
// 			vv[0][0][0][1] = -1 ;
// 			(*shapefunc)[1] = Function(m) ;
			break ;
		}
	case CONSTANT_TIME_QUADRATIC :
		{
			std::cout << "element order not implemented" << std::endl ;
			exit(0) ;
// 			shapefunc = new std::valarray<Function>(Function(),3) ;
// 			Matrix m(1, 3) ;
// 			std::valarray<Matrix> v(m,1) ; v[0] = m ;
// 			std::valarray<std::valarray<Matrix> > vv(v,1) ;
// 			vv[0] = v ;
// 			
// 			vv[0][0][0][2] = .5 ;
// 			vv[0][0][0][1] = -.5 ;
// 			(*shapefunc)[0] = Function(m) ;
// 			vv[0][0][0][0] = 1 ;
// 			vv[0][0][0][1] = 0 ;
// 			vv[0][0][0][2] = -1 ;
// 			(*shapefunc)[1] = Function(m) ;
// 			vv[0][0][0][0] = 0.5 ;
// 			vv[0][0][0][1] = 0 ;
// 			vv[0][0][0][2] = 0.5 ;
// 			(*shapefunc)[1] = Function(m) ;
			break ;
		}
	case LINEAR_TIME_LINEAR :
		{
			shapefunc = new std::valarray<Function>(Function(),6) ;

			Function z2("0") ;

			Function z1("0") ;
			z1.setNumberOfDerivatives(4) ;
			for(int i = 0 ; i < 4 ; i++)
			{
				z1.setDerivative( (const Variable) i, z2) ;
			}
			
			Function zero("0") ;
			zero.setNumberOfDerivatives(4) ;
			for(int i = 0 ; i < 4 ; i++)
			{
				zero.setDerivative( (const Variable) i, z1) ;
			}	
			Function one = zero +1 ;
			Function mone = zero-1 ;
			
			Function half("0.5") ;
			Function halfm("-0.5") ;
			
			Function t0("0.5 t 0.5 * -") ;
// 			t0 *= 0.5 ;
			t0.setNumberOfDerivatives(4) ;
			t0.setDerivative( XI, zero) ;
			t0.setDerivative( ETA, zero) ;
			t0.setDerivative( ZETA, zero) ;
			t0.setDerivative( TIME_VARIABLE, halfm) ;

// 			Function t0m = t0 * -1. ;
// 			t0m.setNumberOfDerivatives(4) ;
// 			t0m.setDerivative( XI, zero) ;
// 			t0m.setDerivative( ETA, zero) ;
// 			t0m.setDerivative( ZETA, zero) ;
// 			t0m.setDerivative( TIME_VARIABLE, half) ;
			
			Function t1("0.5 t 0.5 * +") ;
// 			t1 *= 0.5 ;
			t1.setNumberOfDerivatives(4) ;
			t1.setDerivative( XI, zero) ;
			t1.setDerivative( ETA, zero) ;
			t1.setDerivative( ZETA, zero) ;
			t1.setDerivative( TIME_VARIABLE, half) ;
			
// 			Function t1m = t1 * -1. ;
// 			t1m.setNumberOfDerivatives(4) ;
// 			t1m.setDerivative( XI, zero) ;
// 			t1m.setDerivative( ETA, zero) ;
// 			t1m.setDerivative( ZETA, zero) ;
// 			t1m.setDerivative( TIME_VARIABLE, halfm) ;
			
				Function s0("y") ;
				Function s1("1 x - y -") ;
				Function s2("x") ;
				s0.setNumberOfDerivatives(4) ;
				s0.setDerivative( XI, zero) ;
				s0.setDerivative( ETA, one) ;
				s0.setDerivative( ZETA, zero) ;
				s0.setDerivative( TIME_VARIABLE, zero) ;
				s1.setNumberOfDerivatives(4) ;
				s1.setDerivative( XI, mone) ;
				s1.setDerivative( ETA, mone) ;
				s1.setDerivative( ZETA, zero) ;
				s1.setDerivative( TIME_VARIABLE, zero) ;
				s2.setNumberOfDerivatives(4) ;
				s2.setDerivative( XI, one) ;
				s2.setDerivative( ETA, zero) ;
				s2.setDerivative( ZETA, zero) ;
				s2.setDerivative( TIME_VARIABLE, zero) ;
			
		//0			
			(*shapefunc)[0] = s0*t0 ;
			(*shapefunc)[1] = s1*t0 ;
			(*shapefunc)[2] = s2*t0 ;
			(*shapefunc)[3] = s0*t1 ;
			(*shapefunc)[4] = s1*t1 ;
			(*shapefunc)[5] = s2*t1 ;
			
// 			std::cout << VirtualMachine().eval((*shapefunc)[3], Point(0,1,0,1) ) << std::endl ;
// 			std::cout << VirtualMachine().deval((*shapefunc)[3], XI, Point(0,1,0,1) ) << std::endl  ;
// 			std::cout << VirtualMachine().deval((*shapefunc)[3], ETA, Point(0,1,0,1) ) << std::endl  ;
// 			std::cout << VirtualMachine().deval((*shapefunc)[3], TIME_VARIABLE, Point(0,1,0,1) ) << std::endl  ;
// 			exit(0) ;
			
// 			(*shapefunc)[0].setNumberOfDerivatives(4) ;
// 			(*shapefunc)[0].setDerivative( XI, zero) ;
// 			(*shapefunc)[0].setDerivative( ETA, t0) ;
		//1
// 			(*shapefunc)[1] = Function("1 x - y - 0.5 0.5 t * - *") ;
// 			(*shapefunc)[1].setNumberOfDerivatives(4) ;
// 			(*shapefunc)[1].setDerivative( XI, t0m) ;
// 			(*shapefunc)[1].setDerivative( ETA, t0m) ;
// 			//2
// 			(*shapefunc)[2] = Function("x 0.5 0.5 t * - *") ;
// 			(*shapefunc)[2].setNumberOfDerivatives(4) ;
// 			(*shapefunc)[2].setDerivative( XI, t0) ;
// 			(*shapefunc)[2].setDerivative( ETA, zero) ;
// 			//3
// 			(*shapefunc)[3] = Function("y 0.5 0.5 t * + *") ;
// 			(*shapefunc)[3].setNumberOfDerivatives(4) ;
// 			(*shapefunc)[3].setDerivative( XI, zero) ;
// 			(*shapefunc)[3].setDerivative( ETA, t1) ;
// 		//4
// 			(*shapefunc)[4] = Function("1 x - y - 0.5 0.5 t * + *") ;
// 			(*shapefunc)[4].setNumberOfDerivatives(4) ;
// 			(*shapefunc)[4].setDerivative( XI, t1m) ;
// 			(*shapefunc)[4].setDerivative( ETA, t1m) ;
// 		//5
// 			(*shapefunc)[5] = Function("x 0.5 0.5 t * + *") ;
// 			(*shapefunc)[5].setNumberOfDerivatives(4) ;
// 			(*shapefunc)[5].setDerivative( XI, t1) ;
// 			(*shapefunc)[5].setDerivative( ETA, zero) ;
			
			
			break ;
		}
		case LINEAR_TIME_QUADRATIC :
		{
			shapefunc = new std::valarray<Function>(Function(),9) ;

			Function z2("0") ;

			Function z1("0") ;
			z1.setNumberOfDerivatives(4) ;
			for(int i = 0 ; i < 4 ; i++)
			{
				z1.setDerivative( (const Variable) i, z2) ;
			}
			
			Function zero("0") ;
			zero.setNumberOfDerivatives(4) ;
			for(int i = 0 ; i < 4 ; i++)
			{
				zero.setDerivative( (const Variable) i, z1) ;
			}	
			Function one = zero +1 ;
			Function mone = zero-1 ;
			
			Function s0("y") ;
			Function s1("1 x - y -") ;
			Function s2("x") ;
			
			Function t0("t 2 ^ t - 0.5 *") ;			
			Function t1("1 t 2 ^ -") ;
			Function t2("t 2 ^ t + 0.5 *") ;

			Function tt0("t 0.5 -") ;
			Function tt1("t -2 *") ; 
			Function tt2("t 0.5 +") ;
			
			Function ttt0 = one ;
			Function ttt1 = mone *2 ;
			Function ttt2 = one ;
			
			s0.setNumberOfDerivatives(4) ;
			s0.setDerivative( XI, zero) ;
			s0.setDerivative( ETA, one) ;
			s0.setDerivative( ZETA, zero) ;
			s0.setDerivative( TIME_VARIABLE, zero) ;
			s1.setNumberOfDerivatives(4) ;
			s1.setDerivative( XI, mone) ;
			s1.setDerivative( ETA, mone) ;
			s1.setDerivative( ZETA, zero) ;
			s1.setDerivative( TIME_VARIABLE, zero) ;
			s2.setNumberOfDerivatives(4) ;
			s2.setDerivative( XI, one) ;
			s2.setDerivative( ETA, zero) ;
			s2.setDerivative( ZETA, zero) ;
			s2.setDerivative( TIME_VARIABLE, zero) ;

			tt0.setNumberOfDerivatives(4) ;
			tt0.setDerivative( XI, zero) ;
			tt0.setDerivative( ETA, zero) ;
			tt0.setDerivative( ZETA, zero) ;
			tt0.setDerivative( TIME_VARIABLE, ttt0 ) ;
			tt1.setNumberOfDerivatives(4) ;
			tt1.setDerivative( XI, zero) ;
			tt1.setDerivative( ETA, zero) ;
			tt1.setDerivative( ZETA, zero) ;
			tt1.setDerivative( TIME_VARIABLE, ttt1 ) ;
			tt2.setNumberOfDerivatives(4) ;
			tt2.setDerivative( XI, zero) ;
			tt2.setDerivative( ETA, zero) ;
			tt2.setDerivative( ZETA, zero) ;
			tt2.setDerivative( TIME_VARIABLE, ttt2 ) ;
			
			t0.setNumberOfDerivatives(4) ;
			t0.setDerivative( XI, zero) ;
			t0.setDerivative( ETA, zero) ;
			t0.setDerivative( ZETA, zero) ;
			t0.setDerivative( TIME_VARIABLE, tt0 ) ;
			t1.setNumberOfDerivatives(4) ;
			t1.setDerivative( XI, zero) ;
			t1.setDerivative( ETA, zero) ;
			t1.setDerivative( ZETA, zero) ;
			t1.setDerivative( TIME_VARIABLE, tt1 ) ;
			t2.setNumberOfDerivatives(4) ;
			t2.setDerivative( XI, zero) ;
			t2.setDerivative( ETA, zero) ;
			t2.setDerivative( ZETA, zero) ;
			t2.setDerivative( TIME_VARIABLE, tt2 ) ;	
			
		//0
			(*shapefunc)[0] = s0*t0 ;
// 			std::cout << VirtualMachine().eval( (*shapefunc)[0].d( TIME_VARIABLE), Point(0,1,0,-1) ) << std::endl ;
//			(*shapefunc)[0].d(TIME_VARIABLE) = Function("y t 0.5 - *") ;
		//1
			(*shapefunc)[1] = s1*t0 ;
//			(*shapefunc)[1].d(TIME_VARIABLE) = Function("1 x - y - t 0.5 - *") ;
		//2
			(*shapefunc)[2] = s2*t0 ;
//			(*shapefunc)[2].d(TIME_VARIABLE) = Function("x t 0.5 - *") ;
		//3
			(*shapefunc)[3] = s0*t1 ;
//			(*shapefunc)[3].d(TIME_VARIABLE) = Function("y t 2 * *") ;
		//4
			(*shapefunc)[4] = s1*t1 ;
//			(*shapefunc)[4].d(TIME_VARIABLE) = Function("1 x - y - t 2 * *") ;
		//5
			(*shapefunc)[5] = s2*t1 ;
//			(*shapefunc)[5].d(TIME_VARIABLE) = Function("x t 2 * *") ;
		//6
			(*shapefunc)[6] = s0*t2 ;
//			(*shapefunc)[6].d(TIME_VARIABLE) = Function("y t 0.5 + *") ;
		//7
			(*shapefunc)[7] = s1*t2 ;
// 			(*shapefunc)[7].d(TIME_VARIABLE) = Function("1 x - y - t 0.5 + *") ;
		//8
			(*shapefunc)[8] = s2*t2 ;
// 			(*shapefunc)[8].d(TIME_VARIABLE) = Function("x t 0.5 + *") ;
			
			break ;
		}
		case QUADRATIC_TIME_LINEAR:
		{
			shapefunc = new std::valarray<Function>(Function(),12) ;
			
			Function zero("0") ;
			Function t0("0.5 t 0.5 * -") ;
			Function half("0.5") ;
			Function halfm("-0.5") ;
// 			t0 *= 0.5 ;
			t0.setNumberOfDerivatives(4) ;
			t0.setDerivative( XI, zero) ;
			t0.setDerivative( ETA, zero) ;
			t0.setDerivative( ZETA, zero) ;
			t0.setDerivative( TIME_VARIABLE, halfm) ;

// 			Function t0m = t0 * -1. ;
// 			t0m.setNumberOfDerivatives(4) ;
// 			t0m.setDerivative( XI, zero) ;
// 			t0m.setDerivative( ETA, zero) ;
// 			t0m.setDerivative( ZETA, zero) ;
// 			t0m.setDerivative( TIME_VARIABLE, half) ;
			
			Function t1("0.5 t 0.5 * +") ;
// 			t1 *= 0.5 ;
			t1.setNumberOfDerivatives(4) ;
			t1.setDerivative( XI, zero) ;
			t1.setDerivative( ETA, zero) ;
			t1.setDerivative( ZETA, zero) ;
			t1.setDerivative( TIME_VARIABLE, half) ;

			(*shapefunc)[0] = Function("y 2 ^ 2 * y -") ;
			(*shapefunc)[0].setNumberOfDerivatives(4) ;
			Function d("4 y * 1 -") ;
      (*shapefunc)[0].setDerivative( ETA, d) ;
			d = Function("0") ;
       (*shapefunc)[0].setDerivative( XI, d) ;
			d = Function("0") ;
       (*shapefunc)[0].setDerivative( ZETA, d) ;
       (*shapefunc)[0].setDerivative( TIME_VARIABLE, d) ;
			 (*shapefunc)[0] =  (*shapefunc)[0]*t0 ;

		//1
			(*shapefunc)[1] = Function("y 4 * y x 4 * * - y y 4 * * -") ;
			(*shapefunc)[1].setNumberOfDerivatives(4) ;
			d = Function("4 x 4 * - y 8 * -") ;
			(*shapefunc)[1].setDerivative( ETA, d) ;
			d = Function("y -4 *") ;
      (*shapefunc)[1].setDerivative( XI, d) ;
				d = Function("0") ;
       (*shapefunc)[1].setDerivative( ZETA, d) ;
       (*shapefunc)[1].setDerivative( TIME_VARIABLE, d) ;
			 (*shapefunc)[1] =  (*shapefunc)[1]*t0 ;

		//2
			(*shapefunc)[2] = Function("1 x 3 * - y 3 * - x y 4 * * + x x 2 * * + y y 2 * * +") ;
			(*shapefunc)[2].setNumberOfDerivatives(4) ;
			d = Function("x 4 * y 4 * + 3 -") ;
			(*shapefunc)[2].setDerivative( ETA, d ) ;
			d = Function("y 4 * x 4 * + 3 -") ;
			(*shapefunc)[2].setDerivative( XI, d) ;
						d = Function("0") ;
       (*shapefunc)[2].setDerivative( ZETA, d) ;
       (*shapefunc)[2].setDerivative( TIME_VARIABLE, d) ;
			 (*shapefunc)[2] =  (*shapefunc)[2]*t0 ;

		//3
			(*shapefunc)[3] = Function("x 4 * x x 4 * * - x y 4 * * -") ;
			(*shapefunc)[3].setNumberOfDerivatives(4) ;
			d = Function("x -4 * ") ;
			(*shapefunc)[3].setDerivative( ETA, d) ;
			d = Function("4 8 x * - y 4 * -") ;
			(*shapefunc)[3].setDerivative( XI, d) ;
						d = Function("0") ;
       (*shapefunc)[3].setDerivative( ZETA, d) ;
       (*shapefunc)[3].setDerivative( TIME_VARIABLE, d) ;
			 (*shapefunc)[3] =  (*shapefunc)[3]*t0 ;

		//4
			(*shapefunc)[4] = Function("x x 2 * * x -") ;
			(*shapefunc)[4].setNumberOfDerivatives(4) ;
			d = Function("0") ;
			(*shapefunc)[4].setDerivative( ETA, d) ;
			d = Function("4 x * 1 -") ;
			(*shapefunc)[4].setDerivative( XI, d) ;
			d = Function("0") ;
       (*shapefunc)[4].setDerivative( ZETA, d) ;
       (*shapefunc)[4].setDerivative( TIME_VARIABLE, d) ;
			 (*shapefunc)[4] =  (*shapefunc)[4]*t0 ;

		//5
			(*shapefunc)[5] = Function("x y 4 * *") ;
			(*shapefunc)[5].setNumberOfDerivatives(4) ;
			d = Function("4 x *") ;
			(*shapefunc)[5].setDerivative( ETA, d) ;
			d = Function("4 y *") ;
			(*shapefunc)[5].setDerivative( XI, d) ;
						d = Function("0") ;
       (*shapefunc)[5].setDerivative( ZETA, d) ;
       (*shapefunc)[5].setDerivative( TIME_VARIABLE, d) ;
			 (*shapefunc)[5] =  (*shapefunc)[5]*t0 ;
			
			(*shapefunc)[6] = Function("y 2 ^ 2 * y -") ;
			(*shapefunc)[6].setNumberOfDerivatives(4) ;
			 d = Function("4 y * 1 -") ;
      (*shapefunc)[6].setDerivative( ETA, d) ;
			d = Function("0") ;
       (*shapefunc)[6].setDerivative( XI, d) ;
			d = Function("0") ;
       (*shapefunc)[6].setDerivative( ZETA, d) ;
       (*shapefunc)[6].setDerivative( TIME_VARIABLE, d) ;
			 (*shapefunc)[6] =  (*shapefunc)[6]*t1 ;

		//1
			(*shapefunc)[7] = Function("y 4 * y x 4 * * - y y 4 * * -") ;
			(*shapefunc)[7].setNumberOfDerivatives(4) ;
			d = Function("4 x 4 * - y 8 * -") ;
			(*shapefunc)[7].setDerivative( ETA, d) ;
			d = Function("y -4 *") ;
      (*shapefunc)[7].setDerivative( XI, d) ;
						d = Function("0") ;
       (*shapefunc)[7].setDerivative( ZETA, d) ;
       (*shapefunc)[7].setDerivative( TIME_VARIABLE, d) ;
			 (*shapefunc)[7] =  (*shapefunc)[7]*t1 ;

		//2
			(*shapefunc)[8] = Function("1 x 3 * - y 3 * - x y 4 * * + x x 2 * * + y y 2 * * +") ;
			(*shapefunc)[8].setNumberOfDerivatives(4) ;
			d = Function("x 4 * y 4 * + 3 -") ;
			(*shapefunc)[8].setDerivative( ETA, d ) ;
			d = Function("y 4 * x 4 * + 3 -") ;
			(*shapefunc)[8].setDerivative( XI, d) ;
						d = Function("0") ;
       (*shapefunc)[8].setDerivative( ZETA, d) ;
       (*shapefunc)[8].setDerivative( TIME_VARIABLE, d) ;
			 (*shapefunc)[8] =  (*shapefunc)[8]*t1 ;

		//3
			(*shapefunc)[9] = Function("x 4 * x x 4 * * - x y 4 * * -") ;
			(*shapefunc)[9].setNumberOfDerivatives(4) ;
			d = Function("x -4 * ") ;
			(*shapefunc)[9].setDerivative( ETA, d) ;
			d = Function("4 8 x * - y 4 * -") ;
			(*shapefunc)[9].setDerivative( XI, d) ;
						d = Function("0") ;
       (*shapefunc)[9].setDerivative( ZETA, d) ;
       (*shapefunc)[9].setDerivative( TIME_VARIABLE, d) ;
			 (*shapefunc)[9] =  (*shapefunc)[9]*t1 ;

		//4
			(*shapefunc)[10] = Function("x x 2 * * x -") ;
			(*shapefunc)[10].setNumberOfDerivatives(4) ;
			d = Function("0") ;
			(*shapefunc)[10].setDerivative( ETA, d) ;
			d = Function("4 x * 1 -") ;
			(*shapefunc)[10].setDerivative( XI, d) ;
			d = Function("0") ;
       (*shapefunc)[10].setDerivative( ZETA, d) ;
       (*shapefunc)[10].setDerivative( TIME_VARIABLE, d) ;
			 (*shapefunc)[10] =  (*shapefunc)[10]*t1 ;

		//5
			(*shapefunc)[11] = Function("x y 4 * *") ;
			(*shapefunc)[11].setNumberOfDerivatives(4) ;
			d = Function("4 x *") ;
			(*shapefunc)[11].setDerivative( ETA, d) ;
			d = Function("4 y *") ;
			(*shapefunc)[11].setDerivative( XI, d) ;
						d = Function("0") ;
       (*shapefunc)[11].setDerivative( ZETA, d) ;
       (*shapefunc)[11].setDerivative( TIME_VARIABLE, d) ;
			 (*shapefunc)[11] =  (*shapefunc)[11]*t1 ;
			

			break ;
		}
	case QUADRATIC_TIME_QUADRATIC:
		{
			shapefunc = new std::valarray<Function>(Function(),18) ;
			Function t0("t 2 ^ t - 0.5 *") ;
			Function t1("1 t 2 ^ -") ;
			Function t2("t 2 ^ t + 0.5 *") ;
			Function x0("y 2 ^ 2 * y -") ;
			x0.setNumberOfDerivatives(2) ;
			Function d("4 y * 1 -") ;
      x0.setDerivative( ETA, d) ;
			d = Function("0") ;
      x0.setDerivative( XI, d) ;
			d = Function("0") ;
      x0.setDerivative( ZETA, d) ;
      x0.setDerivative( TIME_VARIABLE, d) ;

		//1
			Function x1("y 4 * y x 4 * * - y y 4 * * -") ;
			x1.setNumberOfDerivatives(4) ;
			d = Function("4 x 4 * - y 8 * -") ;
			x1.setDerivative( ETA, d) ;
			d = Function("y 4 *") ;
      x1.setDerivative( XI, d) ;
			d = Function("0") ;
      x1.setDerivative( ZETA, d) ;
      x1.setDerivative( TIME_VARIABLE, d) ;

		//2
			Function x2("1 x 3 * - y 3 * - x y 4 * * + x x 2 * * + y y 2 * * +") ;
			x2.setNumberOfDerivatives(4) ;
			d = Function("3 x 4 * + y 4 * +") ;
			x2.setDerivative( ETA, d ) ;
			d = Function("3 y 4 * + x 4 * +") ;
			x2.setDerivative( XI, d) ;
			d = Function("0") ;
      x2.setDerivative( ZETA, d) ;
      x2.setDerivative( TIME_VARIABLE, d) ;

		//3
			Function x3("x 4 * x x 4 * * - x y 4 * * -") ;
			x3.setNumberOfDerivatives(4) ;
			d = Function("x 4 *") ;
			x3.setDerivative( ETA, d) ;
			d = Function("4 8 x * - y 4 * -") ;
			x3.setDerivative( XI, d) ;
			d = Function("0") ;
      x3.setDerivative( ZETA, d) ;
      x3.setDerivative( TIME_VARIABLE, d) ;

		//4
			Function x4("x x 2 * * x -") ;
			x4.setNumberOfDerivatives(4) ;
			d = Function("0") ;
			x4.setDerivative( ETA, d) ;
			d = Function("4 x * 1 -") ;
			x4.setDerivative( XI, d) ;
			d = Function("0") ;
      x4.setDerivative( ZETA, d) ;
      x4.setDerivative( TIME_VARIABLE, d) ;

		//5
			Function x5("x y 4 * *") ;
			x5.setNumberOfDerivatives(4) ;
			d = Function("4 x *") ;
			x5.setDerivative( ETA, d) ;
			d = Function("4 y *") ;
			x5.setDerivative( XI, d) ;
			d = Function("0") ;
      x5.setDerivative( ZETA, d) ;
      x5.setDerivative( TIME_VARIABLE, d) ;

			(*shapefunc)[0] = x0*t0;
			(*shapefunc)[1] = x1*t0;
			(*shapefunc)[2] = x2*t0;
			(*shapefunc)[3] = x3*t0;
			(*shapefunc)[4] = x4*t0;
			(*shapefunc)[5] = x5*t0;

			(*shapefunc)[6] = x0*t1; 
			(*shapefunc)[7] = x1*t1;
			(*shapefunc)[8] = x2*t1;
			(*shapefunc)[9] = x3*t1;
			(*shapefunc)[10] = x4*t1;
			(*shapefunc)[11] = x5*t1;

			(*shapefunc)[12] = x0*t2;
			(*shapefunc)[13] = x1*t2;
			(*shapefunc)[14] = x2*t2;
			(*shapefunc)[15] = x3*t2;
			(*shapefunc)[16] = x4*t2;
			(*shapefunc)[17] = x5*t2;
			break ;
		}
	default:
		{
			assert(false) ;
			break ;
		}
	}
	
	for(size_t i = 0 ; i < this->Triangle::getBoundingPoints().size() ; i++)
		this->Triangle::getBoundingPoint(i).setId(-1)  ;
	
}
	
void TriElement::refresh(const TriElement * parent)
{	
	if(!parent)
		return ;
	setOrder( parent->getOrder() );
	
// 	if(shapefunc && isFather)
// 	{
// 		isFather = false ;
// 		delete shapefunc ;
// 	}
	shapefunc = parent->shapefunc ;
// 	for(size_t i = 0 ; i < getBoundingPoints().size() ; i++)
// 	{
// 		(*parent->shapefunc)[i].setPoint(&getBoundingPoint(i)) ;
// 	}
}
	
std::valarray<std::valarray<Matrix> > & TriElement::getElementaryMatrix() 
{
	return cachedElementaryMatrix ;
}

std::valarray<std::valarray<Matrix> > & TriElement::getViscousElementaryMatrix() 
{
	return cachedViscousElementaryMatrix ;
}

std::valarray<std::valarray<Matrix> > TriElement::getNonLinearElementaryMatrix(Vector * state) 
{
	return std::valarray<std::valarray<Matrix> >() ;
}
	
Function TriElement::jacobian() const 
{
	
	Function xdxi = this->getdXTransform(XI) ;
	Function ydxi = this->getdYTransform(XI) ;
	Function xdeta = this->getdXTransform(ETA) ;
	Function ydeta = this->getdYTransform(ETA) ;
	
	Function ret = ydeta*xdxi - ydxi*xdeta  ;

	if(order >= CONSTANT_TIME_LINEAR)
	{
	    Function tdxi = this->getdTTransform(XI) ;
	    Function tdeta = this->getdTTransform(ETA) ;
	    Function xdtau = this->getdXTransform(TIME_VARIABLE) ;
	    Function ydtau = this->getdYTransform(TIME_VARIABLE) ;
	    Function tdtau = this->getdTTransform(TIME_VARIABLE) ;

	    ret = xdxi*ydeta*tdtau + tdeta*xdtau*ydxi + ydtau*tdxi*xdeta  -
				    xdxi*ydtau*tdeta - xdeta*ydxi*tdtau - xdtau*ydeta*tdxi ;
	}
	
	return ret ;
	
}
	
double  TriElement::jacobianAtPoint(const Amie::Point& p) const 
{
//	p.print() ;
	if(order < CONSTANT_TIME_LINEAR)
	{
		double xdxi = 0 ;//this->getdXTransform(XI,p) ;
		double ydxi = 0 ;//this->getdYTransform(XI,p) ;
		
		double xdeta = 0 ;//this->getdXTransform(ETA,p) ;
		double ydeta = 0 ;//this->getdYTransform(ETA,p) ;

		VirtualMachine vm ;
		TriElement father(order) ;
		for(size_t i = 0 ; i < getBoundingPoints().size() ; i++)
		{
			double dxi = vm.deval(father.getShapeFunction(i), XI, p) ;
			double deta = vm.deval(father.getShapeFunction(i), ETA, p) ;
			
			xdxi += dxi*getBoundingPoint(i).getX() ;
			ydxi += dxi*getBoundingPoint(i).getY() ;

			xdeta += deta*getBoundingPoint(i).getX() ;
			ydeta += deta*getBoundingPoint(i).getY() ;
		}
		
		return ydeta*xdxi - ydxi*xdeta ;
	}
	else
	{
		double xdxi = 0 ;//this->getdXTransform(XI,p) ;
		double ydxi = 0 ;//this->getdYTransform(XI,p) ;
		double tdxi = 0 ;
		
		double xdeta = 0 ;//this->getdXTransform(ETA,p) ;
		double ydeta = 0 ;//this->getdYTransform(ETA,p) ;
		double tdeta = 0 ;

		double xdtau = 0 ;
		double ydtau = 0 ;
		double tdtau = 0 ;

		VirtualMachine vm ;
		TriElement father(order) ;
		for(size_t i = 0 ; i < getBoundingPoints().size() ; i++)
		{
			double dxi = vm.deval(father.getShapeFunction(i), XI, p) ;
			double deta = vm.deval(father.getShapeFunction(i), ETA, p) ;
			double dtau = vm.deval(father.getShapeFunction(i),TIME_VARIABLE,p) ;
			
			xdxi += dxi*getBoundingPoint(i).getX() ;
			ydxi += dxi*getBoundingPoint(i).getY() ;
			tdxi += dxi*getBoundingPoint(i).getT() ;

			xdeta += deta*getBoundingPoint(i).getX() ;
			ydeta += deta*getBoundingPoint(i).getY() ;
			tdeta += deta*getBoundingPoint(i).getT() ;

			xdtau += dtau*getBoundingPoint(i).getX() ;
			ydtau += dtau*getBoundingPoint(i).getY() ;
			tdtau += dtau*getBoundingPoint(i).getT() ;

		}
/*		double zdxi = this->getdTTransform(XI, p) ;
		double zdeta = this->getdTTransform(ETA, p) ;
		double xdzeta = this->getdXTransform(TIME_VARIABLE,p) ;
		double ydzeta = this->getdYTransform(TIME_VARIABLE,p) ;
		double zdzeta = this->getdTTransform(TIME_VARIABLE,p) ;*/
		
		return  xdxi*ydeta*tdtau + tdeta*xdtau*ydxi + ydtau*tdxi*xdeta  -
			xdxi*ydtau*tdeta - xdeta*ydxi*tdtau - xdtau*ydeta*tdxi ;
	}
	
}

int getSubIndexInT2Matrix(int a, int b)
{
	// case dX2 dY2 or dT2
	if(a == b)
		return a ;
	// case dXdY or dYdX
	if(a+b == 1)
		return 3 ;
	// case dYdT or dTdY
	if(a+b == 3)
		return 4 ;
	// case dTdX or dXdT
	if(a+b == 2)
		return 5 ;
	// should not happen
	return 0 ;
}

int getSubIndexInT3Matrix(int a, int b, int c)
{
	// case dX3 or dY3 or dZ3
	if((a == b) && (a == c))
		return a ;
	if((a == b) || (a ==c) || (b == c))
	{
		int k = 3*a+c ;
		if( a == c)
			k = 3*a+b ;
		if(b == c)
			k = 3*b+a ;
		// case dX2dY
		if(k == 1)
			return 3 ;
		// case dY2dX
		if(k == 3)
			return 4 ;
		// case dY2dT
		if(k == 5)
			return 5 ;
		// case dT2dY
		if(k == 7)
			return 6 ;
		// case dT2dX
		if(k == 6)
			return 7 ;
		// case dX2dT
		if(k == 2)
			return 8 ;
	}
	// case dXdYdT
	return 9 ;
}

double coeffDeriv2(int X, int Y)
{
	if(X != Y)
		return 0.5 ;
	return 1. ;
}

double coeffDeriv3(int X, int Y, int T)
{
	if(X == Y && X != T)
		return 0.333333333333333 ;
	if(X == T && X != Y)
		return 0.333333333333333 ;		
	if(Y == T && X != Y)
		return 0.333333333333333 ;
	if(X == Y && Y == T)
		return 0.166666666666666 ;
	return 1. ;
}


	
void TriElement::getSecondJacobianMatrix(const Point &p, Matrix & t1, Matrix & t2)  
{
	#warning second jacobian matrices not implemented for spatial-only elements

	if(order >= CONSTANT_TIME_LINEAR)
	{
		Matrix Jinv(3,3) ;
		this->getInverseJacobianMatrix(p, Jinv) ;
		
		t2.resize(6,6) ;
		t1.resize(6,3) ;
		
		Tensor tn1(1,3) ;
		Tensor tn2(2,3) ;
		Tensor tj1(2,3) ;
		Tensor tjinv(Jinv, true) ;
		Tensor tj2(3,3) ;
		
		std::vector<Variable> var ;
		var.push_back(XI) ;
		var.push_back(ETA) ;
		var.push_back(TIME_VARIABLE) ;

		VirtualMachine vm ;
		for(size_t i = 0 ; i < getBoundingPoints().size() ; i++)
		{
			for(int a = 0 ; a < 3 ; a++)
			{
				tn1(a) = vm.deval( getShapeFunction(i), var[a], p ) ;
				// tensor of point ;
				Tensor tp(getBoundingPoint(i), var) ;
				for(int b = 0 ; b < 3 ; b++)
				{
					tn2(a,b) = vm.ddeval( getShapeFunction(i), var[a], var[b], p, 10*default_derivation_delta) ;
					tj1(a,b) += tn1(a) * tp(b) ;
					for(int c = 0 ; c < 3 ; c++)
					{
						tj2(a,b,c) += tn2(a,b)*tp(c)   ;
					}
				}
			}
		}
		
		Tensor dXXdxx(4,3) ;
		for(int x = 0 ; x < 3 ; x++)
		{
			for(int y = 0 ; y < 3 ; y++)
			{
				for(int X = 0 ; X < 3 ; X++)
				{
					for(int Y = 0 ; Y < 3 ; Y++)
					{
						dXXdxx(X,Y,x,y) = tjinv(X,x)*tjinv(Y,y) ;
// 						for(int z = 0 ; z < 3 ; z++)
// 						{
// 							dxxdX(x,y,X) += tj1(y,Y)*tjinv(Y,z)*tj2(x,z,X) ;
// 						}
					}
				}
			}
		}
		dXXdxx.threshold( POINT_TOLERANCE_3D ) ;
		tj2.threshold( POINT_TOLERANCE_3D ) ;
		
		Matrix c1(6,3) ;
		for(int X = 0 ; X < 3 ; X++)
		{
			for(int Y = 0 ; Y < 3 ; Y++)
			{
				int i = getSubIndexInT2Matrix(X,Y) ;
				for(int x = 0 ; x < 3 ; x++)
				{
					int j = x ;
					c1[i][j] = tj2(X,Y,x) ;
					for(int y = 0 ; y < 3 ; y++)
					{
						j = getSubIndexInT2Matrix(x,y) ;
						t2[ j][ i] = dXXdxx(X,Y,x,y) ;
					}
				}
			}
		}
		
		
		c1 *= -1. ;
		
		t1 = t2 * (c1 * Jinv) ;
	}
}

void TriElement::getThirdJacobianMatrix(const Point &p, Matrix & t1, Matrix & t2, Matrix & t3)  
{
	#warning third jacobian matrices not implemented for spatial-only elements

	if(order >= CONSTANT_TIME_LINEAR)
	{
		Matrix Jinv ;
		this->getInverseJacobianMatrix(p, Jinv) ;
		
		t3.resize(10,10) ;
		t2.resize(10,6) ;
		t1.resize(10,3) ;
		
		Matrix j3(10,3) ;
		Matrix j2(6,3) ;
		Matrix j1(3,3) ;
		
		
		j1 = Jinv ;
		invert3x3Matrix(j1) ;
		
		std::vector<Variable> var ;
		var.push_back(XI) ;
		var.push_back(ETA) ;
		var.push_back(TIME_VARIABLE) ;

		// tensors of derivative
		Tensor tn1(1,3) ;
		Tensor tn2(2,3) ;
		Tensor tn3(3,3) ;

		// tensor of jacobians
		Tensor tj1(2,3) ;
		Tensor tjinv(Jinv, false) ;
		Tensor tj2(3,3) ;
		Tensor tj3(4,3) ;
		
		
		
		VirtualMachine vm ;
		for(size_t i = 0 ; i < getBoundingPoints().size() ; i++)
		{
			for(int a = 0 ; a < 3 ; a++)
			{
				tn1(a) = vm.deval( getShapeFunction(i), var[a], p ) ;
				// tensor of point ;
				Tensor tp(getBoundingPoint(i), var) ;
				for(int b = 0 ; b < 3 ; b++)
				{
					tn2(a,b) = vm.ddeval( getShapeFunction(i), var[a], var[b], p, 100*default_derivation_delta) ;
					tj1(a,b) += tn1(a) * tp(b) ;
					for(int c = 0 ; c < 3 ; c++)
					{
						tn3(a,b,c) = vm.dddeval( getShapeFunction(i), var[a], var[b], var[c],p.getX(),p.getY(),p.getZ(),p.getT()) ;
						tn3.threshold(POINT_TOLERANCE_3D) ;
						tj2(a,b,c) += tn2(a,b)*tp(c)  ;
						for(int d = 0 ; d < 3 ; d++)
							tj3(a,b,c,d) += tn3(a,b,c)*tp(d) ;
					}
				}
			}
		}
		
		tjinv.threshold( POINT_TOLERANCE_3D ) ;
		tj1.threshold(POINT_TOLERANCE_3D) ;
		tj2.threshold(POINT_TOLERANCE_3D) ;
		tj3.threshold(POINT_TOLERANCE_3D) ;
		
		Matrix T1 ;
		Matrix T2 ;
		this->getSecondJacobianMatrix(p, T1, T2) ;		
		Tensor dXXdx(3,3) ;
		Tensor dXXdxx(4,3) ;

		for(int X = 0 ; X < 3 ; X++)
		{
			for(int Y = 0 ; Y < 3 ; Y++)
			{
				for(int x = 0 ; x < 3 ; x++)
				{
					dXXdx(X,Y,x) = T1[ getSubIndexInT2Matrix(X,Y) ][x] ;
					for(int y = 0 ; y < 3 ; y++)
						dXXdxx(X,Y,x,y) = T2[ getSubIndexInT2Matrix(X,Y) ][ getSubIndexInT2Matrix(x,y) ] ;
				}
			}
		}
		
		Tensor dxxdX(3,3) ;
		Tensor dxxdXX(4,3) ;
		for(int x = 0 ; x < 3 ; x++)
		{
			for(int y = 0 ; y < 3 ; y++)
			{
				for(int X = 0 ; X < 3 ; X++)
				{
					for(int Y = 0 ; Y < 3 ; Y++)
					{
						dxxdXX(x,y,X,Y) = tj1(x,X)*tj1(y,Y) ;
						for(int u = 0 ; u < 3 ; u++)
							dxxdX(x,y,X) += tj1(y,Y)*tjinv(Y,u)*tj2(x,u,X) ;
					}
				}
			}
		}
		
		Tensor dxxxdX(4,3) ;
		Tensor dxxxdXX(5,3) ;
		Tensor dxxxdXXX(6,3) ;
		Tensor dXXXdxxx(6,3) ;
		for(int x = 0 ; x < 3 ; x++)
		{
			for(int y = 0 ; y < 3 ; y++)
			{
				for(int t = 0 ; t < 3 ; t++)
				{
					for(int X = 0 ; X < 3 ; X++)
					{
						for(int Y = 0 ; Y < 3 ; Y++)
						{
							for(int T = 0 ; T < 3 ; T++)
							{
								dxxxdXXX(x,y,t,X,Y,T) = tj1(x,X)*tj1(y,Y)*tj1(t,T) ;
								dXXXdxxx(X,Y,T,x,y,t) = tjinv(X,x)*tjinv(Y,y)*tjinv(T,t) ;
								for(int u = 0 ; u < 3 ; u++)
								{
									dxxxdXX(x,y,t,X,Y) += tj1(t,T)*tjinv(T,u)*tj2(y,u,Y)*tj1(x,X) ; 
									dxxxdXX(x,y,t,X,Y) += tj1(t,T)*tjinv(T,u)*tj2(x,u,X)*tj1(y,Y) ; 
									dxxxdXX(x,y,t,X,T) += tj1(t,T)*tjinv(Y,u)*tj2(x,u,X)*tj1(y,Y) ;
									dxxxdX(x,y,t,X) += tj1(t,T)*tj1(y,Y)*dXXdx(Y,T,u)*tj2(x,u,X) ;
									for(int v = 0 ; v < 3 ; v++)
									{
dxxxdX(x,y,t,X) += tj1(t,T)*tjinv(Y,u)*tj2(x,u,X)*tjinv(T,v)*tj2(y,v,Y) ;
dxxxdX(x,y,t,X) += tj1(t,T)*tj1(y,Y)*dXXdxx(Y,T,u,v)*tj3(x,u,v,X) ;
									}
								}
							}
						}
					}
				}
			}
		}

		Tensor dxxxdx(4,3) ;
		Tensor dxxxdxx(5,3) ;
		for(int x = 0 ; x < 3 ; x++)
		{
			for(int y = 0 ; y < 3 ; y++)
			{
				for(int t = 0 ; t < 3 ; t++)
				{
					for(int X = 0 ; X < 3 ; X++)
					{
						for(int u = 0 ; u < 3 ; u++)
						{
							dxxxdx(x,y,t,u) += dxxxdX(x,y,t,X) * tjinv(X,u) ;
							for(int Y = 0 ; Y < 3 ; Y++)
							{
								dxxxdx(x,y,t,u) += dxxxdXX(x,y,t,X,Y) * dXXdx(X,Y,u) ;
								for(int v = 0 ; v < 3 ; v++)
								{
									dxxxdxx(x,y,t,u,v) += dxxxdXX(x,y,t,X,Y) * dXXdxx(X,Y,u,v) ;
								}
							}
						}
					}
				}
			}
		}
		
		Tensor dXXXdx(4,3) ;
		Tensor dXXXdxx(5,3) ;
		for(int X = 0 ; X < 3 ; X++)
		{
			for(int Y = 0 ; Y < 3 ; Y++)
			{
				for(int T = 0 ; T < 3 ; T++)
				{
					for(int u = 0 ; u < 3 ; u++)
					{
						for(int x = 0 ; x < 3 ; x++)
						{
							for(int y = 0 ; y < 3 ; y++)
							{
								for(int t = 0 ; t < 3 ; t++)
								{
									dXXXdx(X,Y,T,u) += dXXXdxxx(X,Y,T,x,y,t) * dxxxdx(x,y,t,u) ;
								}
							}
						}
						
						for(int v = 0 ; v < 3 ; v++)
						{
							for(int x = 0 ; x < 3 ; x++)
							{
								for(int y = 0 ; y < 3 ; y++)
								{
									for(int t = 0 ; t < 3 ; t++)
									{
											dXXXdxx(X,Y,T,u,v) += dXXXdxxx(X,Y,T,x,y,t) * dxxxdxx(x,y,t,u,v) ;
									}
								}
							}
						}
					}
				}
			}
		}
		
		dXXXdx.threshold( POINT_TOLERANCE_3D ) ;
		dXXXdxx.threshold( POINT_TOLERANCE_3D ) ;
		dXXXdxxx.threshold( POINT_TOLERANCE_3D ) ;
		
		Matrix m1 = dXXXdx.toMatrix(3,1) ;
		Matrix m2 = dXXXdxx.toMatrix(3,2) ;
		Matrix m3inv = dXXXdxxx.toMatrix(3,3) ;
		
// 		m1.print() ;
// 		m2.print() ;
// 		m3inv.print() ;
//		exit(0) ;
		
		t1 = Matrix(10,3) ;
		t2 = Matrix(10,6) ;
		t3 = Matrix(10,10) ;
		
		t1[6][0] = dXXXdx(1,2,2,0) ;
		t1[6][1] = dXXXdx(1,2,2,1) ; 
		t1[6][2] = dXXXdx(1,2,2,2) ;
		t1[7][0] = dXXXdx(0,2,2,0) ;
		t1[7][1] = dXXXdx(0,2,2,1) ; 
		t1[7][2] = dXXXdx(0,2,2,2) ;
		
		t2[6][0] = dXXXdxx(1,2,2,0,0) ;
		t2[6][1] = dXXXdxx(1,2,2,1,1) ;
		t2[6][2] = dXXXdxx(1,2,2,2,2) ;
		t2[6][3] = (dXXXdxx(1,2,2,0,1) + dXXXdxx(1,2,2,1,0)) ;
		t2[6][4] = (dXXXdxx(1,2,2,1,2) + dXXXdxx(1,2,2,2,1)) ;
		t2[6][5] = (dXXXdxx(1,2,2,2,0) + dXXXdxx(1,2,2,0,2)) ;
		t2[7][0] = dXXXdxx(0,2,2,0,0) ;
		t2[7][1] = dXXXdxx(0,2,2,1,1) ;
		t2[7][2] = dXXXdxx(0,2,2,2,2) ;
		t2[7][3] = (dXXXdxx(0,2,2,0,1) + dXXXdxx(0,2,2,1,0)) ;
		t2[7][4] = (dXXXdxx(0,2,2,1,2) + dXXXdxx(0,2,2,2,1)) ;
		t2[7][5] = (dXXXdxx(0,2,2,2,0) + dXXXdxx(0,2,2,0,2)) ;

// 		t3[6][0] = dXXXdxxx(1,2,2,0,0,0) ;
// 		t3[6][1] = dXXXdxxx(1,2,2,1,1,1) ;
// 		t3[6][2] = dXXXdxxx(1,2,2,2,2,2) ;
// 		t3[6][3] = dXXXdxxx(1,2,2,0,0,1) + dXXXdxxx(1,2,2,0,1,0) + dXXXdxxx(1,2,2,1,0,0) ;
// 		t3[6][4] = dXXXdxxx(1,2,2,1,1,0) + dXXXdxxx(1,2,2,0,1,1) + dXXXdxxx(1,2,2,1,0,1)  ;
// 		t3[6][5] = dXXXdxxx(1,2,2,1,1,2) + dXXXdxxx(1,2,2,1,2,1) + dXXXdxxx(1,2,2,2,1,1)  ;
// 		t3[6][6] = dXXXdxxx(1,2,2,2,2,1) + dXXXdxxx(1,2,2,2,1,2) + dXXXdxxx(1,2,2,1,2,2)  ;
// 		t3[6][7] = dXXXdxxx(1,2,2,2,2,0) + dXXXdxxx(1,2,2,0,2,2) + dXXXdxxx(1,2,2,2,0,2)  ;
// 		t3[6][8] = dXXXdxxx(1,2,2,0,0,2) + dXXXdxxx(1,2,2,0,2,0) + dXXXdxxx(1,2,2,2,0,0)  ;
// 		t3[6][9] = dXXXdxxx(1,2,2,0,1,2) + dXXXdxxx(1,2,2,0,2,1) + dXXXdxxx(1,2,2,1,2,0) ;
// 		t3[6][9] += dXXXdxxx(1,2,2,1,0,2) + dXXXdxxx(1,2,2,2,0,1) + dXXXdxxx(1,2,2,2,1,0) ;
// 		t3[7][0] = dXXXdxxx(0,2,2,0,0,0) ;
// 		t3[7][1] = dXXXdxxx(0,2,2,1,1,1) ;
// 		t3[7][2] = dXXXdxxx(0,2,2,2,2,2) ;
// 		t3[7][3] = dXXXdxxx(0,2,2,0,0,1) + dXXXdxxx(0,2,2,0,1,0) + dXXXdxxx(0,2,2,1,0,0) ;
// 		t3[7][4] = dXXXdxxx(0,2,2,1,1,0) + dXXXdxxx(0,2,2,0,1,1) + dXXXdxxx(0,2,2,1,0,1);
// 		t3[7][5] = dXXXdxxx(0,2,2,1,1,2) + dXXXdxxx(0,2,2,1,2,1) + dXXXdxxx(0,2,2,2,1,1) ;
// 		t3[7][6] = dXXXdxxx(0,2,2,2,2,1) + dXXXdxxx(0,2,2,2,1,2) + dXXXdxxx(0,2,2,1,2,2) ;
// 		t3[7][7] = dXXXdxxx(0,2,2,2,2,0) + dXXXdxxx(0,2,2,0,2,2) + dXXXdxxx(0,2,2,2,0,2)  ;
// 		t3[7][8] = dXXXdxxx(0,2,2,0,0,2) + dXXXdxxx(0,2,2,0,2,0) + dXXXdxxx(0,2,2,2,0,0);
// 		t3[7][9] = dXXXdxxx(0,2,2,0,1,2) + dXXXdxxx(0,2,2,0,2,1) + dXXXdxxx(0,2,2,1,2,0) ;
// 		t3[7][9] += dXXXdxxx(0,2,2,1,0,2) + dXXXdxxx(0,2,2,2,0,1) + dXXXdxxx(0,2,2,2,1,0) ;
		
// 		t1 *= 0. ;
// 		t2 *= 0. ;
// 		for(int i = 0 ; i < 10 ; i++)
// 		{
// 			if(i < 3)
// 			{
// 				t1[6][i] = -m1[8][i] ;
// 				t1[7][i] = -m1[17][i] ;
// // 				//dXXX
// // 				t1[0][i] = d1[0][i] ;
// // 				t1[1][i] = d1[13][i] ;
// // 				t1[2][i] = d1[26][i] ;
// // 				//dXXY
// // 				t1[3][i] = (d1[1][i] + d1[3][i] + d1[9][i])/3. ;
// // 				t1[4][i] = (d1[4][i] + d1[10][i] + d1[12][i])/3. ;
// // 				t1[5][i] = (d1[14][i] + d1[16][i] + d1[22][i])/3. ;
// // 				t1[6][i] = (d1[17][i] + d1[23][i] + d1[25][i])/3. ;
// // 				t1[7][i] = (d1[8][i] + d1[20][i] + d1[24][i])/3. ;
// // 				t1[8][i] = (d1[2][i] + d1[6][i] + d1[18][i])/3. ;
// // 				//dXYT
// // 				t1[9][i] = (d1[5][i] + d1[7][i] + d1[11][i] + d1[15][i] + d1[19][i] + d1[21][i])/6. ;
// 			}
// 			if(i < 6)
// 			{
// 				t2[6][i] = -m2[8][i] ;
// 				t2[7][i] = -m2[17][i] ;
// // 				t2[0][i] = d2[0][i] ;
// // 				t2[1][i] = d2[13][i] ;
// // 				t2[2][i] = d2[26][i] ;
// // 
// // 				t2[3][i] = (d2[1][i] + d2[3][i] + d2[9][i])/3. ;
// // 				t2[4][i] = (d2[4][i] + d2[10][i] + d2[12][i])/3. ;
// // 				t2[5][i] = (d2[14][i] + d2[16][i] + d2[22][i])/3. ;
// // 				t2[6][i] = (d2[17][i] + d2[23][i] + d2[25][i])/3. ;
// // 				t2[7][i] = (d2[8][i] + d2[20][i] + d2[24][i])/3. ;
// // 				t2[8][i] = (d2[2][i] + d2[6][i] + d2[18][i])/3. ;
// // 			  
// // 				t2[9][i] = (d2[5][i] + d2[7][i] + d2[11][i] + d2[15][i] + d2[19][i] + d2[21][i])/6. ;
// 			}
// 			t3[6][i] = m3inv[8][i] ;
// 			t3[7][i] = m3inv[17][i] ;
// // 			t3[0][i] = m3inv[0][i] ;
// // 			t3[1][i] = m3inv[13][i] ;
// // 			t3[2][i] = m3inv[26][i] ;
// // 
// // 			t3[3][i] = (m3inv[1][i] + m3inv[3][i] + m3inv[9][i])/1. ;
// // 			t3[4][i] = (m3inv[4][i] + m3inv[10][i] + m3inv[12][i])/1. ;
// // 			t3[5][i] = (m3inv[14][i] + m3inv[16][i] + m3inv[22][i])/1. ;
// // 			t3[6][i] = (m3inv[17][i] + m3inv[23][i] + m3inv[25][i])/1. ;
// // 			t3[7][i] = (m3inv[8][i] + m3inv[20][i] + m3inv[24][i])/1. ;
// // 			t3[8][i] = (m3inv[2][i] + m3inv[6][i] + m3inv[18][i])/1. ;
// // 
// // 			t3[9][i] = (m3inv[5][i] + m3inv[7][i] + m3inv[11][i] + m3inv[15][i] + m3inv[19][i] + m3inv[21][i])/6. ;
// // 		  
// 		}
		
		
		
// 		t1 = (Matrix) (t3 * d1) ;//( (Matrix) ( d1*Jinv ) + (Matrix) ( d2*T1 ) ) );
// 		t2 = (Matrix) (t3 * d2) ;//( (Matrix) ( d2*T2) )  );

//  		t1.print() ;
// 		std::cout << std::endl ;
// 		t2.print() ;
// 		std::cout << std::endl ;
// 		t3.print() ;
// 		std::cout << std::endl ;
//		exit(0) ;
		

		
		
/*			double dxidxi = vm.ddeval( getShapeFunction(i), XI, XI, p, 100*default_derivation_delta ) ;
			double dxideta = vm.ddeval( getShapeFunction(i), XI, ETA, p, 100*default_derivation_delta ) ;
			double dxidtau = vm.ddeval( getShapeFunction(i), XI, TIME_VARIABLE, p, 100*default_derivation_delta ) ;
			double detadxi = vm.ddeval( getShapeFunction(i), ETA, XI, p, 100*default_derivation_delta ) ;
			double detadeta = vm.ddeval( getShapeFunction(i), ETA, ETA, p, 100*default_derivation_delta ) ;
			double detadtau = vm.ddeval( getShapeFunction(i), ETA, TIME_VARIABLE, p, 100*default_derivation_delta ) ;
			double dtaudxi = vm.ddeval( getShapeFunction(i), TIME_VARIABLE, XI, p, 100*default_derivation_delta ) ;
			double dtaudeta = vm.ddeval( getShapeFunction(i), TIME_VARIABLE, ETA, p, 100*default_derivation_delta ) ;
			double dtaudtau = vm.ddeval( getShapeFunction(i), TIME_VARIABLE, TIME_VARIABLE, p, 100*default_derivation_delta ) ;
			
			double dxidxidxi = vm.dddeval( getShapeFunction(i), XI, XI, XI, p.getX(),p.getY(),p.getZ(),p.getT()) ;
			double detadetadeta = vm.dddeval( getShapeFunction(i), ETA, ETA, ETA, p.getX(),p.getY(),p.getZ(),p.getT()) ;
			double dtaudtaudtau = vm.dddeval( getShapeFunction(i), TIME_VARIABLE, TIME_VARIABLE, TIME_VARIABLE, p.getX(),p.getY(),p.getZ(),p.getT()) ;
			
			double dxidxideta = vm.dddeval( getShapeFunction(i), XI, XI, ETA, p.getX(),p.getY(),p.getZ(),p.getT()) ;
			double dxidetadxi = vm.dddeval( getShapeFunction(i), XI, ETA, XI, p.getX(),p.getY(),p.getZ(),p.getT()) ;
			double detadxidxi = vm.dddeval( getShapeFunction(i), ETA, XI, XI, p.getX(),p.getY(),p.getZ(),p.getT()) ;
			double dxidxidtau = vm.dddeval( getShapeFunction(i), XI, XI, TIME_VARIABLE, p.getX(),p.getY(),p.getZ(),p.getT()) ;
			double dxidtaudxi = vm.dddeval( getShapeFunction(i), XI, TIME_VARIABLE, XI, p.getX(),p.getY(),p.getZ(),p.getT()) ;
			double dtaudxidxi = vm.dddeval( getShapeFunction(i), TIME_VARIABLE, XI, XI, p.getX(),p.getY(),p.getZ(),p.getT()) ;

			double detadetadxi = vm.dddeval( getShapeFunction(i), ETA, ETA, XI, p.getX(),p.getY(),p.getZ(),p.getT()) ;
			double detadxideta = vm.dddeval( getShapeFunction(i), ETA, XI, ETA, p.getX(),p.getY(),p.getZ(),p.getT()) ;
			double dxidetadeta = vm.dddeval( getShapeFunction(i), XI, ETA, ETA, p.getX(),p.getY(),p.getZ(),p.getT()) ;
			double detadetadtau = vm.dddeval( getShapeFunction(i), ETA, ETA, TIME_VARIABLE, p.getX(),p.getY(),p.getZ(),p.getT()) ;
			double detadtaudeta = vm.dddeval( getShapeFunction(i), ETA, TIME_VARIABLE, ETA, p.getX(),p.getY(),p.getZ(),p.getT()) ;
			double dtaudetadeta = vm.dddeval( getShapeFunction(i), TIME_VARIABLE, ETA, ETA, p.getX(),p.getY(),p.getZ(),p.getT()) ;

			double dtaudtaudxi = vm.dddeval( getShapeFunction(i), TIME_VARIABLE, TIME_VARIABLE, XI, p.getX(),p.getY(),p.getZ(),p.getT()) ;
			double dtaudxidtau = vm.dddeval( getShapeFunction(i), TIME_VARIABLE, XI, TIME_VARIABLE, p.getX(),p.getY(),p.getZ(),p.getT()) ;
			double dxidtaudtau = vm.dddeval( getShapeFunction(i), XI, TIME_VARIABLE, TIME_VARIABLE, p.getX(),p.getY(),p.getZ(),p.getT()) ;
			double dtaudtaudeta = vm.dddeval( getShapeFunction(i), TIME_VARIABLE, TIME_VARIABLE, ETA, p.getX(),p.getY(),p.getZ(),p.getT()) ;
			double dtaudetadtau = vm.dddeval( getShapeFunction(i), TIME_VARIABLE, ETA, TIME_VARIABLE, p.getX(),p.getY(),p.getZ(),p.getT()) ;
			double detadtaudtau = vm.dddeval( getShapeFunction(i), ETA, TIME_VARIABLE, TIME_VARIABLE, p.getX(),p.getY(),p.getZ(),p.getT()) ;
			
			double dxidetadtau = vm.dddeval( getShapeFunction(i), XI, ETA, TIME_VARIABLE, p.getX(),p.getY(),p.getZ(),p.getT()) ;
			double dxidtaudeta = vm.dddeval( getShapeFunction(i), XI, TIME_VARIABLE, ETA, p.getX(),p.getY(),p.getZ(),p.getT()) ;
			double detadxidtau = vm.dddeval( getShapeFunction(i), ETA, XI, TIME_VARIABLE, p.getX(),p.getY(),p.getZ(),p.getT()) ;
			double detadtaudxi = vm.dddeval( getShapeFunction(i), ETA, TIME_VARIABLE, XI, p.getX(),p.getY(),p.getZ(),p.getT()) ;
			double dtaudxideta = vm.dddeval( getShapeFunction(i), TIME_VARIABLE, XI, ETA, p.getX(),p.getY(),p.getZ(),p.getT()) ;
			double dtaudetadxi = vm.dddeval( getShapeFunction(i), TIME_VARIABLE, ETA, XI, p.getX(),p.getY(),p.getZ(),p.getT()) ;

			
			j2[0][0] += dxidxi * getBoundingPoint(i).getX() ;
			j2[0][1] += dxidxi * getBoundingPoint(i).getY() ;
			j2[0][2] += dxidxi * getBoundingPoint(i).getT() ;
			
			j2[1][0] += detadeta * getBoundingPoint(i).getX() ;
			j2[1][1] += detadeta * getBoundingPoint(i).getY() ;
			j2[1][2] += detadeta * getBoundingPoint(i).getT() ;
			
			j2[2][0] += dtaudtau * getBoundingPoint(i).getX() ;
			j2[2][1] += dtaudtau * getBoundingPoint(i).getY() ;
			j2[2][2] += dtaudtau * getBoundingPoint(i).getT() ;
			
			j2[3][0] += 0.5 * (detadxi + dxideta) * getBoundingPoint(i).getX() ;
			j2[3][1] += 0.5 * (detadxi + dxideta) * getBoundingPoint(i).getY() ;
			j2[3][2] += 0.5 * (detadxi + dxideta) * getBoundingPoint(i).getT() ;

			j2[4][0] += 0.5 * (dtaudeta + detadtau) * getBoundingPoint(i).getX() ;
			j2[4][1] += 0.5 * (dtaudeta + detadtau) * getBoundingPoint(i).getY() ;
			j2[4][2] += 0.5 * (dtaudeta + detadtau) * getBoundingPoint(i).getT() ;
			
			j2[5][0] += 0.5 * (dxidtau + dtaudxi) * getBoundingPoint(i).getX() ;
			j2[5][1] += 0.5 * (dxidtau + dtaudxi) * getBoundingPoint(i).getY() ;
			j2[5][2] += 0.5 * (dxidtau + dtaudxi) * getBoundingPoint(i).getT() ;
			
			
			
			
			
			
			j3[0][0] += dxidxidxi * getBoundingPoint(i).getX() ;
			j3[0][1] += dxidxidxi * getBoundingPoint(i).getY() ;
			j3[0][2] += dxidxidxi * getBoundingPoint(i).getT() ;
			
			j3[1][0] += detadetadeta * getBoundingPoint(i).getX() ;
			j3[1][1] += detadetadeta * getBoundingPoint(i).getY() ;
			j3[1][2] += detadetadeta * getBoundingPoint(i).getT() ;

			j3[2][0] += dtaudtaudtau * getBoundingPoint(i).getX() ;
			j3[2][1] += dtaudtaudtau * getBoundingPoint(i).getY() ;
			j3[2][2] += dtaudtaudtau * getBoundingPoint(i).getT() ;
			
			j3[3][0] += (dxidxideta + dxidetadxi + detadxidxi) /3. * getBoundingPoint(i).getX() ;
			j3[3][1] += (dxidxideta + dxidetadxi + detadxidxi) /3. * getBoundingPoint(i).getY() ;
			j3[3][2] += (dxidxideta + dxidetadxi + detadxidxi) /3. * getBoundingPoint(i).getT() ;
			
			j3[4][0] += (detadetadxi + detadxideta + dxidetadeta) /3. * getBoundingPoint(i).getX() ;
			j3[4][1] += (detadetadxi + detadxideta + dxidetadeta)  /3. * getBoundingPoint(i).getY() ;
			j3[4][2] += (detadetadxi + detadxideta + dxidetadeta)  /3. * getBoundingPoint(i).getT() ;
			
			j3[5][0] += (detadetadtau + detadtaudeta + dtaudetadeta) /3. * getBoundingPoint(i).getX() ;
			j3[5][1] += (detadetadtau + detadtaudeta + dtaudetadeta) /3. * getBoundingPoint(i).getY() ;
			j3[5][2] += (detadetadtau + detadtaudeta + dtaudetadeta) /3. * getBoundingPoint(i).getT() ;
			
			j3[6][0] += (dtaudtaudeta + dtaudetadtau + detadtaudtau) /3. * getBoundingPoint(i).getX() ;
			j3[6][1] += (dtaudtaudeta + dtaudetadtau + detadtaudtau)  /3. * getBoundingPoint(i).getY() ;
			j3[6][2] += (dtaudtaudeta + dtaudetadtau + detadtaudtau)  /3. * getBoundingPoint(i).getT() ;
			
			j3[7][0] += (dtaudtaudxi + dtaudxidtau + dxidtaudtau) /3. * getBoundingPoint(i).getX() ;
			j3[7][1] += (dtaudtaudxi + dtaudxidtau + dxidtaudtau) /3. * getBoundingPoint(i).getY() ;
			j3[7][2] += (dtaudtaudxi + dtaudxidtau + dxidtaudtau) /3. * getBoundingPoint(i).getT() ;
			
			j3[8][0] += (dxidxidtau + dxidtaudxi + dtaudxidxi) /3. * getBoundingPoint(i).getX() ;
			j3[8][1] += (dxidxidtau + dxidtaudxi + dtaudxidxi)  /3. * getBoundingPoint(i).getY() ;
			j3[8][2] += (dxidxidtau + dxidtaudxi + dtaudxidxi)  /3. * getBoundingPoint(i).getT() ;
			
			j3[9][0] += (dxidetadtau + dxidtaudeta + detadxidtau + detadtaudxi + dtaudxideta + dtaudetadxi ) /6. * getBoundingPoint(i).getX() ;
			j3[9][1] += (dxidetadtau + dxidtaudeta + detadxidtau + detadtaudxi + dtaudxideta + dtaudetadxi ) /6. * getBoundingPoint(i).getY() ;
			j3[9][2] += (dxidetadtau + dxidtaudeta + detadxidtau + detadtaudxi + dtaudxideta + dtaudetadxi ) /6. * getBoundingPoint(i).getT() ;*/
			
			
		
/*		for(size_t i = 0 ; i < j2.numRows() ; i++)
		{
			for(size_t j = 0 ; j < j2.numCols() ; j++)
			{
				if( std::abs(j2[i][j]) < POINT_TOLERANCE_3D)
					j2[i][j] = 0 ;
			}
		}

		for(size_t i = 0 ; i < j3.numRows() ; i++)
		{
			for(size_t j = 0 ; j < j3.numCols() ; j++)
			{
				if( std::abs(j3[i][j]) < POINT_TOLERANCE_3D)
					j3[i][j] = 0 ;
			}
		}*/

/*		Matrix d2(10,6) ;
		for(size_t i = 0 ; i < 3 ; i++)
		{
			for(size_t j = 0 ; j < 3 ; j++)
			{
				d2[i][j] = 3. * j2[i][j]*j1[i][j] ;		
				d2[i][j+3] = 3. * (j2[i][j]*j1[i][(j+1)%3]+j2[i][(j+1)%3]*j1[i][j]) ;
				
				d2[3+i*2][j] = 3. * (j2[i][j]*j1[(i+1)%3][j] + 2*(j2[3+i][j]*j1[i][j])) ;
				d2[3+i*2][j+3] = 3. * (j2[i][j]*j1[(i+1)%3][(j+1)%3] + 2*(j2[3+i][j]*j1[i][(j+1)%3]) + j2[i][(j+1)%3]*j1[(i+1)%3][j] + 2*(j2[3+i][(j+1)%3]*j1[i][j])) ;
				
				d2[3+i*2+1][j] = 3. * (j2[(i+1)%3][j]*j1[i][j] + 2*(j2[3+i][j]*j1[(i+1)%3][j])) ;
				d2[3+i*2+1][j+3] = 3. * ( j2[(i+1)%3][j]*j1[i][(j+1)%3] + 2*(j2[3+i][j]*j1[(i+1)%3][(j+1)%3]) + j2[(i+1)%3][(j+1)%3]*j1[i][j] + 2*(j2[3+i][j]*j1[(i+1)%3][j])) ;
				
				d2[9][j] += 2*( 2*j2[3+i][j]*j1[(i+2)%3][j] + j2[(i+2)%3+3][j]*j1[(i+1)%3][j] ) ;
				d2[9][j+3] += 2*( 2*j2[3+i][(j+1)%3]*j1[(i+2)%3][j] + 2*j2[3+i][j]*j1[(i+2)%3][(j+1)%3] + j2[(i+2)%3+3][(j+1)%3]*j1[(i+1)%3][j] + j2[(i+2)%3+3][j]*j1[(i+1)%3][(j+1)%3] ) ;
			}
		}

// 		j1.print() ;
// 		std::cout << std::endl ;
// 		j2.print() ;
// 		std::cout << std::endl ;
// 		d2.print() ;
// 		std::cout << std::endl ;

		
		Matrix c1 ;
		Matrix c2 ;		
		this->getSecondJacobianMatrix( p, c1, c2) ;
		
		Matrix d3(10,10) ;
		for(size_t i = 0 ; i < 3 ; i++)
		{
			for(size_t j = 0 ; j < 3 ; j++)
			{
				d3[i][j] = Jinv[i][j]*c2[i][j] ;
				d3[i][3+j*2] = Jinv[i][(j+1)%3]*c2[i][j] + Jinv[i][j]*c2[i][j+3] ;
				d3[i][3+j*2+1] = Jinv[i][j]*c2[i][(j+1)%3] + Jinv[i][(j+1)%3]*c2[i][j+3] ;
				d3[i][9] += Jinv[i][(j+2)%3]*c2[i][j+3] ;
				
				d3[3+i*2][j] = Jinv[(i+1)%3][j]*c2[i][j] + Jinv[i][j]*c2[i+3][j] ;
				d3[3+i*2][3+j*2] = Jinv[(i+1)%3][(j+1)%3]*c2[i][j] + Jinv[(i+1)%3][j]*c2[i][j+3] + Jinv[i][(j+1)%3]*c2[i+3][j] + Jinv[i][j] * c2[i+3][j+3] ;
				d3[3+i*2][3+j*2+1] = Jinv[(i+1)%3][j]*c2[i][(j+1)%3] + Jinv[(i+1)%3][(j+1)%3]*c2[i][j+3] + Jinv[i][j]*c2[i+3][(j+1)%3] + Jinv[i][(j+1)%3] * c2[i+3][j+3] ;
				d3[3+i*2][9] += Jinv[(i+1)%3][(j+2)%3]*c2[i][j+3] + Jinv[i][(j+2)%3]*c2[i+3][j+3] ;
 			  
				d3[3+i*2+1][j] = Jinv[i][j]*c2[(i+1)%3][j] + Jinv[(i+1)%3][j]*c2[i+3][j] ;
				d3[3+i*2+1][3+j*2] = Jinv[i][(j+1)%3]*c2[(i+1)%3][j] + Jinv[i][j]*c2[(i+1)%3][j+3] + Jinv[(i+1)%3][(j+1)%3]*c2[i+3][j] + Jinv[(i+1)%3][j] * c2[i+3][j+3] ;
				d3[3+i*2+1][3+j*2+1] = Jinv[i][j]*c2[(i+1)%3][(j+1)%3] + Jinv[i][(j+1)%3]*c2[(i+1)%3][j+3] + Jinv[(i+1)%3][j]*c2[i+3][(j+1)%3] + Jinv[(i+1)%3][(j+1)%3] * c2[i+3][j+3] ;
				d3[3+i*2+1][9] += Jinv[i][(j+2)%3]*c2[(i+1)%3][j+3] + Jinv[(i+1)%3][(j+2)%3]*c2[i+3][j+3] ;
				
				d3[9][j] += Jinv[i][j]*c2[(i+1)%3+3][j] ;
				d3[9][3+j*2] += Jinv[i][(j+1)%3]*c2[(i+1)%3+3][j] + Jinv[i][j]*c2[(i+1)%3+3][j+3] ;
				d3[9][3+j*2+1] += Jinv[i][j]*c2[(i+1)%3+3][(j+1)%3] + Jinv[i][(j+1)%3]*c2[(i+1)%3+3][j+3] ;
				d3[9][9] += Jinv[(i+2)%3][j]*c2[i+3][3+(j+1)%3] ;
			}
		} 
		
		Matrix d1 = j3 ;
		
	*/	
	}
}

	
void TriElement::getInverseJacobianMatrix(const Point & p, Matrix & ret) 
{
	if(order < CONSTANT_TIME_LINEAR)
	{
		if(!isMoved() && !cachedJinv.empty())
		{
			if(ret.isNull())
				ret.resize(2,2) ;
			ret.array() = cachedJinv[0].array() ;
			return ;
		}
		
		if(ret.isNull())
			ret.resize(2,2) ;
		
		double xdxi = 0 ;//this->getdXTransform(XI,p) ;
		double ydxi = 0 ;//this->getdYTransform(XI,p) ;
		
		double xdeta = 0 ;//this->getdXTransform(ETA,p) ;
		double ydeta = 0 ;//this->getdYTransform(ETA,p) ;

		VirtualMachine vm ;
		TriElement * father = nullptr ;

		std::valarray<Function> * functions = shapefunc ;
		if(!shapefunc)
		{
			father = new TriElement(order) ;
			functions = father->shapefunc ;
		}
		
		for(size_t i = 0 ; i < getBoundingPoints().size() ; i++)
		{
			double dxi = vm.deval((*functions)[i], XI, p) ;
			double deta = vm.deval((*functions)[i], ETA, p) ;
			
			xdxi += dxi*getBoundingPoint(i).getX() ;
			ydxi += dxi*getBoundingPoint(i).getY() ;

			xdeta += deta*getBoundingPoint(i).getX() ;
			ydeta += deta*getBoundingPoint(i).getY() ;
		}
		ret[0][0] = xdxi ; ret[0][1] = ydxi ; 
		ret[1][0] = xdeta ; ret[1][1] = ydeta ;
		invert2x2Matrix(ret) ;
		if(cachedJinv.empty() && !isMoved())
			cachedJinv.push_back(ret) ;
		delete father ;
// 		ret.print() ;
// 		exit(0) ;
	}
	else
	{
/* 		if(!isMoved() && !cachedJinv.empty())
 		{
 			if(ret.isNull()  || ret.size() != 9)
 				ret.resize(3,3) ;
 			ret.array() = cachedJinv[0].array() ;
//			std::cout << ret.numCols() << "\t" << ret.numRows() << std::endl ;
//			ret.print() ;
 			return ;
 		}*/
		
		if(ret.isNull() || ret.size() != 9)
			ret.resize(3,3) ;

		double xdxi = 0 ;//this->getdXTransform(XI,p) ;
		double ydxi = 0 ;//this->getdYTransform(XI,p) ;
		double tdxi = 0 ;

		double xdeta = 0 ;//this->getdXTransform(ETA,p) ;
		double ydeta = 0 ;//this->getdYTransform(ETA,p) ;
		double tdeta = 0 ;

		double xdtau = 0 ;
		double ydtau = 0 ;
		double tdtau = 0 ;

		VirtualMachine vm ;
		
		Point dummy ;
		for(size_t i = 0 ; i < getBoundingPoints().size() ; i++)
		{
			double dxi = vm.deval(getShapeFunction(i), XI, p, dummy, 10.*default_derivation_delta) ;
			double deta = vm.deval(getShapeFunction(i), ETA, p, dummy,  10.*default_derivation_delta) ;
			double dtau = vm.deval(getShapeFunction(i),TIME_VARIABLE,p, dummy, 10.*default_derivation_delta) ;
			
 			xdxi += dxi*getBoundingPoint(i).getX() ;
 			ydxi += dxi*getBoundingPoint(i).getY() ;
			tdxi += dxi*getBoundingPoint(i).getT()  ;

 			xdeta += deta*getBoundingPoint(i).getX() ;
 			ydeta += deta*getBoundingPoint(i).getY() ;
			tdeta += deta*getBoundingPoint(i).getT()  ;

			xdtau += dtau*getBoundingPoint(i).getX() ;
			ydtau += dtau*getBoundingPoint(i).getY() ;
			tdtau += dtau*getBoundingPoint(i).getT()  ;
			
		}
		
		ret[0][0] = xdxi ; ret[0][1] = ydxi ; ret[0][2] = tdxi ;
		ret[1][0] = xdeta ; ret[1][1] = ydeta ; ret[1][2] = tdeta ;
		ret[2][0] = xdtau ;  ret[2][1] = ydtau ; ret[2][2] = tdtau;

		invert3x3Matrix(ret) ;
		if(cachedJinv.empty() && !isMoved())
			cachedJinv.push_back(ret) ;
//		ret.print() ;
	}
}
	
const GaussPointArray & TriElement::getGaussPoints()
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
	S[0][0] = getBoundingPoint(0).getX() ; S[0][1] = getBoundingPoint(factor).getX() ;  S[0][2] = getBoundingPoint(factor*2).getX() ; 
	S[1][0] = getBoundingPoint(0).getY() ; S[1][1] = getBoundingPoint(factor).getY() ;  S[1][2] = getBoundingPoint(factor*2).getY() ; 
	S[2][0] = 1 ; S[2][1] = 1 ;  S[2][2] = 1 ; 
	
	Vector v(3) ; 
	v[0] = p.getX() ;
	v[1] = p.getY() ;
	v[2] = 1 ;
	
	Vector coeff = inverse3x3Matrix( S) * v ;
	

	double tmin = getBoundingPoint(0).getT() ;
	double tmax = getBoundingPoint(getBoundingPoints().size()-1).getT() ;
	
	Point ret ;
	ret += Point(0,1,0,0)*coeff[0] ;
	ret += Point(0,0,0,0)*coeff[1] ;
	ret += Point(1,0,0,0)*coeff[2] ;	
	
	
	std::vector<double> instants ;
	for(size_t i = 0 ; i < timePlanes() ; i++)
	{
		instants.push_back(getBoundingPoint( i*getBoundingPoints().size() / timePlanes() ).getT()) ;
	}
	if(instants.size() > 1)
	{
		if(p.getT() < instants[1])
		{
			ret.getT() = -1. + (2./(timePlanes()-1))*(p.getT()-instants[0])/(instants[1]-instants[0]) ;
			return ret ;
		}
		else if(p.getT() > instants[instants.size()-2])
		{
			ret.getT() = 1. - (2./(timePlanes()-1))*(instants[instants.size()-1]-p.getT())/(instants[instants.size()-1]-instants[instants.size()-2]) ;
			return ret ;
		}
		else
		{
			for(size_t i = 1 ; i < instants.size() - 2 ; i++)
			{
				if(p.getT() > instants[i] && p.getT() < instants[i+1])
				{
					ret.getT() = -1. + 2.*((double) i )/ (timePlanes()-1) + (2./(timePlanes()-1))*(p.getT()-instants[i])/(instants[i+1]-instants[i]) ;
					return ret ;
				}
			}
		  
		}
	}
	
	return ret ;
}


std::valarray<std::valarray<Matrix> > TriElement::getNonLinearElementaryMatrix() 
{
	return std::valarray<std::valarray<Matrix> >() ;
}

Vector TriElement::getNonLinearForces()  
{
	
	return Vector(0) ;
}


const GaussPointArray & TetrahedralElement::genGaussPoints() 
{
	if(getCachedGaussPoints())
	{
//		std::cout <<getCachedGaussPoints()->gaussPoints.size() << std::endl ;
		return *getCachedGaussPoints() ;
	}
		
	size_t ordre=0 ;
	if(order == LINEAR)
		ordre = 1 ;
	else if(order == LINEAR_TIME_QUADRATIC  || order == LINEAR_TIME_LINEAR)
		ordre = 3 ;
	else if (order == CUBIC || order == QUADRATIC || order == CUBIC_TIME_LINEAR || order == QUADRATIC_TIME_LINEAR)
		ordre = 5 ;
	else if (order == QUADRIC || order == QUINTIC )
		ordre = 17 ;
	else
		ordre = 10 ;
	
	std::valarray< std::pair<Point, double> > fin(ordre);
	
	if(order == LINEAR )
	{
		fin[0] = std::pair<Point, double>(Point(0.25, 0.25, 0.25), 0.1666666666666667) ;
	}
	else if (order == CUBIC || order == QUADRATIC )
	{
		fin[0] = std::pair<Point, double>(Point(0.25, 0.25, 0.25), 0.1666666666666667*-.8) ;
		fin[1] = std::pair<Point, double>(Point(0.166666666666667, 0.166666666666667, 0.166666666666667), 0.1666666666666667*0.45) ;
		fin[2] = std::pair<Point, double>(Point(0.5, 0.166666666666667, 0.166666666666667), 0.1666666666666667*0.45) ;
		fin[3] = std::pair<Point, double>(Point(0.166666666666667, 0.5, 0.166666666666667), 0.1666666666666667*0.45) ;
		fin[4] = std::pair<Point, double>(Point(0.166666666666667, 0.166666666666667, 0.5), 0.1666666666666667*0.45) ;
	}
	else if (order == QUADRIC || order == QUINTIC )
	{
		fin[0] = std::pair<Point, double>(Point(0.25, 0.25, 0.25), 0.031403092789502) ;
		fin[1] = std::pair<Point, double>(Point(0.731636907957618, 0.089454364014127, 0.089454364014127), 0.011173097287728) ;
		fin[2] = std::pair<Point, double>(Point(0.089454364014127, 0.089454364014127, 0.089454364014127), 0.011173097287728) ;
		fin[3] = std::pair<Point, double>(Point(0.089454364014127, 0.731636907957618, 0.089454364014127), 0.011173097287728) ;
		fin[4] = std::pair<Point, double>(Point(0.089454364014127, 0.089454364014127, 0.731636907957618), 0.011173097287728) ;		
		fin[5] = std::pair<Point, double>(Point(0.132581099938466, 0.024540037929030, 0.421439431066252), 0.007547598727224) ;
		fin[6] = std::pair<Point, double>(Point(0.132581099938466, 0.421439431066252, 0.024540037929030), 0.007547598727224) ;
		fin[7] = std::pair<Point, double>(Point(0.132581099938466, 0.421439431066252, 0.421439431066252), 0.007547598727224) ;
		fin[8] = std::pair<Point, double>(Point(0.024540037929030, 0.132581099938466, 0.421439431066252), 0.007547598727224) ;
		fin[9] = std::pair<Point, double>(Point(0.024540037929030, 0.421439431066252, 0.132581099938466), 0.007547598727224) ;
		fin[10] = std::pair<Point, double>(Point(0.024540037929030, 0.421439431066252, 0.421439431066252), 0.007547598727224) ;
		fin[11] = std::pair<Point, double>(Point(0.421439431066252, 0.132581099938466, 0.024540037929030), 0.007547598727224) ;
		fin[12] = std::pair<Point, double>(Point(0.421439431066252, 0.132581099938466, 0.421439431066252), 0.007547598727224) ;
		fin[13] = std::pair<Point, double>(Point(0.421439431066252, 0.024540037929030, 0.132581099938466), 0.007547598727224) ;
		fin[14] = std::pair<Point, double>(Point(0.421439431066252, 0.024540037929030, 0.421439431066252), 0.007547598727224) ;
		fin[15] = std::pair<Point, double>(Point(0.421439431066252, 0.421439431066252, 0.132581099938466), 0.007547598727224) ;
		fin[16] = std::pair<Point, double>(Point(0.421439431066252, 0.421439431066252, 0.024540037929030), 0.007547598727224) ;
	}
	else if(order == LINEAR_TIME_LINEAR )
	{
		fin[0] = std::pair<Point, double>(Point(0.25, 0.25, 0.25,-std::sqrt(0.6)), 0.1666666666666667*5./9) ;
		fin[1] = std::pair<Point, double>(Point(0.25, 0.25, 0.25,0.), 0.1666666666666667*8./9) ;
		fin[2] = std::pair<Point, double>(Point(0.25, 0.25, 0.25,std::sqrt(0.6)), 0.1666666666666667*5./9) ;
	}
	else if(order == LINEAR_TIME_QUADRATIC )
	{
		fin[0] = std::pair<Point, double>(Point(0.25, 0.25, 0.25,-std::sqrt(0.6)), 0.1666666666666667*5./9) ;
		fin[1] = std::pair<Point, double>(Point(0.25, 0.25, 0.25,0.), 0.1666666666666667*8./9) ;
		fin[2] = std::pair<Point, double>(Point(0.25, 0.25, 0.25,std::sqrt(0.6)), 0.1666666666666667*5./9) ;
	}
	else if (order == CUBIC_TIME_LINEAR || order == QUADRATIC_TIME_LINEAR )
	{
		fin[0] = std::pair<Point, double>(Point(0.25, 0.25, 0.25), -0.5333333333333333333) ;
		fin[1] = std::pair<Point, double>(Point(0.16666666666666666667, 0.166666666666667, 0.1666666666666666667), 0.3) ;
		fin[2] = std::pair<Point, double>(Point(0.5, 0.166666666666666666667, 0.166666666666666666667), 0.3) ;
		fin[3] = std::pair<Point, double>(Point(0.16666666666666666667, 0.5, 0.166666666666667), 0.3) ;
		fin[4] = std::pair<Point, double>(Point(0.16666666666666666667, 0.1666666666666666667, 0.5), 0.3) ;
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
// 	else
// 	{
// 		std::cout << "this set of Gauss points is not implemented" << std::endl ;
// 		assert(false) ;
// 	}
	
	if( !isFather && isMoved())
	{
		for(size_t i = 0 ; i < fin.size() ; i++)
		{
			fin[i].second *= jacobianAtPoint(fin[i].first) ;
		}
	}
	else
	{
		Matrix J ;
// 		getInverseJacobianMatrix(fin[0].first, J);
		double j = volume()*3. ;
                if(order < CONSTANT_TIME_LINEAR)
                    j *= 2. ;
		for(size_t i = 0 ; i < fin.size() ; i++)
		{
			fin[i].second *= j;
		}

// 		double j = volume()/0.1666666666666666 ;
// 
// 		
// 		for(size_t i = 0 ; i < fin.size() ; i++)
// 		{
// 			fin[i].second*=j ;
// 		}
	}
	
	setCachedGaussPoints( new GaussPointArray(fin, order)) ;
	return *getCachedGaussPoints() ;
}

TetrahedralElement::TetrahedralElement( Point * p0,  Point * p1,  Point * p2, Point * p3 ) : Tetrahedron(p0, p1, p2, p3), moved(false)
{
	isFather = false ;
	this->order = LINEAR ;
	shapefunc = nullptr ;
}

TetrahedralElement::TetrahedralElement( Point * p0,  Point * p1,  Point * p2, Point * p3,  Point * p4,  Point * p5,  Point * p6, Point * p7 ) : Tetrahedron(p0, p1, p2, p3, p4, p5, p6, p7), moved(false) 
{
	isFather = false ;
	this->order = QUADRATIC ;
	shapefunc = nullptr ;
}

TetrahedralElement::TetrahedralElement(Order order ): moved(false)
{
	isFather = true ;
	setOrder( order );
	
	if(order == LINEAR)
	{
		
		
		shapefunc = new std::valarray<Function>(4) ;
// 		Matrix zero(2,2) ;
// 		std::valarray<Matrix> xi(zero, 2) ;
// 		xi[1][0][0] = 1 ;
// 		std::valarray<Matrix> eta(zero,2) ;
// 		eta[0][1][0] = 1 ;
// 		std::valarray<Matrix> zeta(zero,2) ;
// 		zeta[0][0][1] = 1 ;
// 		std::valarray<Matrix> f(zero,2) ;
// 		f[0][0][0] = 1 ;
// 		f[1][0][0] = -1 ;
// 		f[0][1][0] = -1 ;
// 		f[0][0][1] = -1 ;
			
			Function zero("0") ;
			Function one("1") ;
			Function mone("-1") ;
// 		//0
// 			(*shapefunc)[0] = Function("y") ;
// 			(*shapefunc)[0].setNumberOfDerivatives(2) ;
// 			(*shapefunc)[0].setDerivative( XI, zero) ;
// 			(*shapefunc)[0].setDerivative( ETA, one) ;
// 		//1
// 			(*shapefunc)[1] = Function("1 x - y -") ;
// 			(*shapefunc)[1].setNumberOfDerivatives(2) ;
// 			(*shapefunc)[1].setDerivative( XI, mone) ;
// 			(*shapefunc)[1].setDerivative( ETA, mone) ;
// 		//2
// 			(*shapefunc)[2] = Function("x") ;
// 			(*shapefunc)[2].setNumberOfDerivatives(2) ;
// 			(*shapefunc)[2].setDerivative( ETA, zero) ;
// 			(*shapefunc)[2].setDerivative( XI, one) ;
		
			//0
		(*shapefunc)[0] = Function("z") ;
		(*shapefunc)[0].setNumberOfDerivatives(3) ;
		(*shapefunc)[0].setDerivative( XI, zero) ;
		(*shapefunc)[0].setDerivative( ETA, zero) ;
		(*shapefunc)[0].setDerivative( ZETA, one) ;
			//1
		(*shapefunc)[1] = Function("1 x - y - z -") ;
		(*shapefunc)[1].setNumberOfDerivatives(3) ;
		(*shapefunc)[1].setDerivative( XI, mone) ;
		(*shapefunc)[1].setDerivative( ETA, mone) ;
		(*shapefunc)[1].setDerivative( ZETA, mone) ;
			//2
		(*shapefunc)[2] = Function("x") ;
		(*shapefunc)[2].setNumberOfDerivatives(3) ;
		(*shapefunc)[2].setDerivative( XI, one) ;
		(*shapefunc)[2].setDerivative( ETA, zero) ;
		(*shapefunc)[2].setDerivative( ZETA, zero) ;
			//3
		(*shapefunc)[3] = Function("y") ;
		(*shapefunc)[3].setNumberOfDerivatives(3) ;
		(*shapefunc)[3].setDerivative( XI, zero) ;
		(*shapefunc)[3].setDerivative( ETA, one) ;
		(*shapefunc)[3].setDerivative( ZETA, zero) ;
	}
	else if(order == QUADRATIC)
	{
		
		shapefunc = new std::valarray<Function>(10) ;
			
		(*shapefunc)[0] = Function("z 2 ^ 2 * z -") ;//z- z*2*(one-x-y-z) - x*z*2 - y*z*2 ;  // z
		(*shapefunc)[0].setNumberOfDerivatives(3) ;
		
		Function d("0") ; (*shapefunc)[0].setDerivative(XI, d) ;
		d = Function("0") ; (*shapefunc)[0].setDerivative(ETA, d) ;
		d = Function("4 z * 1 -") ; (*shapefunc)[0].setDerivative(ZETA, d) ;
		
		
		
		(*shapefunc)[1] = Function("z 4 * z 2 ^ 4 * - z x * 4 * - z y * 4 * -") ; //z*4*(one-x-y-z) ;
		(*shapefunc)[1].setNumberOfDerivatives(3) ;
		d = Function("-4 z *") ; (*shapefunc)[1].setDerivative(XI, d) ;
		d = Function("-4 z *") ; (*shapefunc)[1].setDerivative(ETA, d) ;
		d = Function("4 8 z * - x 4 * - y 4 * -") ; (*shapefunc)[1].setDerivative(ZETA, d) ;
		
		(*shapefunc)[2] = Function("1 z 3 * - 2 z 2 ^ * + z x * 4 * + z y * 4 * + x 3 * - x 2 ^ 2 * + x y * 4 * + y 3 * - y 2 ^ 2 * +") ; //one-x-y-z-(one-x-y-z)*(x+y+z)*2 ; // 0
		(*shapefunc)[2].setNumberOfDerivatives(3) ;
		d = Function("4 z * 3 - 4 x * + 4 y * +") ;(*shapefunc)[2].setDerivative(XI, d) ;
		d = Function("4 z * 4 x * + 3 - 4 y * + ") ;(*shapefunc)[2].setDerivative(ETA, d) ;
		d = Function("-3 4 z * + 4 x * + 4 y * +") ;(*shapefunc)[2].setDerivative(ZETA, d) ;
		
		
		(*shapefunc)[3] = Function("x 4 * x 2 ^ 4 * - z x * 4 * - x y * 4 * -") ; //x*4*(one-x-y-z) ;
		(*shapefunc)[3].setNumberOfDerivatives(3) ;
		d = Function("4 8 x * - z 4 * - 4 y * -") ;(*shapefunc)[3].setDerivative(XI, d) ;
		d = Function("-4 x *") ;                   (*shapefunc)[3].setDerivative(ETA, d) ;
		d = Function("-4 x *") ;                   (*shapefunc)[3].setDerivative(ZETA, d) ;
		
		(*shapefunc)[4] = Function("x 2 ^ 2 * x -") ; //x- x*2*(one-x-y-z) - x*z*2 - y*x*2 ; //x
		(*shapefunc)[4].setNumberOfDerivatives(3) ;
		d = Function("4 x * 1 -") ;(*shapefunc)[4].setDerivative(XI, d) ;
		d = Function("0") ;        (*shapefunc)[4].setDerivative(ETA, d) ;
		d = Function("0") ;        (*shapefunc)[4].setDerivative(ZETA, d) ;
		
		(*shapefunc)[5] = Function("x y * 4 *") ; //x*y*4 ; 
		(*shapefunc)[5].setNumberOfDerivatives(3) ;
		d = Function("y 4 *") ;(*shapefunc)[5].setDerivative(XI, d) ;
		d = Function("x 4 *") ;(*shapefunc)[5].setDerivative(ETA, d) ;
		d = Function("0") ;    (*shapefunc)[5].setDerivative(ZETA, d) ;
		
		(*shapefunc)[6] = Function("y 2 ^ 2 * y -") ; //y- y*2*(one-x-y-z) - y*z*2 - y*x*2 ;  //y
		(*shapefunc)[6].setNumberOfDerivatives(3) ;
		d = Function("0") ;(*shapefunc)[6].setDerivative(XI, d) ;
		d = Function("4 y * 1 -") ;(*shapefunc)[6].setDerivative(ETA, d) ;
		d = Function("0") ;(*shapefunc)[6].setDerivative(ZETA, d) ;
		
		(*shapefunc)[7] = Function("y z * 4 *") ; //y*z*4 ;
		(*shapefunc)[7].setNumberOfDerivatives(3) ;
		d = Function("0") ;(*shapefunc)[7].setDerivative(XI, d) ;
		d = Function("z 4 *") ;(*shapefunc)[7].setDerivative(ETA, d) ;
		d = Function("y 4 *") ;(*shapefunc)[7].setDerivative(ZETA, d) ;
		
		(*shapefunc)[8] = Function("y 4 * y 2 ^ 4 * - z y * 4 * - x y * 4 * -") ; //y*4*(one-x-y-z) ;
		(*shapefunc)[8].setNumberOfDerivatives(3) ;
		d = Function("-4 y *") ;(*shapefunc)[8].setDerivative(XI, d) ;
		d = Function("4 8 y * - 4 z * - 4 x * -") ;(*shapefunc)[8].setDerivative(ETA, d) ;
		d = Function("-4 y *") ;(*shapefunc)[8].setDerivative(ZETA, d) ;
		
		
		(*shapefunc)[9] = Function("x z * 4 *") ; //x*z*4 ;
		(*shapefunc)[9].setNumberOfDerivatives(3) ;
		d = Function("z 4 *") ;(*shapefunc)[9].setDerivative(XI, d) ;
		d = Function("0") ;(*shapefunc)[9].setDerivative(ETA, d) ;
		d = Function("x 4 *") ;(*shapefunc)[9].setDerivative(ZETA, d) ;
	}
	else if(order == LINEAR_TIME_LINEAR)
	{
		shapefunc = new std::valarray<Function>(8) ;
			Function z2("0") ;

			Function z1("0") ;

			z1.setNumberOfDerivatives(4) ;
			for(int i = 0 ; i < 4 ; i++)
			{
				z1.setDerivative( (const Variable) i, z2) ;
			}
			
			Function zero("0") ;
			zero.setNumberOfDerivatives(4) ;
			for(int i = 0 ; i < 4 ; i++)
			{
				zero.setDerivative( (const Variable) i, z1) ;
			}	
			Function one = zero +1 ;
			Function mone = zero-1 ;
			
			Function half("0.5") ;
			half.setNumberOfDerivatives(4) ;
			Function halfm("-0.5") ;
			halfm.setNumberOfDerivatives(4) ;
			for(int i = 0 ; i < 4 ; i++)
			{
				halfm.setDerivative( (const Variable) i, z1) ;
				half.setDerivative( (const Variable) i, z1) ;
			}	
			
			Function t0("0.5 t 0.5 * -") ;
			t0.setNumberOfDerivatives(4) ;
			t0.setDerivative( XI, zero) ;
			t0.setDerivative( ETA, zero) ;
			t0.setDerivative( ZETA, zero) ;
			t0.setDerivative( TIME_VARIABLE, halfm) ;

			Function t1("0.5 t 0.5 * +") ;
			t1.setNumberOfDerivatives(4) ;
			t1.setDerivative( XI, zero) ;
			t1.setDerivative( ETA, zero) ;
			t1.setDerivative( ZETA, zero) ;
			t1.setDerivative( TIME_VARIABLE, half) ;
			
			Function s0("z") ;
			Function s1("1 x - y - z -") ;
			Function s2("x") ;
			Function s3("y") ;
			s0.setNumberOfDerivatives(4) ;
			s0.setDerivative( XI, zero) ;
			s0.setDerivative( ETA, zero) ;
			s0.setDerivative( ZETA, one) ;
			s0.setDerivative( TIME_VARIABLE, zero) ;
			s1.setNumberOfDerivatives(4) ;
			s1.setDerivative( XI, mone) ;
			s1.setDerivative( ETA, mone) ;
			s1.setDerivative( ZETA, mone) ;
			s1.setDerivative( TIME_VARIABLE, zero) ;
			s2.setNumberOfDerivatives(4) ;
			s2.setDerivative( XI, one) ;
			s2.setDerivative( ETA, zero) ;
			s2.setDerivative( ZETA, zero) ;
			s2.setDerivative( TIME_VARIABLE, zero) ;
			s3.setNumberOfDerivatives(4) ;
			s3.setDerivative( XI, zero) ;
			s3.setDerivative( ETA, one) ;
			s3.setDerivative( ZETA, zero) ;
			s3.setDerivative( TIME_VARIABLE, zero) ;
			
		//0			
			(*shapefunc)[0] = s0*t0 ;
			(*shapefunc)[1] = s1*t0 ;
			(*shapefunc)[2] = s2*t0 ;
			(*shapefunc)[3] = s3*t0 ;
			(*shapefunc)[4] = s0*t1 ;
			(*shapefunc)[5] = s1*t1 ;
			(*shapefunc)[6] = s2*t1 ;
			(*shapefunc)[7] = s3*t1 ;
	}
	else if(order == LINEAR_TIME_QUADRATIC)
	{

// 		std::cout << "element order not implemented" << std::endl ;
// 		exit(0) ;

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
// 		std::cout << "element order not implemented" << std::endl ;
// 		exit(0) ;
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
		(*shapefunc)[12] = Function("1 z y x + + 3 * - z z * x x * y y * + + 2 * + z x *  z y * x y * + + 4 * + 1 t + 0.5 * *") ;
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
// 		std::cout << "element order not implemented" << std::endl ;
// 		exit(0) ;
		shapefunc = new std::valarray<Function>(30) ;
		
			//0
		(*shapefunc)[0] = Function("2 z * 1 - z * t 1 - t * 0.5 * *") ;
			//1
		(*shapefunc)[1] = Function("1 x - y - z - z * 4 * t 1 - t * 0.5 * *") ;
			//2
		(*shapefunc)[2] = Function("1 z y x + + 3 * - z z * x x * y y * + + 2 * + z x *  z y * x y * + + 4 * + t 1 - t * 0.5 * *") ;
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
		(*shapefunc)[12] = Function("1 z y x + + 3 * - z z * x x * y y * + + 2 * + z x *  z y * x y * + + 4 * + 1 t t * - *") ;
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
		(*shapefunc)[22] = Function("1 z y x + + 3 * - z z * x x * y y * + + 2 * + z x *  z y * x y * + + 4 * + t 1 + t * 0.5 * *") ;
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
	assert(this->jacobianAtPoint(Tetrahedron::getCenter()) > 0) ;
}
	
TetrahedralElement::TetrahedralElement(TetrahedralElement * parent, Tetrahedron * t)
{
	order =  parent->getOrder() ;
	
	shapefunc = parent->shapefunc ;
// 		std::copy(&parent->getShapeFunctions()[0], &parent->getShapeFunctions()[parent->getShapeFunctions().size()], &this->shapefunc[0]) ;
	for(size_t i =  0  ; i < this->size() ; i++)
	{
		delete &getPoint(i) ;
	}
	
	getInPoints().resize(t->getInPoints().size()) ;
	getBoundingPoints().resize(t->getBoundingPoints().size()) ;
	
	
	
	std::copy(&t->getInPoints()[0], &t->getInPoints()[t->getInPoints().size()],&getInPoints()[0] ) ;
	std::copy(&t->getBoundingPoints()[0], &t->getBoundingPoints()[t->getBoundingPoints().size()],&getBoundingPoints()[0] ) ;
}
	
std::valarray<std::valarray<Matrix> > & TetrahedralElement::getElementaryMatrix() 
{
	return cachedElementaryMatrix ;
}

std::valarray<std::valarray<Matrix> > & TetrahedralElement::getViscousElementaryMatrix() 
{
	return cachedViscousElementaryMatrix ;
}

std::valarray<std::valarray<Matrix> > TetrahedralElement::getNonLinearElementaryMatrix() 
{
	return std::valarray<std::valarray<Matrix> >() ;
}

Function TriElement::getXTransformAtCentralNodalTime() const
{
	if(shapefunc->size() == getBoundingPoints().size())
	{
		TriElement * tmp = new TriElement(LINEAR) ;
//		PointArray nodes = this->getBoundingPoints()[std::slice(0,3,1)] ;
		Function f = XTransform( getBoundingPoints()[std::slice(0,3,1)], tmp->getShapeFunctions()) ;
		delete tmp ;
		return f ;
	}
	return getXTransform() ;
}

Function TriElement::getYTransformAtCentralNodalTime() const
{
//	if(shapefunc->size() == getBoundingPoints().size())
//		return YTransform( this->getBoundingPoints()[std::slice(0,3,1)], getShapeFunctions()) ;
	if(shapefunc->size() == getBoundingPoints().size())
	{
		TriElement * tmp = new TriElement(LINEAR) ;
//		PointArray nodes = this->getBoundingPoints()[std::slice(0,3,1)] ;
		Function f = YTransform( getBoundingPoints()[std::slice(0,3,1)], tmp->getShapeFunctions()) ;
		delete tmp ;
		return f ;
	}
	return getYTransform() ;
}

Function TriElement::getXTransform() const
{
	if(shapefunc->size() == getBoundingPoints().size())
		return XTransform( this->getBoundingPoints(), getShapeFunctions()) ;
	return XTransform( this->getBoundingPoints(), TriElement(getOrder()).getShapeFunctions()) ;
}

Function TriElement::getYTransform() const
{
if(shapefunc->size() == getBoundingPoints().size())
		return YTransform( this->getBoundingPoints(), getShapeFunctions()) ;
	return YTransform( this->getBoundingPoints(), TriElement(getOrder()).getShapeFunctions()) ;
}


Function HexahedralElement::getXTransform() const
{
	return XTransform( this->getBoundingPoints(), HexahedralElement(getOrder()).getShapeFunctions()) ;
}

Function HexahedralElement::getYTransform() const
{
// 	if(shapefunc)
// 		return YTransform( this->getBoundingPoints(), getShapeFunctions()) ;
	return YTransform( this->getBoundingPoints(), HexahedralElement(getOrder()).getShapeFunctions()) ;
}

Function HexahedralElement::getZTransform() const
{
// 	if(shapefunc)
// 		return ZTransform( this->getBoundingPoints(), getShapeFunctions()) ;
	return ZTransform( this->getBoundingPoints(), HexahedralElement(getOrder()).getShapeFunctions()) ;
}

Vector TetrahedralElement::getNonLinearForces() 
{
	return Vector(0) ;
}
	
void TetrahedralElement::refresh(const TetrahedralElement * parent)
{
	if(!parent)
		return ;
	setOrder( parent->getOrder() );
// 	
// 	if(shapefunc)
// 	{
// 		isFather = false ;
// 		delete shapefunc ;
// 	}
	this->shapefunc = parent->shapefunc ;
}

void TetrahedralElement::print()  const
{
	std::cout << "I am a tetrahedron..." << std::endl ;
}

Function ElementarySurface::getXTransform() const
{
	return XTransform( this->getBoundingPoints(), this->getShapeFunctions()) ;
}

Function ElementarySurface::getYTransform() const
{
	return YTransform( this->getBoundingPoints(), this->getShapeFunctions()) ;
}

Function ElementarySurface::getTTransform() const
{
	return TTransform( this->getBoundingPoints(), this->getShapeFunctions()) ;
}

Function ElementarySurface::getdXTransform(Variable v) const
{
	return dXTransform( this->getBoundingPoints(), this->getShapeFunctions(),v) ;
}

Function ElementarySurface::getdYTransform(Variable v) const
{
	return dYTransform( this->getBoundingPoints(), this->getShapeFunctions(),v) ;
}

Function ElementarySurface::getdTTransform(Variable v) const
{
	return dTTransform( this->getBoundingPoints(), this->getShapeFunctions(),v) ;
}

double ElementarySurface::getdXTransform(Variable v, const Point p) const
{
	return dXTransform( this->getBoundingPoints(), this->getShapeFunctions(),v, p) ;
}

double ElementarySurface::getdYTransform(Variable v, const Point p) const
{
	return dYTransform( this->getBoundingPoints(), this->getShapeFunctions(),v, p) ;
}

double ElementarySurface::getdTTransform(Variable v, const Point p) const
{
	return dTTransform( this->getBoundingPoints(), this->getShapeFunctions(),v, p) ;
}

void ElementarySurface::setEnrichment( const Function & p, const Geometry * g)
{
	bool unique = true ;
	for(size_t i = 0 ;  i < enrichfunc.size() ; i++)
	{
		if (getEnrichmentFunction(i).getDofID() == p.getDofID())
		{
			unique = false ;
			break ;
		}
	}
	if(unique)
	{
		enrichmentUpdated = true ;
		enrichfunc.push_back(p) ;
		enrichmentSource.push_back(g) ;
	}
}

const  Function & ElementarySurface::getEnrichmentFunction(size_t i) const
{
	return this->enrichfunc[i];
}

//  Function & ElementarySurface::getEnrichmentFunction(size_t i) 
// {
// 	return this->enrichfunc[i];
// }


const std::valarray< Function >  & TriElement::getShapeFunctions() const
{
	return *shapefunc   ;
}

std::valarray< Function >  & TriElement::getShapeFunctions()
{
	return *shapefunc   ;
}


const std::vector<Function> & ElementarySurface::getEnrichmentFunctions() const
{
	return this->enrichfunc;
}


Point TetrahedralElement::inLocalCoordinates(const Point & p) const
{
	if(timePlanes() == 1)
	{

		if(order < QUADRATIC)
		{
			Matrix S(4,4) ;
			S[0][0] = this->getBoundingPoint(2).getX();
			S[0][1] = this->getBoundingPoint(3).getX(); 
			S[0][2] = this->getBoundingPoint(0).getX() ; 
			S[0][3] = this->getBoundingPoint(1).getX() ;  
		
			S[1][0] = this->getBoundingPoint(2).getY() ;
			S[1][1] = this->getBoundingPoint(3).getY();
			S[1][2] = this->getBoundingPoint(0).getY() ; 
			S[1][3] = this->getBoundingPoint(1).getY() ;  

			S[2][0] = this->getBoundingPoint(2).getZ() ;
			S[2][1] = this->getBoundingPoint(3).getZ();
			S[2][2] = this->getBoundingPoint(0).getZ() ; 
			S[2][3] = this->getBoundingPoint(1).getZ() ;  
		
			S[3][0] = 1 ; S[3][1] = 1 ;  S[3][2] = 1 ; S[3][3]= 1;
		
			Vector v(4) ; 
			v[0] = p.getX() ;
			v[1] = p.getY() ;
			v[2] = p.getZ() ;
			v[3] = 1 ;

			Vector coeff = inverse4x4Matrix(S) * v ;
		
		// 	VirtualMachine vm ;

		// 	Point t = ( Point(1.,0.,0.)*coeff[2] + Point(0.,1.,0.)*coeff[1] + Point(0.,0.,1.)*coeff[0] + Point(0.,0.,0.,p.getT())) ;
		// 	t.print();
		// 	Point test = Point(vm.eval(getXTransform(), t), vm.eval(getYTransform(),  t), vm.eval(getZTransform(),  t)) ;
		// 	test.print() ;
		// 	inLocalCoordinates(test) ;
		// 	std::cout << std::endl ;
			return Point(1.,0.,0.)*coeff[0] + Point(0.,1.,0.)*coeff[1] + Point(0.,0.,1.)*coeff[2] + Point(0.,0.,0.,p.getT()); 
		}
		else
		{
			Matrix S(4,4) ;
			S[0][0] = this->getBoundingPoint(4).getX();
			S[0][1] = this->getBoundingPoint(6).getX(); 
			S[0][2] = this->getBoundingPoint(0).getX() ; 
			S[0][3] = this->getBoundingPoint(2).getX() ;  
		
			S[1][0] = this->getBoundingPoint(4).getY() ;
			S[1][1] = this->getBoundingPoint(6).getY();
			S[1][2] = this->getBoundingPoint(0).getY() ; 
			S[1][3] = this->getBoundingPoint(2).getY() ;  

			S[2][0] = this->getBoundingPoint(4).getZ() ;
			S[2][1] = this->getBoundingPoint(6).getZ();
			S[2][2] = this->getBoundingPoint(0).getZ() ; 
			S[2][3] = this->getBoundingPoint(2).getZ() ;  
		
			S[3][0] = 1 ; S[3][1] = 1 ;  S[3][2] = 1 ; S[3][3]= 1;
		
			Vector v(4) ; 
			v[0] = p.getX() ;
			v[1] = p.getY() ;
			v[2] = p.getZ() ;
			v[3] = 1 ;

			Vector coeff = inverse4x4Matrix(S) * v ;
		
			double time = p.getT() ;
			if(timePlanes() > 1)
			{
				double t0 = this->getBoundingPoint(0).getT() ;
				double t1 = this->getBoundingPoint(this->getBoundingPoints().size()-1).getT() ;
				time = -1 + 2.*( p.getT()-t0)/(t1-t0) ;
			}	

		// 	VirtualMachine vm ;

		// 	Point t = ( Point(1.,0.,0.)*coeff[2] + Point(0.,1.,0.)*coeff[1] + Point(0.,0.,1.)*coeff[0] + Point(0.,0.,0.,p.getT())) ;
		// 	t.print();
		// 	Point test = Point(vm.eval(getXTransform(), t), vm.eval(getYTransform(),  t), vm.eval(getZTransform(),  t)) ;
		// 	test.print() ;
		// 	inLocalCoordinates(test) ;
		// 	std::cout << std::endl ;
			return Point(1.,0.,0.)*coeff[0] + Point(0.,1.,0.)*coeff[1] + Point(0.,0.,1.)*coeff[2] + Point(0.,0.,0.,time); 
		}
	}
	else
	{
			Matrix S(4,4) ;
			S[0][0] = this->getBoundingPoint(2).getX();
			S[0][1] = this->getBoundingPoint(3).getX(); 
			S[0][2] = this->getBoundingPoint(0).getX() ; 
			S[0][3] = this->getBoundingPoint(1).getX() ;  
		
			S[1][0] = this->getBoundingPoint(2).getY() ;
			S[1][1] = this->getBoundingPoint(3).getY();
			S[1][2] = this->getBoundingPoint(0).getY() ; 
			S[1][3] = this->getBoundingPoint(1).getY() ;  

			S[2][0] = this->getBoundingPoint(2).getZ() ;
			S[2][1] = this->getBoundingPoint(3).getZ();
			S[2][2] = this->getBoundingPoint(0).getZ() ; 
			S[2][3] = this->getBoundingPoint(1).getZ() ;  
		
			S[3][0] = 1 ; S[3][1] = 1 ;  S[3][2] = 1 ; S[3][3]= 1;
		
			Vector v(4) ; 
			v[0] = p.getX() ;
			v[1] = p.getY() ;
			v[2] = p.getZ() ;
			v[3] = 1 ;

			Vector coeff = inverse4x4Matrix(S) * v ;

			double time = p.getT() ;
			if(timePlanes() > 1)
			{
				double t0 = this->getBoundingPoint(0).getT() ;
				double t1 = this->getBoundingPoint(this->getBoundingPoints().size()-1).getT() ;
				time = -1 + 2.*( p.getT()-t0)/(t1-t0) ;
			}	

			return Point(1.,0.,0.)*coeff[0] + Point(0.,1.,0.)*coeff[1] + Point(0.,0.,1.)*coeff[2] + Point(0.,0.,0.,time); 

	}
}

Point HexahedralElement::inLocalCoordinates(const Point& p) const
{
	
	
	Matrix S(4,4) ;
	S[0][0] = this->getBoundingPoint(0).getX();
	S[0][1] = this->getBoundingPoint(4).getX();
	S[0][2] = this->getBoundingPoint(2).getX();
	S[0][3]=  this->getBoundingPoint(1).getX();
	
	
	S[1][0] = this->getBoundingPoint(0).getY();
	S[1][1] = this->getBoundingPoint(4).getY();
	S[1][2] = this->getBoundingPoint(2).getY();
	S[1][3]=  this->getBoundingPoint(1).getY();
	
	
	S[2][0] = this->getBoundingPoint(0).getZ();
	S[2][1] = this->getBoundingPoint(4).getZ();
	S[2][2] = this->getBoundingPoint(2).getZ();
	S[2][3]=  this->getBoundingPoint(1).getZ();
	
	
	S[3][0] = 1 ; S[3][1] = 1 ;  S[3][2] = 1 ; S[3][3]=1;
	
	Vector v(4) ; 
	v[0] = p.getX() ;
	v[1] = p.getY() ;
	v[2] = p.getZ() ;
	v[3] = 1 ;
	
	Vector coeff = inverse4x4Matrix(S) * v ;
	
	return Point(-1,-1,-1)*coeff[0] + Point(1,-1,-1)*coeff[1] + Point(-1,1,-1)*coeff[2]+  Point(-1,-1,1)*coeff[3];

}


const Function &  ElementarySurface::getShapeFunction(size_t i) const
{
	return getShapeFunctions()[i] ;
}

Function &  ElementarySurface::getShapeFunction(size_t i) 
{
	return getShapeFunctions()[i] ;
}

void ElementarySurface::compileAndPrecalculate()
{
	for(size_t k = 0 ; k < getShapeFunctions().size() ; k++)
	{
		std::vector<Variable> vars ;
		vars.push_back(XI) ;
		vars.push_back(ETA) ;
		if(order > CONSTANT_TIME_LINEAR)
			vars.push_back(TIME_VARIABLE) ;
		(*static_cast<TriElement *>(this)->shapefunc)[k].preCalculate(getGaussPoints(), vars) ;
	}
}

void ElementaryVolume::compileAndPrecalculate()
{
	for(size_t k = 0 ; k < getShapeFunctions().size() ; k++)
	{
		std::vector<Variable> vars ;
		vars.push_back(XI) ;
		vars.push_back(ETA) ;
		vars.push_back(ZETA) ;
		if(order > CONSTANT_TIME_LINEAR)
			vars.push_back(TIME_VARIABLE) ;
		(*static_cast<TetrahedralElement *>(this)->shapefunc)[k].preCalculate(getGaussPoints(), vars) ;
	}
}

// Function & ElementarySurface::getShapeFunction(size_t i) 
// {
// 	return (*shapefunc)[i] ;
// }

Order ElementarySurface::getOrder() const
{
	return order ;
}

Function ElementaryVolume::getdXTransform(Variable v) const
{
	return dXTransform( getBoundingPoints(), getShapeFunctions(),v) ;
}

Function ElementaryVolume::getdYTransform(Variable v) const
{
	return dYTransform( getBoundingPoints(), getShapeFunctions(),v) ;
}

Function ElementaryVolume::getdZTransform(Variable v) const
{
	return dZTransform( getBoundingPoints(), getShapeFunctions(),v) ;
}

Function ElementaryVolume::getdTTransform(Variable v) const
{
	return dTTransform( getBoundingPoints(), getShapeFunctions(),v) ;
}


double ElementaryVolume::getdXTransform(Variable v, const Point & p) const
{
	
	return dXTransform( getBoundingPoints(),getShapeFunctions(),v, p) ;
}

double ElementaryVolume::getdYTransform(Variable v, const Point & p) const
{
	return dYTransform( getBoundingPoints(), getShapeFunctions(),v, p) ;
}

double ElementaryVolume::getdZTransform(Variable v, const Point & p) const
{
	return dZTransform( getBoundingPoints(), getShapeFunctions(),v, p) ;
}

double ElementaryVolume::getdTTransform(Variable v, const Point & p) const
{
	return dTTransform( getBoundingPoints(), getShapeFunctions(),v, p) ;
}


Function TetrahedralElement::getdXTransform(Variable v) const
{
	if(getShapeFunctions().size() != getBoundingPoints().size())
		return dXTransform( getBoundingPoints(), TetrahedralElement(getOrder()).getShapeFunctions(),v) ;
	
	return dXTransform( getBoundingPoints(), getShapeFunctions(),v) ;
}

Function TetrahedralElement::getdYTransform(Variable v) const
{
	if(getShapeFunctions().size() != getBoundingPoints().size())
		return dYTransform( getBoundingPoints(), TetrahedralElement(getOrder()).getShapeFunctions(),v) ;
	
	return dYTransform( getBoundingPoints(), getShapeFunctions(),v) ;
}

Function TetrahedralElement::getdZTransform(Variable v) const
{
	if(getShapeFunctions().size() != getBoundingPoints().size())
		return dZTransform( getBoundingPoints(), TetrahedralElement(getOrder()).getShapeFunctions(),v) ;
	
	return dZTransform( getBoundingPoints(), getShapeFunctions(),v) ;
}

Function TetrahedralElement::getdTTransform(Variable v) const
{
	if(getShapeFunctions().size() != getBoundingPoints().size())
		return dTTransform( getBoundingPoints(), TetrahedralElement(getOrder()).getShapeFunctions(),v) ;
	
	return dTTransform( getBoundingPoints(), getShapeFunctions(),v) ;
}


double TetrahedralElement::getdXTransform(Variable v, const Point & p) const
{
	if(getShapeFunctions().size() != getBoundingPoints().size())
		return dXTransform(getBoundingPoints(), TetrahedralElement(getOrder()).getShapeFunctions(),v, p) ;
	
	return dXTransform( getBoundingPoints(), getShapeFunctions(),v, p) ;
}

double TetrahedralElement::getdYTransform(Variable v, const Point & p) const
{
	if(getShapeFunctions().size() != getBoundingPoints().size())
		return dYTransform( getBoundingPoints(), TetrahedralElement(getOrder()).getShapeFunctions(),v, p) ;
	
	return dYTransform( getBoundingPoints(), getShapeFunctions(),v, p) ;
}

double TetrahedralElement::getdZTransform(Variable v, const Point & p) const
{
	if(getShapeFunctions().size() != getBoundingPoints().size())
		return dZTransform( getBoundingPoints(), TetrahedralElement(getOrder()).getShapeFunctions(),v, p) ;
	
	return dZTransform( getBoundingPoints(),getShapeFunctions(),v, p) ;
}

double TetrahedralElement::getdTTransform(Variable v, const Point & p) const
{
	if(getShapeFunctions().size() != getBoundingPoints().size())
		return dTTransform( getBoundingPoints(), TetrahedralElement(getOrder()).getShapeFunctions(),v, p) ;
	
	return dTTransform( getBoundingPoints(), getShapeFunctions(),v, p) ;
}

std::vector<size_t> ElementaryVolume::clearEnrichment(const Geometry * g)
{
	std::vector<size_t> ret ;
	std::vector<Function> newFunc ;
	std::vector<Geometry *> newSource ;

	for(size_t i = 0 ; i < enrichmentSource.size() ; i++)
	{
		if(enrichmentSource[i] != g)
		{
			
			newFunc.push_back(enrichfunc[i]) ;
			newSource.push_back(enrichmentSource[i]) ;
		}
		else
			ret.push_back(enrichfunc[i].getDofID());
	}
	blendfunc.clear();
	enrichfunc = newFunc ;
	enrichmentSource = newSource ;
	
	return ret ;
}

std::vector<size_t> ElementarySurface::clearAllEnrichment()
{
	std::vector<size_t> ret ;
	std::vector<Function> newFunc ;
	std::vector<const Geometry *> newSource ;

	for(size_t i = 0 ; i < enrichmentSource.size() ; i++)
	{
		ret.push_back(enrichfunc[i].getDofID());
	}

	enrichfunc = newFunc ;
	enrichmentSource = newSource ;
	
	return ret ;
}

std::vector<size_t> ElementaryVolume::clearAllEnrichment()
{
	std::vector<size_t> ret ;
	std::vector<Function> newFunc ;
	std::vector<Geometry *> newSource ;

	for(size_t i = 0 ; i < enrichmentSource.size() ; i++)
	{
		ret.push_back(enrichfunc[i].getDofID());
	}

	enrichfunc = newFunc ;
	enrichmentSource = newSource ;
	
	return ret ;
}

std::vector<size_t> ElementarySurface::clearEnrichment(const Geometry * g)
{
	std::vector<size_t> ret ;
	std::vector<Function> newFunc ;
	std::vector<const Geometry *> newSource ;

	for(size_t i = 0 ; i < enrichmentSource.size() ; i++)
	{
		if(enrichmentSource[i] != g)
		{
			
			newFunc.push_back(enrichfunc[i]) ;
			newSource.push_back(enrichmentSource[i]) ;
		}
		else
			ret.push_back(enrichfunc[i].getDofID());
	}
	enrichfunc.clear() ;
	blendfunc.clear();
	enrichfunc = newFunc ;
	enrichmentSource = newSource ;
	
	return ret ;
}

ElementaryVolume::ElementaryVolume(bool f ) 
{
	this->behaviour = nullptr ;
	this->nonlinbehaviour = nullptr ;
	enrichmentUpdated = true ;
	behaviourUpdated = true ;
}
////

ElementaryVolume::~ElementaryVolume()
{
	delete behaviour ;
}

Function ElementaryVolume::jacobian() const 
{
	
	Function xdxi = getXTransform().d(XI) ;
	Function ydxi = getYTransform().d(XI) ;
	Function zdxi = getZTransform().d(XI) ;
	Function xdeta = getXTransform().d(ETA) ;
	Function ydeta = getYTransform().d(ETA) ;
	Function zdeta = getZTransform().d(ETA) ;
	Function xdzeta = getXTransform().d(ZETA) ;
	Function ydzeta = getYTransform().d(ZETA) ;
	Function zdzeta = getZTransform().d(ZETA) ;
	
	Function ret = xdxi*ydeta*zdzeta + zdeta*xdzeta*ydxi + ydzeta*zdxi*xdeta  -
						xdxi*ydzeta*zdeta - xdeta*ydxi*zdzeta - xdzeta*ydeta*zdxi ;
	if(order < CONSTANT_TIME_LINEAR)
		return ret ;	

	return ret * getTTransform().d(TIME_VARIABLE) ;
}

double ElementaryVolume::jacobianAtPoint(const Point & p) const 
{
	
	if(order < CONSTANT_TIME_LINEAR)
	{

		double xdxi = this->getdXTransform(XI,p) ;
		double ydxi = this->getdYTransform(XI,p) ;
		double zdxi = this->getdZTransform(XI,p) ;
		
		double xdeta = this->getdXTransform(ETA,p) ;
		double ydeta = this->getdYTransform(ETA,p) ;
		double zdeta = this->getdZTransform(ETA,p) ;
		
		double xdzeta = this->getdXTransform(ZETA,p) ;
		double ydzeta = this->getdYTransform(ZETA,p) ;
		double zdzeta = this->getdZTransform(ZETA,p) ;
// 		VirtualMachine vm ;
// 		TetrahedralElement father(order) ;
// 		for(size_t i = 0 ; i < getBoundingPoints().size() ; i++)
// 		{
// 			double dxi = vm.deval(father.getShapeFunction(i), XI, p.getX(), p.getY(). p.getZ()) ;
// 			double deta = vm.deval(father.getShapeFunction(i), ETA, p.getX(), p.getY(). p.getZ()) ;
// 			double dzeta = vm.deval(father.getShapeFunction(i), ZETA, p.getX(), p.getY(). p.getZ()) ;
// 			xdxi += dxi*getBoundingPoint(i).getX() ;
// 			ydxi += dxi*getBoundingPoint(i).getY() ;
// 			zdxi += dxi*getBoundingPoint(i).getZ() ;
// 
// 			xdeta += deta*getBoundingPoint(i).getX() ;
// 			ydeta += deta*getBoundingPoint(i).getY() ;
// 			zdeta += deta*getBoundingPoint(i).getZ() ;
// 
// 			xdzeta += dzeta*getBoundingPoint(i).getX() ;
// 			ydzeta += dzeta*getBoundingPoint(i).getY() ;
// 			zdzeta += dzeta*getBoundingPoint(i).getZ() ;
// 		}
		
		return xdxi*ydeta*zdzeta + zdeta*xdzeta*ydxi + ydzeta*zdxi*xdeta  -
			xdxi*ydzeta*zdeta - xdeta*ydxi*zdzeta - xdzeta*ydeta*zdxi ;
	}
	else
	{
		Matrix  J0(4,4) ;
		
		double xdxi = getdXTransform(XI,p) ;
		double ydxi = getdYTransform(XI,p) ;
		double zdxi = getdZTransform(XI,p) ;
// 		double tdxi = this->getdTTransform(XI,p) ;
		
		double xdeta = getdXTransform(ETA,p) ;
		double ydeta = getdYTransform(ETA,p) ;
		double zdeta = getdZTransform(ETA,p) ;
// 		double tdeta = this->getdTTransform(ETA,p) ;
		
		double xdzeta = getdXTransform(ZETA,p) ;
		double ydzeta = getdYTransform(ZETA,p) ;
		double zdzeta = getdZTransform(ZETA,p) ;
// 		double tdzeta = this->getdTTransform(ZETA,p) ;
		
// 		double xdtheta = this->getdXTransform(TIME_VARIABLE,p) ;
// 		double ydtheta = this->getdYTransform(TIME_VARIABLE,p) ;
// 		double zdtheta = this->getdZTransform(TIME_VARIABLE,p) ;
		double tdtheta = getdTTransform(TIME_VARIABLE,p) ;
		
		J0[0][0] = xdxi ; J0[0][1] = ydxi ; J0[0][2] = zdxi ; J0[0][3] = 0; 
		J0[1][0] = xdeta ; J0[1][1] = ydeta ; J0[1][2] = zdeta ; J0[1][3] = 0;
		J0[2][0] = xdzeta ; J0[2][1] = ydzeta ; J0[2][2] = zdzeta ; J0[2][3] = 0;
		J0[3][0] = 0 ; J0[3][1] = 0 ; J0[3][2] = 0 ; J0[3][3] = tdtheta;
		
		return det(J0) ;
	}
	
	
	VirtualMachine vm ;
	return vm.eval(jacobian(), p) ;
}

Function ElementaryVolume::getXTransform() const
{
	return XTransform( getBoundingPoints(), getShapeFunctions()) ;
}

Function  ElementaryVolume::getYTransform() const
{
	return YTransform( getBoundingPoints(), getShapeFunctions()) ;
}

Function  ElementaryVolume::getZTransform() const
{
	return ZTransform( getBoundingPoints(), getShapeFunctions()) ;
}

Function  ElementaryVolume::getTTransform() const
{
	return TTransform( getBoundingPoints(), getShapeFunctions()) ;
}

Function TetrahedralElement::getXTransform() const
{
	if(getBoundingPoints().size() != getShapeFunctions().size())
		return XTransform( getBoundingPoints(), TetrahedralElement(getOrder()).getShapeFunctions()) ;
	
	return XTransform( getBoundingPoints(), getShapeFunctions()) ;
}

Function  TetrahedralElement::getYTransform() const
{
	if(getBoundingPoints().size() != getShapeFunctions().size())
		return YTransform( getBoundingPoints(), TetrahedralElement(getOrder()).getShapeFunctions()) ;
	
	return YTransform( getBoundingPoints(), getShapeFunctions()) ;
}

Function  TetrahedralElement::getZTransform() const
{
	if(getBoundingPoints().size() != getShapeFunctions().size())
		return ZTransform( getBoundingPoints(), TetrahedralElement(getOrder()).getShapeFunctions()) ;
	
	return ZTransform( getBoundingPoints(), getShapeFunctions()) ;
}

Function  TetrahedralElement::getTTransform() const
{
	if(getBoundingPoints().size() != getShapeFunctions().size())
		return TTransform( getBoundingPoints(), TetrahedralElement(getOrder()).getShapeFunctions()) ;
	
	return TTransform( getBoundingPoints(), getShapeFunctions()) ;
}

void ElementaryVolume::setEnrichment(const Function & p, Geometry * g)
{
	bool unique = true ;
	for(size_t i = 0 ;  i < enrichfunc.size() ; i++)
	{
		if (getEnrichmentFunction(i).getDofID() == p.getDofID())
		{
			unique = false ;
			break ;
		}
	}
	if(unique)
	{
		enrichfunc.push_back(p) ;
		enrichmentSource.push_back(g) ;
	}
	
	
}

const Function& ElementaryVolume::getEnrichmentFunction(size_t i)  const
{
	return enrichfunc[i] ;
}

// Function& ElementaryVolume::getEnrichmentFunction(size_t i) 
// {
// 	return this->enrichfunc[i] ;
// }

const std::valarray<Function  >  & TetrahedralElement::getShapeFunctions() const
{
	return *shapefunc ;
}

std::valarray<Function  >  & TetrahedralElement::getShapeFunctions()
{
	return *shapefunc ;
}

const std::valarray<Function  >  & HexahedralElement::getShapeFunctions() const
{
	return *shapefunc ;
}

std::valarray<Function  >  & HexahedralElement::getShapeFunctions()
{
	return *shapefunc ;
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
		if(getBoundingPoint(i).getId() >= 0)
			ret.push_back(getBoundingPoint(i).getId()) ;
		else
			std::cout << "negative ID, check numbering !" << std::endl ;
	}
	for (size_t i = 0 ; i < getEnrichmentFunctions().size() ; i++)
	{
		if(getEnrichmentFunction(i).getDofID() >= 0)
			ret.push_back(getEnrichmentFunction(i).getDofID()) ;
		else
			std::cout << "negative ID, check numbering !" << std::endl ;
	}
	
	return ret ;
}

const std::vector< Function> & ElementaryVolume::getEnrichmentFunctions() const
{
	return this->enrichfunc;
}

// std::vector< Function> & ElementaryVolume::getEnrichmentFunctions() 
// {
// 	return this->enrichfunc;
// }


std::valarray<std::valarray<Matrix> > & HexahedralElement::getElementaryMatrix() 
{
	if(!cachedElementaryMatrix.size() == 0)
		return cachedElementaryMatrix ;

		GaussPointArray gp  = this->getGaussPoints(); 
		std::valarray<Matrix> Jinv ;
		std::vector<size_t> dofs = getDofIds() ;


		Matrix J(3,3) ;
		getInverseJacobianMatrix( gp.gaussPoints[0].first, J ) ;
		Jinv.resize( gp.gaussPoints.size(), J) ;


		

		std::valarray< Matrix > v_j(Matrix(3,3), dofs.size()) ;
		cachedElementaryMatrix.resize(dofs.size(), v_j);

	
	VirtualMachine vm ;
	for(size_t i = 0 ; i < getShapeFunctions().size() ; i++)
	{
		 behaviour->apply(getShapeFunction(i), getShapeFunction(i),gp, Jinv,cachedElementaryMatrix[i][i], &vm) ;
		
		for(size_t j = i+1 ; j < getShapeFunctions().size() ; j++)
		{
			behaviour->apply(getShapeFunction(i), getShapeFunction(j),gp, Jinv,cachedElementaryMatrix[i][j], &vm) ;
			behaviour->apply(getShapeFunction(j), getShapeFunction(i),gp, Jinv,cachedElementaryMatrix[j][i], &vm) ;
		}
		for(size_t j = 0 ; j < getEnrichmentFunctions().size() ; j++)
		{
			behaviour->apply(getShapeFunction(i), getEnrichmentFunction(j),gp, Jinv,cachedElementaryMatrix[i][j+getShapeFunctions().size()], &vm) ;
			behaviour->apply(getEnrichmentFunction(j), getShapeFunction(i),gp, Jinv, cachedElementaryMatrix[j+getShapeFunctions().size()][i], &vm) ;
		}
	}
	for(size_t i = 0 ; i < getEnrichmentFunctions().size() ; i++)
	{
		 behaviour->apply(getEnrichmentFunction(i), getEnrichmentFunction(i),gp, Jinv, cachedElementaryMatrix[i+getShapeFunctions().size()][i+getShapeFunctions().size()], &vm) ;
		
		for(size_t j = 0 ; j < getEnrichmentFunctions().size() ; j++)
		{
			behaviour->apply(getEnrichmentFunction(i), getEnrichmentFunction(j),gp, Jinv, cachedElementaryMatrix[i+getShapeFunctions().size()][j+getShapeFunctions().size()], &vm) ;
			behaviour->apply(getEnrichmentFunction(j), getEnrichmentFunction(i),gp, Jinv, cachedElementaryMatrix[j+getShapeFunctions().size()][i+getShapeFunctions().size()], &vm) ;
		}
	}

// 		}
// 	}
		return cachedElementaryMatrix ;
// 		}
		

	}


std::valarray<std::valarray<Matrix> > & HexahedralElement::getViscousElementaryMatrix() 
{
	return cachedViscousElementaryMatrix ;
}


const Function  & ElementaryVolume::getShapeFunction(size_t i) const
{
	return getShapeFunctions()[i] ;
}

Function  & ElementaryVolume::getShapeFunction(size_t i)
{
	return getShapeFunctions()[i] ;
}


void ElementaryVolume::getInverseJacobianMatrix(const Point & p, Matrix & ret) 
{
	if(order < CONSTANT_TIME_LINEAR)
	{
		if(ret.isNull())
			ret.resize(3,3) ;
		
		double xdxi = this->getdXTransform(XI,p) ;
		double ydxi = this->getdYTransform(XI,p) ;
		double zdxi = this->getdZTransform(XI,p) ;
		
		double xdeta = this->getdXTransform(ETA,p) ;
		double ydeta = this->getdYTransform(ETA,p) ;
		double zdeta = this->getdZTransform(ETA,p) ;
		
		double xdzeta = this->getdXTransform(ZETA,p) ;
		double ydzeta = this->getdYTransform(ZETA,p) ;
		double zdzeta = this->getdZTransform(ZETA,p) ;
		
		ret[0][0] = xdxi ; ret[0][1] = ydxi ; ret[0][2] = zdxi ; 
		ret[1][0] = xdeta ; ret[1][1] = ydeta ; ret[1][2] = zdeta ;
		ret[2][0] = xdzeta ; ret[2][1] = ydzeta ; ret[2][2] = zdzeta ;
		invert3x3Matrix(ret) ;
	}
	else
	{
		if(ret.isNull())
			ret.resize(4,4) ;
		
		double xdxi = this->getdXTransform(XI,p) ;
		double ydxi = this->getdYTransform(XI,p) ;
		double zdxi = this->getdZTransform(XI,p) ;
// 		double tdxi = this->getdTTransform(XI,p) ;
		
		double xdeta = this->getdXTransform(ETA,p) ;
		double ydeta = this->getdYTransform(ETA,p) ;
		double zdeta = this->getdZTransform(ETA,p) ;
// 		double tdeta = this->getdTTransform(ETA,p) ;
		
		double xdzeta = this->getdXTransform(ZETA,p) ;
		double ydzeta = this->getdYTransform(ZETA,p) ;
		double zdzeta = this->getdZTransform(ZETA,p) ;
// 		double tdzeta = this->getdTTransform(ZETA,p) ;
		
// 		double xdtheta = this->getdXTransform(TIME_VARIABLE,p) ;
// 		double ydtheta = this->getdYTransform(TIME_VARIABLE,p) ;
// 		double zdtheta = this->getdZTransform(TIME_VARIABLE,p) ;
		double tdtheta = this->getdTTransform(TIME_VARIABLE,p) ;
		
		ret[0][0] = xdxi ; ret[0][1] = ydxi ; ret[0][2] = zdxi ; ret[0][3] = 0; 
		ret[1][0] = xdeta ; ret[1][1] = ydeta ; ret[1][2] = zdeta ; ret[1][3] = 0;
		ret[2][0] = xdzeta ; ret[2][1] = ydzeta ; ret[2][2] = zdzeta ; ret[2][3] = 0;
		ret[3][0] = 0 ; ret[3][1] = 0 ; ret[3][2] = 0 ; ret[3][3] = tdtheta;


		ret = inverse4x4Matrix(ret) ;
	}

}

void TetrahedralElement::getInverseJacobianMatrix(const Point & p, Matrix & ret) 
{
    if(order < CONSTANT_TIME_LINEAR)
    {

        if(ret.isNull()|| ret.size() != 9)
            ret.resize(3,3) ;
        
        double xdxi = 0 ;//this->getdXTransform(XI,p) ;
        double ydxi = 0 ;//this->getdYTransform(XI,p) ;
        double zdxi = 0 ;//this->getdYTransform(XI,p) ;
        
        double xdeta = 0 ;//this->getdXTransform(ETA,p) ;
        double ydeta = 0 ;//this->getdYTransform(ETA,p) ;
        double zdeta = 0 ;//this->getdYTransform(ETA,p) ;

        double xdzeta = 0 ;//this->getdXTransform(ETA,p) ;
        double ydzeta = 0 ;//this->getdYTransform(ETA,p) ;
        double zdzeta = 0 ;//this->getdYTransform(ETA,p) ;
        
        VirtualMachine vm ;
        TetrahedralElement * father = nullptr ;

        std::valarray<Function> * functions = shapefunc ;
        if(!shapefunc)
        {
            father = new TetrahedralElement(order) ;
            functions = father->shapefunc ;
        }
        
        for(size_t i = 0 ; i < getBoundingPoints().size() ; i++)
        {
            double dxi = vm.deval((*functions)[i], XI, p) ;
            double deta = vm.deval((*functions)[i], ETA, p) ;
            double dzeta = vm.deval((*functions)[i], ZETA, p) ;
            
            xdxi += dxi*getBoundingPoint(i).getX() ;
            ydxi += dxi*getBoundingPoint(i).getY() ;
            zdxi += dxi*getBoundingPoint(i).getZ() ;

            xdeta += deta*getBoundingPoint(i).getX() ;
            ydeta += deta*getBoundingPoint(i).getY() ;
            zdeta += deta*getBoundingPoint(i).getZ() ;
            
            xdzeta += dzeta*getBoundingPoint(i).getX() ;
            ydzeta += dzeta*getBoundingPoint(i).getY() ;
            zdzeta += dzeta*getBoundingPoint(i).getZ() ;
        }
        ret[0][0] = xdxi ; ret[0][1] = ydxi ;  ret[0][2] = zdxi ; 
        ret[1][0] = xdeta ; ret[1][1] = ydeta ; ret[1][2] = zdeta ;
        ret[2][0] = xdzeta ; ret[2][1] = ydzeta ; ret[2][2] = zdzeta ;
        invert3x3Matrix(ret) ;
        delete father ;
//      ret.print() ;
//      exit(0) ;
    }
    else
    {
/*      if(!isMoved() && !cachedJinv.empty())
        {
            if(ret.isNull()  || ret.size() != 9)
                ret.resize(3,3) ;
            ret.array() = cachedJinv[0].array() ;
//          std::cout << ret.numCols() << "\t" << ret.numRows() << std::endl ;
//          ret.print() ;
            return ;
        }*/
        
        if(ret.isNull() || ret.size() != 16)
            ret.resize(4,4) ;

        double xdxi = 0 ;//this->getdXTransform(XI,p) ;
        double ydxi = 0 ;//this->getdYTransform(XI,p) ;
        double zdxi = 0 ;//this->getdYTransform(XI,p) ;
        double tdxi = 0 ;

        double xdeta = 0 ;//this->getdXTransform(ETA,p) ;
        double ydeta = 0 ;//this->getdYTransform(ETA,p) ;
        double zdeta = 0 ;//this->getdYTransform(ETA,p) ;
        double tdeta = 0 ;
        
        double xdzeta = 0 ;//this->getdXTransform(ETA,p) ;
        double ydzeta = 0 ;//this->getdYTransform(ETA,p) ;
        double zdzeta = 0 ;//this->getdYTransform(ETA,p) ;
        double tdzeta = 0 ;


        double xdtau = 0 ;
        double ydtau = 0 ;
        double zdtau = 0 ;
        double tdtau = 0 ;

        VirtualMachine vm ;
        
        Point dummy ;
        for(size_t i = 0 ; i < getBoundingPoints().size() ; i++)
        {
            double dxi   = vm.deval(getShapeFunction(i), XI, p, dummy, 10.*default_derivation_delta) ;
            double deta  = vm.deval(getShapeFunction(i), ETA, p, dummy,  10.*default_derivation_delta) ;
            double dzeta = vm.deval(getShapeFunction(i), ZETA, p, dummy,  10.*default_derivation_delta) ;
            double dtau  = vm.deval(getShapeFunction(i),TIME_VARIABLE,p, dummy, 10.*default_derivation_delta) ;
            
            xdxi += dxi*getBoundingPoint(i).getX() ;
            ydxi += dxi*getBoundingPoint(i).getY() ;
            zdxi += dxi*getBoundingPoint(i).getZ() ;
            tdxi += dxi*getBoundingPoint(i).getT()  ;

            xdeta += deta*getBoundingPoint(i).getX() ;
            ydeta += deta*getBoundingPoint(i).getY() ;
            zdeta += deta*getBoundingPoint(i).getZ() ;
            tdeta += deta*getBoundingPoint(i).getT()  ;
            
            xdzeta += dzeta*getBoundingPoint(i).getX() ;
            ydzeta += dzeta*getBoundingPoint(i).getY() ;
            zdzeta += dzeta*getBoundingPoint(i).getZ() ;
            tdzeta += dzeta*getBoundingPoint(i).getT()  ;

            xdtau += dtau*getBoundingPoint(i).getX() ;
            ydtau += dtau*getBoundingPoint(i).getY() ;
            zdtau += dtau*getBoundingPoint(i).getZ() ;
            tdtau += dtau*getBoundingPoint(i).getT()  ;
            
        }
        
        ret[0][0] = xdxi ; ret[0][1] = ydxi ;  ret[0][2] = zdxi ; ret[0][3] = tdxi ;
        ret[1][0] = xdeta ; ret[1][1] = ydeta ;  ret[1][2] = zdeta ;ret[1][3] = tdeta ;
        ret[2][0] = xdzeta ; ret[2][1] = ydzeta ; ret[2][2] = zdzeta ;ret[2][3] = tdzeta ;
        ret[3][0] = xdtau ;  ret[3][1] = ydtau ;  ret[3][2] = zdtau ;ret[3][3] = tdtau;

        ret = inverse4x4Matrix(ret) ;

    }
    
// 	if(getOrder() < CONSTANT_TIME_LINEAR)
// 	{
// 		if(ret.isNull())
// 			ret.resize(3,3) ;
// 		
// 		double xdxi = 0 ;//this->getdXTransform(XI,p) ;
// 		double ydxi = 0 ;//this->getdYTransform(XI,p) ;
// 		double zdxi = 0 ;//this->getdZTransform(XI,p) ;
// 		
// 		double xdeta = 0 ;//this->getdXTransform(ETA,p) ;
// 		double ydeta = 0 ;//this->getdYTransform(ETA,p) ;
// 		double zdeta = 0 ;//this->getdZTransform(ETA,p) ;
// 		
// 		double xdzeta = 0 ;//this->getdXTransform(ZETA,p) ;
// 		double ydzeta = 0 ;//this->getdYTransform(ZETA,p) ;
// 		double zdzeta = 0 ;//this->getdZTransform(ZETA,p) ;
// 		VirtualMachine vm ;
// 		TetrahedralElement * father = nullptr ;
// 
// 		std::valarray<Function> * functions = shapefunc ;
// 		if(!shapefunc)
// 		{
// 			father = new TetrahedralElement(order) ;
// 			functions = father->shapefunc ;
// 		}
// 		
// 		for(size_t i = 0 ; i < getBoundingPoints().size() ; i++)
// 		{
// // 			std::cout << i << "  "<< shapefunc << std::endl ;
// 			double dxi = vm.deval((*functions)[i], XI, p) ;
// 			double deta = vm.deval((*functions)[i], ETA, p) ;
// 			double dzeta = vm.deval((*functions)[i], ZETA, p) ;
// 			xdxi += dxi*getBoundingPoint(i).getX() ;
// 			ydxi += dxi*getBoundingPoint(i).getY() ;
// 			zdxi += dxi*getBoundingPoint(i).getZ() ;
// 
// 			xdeta += deta*getBoundingPoint(i).getX() ;
// 			ydeta += deta*getBoundingPoint(i).getY() ;
// 			zdeta += deta*getBoundingPoint(i).getZ() ;
// 
// 			xdzeta += dzeta*getBoundingPoint(i).getX() ;
// 			ydzeta += dzeta*getBoundingPoint(i).getY() ;
// 			zdzeta += dzeta*getBoundingPoint(i).getZ() ;
// 		}
// 
// 		ret[0][0] = xdxi ; ret[0][1] = ydxi ; ret[0][2] = zdxi ; 
// 		ret[1][0] = xdeta ; ret[1][1] = ydeta ; ret[1][2] = zdeta ;
// 		ret[2][0] = xdzeta ; ret[2][1] = ydzeta ; ret[2][2] = zdzeta ;
// 		invert3x3Matrix(ret) ;
// 		delete father ;
// 	}
// 	else
// 	{
// 		if(ret.numRows() != 4)
// 			ret.resize(4,4) ;
// 		
// // 		double xdxi = this->getdXTransform(XI,p) ;
// // 		double ydxi = this->getdYTransform(XI,p) ;
// // 		double zdxi = this->getdZTransform(XI,p) ;
// // // 		double tdxi = this->getdTTransform(XI,p) ;
// // 		
// // 		double xdeta = this->getdXTransform(ETA,p) ;
// // 		double ydeta = this->getdYTransform(ETA,p) ;
// // 		double zdeta = this->getdZTransform(ETA,p) ;
// // // 		double tdeta = this->getdTTransform(ETA,p) ;
// // 		
// // 		double xdzeta = this->getdXTransform(ZETA,p) ;
// // 		double ydzeta = this->getdYTransform(ZETA,p) ;
// // 		double zdzeta = this->getdZTransform(ZETA,p) ;
// // // 		double tdzeta = this->getdTTransform(ZETA,p) ;
// // 		
// // // 		double xdtheta = this->getdXTransform(TIME_VARIABLE,p) ;
// // // 		double ydtheta = this->getdYTransform(TIME_VARIABLE,p) ;
// // // 		double zdtheta = this->getdZTransform(TIME_VARIABLE,p) ;
// // 		double tdtheta = this->getdTTransform(TIME_VARIABLE,p) ;
// 
// 		double xdxi = 0 ;//this->getdXTransform(XI,p) ;
// 		double ydxi = 0 ;//this->getdYTransform(XI,p) ;
// 		double zdxi = 0 ;//this->getdZTransform(XI,p) ;
// 		double tdxi = 0 ;
// 		
// 		double xdeta = 0 ;//this->getdXTransform(ETA,p) ;
// 		double ydeta = 0 ;//this->getdYTransform(ETA,p) ;
// 		double zdeta = 0 ;//this->getdZTransform(ETA,p) ;
// 		double tdeta = 0 ;
// 		
// 		double xdzeta = 0 ;//this->getdXTransform(ZETA,p) ;
// 		double ydzeta = 0 ;//this->getdYTransform(ZETA,p) ;
// 		double zdzeta = 0 ;//this->getdZTransform(ZETA,p) ;
// 		double tdzeta = 0 ;
// 
// 		double xdtau = 0 ;//this->getdXTransform(ZETA,p) ;
// 		double ydtau = 0 ;//this->getdYTransform(ZETA,p) ;
// 		double zdtau = 0 ;//this->getdZTransform(ZETA,p) ;
// 		double tdtau = 0 ;
// 		
// 		VirtualMachine vm ;
// 		for(size_t i = 0 ; i < getBoundingPoints().size() ; i++)
// 		{
// 			double dxi = vm.deval(getShapeFunction(i), XI, p) ;
// 			double deta = vm.deval(getShapeFunction(i), ETA, p) ;
// 			double dzeta = vm.deval(getShapeFunction(i), ZETA, p) ;
// 			double dtau = vm.deval(getShapeFunction(i), TIME_VARIABLE, p) ;
// 			
// 			xdxi += dxi*getBoundingPoint(i).getX() ;
// 			ydxi += dxi*getBoundingPoint(i).getY() ;
// 			zdxi += dxi*getBoundingPoint(i).getZ() ;
// 			tdxi += dxi*getBoundingPoint(i).getT() ;
// 
// 			xdeta += deta*getBoundingPoint(i).getX() ;
// 			ydeta += deta*getBoundingPoint(i).getY() ;
// 			zdeta += deta*getBoundingPoint(i).getZ() ;
// 			tdeta += deta*getBoundingPoint(i).getT() ;
// 
// 			xdzeta += dzeta*getBoundingPoint(i).getX() ;
// 			ydzeta += dzeta*getBoundingPoint(i).getY() ;
// 			zdzeta += dzeta*getBoundingPoint(i).getZ() ;
// 			tdzeta += dzeta*getBoundingPoint(i).getT() ;
// 			
// 			xdtau += dtau*getBoundingPoint(i).getX() ;
// 			ydtau += dtau*getBoundingPoint(i).getY() ;
// 			zdtau += dtau*getBoundingPoint(i).getZ() ;
// 			tdtau += dtau*getBoundingPoint(i).getT() ;
// 			
// 		}
// 
// 		ret[0][0] = xdxi ; ret[0][1] = ydxi ; ret[0][2] = zdxi ; ret[0][3] = 0; 
// 		ret[1][0] = xdeta ; ret[1][1] = ydeta ; ret[1][2] = zdeta ; ret[1][3] = 0;
// 		ret[2][0] = xdzeta ; ret[2][1] = ydzeta ; ret[2][2] = zdzeta ; ret[2][3] = 0;
// 		ret[3][0] = 0 ; ret[3][1] = 0 ; ret[3][2] = 0 ; ret[3][3] = tdtau;
// 
// 		ret = inverse4x4Matrix(ret) ;
// 	}

}

 Vector HexahedralElement::getNonLinearForces() 
{
	return Vector(0);
}

 std::valarray<std::valarray<Amie::Matrix> >  HexahedralElement::getNonLinearElementaryMatrix()  
{
	return std::valarray<std::valarray<Amie::Matrix> >(0);
}

void ElementaryVolume::setBehaviour( Mesh< DelaunayTetrahedron, DelaunayTreeItem3D >* msh, Form* f )
{	
	bool init = false ;
	if(state)
	{	
 		init = this->getState().getDisplacements().size() > 0 ;
		delete state ;
	}
	state = f->createElementState( this ) ;
//	Form * old = behaviour ;
	behaviour = f ;
	if(init)
	{
		state->initialize(msh) ;
	}
}

Order ElementaryVolume::getOrder() const
{
	return order ;
}

 void ElementaryVolume::setOrder(Order o ) 
{

	order = o;
	
	switch(order)
	{
	case CONSTANT:
		{
			break ;
		}
	case LINEAR :
		{
			break ;
		}
	case QUADRATIC :
		{
			break ;
		}
	case CUBIC :
		{
			break ;
		}
	case QUADRIC :
		{
			break ;
		}
	case QUINTIC :
		{
			break ;
		}
	case QUADTREE_REFINED :
		{
			break ;
		}
	case CONSTANT_TIME_LINEAR :
		{
			timePlanes() = 2 ;
			break ;
		}
	case CONSTANT_TIME_QUADRATIC :
		{
			timePlanes() = 3 ;
			break ;
		}
	case LINEAR_TIME_LINEAR :
		{
			timePlanes() = 2 ;
			break ;
		}
	case LINEAR_TIME_QUADRATIC :
		{
			timePlanes() = 3 ;
			break ;
		}
	case QUADRATIC_TIME_LINEAR :
		{
			timePlanes() = 2 ;
			break ;
		}
	case QUADRATIC_TIME_QUADRATIC :
		{
			timePlanes() = 3 ;
			break ;
		}
	case CUBIC_TIME_LINEAR :
		{
			timePlanes() = 3 ;
			break ;
		}
	case CUBIC_TIME_QUADRATIC	 :
		{
			timePlanes() = 3 ;
			break ;
		}
	case QUADRIC_TIME_LINEAR :
		{
			timePlanes() = 2 ;
			break ;
		}
	case QUADRIC_TIME_QUADRATIC :
		{
			timePlanes() = 3 ;
			break ;
		}
	case QUINTIC_TIME_LINEAR :
		{
			timePlanes() = 2 ;
			break ;
		}
	case QUINTIC_TIME_QUADRATIC :
		{
			timePlanes() = 3 ;
			break ;
		}
	default:
		{
			break ;
		}
	}
}

Function zeroes(int numberOfDerivatives, int derivationDepth)
{
	Function ret("0") ;
	if(derivationDepth > 0)
	{
		ret.setNumberOfDerivatives( numberOfDerivatives ) ;
		Function d = zeroes( numberOfDerivatives, derivationDepth-1 ) ;
		for(size_t i = 0 ; i < numberOfDerivatives ; i++)
		{
			ret.setDerivative( (const Variable) i, d ) ;
		}
	}
	return ret ;
}

Function XTransform(const std::valarray<Amie::Point*> & points, const std::valarray<Function > & basis)
{
//	Function ret = zeroes( 4, basis[0].derivationDepth()-1 ) ;
 	Function ret("0") ;
//  	Function zero("0") ;
//  	ret.setNumberOfDerivatives(4);
//  	ret.setDerivative(XI, zero);
//  	ret.setDerivative(ETA, zero);
//  	ret.setDerivative(ZETA, zero);
//  	ret.setDerivative(TIME_VARIABLE, zero);
	
	for(size_t i = 0 ; i < points.size() ; i++)
	{
		ret += basis[i]*points[i]->getX() ;
	}
	return ret ;
}


double xTransform(const Point & p, const std::valarray<Amie::Point*> & points, const std::valarray<Function > & basis)
{
	double ret = 0;
	
	VirtualMachine vm ;
	for(size_t i = 0 ; i < points.size() ; i++)
	{
		ret += vm.eval(basis[i], p)*points[i]->getX() ;
	}
	return ret ;
}

Function YTransform(const std::valarray<Amie::Point*> & points, const std::valarray<Function > & basis)
{
	Function ret("0") ;
	Function zero("0") ;
	
// 	ret.setNumberOfDerivatives(4);
// 	ret.setDerivative(XI, zero);
// 	ret.setDerivative(ETA, zero);
// 	ret.setDerivative(ZETA, zero);
// 	ret.setDerivative(TIME_VARIABLE, zero);
	
	assert(points.size() == basis.size()) ;
	for(size_t i = 0 ; i < points.size() ; i++)
	{
		ret += basis[i]*points[i]->getY() ;
	}
	
	return ret ;
}


double yTransform(const Point & p, const std::valarray<Amie::Point*> & points, const std::valarray<Function > & basis)
{
	double ret = 0;
	
	VirtualMachine vm ;
	for(size_t i = 0 ; i < points.size() ; i++)
	{
		ret += vm.eval(basis[i], p)*points[i]->getY() ;
	}
	return ret ;
}

Function ZTransform(const std::valarray<Amie::Point*> & points, const std::valarray<Function > & basis)
{
	Function ret("0") ;
	Function zero("0") ;
// 	ret.setNumberOfDerivatives(4);
// 	ret.setDerivative(XI, zero);
// 	ret.setDerivative(ETA, zero);
// 	ret.setDerivative(ZETA, zero);
// 	ret.setDerivative(TIME_VARIABLE, zero);

	assert(points.size() == basis.size()) ;

	for(size_t i = 0 ; i < points.size() ; i++)
	{
		ret += basis[i]*points[i]->getZ() ;

	}

	return ret ;
}

double zTransform(const Point & p, const std::valarray<Amie::Point*> & points, const std::valarray<Function > & basis)
{
	double ret = 0;
	
	VirtualMachine vm ;
	for(size_t i = 0 ; i < points.size() ; i++)
	{
		ret += vm.eval(basis[i], p)*points[i]->getZ() ;
	}
	return ret ;
}

Function TTransform(const std::valarray<Amie::Point*> & points, const std::valarray<Function > & basis)
{
	Function ret("0") ;
	Function zero("0") ;
// 	ret.setNumberOfDerivatives(4);
// 	ret.setDerivative(XI, zero);
// 	ret.setDerivative(ETA, zero);
// 	ret.setDerivative(ZETA, zero);
// 	ret.setDerivative(TIME_VARIABLE, zero);
	
	assert(points.size() == basis.size()) ;
	for(size_t i = 0 ; i < points.size() ; i++)
	{
		ret += basis[i]*points[i]->getT() ;
		
	}

	return ret ;
}

double tTransform(const Point & p, const std::valarray<Amie::Point*> & points, const std::valarray<Function > & basis)
{
	double ret = 0;
	
	VirtualMachine vm ;
	for(size_t i = 0 ; i < points.size() ; i++)
	{
		ret += vm.eval(basis[i], p)*points[i]->getT() ;
	}
	return ret ;
}

Point coordinateTransform(const Point & p, const std::valarray<Amie::Point*> & points, const std::valarray<Function > & basis)
{
	Point ret;
	
	VirtualMachine vm ;
	for(size_t i = 0 ; i < points.size() ; i++)
	{
		double v = vm.eval(basis[i], p) ;
		ret.getX() += v*points[i]->getX() ;
		ret.getY() += v*points[i]->getY() ;
		ret.getZ() += v*points[i]->getZ() ;
		ret.getT() += v*points[i]->getT() ;
	}
	return ret ;
}

Amie::Function dXTransform(const std::valarray<Amie::Point*> & points ,const std::valarray< Amie::Function> &basis, Amie::Variable v)
{
	switch(v)
	{
	case XI:
		{
			Function der_x ;
			
			assert(points.size() == basis.size()) ;
			for(size_t i = 0 ; i < points.size() ; i++)
			{
				der_x += basis[i].d(XI)*points[i]->getX() ;
			}

			return der_x ;
		}
	case ETA:
		{
			Function der_y ;
			
			assert(points.size() == basis.size()) ;
			for(size_t i = 0 ; i < points.size() ; i++)
			{
				der_y += basis[i].d(ETA)*points[i]->getX() ;
			}
			
			return der_y ;
		}
	case ZETA:
		{
			Function der_z ;
			
			assert(points.size() == basis.size()) ;
			for(size_t i = 0 ; i < points.size() ; i++)
			{
				der_z += basis[i].d(ZETA)*points[i]->getX() ;
			}
			
			return der_z ;
		}
	case TIME_VARIABLE:
		{
			Function der_t ;
			
			assert(points.size() == basis.size()) ;
			for(size_t i = 0 ; i < points.size() ; i++)
			{
				der_t += basis[i].d(TIME_VARIABLE)*points[i]->getX() ;
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

Amie::Function dYTransform(const std::valarray<Amie::Point*> & points ,const std::valarray< Amie::Function> &basis, Amie::Variable v)
{
	switch(v)
	{
	case XI:
		{
			Function der_x ;
			
			assert(points.size() == basis.size()) ;
			for(size_t i = 0 ; i < points.size() ; i++)
			{
				der_x += basis[i].d(XI)*points[i]->getY() ;
			}
			
			return der_x ;
		}
	case ETA:
		{
			Function der_y ;
			
			assert(points.size() == basis.size()) ;
			for(size_t i = 0 ; i < points.size() ; i++)
			{
				der_y += basis[i].d(ETA)*points[i]->getY() ;
			}
			
			return der_y ;
		}
	case ZETA:
		{
			Function der_z ;
			
			assert(points.size() == basis.size()) ;
			for(size_t i = 0 ; i < points.size() ; i++)
			{
				der_z += basis[i].d(ZETA)*points[i]->getY() ;
			}
			
			return der_z ;
		}
	case TIME_VARIABLE:
		{
			Function der_t ;
			
			assert(points.size() == basis.size()) ;
			for(size_t i = 0 ; i < points.size() ; i++)
			{
				der_t += basis[i].d(TIME_VARIABLE)*points[i]->getY() ;
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

Amie::Function dZTransform(const std::valarray<Amie::Point*> & points ,const std::valarray< Amie::Function> &basis, Amie::Variable v)
{
	switch(v)
	{
	case XI:
		{
			Function der_x ;
			
			assert(points.size() == basis.size()) ;
			for(size_t i = 0 ; i < points.size() ; i++)
			{
				der_x += basis[i].d(XI)*points[i]->getZ() ;
			}
			
			return der_x ;
		}
	case ETA:
		{
			Function der_y ;
			
			assert(points.size() == basis.size()) ;
			for(size_t i = 0 ; i < points.size() ; i++)
			{
				der_y += basis[i].d(ETA)*points[i]->getZ() ;
			}
			
			return der_y ;
		}
	case ZETA:
		{
			Function der_z ;
			
			assert(points.size() == basis.size()) ;
			for(size_t i = 0 ; i < points.size() ; i++)
			{
				der_z += basis[i].d(ZETA)*points[i]->getZ() ;
			}
			
			return der_z ;
		}
	case TIME_VARIABLE:
		{
			Function der_t ;
			
			assert(points.size() == basis.size()) ;
			for(size_t i = 0 ; i < points.size() ; i++)
			{
				der_t += basis[i].d(TIME_VARIABLE)*points[i]->getT() ;
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

Amie::Function dTTransform(const std::valarray<Amie::Point*> & points ,const std::valarray< Amie::Function> &basis, Amie::Variable v)
{
	switch(v)
	{
	case XI:
		{
			Function der_x ;
			
			assert(points.size() == basis.size()) ;
			for(size_t i = 0 ; i < points.size() ; i++)
			{
				der_x += basis[i].d(XI)*points[i]->getT() ;
			}
			
			return der_x ;
		}
	case ETA:
		{
			Function der_y ;
			
			assert(points.size() == basis.size()) ;
			for(size_t i = 0 ; i < points.size() ; i++)
			{
				der_y += basis[i].d(ETA)*points[i]->getT() ;
			}
			
			return der_y ;
		}
	case ZETA:
		{
			Function der_z ;
			
			assert(points.size() == basis.size()) ;
			for(size_t i = 0 ; i < points.size() ; i++)
			{
				der_z += basis[i].d(ZETA)*points[i]->getT() ;
			}
			
			return der_z ;
		}
	case TIME_VARIABLE:
		{
			Function der_t ;
			
			assert(points.size() == basis.size()) ;
			for(size_t i = 0 ; i < points.size() ; i++)
			{
				der_t += basis[i].d(ZETA)*points[i]->getT() ;
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

double dXTransform(const std::valarray<Amie::Point*> & points ,const std::valarray< Amie::Function> &basis, Amie::Variable v, const Point & p)
{
	switch(v)
	{
	case XI:
		{
			VirtualMachine vm ;
			double der_x = 0;
			
			for(size_t i = 0 ; i < basis.size() ; i++)
			{
				der_x += vm.deval(basis[i],XI,p)*points[i]->getX() ;
			}
			return der_x ;
		}
	case ETA:
		{
			VirtualMachine vm ;
			double der_y = 0 ;

			for(size_t i = 0 ; i < basis.size() ; i++)
			{
				der_y += vm.deval(basis[i],ETA, p)*points[i]->getX() ;
			}
			return der_y ;
		}
	case ZETA:
		{
			VirtualMachine vm ;
			double der_z = 0 ;
			
			for(size_t i = 0 ; i < basis.size() ; i++)
			{
				der_z += vm.deval(basis[i],ZETA,p)*points[i]->getX() ;
			}
			
			return der_z ;
		}
	case TIME_VARIABLE:
		{
			VirtualMachine vm ;
			double der_t = 0 ;
			
			for(size_t i = 0 ; i < basis.size() ; i++)
			{
				der_t += vm.deval(basis[i],TIME_VARIABLE,p)*points[i]->getX() ;
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

double dYTransform(const std::valarray<Amie::Point*> & points ,const std::valarray< Amie::Function> &basis, Amie::Variable v, const Point & p)
{
	switch(v)
	{
	case XI:
		{
			
			VirtualMachine vm ;
			double der_x = 0;

			for(size_t i = 0 ; i < basis.size() ; i++)
			{
				der_x += vm.deval(basis[i],XI,p)*points[i]->getY() ;
			}
			return der_x ;
		}
	case ETA:
		{
			VirtualMachine vm ;
			double der_y = 0 ;
			
			for(size_t i = 0 ; i < basis.size() ; i++)
			{
				der_y += vm.deval(basis[i],ETA, p)*points[i]->getY() ;
			}
			
			return der_y ;
		}
	case ZETA:
		{
			VirtualMachine vm ;
			double der_z = 0 ;
			
			for(size_t i = 0 ; i < basis.size() ; i++)
			{
				der_z += vm.deval(basis[i],ZETA,p)*points[i]->getY() ;
			}
			
			return der_z ;
		}
	case TIME_VARIABLE:
		{
			VirtualMachine vm ;
			double der_t = 0 ;
			

			for(size_t i = 0 ; i < basis.size() ; i++)
			{
				der_t += vm.deval(basis[i],TIME_VARIABLE,p)*points[i]->getY() ;
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

double dZTransform(const std::valarray<Amie::Point*> & points ,const std::valarray< Amie::Function> &basis, Amie::Variable v, const Point & p)
{
	switch(v)
	{
	case XI:
		{
			VirtualMachine vm ;
			double der_x = 0;
			
			for(size_t i = 0 ; i < points.size() ; i++)
			{
				der_x += vm.deval(basis[i],XI,p)*points[i]->getZ() ;
			}
			
			return der_x ;
		}
	case ETA:
		{
			VirtualMachine vm ;
			double der_y = 0 ;
			
			for(size_t i = 0 ; i < basis.size() ; i++)
			{
				der_y += vm.deval(basis[i],ETA, p)*points[i]->getZ() ;
			}
			
			return der_y ;
		}
	case ZETA:
		{
			VirtualMachine vm ;
			double der_z = 0 ;
			
			for(size_t i = 0 ; i < points.size() ; i++)
			{
				der_z += vm.deval(basis[i],ZETA,p)*points[i]->getZ() ;
			}
			
			return der_z ;
		}
	case TIME_VARIABLE:
		{
			VirtualMachine vm ;
			double der_t = 0 ;
			
			for(size_t i = 0 ; i < points.size() ; i++)
			{
				der_t += vm.deval(basis[i],TIME_VARIABLE,p)*points[i]->getZ() ;
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

double dTTransform(const std::valarray<Amie::Point*> & points ,const std::valarray< Amie::Function> &basis, Amie::Variable v, const Point & p)
{
	switch(v)
	{
	case XI:
		{
			VirtualMachine vm ;
			double der_x = 0;
			for(size_t i = 0 ; i < basis.size() ; i++)
			{
				der_x += vm.deval(basis[i],XI,p)*points[i]->getT() ;
			}
			return der_x ;
		}
	case ETA:
		{
			VirtualMachine vm ;
			double der_y = 0 ;
			
			for(size_t i = 0 ; i < points.size() ; i++)
			{
				der_y += vm.deval(basis[i],ETA, p)*points[i]->getT() ;
			}
			
			return der_y ;
		}
	case ZETA:
		{
			VirtualMachine vm ;
			double der_z = 0 ;
			
			for(size_t i = 0 ; i < points.size() ; i++)
			{
				der_z += vm.deval(basis[i],ZETA,p)*points[i]->getT() ;
			}
			
			return der_z ;
		}
	case TIME_VARIABLE:
		{
			VirtualMachine vm ;
			double der_t = 0 ;
// 			std::cout << "-----" << std::endl ;
// 			std::cout << der_t << std::endl ;
			for(size_t i = 0 ; i < points.size() ; i++)
			{
// 				vm.print(basis[i]) ;
				der_t += vm.deval(basis[i],TIME_VARIABLE,p)*points[i]->getT() ;
// 				std::cout << der_t << "   " << points[i]->getT() << std::endl ;
			}
// 			std::cout << "-----" << std::endl ;
			return der_t ;
		}
	default:
		{
			std::cout << "*** dTTransform error ***" << std::endl ;
			return 0 ;
		}
	}
}

const GaussPointArray & HexahedralElement::genGaussPoints()
{
	if(getCachedGaussPoints())
		return *getCachedGaussPoints() ;
	
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
		fin[0] = std::pair<Point, double>(Point(0, 0, 0), 8.0) ;
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
// 	Function  J = this->jacobian() ;
// 	
// 	VirtualMachine vm ;
	for(size_t i = 0 ; i < fin.size() ; i++)
	{
		fin[i].second*= jacobianAtPoint(fin[i].first) ;
	}
	
	setCachedGaussPoints(new GaussPointArray(fin, order)) ;
	return *getCachedGaussPoints() ;
}
	
HexahedralElement::HexahedralElement(Order order, bool f ) : ElementaryVolume( f)
{
	setOrder(order) ;
	visited = false ;

	if(order == LINEAR)
	{
		assert(false) ;
// 		this->Hexahedron::sampleSurface(8) ;
// 		shapefunc = new std::valarray<Function>(8) ;
// 		Matrix zero ;
// 		std::valarray<Matrix> f1(zero, 2) ;
// 		f1[1][0][0] = -0.125 ;
// 		f1[0][1][0] = -0.125 ;
// 		f1[0][0][1] = -0.125 ;
// 		f1[0][0][0] = 0.125 ;
// 		f1[1][1][0] = 0.125 ;
// 		f1[0][1][1] = 0.125 ;
// 		f1[1][0][1] = 0.125 ;
// 		f1[1][1][1] =-0.125 ;
// 		
// 		std::valarray<Matrix> f2(zero,2) ;
// 		f2[1][0][0] = 0.125 ;
// 		f2[0][1][0] = -0.125 ;
// 		f2[0][0][1] = -0.125 ;
// 		f2[0][0][0] = 0.125 ;
// 		f2[1][1][0] = -0.125 ;
// 		f2[0][1][1] = 0.125 ;
// 		f2[1][0][1] = -0.125 ;
// 		f2[1][1][1] = 0.125 ;
// 	
// 		std::valarray<Matrix> f3(zero,2) ;
// 		f3[1][0][0] = 0.125 ;
// 		f3[0][1][0] = 0.125 ;
// 		f3[0][0][1] = -0.125 ;
// 		f3[0][0][0] = 0.125 ;
// 		f3[1][1][0] = 0.125 ;
// 		f3[0][1][1] = -0.125 ;
// 		f3[1][0][1] = -0.125 ;
// 		f3[1][1][1] = -0.125 ;
// 		
// 		std::valarray<Matrix> f4(zero,2) ;
// 		f4[1][0][0] = -0.125 ;
// 		f4[0][1][0] = 0.125 ;
// 		f4[0][0][1] = -0.125 ;
// 		f4[0][0][0] = 0.125 ;
// 		f4[1][1][0] = -0.125 ;
// 		f4[0][1][1] = -0.125 ;
// 		f4[1][0][1] = 0.125 ;
// 		f4[1][1][1] = 0.125 ;
// 			
// 		std::valarray<Matrix> f5(zero,2) ;
// 		f5[1][0][0] = -0.125 ;
// 		f5[0][1][0] = -0.125 ;
// 		f5[0][0][1] = 0.125 ;
// 		f5[0][0][0] = 0.125 ;
// 		f5[1][1][0] = 0.125 ;
// 		f5[0][1][1] = -0.125 ;
// 		f5[1][0][1] = -0.125 ;
// 		f5[1][1][1] = 0.125 ;
// 		
// 		std::valarray<Matrix> f6(zero,2) ;
// 		f6[1][0][0] = 0.125 ;
// 		f6[0][1][0] = -0.125 ;
// 		f6[0][0][1] = 0.125 ;
// 		f6[0][0][0] = 0.125 ;
// 		f6[1][1][0] = -0.125 ;
// 		f6[0][1][1] = -0.125 ;
// 		f6[1][0][1] = 0.125 ;
// 		f6[1][1][1] = -0.125 ;
// 			
// 		std::valarray<Matrix> f7(zero,2) ;
// 		f7[1][0][0] = 0.125 ;
// 		f7[0][1][0] = 0.125 ;
// 		f7[0][0][1] = 0.125 ;
// 		f7[0][0][0] = 0.125 ;
// 		f7[1][1][0] = 0.125 ;
// 		f7[0][1][1] = 0.125 ;
// 		f7[1][0][1] = 0.125 ;
// 		f7[1][1][1] = 0.125 ;
// 		
// 		std::valarray<Matrix> f8(zero,2) ;
// 		f8[1][0][0] = -0.125 ;
// 		f8[0][1][0] = 0.125 ;
// 		f8[0][0][1] = 0.125 ;
// 		f8[0][0][0] = 0.125 ;
// 		f8[1][1][0] = -0.125 ;
// 		f8[0][1][1] = 0.125 ;
// 		f8[1][0][1] = -0.125 ;
// 		f8[1][1][1] = -0.125 ;
// 			
// 			//0
// 		(*shapefunc)[0] = Function(f1) ;
// 			//1
// 		(*shapefunc)[1] = Function(f5) ;
// 			//2
// 		(*shapefunc)[2] = Function(f4) ;
// 			//3
// 		(*shapefunc)[3] = Function(f8) ;
// 			
// 			//4
// 		(*shapefunc)[4] = Function(f2) ;
// 			//5
// 		(*shapefunc)[5] = Function(f6) ;
// 			//6
// 		(*shapefunc)[6] = Function(f3) ;
// 			//7
// 		(*shapefunc)[7] = Function(f7) ;
// 		
		
	}
	else if(order == QUADRATIC)
	{
		assert(false) ;
	}
	else
	{
		assert(false) ;
	}
	assert(this->jacobianAtPoint(Hexahedron::getCenter()) > 0) ;
	
	this->size_x  = 2 ;
	this->size_y  = 2 ;
	this->size_z  = 2 ;
}
	
HexahedralElement::HexahedralElement(HexahedralElement * parent,Hexahedron * t)
{
	setOrder( parent->getOrder() );
	visited = false ;

	this->shapefunc = parent->shapefunc ;

	for(size_t i =  0  ; i < this->size() ; i++)
	{
		delete &this->getPoint(i) ;
	}
	
	this->Hexahedron::center = t->getCenter() ;
	
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
	
	setOrder(parent->getOrder()) ;
	this->shapefunc = parent->shapefunc ;
}

void HexahedralElement::print()  const
{
	std::cout << "I am an hexahedron..." << std::endl ;
}

std::vector<int> validNeighbours(int i, int j, int k, int length)
{
	std::vector<int> ret ;
	
	if(k > 0)
		ret.push_back(i + j*length + (k-1)*length*length) ;
	if(k < length-1)
		ret.push_back(i + j*length + (k+1)*length*length) ;
	if(i > 0)
		ret.push_back(i-1 + j*length + (k)*length*length) ;
	if(i < length-1)
		ret.push_back(i+1 + j*length + (k)*length*length) ;
	if(j > 0)
		ret.push_back(i + (j-1)*length + (k)*length*length) ;
	if(j < length-1)
		ret.push_back(i + (j+1)*length + (k)*length*length) ;
	
// 	for(int ii = i-1 ; ii < i+2 ; ii++)
// 	{
// 		for(int jj = j-1 ; jj < j+2 ; jj++)
// 		{
// 			for(int kk = k-1 ; kk < k+2 ; kk++)
// 			{
// 				if((ii >= 0 && ii < length)
// 				   && (jj >= 0 && jj < length)
// 				   && (kk >= 0 && kk < length)
// 				   && !( ii == i && jj== j && kk==k)
// 				   && ( ii == i || jj== j || kk==k)
// 				  )
// 					ret.push_back(kk + jj*length + ii*length*length) ;
// 			}
// 		}
// 	}
	return ret ;
}

void computeNeighbourhoodForStructuredHexahedralMesh(std::vector<Amie::HexahedralElement *> & vec)
{
	int length = (int)round(std::pow(vec.size(), 1./3.)) ;


	for(int i = 0 ; i < length ; i++)
	{
		for(int j = 0 ; j < length ; j++)
		{
			for(int k = 0 ; k < length ; k++)
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

void burn(std::vector<Amie::HexahedralElement *> & vec)
{
	int length = (int)round(std::pow(vec.size(), 1./3.)) ;

	
	std::vector<Amie::HexahedralElement *> to_check ;
	std::vector<Amie::HexahedralElement *> connected ;
	for(int i = 0 ; i < length ; i++)
	{
		for(int j = 0 ; j < length ; j++)
		{

			int current_index = i + j*length ;
			
			if(vec[current_index]->getBehaviour()->type != VOID_BEHAVIOUR)	
			{
				to_check.push_back(vec[current_index]) ;
			}
				
		}
	}
	
	while(!to_check.empty())
	{
		std::vector<Amie::HexahedralElement *> temp ;
		
		for(size_t i = 0 ; i< to_check.size() ; i++)
		{
			to_check[i]->visited = true ;
			connected.push_back(to_check[i]) ;
			
			for(size_t l = 0 ; l < to_check[i]->neighbourhood.size() ; l++)
			{
				if(!to_check[i]->neighbourhood[l]->visited 
				   && to_check[i]->neighbourhood[l]->getBehaviour()->type != VOID_BEHAVIOUR
				  )
					temp.push_back(to_check[i]->neighbourhood[l]) ;
			}
			
		}
		
		std::stable_sort(temp.begin(), temp.end()) ;
		auto e = std::unique(temp.begin(), temp.end()) ;
		temp.erase(e, temp.end()) ;
		
		to_check = temp ;
	}
	
	std::stable_sort(connected.begin(), connected.end()) ;
	auto e = std::unique(connected.begin(), connected.end()) ;
	connected.erase(e, connected.end()) ;
	
	for(size_t i = 0 ; i < vec.size() ; i++)
		vec[i]->visited = false ;
	for(size_t i = 0 ; i < connected.size() ; i++)
		connected[i]->visited = true ;

}

}