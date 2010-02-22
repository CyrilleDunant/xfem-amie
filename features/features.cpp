// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "features.h"
#include "crackinitiation.h"
#include "layeredinclusion.h"
#include "../physics/void_form.h"
#ifdef HAVE_OPENMP
#include <omp.h>
#endif


using namespace Mu ;


void Feature::setBoundary(Geometry * g)
{
	delete boundary ;
	this->boundary = g ;
}

Geometry * Feature::getBoundary()
{
	return this->boundary ;
}

const Geometry * Feature::getBoundary() const
{
	return this->boundary ;
}

Feature::Feature(Feature * father)
{
	this->isEnrichmentFeature = false ;
	this->isCompositeFeature = false ;
	this->isVirtualFeature = false ;
		
	this->behaviour = new VoidForm() ;

	
	m_f = father ;
	if(father != NULL)
		father->addChild(this) ;
	infRad = DEFAULT_BOUNDARY ;
	this->boundary = NULL ;
	this->boundary2 = NULL ;
}

Feature::Feature()
{
	this->isEnrichmentFeature = false ;
	this->isCompositeFeature = false ;
	this->isVirtualFeature = false ;

	this->behaviour = new VoidForm() ;
	
	m_f = NULL ;
	infRad = DEFAULT_BOUNDARY ;
	this->boundary = NULL ;
	this->boundary2 = NULL ;
}

Feature::Feature(Feature *father, Geometry * b) 
{ 
	boundary = b ; 
	m_f = father ; 
	this->isEnrichmentFeature = false ;
	this->isCompositeFeature = false ;
	this->isVirtualFeature = false ;

	this->behaviour = new VoidForm() ;
	
	if(father != NULL)
		father->addChild(this) ;
	
	infRad = DEFAULT_BOUNDARY ;
	this->boundary2 = NULL ;
}

 bool Feature::inBoundary(const Point & v) const 
{
	bool ret(false) ;
	if(boundary)
		ret = getBoundary()->in(v) ;

	for(size_t i = 0 ;  i < this->m_c.size() ; i++)
		ret = ret || m_c[i]->inBoundary(v) ;
	
	return ret ;
}


std::vector<Point *> Feature::doubleSurfaceSampling()
{
	std::vector<Point *> ret ;
	std::valarray<Point *> newboundingPoints(this->getBoundingPoints().size()*2) ;
	for(size_t i = 0 ; i < this->getBoundingPoints().size()-1 ; i++)
	{
		newboundingPoints[i*2] = &getBoundingPoint(i) ;
		newboundingPoints[i*2+1] = new Point(getBoundingPoint(i)*0.5 + getBoundingPoint(i+1)*0.5) ;
		this->project(newboundingPoints[i*2+1]) ;
		newboundingPoints[i*2+1]->id = -1 ;
		ret.push_back(newboundingPoints[i*2+1]) ;
	}
	
	newboundingPoints[(getBoundingPoints().size()-1)*2] = &getBoundingPoint(getBoundingPoints().size()-1) ;
	newboundingPoints[getBoundingPoints().size()*2-1] = new Point(getBoundingPoint(getBoundingPoints().size()-1)*0.5 + getBoundingPoint(0)*0.5) ;
	newboundingPoints[getBoundingPoints().size()*2-1]->id = -1 ;
	this->project(newboundingPoints[getBoundingPoints().size()*2-1]) ;
	ret.push_back(newboundingPoints[getBoundingPoints().size()*2-1]) ;
	
	dynamic_cast<Geometry *>(this)->setBoundingPoints(newboundingPoints) ;
	
	return ret ;
}

bool Feature::inBoundary(const Point *v) const
{
	bool ret = boundary->in(*v) ;
	for(size_t i = 0 ;  i < this->m_c.size() ; i++)
		ret =  ret || m_c[i]->inBoundary(v) ;
	
	return ret ;
}

bool Feature::inBoundaryLayer(const Point *v) const
{

	if(boundary2)
		return  boundary->in(*v) && !(boundary2->in(*v)) ;
	else
		return boundary->in(*v) ;

}

void  Feature::addChild(Feature *f)
{
	if(std::find(m_c.begin(), m_c.end(), f) == m_c.end())
	{
		m_c.push_back(f) ;
	}
}

void  Feature::removeChild(Feature *f)
{
	std::vector<Feature *>::iterator c = std::find(m_c.begin(), m_c.end(), f) ;
	if(c != m_c.end())
		m_c.erase(c) ;
}

void Feature::setBehaviour(Form * f)
{
	delete this->behaviour ;
	this->behaviour = f ;
}



BoundingBoxDefinedBoundaryCondition::BoundingBoxDefinedBoundaryCondition(LagrangeMultiplierType t, BoundingBoxPosition pos, double d) :BoundaryCondition(t, d), pos(pos) {} ;

void apply2DBC(ElementarySurface *e,  std::vector<Point> & id, LagrangeMultiplierType condition, double data, Assembly * a)
{
	if(e->getBehaviour()->type == VOID_BEHAVIOUR)
		return ;
	for(size_t i = 0 ; i < id.size() ; i++)
	{
		switch(condition)
		{
			case GENERAL :
				std::cout << "I don't know how to form a General Lagrange Multiplier from the data" << std::endl ;
				break ;
			case FIX_ALONG_XI:
				a->setPointAlong(XI, 0, id[i].id) ;
				break ;
			case SET_ALONG_XI:
				a->setPointAlong(XI, data, id[i].id) ;
				break ;
			case FIX_ALONG_ETA:
				a->setPointAlong(ETA, 0, id[i].id) ;
				break ;
			case SET_ALONG_ETA:
				a->setPointAlong(ETA, data, id[i].id) ;
				break ;
			case SET_FORCE_XI:
				a->setForceOn(XI, data, id[i].id) ;
				break ;
			case SET_FORCE_ETA:
				a->setForceOn(ETA, data, id[i].id) ;
				break ;
			case SET_STRESS_XI:
			{
				std::vector<Function> shapeFunctions ;
				for(size_t j = 0 ; j < id.size() ; j++)
				{
					for(size_t i = 0 ; i < e->getBoundingPoints().size() ; i++)
					{
						if(id[j].id == e->getBoundingPoint(i).id)
							shapeFunctions.push_back(e->getShapeFunction(i)) ;
					}
				}
				std::vector<Variable> v(2) ;
				v[0] = XI ; v[1] = ETA ;
				Vector imposed(2) ;
				imposed[0] = data ;
				imposed[1] = 0 ;
				imposed[2] = 0 ;
				GaussPointArray gp = e->getGaussPoints() ;
				std::valarray<Matrix> Jinv( Matrix(), gp.gaussPoints.size()) ;
				for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
				{
					e->getInverseJacobianMatrix(gp.gaussPoints[i].first, Jinv[i]) ;
				}
				
				for(size_t i = 0 ; i < shapeFunctions.size() ; ++i)
				{
					Vector forces =  VirtualMachine().ieval(Gradient(shapeFunctions[i]) * (imposed), gp, Jinv, v) ;
					a->addForceOn(XI,forces[0], id[i].id) ;
					a->addForceOn(ETA,forces[1], id[i].id) ;
				}
				return ;
			}
			case SET_STRESS_ETA:
			{
				std::vector<Function> shapeFunctions ;
				for(size_t j = 0 ; j < id.size() ; j++)
				{
					for(size_t i = 0 ; i < e->getBoundingPoints().size() ; i++)
					{
						if(id[j].id == e->getBoundingPoint(i).id)
							shapeFunctions.push_back(e->getShapeFunction(i)) ;
					}
				}
				std::vector<Variable> v(2) ;
				v[0] = XI ; v[1] = ETA ;
				Vector imposed(3) ;
				imposed[0] = 0 ;
				imposed[1] = data ;
				imposed[2] = 0 ;
				GaussPointArray gp = e->getGaussPoints() ;
				std::valarray<Matrix> Jinv( Matrix(), gp.gaussPoints.size()) ;
				for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
				{
					e->getInverseJacobianMatrix(gp.gaussPoints[i].first, Jinv[i]) ;
				}
				
				for(size_t i = 0 ; i < shapeFunctions.size() ; ++i)
				{
					Vector forces =  VirtualMachine().ieval(Gradient(shapeFunctions[i]) * (imposed), gp, Jinv, v) ;
					a->addForceOn(XI,forces[0], id[i].id) ;
					a->addForceOn(ETA,forces[1], id[i].id) ;
				}
				return ;
			}
			case SET_STRESS_XI_ETA:
			{
				std::vector<Function> shapeFunctions ;
				for(size_t j = 0 ; j < id.size() ; j++)
				{
					for(size_t i = 0 ; i < e->getBoundingPoints().size() ; i++)
					{
						if(id[j].id == e->getBoundingPoint(i).id)
							shapeFunctions.push_back(e->getShapeFunction(i)) ;
					}
				}
				std::vector<Variable> v(2) ;
				v[0] = XI ; v[1] = ETA ;
				Vector imposed(3) ;
				imposed[0] = 0 ;
				imposed[1] = 0 ;
				imposed[2] = data ;
				GaussPointArray gp = e->getGaussPoints() ;
				std::valarray<Matrix> Jinv( Matrix(), gp.gaussPoints.size()) ;
				for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
				{
					e->getInverseJacobianMatrix(gp.gaussPoints[i].first, Jinv[i]) ;
				}
				
				for(size_t i = 0 ; i < shapeFunctions.size() ; ++i)
				{
					Vector forces =  VirtualMachine().ieval(Gradient(shapeFunctions[i]) * (imposed), gp, Jinv, v) ;
					a->addForceOn(XI,forces[0], id[i].id) ;
					a->addForceOn(ETA,forces[1], id[i].id) ;
				}
				return ;
			}
			default:
				break;
		}
	}
}

void apply3DBC(ElementaryVolume *e,  const std::vector<Point> & id, LagrangeMultiplierType condition, double data, Assembly * a)
{
	if(e->getBehaviour()->type == VOID_BEHAVIOUR)
		return ;
	for(size_t i = 0 ; i < id.size() ; i++)
	{
		switch(condition)
		{
			case GENERAL :
				std::cout << "I don't know how to form a General Lagrange Multiplier from the data" << std::endl ;
				break ;
			case FIX_ALONG_XI:
				a->setPointAlong(XI, 0, id[i].id) ;
				break ;
			case SET_ALONG_XI:
				a->setPointAlong(XI, data, id[i].id) ;
				break ;
			case FIX_ALONG_ETA:
				a->setPointAlong(ETA, 0, id[i].id) ;
				break ;
			case SET_ALONG_ETA:
				a->setPointAlong(ETA, data, id[i].id) ;
				break ;
			case FIX_ALONG_ZETA:
				a->setPointAlong(ZETA, 0, id[i].id) ;
				break ;
			case SET_ALONG_ZETA:
				a->setPointAlong(ZETA, data, id[i].id) ;
				break ;
			case SET_FORCE_XI:
				a->setForceOn(XI, data, id[i].id) ;
				break ;
			case SET_FORCE_ETA:
				a->setForceOn(ETA, data, id[i].id) ;
				break ;
			case SET_FORCE_ZETA:
				a->setForceOn(ZETA, data, id[i].id) ;
				break ;
			case SET_STRESS_XI:
			{
				std::vector<Function> shapeFunctions ;
				for(size_t j = 0 ; j < id.size() ; j++)
				{
					for(size_t i = 0 ; i < e->getBoundingPoints().size() ; i++)
					{
						if(id[j].id == e->getBoundingPoint(i).id)
							shapeFunctions.push_back(e->getShapeFunction(i)) ;
					}
				}
				std::vector<Variable> v(3) ;
				v[0] = XI ; v[1] = ETA ; v[1] = ZETA ;
				Vector imposed(6) ;
				imposed[0] = data ;
				imposed[1] = 0 ;
				imposed[2] = 0 ;
				imposed[3] = 0 ;
				imposed[4] = 0 ;
				imposed[5] = 0 ;
				GaussPointArray gp = e->getGaussPoints() ;
				std::valarray<Matrix> Jinv( Matrix(), gp.gaussPoints.size()) ;
				for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
				{
					e->getInverseJacobianMatrix(gp.gaussPoints[i].first, Jinv[i]) ;
				}
				
				for(size_t i = 0 ; i < shapeFunctions.size() ; ++i)
				{
					Vector forces =  VirtualMachine().ieval(Gradient(shapeFunctions[i],true) * (imposed), gp, Jinv, v) ;
					a->addForceOn(XI,forces[0], id[i].id) ;
					a->addForceOn(ETA,forces[1], id[i].id) ;
					a->addForceOn(ZETA,forces[2], id[i].id) ;
				}
				return ;
			}
			case SET_STRESS_ETA:
			{
				std::vector<Function> shapeFunctions ;
				for(size_t j = 0 ; j < id.size() ; j++)
				{
					for(size_t i = 0 ; i < e->getBoundingPoints().size() ; i++)
					{
						if(id[j].id == e->getBoundingPoint(i).id)
							shapeFunctions.push_back(e->getShapeFunction(i)) ;
					}
				}
				std::vector<Variable> v(3) ;
				v[0] = XI ; v[1] = ETA ; v[1] = ZETA ;
				Vector imposed(6) ;
				imposed[0] = 0 ;
				imposed[1] = data ;
				imposed[2] = 0 ;
				imposed[3] = 0 ;
				imposed[4] = 0 ;
				imposed[5] = 0 ;
				GaussPointArray gp = e->getGaussPoints() ;
				std::valarray<Matrix> Jinv( Matrix(), gp.gaussPoints.size()) ;
				for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
				{
					e->getInverseJacobianMatrix(gp.gaussPoints[i].first, Jinv[i]) ;
				}
				
				for(size_t i = 0 ; i < shapeFunctions.size() ; ++i)
				{
					Vector forces =  VirtualMachine().ieval(Gradient(shapeFunctions[i],true) * (imposed), gp, Jinv, v) ;
					a->addForceOn(XI,forces[0], id[i].id) ;
					a->addForceOn(ETA,forces[1], id[i].id) ;
					a->addForceOn(ZETA,forces[2], id[i].id) ;
				}
				return ;
			}
			case SET_STRESS_ZETA:
			{
				std::vector<Function> shapeFunctions ;
				for(size_t j = 0 ; j < id.size() ; j++)
				{
					for(size_t i = 0 ; i < e->getBoundingPoints().size() ; i++)
					{
						if(id[j].id == e->getBoundingPoint(i).id)
							shapeFunctions.push_back(e->getShapeFunction(i)) ;
					}
				}
				std::vector<Variable> v(3) ;
				v[0] = XI ; v[1] = ETA ; v[1] = ZETA ;
				Vector imposed(6) ;
				imposed[0] = 0 ;
				imposed[1] = 0 ;
				imposed[2] = data ;
				imposed[3] = 0 ;
				imposed[4] = 0 ;
				imposed[5] = 0 ;
				GaussPointArray gp = e->getGaussPoints() ;
				std::valarray<Matrix> Jinv( Matrix(), gp.gaussPoints.size()) ;
				for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
				{
					e->getInverseJacobianMatrix(gp.gaussPoints[i].first, Jinv[i]) ;
				}
				
				for(size_t i = 0 ; i < shapeFunctions.size() ; ++i)
				{
					Vector forces =  VirtualMachine().ieval(Gradient(shapeFunctions[i],true) * (imposed), gp, Jinv, v) ;
					a->addForceOn(XI,forces[0], id[i].id) ;
					a->addForceOn(ETA,forces[1], id[i].id) ;
					a->addForceOn(ZETA,forces[2], id[i].id) ;
				}
				return ;
			}
			case SET_STRESS_XI_ETA:
			{
				std::vector<Function> shapeFunctions ;
				for(size_t j = 0 ; j < id.size() ; j++)
				{
					for(size_t i = 0 ; i < e->getBoundingPoints().size() ; i++)
					{
						if(id[j].id == e->getBoundingPoint(i).id)
							shapeFunctions.push_back(e->getShapeFunction(i)) ;
					}
				}
				std::vector<Variable> v(3) ;
				v[0] = XI ; v[1] = ETA ; v[1] = ZETA ;
				Vector imposed(6) ;
				imposed[0] = 0 ;
				imposed[1] = 0 ;
				imposed[2] = 0 ;
				imposed[3] = data ;
				imposed[4] = 0 ;
				imposed[5] = 0 ;
				GaussPointArray gp = e->getGaussPoints() ;
				std::valarray<Matrix> Jinv( Matrix(), gp.gaussPoints.size()) ;
				for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
				{
					e->getInverseJacobianMatrix(gp.gaussPoints[i].first, Jinv[i]) ;
				}
				
				for(size_t i = 0 ; i < shapeFunctions.size() ; ++i)
				{
					Vector forces =  VirtualMachine().ieval(Gradient(shapeFunctions[i],true) * (imposed), gp, Jinv, v) ;
					a->addForceOn(XI,forces[0], id[i].id) ;
					a->addForceOn(ETA,forces[1], id[i].id) ;
					a->addForceOn(ZETA,forces[2], id[i].id) ;
				}
				return ;
			}
			case SET_STRESS_XI_ZETA:
			{
				std::vector<Function> shapeFunctions ;
				for(size_t j = 0 ; j < id.size() ; j++)
				{
					for(size_t i = 0 ; i < e->getBoundingPoints().size() ; i++)
					{
						if(id[j].id == e->getBoundingPoint(i).id)
							shapeFunctions.push_back(e->getShapeFunction(i)) ;
					}
				}
				std::vector<Variable> v(3) ;
				v[0] = XI ; v[1] = ETA ; v[1] = ZETA ;
				Vector imposed(6) ;
				imposed[0] = 0 ;
				imposed[1] = 0 ;
				imposed[2] = 0 ;
				imposed[3] = 0 ;
				imposed[4] = data ;
				imposed[5] = 0 ;
				GaussPointArray gp = e->getGaussPoints() ;
				std::valarray<Matrix> Jinv( Matrix(), gp.gaussPoints.size()) ;
				for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
				{
					e->getInverseJacobianMatrix(gp.gaussPoints[i].first, Jinv[i]) ;
				}
				
				for(size_t i = 0 ; i < shapeFunctions.size() ; ++i)
				{
					Vector forces =  VirtualMachine().ieval(Gradient(shapeFunctions[i],true) * (imposed), gp, Jinv, v) ;
					a->addForceOn(XI,forces[0], id[i].id) ;
					a->addForceOn(ETA,forces[1], id[i].id) ;
					a->addForceOn(ZETA,forces[2], id[i].id) ;
				}
				return ;
			}
			case SET_STRESS_ETA_ZETA:
			{
				std::vector<Function> shapeFunctions ;
				for(size_t j = 0 ; j < id.size() ; j++)
				{
					for(size_t i = 0 ; i < e->getBoundingPoints().size() ; i++)
					{
						if(id[j].id == e->getBoundingPoint(i).id)
							shapeFunctions.push_back(e->getShapeFunction(i)) ;
					}
				}
				std::vector<Variable> v(3) ;
				v[0] = XI ; v[1] = ETA ; v[1] = ZETA ;
				Vector imposed(6) ;
				imposed[0] = 0 ;
				imposed[1] = 0 ;
				imposed[2] = 0 ;
				imposed[3] = 0 ;
				imposed[4] = 0 ;
				imposed[5] = data ;
				GaussPointArray gp = e->getGaussPoints() ;
				std::valarray<Matrix> Jinv( Matrix(), gp.gaussPoints.size()) ;
				for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
				{
					e->getInverseJacobianMatrix(gp.gaussPoints[i].first, Jinv[i]) ;
				}
				
				for(size_t i = 0 ; i < shapeFunctions.size() ; ++i)
				{
					Vector forces =  VirtualMachine().ieval(Gradient(shapeFunctions[i],true) * (imposed), gp, Jinv, v) ;
					a->addForceOn(XI,forces[0], id[i].id) ;
					a->addForceOn(ETA,forces[1], id[i].id) ;
					a->addForceOn(ZETA,forces[2], id[i].id) ;
				}
				return ;
			}
			default:
				break;
		}
	}
}

GeometryDefinedBoundaryCondition::GeometryDefinedBoundaryCondition(LagrangeMultiplierType t, Geometry * source, double d) : BoundaryCondition(t, d), domain(source) { };

void GeometryDefinedBoundaryCondition::apply(Assembly * a, DelaunayTree * t) const
{
	std::vector<ElementarySurface *> & elements = a->getElements2d() ;
	double tol = domain->getRadius()*.0001 ;
	

	for(size_t i = 0 ; i < elements.size() ; ++i)
	{
		std::vector<Point> id  ;
		for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
		{
			Circle c(tol, elements[i]->getBoundingPoint(j)) ; 
			if(domain->intersects(&c))
			{
				id.push_back(elements[i]->getBoundingPoint(j)) ;
			}
		}
		apply2DBC(elements[i], id, condition, data, a) ;
	}
}
void GeometryDefinedBoundaryCondition::apply(Assembly * a, DelaunayTree3D * t)  const
{
	std::vector<ElementaryVolume *> & elements = a->getElements3d() ;
	double tol = domain->getRadius()*.0001 ;
	

	for(size_t i = 0 ; i < elements.size() ; ++i)
	{
		std::vector<Point> id  ;
		for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
		{
			Sphere c(tol, elements[i]->getBoundingPoint(j)) ; 
			if(domain->intersects(&c))
			{
				id.push_back(elements[i]->getBoundingPoint(j)) ;
			}
		}
		apply3DBC(elements[i], id, condition, data, a) ;
	}
}


void BoundingBoxDefinedBoundaryCondition::apply(Assembly * a, DelaunayTree * t) const
{
	std::vector<ElementarySurface *> & elements = a->getElements2d() ;
	double minx = elements.front()->getBoundingPoint(0).x ;
	double miny = elements.front()->getBoundingPoint(0).y ;
	double maxx = elements.front()->getBoundingPoint(0).x ;
	double maxy = elements.front()->getBoundingPoint(0).y ; 
	for(size_t i = 0 ; i < elements.size() ; ++i)
	{
		for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
		{
			if(elements[i]->getBoundingPoint(j).x < minx)
				minx = elements[i]->getBoundingPoint(j).x ;
			if(elements[i]->getBoundingPoint(j).x > maxx)
				maxx = elements[i]->getBoundingPoint(j).x ;
			if(elements[i]->getBoundingPoint(j).y < miny)
				miny = elements[i]->getBoundingPoint(j).y ;
			if(elements[i]->getBoundingPoint(j).y > maxy)
				maxy = elements[i]->getBoundingPoint(j).y ;
		}
	}
	
	double tol = std::min(maxx-minx, maxy-miny)*.0001 ;
	
	switch(pos)
	{
		case TOP:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				std::vector<Point> id  ;
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).y-maxy) < tol)
					{
						id.push_back(elements[i]->getBoundingPoint(j)) ;
					}
				}
				apply2DBC(elements[i], id, condition, data, a) ;
			}
			break ;
		}
		case LEFT:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				std::vector<Point> id  ;
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).x-minx) < tol)
					{
						id.push_back(elements[i]->getBoundingPoint(j)) ;
					}
				}
				apply2DBC(elements[i], id, condition, data, a) ;
			}
			break ;
		}
		case BOTTOM:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				std::vector<Point> id  ;
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).y-miny) < tol)
					{
						id.push_back(elements[i]->getBoundingPoint(j)) ;
					}
				}
				apply2DBC(elements[i], id, condition, data, a) ;
			}
			break ;
		}
		case RIGHT:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				std::vector<Point> id  ;
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).x-maxx) < tol)
					{
						id.push_back(elements[i]->getBoundingPoint(j)) ;
					}
				}
				apply2DBC(elements[i], id, condition, data, a) ;
			}
			break ;
		}
		case TOP_LEFT:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				std::vector<Point> id  ;
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).x-minx) < tol && std::abs(elements[i]->getBoundingPoint(j).y-maxy) < tol)
					{
						id.push_back(elements[i]->getBoundingPoint(j)) ;
					}
				}
				apply2DBC(elements[i], id, condition, data, a) ;
			}
			break ;
		}
		case TOP_RIGHT:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				std::vector<Point> id  ;
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).x-maxx) < tol && std::abs(elements[i]->getBoundingPoint(j).y-maxy) < tol)
					{
						id.push_back(elements[i]->getBoundingPoint(j)) ;
					}
				}
				apply2DBC(elements[i], id, condition, data, a) ;
			}
			break ;
		}
		case BOTTOM_LEFT:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				std::vector<Point> id  ;
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).x-minx) < tol && std::abs(elements[i]->getBoundingPoint(j).y-miny) < tol)
					{
						id.push_back(elements[i]->getBoundingPoint(j)) ;
					}
				}
				apply2DBC(elements[i], id, condition, data, a) ;
			}
			break ;
		}
		case BOTTOM_RIGHT:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				std::vector<Point> id  ;
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).x-maxx) < tol && std::abs(elements[i]->getBoundingPoint(j).y-miny) < tol)
					{
						id.push_back(elements[i]->getBoundingPoint(j)) ;
					}
				}
				apply2DBC(elements[i], id, condition, data, a) ;
			}
			break ;
		}
		default:
		{
			break;
		}
	}
}
void BoundingBoxDefinedBoundaryCondition::apply(Assembly * a, DelaunayTree3D * t)  const 
{
	std::vector<ElementaryVolume *> & elements = a->getElements3d() ;
	double minx = elements.front()->getBoundingPoint(0).x ;
	double miny = elements.front()->getBoundingPoint(0).y ;
	double minz = elements.front()->getBoundingPoint(0).z ;
	double maxx = elements.front()->getBoundingPoint(0).x ;
	double maxy = elements.front()->getBoundingPoint(0).y ; 
	double maxz = elements.front()->getBoundingPoint(0).z ; 
	for(size_t i = 0 ; i < elements.size() ; ++i)
	{
		for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
		{
			if(elements[i]->getBoundingPoint(j).x < minx)
				minx = elements[i]->getBoundingPoint(j).x ;
			if(elements[i]->getBoundingPoint(j).x > maxx)
				maxx = elements[i]->getBoundingPoint(j).x ;
			if(elements[i]->getBoundingPoint(j).y < miny)
				miny = elements[i]->getBoundingPoint(j).y ;
			if(elements[i]->getBoundingPoint(j).y > maxy)
				maxy = elements[i]->getBoundingPoint(j).y ;
			if(elements[i]->getBoundingPoint(j).y < minz)
				minz = elements[i]->getBoundingPoint(j).z ;
			if(elements[i]->getBoundingPoint(j).y > maxz)
				maxz = elements[i]->getBoundingPoint(j).z ;
		}
	}
	
	double tol = std::min(std::min(maxx-minx, maxy-miny), maxz-minz)*.0001 ;
	
	switch(pos)
	{
		case TOP:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				std::vector<Point> id  ;
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).y-maxy) < tol)
					{
						id.push_back(elements[i]->getBoundingPoint(j)) ;
					}
				}
				apply3DBC(elements[i], id, condition, data, a) ;
			}
			break ;
		}
		case LEFT:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				std::vector<Point> id  ;
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).x-minx) < tol)
					{
						id.push_back(elements[i]->getBoundingPoint(j)) ;
					}
				}
				apply3DBC(elements[i], id, condition, data, a) ;
			}
			break ;
		}
		case BOTTOM:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				std::vector<Point> id  ;
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).y-miny) < tol)
					{
						id.push_back(elements[i]->getBoundingPoint(j)) ;
					}
				}
				apply3DBC(elements[i], id, condition, data, a) ;
			}
			break ;
		}
		case RIGHT:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				std::vector<Point> id  ;
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).x-maxx) < tol)
					{
						id.push_back(elements[i]->getBoundingPoint(j)) ;
					}
				}
				apply3DBC(elements[i], id, condition, data, a) ;
			}
			break ;
		}
		case FRONT:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				std::vector<Point> id  ;
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).z-maxz) < tol)
					{
						id.push_back(elements[i]->getBoundingPoint(j)) ;
					}
				}
				apply3DBC(elements[i], id, condition, data, a) ;
			}
			break ;
		}
		case BACK:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				std::vector<Point> id  ;
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).z-minz) < tol)
					{
						id.push_back(elements[i]->getBoundingPoint(j)) ;
					}
				}
				apply3DBC(elements[i], id, condition, data, a) ;
			}
			break ;
		}
		case TOP_LEFT:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				std::vector<Point> id  ;
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).x-minx) < tol && std::abs(elements[i]->getBoundingPoint(j).y-maxy) < tol)
					{
						id.push_back(elements[i]->getBoundingPoint(j)) ;
					}
				}
				apply3DBC(elements[i], id, condition, data, a) ;
			}
			break ;
		}
		case TOP_RIGHT:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				std::vector<Point> id  ;
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).x-maxx) < tol && std::abs(elements[i]->getBoundingPoint(j).y-maxy) < tol)
					{
						id.push_back(elements[i]->getBoundingPoint(j)) ;
					}
				}
				apply3DBC(elements[i], id, condition, data, a) ;
			}
			break ;
		}
		case BOTTOM_LEFT:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				std::vector<Point> id  ;
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).x-minx) < tol && std::abs(elements[i]->getBoundingPoint(j).y-miny) < tol)
					{
						id.push_back(elements[i]->getBoundingPoint(j)) ;
					}
				}
				apply3DBC(elements[i], id, condition, data, a) ;
			}
			break ;
		}
		case BOTTOM_RIGHT:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				std::vector<Point> id  ;
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).x-maxx) < tol && std::abs(elements[i]->getBoundingPoint(j).y-miny) < tol)
					{
						id.push_back(elements[i]->getBoundingPoint(j)) ;
					}
				}
				apply3DBC(elements[i], id, condition, data, a) ;
			}
			break ;
		}
		case FRONT_LEFT:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				std::vector<Point> id  ;
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).x-minx) < tol && std::abs(elements[i]->getBoundingPoint(j).z-maxz) < tol)
					{
						id.push_back(elements[i]->getBoundingPoint(j)) ;
					}
				}
				apply3DBC(elements[i], id, condition, data, a) ;
			}
			break ;
		}
		case FRONT_RIGHT:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				std::vector<Point> id  ;
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).x-maxx) < tol && std::abs(elements[i]->getBoundingPoint(j).z-maxz) < tol)
					{
						id.push_back(elements[i]->getBoundingPoint(j)) ;
					}
				}
				apply3DBC(elements[i], id, condition, data, a) ;
			}
			break ;
		}
		case BACK_LEFT:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				std::vector<Point> id  ;
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).x-minx) < tol && std::abs(elements[i]->getBoundingPoint(j).z-minz) < tol)
					{
						id.push_back(elements[i]->getBoundingPoint(j)) ;
					}
				}
				apply3DBC(elements[i], id, condition, data, a) ;
			}
			break ;
		}
		case BACK_RIGHT:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				std::vector<Point> id  ;
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).x-maxx) < tol && std::abs(elements[i]->getBoundingPoint(j).z-minz) < tol)
					{
						id.push_back(elements[i]->getBoundingPoint(j)) ;
					}
				}
				apply3DBC(elements[i], id, condition, data, a) ;
			}
			break ;
		}
		case FRONT_TOP:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				std::vector<Point> id  ;
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).y-maxy) < tol && std::abs(elements[i]->getBoundingPoint(j).z-maxz) < tol)
					{
						id.push_back(elements[i]->getBoundingPoint(j)) ;
					}
				}
				apply3DBC(elements[i], id, condition, data, a) ;
			}
			break ;
		}
		case FRONT_BOTTOM:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				std::vector<Point> id  ;
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).y-miny) < tol && std::abs(elements[i]->getBoundingPoint(j).z-maxz) < tol)
					{
						id.push_back(elements[i]->getBoundingPoint(j)) ;
					}
				}
				apply3DBC(elements[i], id, condition, data, a) ;
			}
			break ;
		}
		case TOP_LEFT_FRONT:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				std::vector<Point> id  ;
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).x-minx) < tol 
						&& std::abs(elements[i]->getBoundingPoint(j).y-maxy) < tol
						&& std::abs(elements[i]->getBoundingPoint(j).z-maxz) < tol)
					{
						id.push_back(elements[i]->getBoundingPoint(j)) ;
					}
				}
				apply3DBC(elements[i], id, condition, data, a) ;
			}
			break ;
		}
		case TOP_LEFT_BACK:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				std::vector<Point> id  ;
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).x-minx) < tol 
						&& std::abs(elements[i]->getBoundingPoint(j).y-maxy) < tol
						&& std::abs(elements[i]->getBoundingPoint(j).z-minz) < tol)
					{
						id.push_back(elements[i]->getBoundingPoint(j)) ;
					}
				}
				apply3DBC(elements[i], id, condition, data, a) ;
			}
			break ;
		}
		case BOTTOM_LEFT_FRONT:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				std::vector<Point> id  ;
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).x-minx) < tol 
						&& std::abs(elements[i]->getBoundingPoint(j).y-miny) < tol
						&& std::abs(elements[i]->getBoundingPoint(j).z-maxz) < tol)
					{
						id.push_back(elements[i]->getBoundingPoint(j)) ;
					}
				}
				apply3DBC(elements[i], id, condition, data, a) ;
			}
			break ;
		}
		case BOTTOM_LEFT_BACK:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				std::vector<Point> id  ;
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).x-minx) < tol 
						&& std::abs(elements[i]->getBoundingPoint(j).y-miny) < tol
						&& std::abs(elements[i]->getBoundingPoint(j).z-minz) < tol)
					{
						id.push_back(elements[i]->getBoundingPoint(j)) ;
					}
				}
				apply3DBC(elements[i], id, condition, data, a) ;
			}
			break ;
		}
		case TOP_RIGHT_FRONT:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				std::vector<Point> id  ;
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).x-maxx) < tol 
						&& std::abs(elements[i]->getBoundingPoint(j).y-maxy) < tol
						&& std::abs(elements[i]->getBoundingPoint(j).z-minz) < tol)
					{
						id.push_back(elements[i]->getBoundingPoint(j)) ;
					}
				}
				apply3DBC(elements[i], id, condition, data, a) ;
			}
			break ;
		}
		case TOP_RIGHT_BACK:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				std::vector<Point> id  ;
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).x-maxx) < tol 
						&& std::abs(elements[i]->getBoundingPoint(j).y-maxy) < tol
						&& std::abs(elements[i]->getBoundingPoint(j).z-maxz) < tol)
					{
						id.push_back(elements[i]->getBoundingPoint(j)) ;
					}
				}
				apply3DBC(elements[i], id, condition, data, a) ;
			}
			break ;
		}
		case BOTTOM_RIGHT_FRONT:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				std::vector<Point> id  ;
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).x-maxx) < tol 
						&& std::abs(elements[i]->getBoundingPoint(j).y-miny) < tol
						&& std::abs(elements[i]->getBoundingPoint(j).z-minz) < tol)
					{
						id.push_back(elements[i]->getBoundingPoint(j)) ;
					}
				}
				apply3DBC(elements[i], id, condition, data, a) ;
			}
			break ;
		}
		case BOTTOM_RIGHT_BACK:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				std::vector<Point> id  ;
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).x-maxx) < tol 
						&& std::abs(elements[i]->getBoundingPoint(j).y-miny) < tol
						&& std::abs(elements[i]->getBoundingPoint(j).z-maxz) < tol)
					{
						id.push_back(elements[i]->getBoundingPoint(j)) ;
					}
				}
				apply3DBC(elements[i], id, condition, data, a) ;
			}
			break ;
		}
		default:
		{
			break;
		}
	}
}

Form * Feature::getBehaviour( const Point & p)
{
	return this->behaviour ;
}

Feature * Feature::getChild(size_t i) const
{
	return m_c[i] ;
}

const std::vector<Feature *> & Feature::getChildren() const
{
	return m_c ;
}

std::vector<Feature *> & Feature::getChildren()
{
	return m_c ;
}

std::vector<Feature *> Feature::getDescendants() const
{
	std::vector<Feature *> childrenToCheck = m_c ;
	std::vector<Feature *> ret ;
	
	while(!childrenToCheck.empty())
	{
		std::vector<Feature *> newChildren ;
		
		for(size_t i = 0 ; i< childrenToCheck.size() ;i++)
		{
			ret.push_back(childrenToCheck[i]) ;
			newChildren.insert(newChildren.end(), childrenToCheck[i]->getChildren().begin(), childrenToCheck[i]->getChildren().end()) ;
		}
		
		childrenToCheck = newChildren ;
	}
	
	return ret ;
}

void  Feature::setFather(Feature *f)
{
	m_f = f ;
}

 void Feature::setInfluenceRadius(double r)
{
	infRad = r ;
}

CompositeFeature::~CompositeFeature()
{
	for(size_t i = 0 ; i < components.size() ; i++)
		delete components[i] ;
}

std::vector<VirtualFeature *> & CompositeFeature::getComponents()
{
	return components ;
}

const std::vector<VirtualFeature *> & CompositeFeature::getComponents() const
{
	return components ;
}

Feature::~Feature() 
{ 
// 	for(size_t i = 0 ; i < this->boundary->size() ; i++)
// 		delete this->boundary->getPoint(i) ;
	
	delete this->boundary ;
	delete this->boundary2 ;
	
	delete this->behaviour ;
}


DelaunayTree * FeatureTree::getDelaunayTree()
{
	if(this->dtree == NULL && this->dtree3D == NULL)
		this->generateElements() ;
	return this->dtree ;
}

DelaunayTree3D * FeatureTree::getDelaunayTree3D()
{
	if(this->dtree3D == NULL && this->dtree == NULL)
		this->generateElements() ;
	return this->dtree3D ;
}

std::vector<DelaunayTriangle *> FeatureTree::getBoundingTriangles(Feature * f )
{
	if(f == NULL)
		return tree[0]->getBoundingTriangles(dtree) ;
	else
		return f->getBoundingTriangles(dtree) ;
}

FeatureTree::FeatureTree(Feature *first) : grid(NULL), grid3d(NULL)
{



	this->dtree = NULL ;
	this->dtree3D = NULL ;
	if(first)
		this->addFeature(NULL, first) ;

	if(!is3D())
		grid = new Grid(first->getBoundingBox()[1].x-first->getBoundingBox()[0].x,
		                first->getBoundingBox()[1].y-first->getBoundingBox()[2].y, 1,
		                Point((first->getBoundingBox()[1].x+first->getBoundingBox()[0].x)*.5, 
		                      (first->getBoundingBox()[1].y+first->getBoundingBox()[2].y)*.5
		                     )) ;
	if(is3D())
		grid3d =new Grid3D(first->getBoundingBox()[7].x-first->getBoundingBox()[0].x,
		                   first->getBoundingBox()[7].y-first->getBoundingBox()[0].y,
		                   first->getBoundingBox()[7].z-first->getBoundingBox()[0].z, 1, (first->getBoundingBox()[7]+first->getBoundingBox()[0])*.5);
	this->father3D = NULL;
	this->father2D = NULL ;
	this->elemOrder = LINEAR ;
	this->stitched = false ;
	this->renumbered = false ;
	this->needAssembly = true ;
	this->initialized = false ;
	this->setBehaviours = false ;
	meshChange = true ;
	solverConvergence = false ;
	enrichmentChange = true ;

	K = new Assembly() ;
	if(is2D())
		K->set2D() ;
	else
		K->set3D() ;

	crackedVolume = 0 ;
	damagedVolume = 0 ;

}

void FeatureTree::twineFeature(CompositeFeature * father, CompositeFeature * f)
{
	std::vector<Feature *> chain = father->getDescendants() ;
	std::vector<VirtualFeature *> fatherComponents = father->getComponents() ;
	std::vector<VirtualFeature *> childComponents = f->getComponents() ;
	
	//we look for the father of the second descendant of the "father" feature
	Feature * headerEnd =(*std::find(chain.begin(), chain.end(), fatherComponents[0])) ;
	Feature * headerEndChild = (*headerEnd->getChildren().rbegin()) ;
	
	//we attach the first component of the twinee to the twinee
	addFeature(f, childComponents[0]) ;
	
	//the header end loses its child and gets a new one
	headerEnd->removeChild(*headerEnd->getChildren().rbegin());
	addFeature(headerEnd, f) ;
	
	//we re-attach the new header to the rest
	childComponents[0]->addChild(headerEndChild) ;
	headerEndChild->setFather(childComponents[0]) ;
	
	//now that the header is complete, we find each father component in the chain, and insert at that
	//point a new link: the corresponding child component
	
	for(size_t i = 1 ; i < fatherComponents.size() ;i++)
	{
		Feature * attachPoint = (*std::find(chain.begin(), chain.end(), fatherComponents[i])) ;
		Feature * attachPointChild = NULL;
		if(!attachPoint->getChildren().empty())
			attachPointChild = (*attachPoint->getChildren().rbegin()) ;
		
		//we detach the attach point from its father
		attachPoint->removeChild(attachPointChild) ;
		
		//we attach the corresponding child component to the loose end
		addFeature(attachPoint, childComponents[i]) ;

		//we re-attach the attach point to the new feature
		if(attachPointChild)
		{
			childComponents[i]->addChild(attachPointChild) ;
			attachPointChild->setFather(childComponents[i]) ;

		}
		
	}
}


void FeatureTree::addFeature(Feature * father, Feature * f)
{
	if(!tree.empty() && f->spaceDimensions() == SPACE_TWO_DIMENSIONAL && !f->isEnrichmentFeature)
		grid->forceAdd(f) ;
	else if(!tree.empty()&& !f->isEnrichmentFeature)
		grid3d->forceAdd(f) ;
	
	if( f->isCompositeFeature && father && !father->isCompositeFeature)
	{
		std::vector<VirtualFeature *> pile = dynamic_cast<CompositeFeature *>(f)->getComponents();
		f->setFather(father) ;
		if(father != NULL)
			father->addChild(f) ;
		this->tree.push_back(f) ;
		addFeature(f, pile[0]) ;
		for(size_t i = 0 ; i < pile.size()-1 ; i++)
		{
			addFeature(pile[i], pile[i+1]) ;
		}
		return ;
		
	}
	
	f->setFather(father) ;
	if(father != NULL)
		father->addChild(f) ;
	this->tree.push_back(f) ;

}

FeatureTree::~FeatureTree()
{
	delete father3D ;
	delete father2D ;
	delete grid ;
	delete grid3d ;
	delete this->dtree ;
	delete this->dtree3D ;
	delete this->K ;

 	for(size_t i = 0 ; i < additionalPoints.size() ; i++)
		delete additionalPoints[i] ;
	
	for(size_t i = 0 ; i < boundaryCondition.size() ; ++i)
		delete boundaryCondition[i] ;
}

void FeatureTree::addBoundaryCondition(BoundaryCondition * bc)
{
	boundaryCondition.push_back(bc) ;
}

void FeatureTree::removeBoundaryCondition(BoundaryCondition * bc)
{
	std::vector<BoundaryCondition *>::iterator toDelete = std::find(boundaryCondition.begin(), boundaryCondition.end(), bc) ;
	boundaryCondition.erase(toDelete) ;
}


void FeatureTree::setOrder(Order ord)
{
	this->elemOrder = ord ;
	
	if(father3D)
		delete father3D ;
	
	father3D = new TetrahedralElement(this->elemOrder) ;
	father3D->compileAndPrecalculate() ;

	
	if(father2D)
		delete father2D ;
	
	father2D = new TriElement(this->elemOrder) ;
	father2D->compileAndPrecalculate() ;
}

void FeatureTree::renumber()
{

	if(is2D())
	{
		std::vector<DelaunayTriangle *> triangles = this->dtree->getTriangles() ;
		size_t count = 0 ;
		std::cerr << " renumbering... " << std::flush ;

		for(std::vector<DelaunayTriangle *>::iterator i = triangles.begin() ; i != triangles.end() ; ++i)
		{
			for(size_t j = 0 ; j < (*i)->getBoundingPoints().size() ; j++)
			{
				(*i)->getBoundingPoint(j).id = -1 ;
			}
		}
		
		
		for(std::vector<DelaunayTriangle *>::iterator i = triangles.begin() ; i != triangles.end() ; ++i)
		{
			if((*i)->getBehaviour()->type != VOID_BEHAVIOUR)
			{
				for(size_t j = 0 ; j < (*i)->getBoundingPoints().size() ; j++)
				{
					if((*i)->getBoundingPoint(j).id == -1)
						(*i)->getBoundingPoint(j).id = count++ ;
				}
			}
		}

		
		this->dtree->global_counter = count ;
		
		std::cerr << count*2 << " ...done " << std::endl ;

	}
	else if (is3D())
	{
		std::vector<DelaunayTetrahedron *> tets = this->dtree3D->getTetrahedrons() ;
		size_t count = 0 ;
		std::cerr << " renumbering... " << std::flush ;
		
		for(std::vector<DelaunayTetrahedron *>::iterator i = tets.begin() ; i != tets.end() ; ++i)
		{
			for(size_t j = 0 ; j < (*i)->getBoundingPoints().size() ; j++)
			{
				(*i)->getBoundingPoint(j).id = -1 ;
			}
		}
		

		for(std::vector<DelaunayTetrahedron *>::iterator i = tets.begin() ; i != tets.end() ; ++i)
		{
			if((*i)->getBehaviour()->type != VOID_BEHAVIOUR)
			{
				for(size_t j = 0 ; j < (*i)->getBoundingPoints().size() ; j++)
				{
					if((*i)->getBoundingPoint(j).id == -1)
						(*i)->getBoundingPoint(j).id = count++ ;
				}
			}
		}
		this->dtree3D->global_counter = count ;
		
		std::cerr << count*3 << " ...done " << std::endl ;
	}
	
	renumbered = true ;
	
	
}

bool FeatureTree::inRoot(const Point &p) const
{
	return this->tree[0]->in(p) ;
}

void FeatureTree::stitch()
{

	size_t count = 0 ; 
	size_t pd = 0 ;
	if(is2D())
	{
		if((int)elemOrder-1 > 0 )
		{
			switch(elemOrder)
			{
			case CONSTANT:
				{
					break ;
				}
			case LINEAR:
				{
					break ;
				}
			case QUADRATIC:
				{
					this->dtree->addSharedNodes(1,1,0) ;
					break ;
				}
			case CUBIC:
				{
					this->dtree->addSharedNodes(2,1,0) ;
					break ;
				}
			case QUADRIC:
				{
					this->dtree->addSharedNodes(3,1,0) ;
					break ;
				}
			case QUINTIC:
				{
					this->dtree->addSharedNodes(3,1,0) ;
					break ;
				}
			case CONSTANT_TIME_LINEAR:
				{
					this->dtree->addSharedNodes(0,2,2) ;
					break ;
				}
			case CONSTANT_TIME_QUADRATIC:
				{
					this->dtree->addSharedNodes(0,3,2) ;
					break ;
				}
			case LINEAR_TIME_LINEAR:
				{
					this->dtree->addSharedNodes(0,2,2) ;
					break ;
				}
			case LINEAR_TIME_QUADRATIC:
				{
					this->dtree->addSharedNodes(0,3,2) ;
					break ;
				}
			case QUADRATIC_TIME_LINEAR:
				{
					this->dtree->addSharedNodes(1,2,2) ;
					break ;
				}
			case QUADRATIC_TIME_QUADRATIC:
				{
					this->dtree->addSharedNodes(1,3,2) ;
					break ;
				}
			case CUBIC_TIME_LINEAR:
				{
					this->dtree->addSharedNodes(2,2,2) ;
					break ;
				}
			case CUBIC_TIME_QUADRATIC:
				{
					this->dtree->addSharedNodes(2,3,2) ;
					break ;
				}
			case QUADRIC_TIME_LINEAR:
				{
					this->dtree->addSharedNodes(3,2,2) ;
					break ;
				}
			case QUADRIC_TIME_QUADRATIC:
				{
					this->dtree->addSharedNodes(3,3,2) ;
					break ;
				}
			case QUINTIC_TIME_LINEAR:
				{
					this->dtree->addSharedNodes(3,2,2) ;
					break ;
				}
			case QUINTIC_TIME_QUADRATIC:
				{
					this->dtree->addSharedNodes(3,3,2) ;
					break ;
				}
			default:
				break ;
				
			}
			stitched  = true ;	
// 			return ;
			for(size_t j = 1 ; j < this->tree.size() ; j++)
			{
				if(!tree[j]->isEnrichmentFeature)
				{
					
					std::vector<DelaunayTriangle *> triangles = this->tree[j]->getTriangles(dtree) ;
	
					for(size_t i = 0 ; i < triangles.size() ; i++)
					{
						if(true /*|| triangles[i]->Triangle::intersects(tree[j]->getPrimitive())*/)
						{
							Point proj_0(*triangles[i]->first) ;
							tree[j]->project(&proj_0) ;
							Point proj_1(*triangles[i]->second) ;
							tree[j]->project(&proj_1) ;
							Point proj_2(*triangles[i]->third) ;
							tree[j]->project(&proj_2) ;
							bool changed  = true;
							
							if(squareDist2D(&proj_0 , triangles[i]->first ) < POINT_TOLERANCE && 
								squareDist2D(&proj_1 , triangles[i]->second) < POINT_TOLERANCE )
							{
								count+=changed ; 
								changed = false ;
								Point test = triangles[i]->getBoundingPoint(1) ;
								tree[j]->project(&test) ;
								if (inRoot(test)/* && triangles[i]->in(test)*/)
								{
									tree[j]->project(&triangles[i]->getBoundingPoint(1)) ;
									if(elemOrder >= CONSTANT_TIME_LINEAR)
										tree[j]->project(&triangles[i]->getBoundingPoint(7)) ;
									triangles[i]->moved = true ;
									
									for(size_t j = 0 ; j < 3 ; j++)
									{
										if(triangles[i]->getNeighbour(j)->isTriangle)
										{
											dynamic_cast<DelaunayTriangle *>(triangles[i]->getNeighbour(j))->moved = true ;
										}
									}
								}
		// 						std::cerr << "--> " << (*triangles)[i]->getBoundingPoint(1)->x << ", " << (*triangles)[i]->getBoundingPoint(1)->y << std::endl ;
							}
							if(squareDist2D(&proj_1 , triangles[i]->second) < POINT_TOLERANCE && 
								squareDist2D(&proj_2 , triangles[i]->third) < POINT_TOLERANCE )
							{
								count+=changed ; 
								changed = false ;								
								Point test = triangles[i]->getBoundingPoint(3) ;
								tree[j]->project(&test) ;
								if (inRoot(test)/* && triangles[i]->in(test)*/)
								{
									tree[j]->project(&triangles[i]->getBoundingPoint(3)) ;
									if(elemOrder >= CONSTANT_TIME_LINEAR)
										tree[j]->project(&triangles[i]->getBoundingPoint(9)) ;
									triangles[i]->moved = true ;
									for(size_t j = 0 ; j < 3 ; j++)
									{
										if(triangles[i]->getNeighbour(j)->isTriangle)
										{
											dynamic_cast<DelaunayTriangle *>(triangles[i]->getNeighbour(j))->moved = true ;
										}
									}
								}
								
		// 						std::cerr << "--> " << (*triangles)[i]->getBoundingPoint(3)->x << ", " << (*triangles)[i]->getBoundingPoint(3)->y << std::endl ;
							}
							if(squareDist2D(&proj_2 , triangles[i]->third) < POINT_TOLERANCE && 
							   squareDist2D(&proj_0, triangles[i]->first) < POINT_TOLERANCE) 
							{
								count+=changed ; 
								changed = false ;								
								Point test = triangles[i]->getBoundingPoint(5) ;
								tree[j]->project(&test) ;
								if (inRoot(test)/* && triangles[i]->in(test)*/)
								{
									tree[j]->project(&triangles[i]->getBoundingPoint(5)) ;
									if(elemOrder >= CONSTANT_TIME_LINEAR)
										tree[j]->project(&triangles[i]->getBoundingPoint(11)) ;
									triangles[i]->moved = true ;
									for(size_t j = 0 ; j < 3 ; j++)
									{
										if(triangles[i]->getNeighbour(j)->isTriangle)
										{
											dynamic_cast<DelaunayTriangle *>(triangles[i]->getNeighbour(j))->moved = true ;
										}
									}
								
		// 						std::cerr << "--> " << (*triangles)[i]->getBoundingPoint(5)->x << ", " << (*triangles)[i]->getBoundingPoint(5)->y << std::endl ;
									for(size_t j = 0 ; j < 3 ; j++)
									{
										if(triangles[i]->getNeighbour(j)->isTriangle)
										{
											dynamic_cast<DelaunayTriangle *>(triangles[i]->getNeighbour(j))->moved = true ;
										}
									}
								}
							}
							
						}
						if(count % 1000 == 0)
							std::cerr << "\r projecting points on boundaries... triangle " << count << "/" << triangles.size() << " feature " << i << std::flush ; 
						
					}
				}
			}
		}
	}
	else if (is3D())
	{
		if((int)elemOrder-1 > 0 )
		{
			
			switch(elemOrder)
			{
			case CONSTANT:
				{
					break ;
				}
			case LINEAR:
				{
					break ;
				}
			case QUADRATIC:
				{
					this->dtree3D->addSharedNodes(1,1,0) ;
					break ;
				}
			case CUBIC:
				{
					this->dtree3D->addSharedNodes(2,1,0) ;
					break ;
				}
			case QUADRIC:
				{
					this->dtree3D->addSharedNodes(3,1,0) ;
					break ;
				}
			case QUINTIC:
				{
					this->dtree3D->addSharedNodes(3,1,0) ;
					break ;
				}
			case CONSTANT_TIME_LINEAR:
				{
					this->dtree3D->addSharedNodes(0,2,2) ;
					break ;
				}
			case CONSTANT_TIME_QUADRATIC:
				{
					this->dtree3D->addSharedNodes(0,3,2) ;
					break ;
				}
			case LINEAR_TIME_LINEAR:
				{
					this->dtree3D->addSharedNodes(0,2,2) ;
					break ;
				}
			case LINEAR_TIME_QUADRATIC:
				{
					this->dtree3D->addSharedNodes(0,3,2) ;
					break ;
				}
			case QUADRATIC_TIME_LINEAR:
				{
					this->dtree3D->addSharedNodes(1,2,2) ;
					break ;
				}
			case QUADRATIC_TIME_QUADRATIC:
				{
					this->dtree3D->addSharedNodes(1,3,2) ;
					break ;
				}
			case CUBIC_TIME_LINEAR:
				{
					this->dtree3D->addSharedNodes(2,2,2) ;
					break ;
				}
			case CUBIC_TIME_QUADRATIC:
				{
					this->dtree3D->addSharedNodes(2,3,2) ;
					break ;
				}
			case QUADRIC_TIME_LINEAR:
				{
					this->dtree3D->addSharedNodes(3,2,2) ;
					break ;
				}
			case QUADRIC_TIME_QUADRATIC:
				{
					this->dtree3D->addSharedNodes(3,3,2) ;
					break ;
				}
			case QUINTIC_TIME_LINEAR:
				{
					this->dtree3D->addSharedNodes(3,2,2) ;
					break ;
				}
			case QUINTIC_TIME_QUADRATIC:
				{
					this->dtree3D->addSharedNodes(3,3,2) ;
					break ;
				}
			default:
				break ;
				
			}
			stitched = true ;
			return ;
			std::vector<DelaunayTetrahedron *> tets = this->dtree3D->getTetrahedrons() ;
			for(size_t j = 1 ; j < this->tree.size() ; j++)
			{
				if(!tree[j]->isEnrichmentFeature)
				{
					//In two pass
					
					for(size_t i = 0 ; i < tets.size() ; i++)
					{
						Point proj_0(*tets[i]->first) ;
						tree[j]->project(&proj_0) ;
						Point proj_1(*tets[i]->second) ;
						tree[j]->project(&proj_1) ;
						Point proj_2(*tets[i]->third) ;
						tree[j]->project(&proj_2) ;
						Point proj_3(*tets[i]->fourth) ;
						tree[j]->project(&proj_3) ;

						
						
						if(
						    squareDist3D(&proj_0 , tets[i]->first ) < POINT_TOLERANCE*POINT_TOLERANCE && 
						    squareDist3D(&proj_1 , tets[i]->second) < POINT_TOLERANCE*POINT_TOLERANCE 
						  )
						{
							count +=1; 
							Point test = tets[i]->getBoundingPoint(1) ;
							tree[j]->project(&test) ;
							if (inRoot(test) && tets[i]->in(test))
							{
								tree[j]->project(&tets[i]->getBoundingPoint(1)) ;
								if(elemOrder >= CONSTANT_TIME_LINEAR)
									tree[j]->project(&tets[i]->getBoundingPoint(11)) ;
								tets[i]->moved = true ;
							}
						}
						if(
						    squareDist3D(&proj_0 , tets[i]->first ) < POINT_TOLERANCE*POINT_TOLERANCE && 
						    squareDist3D(&proj_2 , tets[i]->third) < POINT_TOLERANCE*POINT_TOLERANCE 
						  )
						{
							count +=1; 
							Point test = tets[i]->getBoundingPoint(9) ;
							tree[j]->project(&test) ;
							if (inRoot(test) && tets[i]->in(test))
							{
								tree[j]->project(&tets[i]->getBoundingPoint(9)) ;
								if(elemOrder >= CONSTANT_TIME_LINEAR)
									tree[j]->project(&tets[i]->getBoundingPoint(19)) ;
								tets[i]->moved = true ;
							}
						}
						if(
						    squareDist3D(&proj_0 , tets[i]->first ) < POINT_TOLERANCE*POINT_TOLERANCE && 
						    squareDist3D(&proj_3 , tets[i]->fourth) < POINT_TOLERANCE*POINT_TOLERANCE 
						  )
						{
							count +=1; 
							Point test = tets[i]->getBoundingPoint(7) ;
							tree[j]->project(&test) ;
							if (inRoot(test) && tets[i]->in(test))
							{
								tree[j]->project(&tets[i]->getBoundingPoint(7)) ;
								if(elemOrder >= CONSTANT_TIME_LINEAR)
									tree[j]->project(&tets[i]->getBoundingPoint(17)) ;
								tets[i]->moved = true ;
							}
						}
						if(
						    squareDist3D(&proj_1 , tets[i]->second ) < POINT_TOLERANCE*POINT_TOLERANCE && 
						    squareDist3D(&proj_3 , tets[i]->fourth) < POINT_TOLERANCE*POINT_TOLERANCE 
						  )
						{
							count +=1; 
							Point test = tets[i]->getBoundingPoint(8) ;
							tree[j]->project(&test) ;
							if (inRoot(test) && tets[i]->in(test))
							{
								tree[j]->project(&tets[i]->getBoundingPoint(8)) ;
								if(elemOrder >= CONSTANT_TIME_LINEAR)
									tree[j]->project(&tets[i]->getBoundingPoint(18)) ;
								tets[i]->moved = true ;
							}
						}
						if(
						    squareDist3D(&proj_1 , tets[i]->second ) < POINT_TOLERANCE*POINT_TOLERANCE && 
						    squareDist3D(&proj_2 , tets[i]->third) < POINT_TOLERANCE*POINT_TOLERANCE
						  )
						{
							count +=1; 
							Point test = tets[i]->getBoundingPoint(3) ;
							tree[j]->project(&test) ;
							if (inRoot(test) && tets[i]->in(test))
							{
								tree[j]->project(&tets[i]->getBoundingPoint(3)) ;
								if(elemOrder >= CONSTANT_TIME_LINEAR)
									tree[j]->project(&tets[i]->getBoundingPoint(13)) ;
								tets[i]->moved = true ;
							}
						}
						if(
						    squareDist3D(&proj_3 , tets[i]->fourth ) < POINT_TOLERANCE*POINT_TOLERANCE && 
						    squareDist3D(&proj_2 , tets[i]->third) < POINT_TOLERANCE*POINT_TOLERANCE 
						  )
						{
							count +=1; 
							Point test = tets[i]->getBoundingPoint(5) ;
							tree[j]->project(&test) ;
							if (inRoot(test) && tets[i]->in(test))
							{
								tree[j]->project(&tets[i]->getBoundingPoint(5)) ;
								if(elemOrder >= CONSTANT_TIME_LINEAR)
									tree[j]->project(&tets[i]->getBoundingPoint(15)) ;
								tets[i]->moved = true ;
							}
						}

						if(count % 1000 == 0)
							std::cerr << "\r projecting points on boundaries... point " << count << "/" << ++pd << " feature " << i << std::flush ; 
					}
				}
			}
		}
	}
	stitched = true ;
	std::cerr << " ...done."<< std::endl ;
}

void FeatureTree::sample(size_t n)
{
	if(is2D())
	{
		std::cerr << "2D features" << std::endl ;
		double total_area = tree[0]->area() ;
		tree[0]->sample(2*n) ;
		for(size_t i  = 1 ; i < this->tree.size() ; i ++)
		{
//			std::cout << i << std::endl ;
			double shape_factor = (sqrt(tree[0]->area())/(2.*M_PI*tree[0]->getRadius()))/(sqrt(tree[i]->area())/(2.*M_PI*tree[i]->getRadius()));
			size_t npoints = std::max((size_t)((double)n*sqrt(tree[i]->area()/(total_area*shape_factor))),(size_t)8) ;
			if(npoints < n)
				tree[i]->sample(npoints) ;
		}
	}
	else if (is3D())
	{
		std::cerr << " 3D features..." << std::flush ;
		double total_area = tree[0]->area() ;
// 		total_area *= tree[0]->area()/(4.*M_PI*tree[0]->getRadius()*tree[0]->getRadius()) ;
		tree[0]->sample(n) ;
		int count = 0 ;
#ifdef HAVE_OPENMP
		omp_set_num_threads(8) ;
#endif
#pragma omp parallel for
		for(int i  = 1 ; i < (int)tree.size() ; i++)
		{
// 			if(i%100 == 0)
			
			std::cerr << "\r 3D features... sampling feature "<< count << "/" << this->tree.size() << "          " << std::flush ;
			if(inRoot(this->tree[i]->getCenter()))
			{
				double a =pow(tree[i]->area()/total_area, .7) ;
				tree[i]->sample(a*n) ;
				count++ ;
				
			}
			else
			{
				double a =tree[i]->area() ;
				tree[i]->sample(a*n) ;
				count++ ;
			}
		}
		std::cerr << "...done" << std::endl ;
	}
}

void FeatureTree::refine(size_t nit, SamplingCriterion *cri)
{
	for(size_t t = 0 ; t < 512 ; t++)
	{
		bool corrected = false ;
		std::vector <DelaunayTriangle *> triangles  =  dtree->getTriangles(true) ;
		
		int count = 0 ;
		for(size_t j = 0;  j < triangles.size() ; j++)
		{
			if(!cri->meetsCriterion(triangles[j]))
			{
				count++ ;
			}
		}
		
		std::cout << count << " non-conformant triangles " << std::endl ;
		for(size_t j = 0;  j < triangles.size() ; j++)
		{
			if(!cri->meetsCriterion(triangles[j]))
			{
				std::vector<Point> temp = cri->suggest(triangles[j]) ;
				if( !temp.empty())
				{
					std::random_shuffle(temp.begin(), temp.end()) ;
					std::cout << "inserting " << temp.size() << " points" << std::endl ;
					for(size_t i = 0 ; i< temp.size() ; i++)
					{
						dtree->insert(new Point(temp[i])) ;
						corrected = true ;
					}
					break ;
				}
			}
		}
		
		if(!corrected)
			break ;
	}
}

void FeatureTree::refine( size_t level )
{
	if(this->dtree == NULL && this->dtree3D == NULL)
		this->generateElements() ;
	
	if(level < 1)
		return ;
	
	if(this->dtree == NULL && this->dtree3D != NULL)
	{
		std::vector<std::pair<std::vector<Geometry *>, Feature *> >zonesVec ;
		
		for(size_t j = 1;  j < this->tree.size() ; j++)
		{
			zonesVec.push_back(std::pair<std::vector<Geometry *>, Feature *>( this->tree[j]->getRefinementZones(level),  this->tree[j])) ;
		}
		
		std::vector<Feature *> enrichmentFeature ;
		for(size_t i  = 0 ; i < this->tree.size() ; i++)
		{
			if(tree[i]->isEnrichmentFeature)
			{
				enrichmentFeature.push_back(tree[i]) ;
			}
		}
		
		size_t points_added = 0 ;
			
		for(size_t i = 0 ; i < zonesVec.size() ; i++)
		{
			for(size_t j = 0 ; j < zonesVec[i].first.size() ; j++)
			{
				if(!zonesVec[i].second->isEnrichmentFeature)
				{

				
					std::vector<Point> toAdd ;
					std::vector<DelaunayTetrahedron *>  tet = this->dtree3D->conflicts(zonesVec[i].first[j]) ;
	// 			std::vector<DelaunayTriangle *>  * tri_in = this->getBoundingTriangles(zonesVec[i].second) ;
				
					for(size_t k = 0 ; k < tet.size() ; k++)
					{
						size_t count_0 = 0 ;
						size_t count_1 = 0 ;
						size_t count_2 = 0 ;
						size_t count_3 = 0 ;
						size_t count_4 = 0 ;
						size_t count_5 = 0 ;
	
						Point p0  = *tet[k]->first*(0.5) + *tet[k]->second*(0.5) ;
						p0.id = -1 ;
						Point p1  = *tet[k]->first*(0.5) + *tet[k]->third*(0.5) ;
						p1.id = -1 ;
						Point p2  = *tet[k]->first*(0.5) + *tet[k]->fourth*(0.5) ;
						p2.id = -1 ;
						Point p3  = *tet[k]->third*(0.5) + *tet[k]->second*(0.5) ;
						p3.id = -1 ;
						Point p4  = *tet[k]->fourth*(0.5) + *tet[k]->third*(0.5) ;
						p4.id = -1 ;
						Point p5  = *tet[k]->second*(0.5) + *tet[k]->fourth*(0.5) ;
						p5.id = -1 ;
						
		
						for(size_t l = 0 ;  l < zonesVec[i].second->getFather()->getChildren().size() ; l++)
						{
							if(!zonesVec[i].second->getFather()->getChild(l)->isEnrichmentFeature )
							{
								if(zonesVec[i].second->getFather()->getChild(l)->inBoundary(&p0))
								{
									count_0++ ;
								}
								if(zonesVec[i].second->getFather()->getChild(l)->inBoundary(&p1))
								{
									count_1++ ;
								}
								if(zonesVec[i].second->getFather()->getChild(l)->inBoundary(&p1))
								{
									count_2++ ;
								}
								if(zonesVec[i].second->getFather()->getChild(l)->inBoundary(&p2))
								{
									count_2++ ;
								}
								if(zonesVec[i].second->getFather()->getChild(l)->inBoundary(&p3))
								{
									count_3++ ;
								}
								if(zonesVec[i].second->getFather()->getChild(l)->inBoundary(&p4))
								{
									count_4++ ;
								}
							}
						}
						
						
	// 					for(size_t m = 0 ; m < enrichmentFeature.size() ; m++)
	// 					{
	// 						if(enrichmentFeature[m]->inBoundary(&p0))
	// 							count_0++ ;
	// 						if(enrichmentFeature[m]->inBoundary(&p1))
	// 							count_1++ ;
	// 						if(enrichmentFeature[m]->inBoundary(&p2))
	// 							count_2++ ;
	// 						if(enrichmentFeature[m]->inBoundary(&p3))
	// 							count_3++ ;
	// 						if(enrichmentFeature[m]->inBoundary(&p4))
	// 							count_4++ ;
	// 						if(enrichmentFeature[m]->inBoundary(&p5))
	// 							count_5++ ;
	// 					}
						
						if(count_0 == 0 )
						{
							toAdd.push_back(p0) ;
						}
						if(count_1 == 0 )
						{
							toAdd.push_back(p1) ;
						}
						if(count_2 == 0 )
						{
							toAdd.push_back(p2) ;
						}
						if(count_3 == 0 )
						{
							toAdd.push_back(p3) ;
						}
						if(count_4 == 0 )
						{
							toAdd.push_back(p4) ;
						}
						if(count_5 == 0 )
						{
							toAdd.push_back(p5) ;
						}
						
					}
		
					
					std::sort(toAdd.begin(), toAdd.end()) ;
					std::vector<Point>::iterator e = std::unique(toAdd.begin(), toAdd.end());
					toAdd.erase(e, toAdd.end()) ;
				
	// 			std::cerr << "we have " << toAdd.size() << " points for refinement" << std::endl ;

					std::random_shuffle(toAdd.begin(), toAdd.end()) ;
					for(size_t k = 0 ; k< toAdd.size() ;k++)
					{
						std::cerr << "\r refining feature " << i+1 << "/" << zonesVec.size() << "...added " << points_added << " points"<< std::flush ;
						Point * p = new Point(toAdd[k]) ;
						points_added++ ;
						this->meshPoints.push_back(std::pair<Point*, Feature *>(p, zonesVec[i].second->getFather())) ;
						this->dtree3D->insert(p) ;
					}
				}
			}
		}
		std::cerr << "...done"<< std::endl ;
	}
	else
	{
			std::vector<std::pair<std::vector<Geometry *>, Feature *> >zonesVec ;
	
	for(size_t j = 1;  j < this->tree.size() ; j++)
	{
		zonesVec.push_back(std::pair<std::vector<Geometry *>, Feature *>( this->tree[j]->getRefinementZones(level),  this->tree[j])) ;
	}
	
	std::vector<Feature *> enrichmentFeature ;
	for(size_t i  = 0 ; i < this->tree.size() ; i++)
	{
		if(tree[i]->isEnrichmentFeature)
		{
			enrichmentFeature.push_back(tree[i]) ;
		}
	}
	
	size_t points_added = 0 ;
		
	for(size_t i = 0 ; i < zonesVec.size() ; i++)
	{

		
		for(size_t j = 0 ; j < zonesVec[i].first.size() ; j++)
		{
			
			if(!zonesVec[i].second->isEnrichmentFeature)
			{
				std::vector<Point *> sample = zonesVec[i].second->doubleSurfaceSampling() ;
				
				for(size_t k = 0 ; k < sample.size() ; k++)
				{
					if(tree[0]->in(*sample[k]))
					{
						bool yes =true ;
						for(size_t l = 0 ; l < enrichmentFeature.size() ; l++)
						{
							if(enrichmentFeature[l]->inBoundary(sample[k]))
							{
								yes = false ;
								break ;
							}
						}
						for(size_t l = 0 ; l < zonesVec[i].second->getChildren().size() ; l++)
						{
							if(zonesVec[i].second->getChild(l)->inBoundary(sample[k]))
							{
								yes = false ;
								break ;
							}
						}
						if(yes)
						{
							dtree->insert(sample[k]) ;
							this->meshPoints.push_back(std::make_pair<Point *, Feature *>(sample[k],zonesVec[i].second )) ;
						}
					}
				}
		}
			
			std::vector<Point> toAdd ;

			std::vector<DelaunayTriangle *>  tri = this->dtree->conflicts(zonesVec[i].first[j]) ;
// 			std::vector<DelaunayTriangle *>  * tri_in = this->getBoundingTriangles(zonesVec[i].second) ;
			
			for(size_t k = 0 ; k < tri.size() ; k++)
			{
// 				if((*tri)[k]->area() > 4e-4)
// 				{
					size_t count_0 = 0 ;
					size_t count_1 = 0 ;
					size_t count_2 = 0 ;
					double rand0 = 0 ;/*((2.*rand()/(RAND_MAX+1.0))-1.)*0.1 ;*/
					double rand1 = 0 ;/*((2.*rand()/(RAND_MAX+1.0))-1.)*0.1 ;*/
					double rand2 = 0 ;/*((2.*rand()/(RAND_MAX+1.0))-1.)*0.1 ;*/
						
					Point p0  = *tri[k]->first*(0.5+rand0) + *tri[k]->second*(0.5-rand0) ;
					p0.id = -1 ;
					Point p1  = *tri[k]->first*(0.5+rand1) + *tri[k]->third*(0.5-rand1) ;
					p1.id = -1 ;
					Point p2  = *tri[k]->second*(0.5+rand2) + *tri[k]->third*(0.5-rand2) ;
					p2.id = -1 ;
// 					
	
					for(size_t l = 0 ;  l < zonesVec[i].second->getFather()->getChildren().size() ; l++)
					{
						if(!zonesVec[i].second->getFather()->getChild(l)->isEnrichmentFeature )
						{
							if(zonesVec[i].second->getFather()->getChild(l)->inBoundary(&p0))
							{
								count_0++ ;
							}
							if(zonesVec[i].second->getFather()->getChild(l)->inBoundary(&p1))
							{
								count_1++ ;
							}
							if(zonesVec[i].second->getFather()->getChild(l)->inBoundary(&p1))
							{
								count_2++ ;
							}
						}
					}
				
			
			
				for(size_t m = 0 ; m < enrichmentFeature.size() ; m++)
				{
					if(enrichmentFeature[m]->inBoundary(&p0))
						count_0++ ;
					if(enrichmentFeature[m]->inBoundary(&p1))
						count_1++ ;
					if(enrichmentFeature[m]->inBoundary(&p2))
						count_2++ ;
				}
				
					if(count_0 == 0 && zonesVec[i].first[j]->in(*tri[k]->first) && zonesVec[i].first[j]->in(*tri[k]->second))
					{
						toAdd.push_back(p0) ;
					}
					if(count_1 == 0 && zonesVec[i].first[j]->in(*tri[k]->first) && zonesVec[i].first[j]->in(*tri[k]->third))
					{
						toAdd.push_back(p1) ;
					}
					if(count_2 == 0 && zonesVec[i].first[j]->in(*tri[k]->second) && zonesVec[i].first[j]->in(*tri[k]->third))
					{
						toAdd.push_back(p2) ;
					}
// 				}
				
			}

			
			std::sort(toAdd.begin(), toAdd.end()) ;
			std::vector<Point>::iterator e = std::unique(toAdd.begin(), toAdd.end());
			toAdd.erase(e, toAdd.end()) ;
			
// 			std::cerr << "we have " << toAdd.size() << " points for refinement" << std::endl ;
			std::random_shuffle(toAdd.begin(), toAdd.end()) ;
			for(size_t k = 0 ; k< toAdd.size() ;k++)
			{
				std::cerr << "\r refining feature " << i+1 << "/" << zonesVec.size() << "...added " << points_added << " points"<< std::flush ;
				Point * p = new Point(toAdd[k]) ;
				points_added++ ;
				this->meshPoints.push_back(std::pair<Point*, Feature *>(p, zonesVec[i].second->getFather())) ;
				this->dtree->insert(p) ;
			}
		}
	}
	std::cerr << "...done"<< std::endl ;
	}
}

Feature * Feature::getFather() const
{
	return this->m_f ;
}

Form * FeatureTree::getElementBehaviour(const DelaunayTriangle * t) const
{

	if(!inRoot(t->getCenter())) 
		return new VoidForm() ;

	if(t->getBoundingPoints().size()%3 != 0)
		return new VoidForm() ;
	
	for(size_t i = 0 ; i < t->getBoundingPoints().size() ; i++)
		if( t->getBoundingPoint(i).id == -1)
		{
			t->getCenter().print() ;
			return new VoidForm() ;
		}
	
	std::vector<Feature *> targets = grid->coOccur(static_cast<const Triangle *>(t)) ;
	if(!targets.empty())
	{
		
		for(int i = targets.size()-1 ; i >=0  ; i--)
		{
			if (!targets[i]->isEnrichmentFeature && targets[i]->in(t->getCenter()))
			{
				
				bool notInChildren  = true ;
				
				std::vector<Feature *> descendants = targets[i]->getDescendants() ;
				
				for(size_t j = 0 ; j < descendants.size() ; j++)
				{
					if(!descendants[j]->isEnrichmentFeature && descendants[j]->in(t->getCenter()))
					{
						notInChildren = false ;
						break ;
					}
				}
				
				if(notInChildren)
				{
					if(targets[i]->getBehaviour(t->getCenter())->timeDependent())
					{
						if( !targets[i]->getBehaviour(t->getCenter())->spaceDependent())
							return targets[i]->getBehaviour(t->getCenter())->getCopy() ;
						else
						{
							Form * b = targets[i]->getBehaviour(t->getCenter())->getCopy() ;
							b->transform(t->getXTransform(), t->getYTransform()) ;
							
							return b ;
						}
					}
					else if(!targets[i]->getBehaviour(t->getCenter())->spaceDependent())
						return targets[i]->getBehaviour(t->getCenter())->getCopy() ;
					else
					{
						Form * b = targets[i]->getBehaviour(t->getCenter())->getCopy() ;
						b->transform(t->getXTransform(), t->getYTransform()) ;
						
						return b ;
					}
					
					return targets[i]->getBehaviour(t->getCenter())->getCopy() ;
				}
			}
		}
	}

	
	if(tree[0]->getBehaviour(t->getCenter())->timeDependent())
	{
		if( !tree[0]->getBehaviour(t->getCenter())->spaceDependent())
			return tree[0]->getBehaviour(t->getCenter())->getCopy() ;
		else
		{
			Form * b = tree[0]->getBehaviour(t->getCenter())->getCopy() ;
			b->transform(t->getXTransform(), t->getYTransform()) ;
			return b ;
		}
	}
	else if(!tree[0]->getBehaviour(t->getCenter())->spaceDependent())
		return tree[0]->getBehaviour(t->getCenter())->getCopy() ;
	else
	{
		Form * b = tree[0]->getBehaviour(t->getCenter())->getCopy() ;
		b->transform(t->getXTransform(), t->getYTransform()) ;
		return b ;
	}
	
	return tree[0]->getBehaviour(t->getCenter())->getCopy() ;

}

Form * FeatureTree::getElementBehaviour(const DelaunayTetrahedron * t) const
{
	
	if(!inRoot(t->getCenter())) 
		return new VoidForm() ;
	
// 	for(size_t i = 0 ; i < t->getBoundingPoints().size() ; i++)
// 		if( t->getBoundingPoint(i).id == -1)
// 		{
// 			return new VoidForm() ;
// 		}
	
	std::vector<Feature *> targets = grid3d->coOccur(t->getCenter()) ;
	if(!targets.empty())
	{
		for(int i = targets.size()-1 ; i >=0  ; i--)
		{
			if (!targets[i]->isEnrichmentFeature && targets[i]->in(t->getCenter()))
			{
				
				
				bool notInChildren  = true ;
				
				std::vector<Feature *> descendants = targets[i]->getDescendants() ;
				
				for(size_t j = 0 ; j < descendants.size() ; j++)
				{
					if(!descendants[j]->isEnrichmentFeature && descendants[j]->in(t->getCenter()))
					{
						notInChildren = false ;
						break ;
					}
				}
				
				if(notInChildren)
				{
					if(targets[i]->getBehaviour(t->getCenter())->timeDependent())
					{
						if( !targets[i]->getBehaviour(t->getCenter())->spaceDependent())
							return targets[i]->getBehaviour(t->getCenter())->getCopy() ;
						else
						{
							Form * b = targets[i]->getBehaviour(t->getCenter())->getCopy() ;
							b->transform(t->getXTransform(), t->getYTransform(), t->getZTransform()) ;
							return b ;
						}
					}
					else if(!targets[i]->getBehaviour(t->getCenter())->spaceDependent())
						return targets[i]->getBehaviour(t->getCenter())->getCopy() ;
					else
					{
						Form * b = targets[i]->getBehaviour(t->getCenter())->getCopy() ;
						b->transform(t->getXTransform(), t->getYTransform(), t->getZTransform()) ;
						return b ;
					}
					
					return targets[i]->getBehaviour(t->getCenter())->getCopy() ;
				}
			}
		}
	}
	else
	{
		for(int i = tree.size()-1 ; i >=0  ; i--)
		{
			if (!tree[i]->isEnrichmentFeature && tree[i]->in(t->getCenter()))
			{
				
				
				bool notInChildren  = true ;
				
				std::vector<Feature *> descendants = tree[i]->getDescendants() ;
				
				for(size_t j = 0 ; j < descendants.size() ; j++)
				{
					if(!descendants[j]->isEnrichmentFeature && descendants[j]->in(t->getCenter()))
					{
						notInChildren = false ;
						break ;
					}
				}
				
				if(notInChildren)
				{
					if(tree[i]->getBehaviour(t->getCenter())->timeDependent())
					{
						if( !tree[i]->getBehaviour(t->getCenter())->spaceDependent())
							return targets[i]->getBehaviour(t->getCenter())->getCopy() ;
						else
						{
							Form * b = tree[i]->getBehaviour(t->getCenter())->getCopy() ;
							b->transform(t->getXTransform(), t->getYTransform(), t->getZTransform()) ;
							return b ;
						}
					}
					else if(!tree[i]->getBehaviour(t->getCenter())->spaceDependent())
						return tree[i]->getBehaviour(t->getCenter())->getCopy() ;
					else
					{
						Form * b = tree[i]->getBehaviour(t->getCenter())->getCopy() ;
						b->transform(t->getXTransform(), t->getYTransform(), t->getZTransform()) ;
						return b ;
					}
					
					return tree[i]->getBehaviour(t->getCenter())->getCopy() ;
				}
			}
		}
	}
	
	if(tree[0]->getBehaviour(t->getCenter())->timeDependent())
	{
		if( !tree[0]->getBehaviour(t->getCenter())->spaceDependent())
			return tree[0]->getBehaviour(t->getCenter())->getCopy() ;
		else
		{
			Form * b = tree[0]->getBehaviour(t->getCenter())->getCopy() ;
			b->transform(t->getXTransform(), t->getYTransform(), t->getZTransform()) ;
			return b ;
		}
	}
	else if(!tree[0]->getBehaviour(t->getCenter())->spaceDependent())
		return tree[0]->getBehaviour(t->getCenter())->getCopy() ;
	else
	{
		Form * b = tree[0]->getBehaviour(t->getCenter())->getCopy() ;
		b->transform(t->getXTransform(), t->getYTransform(), t->getZTransform()) ;
		return b ;
	}
	
	return tree[0]->getBehaviour(t->getCenter())->getCopy() ;

}

// Point * FeatureTree::checkElement( const DelaunayTetrahedron * t ) const
// {
// 		
// 	if(!inRoot(t->getCenter())) 
// 		return NULL;
// 		
// 	for(int i = tree.size()-1 ; i >= 0 ; i--)
// 	{
// 		if (tree[i]->in(t->getCenter()) && (inRoot(t->getCenter())))
// 		{
// 			bool inChild = false ;
// 			
// 			std::vector<Feature *> tocheck = tree[i]->getDescendants();
// 			std::vector<Feature *> tocheckNew =  tocheck;
// 			
// 			for(size_t j = 0 ; j < tocheck.size() ; j++)
// 			{	
// 				if(tocheck[j]->in(t->getCenter()) )
// 				{
// 					inChild = true ;
// 					break ;
// 				}
// 			}
// 				
// 			
// 			if(!inChild)
// 			{
// 				
// 				size_t count_in = 0 ;
// 				
// 				count_in += tree[i]->inBoundary(t->first) ;
// 				count_in += tree[i]->inBoundary(t->second) ;
// 				count_in += tree[i]->inBoundary(t->third) ;
// 				count_in += tree[i]->inBoundary(t->fourth) ;
// 				
// 				if(count_in == 4 && tree[i]->in(t->getCenter()))
// 				{
// 					return NULL;
// 				}
// 				
// 				else 
// 				{
// 					Point *p = new Point(*t->getCircumCenter()) ;
// 					tree[i]->project(p);
// 					return p;
// 				}
// 					
// 			}
// 		}
// 	}
// 	return NULL;
// }


Point * FeatureTree::checkElement( const DelaunayTetrahedron * t ) const
{
	
	if(!inRoot(t->getCenter())) 
		return NULL;
	
	for(int i = tree.size()-1 ; i >= 0 ; i--)
	{
		int inCount = tree[i]->in(*t->first) 
					+ tree[i]->in(*t->second) 
					+ tree[i]->in(*t->third) 
					+ tree[i]->in(*t->fourth) 
					+ tree[i]->in(t->getCenter());
		if (inCount > 2 && inRoot(t->getCenter()))
		{
			bool inChild = false ;
			
			std::vector<Feature *> tocheck = tree[i]->getDescendants();
			
			for(size_t j = 0 ; j < tocheck.size() ; j++)
			{	
				if(tocheck[j]->in(t->getCenter()) )
				{
					inChild = true ;
					break ;
				}
			}
			
			
			if(!inChild)
			{
				
				size_t count_in = 0 ;
				
				count_in += tree[i]->inBoundary(t->first) ;
				count_in += tree[i]->inBoundary(t->second) ;
				count_in += tree[i]->inBoundary(t->third) ;
				count_in += tree[i]->inBoundary(t->fourth) ;
				
				if(count_in == 4 && tree[i]->in(t->getCenter()))
				{
					return NULL;
				}
				
				else 
				{
					Point p1(t->getCenter()) ;
					Point p0(*t->getCircumCenter()) ;
					tree[i]->project(&p1);
					tree[i]->project(&p0);
					double d0 = std::min(dist(&p0, t->first), std::min(dist(&p0, t->second),dist(&p0, t->third))) ;
					double d1 = std::min(dist(&p1, t->first), std::min(dist(&p1, t->second),dist(&p1, t->third))) ;
					if(t->inCircumSphere(p0) 
					   && d0 > 1e-8 
					   && d0>d1
					  ) 
						return new Point(p0);
					else if(t->inCircumSphere(p1) 
					        && d1 > 1e-8 
					        && d1>d0
					       ) 
						return new Point(p1);
					else
						return NULL ;
				}
				
			}
		}
	}
	return NULL;
}



Point * FeatureTree::checkElement( const DelaunayTriangle * t ) const
{
	
	if(!inRoot(t->getCenter())) 
		return NULL;
	
	for(int i = tree.size()-1 ; i >= 0 ; i--)
	{
		int inCount = tree[i]->in(*t->first) 
		         + tree[i]->in(*t->second) 
			+ tree[i]->in(*t->third) +tree[i]->in(t->getCenter());
		if (inCount > 1 && inRoot(t->getCenter()))
		{
			bool inChild = false ;
			
			std::vector<Feature *> tocheck = tree[i]->getDescendants();
			
			for(size_t j = 0 ; j < tocheck.size() ; j++)
			{	
				if(tocheck[j]->in(t->getCenter()) )
				{
					inChild = true ;
					break ;
				}
			}
			
			
			if(!inChild)
			{
				
				size_t count_in = 0 ;
				
				count_in += tree[i]->inBoundary(t->first) ;
				count_in += tree[i]->inBoundary(t->second) ;
				count_in += tree[i]->inBoundary(t->third) ;
				
				if(count_in == 3 && tree[i]->in(t->getCenter()))
				{
					return NULL;
				}
				
				else 
				{
					Point p1(t->getCenter()) ;
					Point p0(t->getCircumCenter()) ;
					tree[i]->project(&p1);
					tree[i]->project(&p0);
					double d0 = std::min(dist(&p0, t->first), std::min(dist(&p0, t->second),dist(&p0, t->third))) ;
					double d1 = std::min(dist(&p1, t->first), std::min(dist(&p1, t->second),dist(&p1, t->third))) ;
					if(t->inCircumCircle(p0) 
					   && d0 > 1e-4 
					   && d0>d1
					  ) 
						return new Point(p0);
					else if(t->inCircumCircle(p1) 
					        && d1 > 1e-4 
					        && d1>d0
					       ) 
						return new Point(p1);
					else
						return NULL ;
				}
				
			}
		}
	}
	return NULL;
}


Feature * FeatureTree::getFeatForTetra( const DelaunayTetrahedron * t ) const
{
			
	if(!inRoot(t->getCenter())) 
		return NULL;
		
	for(int i = tree.size()-1 ; i >= 0 ; i--)
	{
		if (tree[i]->in(t->getCenter()) && (inRoot(t->getCenter())))
		{
			bool inChild = false ;
			
			std::vector<Feature *> tocheck = tree[i]->getChildren();
			std::vector<Feature *> tocheckNew =  tocheck;
			
			while(!tocheckNew.empty())
			{
				std::vector<Feature *> tocheckTemp ;
				for(size_t k = 0 ; k < tocheckNew.size() ; k++)
				{
					tocheckTemp.insert(
									   tocheckTemp.end(), 
							tocheckNew[k]->getChildren().begin(), 
									tocheckNew[k]->getChildren().end()
									  ) ;
				}
				tocheck.insert(tocheck.end(), tocheckTemp.begin(), tocheckTemp.end()) ;
				tocheckNew = tocheckTemp ;
			}
			
			for(size_t j = 0 ; j < tocheck.size() ; j++)
			{	
				if(tocheck[j]->in(t->getCenter()) )
				{
					inChild = true ;
					break ;
				}
			}
			return tree[i];	
		}
	}
	return NULL;
}

void FeatureTree::setElementBehaviours()
{
	if(this->dtree3D == NULL && this->dtree !=NULL)
	{
		std::vector<DelaunayTriangle *> triangles = this->dtree->getTriangles() ;
			
		std::cerr << " refreshing..." << std::flush ;
		this->dtree->refresh(father2D) ;
		std::cerr << " ...done" << std::endl ;
				
		std::cerr << " setting behaviours..." << std::flush ;
		for(size_t i = 0 ; i < triangles.size() ;i++)
		{
			if (i%1000 == 0)
				std::cerr << "\r setting behaviours... triangle " << i << "/" << triangles.size() << std::flush ;
			
			if(!triangles[i]->getBehaviour())
				triangles[i]->setBehaviour(getElementBehaviour(triangles[i])) ;
			
		}
		std::cerr << " ...done" << std::endl ;
		
		setBehaviours = true ;
	}
	else
	{
		std::cerr << " refreshing..." << std::flush ;
		this->dtree3D->refresh(father3D) ;
		std::cerr << " ...done" << std::endl ;
				
		std::vector<DelaunayTetrahedron *> tetrahedrons = this->dtree3D->getTetrahedrons() ;
		
		std::cerr << " setting behaviours..." << std::flush ;
		for(size_t i = 0 ; i < tetrahedrons.size() ;i++)
		{
			if (i%1000 == 0)
				std::cerr << "\r setting behaviours... tet " << i << "/" << tetrahedrons.size() << std::flush ;
			
			if(!tetrahedrons[i]->getBehaviour())
				tetrahedrons[i]->setBehaviour(getElementBehaviour(tetrahedrons[i])) ;
		}
		
		std::cerr << " ...done" << std::endl ;
		
		setBehaviours = true ;
	}
}

void FeatureTree::assemble()
{	
	std::vector<DelaunayTriangle *> triangles ; 
	std::vector<DelaunayTetrahedron *> tetrahedrons ; 
	
	if(!father3D)
		father3D = new TetrahedralElement(this->elemOrder) ;
	father3D->compileAndPrecalculate() ;
	
	if(!father2D)
		father2D = new TriElement(this->elemOrder) ;
	father2D->compileAndPrecalculate() ;
	
	if( this->dtree == NULL && this->dtree3D == NULL)
		this->generateElements() ;
	
	 
	if(!initialized)
	{
		if(this->dtree3D == NULL && this->dtree !=NULL)
		{
			if(!stitched)
				stitch() ;
			
			if(!setBehaviours)
				setElementBehaviours() ;
			
			if(!renumbered)
				renumber() ;
			
			initializeElements() ;
			
		}
		else if(this->dtree3D != NULL && this->dtree ==NULL)
		{
			
			if(!stitched)
				stitch() ;
			
			if(!setBehaviours)
				setElementBehaviours() ;
			
			if(!renumbered)
				renumber() ;

			initializeElements() ;
		}
	}
	
	size_t en_counter =  0 ;
	if(dtree3D)
		en_counter =  this->dtree3D->global_counter ;
	else
		en_counter =  this->dtree->global_counter ;
	
	std::cerr << " enriching..." << std::flush ;		
	
	if(enrichmentChange)
	{
		if(this->dtree3D == NULL && this->dtree !=NULL)
		{
			en_counter = this->dtree->global_counter ;
			initializeElements() ;
		}
		else
		{
			en_counter = this->dtree3D->global_counter ;
			initializeElements() ;

		}
		
		for(size_t i = 1 ; i < this->tree.size() ; i++)
		{
			if(is3D())
			{
				if(this->tree[i]->isEnrichmentFeature && dynamic_cast<EnrichmentFeature *>(this->tree[i])->moved())
				{
					dynamic_cast<EnrichmentFeature *>(this->tree[i])->enrich(en_counter, this->dtree3D) ;
				}
			
				if(i%10 == 0)
					std::cerr << "\r enriching... feature " << i+1 <<"/" << this->tree.size() << std::flush ;
			}
			else
			{
				if(this->tree[i]->isEnrichmentFeature && dynamic_cast<EnrichmentFeature *>(this->tree[i])->moved())
				{
					dynamic_cast<EnrichmentFeature *>(this->tree[i])->enrich(en_counter, this->dtree) ;
				}
			
				if(i%10 == 0)
					std::cerr << "\r enriching... feature " << i+1 <<"/" << this->tree.size() << std::flush ;
			}
		}

		std::cerr << " ...done" << std::endl ;
	}
	
	
	this->numdofs = en_counter ;
	
	if(this->dtree != NULL)
	{
		triangles = this->dtree->getTriangles() ;
		
		for(size_t j = 0 ; j < triangles.size() ; j++)
		{
			if(	triangles[j]->getBehaviour()->type != VOID_BEHAVIOUR)
			{

				std::cerr << "\r assembling stiffness matrix... triangle " << j+1 << "/" << triangles.size() << std::flush ;
				triangles[j]->refresh(father2D) ;
				K->add(triangles[j]) ;
			}
		}
		std::cerr << " ...done." << std::endl ;
	}
	else
	{
		std::vector<DelaunayTetrahedron *> tets = this->dtree3D->getTetrahedrons() ;
		
		for(size_t j = 0 ; j < tets.size() ; j++)
		{
			if(	tets[j]->getBehaviour()->type != VOID_BEHAVIOUR)
			{
				
				if(j%1000 == 0)
					std::cerr << "\r assembling stiffness matrix... tetrahedron " << j+1 << "/" << tets.size() << std::flush ;
				
// 				for(size_t k = 0 ; k < tets[j]->getEnrichmentFunctions().size() ; k++)
// 				{
// 					tets[j]->getEnrichmentFunction(k).compile() ;
// 				}
				tets[j]->refresh(father3D) ;
				K->add(tets[j]) ;
// 				elements3D.push_back(tets[j]) ;
			}
		}
		
		std::cerr << " ...done." << std::endl ;
	}
}


std::vector<DelaunayTriangle> FeatureTree::getSnapshot2D() const
{
	std::vector<DelaunayTriangle> copy ;
	
	std::vector<DelaunayTriangle *> tris = dtree->getTriangles() ;
	
	for(size_t i = 0 ; i < tris.size() ; i++)
	{
		copy.push_back(*tris[i]) ;
		copy.back().setBehaviour(tris[i]->getBehaviour()->getCopy()) ;
		copy.back().getState().initialize() ;
	}
	
	return copy ;
}

Vector FeatureTree::stressFromDisplacements() const
{
	
	if(dtree != NULL)
	{
		std::vector<DelaunayTriangle *> elements = dtree->getTriangles() ;
		Vector stress(0.f, 3*3*elements.size()) ;
	
		for(size_t i  = 0 ; i < elements.size() ; i++)
		{
			if(elements[i]->getBehaviour()->type !=VOID_BEHAVIOUR)
			{
				std::valarray<Point *> pts(3) ;
				pts[0] =  elements[i]->first ;
				pts[1] =  elements[i]->second ;
				pts[2] =  elements[i]->third ;
				
				Vector str = elements[i]->getState().getStress(pts) ;
				for(size_t j = 0 ; j < 9 ; j++)
					stress[i*3*3+j] = str[j] ;
		
				std::cerr << "\r computing stress... element " << i+1 << "/" << elements.size() << std::flush ;
			}
			
		}
		std::cerr << " ...done." << std::endl ;
		return stress ;
	}
	else
	{
		Vector stress(0., 4*6*elements3D.size()) ;
		
		for(size_t i  = 0 ; i < elements3D.size() ; i++)
		{
			std::valarray<Point *> pts(4) ;
			pts[0] =  elements3D[i]->first ;
			pts[1] =  elements3D[i]->second ;
			pts[2] =  elements3D[i]->third ;
			pts[2] =  elements3D[i]->fourth ;
			
			Vector str = elements3D[i]->getState().getStress(pts) ;
			for(size_t j = 0 ; j < 9 ; j++)
				stress[i*4*6+j] = str[j] ;
			
			std::cerr << "\r computing stress... element " << i+1 << "/" << elements3D.size() << std::flush ;
			
		}
		std::cerr << " ...done." << std::endl ;
		return stress ;
	}
}

Vector FeatureTree::getDisplacements() const
{
	return K->getDisplacements() ;
}

std::pair<Vector , Vector > FeatureTree::getStressAndStrain()
{
	if(dtree != NULL)
	{
		std::vector<DelaunayTriangle *> elements = dtree->getTriangles() ;
		std::pair<Vector , Vector > stress_strain(Vector(0., elements[0]->getBoundingPoints().size()*3*elements.size()), Vector(0., elements[0]->getBoundingPoints().size()*3*elements.size())) ;
		for(size_t i  = 0 ; i < elements.size() ; i++)
		{
			if(elements[i]->getBehaviour()->type != VOID_BEHAVIOUR)
			{
// 				std::valarray<Point *> pts(3) ;
// 				pts[0] =  elements[i]->first ;
// 				pts[1] =  elements[i]->second ;
// 				pts[2] =  elements[i]->third ;
				
				std::pair<Vector , Vector > str = elements[i]->getState().getStressAndStrain(elements[i]->getBoundingPoints()) ;
				
				for(size_t j = 0 ; j < elements[0]->getBoundingPoints().size()*3 ; j++)
				{
					stress_strain.first[i*elements[0]->getBoundingPoints().size()*3+j] = str.first[j] ;
					stress_strain.second[i*elements[0]->getBoundingPoints().size()*3+j] = str.second[j] ;
				}
				std::cerr << "\r computing strain+stress... element " << i+1 << "/" << elements.size() << std::flush ;
			}
		}
		std::cerr << " ...done." << std::endl ;
		return stress_strain ;
	}
	else
	{
		std::vector<DelaunayTetrahedron *> tets = dtree3D->getTetrahedrons() ;
		std::pair<Vector , Vector > stress_strain(Vector(0.f, 4*6*tets.size()), Vector(0.f, 4*6*tets.size())) ;
		
		for(size_t i  = 0 ; i < tets.size() ; i++)
		{
			std::valarray<Point *> pts(4) ;
			pts[0] =  tets[i]->first ;
			pts[1] =  tets[i]->second ;
			pts[2] =  tets[i]->third ;
			pts[3] =  tets[i]->fourth ;
			
			std::pair<Vector , Vector > str = tets[i]->getState().getStressAndStrain(pts) ;
			for(size_t j = 0 ; j < 24 ; j++)
			{
				stress_strain.first[i*4*6+j] = str.first[j] ;
				stress_strain.second[i*4*6+j] = str.second[j] ;
			}
			std::cerr << "\r computing strain+stress... element " << i+1 << "/" << tets.size() << std::flush ;
		}
		std::cerr << " ...done." << std::endl ;
		return stress_strain ;
	}
}

std::pair<Vector , Vector > FeatureTree::getStressAndStrain(const std::vector<DelaunayTetrahedron *> & tets)
{
		std::pair<Vector , Vector > stress_strain(Vector(4*6*tets.size()), Vector(4*6*tets.size())) ;
		
		for(size_t i  = 0 ; i < tets.size() ; i++)
		{
			std::valarray<Point *> pts(4) ;
			pts[0] =  tets[i]->first ;
			pts[1] =  tets[i]->second ;
			pts[2] =  tets[i]->third ;
			pts[3] =  tets[i]->fourth ;
			
			std::pair<Vector , Vector > str ;
			str.first.resize(24) ;
			str.second.resize(24) ;
			str = tets[i]->getState().getStressAndStrain(pts) ;
			for(size_t j = 0 ; j < 4 ; j++)
			{
				for(size_t k = 0 ; k < 6 ; k++)
				{
					stress_strain.first[i*4*6+j*6+k] = str.first[j*6+k] ;
					stress_strain.second[i*4*6+j*6+k] = str.second[j*6+k] ;
				}
			}
			std::cerr << "\r computing strain+stress... element " << i+1 << "/" << tets.size() << std::flush ;
		}
		std::cerr << " ...done." << std::endl ;
		return stress_strain ;
}

Vector FeatureTree::strainFromDisplacements() const
{
	if(dtree != NULL)
	{
		std::vector<DelaunayTriangle *> elements = dtree->getTriangles() ;
		Vector strain(0.f, 3*3*elements.size()) ;
		
		for(size_t i  = 0 ; i < elements.size() ; i++)
		{
			if(elements[i]->getBehaviour()->type !=VOID_BEHAVIOUR)
			{
				std::valarray<Point *> pts(3) ;
				pts[0] =  elements[i]->first ;
				pts[1] =  elements[i]->second ;
				pts[2] =  elements[i]->third ;
				
				Vector str = elements[i]->getState().getStrain(pts) ;
				
				for(size_t j = 0 ; j < 9 ; j++)
					strain[i*3*3+j] = str[j] ;
				std::cerr << "\r computing strain... element " << i+1 << "/" << elements.size() << std::flush ;
			}
		}
		std::cerr << " ...done." << std::endl ;
		return strain ;
	}
	else
	{
		Vector strain(0., 4*6*elements3D.size()) ;
		
		for(size_t i  = 0 ; i < elements3D.size() ; i++)
		{
			std::valarray<Point *>  pts(4) ;
			pts[0] =  elements3D[i]->first ;
			pts[1] =  elements3D[i]->second ;
			pts[2] =  elements3D[i]->third ;
			pts[3] =  elements3D[i]->fourth ;
			
			
			Vector str = elements3D[i]->getState().getStrain(pts) ;
			
			for(size_t j = 0 ; j < 24 ; j++)
				strain[i*4*6+j] = str[j] ;
			std::cerr << "\r computing strain... element " << i+1 << "/" << elements3D.size() << std::flush ;
		}
		std::cerr << " ...done." << std::endl ;
		return strain ;
	}
	
}

Assembly * FeatureTree::getAssembly()
{
	if(!stitched)
		stitch() ;
	
	if(!setBehaviours)
		setElementBehaviours() ;
	
	if(!renumbered)
		renumber() ;
	
	return K ;
}

void FeatureTree::insert(Point * p )
{
	Feature * mother = NULL;
	
	for(size_t i  = 0 ; i < this->tree.size() ; i ++)
	{
		if(this->tree[i]->in((*p)))
		{
			bool yes = true ;
			
			for(size_t k  =  0 ; k < this->tree[i]->getChildren().size() && yes; k++)
			{
				if(this->tree[i]->getChild(k)->in((*p)))
					yes = false ;
			}
			
			if(yes)
			{
				mother = this->tree[i] ;
				break ;
			}
		}
	}
	
	bool yes = true ;
	
	for(size_t k  =  0 ; k <  mother->getChildren().size() && yes; k++)
	{
		if( mother->getChild(k)->inBoundary(p))
			yes = false ;
	}
	
	if(yes)
	{
		this->meshPoints.push_back(std::pair<Point*, Feature *>(p, mother)) ;
		if(dtree != NULL)
			this->dtree->insert(p) ;
		else
			this->dtree3D->insert(p) ;
	}
}


void FeatureTree::stepBack()
{
	if(is2D())
	{
		std::vector<DelaunayTriangle *> elements = dtree->getTriangles() ;
		
		for(size_t i = 0 ; i < elements.size() ;i++)
		{	
			if(i%1000 == 0)
				std::cerr << "\r stepping through elements... " << i << "/" << elements.size() << std::flush ;
			elements[i]->stepBack() ;
		}
		std::cerr << " ...done" << std::endl ;
	}
	else
	{
		std::vector<DelaunayTetrahedron *> elements = dtree3D->getTetrahedrons() ;
		
		for(size_t i = 0 ; i < elements.size() ;i++)
		{	
			if(i%1000 == 0)
				std::cerr << "\r stepping through elements... " << i << "/" << elements.size() << std::flush ;
			elements[i]->stepBack() ;
		}
		std::cerr << " ...done" << std::endl ;
	}
}

bool FeatureTree::solverConverged() const
{
	return solverConvergence ;
}

bool FeatureTree::meshChanged() const
{
	return meshChange ;
}

bool FeatureTree::enrichmentChanged() const
{
	return enrichmentChange ;
}

void FeatureTree::elasticStep()
{
	Vector lastx(K->getDisplacements()) ;
	this->K->clear() ;
	assemble() ;

	this->K->cgsolve(lastx) ;
	if(is2D())
	{
		std::vector<DelaunayTriangle *> elements = dtree->getTriangles() ;
		
		//this will update the state of all elements. This is necessary as 
		//the behaviour updates might depend on the global state of the 
		//simulation.
		std::cerr << " stepping through elements... " << std::flush ;
		for(size_t i = 0 ; i < elements.size() ;i++)
		{	
			if(i%1000 == 0)
				std::cerr << "\r stepping through elements... " << i << "/" << elements.size() << std::flush ;
			elements[i]->step(0., &K->getDisplacements()) ;
		}
		std::cerr << " ...done" << std::endl ;

		
	}
	else if(is3D())
	{
		
		std::vector<DelaunayTetrahedron *> elements = dtree3D->getTetrahedrons() ;
		std::cerr << " stepping through elements... " << std::flush ;
		
		for(size_t i = 0 ; i < elements.size() ;i++)
		{	
			if(i%1000 == 0)
				std::cerr << "\r stepping through elements... " << i << "/" << elements.size() << std::flush ;
			elements[i]->step(0., &K->getDisplacements()) ;
			std::cerr << " ...done" << std::endl ;
		}
	}
	
}

bool FeatureTree::step(double dt)
{
	Vector lastx(K->getDisplacements()) ;
	bool ret = true ;
	if(enrichmentChange)
		this->K->clear() ;

	assemble() ;
	for(size_t i = 0 ; i < boundaryCondition.size() ; ++i)
	{
		if(dtree)
			boundaryCondition[i]->apply(K, dtree) ;
		if(dtree3D)
			boundaryCondition[i]->apply(K, dtree3D) ;
	}
	
	needAssembly = true ;
	
	meshChange = false ;
	if(solverConvergence)
		solverConvergence = this->K->cgsolve(lastx) ;
	else
	{
		lastx = 0 ;
		solverConvergence = this->K->cgsolve(lastx) ;
	}
	enrichmentChange = false ;

	if(is2D())
	{
		std::vector<DelaunayTriangle *> elements = dtree->getTriangles() ;
		
		double volume = 0;	
		crackedVolume = 0 ;	
		damagedVolume = 0 ;	
		//this will update the state of all elements. This is necessary as 
		//the behaviour updates might depend on the global state of the 
		//simulation.
		std::cerr << " stepping through elements... " << std::flush ;
		for(size_t i = 0 ; i < elements.size() ;i++)
		{	
			if(i%1000 == 0)
				std::cerr << "\r stepping through elements... " << i << "/" << elements.size() << std::flush ;
			elements[i]->step(dt, &K->getDisplacements()) ;
		}
		std::cerr << " ...done" << std::endl ;
		int fracturedCount = 0 ;
		
		for(size_t i = 0 ; i < elements.size() ;i++)
		{	
			if(i%1000 == 0)
				std::cerr << "\r checking for fractures... " << i << "/" << elements.size() << std::flush ;
			
			if(elements[i]->getBehaviour()->type !=VOID_BEHAVIOUR )
			{
				volume += elements[i]->area() ;
				
				elements[i]->getBehaviour()->step(dt, elements[i]->getState()) ;
				
				if(elements[i]->getBehaviour()->changed())
				{
					needAssembly = true ;
					meshChange = true ;
					ret = false ;
				}
				
				if(elements[i]->getBehaviour()->fractured())
				{
					fracturedCount++ ;
					crackedVolume +=  elements[i]->area() ;
				}
				else if(std::abs(elements[i]->getBehaviour()->getTensor(Point(.3,.3))[0][0] - elements[i]->getBehaviour()->param[0][0]) > 1e-12*elements[i]->getBehaviour()->getTensor(Point(.3,.3))[0][0] )
				{
					damagedVolume +=  elements[i]->area() ;
				}
			}
			else if (elements[i]->getBehaviour()->fractured())
				crackedVolume +=  elements[i]->area() ;
		}
		std::cerr << " ...done" << std::endl ;
		for(size_t i = 0 ; i < elements.size() ;i++)
			elements[i]->clearVisited() ;
// 		std::cout << " Fractured " << fracturedCount << " Elements" << std::endl ;
// 		std::cout << " Fractured Fraction " <<  crackedVolume / volume << std::endl ;
		
		for(size_t i = 0 ; i< tree.size() ; i++)
		{
			if(tree[i]->isEnrichmentFeature)
			{
				dynamic_cast<EnrichmentFeature *>(tree[i])->step(dt, &K->getForces(), dtree) ;
				enrichmentChange = enrichmentChange || dynamic_cast<EnrichmentFeature *>(tree[i])->moved() ;
				needAssembly = true ;
			}
		}

// 		CrackInitiation().step(.01, dtree) ;
		
	}
	else if(is3D())
	{
		
		std::vector<DelaunayTetrahedron *> elements = dtree3D->getTetrahedrons() ;
		
		//this will update the state of all elements. This is necessary as 
		//the behaviour updates might depend on the global state of the 
		//simulation.
		double volume = 0;	
		crackedVolume = 0 ;	
		damagedVolume = 0 ;	
		//this will update the state of all elements. This is necessary as 
		//the behaviour updates might depend on the global state of the 
		//simulation.
		std::cerr << " stepping through elements... " << std::flush ;
		for(size_t i = 0 ; i < elements.size() ;i++)
		{	
			if(i%1000 == 0)
				std::cerr << "\r stepping through elements... " << i << "/" << elements.size() << std::flush ;
			elements[i]->step(dt, &K->getDisplacements()) ;
		}
		std::cerr << " ...done" << std::endl ;
		int fracturedCount = 0 ;
		
		for(size_t i = 0 ; i < elements.size() ;i++)
		{	
			if(i%1000 == 0)
				std::cerr << "\r checking for fractures... " << i << "/" << elements.size() << std::flush ;
			
			if(elements[i]->getBehaviour()->type !=VOID_BEHAVIOUR )
			{
				volume += elements[i]->area() ;
				
				elements[i]->getBehaviour()->step(dt, elements[i]->getState()) ;
				
				if(elements[i]->getBehaviour()->changed())
				{
					needAssembly = true ;
					meshChange = true ;
					ret = false ;
				}
				
				if(elements[i]->getBehaviour()->fractured())
				{
					fracturedCount++ ;
					crackedVolume +=  elements[i]->area() ;
				}
				else if(std::abs(elements[i]->getBehaviour()->getTensor(Point(.25,.25,.25 ))[0][0] - elements[i]->getBehaviour()->param[0][0]) > 1e-12*elements[i]->getBehaviour()->getTensor(Point(.25,.25,.25))[0][0] )
				{
					damagedVolume +=  elements[i]->area() ;
				}
			}
			else if (elements[i]->getBehaviour()->fractured())
				crackedVolume +=  elements[i]->area() ;
		}
		std::cerr << " ...done" << std::endl ;
		for(size_t i = 0 ; i < elements.size() ;i++)
			elements[i]->clearVisited() ;
		// 		std::cout << " Fractured " << fracturedCount << " Elements" << std::endl ;
		// 		std::cout << " Fractured Fraction " <<  crackedVolume / volume << std::endl ;
		
		for(size_t i = 0 ; i< tree.size() ; i++)
		{
			if(tree[i]->isEnrichmentFeature)
			{
				dynamic_cast<EnrichmentFeature *>(tree[i])->step(dt, &K->getForces(), dtree) ;
				enrichmentChange = enrichmentChange || dynamic_cast<EnrichmentFeature *>(tree[i])->moved() ;
				needAssembly = true ;
			}
		}
	}
	return ret ;
	
}

double FeatureTree::getMaximumDisplacement() const
{
	if(is2D())
	{
		std::vector<DelaunayTriangle *> tri = dtree->getTriangles() ;
		
		double max = 0 ;
		
		for(std::vector<DelaunayTriangle *>::const_iterator i = tri.begin() ; i != tri.end() ; ++i)
		{
			if((*i)->getBehaviour()->type != VOID_BEHAVIOUR)
				max = std::max(max, (*i)->getState().getDisplacements().max()) ;
		}
		
		return max ;
	}
	else if(is3D())
	{
		std::vector<DelaunayTetrahedron *> tets = dtree3D->getTetrahedrons() ;
		
		double max = 0 ;
		
		for(std::vector<DelaunayTetrahedron *>::const_iterator i = tets.begin() ; i != tets.end() ; ++i)
		{
			if((*i)->getBehaviour()->type != VOID_BEHAVIOUR)
				max = std::max(max, (*i)->getState().getDisplacements().max()) ;
		}
		
		return max ;
	}
	
	return 0 ;
}

double FeatureTree::getMinimumDisplacement() const
{
	if(is2D())
	{
		std::vector<DelaunayTriangle *> tri = dtree->getTriangles() ;
		
		double max = 0 ;
		
		for(std::vector<DelaunayTriangle *>::const_iterator i = tri.begin() ; i != tri.end() ; ++i)
		{
			if((*i)->getBehaviour()->type != VOID_BEHAVIOUR)
				max = std::min(max, (*i)->getState().getDisplacements().min()) ;
		}
		
		return max ;
	}
	else if(is3D())
	{
		std::vector<DelaunayTetrahedron *> tets = dtree3D->getTetrahedrons() ;
		
		double max = 0 ;
		
		for(std::vector<DelaunayTetrahedron *>::const_iterator i = tets.begin() ; i != tets.end() ; ++i)
		{
			if((*i)->getBehaviour()->type != VOID_BEHAVIOUR)
				max = std::min(max, (*i)->getState().getDisplacements().min()) ;
		}
		
		return max ;
	}
	
	return 0 ;
}

size_t FeatureTree::numPoints() const
{
	return this->numdofs ;
}

std::deque<std::pair<Point *, Feature *> >::iterator FeatureTree::begin()
{
	return this->meshPoints.begin() ;
}

std::deque<std::pair<Point *, Feature *> >::iterator FeatureTree::end()
{
	return this->meshPoints.end() ;
}


void FeatureTree::print() const
{
	printForFeature(tree[0]);
	
}

void FeatureTree::printForFeature(const Feature *f) const
{
	f->print();
	std::vector<Feature *> children = f->getChildren();
	for (size_t i = 0; i != children.size(); ++i)
	{
// 		if ( !(*children)[i]->getChildren().empty()) 
			printForFeature(children[i]);
	}

}

bool FeatureTree::is3D() const
{
	return tree[0]->spaceDimensions() == SPACE_THREE_DIMENSIONAL ;
}

bool FeatureTree::is2D() const
{
	return tree[0]->spaceDimensions() == SPACE_TWO_DIMENSIONAL ;
}

void FeatureTree::initializeElements() 
{
	if(this->dtree3D == NULL && this->dtree !=NULL)
	{
		std::cerr << " refreshing..." << std::flush ;
		this->dtree->refresh(father2D) ;
		std::cerr << " ...done" << std::endl ;
	}
	else if(this->dtree3D != NULL && this->dtree ==NULL)
	{
		std::cerr << " refreshing..." << std::flush ;
		this->dtree3D->refresh(father3D) ;
		std::cerr << " ...done" << std::endl ;
	}
	
	if(dtree)
	{
		std::vector<DelaunayTriangle *> triangles = this->dtree->getTriangles() ;
		std::cout << " initialising..." ;
		for(size_t i = 0 ; i < triangles.size() ;i++)
		{
			if(i%1000 == 0)
				std::cout << "\r initialising... element " << i << "/" << triangles.size() << std::flush ;
			triangles[i]->getState().initialize() ;
		}
		std::cout << "\r initialising... element " << triangles.size() << "/" << triangles.size() << std::endl ;
	}
	
	if(dtree3D)
	{
		std::vector<DelaunayTetrahedron *> tets = this->dtree3D->getTetrahedrons() ;
		std::cout << " initialising..." ;
		for(size_t i = 0 ; i < tets.size() ;i++)
		{
			if(i%1000 == 0)
				std::cout << "\r initialising... element " << i << "/" << tets.size() << std::flush ;
			tets[i]->getState().initialize() ;
		}
		std::cout << "\r initialising... element " << tets.size() << "/" << tets.size() << std::endl ;
	}
	
	initialized = true ;
}

void FeatureTree::generateElements( size_t correctionSteps, bool computeIntersections) 
{

	std::valarray<Point> bbox(8) ;
	double min_x = 0, min_y = 0, max_x = 0, max_y = 0, max_z = 0, min_z = 0;

	for(size_t j  =  0 ; j <  this->tree[0]->getBoundingPoints().size() ; j++)
	{
		if(this->tree[0]->getBoundingPoint(j).y < min_y)
			min_y = this->tree[0]->getBoundingPoint(j).y ;
		if(this->tree[0]->getBoundingPoint(j).y > max_y)
			max_y = this->tree[0]->getBoundingPoint(j).y ;
		
		if(this->tree[0]->getBoundingPoint(j).x < min_x)
			min_x = this->tree[0]->getBoundingPoint(j).x ;
		if(this->tree[0]->getBoundingPoint(j).x > max_x)
			max_x = this->tree[0]->getBoundingPoint(j).x ;
		
		if(this->tree[0]->getBoundingPoint(j).z < min_z)
			min_z = this->tree[0]->getBoundingPoint(j).z ;
		if(this->tree[0]->getBoundingPoint(j).z > max_z)
			max_z = this->tree[0]->getBoundingPoint(j).z ;
	}

	bbox[0] = Point(min_x, min_y, min_z) ;
	bbox[1] = Point(min_x, min_y, max_z) ;
	bbox[2] = Point(min_x, max_y, min_z) ;
	bbox[3] = Point(min_x, max_y, max_z) ;
	bbox[4] = Point(max_x, min_y, min_z) ;
	bbox[5] = Point(max_x, min_y, max_z) ;
	bbox[6] = Point(max_x, max_y, min_z) ;
	bbox[7] = Point(max_x, max_y, max_z) ;

	std::vector<Feature *> enrichmentFeature ;
	for(size_t i  = 0 ; i < this->tree.size() ; i++)
	{
		if(tree[i]->isEnrichmentFeature)
		{
			enrichmentFeature.push_back(tree[i]) ;
		}
	}
	int bpcount = 0 ;
	size_t basepoints = 0 ;
	std::cerr << " getting mesh points..." << std::flush ;
	for(size_t i  = 0 ; i < tree.size() ; i++)
	{
		std::cerr << "\r getting mesh points... feature " << i << "/"<< tree.size() << std::flush ;
		if(!tree[i]->isEnrichmentFeature && !tree[i]->isVirtualFeature)
		{
			std::vector<Feature *> descendants = tree[i]->getDescendants() ;
			std::stable_sort(descendants.begin(), descendants.end()) ;
			for(size_t j  =  0 ; j <  tree[i]->getBoundingPoints().size() ; j++)
			{
				
				bool isIn = false ;
				
				std::vector<Feature *> potentialFeatures ;
				
				if(tree[0]->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
					potentialFeatures = grid->coOccur(tree[i]->getBoundingPoint(j)) ;
				else
				{
					potentialFeatures = grid3d->coOccur(tree[i]->getBoundingPoint(j)) ;
				}
				
				std::vector<Feature *> potentialChildren ;
				for(size_t l = 0 ; l < potentialFeatures.size() ; l++)
				{
					if(!potentialFeatures[l]->isEnrichmentFeature 
						&& std::binary_search(descendants.begin(), descendants.end(), potentialFeatures[l]))
						potentialChildren.push_back(potentialFeatures[l]) ;
				}
// 				for(size_t k = 0 ; k < potentialFeatures.size() ; k++)
// 				{
// 					for(size_t l = 0 ; l < descendants.size() ; l++)
// 					{
// 						if(potentialFeatures[k] == descendants[l])
// 						{
// 							
// 							break ;
// 						}
// 					}
// 				}
				
				for(size_t k  =  0 ; k <  potentialChildren.size() ; k++)
				{
					if((!potentialChildren[k]->isVirtualFeature 
					    && potentialChildren[k]->getBoundary()->in(tree[i]->getBoundingPoint(j))) 
					   || (potentialChildren[k]->isVirtualFeature 
					       && tree[i]->isVirtualFeature 
					       && (dynamic_cast<VirtualFeature *>(potentialChildren[k])->getSource() 
					           != dynamic_cast<VirtualFeature *>(tree[i])->getSource())
					       && potentialChildren[k]->getBoundary()->in(tree[i]->getBoundingPoint(j))
					      )
					   || (potentialChildren[k]->isVirtualFeature 
					       && potentialChildren[k]->in(tree[i]->getBoundingPoint(j)))
					  )
					{
						isIn = true ;
						break ;
					}
// 					else if(potentialChildren[k]->isVirtualFeature 
// 					        && potentialChildren[k]->in(tree[i]->getBoundingPoint(j)))
// 					{
// 						isIn = true ;
// 						break ;
// 					}
				}
				

				if(i != 0 && !inRoot(tree[i]->getBoundingPoint(j)))
					isIn = true ;
				if(tree[i]->isVirtualFeature && !tree[i]->in(tree[i]->getBoundingPoint(j)))
					isIn = true ;
				
				if(!isIn)
				{
					bpcount++ ;
					meshPoints.push_back(std::pair<Point *, Feature *>(&tree[i]->getBoundingPoint(j), this->tree[i])) ;
					if(i == 0)
						basepoints++ ;
				}
			}

			for(size_t j  =  0 ; j <  tree[i]->getInPoints().size() ; j++)
			{
				bool isIn = false ;
				std::vector<Feature *> potentialFeatures  ;
				if(tree[0]->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
					potentialFeatures = grid->coOccur(tree[i]->getInPoint(j)) ;
				else
					potentialFeatures = grid3d->coOccur(tree[i]->getInPoint(j)) ;
				std::vector<Feature *> potentialChildren ;
				
				for(size_t l = 0 ; l < potentialFeatures.size() ; l++)
				{
					if(!potentialFeatures[l]->isVirtualFeature
					   &&!potentialFeatures[l]->isEnrichmentFeature 
					   && std::binary_search(descendants.begin(), descendants.end(), potentialFeatures[l] ))
						potentialChildren.push_back(potentialFeatures[l]) ;
				}
// 				for(size_t k = 0 ; k < potentialFeatures.size() ; k++)
// 				{
// 					for(size_t l = 0 ; l < descendants.size() ; l++)
// 					{
// 						if(potentialFeatures[k] == descendants[l])
// 						{
// 							potentialChildren.push_back(potentialFeatures[k]) ;
// 							break ;
// 						}
// 					}
// 				}
				for(size_t k  =  0 ; k <  potentialChildren.size() ; k++)
				{
					if(
					    (
					      !potentialChildren[k]->isVirtualFeature 
					      && potentialChildren[k]->getBoundary()->in(tree[i]->getInPoint(j))
					    ) 
					   || 
					    (
					      (
					        potentialChildren[k]->isVirtualFeature 
					        && tree[i]->isVirtualFeature 
					        && 
					        (
					          dynamic_cast<VirtualFeature *>(potentialChildren[k])->getSource() 
					          != dynamic_cast<VirtualFeature *>(tree[i])->getSource()
					        )
					        && potentialChildren[k]->getBoundary()->in(tree[i]->getInPoint(j))
					      )
					      || 
					      (
					        potentialChildren[k]->isVirtualFeature 
					        && potentialChildren[k]->in(tree[i]->getInPoint(j))
					      )
					    )
					  )
					{
						isIn = true ;
						break ;
					}
// 					else if(potentialChildren[k]->isVirtualFeature 
// 					        && potentialChildren[k]->in(tree[i]->getInPoint(j)))
// 					{
// 						isIn = true ;
// 						break ;
// 					}

				}
				
				
				if(i != 0 && !inRoot(tree[i]->getInPoint(j)))
					isIn = true ;
				if(tree[i]->isVirtualFeature && !tree[i]->in(tree[i]->getInPoint(j)))
					isIn = true ;
					
				if(!isIn)
				{
					meshPoints.push_back(std::pair<Point *, Feature *>(&tree[i]->getInPoint(j), tree[i])) ;	
					if(i == 0)
						basepoints++ ;
				}
			}
			
		}
	}
	

	std::cerr << "...done" << std::endl ;
	
	size_t count  = 0 ;

	if(computeIntersections)
	{
		for(size_t i = 1 ;  i < tree.size() ; i++)
		{
			if(!tree[i]->isEnrichmentFeature && !tree[i]->isVirtualFeature)
			{
				std::vector<Feature *> coOccuringFeatures ;
				std::vector<Feature *> descendants = tree[i]->getDescendants() ;

				if(is3D())
					coOccuringFeatures = grid3d->coOccur(tree[i]) ;
				else
					coOccuringFeatures = grid->coOccur(tree[i]) ;
				for(size_t j  = 0 ; j < coOccuringFeatures.size() ; j++)
				{
					if(!coOccuringFeatures[j]->isEnrichmentFeature 
					   && !coOccuringFeatures[j]->isVirtualFeature 
					   && tree[i] != coOccuringFeatures[j] 
					   && tree[i]->intersects(coOccuringFeatures[j]))
					{
						std::vector<Point> inter = tree[i]->intersection(coOccuringFeatures[j]) ;
						for(size_t k = 0 ;  k < inter.size() ; k++)
						{
							bool indescendants = false ;
							for(size_t l = 0 ; l < descendants.size() ; l++)
							{
								if(descendants[l]->inBoundary(inter[k]))
								{
									indescendants = true ;
									break ;
								}
							}
	
							if(!indescendants)
							{
								if(inRoot(inter[k]))
								{
									Point *p = new Point(inter[k]) ;
									additionalPoints.push_back(p) ;
									++count ;
									meshPoints.push_back(std::make_pair(p, tree[i])) ;
								}
							}
							
							if(count%100 == 0)
								std::cerr << "\r adding intersection points... " << count << std::flush ;
						}
					}
				}
			}
		}

		for(size_t i = 1 ;  i < tree.size() ; i++)
		{
			
			if(!tree[i]->isEnrichmentFeature && !tree[i]->isVirtualFeature && tree[0]->intersects(tree[i]))
			{
				std::vector<Point> inter = tree[0]->intersection(tree[i]) ;
				std::vector<Feature *> descendants = tree[i]->getDescendants() ;
				for(size_t k = 0 ;  k < inter.size() ; k++)
				{
					bool indescendants = false ;
					for(size_t l = 0 ; l < descendants.size() ; l++)
					{
						if(descendants[l]->inBoundary(inter[k]))
						{
							indescendants = true ;
							break ;
						}
					}
					
					if(!indescendants)
					{
						Point *p = new Point(inter[k]) ;
						additionalPoints.push_back(p) ;
						++count ;
						meshPoints.push_back(std::make_pair(p, tree[i])) ;
					}

				}
					
				if(count%100 == 0)
					std::cerr << "\r adding intersection points... " << count << std::flush ;
			}
	
		}
	}
	count = 0 ;
	std::cerr << " ...done." << std::endl ;
	
	//let us make sure we have no overlap
	std::stable_sort(meshPoints.begin(), meshPoints.end(), PairPointFeatureLess_Than_x()) ;
	std::stable_sort(meshPoints.begin(), meshPoints.end(), PairPointFeatureLess_Than_y()) ;
	std::stable_sort(meshPoints.begin(), meshPoints.end(), PairPointFeatureLess_Than_z()) ;
	std::deque<std::pair<Point *, Feature *> > ::iterator e = std::unique(meshPoints.begin(), meshPoints.end(), PairPointFeatureEqual());
	meshPoints.erase(e, meshPoints.end()) ;

	//shuffle for efficiency
	std::random_shuffle(meshPoints.begin(),meshPoints.end()) ;
// 	std::random_shuffle(meshPoints.begin(),meshPoints.begin()+bpcount) ;
// 	std::random_shuffle(meshPoints.begin()+bpcount, meshPoints.end()) ;
	
//	for(size_t i = 0 ; i < meshPoints.size() ; i++)
//		meshPoints[i].first->print() ;

	if(is2D())
	{	
		additionalPoints.push_back(new Point(bbox[0])) ;
		additionalPoints.push_back(new Point(bbox[2])) ;
		additionalPoints.push_back(new Point(bbox[4])) ;
		additionalPoints.push_back(new Point(bbox[6])) ;

		meshPoints.push_front(std::make_pair(additionalPoints[additionalPoints.size()-4], tree[0])) ;
		meshPoints.push_front(std::make_pair(additionalPoints[additionalPoints.size()-3], tree[0])) ;
		meshPoints.push_front(std::make_pair(additionalPoints[additionalPoints.size()-2], tree[0])) ;
		meshPoints.push_front(std::make_pair(additionalPoints[additionalPoints.size()-1], tree[0])) ;
			
		this->dtree = new DelaunayTree( meshPoints[0].first, meshPoints[1].first, meshPoints[2].first) ;
		this->dtree->insert(meshPoints[3].first) ;
		
		for( std::deque<std::pair<Point *, Feature *> >::iterator i = meshPoints.begin()+4 ; i != meshPoints.end(); ++i)
		{
//			if( (i - meshPoints.begin())%1000 == 0)
//				std::cerr << "\r generating triangles... point " << count << "/" << meshPoints.size() << std::flush ;
			
			++count ;
			
			if(*i->first != bbox[0] &&
			   *i->first != bbox[2] &&
			   *i->first != bbox[4] &&
			   *i->first != bbox[6] && inRoot(*i->first)
			  )
			{
				dtree->insert(i->first) ;
			}
		}
		
//		std::cerr << "\r generating triangles.... point " << meshPoints.size()-3 << "/" << meshPoints.size()-4 << " ...done" << std::endl ;
		
		bool correct = false ;
		int tries = correctionSteps ;
		
		while(!correct && tries)
		{
			std::vector< DelaunayTriangle * > tets = dtree->getTriangles();
			std::vector< Point *> to_insert ;
			
			for(size_t i = 0 ; i < tets.size() ;i++)
			{
				Point *test = checkElement(tets[i]);
				if(test)
				{
					to_insert.push_back(test);
				}
			}
			
			if(to_insert.empty())
				correct = true ;
			
			for(size_t i = 0 ; i < to_insert.size() ;i++)
			{
				std::cerr << "\r generating triangles.... point " << ++count << "/" << meshPoints.size()-3 << std::flush ;
				if(*to_insert[i] != bbox[0] &&
				   *to_insert[i] != bbox[2] &&
				   *to_insert[i] != bbox[4] &&
				   *to_insert[i] != bbox[6] &&
				   inRoot( *to_insert[i])
				  )
					dtree->insert(to_insert[i]) ;
				if(to_insert[i]->id == -1)
					delete to_insert[i] ;
			}
			
			tries-- ;
			
		}
		std::cerr << " ...done."<< std::endl ;
		
	}
	else if (is3D())
	{
		additionalPoints.push_back(new Point(bbox[0])) ;
		additionalPoints.push_back(new Point(bbox[1])) ;
		additionalPoints.push_back(new Point(bbox[7])) ;
		additionalPoints.push_back(new Point(bbox[2])) ;
		additionalPoints.push_back(new Point(bbox[3])) ;
		additionalPoints.push_back(new Point(bbox[4])) ;
		additionalPoints.push_back(new Point(bbox[5])) ;
		additionalPoints.push_back(new Point(bbox[6])) ;
		
		meshPoints.push_front(std::make_pair(additionalPoints[additionalPoints.size()-8], tree[0])) ;
		meshPoints.push_front(std::make_pair(additionalPoints[additionalPoints.size()-7], tree[0])) ;
		meshPoints.push_front(std::make_pair(additionalPoints[additionalPoints.size()-6], tree[0])) ;
		meshPoints.push_front(std::make_pair(additionalPoints[additionalPoints.size()-5], tree[0])) ;
		meshPoints.push_front(std::make_pair(additionalPoints[additionalPoints.size()-4], tree[0])) ;
		meshPoints.push_front(std::make_pair(additionalPoints[additionalPoints.size()-3], tree[0])) ;
		meshPoints.push_front(std::make_pair(additionalPoints[additionalPoints.size()-2], tree[0])) ;
		meshPoints.push_front(std::make_pair(additionalPoints[additionalPoints.size()-1], tree[0])) ;
			
		this->dtree3D = new DelaunayTree3D( meshPoints[0].first, meshPoints[1].first, meshPoints[2].first, meshPoints[3].first) ;
		this->dtree3D->insert(meshPoints[4].first) ;
		this->dtree3D->insert(meshPoints[5].first) ;
		this->dtree3D->insert(meshPoints[6].first) ;
		this->dtree3D->insert(meshPoints[7].first) ;

		assert(meshPoints[0].first->id > -1 ) ;
		assert(meshPoints[1].first->id > -1 ) ;
		assert(meshPoints[2].first->id > -1 ) ;
		assert(meshPoints[3].first->id > -1 ) ;
		
		
		for( std::deque<std::pair<Point *, Feature *> >::iterator i = meshPoints.begin()+8 ; i != this->meshPoints.end(); ++i)
		{
			if((i - meshPoints.begin())%1000 == 0)
				std::cerr << "\r generating tetrahedrons... point " << ++count*1000 << "/" << meshPoints.size()-8 << " tets: "<< dtree3D->tree.size() <<  std::flush ;
			if(*i->first != bbox[0] &&
			   *i->first != bbox[1] &&
			   *i->first != bbox[2] &&
			   *i->first != bbox[3] &&
			   *i->first != bbox[4] &&
			   *i->first != bbox[5] &&
			   *i->first != bbox[6] &&
			   *i->first != bbox[7]
			  )
			{
				dtree3D->insert(i->first) ; 
				if(i->first->id == -1)
				{
					std::cout << "insertion failed" << std::endl ;
				}
			}
		}
		
		for(size_t k  =  0 ; k <  enrichmentFeature.size() ; k++)
		{
			std::vector<Point *> pts = dynamic_cast<EnrichmentFeature *>(enrichmentFeature[k])->getSamplingPoints() ;
			
			for(size_t i = 0 ; i < pts.size() ;i++)
			{
				if(inRoot(*pts[i]))
				{
					meshPoints.push_back(std::pair<Point *, Feature *>(pts[i], this->tree[0])) ;
					std::cerr << "\r generating tetrahedrons.... point " << ++count << "/" << meshPoints.size()-3 << std::flush ;
					dtree3D->insert(pts[i]) ;
				}
			}
		}
		
		bool correct = false ;
		int tries = correctionSteps ;
		
		while(!correct && tries)
		{
			std::vector< DelaunayTetrahedron * > tets = dtree3D->getTetrahedrons();
			std::vector< Point *> to_insert ;
			
			for(size_t i = 0 ; i < tets.size() ;i++)
			{
				Point *test = checkElement(tets[i]);
				if(test)
				{
					to_insert.push_back(test);
				}
			}
			
			if(to_insert.empty())
				correct = true ;
			
			for(size_t i = 0 ; i < to_insert.size() ;i++)
			{
				std::cerr << "\r generating tetrahedrons.... point " << ++count << "/" << meshPoints.size()-3 << std::flush ;
				if(*to_insert[i] != bbox[0] &&
							*to_insert[i] != bbox[1] &&
							*to_insert[i] != bbox[2] &&
							*to_insert[i] != bbox[3] &&
							*to_insert[i] != bbox[4] &&
							*to_insert[i] != bbox[5] &&
							*to_insert[i] != bbox[6] &&
							*to_insert[i] != bbox[7] &&
							inRoot( *to_insert[i])
				)
					dtree3D->insert(to_insert[i]) ;
				
				if(to_insert[i]->id == -1)
					delete to_insert[i] ;
			}
			
			tries-- ;
			
		}
		std::cerr << " ...done."<< std::endl ;
// 		dtree3D->purge() ;
		
	}
}

BoundaryCondition::BoundaryCondition(LagrangeMultiplierType t, const double & d) : condition(t), data(d) { } ;



ProjectionDefinedBoundaryCondition::ProjectionDefinedBoundaryCondition(LagrangeMultiplierType t, const Point & dir, double d) : BoundaryCondition(t,d), direction(dir) { }

void ProjectionDefinedBoundaryCondition::apply(Assembly * a, DelaunayTree * t) const
{
	std::vector<DelaunayTriangle *> tris = t->getTriangles() ;
	for(size_t i = 0 ; i < tris.size() ; i++)
	{
		DelaunayTreeItem * VoidNeighbour ;
		bool border = false ;
		for(size_t j = 0 ; j < tris[i]->neighbour.size() ; j++)
		{
			bool voidNeighbour = (tris[i]->getNeighbour(j)->isTriangle 
			                    && dynamic_cast<DelaunayTriangle *>(tris[i]->getNeighbour(j))->getBehaviour()->type == VOID_BEHAVIOUR) ;
			border = border || tris[i]->getNeighbour(j)->isPlane
			                || voidNeighbour ;
			if(voidNeighbour)
				VoidNeighbour = tris[i]->getNeighbour(j) ;
			if(tris[i]->getNeighbour(j)->isPlane)
				VoidNeighbour = tris[i]->getNeighbour(j) ;
		}
		
		if(tris[i]->getBehaviour()->type == VOID_BEHAVIOUR)
			border = false ;
		
		if(border)
		{
			Segment ray((tris[i]->getCenter()), (tris[i]->getCenter())-direction*(tris[i]->getRadius())) ;
			bool isOnTheRightSide = true ;
			for(size_t j = 0 ; j < tris[i]->neighbour.size() ; j++)
			{
				if(tris[i]->getNeighbour(j)->isTriangle && !(dynamic_cast<DelaunayTriangle *>(tris[i]->getNeighbour(j))->getBehaviour()->type == VOID_BEHAVIOUR))
				{
					isOnTheRightSide = isOnTheRightSide && !(ray.intersects(dynamic_cast<Triangle *>(tris[i]->getNeighbour(j)))) ;
				}
			}
			
			if(isOnTheRightSide)
			{
				std::vector<Point> id ;
				for(size_t j = 0 ; j < tris[i]->getBoundingPoints().size() ; j++)
				{
					
					Segment side(tris[i]->getBoundingPoint(j), tris[i]->getBoundingPoint((j+1)%tris[i]->getBoundingPoints().size())) ;
					if(side.intersects(ray))
					{
						id.push_back(tris[i]->getBoundingPoint(j)) ;
						id.push_back(tris[i]->getBoundingPoint((j+1)%tris[i]->getBoundingPoints().size())) ;
					}
				}
				if(!id.empty())
				{
					apply2DBC(tris[i], id, condition, data, a) ;
				}
			}
		}
	}
}

void ProjectionDefinedBoundaryCondition::apply(Assembly * a, DelaunayTree3D * t)  const
{
	std::vector<DelaunayTetrahedron *> tris = t->getTetrahedrons() ;
	for(size_t i = 0 ; i < tris.size() ; i++)
	{
		std::vector<DelaunayDemiSpace *> space ;
		for(size_t j = 0 ; j < tris[i]->neighbour.size() ; j++)
		{
			if(tris[i]->getNeighbour(j)->isSpace())
			{
				space.push_back(static_cast<DelaunayDemiSpace *>(tris[i]->getNeighbour(j))) ;
			}
		}
		std::vector<Point> id ;
		for(size_t s = 0 ; s < space.size() ; s++)
		{
			Segment ray(tris[i]->getCenter(), tris[i]->getCenter()-direction*2.*tris[i]->getRadius()) ;
			std::vector<Point *> points = space[s]->commonSurface(tris[i]) ;
			TriPoint surf(points[0], points[1], points[2]) ;
			for(size_t j = 3 ; j < tris[i]->getBoundingPoints().size() ; j++)
			{
				if(isCoplanar(*points[0], *points[1], *points[2], tris[i]->getBoundingPoint(j)))
					points.push_back(&tris[i]->getBoundingPoint(j)) ;
			}
			
			if(ray.intersects(&surf) && !points.empty())
			{
				for(size_t j = 0 ; j < points.size() ; j++)
					id.push_back(*points[j]) ;
			}
		}
		apply3DBC(tris[i], id, condition, data, a) ;
	}
}


std::vector<DelaunayTriangle *> FeatureTree::getTriangles()
{
	if(dtree == NULL && dtree3D == NULL)
		this->generateElements() ;
	
	if(dtree != NULL)
	{
		if(!stitched)
			stitch() ;
		
		if(!setBehaviours)
			setElementBehaviours() ;
		
		if(!renumbered)
			renumber() ;

		return dtree->getTriangles() ;
	}
	else
		return std::vector<DelaunayTriangle *>(0) ;
}

std::vector<DelaunayTetrahedron *> FeatureTree::getTetrahedrons()
{
	if(dtree3D == NULL && dtree == NULL)
		this->generateElements() ;
	
	if(dtree3D != NULL)
	{
		
		if(!stitched)
			stitch() ;

		if(!setBehaviours)
			setElementBehaviours() ;
		
		if(!renumbered)
			renumber() ;
		
		if(!elements3D.empty())
			return elements3D ;
		
		return dtree3D->getTetrahedrons() ;
	}
	else
		return std::vector<DelaunayTetrahedron *>(0) ;
}

std::vector<DelaunayTriangle *> Feature::getBoundingTriangles( DelaunayTree * dt)
{
	std::vector<DelaunayTriangle *> tri = dt->conflicts(dynamic_cast<Geometry *>(this)) ;
	
	std::vector<DelaunayTriangle *> ret  ;
	
	for(size_t i = 0 ; i < tri.size() ; i++)
	{
		if(tri[i]->Triangle::intersects(dynamic_cast<Geometry *>(this)))
			ret.push_back(tri[i]) ;
	}
	
	return ret ;

}

std::vector<DelaunayTetrahedron *> Feature::getBoundingTetrahedrons( DelaunayTree3D * dt)
{
	std::vector<DelaunayTetrahedron *> tri = dt->conflicts(dynamic_cast<Geometry *>(this)) ;
	
	std::vector<DelaunayTetrahedron *> ret  ;
	
	for(size_t i = 0 ; i < tri.size() ; i++)
	{
		if(tri[i]->Tetrahedron::intersects(dynamic_cast<Geometry *>(this)))
			ret.push_back(tri[i]) ;
	}
	
	return ret ;
	
}


Voxel::Voxel()
{
	tlf.print() ;
	filled = false ;
}

Voxel::Voxel(double x, double y, double z ,double s) : tlf(x-s*.5, y+s*.5, z+s*.5), trf(x+s*.5, y+s*.5, z+s*.5), blf(x-s*.5, y-s*.5, z+s*.5), brf(x+s*.5, y-s*.5, z+s*.5),tlb(x-s*.5, y+s*.5, z-s*.5), trb(x+s*.5, y+s*.5, z-s*.5), blb(x-s*.5, y-s*.5, z-s*.5), brb(x+s*.5, y-s*.5, z-s*.5), filled(false)
{
}

const std::vector<Feature *> & Voxel::getFeatures() const
{
	return this->features ;
}

std::vector<Feature *> & Voxel::getFeatures()
{
	return this->features ;
}

bool Voxel::in(const Point & p) const
{
	return (p.x >= tlf.x)  && (p.x <= brf.x) && (p.y >= brf.y) && (p.y <= tlf.y) && (p.z >= brb.z) && (p.z <= tlf.z);
}

void Voxel::coOccuringFeatures(std::vector<Feature *> &f , const Geometry * inc) const
{
	if(pixels.empty())
	{
		f.insert(f.end(), features.begin(), features.end()) ;
		return ;
	}
	
	for(size_t i = 0 ; i < 8 ; i++)
	{
		if(pixels[i]->coOccur(inc))
		{
			pixels[i]->coOccuringFeatures(f, inc) ;
		}
	}
}

void Voxel::coOccuringFeatures(std::vector<Feature *> &f , const Point & p) const
{
	if(pixels.empty())
	{
		f.insert(f.end(), features.begin(), features.end()) ;
		return ;
	}
	
	std::vector<Feature *> ret ;
	
	for(size_t i = 0 ; i < 8 ; i++)
	{
		if(pixels[i]->coOccur(p))
		{
			pixels[i]->coOccuringFeatures(f, p) ;
		}
	}
}

bool Voxel::coOccur(const Geometry * inc) const
{
	return inc->in(tlf) 
		|| inc->in(trf) 
		|| inc->in(brf) 
		|| inc->in(blf) 
		|| inc->in(tlb) 
		|| inc->in(trb) 
		|| inc->in(brb) 
		|| inc->in(blb) 
		|| in(inc->getCenter()+Point(inc->getRadius(), 0, inc->getRadius()))
		|| in(inc->getCenter()+Point(-inc->getRadius(), 0, inc->getRadius())) 
		|| in(inc->getCenter()+Point(0,inc->getRadius(), inc->getRadius())) 
		|| in(inc->getCenter()+Point(0,-inc->getRadius(), inc->getRadius())) 
		|| in(inc->getCenter()+Point(inc->getRadius(), 0, -inc->getRadius()))
		|| in(inc->getCenter()+Point(-inc->getRadius(), 0, -inc->getRadius())) 
		|| in(inc->getCenter()+Point(0,inc->getRadius(), -inc->getRadius())) 
		|| in(inc->getCenter()+Point(0,-inc->getRadius(), -inc->getRadius()))
		|| in(inc->getCenter()+Point(inc->getRadius(),inc->getRadius(), inc->getRadius())) 
		|| in(inc->getCenter()+Point(inc->getRadius(),inc->getRadius(), -inc->getRadius())) 
		|| in(inc->getCenter()+Point(inc->getRadius(), -inc->getRadius(), inc->getRadius()))
		|| in(inc->getCenter()+Point(inc->getRadius(), -inc->getRadius(), -inc->getRadius())) 
		|| in(inc->getCenter()+Point(-inc->getRadius(),inc->getRadius(), inc->getRadius())) 
		|| in(inc->getCenter()+Point(-inc->getRadius(),-inc->getRadius(), -inc->getRadius()))
		|| in(inc->getCenter());
}

void Voxel::remove(Feature * inc)
{

	std::vector<Feature *>::iterator e = std::find(features.begin(), features.end(), inc) ;
	if(e != features.end())
		features.erase(e) ;
	
	if(!pixels.empty())
		for(size_t i = 0 ; i < 8 ; i++)
			pixels[i]->remove(inc) ;
	
	if(features.empty() && !pixels.empty())
	{
		pixels.clear() ;
	}
	
	filled = false ;
}

bool Voxel::isFilled() const
{
	return filled ;
}

Point Voxel::center() const
{
	return (tlf + brb)*.5 ;
}

void Voxel::refine()
{
	if(pixels.empty())
	{
		pixels.push_back(new Voxel(tlf.x+(trf.x-tlf.x)*.25, blf.y+(tlf.y-blf.y)*.75, tlb.z+(trf.z-tlb.z)*.25, (trf.x-tlf.x)*.5)) ;
		pixels.push_back(new Voxel(tlf.x+(trf.x-tlf.x)*.75, blf.y+(tlf.y-blf.y)*.75, tlb.z+(trf.z-tlb.z)*.25, (trf.x-tlf.x)*.5)) ;
		pixels.push_back(new Voxel(tlf.x+(trf.x-tlf.x)*.25, blf.y+(tlf.y-blf.y)*.25, tlb.z+(trf.z-tlb.z)*.25, (trf.x-tlf.x)*.5)) ;
		pixels.push_back(new Voxel(tlf.x+(trf.x-tlf.x)*.75, blf.y+(tlf.y-blf.y)*.25, tlb.z+(trf.z-tlb.z)*.25, (trf.x-tlf.x)*.5)) ;
		pixels.push_back(new Voxel(tlf.x+(trf.x-tlf.x)*.25, blf.y+(tlf.y-blf.y)*.75, tlb.z+(trf.z-tlb.z)*.75, (trf.x-tlf.x)*.5)) ;
		pixels.push_back(new Voxel(tlf.x+(trf.x-tlf.x)*.75, blf.y+(tlf.y-blf.y)*.75, tlb.z+(trf.z-tlb.z)*.75, (trf.x-tlf.x)*.5)) ;
		pixels.push_back(new Voxel(tlf.x+(trf.x-tlf.x)*.25, blf.y+(tlf.y-blf.y)*.25, tlb.z+(trf.z-tlb.z)*.75, (trf.x-tlf.x)*.5)) ;
		pixels.push_back(new Voxel(tlf.x+(trf.x-tlf.x)*.75, blf.y+(tlf.y-blf.y)*.25, tlb.z+(trf.z-tlb.z)*.75, (trf.x-tlf.x)*.5)) ;
	}
	else
		return ;
	
	for(size_t i = 0 ; i < 8 ; i++)
	{
		for(size_t j = 0 ; j < features.size() ; j++)
			if(pixels[i]->coOccur(features[j]))
				pixels[i]->forceAdd(features[j]) ;
	}

	features.clear() ;
}

bool Voxel::add(Feature * inc)
{
	
	if(filled)
		return false;
	
	if(!pixels.empty())
	{
		bool ret = true ;
		#pragma omp parallel for
		for(int i = 0 ; i < 8 ; i++)
			if(pixels[i]->coOccur(inc))
				ret = ret && pixels[i]->add(inc) ;
		
		if(!ret)
			for(size_t i = 0 ; i < 8 ; i++)
				pixels[i]->remove(inc) ;
		
		if(ret)
			forceAdd(inc) ;
		return ret ;
	}
	
	if(!features.empty())
	{
		for(size_t i = 0 ; i < this->features.size() ; i++)
		{
			if(this->features[i]->intersects(inc))
				return false;
		}
		
		this->features.push_back(inc) ;
		
		if(features.size() > 64 )
			refine() ;
		
		return true;
	}
	else
	{
		if(inc->in(tlf) 
			&& inc->in(trf) 
			&& inc->in(brf) 
			&& inc->in(blf)
			&& inc->in(tlb) 
			&& inc->in(trb) 
			&& inc->in(brb) 
			&& inc->in(blb)
			)
			filled = true ;
		
		this->features.push_back(inc) ;
		
		return true ;
	}
}

void Voxel::forceAdd(Feature * inc)
{
	if(pixels.empty())
		this->features.push_back(inc) ;
	
	if(pixels.empty())
		if(inc->in(tlf) 
		&& inc->in(trf) 
		&& inc->in(brf) 
		&& inc->in(blf)
		&& inc->in(tlb) 
		&& inc->in(trb) 
		&& inc->in(brb) 
		&& inc->in(blb)
		)
			filled = true ;
	
	if(!pixels.empty())
	{
		#pragma omp parallel for
		for(int i = 0 ; i < 8 ; i++)
		{
			if(pixels[i]->coOccur(inc))
				pixels[i]->forceAdd(inc) ;
		}
	}
	
	if(features.size() > 64 )
	{
		refine() ;
	}
}

void Voxel::print() const
{
	for(size_t i = 0 ; i < features.size() ; i++)
		features[i]->print() ;
}

Voxel::~Voxel()
{
	if(!pixels.empty())
		for(size_t i = 0 ; i < 8 ; i++)
			delete pixels[i] ;
}

bool Voxel::coOccur(const Point & p) const
{
	return p.x > tlf.x &&  p.x < trf.x && p.y > blf.y &&  p.y < tlf.y && p.z > blb.z &&  p.z < tlf.z;
}

int Voxel::computeFillFactor() const
{
	if(!pixels.empty())
		return 0 ;
	
	int count  = 0 ;
	for(size_t i = 0 ; i < 10 ; i++)
	{
		for(size_t j = 0 ; j < 10 ; j++)
		{
			for(size_t k = 0 ; k < 10 ; k++)
			{
				Point test(blf.x+(double)(i+1)*.1*(brf.x-blf.x), blf.y+(double)(j+1)*.1*(tlf.y-blf.y), blf.y+(double)(k+1)*.1*(tlf.y-blb.y)) ;
				
				for(size_t l = 0 ; l < features.size() ; l++)
				{
					if(features[l]->in(test))
					{
						count++ ;
						break ;
					}
				}
			}
		}
	}
	
	return count ;
}

Pixel::Pixel() : pixels(0)
{
	filled = false ;
}

Pixel::Pixel(double x, double y, double s) : tl(x-s*.5, y+s*.5), tr(x+s*.5, y+s*.5), bl(x-s*.5, y-s*.5), br(x+s*.5, y-s*.5), filled(false) {} ;

const std::vector<Feature *> & Pixel::getFeatures() const
{
	return this->features ;
}

std::vector<Feature *> & Pixel::getFeatures()
{
	return this->features ;
}

bool Pixel::in(const Point & p) const
{
	return (p.x >= tl.x)  && (p.x <= br.x) && (p.y >= br.y) && (p.y <= tl.y);
}

void Pixel::refine()
{
	
	if(pixels.size() == 0)
	{
		pixels.resize(4) ;
		pixels[0] = Pixel(tl.x+(tr.x-tl.x)*.25, bl.y+(tl.y-bl.y)*.75, (tr.x-tl.x)*.5) ;
		pixels[1] = Pixel(tl.x+(tr.x-tl.x)*.75, bl.y+(tl.y-bl.y)*.75, (tr.x-tl.x)*.5) ;
		pixels[2] = Pixel(tl.x+(tr.x-tl.x)*.25, bl.y+(tl.y-bl.y)*.25, (tr.x-tl.x)*.5) ;
		pixels[3] = Pixel(tl.x+(tr.x-tl.x)*.75, bl.y+(tl.y-bl.y)*.25, (tr.x-tl.x)*.5) ;
		
		for(size_t i = 0 ; i < 4 ; i++)
		{
			for(size_t j = 0 ; j < features.size() ; j++)
				if(pixels[i].coOccur(features[j]))
					pixels[i].forceAdd(features[j]) ;
		}
	}
	else
		return ;
}

Pixel::~Pixel()
{
}

bool Pixel::coOccur(const Geometry * inc) const
{
	std::vector<Point> bbox = inc->getBoundingBox() ;
	return inc->in(tl) 
		|| inc->in(tr) 
		|| inc->in(br) 
		|| inc->in(bl) 
		|| inc->in(tr*.5+bl*.5)
		|| in(inc->getCenter()+Point(inc->getRadius(), 0))
		|| in(inc->getCenter())
		|| in(inc->getCenter()+Point(-inc->getRadius(), 0)) 
		|| in(inc->getCenter()+Point(0,inc->getRadius())) 
		|| in(inc->getCenter()+Point(0,-inc->getRadius()))
		|| in(bbox[0])
		|| in(bbox[1])
		|| in(bbox[2]) 
		|| in(bbox[3])  ;
}

bool Pixel::coOccur(const Point & p) const
{
	return p.x > tl.x &&  p.x < tr.x && p.y > bl.y &&  p.y < tl.y ;
}

void Pixel::remove(Feature * inc)
{
	std::vector<Feature *>::iterator e = std::find(features.begin(), features.end(), inc) ;
	if(e != features.end())
		features.erase(e) ;
	
	if(pixels.size())
		for(size_t i = 0 ; i < 4 ; i++)
			pixels[i].remove(inc) ;
	
	if(features.empty() && pixels.size())
	{
		pixels.resize(0) ;
	}
	
	filled = false ;
}

int Pixel::computeFillFactor() const
{
	if(!pixels.size())
		return 0 ;
	
	int count  = 0 ;
	for(size_t i = 0 ; i < 10 ; i++)
	{
		for(size_t j = 0 ; j < 10 ; j++)
		{
			Point test(bl.x+(double)(i+1)*.1*(br.x-bl.x), bl.y+(double)(j+1)*.1*(tl.y-bl.y)) ;
			
			for(size_t k = 0 ; k < features.size() ; k++)
			{
				if(features[k]->in(test))
				{
					count++ ;
					break ;
				}
			}
		}
	}
	
	return count ;
}

void Pixel::coOccuringFeatures(std::vector<Feature *> &f , const Geometry * inc) const
{
	if(!pixels.size() && !features.empty())
	{
		f.insert(f.end(), features.begin(), features.end()) ;
		return ;
	}
	
	std::vector<Feature *> ret ;

	
	if(pixels.size())
	{
		for(size_t i = 0 ; i < 4 ; i++)
		{
			if(pixels[i].coOccur(inc))
			{
				pixels[i].coOccuringFeatures(f, inc) ;
			}
		}
	}
}

void Pixel::coOccuringFeatures(std::vector<Feature *> &f , const Point & p) const
{
	if(!pixels.size() && !features.empty())
	{
		f.insert(f.end(), features.begin(), features.end()) ;
		return ;
	}
	
	std::vector<Feature *> ret ;
	
	if(pixels.size())
	{
		for(size_t i = 0 ; i < 4 ; i++)
		{
			if(pixels[i].coOccur(p))
			{
				pixels[i].coOccuringFeatures(f, p) ;
			}
		}
	}
}

bool Pixel::add(Feature * inc)
{
	if(filled)
		return false;
	
	if(pixels.size())
	{
		bool ret = true ;
		for(size_t i = 0 ; i < 4 ; i++)
		{
			if(pixels[i].coOccur(inc))
				ret = ret && pixels[i].add(inc) ;
		}
		
		if(!ret)
			for(size_t i = 0 ; i < 4 ; i++)
				pixels[i].remove(inc) ;
		
		if(ret)
			forceAdd(inc) ;
		
		return ret ;
	}

	if(!features.empty())
	{
		for(size_t i = 0 ; i < this->features.size() ; i++)
		{
			if(this->features[i]->intersects(inc))
				return false;
		}
		this->features.push_back(inc) ;
		
		if(features.size() > 32 )
			refine() ;
		
		return true;
	}
	else
	{
		if(inc->in(tl) && inc->in(tr) && inc->in(br) && inc->in(bl))
			filled = true ;
		
		this->features.push_back(inc) ;
		
		return true ;
	}
}

void Pixel::forceAdd(Feature * inc)
{	
	this->features.push_back(inc) ;
	
	if(inc->in(tl) && inc->in(tr) && inc->in(br) && inc->in(bl))
		filled = true ;
	
	if(pixels.size())
	{
		for(size_t i = 0 ; i < 4 ; i++)
			if(pixels[i].coOccur(inc))
				pixels[i].forceAdd(inc) ;
	}
	

	
	if(!filled && features.size() > 32 )
	{
		refine() ;
	}
}

void Pixel::print() const
{
	for(size_t i = 0 ; i < features.size() ; i++)
		features[i]->print() ;
}

Grid3D::Grid3D(double sizeX, double sizeY, double sizeZ, int div, const Point & center ): x(std::abs(sizeX)), y(std::abs(sizeY)) , z(std::abs(sizeZ)), c(center)
{
	dirtyCounter = 0;
	if(x < y && x < z)
	{
		lengthX = div ;
		lengthY = div*(y/x) ;
		lengthZ = div*(z/x) ;
		pixels.resize(lengthX,std::valarray<std::valarray<Voxel *> >(std::valarray<Voxel *>((Voxel *)NULL,lengthZ),lengthY)) ;
	}
	else if(y < x && y < z)
	{
		lengthY = div ;
		lengthX = div*(x/y) ;
		lengthZ = div*(z/y) ;
		pixels.resize(lengthX,std::valarray<std::valarray<Voxel *> >(std::valarray<Voxel *>((Voxel *)NULL,lengthZ),lengthY)) ;
	}
	else if(z < x && z < y)
	{
		lengthZ = div ;
		lengthX = div*(x/z) ;
		lengthY = div*(y/z) ;
		pixels.resize(lengthX,std::valarray<std::valarray<Voxel *> >(std::valarray<Voxel *>((Voxel *)NULL,lengthZ),lengthY)) ;
	}
	else
	{
		lengthX = div ;
		lengthY = div ;
		lengthZ = div ;
		pixels.resize(lengthX,std::valarray<std::valarray<Voxel *> >(std::valarray<Voxel *>((Voxel *)NULL,lengthZ),lengthY)) ;
	}
	
	psize = x/lengthX;
	for(size_t i = 0 ; i < lengthX ; i++)
	{
		for(size_t j = 0 ; j < lengthY ; j++)
		{
			for(size_t k = 0 ; k < lengthZ ; k++)
			{
				
// 				pixels[i][j][k] = new Voxel(x*(double)(i)/(double)lengthX+psize*.5+center.x,
// 				                            y*(double)(j)/(double)lengthY+psize*.5+center.y,
// 				                            z*(double)(k)/(double)lengthZ+psize*.5+center.z,  psize) ;
// 				freepixel.push_back(pixels[i][j][k]) ;
				pixels[i][j][k] = new Voxel(x*(double)(i)/(double)lengthX+psize*.5+center.x-x/2.,
				                            y*(double)(j)/(double)lengthY+psize*.5+center.y-y/2.,
				                            z*(double)(k)/(double)lengthZ+psize*.5+center.z-z/2.,  psize) ;
				unfilledpixel.push_back(pixels[i][j][k]) ;
			}
		}
	}
	
// 	std::sort(freepixel.begin(), freepixel.end()) ;
	std::sort(unfilledpixel.begin(), unfilledpixel.end()) ;
}

Point Grid3D::randomFreeCenter() const 
{
	return Point(x*((double)rand()/(RAND_MAX)-.5)+c.x, 
	             y*((double)rand()/(RAND_MAX)-.5)+c.y, 
	             z*((double)rand()/(RAND_MAX)-.5)+c.z) ;
}

Grid3D Grid3D::getGrid(int div) const
{
	Hexahedron all(x, y, z, c) ;
	std::vector<Feature *> features = coOccur(&all) ;
	Grid3D ret(x,y, z,div,c) ;
	for(size_t i = 0 ; i < features.size() ; i++)
	{
		ret.forceAdd(features[i]) ;
	}
	
	return ret ;
}

// std::vector<HexahedralElement *> elementsFromGrid(int div)
// {
// 	Grid3D newGrid(div) ;
// }

Grid3D::~Grid3D()
{
	for(size_t i = 0 ; i < lengthX ; i++)
	{
		for(size_t j = 0 ; j < lengthY ; j++)
		{
			for(size_t k = 0 ; k < lengthZ ; k++)
			{
				delete pixels[i][j][k] ;
			}
		}
	}
}

double Grid3D::fraction() const
{
	return (double)dirtyCounter/(lengthX*lengthY*lengthZ) ;
}

bool Grid3D::add(Feature * inc)
{
	
	std::vector<Feature *> toTest = coOccur(inc);
	for(size_t i = 0 ; i < toTest.size() ; i++)
		if(inc->intersects(toTest[i]))
			return false ;
	
	forceAdd(inc) ;
	
	return true ;
	
// 	bool ret = true ;
// 	std::vector<Voxel *> cleanup ;
// 	
// 	double startX = .5*x + inc->getCenter().x-inc->getRadius() ;
// 	int startI = std::max(0., startX/psize - 2) ;
// 	
// 	double endX =  startX+2.*inc->getRadius();
// 	int endI = std::min(endX/psize + 2, (double)lengthX);
// 	
// 	double startY = .5*y + inc->getCenter().y-inc->getRadius() ;
// 	int startJ = std::max(0., startY/psize - 2) ;
// 	
// 	double endY =  startY+2.*inc->getRadius();
// 	int endJ = std::min(endY/psize + 2, (double)lengthY);
// 	
// 	double startZ = .5*y + inc->getCenter().z-inc->getRadius() ;
// 	int startK = std::max(0., startZ/psize - 2) ;
// 	
// 	double endZ =  startZ+2.*inc->getRadius();
// 	int endK = std::min(endZ/psize + 2, (double)lengthZ);
// 	
// 	for(int i = startI ; i < endI ; i++)
// 	{
// 		for(int j = startJ ; j < endJ ; j++)
// 		{
// 			for(int k = startK ; k < endK ; k++)
// 			{
// 
// 				if(pixels[i][j][k]->coOccur(inc))
// 				{
// 					if(pixels[i][j][k]->add(inc))
// 					{
// 						cleanup.push_back(pixels[i][j][k]) ;
// 					}
// 					else
// 					{
// 						for(size_t l = 0 ; l < cleanup.size() ; l++)
// 						{
// 							cleanup[l]->remove(inc) ;
// 						}
// 						return false ;
// 					}
// 					
// 				}
// 			}
// 		}
// 		
// 	}
// 	
// 	return ret ;
}

void Grid3D::forceAdd(Feature * inc)
{
	double startX = .5*x-c.x + inc->getCenter().x-inc->getRadius() ;
	int startI = std::max(0., startX/psize - 2) ;
	
	double endX =  startX+2.*inc->getRadius();
	int endI = std::min(endX/psize + 2, (double)lengthX);
	
	double startY = .5*y-c.y + inc->getCenter().y-inc->getRadius() ;
	int startJ = std::max(0., startY/psize - 2) ;
	
	double endY =  startY+2.*inc->getRadius();
	int endJ = std::min(endY/psize + 2, (double)lengthY);
	
	double startZ = .5*z-c.z + inc->getCenter().z-inc->getRadius() ;
	int startK = std::max(0., startZ/psize - 2) ;
	
	double endZ =  startZ+2.*inc->getRadius();
	int endK = std::min(endZ/psize + 2, (double)lengthZ);
	std::vector<Voxel *> cleanup ;
	for(int i = startI ; i < endI ; i++)
	{
		for(int j = startJ ; j < endJ ; j++)
		{
			for(int k = startK ; k < endK ; k++)
			{
// 	for(int i = 0 ; i < lengthX ; i++)
// 	{
// 		for(int j = 0 ; j < lengthY ; j++)
// 		{
// 			for(int k = 0 ; k < lengthZ ; k++)
// 			{
				if(pixels[i][j][k]->coOccur(inc))
				{
					pixels[i][j][k]->forceAdd(inc) ;
				}
			}
		}
	}

}

std::vector<Feature *> Grid3D::coOccur(const Geometry * geo) const
{
	std::vector<Feature *> ret ;
	double startX = .5*x-c.x + geo->getCenter().x-geo->getRadius() ;
	int startI = std::max(0., startX/psize - 2) ;
	
	double endX =  startX+2.*geo->getRadius();
	int endI = std::min(endX/psize + 2, (double)lengthX);
	
	double startY = .5*y-c.y + geo->getCenter().y-geo->getRadius() ;
	int startJ = std::max(0., startY/psize - 2) ;
	
	double endY =  startY+2.*geo->getRadius();
	int endJ = std::min(endY/psize + 2, (double)lengthY);
	
	double startZ = .5*z-c.z + geo->getCenter().z-geo->getRadius() ;
	int startK = std::max(0., startZ/psize - 2) ;
	
	double endZ =  startZ+2.*geo->getRadius();
	int endK = std::min(endZ/psize + 2, (double)lengthZ);

	if(geo->getGeometryType() == TETRAHEDRON)
	{
		const Tetrahedron * t = dynamic_cast<const Tetrahedron *>(geo) ;
		startX = .5*x + t->getCircumCenter()->x-t->getRadius()*1.1 ;
		startI = std::max(0., startX/psize - 2) ;
		
		endX =  startX+2.2*geo->getRadius();
		endI = std::min(endX/psize + 2, (double)lengthX);
		
		startY = .5*y + t->getCircumCenter()->y-t->getRadius()*1.1 ;
		startJ = std::max(0., startY/psize - 2) ;
		
		endY =  startY+2.2*t->getRadius();
		endJ = std::min(endY/psize + 2, (double)lengthY);
		
		startZ = .5*z + t->getCircumCenter()->z-t->getRadius()*1.1 ;
		startK = std::max(0., startZ/psize - 2) ;
		
		endZ =  startZ+2.2*t->getRadius();
		endK = std::min(endZ/psize + 2, (double)lengthZ);
	}
	
	bool foundPixel = false ;
	for(int i = startI ; i < endI ; i++)
	{
		for(int j = startJ ; j < endJ ; j++)
		{
			for(int k = startK ; k < endK ; k++)
// 	for(int i = 0 ; i < lengthX ; i++)
	{
// 		for(int j = 0 ; j < lengthY ; j++)
// 		{
// 			for(int k = 0 ; k < lengthZ ; k++)
// 			{
				if(pixels[i][j][k]->coOccur(geo))
				{
					foundPixel = true ;
					pixels[i][j][k]->coOccuringFeatures(ret,geo) ;
					std::stable_sort(ret.begin(), ret.end());
					std::vector<Feature *>::iterator e = std::unique(ret.begin(), ret.end()) ;
				}
			}
		}
	}
	
	std::stable_sort(ret.begin(), ret.end());
	std::vector<Feature *>::iterator e = std::unique(ret.begin(), ret.end()) ;
	ret.erase(e, ret.end()) ;
	return ret ;
}

std::vector<Feature *> Grid3D::coOccur(const Point & p) const 
{
	std::vector<Feature *> ret ;
// 	double startX = .5*x + p.x-.05*x ;
// 	int startI = std::max(0., startX/psize - 2) ;
// 	
// 	double endX =  startX+.1*x;
// 	int endI = std::min(endX/psize + 2, (double)lengthX);
// 	
// 	double startY = .5*y + p.y-.05*y ;
// 	int startJ = std::max(0., startY/psize - 2) ;
// 	
// 	double endY =  startY+.1*y;
// 	int endJ = std::min(endY/psize + 2, (double)lengthY);
// 	
// 	double startZ = .5*z + p.z-.05*z ;
// 	int startK = std::max(0., startZ/psize - 2) ;
// 	
// 	double endZ =  startZ+.1*z;
// 	int endK = std::min(endZ/psize + 2, (double)lengthZ);
// 	
	for(int i = std::max((p.x-c.x-psize)/psize, 0.) ; i < std::min((int)((p.x-c.x+psize)/psize+1),(int)pixels.size()) ; i++)
	{
		for(int j = std::max((p.y-c.y-psize)/psize,0.) ; 
		    j < std::min((int)((p.y-c.y+psize)/psize+1),(int)pixels[i].size()) ; 
		    j++)
		{
			for(int k = std::max((p.z-c.z-psize)/psize,0.) ; 
			    k <std::min((int)((p.z-c.z+psize)/psize+1),(int)pixels[i][j].size()) ;
			    k++)
			{
				if(pixels[i][j][k]->coOccur(p))
				{
					pixels[i][j][k]->coOccuringFeatures(ret,p) ;
				}
			}
		}
	}
	
	std::stable_sort(ret.begin(), ret.end());
	std::vector<Feature *>::iterator e = std::unique(ret.begin(), ret.end()) ;
	ret.erase(e, ret.end()) ;
	return ret ;
}

Grid::~Grid()
{
	for(size_t i = 0 ; i < lengthX ; i++)
	{
		for(size_t j = 0 ; j < lengthY ; j++)
		{
			delete pixels[i][j] ;
		}
	}
}

Grid::Grid(double sizeX, double sizeY, int div, const Point & center ) : x(sizeX), y(sizeY), c(center)
{
	if(x>y)
	{
		
		lengthX = div*(x/y) ;
		lengthY = div ;
		pixels.resize(lengthX,std::valarray<Pixel *>((Pixel *)NULL,lengthY)) ;
	}
	else
	{
		lengthX = div ;
		lengthY = div*(y/x) ;
		pixels.resize(lengthX,std::valarray<Pixel *>((Pixel *)NULL,lengthY)) ;
	}
	
	psize = std::max(std::abs(x/lengthX), std::abs(y/lengthY));
	
	for(size_t i = 0 ; i < lengthX ; i++)
	{
		for(size_t j = 0 ; j < lengthY ; j++)
		{
			pixels[i][j] = new Pixel(x*(double)(i)/(double)lengthX+psize*.5+center.x-.5*x,
			                         y*(double)(j)/(double)lengthY+psize*.5+center.y-.5*y, psize) ;
		}
	}
}

std::vector<Feature *> Grid::coOccur(const Geometry * geo) const
{
	std::vector<Feature *> ret ;
	double startX = x*.5-c.x + geo->getCenter().x-geo->getRadius() ;
	int startI = std::max(0., startX/psize - 2) ;
	
	double endX =  startX+2.*geo->getRadius();
	int endI = std::min(endX/psize + 2, (double)lengthX);
	
	double startY = y*.5-c.y + geo->getCenter().y-geo->getRadius() ;
	int startJ = std::max(0., startY/psize - 2) ;
		
	double endY =  startY+2.*geo->getRadius();
	int endJ = std::min(endY/psize + 2, (double)lengthY);
	for(int i = startI ; i < endI ; i++)
	{
		for(int j = startJ ; j < endJ ; j++)
		{
			if(pixels[i][j]->coOccur(geo))
			{
				pixels[i][j]->coOccuringFeatures(ret,geo) ;
				std::stable_sort(ret.begin(), ret.end());
				std::vector<Feature *>::iterator e = std::unique(ret.begin(), ret.end()) ;
				ret.erase(e, ret.end()) ;
			}
		}
	}

	std::stable_sort(ret.begin(), ret.end());
	std::vector<Feature *>::iterator e = std::unique(ret.begin(), ret.end()) ;
	ret.erase(e, ret.end()) ;
	return ret ;
}

 std::vector<Feature *> Grid::coOccur(const Point & p) const 
{

	std::vector<Feature *> ret ;
	double startX = x*.5-c.x + p.x ;
	int startI = std::max(0., startX/psize - 2) ;
	
	double endX =  startX+.05*x;
	int endI = std::min(endX/psize + 2, (double)lengthX);
	
	double startY = y*.5-c.y + p.y ;
	int startJ = std::max(0., startY/psize - 2) ;
	
	double endY =  startY+.05*y;
	int endJ = std::min(endY/psize + 2, (double)lengthY);
	
	for(int i = startI ; i < endI ; i++)
	{
		for(int j = startJ ; j < endJ ; j++)
		{
			if(pixels[i][j]->coOccur(p))
			{
				pixels[i][j]->coOccuringFeatures(ret,p) ;
			}
		}
	}
	
	std::stable_sort(ret.begin(), ret.end());
	std::vector<Feature *>::iterator e = std::unique(ret.begin(), ret.end()) ;
	ret.erase(e, ret.end()) ;
	return ret ;
}


void Grid::forceAdd(Feature * inc)
{
	
	double startX = x*.5-c.x + inc->getCenter().x-inc->getRadius() ;
	int startI = std::max(0., startX/psize - 2) ;
	
	double endX =  startX+2.*inc->getRadius();
	int endI = std::min(endX/psize + 2, (double)lengthX);
	
	double startY = y*.5-c.y + inc->getCenter().y-inc->getRadius() ;
	int startJ = std::max(0., startY/psize - 2) ;
	
	double endY =  startY+2.*inc->getRadius();
	int endJ = std::min(endY/psize + 2, (double)lengthY);

	for(int i = startI ; i < endI ; i++)
	{
		for(int j = startJ ; j < endJ ; j++)
		{
			if(pixels[i][j]->coOccur(inc))
			{
				pixels[i][j]->forceAdd(inc) ;
			}
		}
		
	}


}

Grid Grid::getGrid(int div) const
{
	Rectangle all(x, y, c) ;
	std::vector<Feature *> features = coOccur(&all) ;
	Grid ret(x,y,div,c) ;
	for(size_t i = 0 ; i < features.size() ; i++)
	{
		ret.forceAdd(features[i]) ;
	}
	
	return ret ;
}

bool Grid::add(Feature * inc)
{
	std::vector<Feature *> toTest = coOccur(inc);
//	std::cout << toTest.size() << std::endl ;
//	c.print() ;
//	inc->getCenter().print() ;
	for(size_t i = 0 ; i < toTest.size() ; i++)
	{
//		toTest[i]->getCenter().print() ;
//		std::cout << inc->intersects(toTest[i]) ;
		if(inc->intersects(toTest[i]))
			return false ;
	}
//	std::cout << std::endl ;
	forceAdd(inc) ;
	
	return true ;
	
}
