// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010

#include "boundarycondition.h"

using namespace Mu ;

BoundingBoxDefinedBoundaryCondition::BoundingBoxDefinedBoundaryCondition(LagrangeMultiplierType t, BoundingBoxPosition pos, double d) :BoundaryCondition(t, d), pos(pos) { } ;

BoundingBoxAndRestrictionDefinedBoundaryCondition::BoundingBoxAndRestrictionDefinedBoundaryCondition(LagrangeMultiplierType t, BoundingBoxPosition pos, double xm, double xp, double ym, double yp, double zm, double zp, double d) : BoundaryCondition(t, d), pos(pos),  xmin(xm), xmax(xp), ymin(ym), ymax(yp), zmin(zm), zmax(zp)
{

}

BoundingBoxNearestNodeDefinedBoundaryCondition::BoundingBoxNearestNodeDefinedBoundaryCondition(LagrangeMultiplierType t, BoundingBoxPosition pos, Point p, double d ) : BoundaryCondition(t, d), pos(pos), nearest(p) {} ;


BoundingBoxAndRestrictionDefinedBoundaryCondition::BoundingBoxAndRestrictionDefinedBoundaryCondition(LagrangeMultiplierType t, BoundingBoxPosition pos, double xm, double xp, double ym, double yp, double d): BoundaryCondition(t, d), pos(pos),  xmin(xm), xmax(xp), ymin(ym), ymax(yp), zmin(0), zmax(0)
{

}

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
				Vector imposed(3) ;
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
				v[0] = XI ; v[1] = ETA ; v[2] = ZETA ;
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
				v[0] = XI ; v[1] = ETA ; v[2] = ZETA ;
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
				v[0] = XI ; v[1] = ETA ; v[2] = ZETA ;
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
				v[0] = XI ; v[1] = ETA ; v[2] = ZETA ;
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
				v[0] = XI ; v[1] = ETA ; v[2] = ZETA ;
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
				v[0] = XI ; v[1] = ETA ; v[2] = ZETA ;
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

void BoundingBoxNearestNodeDefinedBoundaryCondition::apply(Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t) const 
{
	std::vector<ElementarySurface *> & elements = a->getElements2d() ;
	if(elements.empty())
	{
		std::cout << "no elements in assembly" << std::endl ;
		return ;
	}
	double minx = elements.front()->getBoundingPoint(0).x ;
	double miny = elements.front()->getBoundingPoint(0).y ;
	double maxx = elements.front()->getBoundingPoint(0).x ;
	double maxy = elements.front()->getBoundingPoint(0).y ; 
	for(size_t i = 0 ; i < elements.size() ; ++i)
	{
		if(elements[i]->getBehaviour()->type != VOID_BEHAVIOUR)
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
	}
	double tol = std::min(maxx-minx, maxy-miny)*.0001 ;
	std::map<double, std::pair<Point, ElementarySurface*> > id  ;
	switch(pos)
	{
		case TOP:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).y-maxy) < tol )
					{
						id[dist(elements[i]->getBoundingPoint(j), nearest)] = std::make_pair(elements[i]->getBoundingPoint(j), elements[i]) ;
					}
				}
			}
			std::vector<Point> target ;
			target.push_back(id.begin()->second.first) ;
			apply2DBC(id.begin()->second.second, target, condition, data, a) ;
			break ;
		}
		case LEFT:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).x-minx) < tol)
					{
						id[dist(elements[i]->getBoundingPoint(j), nearest)] = std::make_pair(elements[i]->getBoundingPoint(j), elements[i]) ;
					}
				}
			}
			std::vector<Point> target ;
			target.push_back(id.begin()->second.first) ;
			apply2DBC(id.begin()->second.second, target, condition, data, a) ;
			break ;
		}
		case BOTTOM:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).y-miny) < tol)
					{
						id[dist(elements[i]->getBoundingPoint(j), nearest)] = std::make_pair(elements[i]->getBoundingPoint(j), elements[i]) ;
					}
				}
			}
			std::vector<Point> target ;
			target.push_back(id.begin()->second.first) ;
			apply2DBC(id.begin()->second.second, target, condition, data, a) ;
			break ;
		}
		case RIGHT:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).x-maxx) < tol)
					{
						id[dist(elements[i]->getBoundingPoint(j), nearest)] = std::make_pair(elements[i]->getBoundingPoint(j), elements[i]) ;
					}
				}
			}
			std::vector<Point> target ;
			target.push_back(id.begin()->second.first) ;
			apply2DBC(id.begin()->second.second, target, condition, data, a) ;
			break ;
		}
		case TOP_LEFT:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).x-minx) < tol && std::abs(elements[i]->getBoundingPoint(j).y-maxy) < tol)
					{
						id[dist(elements[i]->getBoundingPoint(j), nearest)] = std::make_pair(elements[i]->getBoundingPoint(j), elements[i]) ;
					}
				}
			}
			std::vector<Point> target ;
			target.push_back(id.begin()->second.first) ;
			apply2DBC(id.begin()->second.second, target, condition, data, a) ;
			break ;
		}
		case TOP_RIGHT:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).x-maxx) < tol && std::abs(elements[i]->getBoundingPoint(j).y-maxy) < tol)
					{
						id[dist(elements[i]->getBoundingPoint(j), nearest)] = std::make_pair(elements[i]->getBoundingPoint(j), elements[i]) ;
					}
				}
			}
			std::vector<Point> target ;
			target.push_back(id.begin()->second.first) ;
			apply2DBC(id.begin()->second.second, target, condition, data, a) ;
			break ;
		}
		case BOTTOM_LEFT:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).x-minx) < tol && std::abs(elements[i]->getBoundingPoint(j).y-miny) < tol)
					{
						id[dist(elements[i]->getBoundingPoint(j), nearest)] = std::make_pair(elements[i]->getBoundingPoint(j), elements[i]) ;
					}
				}
			}
			std::vector<Point> target ;
			target.push_back(id.begin()->second.first) ;
			apply2DBC(id.begin()->second.second, target, condition, data, a) ;
			break ;
		}
		case BOTTOM_RIGHT:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).x-maxx) < tol && std::abs(elements[i]->getBoundingPoint(j).y-miny) < tol)
					{
						id[dist(elements[i]->getBoundingPoint(j), nearest)] = std::make_pair(elements[i]->getBoundingPoint(j), elements[i]) ;
					}
				}
			}
			std::vector<Point> target ;
			target.push_back(id.begin()->second.first) ;
			apply2DBC(id.begin()->second.second, target, condition, data, a) ;
			break ;
		}
		default:
		{
			break;
		}
	}
}

void BoundingBoxNearestNodeDefinedBoundaryCondition::apply(Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t)  const 
{
	std::vector<ElementaryVolume *> & elements = a->getElements3d() ;
	if(elements.empty())
	{
		std::cout << "no elements in assembly" << std::endl ;
		return ;
	}
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
	std::map<double, std::pair<Point, ElementaryVolume *> > id  ;
	switch(pos)
	{
		case TOP:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).y-maxy) < tol)
					{
						id[dist(elements[i]->getBoundingPoint(j), nearest)] = std::make_pair(elements[i]->getBoundingPoint(j), elements[i]) ;
					}
				}
			}
			std::vector<Point> target ;
			target.push_back(id.begin()->second.first) ;
			apply3DBC(id.begin()->second.second, target, condition, data, a) ;
			break ;
		}
		case LEFT:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).x-minx) < tol)
					{
						id[dist(elements[i]->getBoundingPoint(j), nearest)] = std::make_pair(elements[i]->getBoundingPoint(j), elements[i]) ;
					}
				}
			}
			std::vector<Point> target ;
			target.push_back(id.begin()->second.first) ;
			apply3DBC(id.begin()->second.second, target, condition, data, a) ;
			break ;
		}
		case BOTTOM:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).y-miny) < tol)
					{
						id[dist(elements[i]->getBoundingPoint(j), nearest)] = std::make_pair(elements[i]->getBoundingPoint(j), elements[i]) ;
					}
				}
			}
			std::vector<Point> target ;
			target.push_back(id.begin()->second.first) ;
			apply3DBC(id.begin()->second.second, target, condition, data, a) ;
			break ;
		}
		case RIGHT:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).x-maxx) < tol)
					{
						id[dist(elements[i]->getBoundingPoint(j), nearest)] = std::make_pair(elements[i]->getBoundingPoint(j), elements[i]) ;
					}
				}
			}
			std::vector<Point> target ;
			target.push_back(id.begin()->second.first) ;
			apply3DBC(id.begin()->second.second, target, condition, data, a) ;
			break ;
		}
		case FRONT:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).z-maxz) < tol)
					{
						id[dist(elements[i]->getBoundingPoint(j), nearest)] = std::make_pair(elements[i]->getBoundingPoint(j), elements[i]) ;
					}
				}
			}
			std::vector<Point> target ;
			target.push_back(id.begin()->second.first) ;
			apply3DBC(id.begin()->second.second, target, condition, data, a) ;
			break ;
		}
		case BACK:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).z-minz) < tol)
					{
						id[dist(elements[i]->getBoundingPoint(j), nearest)] = std::make_pair(elements[i]->getBoundingPoint(j), elements[i]) ;
					}
				}
			}
			std::vector<Point> target ;
			target.push_back(id.begin()->second.first) ;
			apply3DBC(id.begin()->second.second, target, condition, data, a) ;
			break ;
		}
		case TOP_LEFT:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).x-minx) < tol && std::abs(elements[i]->getBoundingPoint(j).y-maxy) < tol)
					{
						id[dist(elements[i]->getBoundingPoint(j), nearest)] = std::make_pair(elements[i]->getBoundingPoint(j), elements[i]) ;
					}
				}
			}
			std::vector<Point> target ;
			target.push_back(id.begin()->second.first) ;
			apply3DBC(id.begin()->second.second, target, condition, data, a) ;
			break ;
		}
		case TOP_RIGHT:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).x-maxx) < tol && std::abs(elements[i]->getBoundingPoint(j).y-maxy) < tol)
					{
						id[dist(elements[i]->getBoundingPoint(j), nearest)] = std::make_pair(elements[i]->getBoundingPoint(j), elements[i]) ;
					}
				}
			}
			std::vector<Point> target ;
			target.push_back(id.begin()->second.first) ;
			apply3DBC(id.begin()->second.second, target, condition, data, a) ;
			break ;
		}
		case BOTTOM_LEFT:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).x-minx) < tol && std::abs(elements[i]->getBoundingPoint(j).y-miny) < tol)
					{
						id[dist(elements[i]->getBoundingPoint(j), nearest)] = std::make_pair(elements[i]->getBoundingPoint(j), elements[i]) ;
					}
				}
			}
			std::vector<Point> target ;
			target.push_back(id.begin()->second.first) ;
			apply3DBC(id.begin()->second.second, target, condition, data, a) ;
			break ;
		}
		case BOTTOM_RIGHT:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).x-maxx) < tol && std::abs(elements[i]->getBoundingPoint(j).y-miny) < tol)
					{
						id[dist(elements[i]->getBoundingPoint(j), nearest)] = std::make_pair(elements[i]->getBoundingPoint(j), elements[i]) ;
					}
				}
			}
			std::vector<Point> target ;
			target.push_back(id.begin()->second.first) ;
			apply3DBC(id.begin()->second.second, target, condition, data, a) ;
			break ;
		}
		case FRONT_LEFT:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).x-minx) < tol && std::abs(elements[i]->getBoundingPoint(j).z-maxz) < tol)
					{
						id[dist(elements[i]->getBoundingPoint(j), nearest)] = std::make_pair(elements[i]->getBoundingPoint(j), elements[i]) ;
					}
				}
			}
			std::vector<Point> target ;
			target.push_back(id.begin()->second.first) ;
			apply3DBC(id.begin()->second.second, target, condition, data, a) ;
			break ;
		}
		case FRONT_RIGHT:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).x-maxx) < tol && std::abs(elements[i]->getBoundingPoint(j).z-maxz) < tol)
					{
						id[dist(elements[i]->getBoundingPoint(j), nearest)] = std::make_pair(elements[i]->getBoundingPoint(j), elements[i]) ;
					}
				}
			}
			std::vector<Point> target ;
			target.push_back(id.begin()->second.first) ;
			apply3DBC(id.begin()->second.second, target, condition, data, a) ;
			break ;
		}
		case BACK_LEFT:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).x-minx) < tol && std::abs(elements[i]->getBoundingPoint(j).z-minz) < tol)
					{
						id[dist(elements[i]->getBoundingPoint(j), nearest)] = std::make_pair(elements[i]->getBoundingPoint(j), elements[i]) ;
					}
				}
			}
			std::vector<Point> target ;
			target.push_back(id.begin()->second.first) ;
			apply3DBC(id.begin()->second.second, target, condition, data, a) ;
			break ;
		}
		case BACK_RIGHT:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).x-maxx) < tol && std::abs(elements[i]->getBoundingPoint(j).z-minz) < tol)
					{
						id[dist(elements[i]->getBoundingPoint(j), nearest)] = std::make_pair(elements[i]->getBoundingPoint(j), elements[i]) ;
					}
				}
			}
			std::vector<Point> target ;
			target.push_back(id.begin()->second.first) ;
			apply3DBC(id.begin()->second.second, target, condition, data, a) ;
			break ;
		}
		case FRONT_TOP:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).y-maxy) < tol && std::abs(elements[i]->getBoundingPoint(j).z-maxz) < tol)
					{
						id[dist(elements[i]->getBoundingPoint(j), nearest)] = std::make_pair(elements[i]->getBoundingPoint(j), elements[i]) ;
					}
				}
			}
			std::vector<Point> target ;
			target.push_back(id.begin()->second.first) ;
			apply3DBC(id.begin()->second.second, target, condition, data, a) ;
			break ;
		}
		case FRONT_BOTTOM:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).y-miny) < tol && std::abs(elements[i]->getBoundingPoint(j).z-maxz) < tol)
					{
						id[dist(elements[i]->getBoundingPoint(j), nearest)] = std::make_pair(elements[i]->getBoundingPoint(j), elements[i]) ;
					}
				}
			}
			std::vector<Point> target ;
			target.push_back(id.begin()->second.first) ;
			apply3DBC(id.begin()->second.second, target, condition, data, a) ;
			break ;
		}
		case TOP_LEFT_FRONT:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).x-minx) < tol 
						&& std::abs(elements[i]->getBoundingPoint(j).y-maxy) < tol
						&& std::abs(elements[i]->getBoundingPoint(j).z-maxz) < tol)
					{
						id[dist(elements[i]->getBoundingPoint(j), nearest)] = std::make_pair(elements[i]->getBoundingPoint(j), elements[i]) ;
					}
				}
			}
			std::vector<Point> target ;
			target.push_back(id.begin()->second.first) ;
			apply3DBC(id.begin()->second.second, target, condition, data, a) ;
			break ;
		}
		case TOP_LEFT_BACK:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).x-minx) < tol 
						&& std::abs(elements[i]->getBoundingPoint(j).y-maxy) < tol
						&& std::abs(elements[i]->getBoundingPoint(j).z-minz) < tol)
					{
						id[dist(elements[i]->getBoundingPoint(j), nearest)] = std::make_pair(elements[i]->getBoundingPoint(j), elements[i]) ;
					}
				}
			}
			std::vector<Point> target ;
			target.push_back(id.begin()->second.first) ;
			apply3DBC(id.begin()->second.second, target, condition, data, a) ;
			break ;
		}
		case BOTTOM_LEFT_FRONT:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).x-minx) < tol 
						&& std::abs(elements[i]->getBoundingPoint(j).y-miny) < tol
						&& std::abs(elements[i]->getBoundingPoint(j).z-maxz) < tol)
					{
						id[dist(elements[i]->getBoundingPoint(j), nearest)] = std::make_pair(elements[i]->getBoundingPoint(j), elements[i]) ;
					}
				}
			}
			std::vector<Point> target ;
			target.push_back(id.begin()->second.first) ;
			apply3DBC(id.begin()->second.second, target, condition, data, a) ;
			break ;
		}
		case BOTTOM_LEFT_BACK:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).x-minx) < tol 
						&& std::abs(elements[i]->getBoundingPoint(j).y-miny) < tol
						&& std::abs(elements[i]->getBoundingPoint(j).z-minz) < tol)
					{
						id[dist(elements[i]->getBoundingPoint(j), nearest)] = std::make_pair(elements[i]->getBoundingPoint(j), elements[i]) ;
					}
				}
			}
			std::vector<Point> target ;
			target.push_back(id.begin()->second.first) ;
			apply3DBC(id.begin()->second.second, target, condition, data, a) ;
			break ;
		}
		case TOP_RIGHT_FRONT:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).x-maxx) < tol 
						&& std::abs(elements[i]->getBoundingPoint(j).y-maxy) < tol
						&& std::abs(elements[i]->getBoundingPoint(j).z-minz) < tol)
					{
						id[dist(elements[i]->getBoundingPoint(j), nearest)] = std::make_pair(elements[i]->getBoundingPoint(j), elements[i]) ;
					}
				}
			}
			std::vector<Point> target ;
			target.push_back(id.begin()->second.first) ;
			apply3DBC(id.begin()->second.second, target, condition, data, a) ;
			break ;
		}
		case TOP_RIGHT_BACK:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).x-maxx) < tol 
						&& std::abs(elements[i]->getBoundingPoint(j).y-maxy) < tol
						&& std::abs(elements[i]->getBoundingPoint(j).z-maxz) < tol)
					{
						id[dist(elements[i]->getBoundingPoint(j), nearest)] = std::make_pair(elements[i]->getBoundingPoint(j), elements[i]) ;
					}
				}
			}
			std::vector<Point> target ;
			target.push_back(id.begin()->second.first) ;
			apply3DBC(id.begin()->second.second, target, condition, data, a) ;
			break ;
		}
		case BOTTOM_RIGHT_FRONT:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).x-maxx) < tol 
						&& std::abs(elements[i]->getBoundingPoint(j).y-miny) < tol
						&& std::abs(elements[i]->getBoundingPoint(j).z-minz) < tol)
					{
						id[dist(elements[i]->getBoundingPoint(j), nearest)] = std::make_pair(elements[i]->getBoundingPoint(j), elements[i]) ;
					}
				}
			}
			std::vector<Point> target ;
			target.push_back(id.begin()->second.first) ;
			apply3DBC(id.begin()->second.second, target, condition, data, a) ;
			break ;
		}
		case BOTTOM_RIGHT_BACK:
		{
			for(size_t i = 0 ; i < elements.size() ; ++i)
			{
				for(size_t j = 0 ;  j< elements[i]->getBoundingPoints().size() ; ++j)
				{
					if(std::abs(elements[i]->getBoundingPoint(j).x-maxx) < tol 
						&& std::abs(elements[i]->getBoundingPoint(j).y-miny) < tol
						&& std::abs(elements[i]->getBoundingPoint(j).z-maxz) < tol)
					{
						id[dist(elements[i]->getBoundingPoint(j), nearest)] = std::make_pair(elements[i]->getBoundingPoint(j), elements[i]) ;
					}
				}
			}
			std::vector<Point> target ;
			target.push_back(id.begin()->second.first) ;
			apply3DBC(id.begin()->second.second, target, condition, data, a) ;
			break ;
		}
		default:
		{
			break;
		}
	}
}


GeometryDefinedBoundaryCondition::GeometryDefinedBoundaryCondition(LagrangeMultiplierType t, Geometry * source, double d) : BoundaryCondition(t, d), domain(source) { };

void GeometryDefinedBoundaryCondition::apply(Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t) const
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
void GeometryDefinedBoundaryCondition::apply(Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t)  const
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
void BoundingBoxAndRestrictionDefinedBoundaryCondition::apply(Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t) const
{
	std::vector<ElementarySurface *> & elements = a->getElements2d() ;
	if(elements.empty())
	{
		std::cout << "no elements in assembly" << std::endl ;
		return ;
	}
	double minx = elements.front()->getBoundingPoint(0).x ;
	double miny = elements.front()->getBoundingPoint(0).y ;
	double maxx = elements.front()->getBoundingPoint(0).x ;
	double maxy = elements.front()->getBoundingPoint(0).y ; 
	for(size_t i = 0 ; i < elements.size() ; ++i)
	{
		if(elements[i]->getBehaviour()->type != VOID_BEHAVIOUR)
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
					if(std::abs(elements[i]->getBoundingPoint(j).y-maxy) < tol 
						&& elements[i]->getBoundingPoint(j).x >= xmin 
						&& elements[i]->getBoundingPoint(j).x <= xmax
						&& elements[i]->getBoundingPoint(j).y >= ymin 
						&& elements[i]->getBoundingPoint(j).y <= ymax
						)
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
					if(std::abs(elements[i]->getBoundingPoint(j).x-minx) < tol
						&& elements[i]->getBoundingPoint(j).x >= xmin 
						&& elements[i]->getBoundingPoint(j).x <= xmax
						&& elements[i]->getBoundingPoint(j).y >= ymin 
						&& elements[i]->getBoundingPoint(j).y <= ymax
						)
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
					if(std::abs(elements[i]->getBoundingPoint(j).y-miny) < tol
						&& elements[i]->getBoundingPoint(j).x >= xmin 
						&& elements[i]->getBoundingPoint(j).x <= xmax
						&& elements[i]->getBoundingPoint(j).y >= ymin 
						&& elements[i]->getBoundingPoint(j).y <= ymax
						)
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
					if(std::abs(elements[i]->getBoundingPoint(j).x-maxx) < tol
						&& elements[i]->getBoundingPoint(j).x >= xmin 
						&& elements[i]->getBoundingPoint(j).x <= xmax
						&& elements[i]->getBoundingPoint(j).y >= ymin 
						&& elements[i]->getBoundingPoint(j).y <= ymax
						)
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
					if(std::abs(elements[i]->getBoundingPoint(j).x-minx) < tol && std::abs(elements[i]->getBoundingPoint(j).y-maxy) < tol
						&& elements[i]->getBoundingPoint(j).x >= xmin 
						&& elements[i]->getBoundingPoint(j).x <= xmax
						&& elements[i]->getBoundingPoint(j).y >= ymin 
						&& elements[i]->getBoundingPoint(j).y <= ymax
						)
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
					if(std::abs(elements[i]->getBoundingPoint(j).x-maxx) < tol && std::abs(elements[i]->getBoundingPoint(j).y-maxy) < tol
						&& elements[i]->getBoundingPoint(j).x >= xmin 
						&& elements[i]->getBoundingPoint(j).x <= xmax
						&& elements[i]->getBoundingPoint(j).y >= ymin 
						&& elements[i]->getBoundingPoint(j).y <= ymax
						)
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
					if(std::abs(elements[i]->getBoundingPoint(j).x-minx) < tol && std::abs(elements[i]->getBoundingPoint(j).y-miny) < tol
						&& elements[i]->getBoundingPoint(j).x >= xmin 
						&& elements[i]->getBoundingPoint(j).x <= xmax
						&& elements[i]->getBoundingPoint(j).y >= ymin 
						&& elements[i]->getBoundingPoint(j).y <= ymax
						)
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
					if(std::abs(elements[i]->getBoundingPoint(j).x-maxx) < tol && std::abs(elements[i]->getBoundingPoint(j).y-miny) < tol
						&& elements[i]->getBoundingPoint(j).x >= xmin 
						&& elements[i]->getBoundingPoint(j).x <= xmax
						&& elements[i]->getBoundingPoint(j).y >= ymin 
						&& elements[i]->getBoundingPoint(j).y <= ymax
						)
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
void BoundingBoxDefinedBoundaryCondition::apply(Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t) const
{
	std::vector<ElementarySurface *> & elements = a->getElements2d() ;
	if(elements.empty())
	{
		std::cout << "no elements in assembly" << std::endl ;
		return ;
	}
	double minx = elements.front()->getBoundingPoint(0).x ;
	double miny = elements.front()->getBoundingPoint(0).y ;
	double maxx = elements.front()->getBoundingPoint(0).x ;
	double maxy = elements.front()->getBoundingPoint(0).y ; 
	for(size_t i = 0 ; i < elements.size() ; ++i)
	{
		if(elements[i]->getBehaviour()->type != VOID_BEHAVIOUR)
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

void BoundingBoxAndRestrictionDefinedBoundaryCondition::apply(Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t)  const 
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
					if(std::abs(elements[i]->getBoundingPoint(j).y-maxy) < tol
						&& elements[i]->getBoundingPoint(j).x >= xmin 
						&& elements[i]->getBoundingPoint(j).x <= xmax
						&& elements[i]->getBoundingPoint(j).y >= ymin 
						&& elements[i]->getBoundingPoint(j).y <= ymax
						&& elements[i]->getBoundingPoint(j).z >= zmin
						&& elements[i]->getBoundingPoint(j).z <= zmax
						)
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
					if(std::abs(elements[i]->getBoundingPoint(j).x-minx) < tol
						&& elements[i]->getBoundingPoint(j).x >= xmin 
						&& elements[i]->getBoundingPoint(j).x <= xmax
						&& elements[i]->getBoundingPoint(j).y >= ymin 
						&& elements[i]->getBoundingPoint(j).y <= ymax
						&& elements[i]->getBoundingPoint(j).z >= zmin
						&& elements[i]->getBoundingPoint(j).z <= zmax
						)
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
					if(std::abs(elements[i]->getBoundingPoint(j).y-miny) < tol
						&& elements[i]->getBoundingPoint(j).x >= xmin 
						&& elements[i]->getBoundingPoint(j).x <= xmax
						&& elements[i]->getBoundingPoint(j).y >= ymin 
						&& elements[i]->getBoundingPoint(j).y <= ymax
						&& elements[i]->getBoundingPoint(j).z >= zmin
						&& elements[i]->getBoundingPoint(j).z <= zmax
						)
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
					if(std::abs(elements[i]->getBoundingPoint(j).x-maxx) < tol
						&& elements[i]->getBoundingPoint(j).x >= xmin 
						&& elements[i]->getBoundingPoint(j).x <= xmax
						&& elements[i]->getBoundingPoint(j).y >= ymin 
						&& elements[i]->getBoundingPoint(j).y <= ymax
						&& elements[i]->getBoundingPoint(j).z >= zmin
						&& elements[i]->getBoundingPoint(j).z <= zmax
						)
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
					if(std::abs(elements[i]->getBoundingPoint(j).z-maxz) < tol
						&& elements[i]->getBoundingPoint(j).x >= xmin 
						&& elements[i]->getBoundingPoint(j).x <= xmax
						&& elements[i]->getBoundingPoint(j).y >= ymin 
						&& elements[i]->getBoundingPoint(j).y <= ymax
						&& elements[i]->getBoundingPoint(j).z >= zmin
						&& elements[i]->getBoundingPoint(j).z <= zmax
						)
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
					if(std::abs(elements[i]->getBoundingPoint(j).z-minz) < tol
						&& elements[i]->getBoundingPoint(j).x >= xmin 
						&& elements[i]->getBoundingPoint(j).x <= xmax
						&& elements[i]->getBoundingPoint(j).y >= ymin 
						&& elements[i]->getBoundingPoint(j).y <= ymax
						&& elements[i]->getBoundingPoint(j).z >= zmin
						&& elements[i]->getBoundingPoint(j).z <= zmax
						)
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
					if(std::abs(elements[i]->getBoundingPoint(j).x-minx) < tol && std::abs(elements[i]->getBoundingPoint(j).y-maxy) < tol
						&& elements[i]->getBoundingPoint(j).x >= xmin 
						&& elements[i]->getBoundingPoint(j).x <= xmax
						&& elements[i]->getBoundingPoint(j).y >= ymin 
						&& elements[i]->getBoundingPoint(j).y <= ymax
						&& elements[i]->getBoundingPoint(j).z >= zmin
						&& elements[i]->getBoundingPoint(j).z <= zmax
						)
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
					if(std::abs(elements[i]->getBoundingPoint(j).x-maxx) < tol && std::abs(elements[i]->getBoundingPoint(j).y-maxy) < tol
						&& elements[i]->getBoundingPoint(j).x >= xmin 
						&& elements[i]->getBoundingPoint(j).x <= xmax
						&& elements[i]->getBoundingPoint(j).y >= ymin 
						&& elements[i]->getBoundingPoint(j).y <= ymax
						&& elements[i]->getBoundingPoint(j).z >= zmin
						&& elements[i]->getBoundingPoint(j).z <= zmax
						)
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
					if(std::abs(elements[i]->getBoundingPoint(j).x-minx) < tol && std::abs(elements[i]->getBoundingPoint(j).y-miny) < tol
						&& elements[i]->getBoundingPoint(j).x >= xmin 
						&& elements[i]->getBoundingPoint(j).x <= xmax
						&& elements[i]->getBoundingPoint(j).y >= ymin 
						&& elements[i]->getBoundingPoint(j).y <= ymax
						&& elements[i]->getBoundingPoint(j).z >= zmin
						&& elements[i]->getBoundingPoint(j).z <= zmax
						)
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
					if(std::abs(elements[i]->getBoundingPoint(j).x-maxx) < tol && std::abs(elements[i]->getBoundingPoint(j).y-miny) < tol
						&& elements[i]->getBoundingPoint(j).x >= xmin 
						&& elements[i]->getBoundingPoint(j).x <= xmax
						&& elements[i]->getBoundingPoint(j).y >= ymin 
						&& elements[i]->getBoundingPoint(j).y <= ymax
						&& elements[i]->getBoundingPoint(j).z >= zmin
						&& elements[i]->getBoundingPoint(j).z <= zmax
						)
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
					if(std::abs(elements[i]->getBoundingPoint(j).x-minx) < tol && std::abs(elements[i]->getBoundingPoint(j).z-maxz) < tol
						&& elements[i]->getBoundingPoint(j).x >= xmin 
						&& elements[i]->getBoundingPoint(j).x <= xmax
						&& elements[i]->getBoundingPoint(j).y >= ymin 
						&& elements[i]->getBoundingPoint(j).y <= ymax
						&& elements[i]->getBoundingPoint(j).z >= zmin
						&& elements[i]->getBoundingPoint(j).z <= zmax
						)
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
					if(std::abs(elements[i]->getBoundingPoint(j).x-maxx) < tol && std::abs(elements[i]->getBoundingPoint(j).z-maxz) < tol
						&& elements[i]->getBoundingPoint(j).x >= xmin 
						&& elements[i]->getBoundingPoint(j).x <= xmax
						&& elements[i]->getBoundingPoint(j).y >= ymin 
						&& elements[i]->getBoundingPoint(j).y <= ymax
						&& elements[i]->getBoundingPoint(j).z >= zmin
						&& elements[i]->getBoundingPoint(j).z <= zmax
						)
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
					if(std::abs(elements[i]->getBoundingPoint(j).x-minx) < tol && std::abs(elements[i]->getBoundingPoint(j).z-minz) < tol
						&& elements[i]->getBoundingPoint(j).x >= xmin 
						&& elements[i]->getBoundingPoint(j).x <= xmax
						&& elements[i]->getBoundingPoint(j).y >= ymin 
						&& elements[i]->getBoundingPoint(j).y <= ymax
						&& elements[i]->getBoundingPoint(j).z >= zmin
						&& elements[i]->getBoundingPoint(j).z <= zmax
						)
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
					if(std::abs(elements[i]->getBoundingPoint(j).x-maxx) < tol && std::abs(elements[i]->getBoundingPoint(j).z-minz) < tol
						&& elements[i]->getBoundingPoint(j).x >= xmin 
						&& elements[i]->getBoundingPoint(j).x <= xmax
						&& elements[i]->getBoundingPoint(j).y >= ymin 
						&& elements[i]->getBoundingPoint(j).y <= ymax
						&& elements[i]->getBoundingPoint(j).z >= zmin
						&& elements[i]->getBoundingPoint(j).z <= zmax
						)
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
					if(std::abs(elements[i]->getBoundingPoint(j).y-maxy) < tol && std::abs(elements[i]->getBoundingPoint(j).z-maxz) < tol
						&& elements[i]->getBoundingPoint(j).x >= xmin 
						&& elements[i]->getBoundingPoint(j).x <= xmax
						&& elements[i]->getBoundingPoint(j).y >= ymin 
						&& elements[i]->getBoundingPoint(j).y <= ymax
						&& elements[i]->getBoundingPoint(j).z >= zmin
						&& elements[i]->getBoundingPoint(j).z <= zmax
						)
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
					if(std::abs(elements[i]->getBoundingPoint(j).y-miny) < tol && std::abs(elements[i]->getBoundingPoint(j).z-maxz) < tol
						&& elements[i]->getBoundingPoint(j).x >= xmin 
						&& elements[i]->getBoundingPoint(j).x <= xmax
						&& elements[i]->getBoundingPoint(j).y >= ymin 
						&& elements[i]->getBoundingPoint(j).y <= ymax
						&& elements[i]->getBoundingPoint(j).z >= zmin
						&& elements[i]->getBoundingPoint(j).z <= zmax
						)
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
						&& std::abs(elements[i]->getBoundingPoint(j).z-maxz) < tol
						&& elements[i]->getBoundingPoint(j).x >= xmin 
						&& elements[i]->getBoundingPoint(j).x <= xmax
						&& elements[i]->getBoundingPoint(j).y >= ymin 
						&& elements[i]->getBoundingPoint(j).y <= ymax
						&& elements[i]->getBoundingPoint(j).z >= zmin
						&& elements[i]->getBoundingPoint(j).z <= zmax
						)
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
						&& std::abs(elements[i]->getBoundingPoint(j).z-minz) < tol
						&& elements[i]->getBoundingPoint(j).x >= xmin 
						&& elements[i]->getBoundingPoint(j).x <= xmax
						&& elements[i]->getBoundingPoint(j).y >= ymin 
						&& elements[i]->getBoundingPoint(j).y <= ymax
						&& elements[i]->getBoundingPoint(j).z >= zmin
						&& elements[i]->getBoundingPoint(j).z <= zmax
						)
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
						&& std::abs(elements[i]->getBoundingPoint(j).z-maxz) < tol
						&& elements[i]->getBoundingPoint(j).x >= xmin 
						&& elements[i]->getBoundingPoint(j).x <= xmax
						&& elements[i]->getBoundingPoint(j).y >= ymin 
						&& elements[i]->getBoundingPoint(j).y <= ymax
						&& elements[i]->getBoundingPoint(j).z >= zmin
						&& elements[i]->getBoundingPoint(j).z <= zmax
						)
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
						&& std::abs(elements[i]->getBoundingPoint(j).z-minz) < tol
						&& elements[i]->getBoundingPoint(j).x >= xmin 
						&& elements[i]->getBoundingPoint(j).x <= xmax
						&& elements[i]->getBoundingPoint(j).y >= ymin 
						&& elements[i]->getBoundingPoint(j).y <= ymax
						&& elements[i]->getBoundingPoint(j).z >= zmin
						&& elements[i]->getBoundingPoint(j).z <= zmax
						)
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
						&& std::abs(elements[i]->getBoundingPoint(j).z-minz) < tol
						&& elements[i]->getBoundingPoint(j).x >= xmin 
						&& elements[i]->getBoundingPoint(j).x <= xmax
						&& elements[i]->getBoundingPoint(j).y >= ymin 
						&& elements[i]->getBoundingPoint(j).y <= ymax
						&& elements[i]->getBoundingPoint(j).z >= zmin
						&& elements[i]->getBoundingPoint(j).z <= zmax
						)
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
						&& std::abs(elements[i]->getBoundingPoint(j).z-maxz) < tol
						&& elements[i]->getBoundingPoint(j).x >= xmin 
						&& elements[i]->getBoundingPoint(j).x <= xmax
						&& elements[i]->getBoundingPoint(j).y >= ymin 
						&& elements[i]->getBoundingPoint(j).y <= ymax
						&& elements[i]->getBoundingPoint(j).z >= zmin
						&& elements[i]->getBoundingPoint(j).z <= zmax
						)
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
						&& std::abs(elements[i]->getBoundingPoint(j).z-minz) < tol
						&& elements[i]->getBoundingPoint(j).x >= xmin 
						&& elements[i]->getBoundingPoint(j).x <= xmax
						&& elements[i]->getBoundingPoint(j).y >= ymin 
						&& elements[i]->getBoundingPoint(j).y <= ymax
						&& elements[i]->getBoundingPoint(j).z >= zmin
						&& elements[i]->getBoundingPoint(j).z <= zmax
						)
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
						&& std::abs(elements[i]->getBoundingPoint(j).z-maxz) < tol
						&& elements[i]->getBoundingPoint(j).x >= xmin 
						&& elements[i]->getBoundingPoint(j).x <= xmax
						&& elements[i]->getBoundingPoint(j).y >= ymin 
						&& elements[i]->getBoundingPoint(j).y <= ymax
						&& elements[i]->getBoundingPoint(j).z >= zmin
						&& elements[i]->getBoundingPoint(j).z <= zmax
						)
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

void BoundingBoxDefinedBoundaryCondition::apply(Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t)  const 
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

BoundaryCondition::BoundaryCondition(LagrangeMultiplierType t, const double & d) : condition(t), data(d) { } ;



ProjectionDefinedBoundaryCondition::ProjectionDefinedBoundaryCondition(LagrangeMultiplierType t, const Point & dir, double d) : BoundaryCondition(t,d), direction(dir) { }

void ProjectionDefinedBoundaryCondition::apply(Assembly * a, Mesh<DelaunayTriangle, DelaunayTreeItem> * t) const
{
	std::vector<DelaunayTriangle *> tris = t->getElements() ;
	for(size_t i = 0 ; i < tris.size() ; i++)
	{
		DelaunayTreeItem * VoidItem ;
		bool border = false ;
		for(size_t j = 0 ; j < tris[i]->neighbour.size() ; j++)
		{
			bool voidNeighbour = (tris[i]->getNeighbour(j)->isTriangle 
			                    && dynamic_cast<DelaunayTriangle *>(tris[i]->getNeighbour(j))->getBehaviour()->type == VOID_BEHAVIOUR) ;
			border = border || tris[i]->getNeighbour(j)->isPlane
			                || voidNeighbour ;
			if(voidNeighbour)
				VoidItem = tris[i]->getNeighbour(j) ;
			if(tris[i]->getNeighbour(j)->isPlane)
				VoidItem = tris[i]->getNeighbour(j) ;
		}
		
		if(tris[i]->getBehaviour()->type == VOID_BEHAVIOUR)
			border = false ;
		
		if(border)
		{
			std::pair<Point *, Point*> commonSurface = tris[i]->commonEdge(VoidItem) ;
			
			Segment ray((tris[i]->getCenter()), (tris[i]->getCenter())-direction*(tris[i]->getRadius())) ;
			bool isOnTheRightSide = ray.intersects(Segment(*commonSurface.first, *commonSurface.second)) ;
			
			if(isOnTheRightSide)
			{
				std::vector<Point> id ;
				for(size_t j = 0 ; j < tris[i]->getBoundingPoints().size() ; j++)
				{
					
					Line side(tris[i]->getBoundingPoint(j), tris[i]->getBoundingPoint(j)-tris[i]->getBoundingPoint((j+1)%tris[i]->getBoundingPoints().size())) ;
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

void ProjectionDefinedBoundaryCondition::apply(Assembly * a, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * t)  const
{
	std::vector<DelaunayTetrahedron *> tris = t->getElements() ;
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
			Plane surf(*points[0], *points[1], *points[2]) ;
			for(size_t j = 3 ; j < tris[i]->getBoundingPoints().size() ; j++)
			{
				if(isCoplanar(*points[0], *points[1], *points[2], tris[i]->getBoundingPoint(j)))
					points.push_back(&tris[i]->getBoundingPoint(j)) ;
			}
			
			if(surf.intersects(ray) && !points.empty())
			{
				for(size_t j = 0 ; j < points.size() ; j++)
					id.push_back(*points[j]) ;
			}
		}
		apply3DBC(tris[i], id, condition, data, a) ;
	}
}

