//
// C++ Implementation: fracturecriterion
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "fracturecriterion.h"
#include "../damagemodels/damagemodel.h"
#include "../../mesher/delaunay.h"
#include "../../mesher/delaunay_3d.h"
#include "../../solvers/assembly.h"
using namespace Mu ;

FractureCriterion::FractureCriterion(MirrorState mirroring, double delta_x, double delta_y, double delta_z) : neighbourhoodradius(.0005), neighbourhoodvolume(-1), physicalCharacteristicRadius(.008), scoreAtState(-1), metInTension(false), metInCompression(false), metAtStep(false), mirroring(mirroring), delta_x(delta_x), delta_y(delta_y), delta_z(delta_z), deltaScoreAtState(0), deltaEnergyAtState(0), energyDamageDifferential(0), criterionDamageDifferential(0), energyIndexed(false), noEnergyUpdate(true), mesh2d(NULL), mesh3d(NULL), stable(true)
{
}

double FractureCriterion::getDeltaEnergy(const ElementState & s, double delta_d)
{
	Assembly K ;
	double originalscore = 0 ;
	double originalenergy = 0 ;
	double score = 0 ;
	double energy = 0 ;
	double volume = 0 ;

	if(s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
	{
		Circle c(getMaterialCharacteristicRadius(), s.getParent()->getCenter()) ;

		std::vector<DelaunayTriangle *> elements ;
		std::vector<LinearForm *> behaviours ;
		TriElement father(LINEAR) ;

		for(size_t i = 0 ; i < cache.size() ; i++)
		{
			DelaunayTriangle * tric = static_cast<DelaunayTriangle *>((*mesh2d)[cache[i]]) ;
			elements.push_back(new DelaunayTriangle(NULL, NULL,   tric->first,  tric->second,   tric->third,  NULL) );
			elements.back()->setBehaviour(tric->getBehaviour()->getCopy()) ;
			elements.back()->refresh(&father);
			elements.back()->getState().initialize() ;
			K.add(elements.back());
			volume += elements.back()->area() ;
			for(size_t j = 0 ; j < tric->neighbour.size() ; j++)
			{
				if(tric->getNeighbour(j)->isTriangle )
				{
					DelaunayTriangle * tri = static_cast<DelaunayTriangle *>(tric->getNeighbour(j)) ;
					if(tri->getBehaviour()->type != VOID_BEHAVIOUR)
					{
						for(size_t k = 0 ; k <  tri->getBoundingPoints().size() ; k++)
						{
							if(!c.in(tri->getBoundingPoint(k)))
							{
								for(size_t l = 0 ; l <  tri->getBoundingPoints().size() ; l++)
								{
									int id = tri->getBoundingPoint(k).id ;
									double ex = tri->getState().getDisplacements()[l*2];
									double ey = tri->getState().getDisplacements()[l*2+1];
									K.setPoint(ex, ey ,id);
									break ;
								}
							}
						}
					}
				}
			}
		}
		K.cgsolve() ;
		for(size_t i = 0 ; i < elements.size() ; i++)
			elements[i]->step(0., &K.getDisplacements()) ;
		for(size_t i = 0 ; i < elements.size() ;i++)
		{
			if(elements[i]->getBehaviour()->getFractureCriterion())
			{
				double a = elements[i]->area() ;
				originalenergy += elements[i]->getState().elasticEnergy()*a ;
			}
		}
		
		std::valarray<Point *> nularray(0) ;
		
		for(size_t i = 0 ; i < elements.size() ;i++)
		{
			elements[i]->setBoundingPoints(nularray) ;
			delete elements[i] ;
		}
		
	}
	else
	{
		Sphere c(getMaterialCharacteristicRadius(), s.getParent()->getCenter()) ;
		
		std::vector<DelaunayTetrahedron *> elements ;
		std::vector<LinearForm *> behaviours ;
		TetrahedralElement father ;
		for(size_t i = 0 ; i < cache.size() ; i++)
		{
			DelaunayTetrahedron * tet = static_cast<DelaunayTetrahedron *>((*mesh3d)[ cache[i]]) ;
			elements.push_back(new DelaunayTetrahedron(NULL, NULL,   tet->first,  tet->second,   tet->third, tet->fourth,  NULL) );
			elements.back()->setBehaviour(tet->getBehaviour()->getCopy()) ;
			
			elements.back()->refresh(&father);
			elements.back()->getState().initialize() ;
			K.add(elements.back());
			volume+= elements.back()->volume() ;
			
			for(size_t j = 0 ; j < tet->neighbour.size() ; j++)
			{
				if(tet->getNeighbour(j)->isTetrahedron() )
				{
					DelaunayTetrahedron * tri = dynamic_cast<DelaunayTetrahedron *>(tet->getNeighbour(j)) ;
					if(tri->getBehaviour()->type != VOID_BEHAVIOUR)
					{
						for(size_t k = 0 ; k <  tri->getBoundingPoints().size() ; k++)
						{
							if(!c.in(tri->getBoundingPoint(k)))
							{
								for(size_t l = 0 ; l <  tri->getBoundingPoints().size() ; l++)
								{
									int id = tri->getBoundingPoint(k).id ;
									double ex = tri->getState().getDisplacements()[l*3];
									double ey = tri->getState().getDisplacements()[l*3+1];
									double ez = tri->getState().getDisplacements()[l*3+2];
									K.setPoint(ex, ey, ez ,id);
									break ;
								}
							}
						}
					}
				}
			}
		
			
		}
	
		K.cgsolve() ;
		for(size_t i = 0 ; i < elements.size() ; i++)
			elements[i]->step(0., &K.getDisplacements()) ;
		
		double energy = 0 ;
		for(size_t i = 0 ; i < elements.size() ;i++)
		{
			if(elements[i]->getBehaviour()->getFractureCriterion())
			{
				double v = elements[i]->volume() ;
				originalenergy += elements[i]->getState().elasticEnergy()*v ;
			}
		}
		for(size_t i = 0 ; i < elements.size() ;i++)
		{
			delete elements[i] ;
		}
	}
	
	
	
	if(s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
	{
		Circle c(getMaterialCharacteristicRadius(), s.getParent()->getCenter()) ;

		std::vector<DelaunayTriangle *> elements ;
		std::vector<LinearForm *> behaviours ;
		TriElement father(LINEAR) ;

		for(size_t i = 0 ; i < cache.size() ; i++)
		{
			DelaunayTriangle * tric = static_cast<DelaunayTriangle *>((*mesh2d)[ cache[i]]) ;
			elements.push_back(new DelaunayTriangle(NULL, NULL,   tric->first,  tric->second,   tric->third,  NULL) );
			elements.back()->setBehaviour(tric->getBehaviour()->getCopy()) ;

			if(tric == s.getParent())
				elements.back()->getBehaviour()->setTensor(elements.back()->getBehaviour()->getTensor(elements.back()->getCenter())*(1.-delta_d)) ;

			elements.back()->getBehaviour()->getFractureCriterion()->setEnergyIndexed(false) ;
			elements.back()->refresh(&father);
			elements.back()->getState().initialize() ;
			K.add(elements.back());
			volume += elements.back()->area() ;
			for(size_t j = 0 ; j < tric->neighbour.size() ; j++)
			{
				if(tric->getNeighbour(j)->isTriangle )
				{
					DelaunayTriangle * tri = dynamic_cast<DelaunayTriangle *>(tric->getNeighbour(j)) ;
					if(tri->getBehaviour()->type != VOID_BEHAVIOUR)
					{
						for(size_t k = 0 ; k <  tri->getBoundingPoints().size() ; k++)
						{
							if(!c.in(tri->getBoundingPoint(k)))
							{
								for(size_t l = 0 ; l <  tri->getBoundingPoints().size() ; l++)
								{
									int id = tri->getBoundingPoint(k).id ;
									double ex = tri->getState().getDisplacements()[l*2];
									double ey = tri->getState().getDisplacements()[l*2+1];
									K.setPoint(ex, ey ,id);
									break ;
								}
							}
						}
					}
				}
			}
		}
	
		K.cgsolve() ;
		for(size_t i = 0 ; i < elements.size() ; i++)
			elements[i]->step(0., &K.getDisplacements()) ;

		for(size_t i = 0 ; i < elements.size() ;i++)
		{
			if(elements[i]->getBehaviour()->getFractureCriterion())
			{
				double a = elements[i]->area() ;
				energy += elements[i]->getState().elasticEnergy()*a ;
			}
		}

		std::valarray<Point *> nularray(0) ;
		
		for(size_t i = 0 ; i < elements.size() ;i++)
		{
			elements[i]->setBoundingPoints(nularray) ;
			delete elements[i] ;
		}
	}
	else
	{
		Sphere c(getMaterialCharacteristicRadius(), s.getParent()->getCenter()) ;
		
		std::vector<DelaunayTetrahedron *> elements ;
		std::vector<LinearForm *> behaviours ;
		TetrahedralElement father ;
		for(size_t i = 0 ; i < cache.size() ; i++)
		{
			DelaunayTetrahedron * tet = static_cast<DelaunayTetrahedron * >((*mesh3d)[cache[i]]) ;
			elements.push_back(new DelaunayTetrahedron(NULL, NULL,   tet->first,  tet->second,   tet->third, tet->fourth,  NULL) );
			elements.back()->setBehaviour(tet->getBehaviour()->getCopy()) ;
			if(tet == s.getParent())
				elements.back()->getBehaviour()->setTensor(elements.back()->getBehaviour()->getTensor(elements.back()->getCenter())*(1.-delta_d)) ;
			elements.back()->getBehaviour()->getFractureCriterion()->setEnergyIndexed(false) ;
			elements.back()->refresh(&father);
			elements.back()->getState().initialize() ;
			K.add(elements.back());
			volume+= elements.back()->volume() ;
			
			for(size_t j = 0 ; j < tet->neighbour.size() ; j++)
			{
				if(tet->getNeighbour(j)->isTetrahedron() )
				{
					DelaunayTetrahedron * tri = dynamic_cast<DelaunayTetrahedron *>(tet->getNeighbour(j)) ;
					if(tri->getBehaviour()->type != VOID_BEHAVIOUR)
					{
						for(size_t k = 0 ; k <  tri->getBoundingPoints().size() ; k++)
						{
							if(!c.in(tri->getBoundingPoint(k)))
							{
								int id = tri->getBoundingPoint(k).id ;
								double ex = tri->getState().getDisplacements()[k*3];
								double ey = tri->getState().getDisplacements()[k*3+1];
								double ez = tri->getState().getDisplacements()[k*3+2];
								K.setPoint(ex, ey, ez ,id);
							}
						}
					}
				}
			}
		
			
		}
	
		K.cgsolve() ;
		for(size_t i = 0 ; i < elements.size() ; i++)
			elements[i]->step(0., &K.getDisplacements()) ;
		
		double energy = 0 ;
		for(size_t i = 0 ; i < elements.size() ;i++)
		{
			if(elements[i]->getBehaviour()->getFractureCriterion())
			{
				energy += elements[i]->getState().elasticEnergy()*elements[i]->volume() ;
			}
		}
		for(size_t i = 0 ; i < elements.size() ;i++)
		{
			delete elements[i] ;
		}
	}
	
	
	return (originalenergy-energy)/(delta_d) ;
}

void FractureCriterion::initialiseCache(const ElementState & s)
{
	DelaunayTriangle * testedTri = dynamic_cast<DelaunayTriangle *>(s.getParent()) ;
	DelaunayTetrahedron * testedTet = dynamic_cast<DelaunayTetrahedron *>(s.getParent()) ;
	if(testedTri)
	{
		if(!cache.empty())
		{
			cache.clear();
			area.clear();
		}
		Circle epsilon(neighbourhoodradius,testedTri->getCenter()) ;
		if(!testedTri->tree)
			return ;
		mesh2d = &testedTri->tree->getTree() ;
		std::vector<DelaunayTriangle *> tempcache = testedTri->tree->getNeighbouringElementsInGeometry(testedTri, &epsilon);
		std::vector<DelaunayTriangle *> neighbourhood ;
		for(size_t i = 0 ; i < testedTri->neighbourhood.size() ; i++)
		{
			if(testedTri->getNeighbourhood(i)->getBehaviour() && 
				testedTri->getNeighbourhood(i)->getBehaviour()->type != VOID_BEHAVIOUR && 
				testedTri->getNeighbourhood(i)->getBehaviour()->getFractureCriterion())
			{
				neighbourhood.push_back(testedTri->getNeighbourhood(i));
				cache.push_back(testedTri->getNeighbourhood(i)->index);
				area.push_back(testedTri->getNeighbourhood(i)->area());
			}
		}
		
		
		for(size_t i = 0 ; i < tempcache.size() ; i++)
		{
			if(tempcache[i]->getBehaviour() && 
				tempcache[i]->getBehaviour()->type != VOID_BEHAVIOUR && 
				tempcache[i]->getBehaviour()->getFractureCriterion())
			{
				bool inNeighbourhood = false ;
				for(size_t j = 0 ; j < neighbourhood.size() ; j++)
				{
					if(neighbourhood[j] == tempcache[i])
					{
						inNeighbourhood = true ;
						break ;
					}
				}
				if(!inNeighbourhood)
				{
					cache.push_back(tempcache[i]->index);
					area.push_back(tempcache[i]->area());
				}
			}
		}
	}
	else if(testedTet)
	{
		if(!cache.empty())
			cache.clear();
		
		Sphere epsilon(neighbourhoodradius,testedTet->getCenter()) ;
		if(!testedTet->tree)
			return ;
		mesh3d = &testedTet->tree->getTree() ;
		std::vector<DelaunayTetrahedron *> tempcache3d = testedTet->tree->getConflictingElements(&epsilon);
		std::vector<DelaunayTetrahedron *> neighbourhood ;
		for(size_t i = 0 ; i < testedTet->neighbourhood.size() ; i++)
		{
			if(testedTet->getNeighbourhood(i)->getBehaviour()
				&& testedTet->getNeighbourhood(i)->getBehaviour()->type != VOID_BEHAVIOUR  
				&& testedTet->getNeighbourhood(i)->getBehaviour()->getFractureCriterion())
			{
				neighbourhood.push_back(testedTet->getNeighbourhood(i));
				cache.push_back(testedTet->getNeighbourhood(i)->index);
				area.push_back(testedTet->getNeighbourhood(i)->volume());
			}
		}
		
		for(size_t i = 0 ; i < tempcache3d.size() ; i++)
		{
			if(tempcache3d[i]->getBehaviour()
				&& tempcache3d[i]->getBehaviour()->type != VOID_BEHAVIOUR 
				&& tempcache3d[i]->getBehaviour()->getFractureCriterion())
			{
				bool inNeighbourhood = false ;
				for(size_t j = 0 ; j < neighbourhood.size() ; j++)
				{
					if(neighbourhood[j] == tempcache3d[i])
					{
						inNeighbourhood = true ;
						break ;
					}
				}
				if(!inNeighbourhood)
				{
					cache.push_back(tempcache3d[i]->index);
					area.push_back(tempcache3d[i]->volume());
				}
			}
		}
	}
}

std::pair<double, double> FractureCriterion::getDeltaEnergyDeltaCriterion(const ElementState & s, double delta_d) const
{
	Assembly K ;
	double originalscore = 0 ;
	double originalenergy = 0 ;
	double score = 0 ;
	double energy = 0 ;
	double volume = 0 ;

	if(s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
	{
		Circle c(getMaterialCharacteristicRadius(), s.getParent()->getCenter()) ;

		std::vector<DelaunayTriangle *> elements ;
		std::vector<LinearForm *> behaviours ;
		TriElement father(LINEAR) ;

		for(size_t i = 0 ; i < cache.size() ; i++)
		{
			DelaunayTriangle * tric = static_cast<DelaunayTriangle * > ((*mesh2d)[cache[i]]) ;
			elements.push_back(new DelaunayTriangle(NULL, NULL,   tric->first,  tric->second,   tric->third,  NULL) );
			elements.back()->setBehaviour(tric->getBehaviour()->getCopy()) ;
			elements.back()->refresh(&father);
			elements.back()->getState().initialize() ;
			K.add(elements.back());
			volume += elements.back()->area() ;
			for(size_t j = 0 ; j < tric->neighbour.size() ; j++)
			{
				if(tric->getNeighbour(j)->isTriangle )
				{
					DelaunayTriangle * tri = dynamic_cast<DelaunayTriangle *>(tric->getNeighbour(j)) ;
					if(tri->getBehaviour()->type != VOID_BEHAVIOUR)
					{
						for(size_t k = 0 ; k <  tri->getBoundingPoints().size() ; k++)
						{
							if(!c.in(tri->getBoundingPoint(k)))
							{
								for(size_t l = 0 ; l <  tri->getBoundingPoints().size() ; l++)
								{
									int id = tri->getBoundingPoint(k).id ;
									double ex = tri->getState().getDisplacements()[l*2];
									double ey = tri->getState().getDisplacements()[l*2+1];
									K.setPoint(ex, ey ,id);
									break ;
								}
							}
						}
					}
				}
			}
		}
		K.cgsolve() ;
		for(size_t i = 0 ; i < elements.size() ; i++)
			elements[i]->step(0., &K.getDisplacements()) ;
		for(size_t i = 0 ; i < elements.size() ;i++)
		{
			if(elements[i]->getBehaviour()->getFractureCriterion())
			{
				double a = elements[i]->area() ;
				double sc = elements[i]->getBehaviour()->getFractureCriterion()->grade(elements[i]->getState()) ;
				if(dist(elements[i]->getCenter(), s.getParent()->getCenter()) < POINT_TOLERANCE_2D)
					originalscore += 1./(1.-sc) ;
				originalenergy += elements[i]->getState().elasticEnergy()*a ;
			}
		}
		
		std::valarray<Point *> nularray(0) ;
		
		for(size_t i = 0 ; i < elements.size() ;i++)
		{
			elements[i]->setBoundingPoints(nularray) ;
			delete elements[i] ;
		}
		
	}
	else
	{
		Sphere c(getMaterialCharacteristicRadius(), s.getParent()->getCenter()) ;
		
		std::vector<DelaunayTetrahedron *> elements ;
		std::vector<LinearForm *> behaviours ;
		TetrahedralElement father ;
		for(size_t i = 0 ; i < cache.size() ; i++)
		{
			DelaunayTetrahedron * tet = dynamic_cast<DelaunayTetrahedron *>((*mesh3d)[cache[i]]) ;
			elements.push_back(new DelaunayTetrahedron(NULL, NULL,   tet->first,  tet->second,   tet->third, tet->fourth,  NULL) );
			elements.back()->setBehaviour(tet->getBehaviour()->getCopy()) ;
			
			elements.back()->refresh(&father);
			elements.back()->getState().initialize() ;
			K.add(elements.back());
			volume+= elements.back()->volume() ;
			
			for(size_t j = 0 ; j < tet->neighbour.size() ; j++)
			{
				if(tet->getNeighbour(j)->isTetrahedron() )
				{
					DelaunayTetrahedron * tri = dynamic_cast<DelaunayTetrahedron *>(tet->getNeighbour(j)) ;
					if(tri->getBehaviour()->type != VOID_BEHAVIOUR)
					{
						for(size_t k = 0 ; k <  tri->getBoundingPoints().size() ; k++)
						{
							if(!c.in(tri->getBoundingPoint(k)))
							{
								for(size_t l = 0 ; l <  tri->getBoundingPoints().size() ; l++)
								{
									int id = tri->getBoundingPoint(k).id ;
									double ex = tri->getState().getDisplacements()[l*3];
									double ey = tri->getState().getDisplacements()[l*3+1];
									double ez = tri->getState().getDisplacements()[l*3+2];
									K.setPoint(ex, ey, ez ,id);
									break ;
								}
							}
						}
					}
				}
			}
		
			
		}
	
		K.cgsolve() ;
		for(size_t i = 0 ; i < elements.size() ; i++)
			elements[i]->step(0., &K.getDisplacements()) ;
		
		double energy = 0 ;
		for(size_t i = 0 ; i < elements.size() ;i++)
		{
			if(elements[i]->getBehaviour()->getFractureCriterion())
			{
				double v = elements[i]->volume() ;
				double sc = elements[i]->getBehaviour()->getFractureCriterion()->getScoreAtState() ;
				if(dist(elements[i]->getCenter(), s.getParent()->getCenter()) < POINT_TOLERANCE_3D)
					originalscore += 1./(1.-sc) ;
				originalenergy += elements[i]->getState().elasticEnergy()*v ;
			}
		}
		for(size_t i = 0 ; i < elements.size() ;i++)
		{
			delete elements[i] ;
		}
	}
	
	
	
	if(s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
	{
		Circle c(getMaterialCharacteristicRadius(), s.getParent()->getCenter()) ;

		std::vector<DelaunayTriangle *> elements ;
		std::vector<LinearForm *> behaviours ;
		TriElement father(LINEAR) ;

		for(size_t i = 0 ; i < cache.size() ; i++)
		{
			DelaunayTriangle * tric = static_cast<DelaunayTriangle * > ((*mesh2d)[cache[i]]) ;
			elements.push_back(new DelaunayTriangle(NULL, NULL,   tric->first,   tric->second,    tric->third,  NULL) );
			elements.back()->setBehaviour(tric->getBehaviour()->getCopy()) ;

			if(tric == s.getParent())
				elements.back()->getBehaviour()->setTensor(elements.back()->getBehaviour()->getTensor(elements.back()->getCenter())*(1.-delta_d)) ;

			elements.back()->getBehaviour()->getFractureCriterion()->setEnergyIndexed(false) ;
			elements.back()->refresh(&father);
			elements.back()->getState().initialize() ;
			K.add(elements.back());
			volume += elements.back()->area() ;
			for(size_t j = 0 ; j <  tric->neighbour.size() ; j++)
			{
				if( tric->getNeighbour(j)->isTriangle )
				{
					DelaunayTriangle * tri = dynamic_cast<DelaunayTriangle *>( tric->getNeighbour(j)) ;
					if(tri->getBehaviour()->type != VOID_BEHAVIOUR)
					{
						for(size_t k = 0 ; k <  tri->getBoundingPoints().size() ; k++)
						{
							if(!c.in(tri->getBoundingPoint(k)))
							{
								for(size_t l = 0 ; l <  tri->getBoundingPoints().size() ; l++)
								{
									int id = tri->getBoundingPoint(k).id ;
									double ex = tri->getState().getDisplacements()[l*2];
									double ey = tri->getState().getDisplacements()[l*2+1];
									K.setPoint(ex, ey ,id);
									break ;
								}
							}
						}
					}
				}
			}
		}
	
		K.cgsolve() ;
		for(size_t i = 0 ; i < elements.size() ; i++)
			elements[i]->step(0., &K.getDisplacements()) ;

		for(size_t i = 0 ; i < elements.size() ;i++)
		{
			if(elements[i]->getBehaviour()->getFractureCriterion())
			{
				double a = elements[i]->area() ;
				double sc = elements[i]->getBehaviour()->getFractureCriterion()->grade(elements[i]->getState()) ;
				if(dist(elements[i]->getCenter(), s.getParent()->getCenter()) < POINT_TOLERANCE_2D)
					score += 1./(1.-sc);
				energy += elements[i]->getState().elasticEnergy()*a ;
			}
		}

		std::valarray<Point *> nularray(0) ;
		
		for(size_t i = 0 ; i < elements.size() ;i++)
		{
			elements[i]->setBoundingPoints(nularray) ;
			delete elements[i] ;
		}
	}
	else
	{
		Sphere c(getMaterialCharacteristicRadius(), s.getParent()->getCenter()) ;
		
		std::vector<DelaunayTetrahedron *> elements ;
		std::vector<LinearForm *> behaviours ;
		TetrahedralElement father ;
		for(size_t i = 0 ; i < cache.size() ; i++)
		{
			DelaunayTetrahedron * tet = static_cast<DelaunayTetrahedron * > ((*mesh3d)[cache[i]]) ;
			elements.push_back(new DelaunayTetrahedron(NULL, NULL,   tet->first,  tet->second,   tet->third, tet->fourth,  NULL) );
			elements.back()->setBehaviour(tet->getBehaviour()->getCopy()) ;
			if(tet == s.getParent())
				elements.back()->getBehaviour()->setTensor(elements.back()->getBehaviour()->getTensor(elements.back()->getCenter())*(1.-delta_d)) ;
			elements.back()->getBehaviour()->getFractureCriterion()->setEnergyIndexed(false) ;
			elements.back()->refresh(&father);
			elements.back()->getState().initialize() ;
			K.add(elements.back());
			volume+= elements.back()->volume() ;
			
			for(size_t j = 0 ; j < tet->neighbour.size() ; j++)
			{
				if(tet->getNeighbour(j)->isTetrahedron() )
				{
					DelaunayTetrahedron * tri = dynamic_cast<DelaunayTetrahedron *>(tet->getNeighbour(j)) ;
					if(tri->getBehaviour()->type != VOID_BEHAVIOUR)
					{
						for(size_t k = 0 ; k <  tri->getBoundingPoints().size() ; k++)
						{
							if(!c.in(tri->getBoundingPoint(k)))
							{
								int id = tri->getBoundingPoint(k).id ;
								double ex = tri->getState().getDisplacements()[k*3];
								double ey = tri->getState().getDisplacements()[k*3+1];
								double ez = tri->getState().getDisplacements()[k*3+2];
								K.setPoint(ex, ey, ez ,id);
							}
						}
					}
				}
			}
		
			
		}
	
		K.cgsolve() ;
		for(size_t i = 0 ; i < elements.size() ; i++)
			elements[i]->step(0., &K.getDisplacements()) ;
		
		double energy = 0 ;
		for(size_t i = 0 ; i < elements.size() ;i++)
		{
			if(elements[i]->getBehaviour()->getFractureCriterion())
			{
				double sc = elements[i]->getBehaviour()->getFractureCriterion()->getScoreAtState() ;
				if(dist(elements[i]->getCenter(), s.getParent()->getCenter()) < POINT_TOLERANCE_3D)
					score += 1./(1.-sc) ;
				energy += elements[i]->getState().elasticEnergy()*elements[i]->volume() ;
			}
		}
		for(size_t i = 0 ; i < elements.size() ;i++)
		{
			delete elements[i] ;
		}
	}
	
	
	return std::make_pair( (energy-originalenergy)/(delta_d),(score-originalscore)/(delta_d)) ;
}

std::pair<bool, bool> FractureCriterion::inSetAndSetChanged(int fractiles, const ElementState &s) 
{
	stable = true ;
	if( s.getParent()->getBehaviour()->getDamageModel() == NULL )
		return std::make_pair(false, false) ;
	
	DelaunayTriangle * testedTri = dynamic_cast<DelaunayTriangle *>(s.getParent()) ;
	DelaunayTetrahedron * testedTet = dynamic_cast<DelaunayTetrahedron *>(s.getParent()) ;
	HexahedralElement * testedHex = dynamic_cast<HexahedralElement *>(s.getParent()) ;
	if(testedTri)
	{
		if(cache.size() == 0)
			std::make_pair(false, false) ;
		
		unsigned int idx = testedTri->index ;
		bool inset = false ;
		
		std::vector<unsigned int> newSet ;
		for(size_t i = 0 ; i< cache.size() ; i++)
		{
			DelaunayTriangle * ci = static_cast<DelaunayTriangle *>((*mesh2d)[cache[i]]) ;
			if( ci->getBehaviour()->getFractureCriterion() 
				&& ci->getBehaviour()->getFractureCriterion()->metAtStep)
			{
				if(ci->index == idx)
					inset = true ;
				
				newSet.push_back(cache[i]) ;
			}
		}
		std::stable_sort(newSet.begin(), newSet.end()) ;
		
		if(s.getDeltaTime() > POINT_TOLERANCE_2D) //new iteration
		{
			damagingSet = newSet ;
			return std::make_pair(inset, true) ;
		}
		else
		{
			if(damagingSet.empty())
				return std::make_pair(false, false) ;
			
			bool identicalSets = (damagingSet.size() == newSet.size()) ;
			bool allConverged = true ;
			for(size_t i = 0 ; i< damagingSet.size(); i++)
			{
				DelaunayTriangle * ci = static_cast<DelaunayTriangle *>((*mesh2d)[damagingSet[i]]) ;
				if(ci->getBehaviour()->getDamageModel() 
					&& !ci->getBehaviour()->getDamageModel()->converged)
				{
					allConverged = false ;
				}
				if(ci->index == idx)
					inset = true ;
			}
			
			if(identicalSets)
			{
				for(size_t i = 0 ; i< newSet.size() ; i++)
				{
					if(newSet[i] != damagingSet[i])
					{
						identicalSets = false ;
						break ;
					}
				}
			}
			
			if(allConverged && !identicalSets) // checkpoint: we work on a new set
			{
				damagingSet = newSet ;
				inset = false ;
				for(size_t i = 0 ; i< damagingSet.size(); i++)
				{
					DelaunayTriangle * ci = static_cast<DelaunayTriangle *>((*mesh2d)[damagingSet[i]]) ;
					if(damagingSet[i] == idx && metAtStep)
					{
						inset = true ;
					}
				}
				
				return std::make_pair(inset, false) ;
			}
			
			if(allConverged && identicalSets)
			{
				std::cout << "loss of stability." << std::endl ;
				stable = false ;
				return std::make_pair(inset, false) ;
			}
			
			
			return std::make_pair(inset, identicalSets) ;
		}
	}
	if(testedTet)
	{
		if(cache.size() == 0)
			std::make_pair(false, false) ;
		
		unsigned int idx = testedTet->index ;
		bool inset = false ;
		
		std::vector<unsigned int> newSet ;
		for(size_t i = 0 ; i< cache.size() ; i++)
		{
			DelaunayTetrahedron * ci = static_cast<DelaunayTetrahedron *>((*mesh3d)[cache[i]]) ;
			if( ci->getBehaviour()->getFractureCriterion() 
				&& ci->getBehaviour()->getFractureCriterion()->metAtStep)
			{
				if(ci->index == idx)
					inset = true ;
				
				newSet.push_back(cache[i]) ;
			}
		}
		std::stable_sort(newSet.begin(), newSet.end()) ;
		
		if(s.getDeltaTime() > POINT_TOLERANCE_3D) //new iteration
		{
			damagingSet = newSet ;
			return std::make_pair(inset, true) ;
		}
		else
		{
			bool identicalSets = (damagingSet.size() == newSet.size()) ;
			bool allConverged = true ;
			for(size_t i = 0 ; i< damagingSet.size(); i++)
			{
				DelaunayTetrahedron * ci = static_cast<DelaunayTetrahedron *>((*mesh3d)[damagingSet[i]]) ;
				if(ci->getBehaviour()->getDamageModel() 
					&& !ci->getBehaviour()->getDamageModel()->converged)
				{
					allConverged = false ;
				}
				if(ci->index == idx)
					inset = true ;
			}
			
			if(identicalSets)
			{
				for(size_t i = 0 ; i< newSet.size() ; i++)
				{
					if(newSet[i] != damagingSet[i])
					{
						identicalSets = false ;
						break ;
					}
				}
			}
			
			if(allConverged && identicalSets) // checkpoint: we work on a new set
			{
				damagingSet = newSet ;
				inset = false ;
				for(size_t i = 0 ; i< damagingSet.size(); i++)
				{
					DelaunayTriangle * ci = static_cast<DelaunayTriangle *>((*mesh2d)[damagingSet[i]]) ;
					if(damagingSet[i] == idx && metAtStep)
					{
						inset = true ;
					}
				}
				
				return std::make_pair(inset, false) ;
			}
			
			return std::make_pair(inset, identicalSets) ;
		}
	}
// 	else if(testedHex)
// 	{
// 		return std::make_pair(false, false) ;
// 	}
	else
	{
		std::cout << " criterion not implemented for this kind of element" << std::endl ;
		return std::make_pair(false, false) ;
	}
	
	//shut up the compiler
	return std::make_pair(false, false) ;
}

void FractureCriterion::step(const ElementState &s)
{
		if(cache.empty() )
				initialiseCache(s) ;
		
	scoreAtState = grade(s) ;
	if(energyIndexed && s.getDeltaTime() > POINT_TOLERANCE_2D )
		noEnergyUpdate = true ;
	
	if(energyIndexed && noEnergyUpdate /*&& scoreAtState > -.5*/)
	{
		noEnergyUpdate = false ;
		currentEnergy = 0 ;
		bool buildvolume = false ;
		if(neighbourhoodvolume < 0)
		{
			neighbourhoodvolume = 0 ;
			buildvolume = true ;
		}
	
		if(mesh2d)
		{
			Circle c(getMaterialCharacteristicRadius(), s.getParent()->getCenter()) ;
			for(size_t i = 0 ; i < getCache().size() ; i++)
			{
				DelaunayTriangle *tri = static_cast<DelaunayTriangle *>( (*mesh2d)[getCache()[i]] ) ;
				if(c.in(*tri->first) || c.in(*tri->second) || c.in(*tri->third))
				{
					currentEnergy += tri->getState().elasticEnergy()*tri->area() ;
					if(buildvolume)
						neighbourhoodvolume += tri->area() ;
				}
			}
		}
		else if(mesh3d)
		{
			Sphere c(getMaterialCharacteristicRadius(), s.getParent()->getCenter()) ;
			
			for(size_t i = 0 ; i < getCache().size() ; i++)
			{
				DelaunayTetrahedron *tet = static_cast<DelaunayTetrahedron *>( (*mesh3d)[getCache()[i]] ) ;
				if(c.in(*tet->first) || c.in(*tet->second) || c.in(*tet->third) || c.in(*tet->fourth))
				{
					currentEnergy += tet->getState().elasticEnergy()*tet->volume() ;
					if(buildvolume)
						neighbourhoodvolume += tet->volume() ;
				}
			}
		}
		deltaEnergyAtState = (currentEnergy-previousEnergy) ;
		
		previousEnergy = currentEnergy ;
		
		return ;
	}
	
	if(energyIndexed /*&& scoreAtState > -.5*/)
	{
		std::pair<double, double> dedc = getDeltaEnergyDeltaCriterion(s, 0.0001) ;
		energyDamageDifferential = dedc.first ;
		criterionDamageDifferential = dedc.second ;
	}

}

void FractureCriterion::computeNonLocalState(const ElementState &s)
{
	if( s.getParent()->getBehaviour()->getDamageModel() == NULL )
	{
		metAtStep = false ;
		return  ;
	}
	if( s.getParent()->getBehaviour()->getDamageModel() && s.getParent()->getBehaviour()->getDamageModel()->fractured())
	{
		metAtStep = false ;
		return  ;
	}
	
	double tol = 1e-6 ;
	DelaunayTriangle * testedTri = dynamic_cast<DelaunayTriangle *>(s.getParent()) ;
	DelaunayTetrahedron * testedTet = dynamic_cast<DelaunayTetrahedron *>(s.getParent()) ;
	HexahedralElement * testedHex = dynamic_cast<HexahedralElement *>(s.getParent()) ;
	if(testedTri)
	{
		if (scoreAtState < 0)
		{
			metAtStep = false ;
			return  ;
		}
		double maxNeighbourhoodScore = 0 ;
		double matchedArea = 0 ;
		std::map<double, DelaunayTriangle *> scores ;
		std::vector<double> unsortedScores ;
		std::map<DelaunayTriangle *, double> areatemp ;
		DelaunayTriangle * maxLocus = NULL;
		double areamax = 0 ;
		if(!cache.empty())
		{
			for(size_t i = 0 ; i< cache.size() ; i++)
			{
				DelaunayTriangle * ci = static_cast<DelaunayTriangle *>((*mesh2d)[cache[i]]) ;

				areamax += area[i] ;

				double s = 0. ;
				if(ci->getBehaviour()->getFractureCriterion() != NULL && !ci->getBehaviour()->fractured())
					s = ci->getBehaviour()->getFractureCriterion()->getSteppedScore() ;

				scores[-s] =  ci;
				unsortedScores.push_back(s);
				if(s > maxNeighbourhoodScore)
				{
					maxNeighbourhoodScore = s ;
					maxLocus = ci ;
				}

				areatemp[ci] = area[i] ;
			}
		}
		
		if(maxLocus == NULL || maxNeighbourhoodScore < 0)
		{
			metAtStep = false ;
			return  ;
		}
		
		bool nearmaxlocus = false;
		
		for(size_t i = 0 ; i< cache.size() ; i++)
		{
			DelaunayTriangle * ci = static_cast<DelaunayTriangle *>((*mesh2d)[cache[i]]) ;
			

			if(ci->getBehaviour()->getFractureCriterion() && (ci->getBehaviour()->getFractureCriterion()->getSteppedScore() >= maxNeighbourhoodScore-tol || ci->getBehaviour()->fractured()))
			{
				if(squareDist2D(ci->getCenter(), s.getParent()->getCenter()) < physicalCharacteristicRadius*physicalCharacteristicRadius)
				{
					nearmaxlocus = true ;
					break ;
				}
			}
		}
		
		bool foundcutoff = false ;
		double thresholdscore = maxNeighbourhoodScore ;
		double avgscore = 0 ;
		double trialarea = std::min(physicalCharacteristicRadius*physicalCharacteristicRadius*M_PI, areamax) ;
		for(auto i = scores.begin() ; i != scores.end() ; ++i)
		{

			double parea = matchedArea ;
			double narea = areatemp[i->second] ;
			matchedArea += narea ;
			double a(narea) ;
			if(mirroring == MIRROR_X && std::abs(i->second->getCenter().x  - delta_x) < physicalCharacteristicRadius) // MIRROR_X
			{
				matchedArea += a ;
				narea += a ;
			}
			if(mirroring == MIRROR_Y &&  std::abs(i->second->getCenter().y  - delta_y) < physicalCharacteristicRadius) // MIRROR_Y
			{
				matchedArea += a ;
				narea += a ;
			}
			if(mirroring == MIRROR_XY &&  std::abs(i->second->getCenter().x  - delta_x) < physicalCharacteristicRadius) // MIRROR_XY
			{
				matchedArea += a ;
				narea += a ;
			}
			if(mirroring == MIRROR_XY &&  std::abs(i->second->getCenter().y  - delta_y) < physicalCharacteristicRadius) // MIRROR_XY
			{
				matchedArea += a ;
				narea += a ;
			}
			avgscore = (parea*avgscore - narea*i->first)/matchedArea ;

			if (avgscore <= 0 && matchedArea >= trialarea)
			{
				thresholdscore = -i->first ;
				foundcutoff = true ;
				break ;
			}
			else if (avgscore < 0 && matchedArea < trialarea)
			{
				foundcutoff = false ;
				break ;
			}
			else if (avgscore > 0 && matchedArea >= trialarea)
			{
				thresholdscore = -i->first ;
				foundcutoff = true ;
				break ;
			}
// 			else
// 				i->second->visited = true ;
		}
		
		if (!foundcutoff && areamax > s.getParent()->area())
		{
			metAtStep = false ;
			return  ;
		}
		if (nearmaxlocus)
		{
			metAtStep = true ;
			return  ;
		}
		metAtStep = false ;
		return  ;

	}
	if(testedTet)
	{

		if(testedTet->visited())
		{
			metAtStep = false ;
			return  ;
		}
				
		if (scoreAtState <= 0)
		{
			metAtStep = false ;
			return  ;
		}

		double maxNeighbourhoodScore = 0 ;
		double matchedArea = 0 ;
		std::map<double, DelaunayTetrahedron *> scores ;
		std::map<DelaunayTetrahedron *, double> areatemp ;
		DelaunayTetrahedron * maxLocus = NULL;
		
		if(!cache.empty())
		{
			for(size_t i = 0 ; i< cache.size() ; i++)
			{
				DelaunayTetrahedron * ci = static_cast<DelaunayTetrahedron *>((*mesh3d)[cache[i]]) ;

				if( !ci->getBehaviour()->fractured())
				{
					double s = ci->getBehaviour()->getFractureCriterion()->getSteppedScore() ;
					scores[-s] =  ci;
					if(s > maxNeighbourhoodScore)
					{
						maxNeighbourhoodScore = s ;
						maxLocus = ci ;
					}
				}
				else if(ci->getBehaviour()->fractured())
				{
					double s = POINT_TOLERANCE_2D ;
					scores[-s] =  ci;
				}
				areatemp[ci] = area[i] ;

			}
		}
		
		if(!maxLocus)
		{
			metAtStep = false ;
			return  ;
		}
		
		std::vector<DelaunayTetrahedron *> maxloci ;
		
		for(size_t i = 0 ; i< cache.size() ; i++)
		{
			DelaunayTetrahedron * ci = static_cast<DelaunayTetrahedron *>((*mesh3d)[cache[i]]) ;
			if(ci->getBehaviour()->getFractureCriterion())
				if(std::abs(ci->getBehaviour()->getFractureCriterion()->getSteppedScore()-maxNeighbourhoodScore) < tol || ci->getBehaviour()->fractured())
					maxloci.push_back(ci) ;
		}
		
		bool foundcutoff = false ;
		double thresholdscore = maxNeighbourhoodScore ;
		
		for(auto i = scores.begin() ; i != scores.end() ; ++i)
		{
			
			if(!foundcutoff)
			{
				if(-i->first > 0 )
				{
					matchedArea += areatemp[i->second] ;
					if(mirroring == MIRROR_X && std::abs(i->second->getCenter().x  - delta_x) < physicalCharacteristicRadius) // MIRROR_X
					  matchedArea += areatemp[i->second] ;
					if(mirroring == MIRROR_Y &&  std::abs(i->second->getCenter().y  - delta_y) < physicalCharacteristicRadius) // MIRROR_Y
					  matchedArea += areatemp[i->second] ;
					if(mirroring == MIRROR_Z &&  std::abs(i->second->getCenter().z  - delta_z) < physicalCharacteristicRadius) // MIRROR_Y
					  matchedArea += areatemp[i->second] ;
					if(mirroring == MIRROR_XY &&  std::abs(i->second->getCenter().x  - delta_x) < physicalCharacteristicRadius) // MIRROR_XY
					  matchedArea += areatemp[i->second] ;
					if(mirroring == MIRROR_XY &&  std::abs(i->second->getCenter().y  - delta_y) < physicalCharacteristicRadius) // MIRROR_XY
					  matchedArea += areatemp[i->second] ;
					if(mirroring == MIRROR_XZ &&  std::abs(i->second->getCenter().x  - delta_x) < physicalCharacteristicRadius) // MIRROR_XY
					  matchedArea += areatemp[i->second] ;
					if(mirroring == MIRROR_XZ &&  std::abs(i->second->getCenter().z  - delta_z) < physicalCharacteristicRadius) // MIRROR_XY
					  matchedArea += areatemp[i->second] ;
					if(mirroring == MIRROR_YZ &&  std::abs(i->second->getCenter().y  - delta_y) < physicalCharacteristicRadius) // MIRROR_XY
					  matchedArea += areatemp[i->second] ;
					if(mirroring == MIRROR_YZ &&  std::abs(i->second->getCenter().z  - delta_z) < physicalCharacteristicRadius) // MIRROR_XY
					  matchedArea += areatemp[i->second] ;
				}
				if (matchedArea > 1.333333333333333333*physicalCharacteristicRadius*physicalCharacteristicRadius*physicalCharacteristicRadius*M_PI)
				{
					thresholdscore = -i->first ;
					foundcutoff  = true ;
					break ;
				}
			}
// 			else
// 				i->second->visited = true ;
		}
		if (!foundcutoff )
		{
			metAtStep = false ;
			return  ;
		}

		for(size_t i = 0 ; i < maxloci.size() ; i++)
			if(squareDist3D(maxloci[i]->getCenter(), s.getParent()->getCenter()) < physicalCharacteristicRadius*physicalCharacteristicRadius)
		{
			metAtStep = true ;
			return  ;
		}
		
		metAtStep = false ;
		return  ;
	}
	else if(testedHex)
	{
		std::set<HexahedralElement *> neighbourhood ;
		std::vector<HexahedralElement *> neighbours = testedHex->neighbourhood ;
		for(size_t i = 0 ; i < neighbours.size() ; i++)
		{
			for(size_t j = 0 ; j <  neighbours[i]->neighbourhood.size() ; j++)
			{
				if(neighbours[i]->neighbourhood[j] != testedHex 
				   && !neighbours[i]->neighbourhood[j]->getBehaviour()->fractured())
					neighbourhood.insert(neighbours[i]->neighbourhood[j]) ;
			}
		}
		double score = grade(s) ;
		double maxNeighbourhoodScore = 0 ;
		if(!neighbourhood.empty())
		{
			for(auto i= neighbourhood.begin() ; i != neighbourhood.end() ; ++i)
			{
				if((*i)->getBehaviour()->getFractureCriterion() 
					&& !(*i)->getBehaviour()->fractured())
					maxNeighbourhoodScore = std::max(maxNeighbourhoodScore,
					                                 (*i)->getBehaviour()->getFractureCriterion()->grade((*i)->getState())) ;
				if((*i)->getBehaviour()->changed())
				{
					maxNeighbourhoodScore = 10.*score ;
					break ;
				}
				
				if (maxNeighbourhoodScore > score)
					break ;
				
			}
		}
		
		if( score > 0 )
		{
			if(score > maxNeighbourhoodScore)
			{
				metAtStep = true ;
				return  ;
			}
		}

		metAtStep = false ;
		return  ;
	}
	else
	{
		std::cout << " criterion not implemented for this kind of element" << std::endl ;
		metAtStep = false ;
		return  ;
	}
	
	//shut up the compiler
		metAtStep = false ;
		return  ;
}

FractureCriterion::~FractureCriterion()
{
}

void FractureCriterion::setNeighbourhoodRadius(double r)
{
	neighbourhoodradius = r ;
	cache.clear() ;
}

void FractureCriterion::setMaterialCharacteristicRadius(double r)
{
	physicalCharacteristicRadius = r ;
}

bool FractureCriterion::met(const ElementState &s)
{
	return metAtStep ;
}

Material FractureCriterion::toMaterial()
{
	Material mat ;
	return mat ;
}


