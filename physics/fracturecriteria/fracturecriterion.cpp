//
// C++ Implementation: fracturecriterion
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
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

FractureCriterion::FractureCriterion(MirrorState mirroring, double delta_x, double delta_y, double delta_z) : neighbourhoodradius(.0005), neighbourhoodvolume(-1), physicalCharacteristicRadius(.008), scoreAtState(0), metInTension(false), metInCompression(false), mirroring(mirroring), delta_x(delta_x), delta_y(delta_y), delta_z(delta_z), deltaScoreAtState(0), deltaEnergyAtState(0), energyDamageDifferential(0), criterionDamageDifferential(0), energyIndexed(false), noEnergyUpdate(true)
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

			elements.push_back(new DelaunayTriangle(NULL, NULL,   cache[i]->first,  cache[i]->second,   cache[i]->third,  NULL) );
			elements.back()->setBehaviour(cache[i]->getBehaviour()->getCopy()) ;
			elements.back()->refresh(&father);
			elements.back()->getState().initialize() ;
			K.add(elements.back());
			volume += elements.back()->area() ;
			for(size_t j = 0 ; j < cache[i]->neighbour.size() ; j++)
			{
				if(cache[i]->getNeighbour(j)->isTriangle )
				{
					DelaunayTriangle * tri = dynamic_cast<DelaunayTriangle *>(cache[i]->getNeighbour(j)) ;
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
		for(size_t i = 0 ; i < cache3d.size() ; i++)
		{

			elements.push_back(new DelaunayTetrahedron(NULL, NULL,   cache3d[i]->first,  cache3d[i]->second,   cache3d[i]->third, cache3d[i]->fourth,  NULL) );
			elements.back()->setBehaviour(cache3d[i]->getBehaviour()->getCopy()) ;
			
			elements.back()->refresh(&father);
			elements.back()->getState().initialize() ;
			K.add(elements.back());
			volume+= elements.back()->volume() ;
			
			for(size_t j = 0 ; j < cache3d[i]->neighbour.size() ; j++)
			{
				if(cache3d[i]->getNeighbour(j)->isTetrahedron() )
				{
					DelaunayTetrahedron * tri = dynamic_cast<DelaunayTetrahedron *>(cache3d[i]->getNeighbour(j)) ;
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

			elements.push_back(new DelaunayTriangle(NULL, NULL,   cache[i]->first,  cache[i]->second,   cache[i]->third,  NULL) );
			elements.back()->setBehaviour(cache[i]->getBehaviour()->getCopy()) ;

			if(cache[i] == s.getParent())
				elements.back()->getBehaviour()->setTensor(elements.back()->getBehaviour()->getTensor(elements.back()->getCenter())*(1.-delta_d)) ;

			elements.back()->getBehaviour()->getFractureCriterion()->setEnergyIndexed(false) ;
			elements.back()->refresh(&father);
			elements.back()->getState().initialize() ;
			K.add(elements.back());
			volume += elements.back()->area() ;
			for(size_t j = 0 ; j < cache[i]->neighbour.size() ; j++)
			{
				if(cache[i]->getNeighbour(j)->isTriangle )
				{
					DelaunayTriangle * tri = dynamic_cast<DelaunayTriangle *>(cache[i]->getNeighbour(j)) ;
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
		for(size_t i = 0 ; i < cache3d.size() ; i++)
		{
			
			elements.push_back(new DelaunayTetrahedron(NULL, NULL,   cache3d[i]->first,  cache3d[i]->second,   cache3d[i]->third, cache3d[i]->fourth,  NULL) );
			elements.back()->setBehaviour(cache3d[i]->getBehaviour()->getCopy()) ;
			if(cache3d[i] == s.getParent())
				elements.back()->getBehaviour()->setTensor(elements.back()->getBehaviour()->getTensor(elements.back()->getCenter())*(1.-delta_d)) ;
			elements.back()->getBehaviour()->getFractureCriterion()->setEnergyIndexed(false) ;
			elements.back()->refresh(&father);
			elements.back()->getState().initialize() ;
			K.add(elements.back());
			volume+= elements.back()->volume() ;
			
			for(size_t j = 0 ; j < cache3d[i]->neighbour.size() ; j++)
			{
				if(cache3d[i]->getNeighbour(j)->isTetrahedron() )
				{
					DelaunayTetrahedron * tri = dynamic_cast<DelaunayTetrahedron *>(cache3d[i]->getNeighbour(j)) ;
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
			return ;
		Circle epsilon(neighbourhoodradius,testedTri->getCenter()) ;
		if(!testedTri->tree)
			return ;
		std::vector<DelaunayTriangle *> tempcache = testedTri->tree->getConflictingElements(&epsilon);
		std::vector<DelaunayTriangle *> neighbourhood ;
		for(size_t i = 0 ; i < testedTri->neighbourhood.size() ; i++)
		{
			if(testedTri->getNeighbourhood(i)->getBehaviour() && testedTri->getNeighbourhood(i)->getBehaviour()->type != VOID_BEHAVIOUR)
			{
				neighbourhood.push_back(testedTri->getNeighbourhood(i));
				cache.push_back(testedTri->getNeighbourhood(i));
			}
		}
		
		
		for(size_t i = 0 ; i < tempcache.size() ; i++)
		{
			if(tempcache[i]->getBehaviour() && tempcache[i]->getBehaviour()->type != VOID_BEHAVIOUR)
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
					cache.push_back(tempcache[i]);
					area.push_back(cache.back()->area());
				}
			}
		}
	}
	else if(testedTet)
	{
		if(!cache3d.empty())
			return ;
		Sphere epsilon(neighbourhoodradius,testedTet->getCenter()) ;
		if(!testedTet->tree)
			return ;
		std::vector<DelaunayTetrahedron *> tempcache3d = testedTet->tree->getConflictingElements(&epsilon);
		std::vector<DelaunayTetrahedron *> neighbourhood ;
		for(size_t i = 0 ; i < testedTet->neighbourhood.size() ; i++)
		{
			if(testedTet->getNeighbourhood(i)->getBehaviour() && testedTet->getNeighbourhood(i)->getBehaviour()->type != VOID_BEHAVIOUR)
			{
				neighbourhood.push_back(testedTet->getNeighbourhood(i));
				cache3d.push_back(testedTet->getNeighbourhood(i));
			}
		}
		
		for(size_t i = 0 ; i < tempcache3d.size() ; i++)
		{
			if(tempcache3d[i]->getBehaviour()&& tempcache3d[i]->getBehaviour()->type != VOID_BEHAVIOUR)
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
					cache3d.push_back(tempcache3d[i]);
					area.push_back(cache3d.back()->volume());
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

			elements.push_back(new DelaunayTriangle(NULL, NULL,   cache[i]->first,  cache[i]->second,   cache[i]->third,  NULL) );
			elements.back()->setBehaviour(cache[i]->getBehaviour()->getCopy()) ;
			elements.back()->refresh(&father);
			elements.back()->getState().initialize() ;
			K.add(elements.back());
			volume += elements.back()->area() ;
			for(size_t j = 0 ; j < cache[i]->neighbour.size() ; j++)
			{
				if(cache[i]->getNeighbour(j)->isTriangle )
				{
					DelaunayTriangle * tri = dynamic_cast<DelaunayTriangle *>(cache[i]->getNeighbour(j)) ;
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
				if(dist(elements[i]->getCenter(), s.getParent()->getCenter()) < POINT_TOLERANCE)
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
		for(size_t i = 0 ; i < cache3d.size() ; i++)
		{

			elements.push_back(new DelaunayTetrahedron(NULL, NULL,   cache3d[i]->first,  cache3d[i]->second,   cache3d[i]->third, cache3d[i]->fourth,  NULL) );
			elements.back()->setBehaviour(cache3d[i]->getBehaviour()->getCopy()) ;
			
			elements.back()->refresh(&father);
			elements.back()->getState().initialize() ;
			K.add(elements.back());
			volume+= elements.back()->volume() ;
			
			for(size_t j = 0 ; j < cache3d[i]->neighbour.size() ; j++)
			{
				if(cache3d[i]->getNeighbour(j)->isTetrahedron() )
				{
					DelaunayTetrahedron * tri = dynamic_cast<DelaunayTetrahedron *>(cache3d[i]->getNeighbour(j)) ;
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
				if(dist(elements[i]->getCenter(), s.getParent()->getCenter()) < POINT_TOLERANCE)
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

			elements.push_back(new DelaunayTriangle(NULL, NULL,   cache[i]->first,  cache[i]->second,   cache[i]->third,  NULL) );
			elements.back()->setBehaviour(cache[i]->getBehaviour()->getCopy()) ;

			if(cache[i] == s.getParent())
				elements.back()->getBehaviour()->setTensor(elements.back()->getBehaviour()->getTensor(elements.back()->getCenter())*(1.-delta_d)) ;

			elements.back()->getBehaviour()->getFractureCriterion()->setEnergyIndexed(false) ;
			elements.back()->refresh(&father);
			elements.back()->getState().initialize() ;
			K.add(elements.back());
			volume += elements.back()->area() ;
			for(size_t j = 0 ; j < cache[i]->neighbour.size() ; j++)
			{
				if(cache[i]->getNeighbour(j)->isTriangle )
				{
					DelaunayTriangle * tri = dynamic_cast<DelaunayTriangle *>(cache[i]->getNeighbour(j)) ;
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
				if(dist(elements[i]->getCenter(), s.getParent()->getCenter()) < POINT_TOLERANCE)
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
		for(size_t i = 0 ; i < cache3d.size() ; i++)
		{
			
			elements.push_back(new DelaunayTetrahedron(NULL, NULL,   cache3d[i]->first,  cache3d[i]->second,   cache3d[i]->third, cache3d[i]->fourth,  NULL) );
			elements.back()->setBehaviour(cache3d[i]->getBehaviour()->getCopy()) ;
			if(cache3d[i] == s.getParent())
				elements.back()->getBehaviour()->setTensor(elements.back()->getBehaviour()->getTensor(elements.back()->getCenter())*(1.-delta_d)) ;
			elements.back()->getBehaviour()->getFractureCriterion()->setEnergyIndexed(false) ;
			elements.back()->refresh(&father);
			elements.back()->getState().initialize() ;
			K.add(elements.back());
			volume+= elements.back()->volume() ;
			
			for(size_t j = 0 ; j < cache3d[i]->neighbour.size() ; j++)
			{
				if(cache3d[i]->getNeighbour(j)->isTetrahedron() )
				{
					DelaunayTetrahedron * tri = dynamic_cast<DelaunayTetrahedron *>(cache3d[i]->getNeighbour(j)) ;
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
				if(dist(elements[i]->getCenter(), s.getParent()->getCenter()) < POINT_TOLERANCE)
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

void FractureCriterion::step(const ElementState &s)
{
	

		if(cache.empty() && cache3d.empty())
				initialiseCache(s) ;
		
	scoreAtState = grade(s) ;
	if(energyIndexed && s.getDeltaTime() > POINT_TOLERANCE )
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
	
		if(!getCache().empty())
		{
			Circle c(getMaterialCharacteristicRadius(), s.getParent()->getCenter()) ;
			for(size_t i = 0 ; i < getCache().size() ; i++)
			{
				if(c.in(*getCache()[i]->first) || c.in(*getCache()[i]->second) || c.in(*getCache()[i]->third))
				{
					currentEnergy += getCache()[i]->getState().elasticEnergy()*getCache()[i]->area() ;
					if(buildvolume)
						neighbourhoodvolume += getCache()[i]->area() ;
				}
			}
		}
		else if(!getCache3d().empty())
		{
			Sphere c(getMaterialCharacteristicRadius(), s.getParent()->getCenter()) ;
			
			for(size_t i = 0 ; i < getCache3d().size() ; i++)
			{
				if(c.in(*getCache3d()[i]->first) || c.in(*getCache3d()[i]->second) || c.in(*getCache3d()[i]->third) || c.in(*getCache3d()[i]->fourth))
				{
					currentEnergy += getCache3d()[i]->getState().elasticEnergy()*getCache3d()[i]->volume() ;
					if(buildvolume)
						neighbourhoodvolume += getCache3d()[i]->volume() ;
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
	if( s.getParent()->getBehaviour()->getDamageModel() == NULL )
		return false ;
	if( s.getParent()->getBehaviour()->getDamageModel() && s.getParent()->getBehaviour()->getDamageModel()->fractured())
		return false ;
	
	double tol = 1e-2 ;
	DelaunayTriangle * testedTri = dynamic_cast<DelaunayTriangle *>(s.getParent()) ;
	DelaunayTetrahedron * testedTet = dynamic_cast<DelaunayTetrahedron *>(s.getParent()) ;
	HexahedralElement * testedHex = dynamic_cast<HexahedralElement *>(s.getParent()) ;
	if(testedTri)
	{

		if(testedTri->visited)
			return false ;
		
		if (scoreAtState <= -tol)
			return false ;

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
				if(cache[i]->getBehaviour()->getFractureCriterion())
				{
					if(!cache[i]->getBehaviour()->fractured())
						areamax += area[i] ;
// 					if( !cache[i]->getBehaviour()->fractured())
// 					{
						double s = cache[i]->getBehaviour()->getFractureCriterion()->getSteppedScore() ;
						scores[-s] =  cache[i];
						unsortedScores.push_back(s);
						if(s > maxNeighbourhoodScore)
						{
							maxNeighbourhoodScore = s ;
							maxLocus = cache[i] ;
						}
// 					}
// 					else if(cache[i]->getBehaviour()->fractured())
// 					{
// 						double s = POINT_TOLERANCE ;
// 						scores[-s] =  cache[i];
// 						unsortedScores.push_back(s);
// 						if(s > maxNeighbourhoodScore)
// 						{
// 							maxNeighbourhoodScore = s ;
// 							maxLocus = cache[i] ;
// 						}
						
// 					}
					areatemp[cache[i]] = area[i] ;
				}
				
// 				if ((maxNeighbourhoodScore*tol) > score)
// 					return false ;
					
			}
		}
		
		if(!maxLocus)
			return false ;
		
// 		std::vector<DelaunayTriangle *> maxloci ;
		bool nearmaxlocus = false;
		
		for(size_t i = 0 ; i< cache.size() ; i++)
		{
			if(cache[i]->getBehaviour()->getFractureCriterion())
			{
				if(maxNeighbourhoodScore-cache[i]->getBehaviour()->getFractureCriterion()->getSteppedScore() < tol)
				{
// 					maxloci.push_back(cache[i]) ;
					if(squareDist2D(cache[i]->getCenter(), s.getParent()->getCenter()) < physicalCharacteristicRadius*physicalCharacteristicRadius)
					{
						nearmaxlocus = true ;
						break ;
					}
				}
			}
		}
		
		bool foundcutoff = false ;
		double thresholdscore = maxNeighbourhoodScore ;
		double avgscore = 0 ;
		double trialarea = std::min(physicalCharacteristicRadius*physicalCharacteristicRadius*M_PI, areamax) ;
		for(std::map<double, DelaunayTriangle *>::iterator i = scores.begin() ; i != scores.end() ; ++i)
		{
			if(!foundcutoff)
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
			}
// 			else
// 				i->second->visited = true ;
		}
		
		if (!foundcutoff && areamax > s.getParent()->area())
			return false ;
		if (nearmaxlocus)
			return true ;
		
		return false ;

	}
	if(testedTet)
	{

		if(testedTri->visited)
			return false ;
				
		if (scoreAtState <= 0)
			return false ;

		double maxNeighbourhoodScore = 0 ;
		double matchedArea = 0 ;
		std::map<double, DelaunayTetrahedron *> scores ;
		std::vector<double> unsortedScores ;
		std::map<DelaunayTetrahedron *, double> areatemp ;
		DelaunayTetrahedron * maxLocus = NULL;
		
		if(!cache.empty())
		{
			for(size_t i = 0 ; i< cache3d.size() ; i++)
			{
				if( cache3d[i]->getBehaviour()->getFractureCriterion())
				{
					if( !cache3d[i]->getBehaviour()->fractured())
					{
						double s = cache3d[i]->getBehaviour()->getFractureCriterion()->getSteppedScore() ;
						scores[-s] =  cache3d[i];
						unsortedScores.push_back(s);
						if(s > maxNeighbourhoodScore)
						{
							maxNeighbourhoodScore = s ;
							maxLocus = cache3d[i] ;
						}
					}
					else if(cache3d[i]->getBehaviour()->fractured())
					{
						double s = POINT_TOLERANCE ;
						scores[-s] =  cache3d[i];
						unsortedScores.push_back(s);
					}
					areatemp[cache3d[i]] = area[i] ;
				}
			}
		}
		
		if(!maxLocus)
			return false ;
		
		std::vector<DelaunayTetrahedron *> maxloci ;
		
		for(size_t i = 0 ; i< cache3d.size() ; i++)
		{
			if(cache3d[i]->getBehaviour()->getFractureCriterion())
				if(std::abs(cache3d[i]->getBehaviour()->getFractureCriterion()->getSteppedScore()-maxNeighbourhoodScore) < tol)
					maxloci.push_back(cache3d[i]) ;
		}
		
		bool foundcutoff = false ;
		double thresholdscore = maxNeighbourhoodScore ;
		
		for(std::map<double, DelaunayTetrahedron *>::iterator i = scores.begin() ; i != scores.end() ; ++i)
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
				if (matchedArea > physicalCharacteristicRadius*physicalCharacteristicRadius*M_PI)
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
			return false ;

		for(size_t i = 0 ; i < maxloci.size() ; i++)
			if(squareDist3D(maxloci[i]->getCenter(), s.getParent()->getCenter()) < physicalCharacteristicRadius*physicalCharacteristicRadius)
				return true ;
		
		return false ;
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
			for(std::set<HexahedralElement *>::iterator i= neighbourhood.begin() ; i != neighbourhood.end() ; ++i)
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
				return true ;
			}
		}

		return false ;
	}
	else
	{
		std::cout << " criterion not implemented for this kind of element" << std::endl ;
		return false ;
	}
	
	//shut up the compiler
	return false ;
}

Material FractureCriterion::toMaterial()
{
	Material mat ;
	return mat ;
}


