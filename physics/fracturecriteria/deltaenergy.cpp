#include "deltaenergy.h"
#include "../../solvers/assembly.h"
#include "../stiffness.h"

using namespace Mu ;

DeltaEnergy::DeltaEnergy(double radius, FractureCriterion * criterion) : radius(radius), criterion(criterion)
{
} ;

DeltaEnergy::~DeltaEnergy()  { } ;

double DeltaEnergy::grade(const ElementState &s)
{
	Assembly K ;
	if(s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
	{
		Circle c(radius, s.getParent()->getCenter()) ;
		if(cache2d.empty())
		{
			cache2d = s.getParent()->get2DMesh()->getConflictingElements(&c) ;
		}
		
		std::vector<DelaunayTriangle> elements ;
		std::vector<LinearForm *> behaviours ;
		TriElement father ;
		double originalenergy = 0 ;
		double area = 0 ;
		for(size_t i = 0 ; i < cache2d.size() ; i++)
		{
			if(cache2d[i]->getBehaviour()->getFractureCriterion())
			{
				DeltaEnergy * de = dynamic_cast<DeltaEnergy *>(cache2d[i]->getBehaviour()->getFractureCriterion()) ;
				if(de)
				{
					de->criterion->step(s);
					originalenergy += de->criterion->getScoreAtState()*cache2d[i]->area() ;
				}
				else
				{
					cache2d[i]->getBehaviour()->getFractureCriterion()->step(s);
					originalenergy += cache2d[i]->getBehaviour()->getFractureCriterion()->getScoreAtState()*cache2d[i]->area() ;
				}
			}
			
			elements.push_back(DelaunayTriangle(NULL, NULL,   cache2d[i]->first,  cache2d[i]->second,   cache2d[i]->third,  NULL) );
			elements.back().setBehaviour(cache2d[i]->getBehaviour()->getCopy()) ;
			if(cache2d[i] == s.getParent())
				elements.back().getBehaviour()->setTensor(elements.back().getBehaviour()->getTensor(elements.back().getCenter())*0.5) ;
			
			elements.back().refresh(&father);
			elements.back().getState().initialize() ;
			K.add(&elements.back());
			area += elements.back().area() ;
			for(size_t j = 0 ; j < cache2d[i]->neighbour.size() ; j++)
			{
				if(cache2d[i]->getNeighbour(j)->isTriangle )
				{
					DelaunayTriangle * tri = dynamic_cast<DelaunayTriangle *>(cache2d[i]->getNeighbour(j)) ;
					if(tri->getBehaviour()->type != VOID_BEHAVIOUR)
					{
						for(size_t k = 0 ; k <  tri->getBoundingPoints().size() ; k++)
						{
							if(!c.in(tri->getBoundingPoint(k)))
							{
								int id = tri->getBoundingPoint(k).id ;
								double ex = tri->getState().getDisplacements()[k*2];
								double ey = tri->getState().getDisplacements()[k*2+1];
								K.setPoint(ex, ey ,id);
							}
						}
					}
				}
			}
		
			
		}
	
		K.cgsolve() ;
		for(size_t i = 0 ; i < elements.size() ; i++)
			elements[i].step(0., &K.getDisplacements()) ;
		
		double energy = 0 ;
		for(size_t i = 0 ; i < elements.size() ;i++)
		{
			if(elements[i].getBehaviour()->getFractureCriterion())
			{
				DeltaEnergy * de = dynamic_cast<DeltaEnergy *>(elements[i].getBehaviour()->getFractureCriterion()) ;
				if(de)
				{
					de->criterion->step(elements[i].getState()) ;
					energy += de->criterion->getScoreAtState()*elements[i].area() ;
				}
				else
				{
					elements[i].getBehaviour()->getFractureCriterion()->step(elements[i].getState()) ;
					energy += elements[i].getBehaviour()->getFractureCriterion()->getScoreAtState()*elements[i].area() ;
				}
			}
		}
		
		return (originalenergy-energy)/area ;
	}
	else
	{
		Sphere c(radius, s.getParent()->getCenter()) ;
		if(cache3d.empty())
		{
			cache3d = s.getParent()->get3DMesh()->getConflictingElements(&c) ;
		}
		
		std::vector<DelaunayTetrahedron *> elements ;
		std::vector<LinearForm *> behaviours ;
		TetrahedralElement father ;
		double originalenergy = 0 ;
		double volume = 0 ;
		for(size_t i = 0 ; i < cache3d.size() ; i++)
		{
			if(cache3d[i]->getBehaviour()->getFractureCriterion())
			{
				DeltaEnergy * de = dynamic_cast<DeltaEnergy *>(cache3d[i]->getBehaviour()->getFractureCriterion()) ;
				if(de)
				{
					de->criterion->step(s);
					originalenergy += de->criterion->getScoreAtState()*cache3d[i]->area() ;
				}
				else
				{
					cache3d[i]->getBehaviour()->getFractureCriterion()->step(s);
					originalenergy += cache3d[i]->getBehaviour()->getFractureCriterion()->getScoreAtState()*cache2d[i]->area() ;
				}
			}
			
			elements.push_back(new DelaunayTetrahedron(NULL, NULL,   cache3d[i]->first,  cache3d[i]->second,   cache3d[i]->third, cache3d[i]->fourth,  NULL) );
			elements.back()->setBehaviour(cache2d[i]->getBehaviour()->getCopy()) ;
			if(cache3d[i] == s.getParent())
				elements.back()->getBehaviour()->setTensor(elements.back()->getBehaviour()->getTensor(elements.back()->getCenter())*0.5) ;
			
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
				DeltaEnergy * de = dynamic_cast<DeltaEnergy *>(elements[i]->getBehaviour()->getFractureCriterion()) ;
				if(de)
				{
					de->criterion->step(elements[i]->getState()) ;
					energy += de->criterion->getScoreAtState()*elements[i]->volume() ;
				}
				else
				{
					elements[i]->getBehaviour()->getFractureCriterion()->step(elements[i]->getState()) ;
					energy += elements[i]->getBehaviour()->getFractureCriterion()->getScoreAtState()*elements[i]->volume() ;
				}
			}
		}
		
		return (originalenergy-energy)/volume ;
	}

	
} ;

void DeltaEnergy::step(const ElementState &s)
{
	criterion->step(s);
	scoreAtState = grade(s) ;
}

FractureCriterion * DeltaEnergy::getCopy() const
{
	return new DeltaEnergy(radius, criterion->getCopy()) ;
}

Material DeltaEnergy::toMaterial()
{
	Material mat ;
	return mat ;
}