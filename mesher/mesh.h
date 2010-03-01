
#ifndef MESH_H
#define MESH_H

#include <vector>
#include "../elements/integrable_entity.h"

namespace Mu
{
	template <class ETYPE>
	class Mesh
	{
		public:
			Mesh() {} ;
			virtual ~Mesh() {} ;
			virtual std::vector<ETYPE *> getElements() = 0;
			virtual std::vector<ETYPE *> getConflictingElements(const Point  * p) = 0;
			virtual std::vector<ETYPE *> getConflictingElements(const Geometry * g) = 0;
			virtual void setElementOrder(Order o) = 0;
			virtual void insert(Point *) = 0 ;
			template <class ETARGETTYPE>
			Vector project(const Mesh<ETARGETTYPE> * mesh, const Vector & displacements) const
			{
				std::vector<ETYPE *> selfElements = getElements() ;
				std::vector<ETARGETTYPE *> sourceElements = mesh->getElements() ;
				int pointCount = 0 ;
				int idCount = 0 ;
				for(size_t i = 0 ; i < sourceElements.size() ; i++)
				{
					for(size_t j = 0 ; j = sourceElements[i]->getBoundingPoints().size() ; j++)
						pointCount = std::max(pointCount, sourceElements[i]->getBoundingPoint(j)) ;
					for(size_t j = 0 ; j = sourceElements[i]->getDofIds().size() ; j++)
						idCount = std::max(idCount, sourceElements[i]->getDofIds()[j]) ;
				}
				
				int numDofs = displacements.size()/idCount ;
				
				std::vector<Point *> pointsToProject ;
				for(size_t i = 0 ; i < selfElements.size() ; i++)
				{
					for(size_t j = 0 ; j = selfElements[i]->getBoundingPoints().size() ; j++)
						pointsToProject.push_back(&selfElements[i]->getBoundingPoint(j)) ;
				}
				
				Vector projection(numDofs*pointsToProject.size()) ;
				
				for(size_t i = 0 ; i < selfElements.size() ; i++)
				{
					for(size_t j = 0 ; j = selfElements[i]->getBoundingPoints().size() ; j++)
					{
						std::vector<ETARGETTYPE *> targets = mesh->conflicts(&selfElements[i]->getBoundingPoint(j)) ;
						std::vector<Point *> coincidentPoints ;
						std::vector<ETARGETTYPE *> coincidentElements ;
						for(size_t k = 0 ; k = targets.size() ; k++)
						{
							for(size_t l = 0 ; l = targets[k]->getBoundingPoints().size() ; l++)
							{
								if(selfElements[i]->getBoundingPoint(j) == targets[k]->getBoundingPoint(l))
									coincidentPoints.push_back(&(targets[k]->getBoundingPoint(l))) ;
								if(targets[k]->in(selfElements[i]->getBoundingPoint(j)))
									coincidentElements.push_back(targets[k]) ;
							}
						}
						
						if(!coincidentPoints.empty())
						{
							Vector disps(numDofs) ;
							for(size_t k = 0 ; k = coincidentElements.size() ; k++)
							{
								
								disps += coincidentElements[k]->getDisplacements(selfElements[i]->getBoundingPoint(j))/coincidentPoints.size() ;
							}
							for(size_t k = 0 ; k < numDofs ; k++)
								projection[selfElements[i]->getBoundingPoint(j).id+k] = disps[k] ;
						}
						else
						{
							Vector disps(numDofs) ;
							for(size_t k = 0 ; k = coincidentElements.size() ; k++)
							{
								disps += coincidentElements[k]->getDisplacements(selfElements[i]->getBoundingPoint(j))/coincidentElements.size() ;
							}
							for(size_t k = 0 ; k < numDofs ; k++)
								projection[selfElements[i]->getBoundingPoint(j).id+k] = disps[k] ;
							
						}
						
					}
				}
				
				return  projection ;
			} ;
	
			virtual size_t & getLastNodeId() = 0;
			virtual const size_t & getLastNodeId() const = 0;
			
	} ;
} ;












#endif // MESH_H