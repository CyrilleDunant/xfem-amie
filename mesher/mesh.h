
#ifndef MESH_H
#define MESH_H

#include <vector>
#include <limits>
#include "../elements/integrable_entity.h"

namespace Mu
{
	template <class ETYPE, class EABSTRACTTYPE>
	class Mesh
	{
		protected:
			std::map<Mesh<ETYPE, EABSTRACTTYPE> *, std::map<ETYPE *, std::vector<ETYPE *> > > cache ;
		public:
			virtual std::vector<EABSTRACTTYPE *> & getTree() = 0;
			virtual const std::vector<EABSTRACTTYPE *> & getTree() const = 0 ;
		
		public:
			Mesh() {} ;
			virtual ~Mesh() {} ;
			virtual std::vector<ETYPE *> getElements() = 0;
			virtual std::vector<ETYPE *> getConflictingElements(const Point  * p) = 0;
			virtual std::vector<ETYPE *> getConflictingElements(const Geometry * g) = 0;
			virtual void setElementOrder(Order o) = 0;
			virtual void insert(Point *) = 0 ;
			template <class ETARGETTYPE>
			/** \brief Return the displacements in source mesh projected on current mesh. 
			*/
			void project(Mesh<ETARGETTYPE, EABSTRACTTYPE> * mesh, Vector & projection)
			{
				if(cache.find(mesh) == cache.end())
				{
					std::map<ETYPE *, std::vector<ETYPE *> > projectionCache ;
					std::vector<ETYPE *> selfElements = getElements() ;
					int pointCount = 0 ;
					size_t idCount = 0 ;
					size_t numDofs = 0 ;
					for(size_t i = 0 ; i < selfElements.size() ; i++)
					{
						if(selfElements[i]->getBehaviour()->type != VOID_BEHAVIOUR)
						{
							numDofs = std::max(selfElements[i]->getBehaviour()->getNumberOfDegreesOfFreedom(),numDofs) ;
							if(numDofs)
								break ;
						}
					}
					
					projection.resize(numDofs*getLastNodeId()) ;
					projection = 0 ;
					
					for(size_t i = 0 ; i < selfElements.size() ; i++)
					{
						if(selfElements[i]->getBehaviour()->type != VOID_BEHAVIOUR)
						{
							std::vector<ETARGETTYPE *> targets = mesh->getConflictingElements(selfElements[i]->getPrimitive()) ;
							if(targets.empty())
							{
								Circle c(selfElements[i]->getRadius()*1.1, selfElements[i]->getCenter()) ;
								targets =  mesh->getConflictingElements(&c) ;
								if(targets.empty())
								{
									std::vector<ETARGETTYPE *> targets = mesh->getConflictingElements(&selfElements[i]->getCenter()) ;
								}
							}
							
							projectionCache[selfElements[i]] = targets ;
							for(size_t j = 0 ; j < selfElements[i]->getBoundingPoints().size() ; j++)
							{
								std::vector<ETARGETTYPE *> coincidentElements ;

								for(size_t k = 0 ; k < targets.size() ; k++)
								{
									for(size_t l = 0 ; l < targets[k]->getBoundingPoints().size() ; l++)
									{
										Point proj(selfElements[i]->getBoundingPoint(j)) ; 
										targets[k]->project(&proj) ;

										if((targets[k]->in(selfElements[i]->getBoundingPoint(j))
											|| squareDist2D(proj, selfElements[i]->getBoundingPoint(j)) < POINT_TOLERANCE) 
											&& targets[k]->getBehaviour()->type != VOID_BEHAVIOUR)
											coincidentElements.push_back(targets[k]) ;
									}
								}

								Vector disps(numDofs) ;
								for(size_t k = 0 ; k < coincidentElements.size() ; k++)
								{
									disps = coincidentElements[k]->getState().getDisplacements(selfElements[i]->getBoundingPoint(j)) ;
								}
								for(size_t k = 0 ; k < numDofs ; k++)
									projection[selfElements[i]->getBoundingPoint(j).id*numDofs+k] = disps[k] ;
								
								
							}
						}
					}
					
					
					cache[mesh] = projectionCache ;
				}
				else
				{
					std::vector<ETYPE *> selfElements = getElements() ;
					size_t numDofs = 0 ;
					for(size_t i = 0 ; i < selfElements.size() ; i++)
					{
						if(selfElements[i]->getBehaviour()->type != VOID_BEHAVIOUR)
						{
							numDofs = std::max(selfElements[i]->getBehaviour()->getNumberOfDegreesOfFreedom(),numDofs) ;
							if(numDofs)
								break ;
						}
					}
					
					projection = 0 ;
					Vector disps(numDofs) ;
					for(size_t i = 0 ; i < selfElements.size() ; i++)
					{
						if(selfElements[i]->getBehaviour()->type != VOID_BEHAVIOUR)
						{
							std::vector<ETARGETTYPE *> targets = cache[mesh][selfElements[i]] ;

							for(size_t j = 0 ; j < selfElements[i]->getBoundingPoints().size() ; j++)
							{
								if(std::abs(std::accumulate(&projection[selfElements[i]->getBoundingPoint(j).id*numDofs], 
									               &projection[selfElements[i]->getBoundingPoint(j).id*numDofs+numDofs], 0.)) < std::numeric_limits<double>::epsilon())
								{
									std::vector<Point *> coincidentPoints ;
									std::vector<ETARGETTYPE *> coincidentElements ;

									for(size_t k = 0 ; k < targets.size() ; k++)
									{
										for(size_t l = 0 ; l < targets[k]->getBoundingPoints().size() ; l++)
										{
											Point proj(selfElements[i]->getBoundingPoint(j)) ; 
											targets[k]->project(&proj) ;
											if((targets[k]->in(selfElements[i]->getBoundingPoint(j))
												|| squareDist2D(proj, selfElements[i]->getBoundingPoint(j)) < POINT_TOLERANCE) 
												&& targets[k]->getBehaviour()->type != VOID_BEHAVIOUR)
												coincidentElements.push_back(targets[k]) ;
										}
									}

									
									for(size_t k = 0 ; k < coincidentElements.size() ; k++)
									{
										disps = coincidentElements[k]->getState().getDisplacements(selfElements[i]->getBoundingPoint(j)) ;
									}
									for(size_t k = 0 ; k < numDofs ; k++)
										projection[selfElements[i]->getBoundingPoint(j).id*numDofs+k] = disps[k] ;
									
								}
							}
						}
					}
				}
			
			} ;
	
			virtual size_t & getLastNodeId() = 0;
			virtual const size_t & getLastNodeId() const = 0;
	} ;
} ;












#endif // MESH_H