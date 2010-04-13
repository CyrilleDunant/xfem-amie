
#ifndef MESH_H
#define MESH_H

#include <vector>
#include <limits>
#include <set>
#include "../elements/integrable_entity.h"

namespace Mu
{
	template <class ETYPE, class EABSTRACTTYPE>
	class Mesh
	{
		protected:
			std::map<const Mesh<ETYPE, EABSTRACTTYPE> *, std::map<Point *, std::pair<ETYPE *, std::vector<double> > > > cache ;
			std::map<const Mesh<ETYPE, EABSTRACTTYPE> *, std::map<Point *, Point * > > pointcache ;
		public:
			virtual std::vector<EABSTRACTTYPE *> & getTree() = 0;
			virtual const std::vector<EABSTRACTTYPE *> & getTree() const = 0 ;
		
		public:
			Mesh() {} ;
			virtual ~Mesh() {} ;
			virtual std::vector<ETYPE *> getElements() = 0;
			virtual std::vector<ETYPE *> getConflictingElements(const Point  * p) const = 0;
			virtual std::vector<ETYPE *> getConflictingElements(const Geometry * g) const = 0;
			virtual void setElementOrder(Order o) = 0;
			virtual void insert(Point *) = 0 ;
			template <class ETARGETTYPE>
			/** \brief Return the displacements in source mesh projected on current mesh. 
			*/
			void project(const Mesh<ETARGETTYPE, EABSTRACTTYPE> * mesh, Vector & projection, const Vector & source, bool fast = false)
			{
				typedef typename std::map<Point *, std::pair<ETARGETTYPE *, std::vector<double> > >::const_iterator MapIterator ;
				typedef typename std::vector<ETARGETTYPE *>::const_iterator VecIterator ;
				if(cache.find(mesh) == cache.end())
				{
					std::map<Point *, std::pair<ETYPE *, std::vector<double> > > projectionCache ;
					std::map<Point *, Point * > projectionPointCache ;
					std::vector<ETYPE *> selfElements = getElements() ;
					int pointCount = 0 ;
					size_t idCount = 0 ;
					size_t numDofs = 0 ;
					double rav = 0 ;
					double ecount  = 0;
					std::set<Point *> points ;
					for(size_t i = 0 ; i < selfElements.size() ; i++)
					{
						if(selfElements[i]->getBehaviour()->type != VOID_BEHAVIOUR)
						{
							numDofs = std::max(selfElements[i]->getBehaviour()->getNumberOfDegreesOfFreedom(),numDofs) ;
							for(size_t j = 0 ; j < selfElements[i]->getBoundingPoints().size() ; j++)
								points.insert(&selfElements[i]->getBoundingPoint(j)) ;
							rav += selfElements[i]->getRadius() ;
							ecount++ ;
						}

					}
					rav /= ecount ;
					if(projection.size() == 0)
					{
						projection.resize(numDofs*getLastNodeId()) ;
						projection = 0 ;
					}
					
					for(std::set<Point *>::iterator i = points.begin() ; i != points.end() ; i++)
					{
						std::vector<ETARGETTYPE *> targets ; 

						std::vector<ETARGETTYPE *> temp = mesh->getConflictingElements(*i) ;
						for(size_t k = 0 ; k < temp.size() ; k++)
						{
							if(temp[k]->getBehaviour()->type != VOID_BEHAVIOUR)
								targets.push_back(temp[k]) ;
						}

						if(!targets.empty())
						{
							std::sort(targets.begin(), targets.end()) ;
							targets.erase(std::unique(targets.begin(), targets.end()),targets.end()) ;
						}
						else
						{
							Circle c(rav*2., *(*i)) ;
							targets =  mesh->getConflictingElements(&c) ;
						}
						
						
						if(targets.empty())
						{
							(*i)->print() ;
							exit(0) ;
						}
						projectionPointCache[(*i)] = NULL ;
						std::map<double, ETARGETTYPE *> coincidentElements;
						for(size_t k = 0 ; k < targets.size() ; k++)
						{
							Point proj(*(*i)) ;
							targets[k]->project(&proj) ;
							if(targets[k]->in(*(*i)) || dist(proj, *(*i)) < 128.*POINT_TOLERANCE)
								coincidentElements[dist(*(*i), targets[k]->getCenter())] = targets[k] ;
						}
						
						if(!coincidentElements.empty())
						{
							Vector disps(0., numDofs) ;
							projectionCache[(*i)] = std::make_pair(coincidentElements.begin()->second, coincidentElements.begin()->second->getState().getInterpolatingFactors(*(*i), false)) ;
							
							for(size_t j = 0 ; j < coincidentElements.begin()->second->getBoundingPoints().size(); j++)
							{
								if(coincidentElements.begin()->second->getBoundingPoint(j) == *(*i))
								{
									projectionPointCache[(*i)] = &coincidentElements.begin()->second->getBoundingPoint(j) ;
									break ;
								}
							}
							
							if(projectionPointCache[(*i)] && source.size())
							{
								projectionCache.erase(projectionCache.find(*i)) ;
								for(size_t k = 0 ; k < numDofs ; k++)
									disps[k] = source[projectionPointCache[(*i)]->id*numDofs+k] ;
							}
							else
							{
								projectionPointCache.erase(projectionPointCache.find(*i)) ;
								if(coincidentElements.begin()->second->getBehaviour()->type != VOID_BEHAVIOUR)
								{
									for(size_t j = 0 ; j < coincidentElements.begin()->second->getBoundingPoints().size() ; j++)
									{
										int id = coincidentElements.begin()->second->getBoundingPoint(j).id ;
										double d = projectionCache[(*i)].second[j] ;
										for(size_t k = 0 ; k < numDofs ; k++)
											disps[k] += source[id*numDofs+k]*d ;
									}
								}
							}
							
							for(size_t k = 0 ; k < numDofs ; k++)
								projection[(*i)->id*numDofs+k] = disps[k] ;
						}
						else
						{
							for(size_t k = 0 ; k < targets.size() ; k++)
							{
								Point proj(*(*i)) ;
								targets[k]->project(&proj) ;
								coincidentElements[dist(proj, *(*i))] = targets[k] ;
							}
							Vector disps(0., numDofs); 
							projectionCache[(*i)] = std::make_pair(coincidentElements.begin()->second, coincidentElements.begin()->second->getState().getInterpolatingFactors(*(*i), false)) ;
							
							for(size_t j = 0 ; j < coincidentElements.begin()->second->getBoundingPoints().size(); j++)
							{
								if(coincidentElements.begin()->second->getBoundingPoint(j) == *(*i))
								{
									projectionPointCache[(*i)] = &coincidentElements.begin()->second->getBoundingPoint(j) ;
									break ;
								}
							}
							
							if(projectionPointCache[(*i)] && source.size())
							{
								projectionCache.erase(projectionCache.find(*i)) ;
								int id = projectionPointCache[(*i)]->id ;
								for(size_t k = 0 ; k < numDofs ; k++)
									disps[k] = source[id*numDofs+k] ;
							}
							else
							{
								projectionPointCache.erase(projectionPointCache.find(*i)) ;
								if(coincidentElements.begin()->second->getBehaviour()->type != VOID_BEHAVIOUR)
								{
									for(size_t j = 0 ; j < coincidentElements.begin()->second->getBoundingPoints().size() ; j++)
									{
										double d = projectionCache[(*i)].second[j] ;
										int id = coincidentElements.begin()->second->getBoundingPoint(j).id ;
										for(size_t k = 0 ; k < numDofs ; k++)
											disps[k] += source[id*numDofs+k]*d ;
									}
								}
							}
							
							for(size_t k = 0 ; k < numDofs ; k++)
									projection[(*i)->id*numDofs+k] = disps[k] ;
						}
							
					}
					pointcache[mesh] = projectionPointCache ;
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
					Vector disps(0.,numDofs) ;
					
					for( MapIterator i = cache[mesh].begin() ; i != cache[mesh].end() ; ++i)
					{
						disps = 0 ;
						if(i->second.first->getBehaviour()->type != VOID_BEHAVIOUR)
						{
							for(size_t j = 0 ; j < i->second.first->getBoundingPoints().size() ; j++)
							{
								double d = i->second.second[j] ;
								int id = i->second.first->getBoundingPoint(j).id ;
								for(size_t k = 0 ; k < numDofs ; k++)
									disps[k] += source[id*numDofs+k]*d ;
							}
						}
						for(size_t k = 0 ; k < numDofs ; k++)
							projection[i->first->id*numDofs+k] = disps[k] ;
						
					}
					
					for(  std::map<Point*, Point*>::const_iterator j = pointcache[mesh].begin() ; j != pointcache[mesh].end() ; ++j)
					{

						for(size_t k = 0 ; k < numDofs ; k++)
							projection[j->first->id*numDofs+k] = source[j->second->id*numDofs+k] ;
						
					}
					
				}
			
			} ;
	
			virtual size_t & getLastNodeId() = 0;
			virtual const size_t & getLastNodeId() const = 0;
	} ;
} ;












#endif // MESH_H
