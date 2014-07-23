/*
    mesh abstract implementation for AMIE
    Copyright (C) 2010  Cyrille Dunant

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

*/


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
// 			virtual std::vector<EABSTRACTTYPE *> & getTree() = 0;
// 			virtual const std::vector<EABSTRACTTYPE *> & getTree() const = 0 ;
			virtual size_t addToTree(EABSTRACTTYPE * toAdd) = 0 ;
			virtual EABSTRACTTYPE * getInTree(int index) = 0 ;
			virtual std::vector<Point * > & getAdditionalPoints() = 0 ;
			virtual const std::vector<Point * > & getAdditionalPoints() const = 0 ;
			virtual void extrude(double dt) = 0 ;
			virtual void extrude(const Vector & dt) = 0 ;
			virtual double getInternalScale() const { return 1. ;} ;
		public:
			Mesh() {} ;
			virtual ~Mesh() {} ;
			virtual std::vector<ETYPE *> getElements() = 0;
			virtual std::vector<ETYPE *> getConflictingElements(const Point  * p)  = 0;
			virtual std::vector<ETYPE *> getConflictingElements(const Geometry * g) = 0;
			
			virtual std::vector<ETYPE *> getNeighbouringElementsInGeometry(ETYPE * start , const Geometry * g) const
			{
				if(!start)
				{
// 					std::cout << "nullptr" << std::endl ;
					return std::vector<ETYPE *>() ;
				}

				std::set<ETYPE *> to_test ;
				std::set<ETYPE *> found ;
				found.insert(start) ;
				for(size_t i = 0 ; i < start->neighbourhood.size() ; i++)
				{
					if(start->getNeighbourhood(i)->timePlanes() > 1)
					{
						if(start->getNeighbourhood(i)->in(g->getCenter()) || g->intersects(start->getNeighbourhood(i)->getPrimitive()) )
						{
							to_test.insert(start->getNeighbourhood(i)) ;
						}
						
						for(size_t j = 0 ; j < start->getNeighbourhood(i)->timePlanes() ; j++)
						{
							Point c = start->getNeighbourhood(i)->getCenter() ;
							c.t = start->getNeighbourhood(i)->getBoundingPoint( start->getNeighbourhood(i)->getBoundingPoints().size() * j / start->getNeighbourhood(i)->timePlanes() ).t ;
							if(g->in(c)) 
							{
								to_test.insert(start->getNeighbourhood(i)) ;
							}

							if(g->getGeometryType() == TIME_DEPENDENT_CIRCLE)
							{	
								for(size_t k = 0 ; k <  start->getNeighbourhood(i)->getBoundingPoints().size() / start->getNeighbourhood(i)->timePlanes()-1 ; k++)
								{
									Point A = start->getNeighbourhood(i)->getBoundingPoint( start->getNeighbourhood(i)->getBoundingPoints().size() * j / start->getNeighbourhood(i)->timePlanes() + k ) ;
									Point B = start->getNeighbourhood(i)->getBoundingPoint( start->getNeighbourhood(i)->getBoundingPoints().size() * j / start->getNeighbourhood(i)->timePlanes() + k + 1) ;
									Segment s(A,B) ;
									if(s.intersects(g) || g->in(A) || g->in(B))
										to_test.insert(start->getNeighbourhood(i)) ;
								}
							}

						}
					}
					else if(g->in(start->getNeighbourhood(i)->getCenter()) || start->getNeighbourhood(i)->in(g->getCenter()) || g->intersects(start->getNeighbourhood(i)->getPrimitive()) )
					{
						to_test.insert(start->getNeighbourhood(i)) ;
					}
				}
				found.insert(to_test.begin(), to_test.end()) ;
				
				while(!to_test.empty())
				{
					std::set<ETYPE *> new_test ;
					for(auto j = to_test.begin() ; j != to_test.end() ; j++)
					{
						for(size_t i = 0 ; i < (*j)->neighbourhood.size() ; i++)
						{
							if(to_test.find((*j)->getNeighbourhood(i)) == to_test.end() 
								&& found.find((*j)->getNeighbourhood(i)) == found.end())
							{

								if((*j)->getNeighbourhood(i)->timePlanes() > 1)
								{
									if((*j)->getNeighbourhood(i)->in(g->getCenter()) || g->intersects((*j)->getNeighbourhood(i)->getPrimitive()) )
									{
										new_test.insert((*j)->getNeighbourhood(i)) ;
									}
									
									for(size_t k = 0 ; k < (*j)->getNeighbourhood(i)->getBoundingPoints().size()-1 ; k++)
									{
										Point A = (*j)->getNeighbourhood(i)->getBoundingPoint(k) ;
										Point B = (*j)->getNeighbourhood(i)->getBoundingPoint(k+1) ;
										Segment s(A,B) ;
										if(g->in(A) || g->in(B) || s.intersects(g))
										{
											new_test.insert((*j)->getNeighbourhood(i)) ;
											break ;
										}

									}
								}
								else if(g->in((*j)->getNeighbourhood(i)->getCenter()) || (*j)->getNeighbourhood(i)->in(g->getCenter()) || g->intersects((*j)->getNeighbourhood(i)->getPrimitive()) )
								{
									new_test.insert((*j)->getNeighbourhood(i)) ;
								}
							  
							}
						}
					}
					to_test = new_test ;
					found.insert(new_test.begin(), new_test.end()) ;
				}
				
				std::vector<ETYPE *> ret(found.begin(),found.end()) ;
				return ret ;
			}

			virtual ETYPE * getUniqueConflictingElement(const Point  * p) 
			{
				std::vector<ETYPE *> elements = getConflictingElements(p) ;
				for(size_t i = 0 ; i < elements.size() ; i++)
				{
					if(elements[i]->in(*p))
						return elements[i] ;
				}
				
				return nullptr ;
			}
			
			virtual void setElementOrder(Order o, double dt = 0.) = 0;
			virtual void insert(Point *) = 0 ;
			template <class ETARGETTYPE>
			/** \brief Return the displacements in source mesh projected on current mesh. 
			*/
			void project( Mesh<ETARGETTYPE, EABSTRACTTYPE> * mesh, Vector & projection, const Vector & source, bool fast = false)
			{
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
					
					for(auto i = points.begin() ; i != points.end() ; i++)
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
							auto e = std::unique(targets.begin(), targets.end()) ;
							targets.erase(e,targets.end()) ;
						}
						else
						{
							Circle c(rav*2., *(*i)) ;
							targets =  mesh->getConflictingElements(&c) ;
						}
						
						
						if(targets.empty())
						{
							std::cout << "failed projection, empty mesh" << std::endl ;
							(*i)->print() ;
							exit(0) ;
						}
						projectionPointCache[(*i)] = nullptr ;
						std::map<double, ETARGETTYPE *> coincidentElements;
						for(size_t k = 0 ; k < targets.size() ; k++)
						{
							Point proj(*(*i)) ;
							targets[k]->project(&proj) ;
							if(targets[k]->in(*(*i)) || dist(proj, *(*i)) < 128.*POINT_TOLERANCE_2D)
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
					
					for( auto i = cache[mesh].begin() ; i != cache[mesh].end() ; ++i)
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
					
					for( auto j = pointcache[mesh].begin() ; j != pointcache[mesh].end() ; ++j)
					{

						for(size_t k = 0 ; k < numDofs ; k++)
							projection[j->first->id*numDofs+k] = source[j->second->id*numDofs+k] ;
						
					}
					
				}
			
			} ;
	
			template <class ETARGETTYPE>
			/** \brief Return the displacements in source mesh projected on current mesh. 
			*/
			void leastSquareProject(const Mesh<ETARGETTYPE, EABSTRACTTYPE> * mesh, Vector & projection, const Vector & source, bool fast = false)
			{
				
			}
			virtual const size_t & getLastNodeId() const = 0;
	} ;
	
	template<class ETYPE, class EABSTRACTTYPE>
	class SingleElementMesh : public Mesh<ETYPE, EABSTRACTTYPE>
	{
	protected:
		ETYPE * element ;
		std::vector<EABSTRACTTYPE *> tree ;
		std::vector<Point *> points ;
		std::map<int *, int> trans ; 
		size_t global_counter ;
		
		void addSharedNodes(size_t nodes_per_side, size_t time_planes, double timestep)
		{
			for(auto  i = tree.begin() ; i != tree.end() ; ++i)
			{
				
				(*i)->visited = true ;
					
				size_t nodes_per_plane = nodes_per_side*3+3 ;
				
				std::valarray<Point *> newPoints(nodes_per_plane*time_planes) ;
				std::valarray<bool> done(false, nodes_per_plane*time_planes) ;
				
				for(size_t plane = 0 ; plane < time_planes ; plane++)
				{
					for(size_t side = 0 ; side < 3 ; side++)
					{
						Point a(static_cast<ETYPE *>(*i)->getBoundingPoint(side)) ;
						Point b(static_cast<ETYPE *>(*i)->getBoundingPoint((side+1)%3)) ;
						
						if(time_planes> 1)
						{
							a.t = (double)plane*(timestep/(double)(time_planes-1))-timestep/2.;
							b.t = (double)plane*(timestep/(double)(time_planes-1))-timestep/2.;
						}
						for(size_t node = 0 ; node < nodes_per_side+1 ; node++)
						{
							double fraction = (double)(node)/((double)nodes_per_side+1) ;
							Point proto = a*(1.-fraction) + b*fraction ;
							Point * foundPoint = nullptr ;
							
							for(size_t j = 0 ; j< static_cast<ETYPE *>(*i)->getBoundingPoints().size() ; j++)
							{
								if(static_cast<ETYPE *>(*i)->getBoundingPoint(j) == proto)
								{
									foundPoint = &static_cast<ETYPE *>(*i)->getBoundingPoint(j) ;
									break ;
								}
							}
							
							if(!foundPoint)
							{
								for(size_t j = 0 ; j < static_cast<ETYPE *>(*i)->neighbourhood.size() ; j++)
								{
									if(static_cast<ETYPE *>(*i)->getNeighbourhood(j)->visited)
									{
										ETYPE * n = static_cast<ETYPE *>(*i)->getNeighbourhood(j) ;
										for(size_t k = 0 ; k < n->getBoundingPoints().size();k++)
										{
											if(n->getBoundingPoint(k) == proto)
											{
												foundPoint = &n->getBoundingPoint(k) ;
												break ;
											}
										}
										
										if(foundPoint)
											break ;
									}
								}
							}
							
							if(!done[nodes_per_plane*plane+side*(nodes_per_side+1)+node])
							{
								if(foundPoint)
								{
									newPoints[nodes_per_plane*plane+side*(nodes_per_side+1)+node]  = foundPoint ;
								}
								else
								{
									points.push_back(new Point(proto) ) ;
									newPoints[nodes_per_plane*plane+side*(nodes_per_side+1)+node]  = points.back();
									newPoints[nodes_per_plane*plane+side*(nodes_per_side+1)+node]->id = global_counter++ ;
								}
								
								done[nodes_per_plane*plane+side*(nodes_per_side+1)+node] = true ;
							}
						}
					}
				}
				
				static_cast<ETYPE *>(*i)->setBoundingPoints(newPoints) ;
			}
					
			
			for(auto i = tree.begin() ; i != tree.end() ; ++i)
			{
				(*i)->clearVisited() ;
			}

		}

		
	public:
		
		int elementLayer(ETYPE * e) const { return -1 ;}
		void addElementToLayer(const ETYPE * element, int layer) { } ;
		
		SingleElementMesh(ETYPE * element) : element(element), global_counter(0)
		{
			tree.push_back(element) ;
			for(size_t i = 0 ; i < element->getBoundingPoints().size() ; ++i)
			{
				if(element->getBoundingPoint(i).id < 0)
					element->getBoundingPoint(i).id = global_counter ;
				trans[& (element->getBoundingPoint(i).id)] = global_counter ;
				global_counter++ ;
			}
		} ;
		
		virtual ~SingleElementMesh() 
		{
			for(int i = 0 ; i < points.size() ; i++)
				delete points[i] ;
		} ;
		virtual std::vector<ETYPE *> getElements() 
		{
			std::vector<ETYPE *> ret ; 
			ret.push_back(element) ; 
			return ret ; 
		}
		virtual std::vector<ETYPE *> getConflictingElements(const Point  * p) 
		{
			std::vector<ETYPE *> ret ; 
			if(element->in(*p))	
			{
				ret.push_back(element) ; 
			}
			return ret ; 
		}
		
		virtual void extrude(double dt) { std::cout << "should extrude.." << std::endl  ;} ;
		virtual void extrude(const Vector & dt) { std::cout << "should extrude.." << std::endl  ;} ;
		
		virtual std::vector<Point * > & getAdditionalPoints() {return points ; };
		virtual const std::vector<Point * > & getAdditionalPoints() const { return points ;};
		
		virtual std::vector<ETYPE *> getConflictingElements(const Geometry * g) 
		{
			std::vector<ETYPE *> ret ; 
			if(element->intersects(g) || g->in(element->getCenter()) || element->in(g->getCenter()))
				ret.push_back(element) ; 
			
			return ret ; 
		}
		/** Does nothing as this is a special-purpose mesh*/
		virtual void setElementOrder(Order o, double dt = 0.)
		{
			switch(o)
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
					addSharedNodes(1,1,0) ;
					break ;
				}
			case CUBIC:
				{
					addSharedNodes(2,1,0) ;
					break ;
				}
			case QUADRIC:
				{
					addSharedNodes(3,1,0) ;
					break ;
				}
			case QUINTIC:
				{
					addSharedNodes(3,1,0) ;
					break ;
				}
			case CONSTANT_TIME_LINEAR:
				{
					addSharedNodes(0,2,dt) ;
					break ;
				}
			case CONSTANT_TIME_QUADRATIC:
				{
					addSharedNodes(0,3,dt) ;
					break ;
				}
			case LINEAR_TIME_LINEAR:
				{
					addSharedNodes(0,2,dt) ;
					break ;
				}
			case LINEAR_TIME_QUADRATIC:
				{
					addSharedNodes(0,3,dt) ;
					break ;
				}
			case QUADRATIC_TIME_LINEAR:
				{
					addSharedNodes(1,2,dt) ;
					break ;
				}
			case QUADRATIC_TIME_QUADRATIC:
				{
					addSharedNodes(1,3,dt) ;
					break ;
				}
			case CUBIC_TIME_LINEAR:
				{
					addSharedNodes(2,2,dt) ;
					break ;
				}
			case CUBIC_TIME_QUADRATIC:
				{
					addSharedNodes(2,3,dt) ;
					break ;
				}
			case QUADRIC_TIME_LINEAR:
				{
					addSharedNodes(3,2,dt) ;
					break ;
				}
			case QUADRIC_TIME_QUADRATIC:
				{
					addSharedNodes(3,3,dt) ;
					break ;
				}
			case QUINTIC_TIME_LINEAR:
				{
					addSharedNodes(3,2,dt) ;
					break ;
				}
			case QUINTIC_TIME_QUADRATIC:
				{
					addSharedNodes(3,3,dt) ;
					break ;
				}
			default:
				break ;
				
			}
		}
		
		/** Does nothing as this is a special-purpose mesh*/
		virtual void insert(Point *)
		{
		}
		
		virtual const size_t & getLastNodeId() const { return global_counter ; }
		virtual size_t addToTree(EABSTRACTTYPE * toAdd)
		{
			tree.push_back(toAdd) ;
			return tree.size() -1 ;
		}
		
		virtual EABSTRACTTYPE * getInTree(int index) 
		{
			return tree[std::abs(index)] ;
		}
		
// 		virtual std::vector<EABSTRACTTYPE *> & getTree() {return tree ; }
// 		virtual const std::vector<EABSTRACTTYPE *> & getTree() const {return tree ; }
	} ;
} ;












#endif // MESH_H
