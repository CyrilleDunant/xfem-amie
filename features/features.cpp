
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "features.h"
#include "crackinitiation.h"
#include "layeredinclusion.h"
#include "sample.h"
#include "sample3d.h"
#include "../physics/void_form.h"
#include "../physics/fracturecriteria/fracturecriterion.h"
#ifdef HAVE_OPENMP
#include <omp.h>
#endif
#include "../physics/homogeneised_behaviour.h"
#include "../solvers/multigrid.h"
#include "../solvers/multigridstep.h"



using namespace Mu ;



Mesh<DelaunayTriangle, DelaunayTreeItem> * FeatureTree::get2DMesh(int g)
{
	if(this->dtree == NULL && this->dtree3D == NULL)
		this->generateElements() ;
	if(g == -1)
		return this->dtree ;
	else
		return coarseTrees[g] ;
}

Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * FeatureTree::get3DMesh(int g)
{
	if(this->dtree3D == NULL && this->dtree == NULL)
		this->generateElements() ;
	if(g == -1)
		return this->dtree3D ;
	else
		return coarseTrees3D[g] ;
}

std::vector<DelaunayTriangle *> FeatureTree::getBoundingTriangles(Feature * f )
{
	if(f == NULL)
		return tree[0]->getBoundingElements(dtree) ;
	else
		return f->getBoundingElements(dtree) ;
}

FeatureTree::FeatureTree(Feature *first) : grid(NULL), grid3d(NULL)
{
	reuseDisplacements = false ; 
	useMultigrid = false ;
	this->dtree = NULL ;
	this->dtree3D = NULL ;
	if(first)
		this->addFeature(NULL, first) ;

	if(is2D())
		grid = new Grid((first->getBoundingBox()[1].x-first->getBoundingBox()[0].x)*1.1,
		                (first->getBoundingBox()[1].y-first->getBoundingBox()[2].y)*1.1, 30,
		                Point((first->getBoundingBox()[1].x+first->getBoundingBox()[0].x)*.5, 
		                      (first->getBoundingBox()[1].y+first->getBoundingBox()[2].y)*.5
		                     )) ;
 
	if(is3D())
		grid3d =new Grid3D((first->getBoundingBox()[7].x-first->getBoundingBox()[0].x)*1.1,
		                   (first->getBoundingBox()[7].y-first->getBoundingBox()[0].y)*1.1,
		                   (first->getBoundingBox()[7].z-first->getBoundingBox()[0].z)*1.1, 10, (first->getBoundingBox()[7]+first->getBoundingBox()[0])*.5);
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
	residualError = 10000 ;

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
	for(size_t i = 0 ; i < coarseTrees.size() ;i++)
		delete coarseTrees[i] ;
	for(size_t i = 0 ; i < coarseTrees3D.size() ;i++)
		delete coarseTrees3D[i] ;
	for(size_t i = 0 ; i < coarseAssemblies.size() ;i++)
		delete coarseAssemblies[i] ;
	
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
		std::vector<DelaunayTriangle *> triangles = this->dtree->getElements() ;
		size_t count = 0 ;
		std::cerr << " renumbering... " << std::flush ;

		for(std::vector<DelaunayTriangle *>::iterator i = triangles.begin() ; i != triangles.end() ; ++i)
		{
			for(size_t j = 0 ; j < (*i)->getBoundingPoints().size() ; j++)
			{
				(*i)->getBoundingPoint(j).id = -1 ;
			}
		}

		
		Grid tmpgrid =this->grid->getGrid(std::max((size_t)round(sqrt(triangles.size()))/1024, (size_t)1)) ; //magic number such that the cache is full, but not too much
		for(std::vector<DelaunayTriangle *>::iterator i = triangles.begin() ; i != triangles.end() ; ++i)
			tmpgrid.forceAdd((*i)->getPrimitive()) ;
		
		std::vector<Geometry *> sortedElements ;
		std::set<Geometry *> placedElements ;
		for(size_t i = 0 ; i < tmpgrid.pixels.size() ; ++i)
		{
			for(size_t j = 0 ; j < tmpgrid.pixels[i].size() ; ++j)
			{
				for(size_t k = 0 ; k < tmpgrid.pixels[i][j]->getFeatures().size() ; ++k)
				{
					if(placedElements.find( tmpgrid.pixels[i][j]->getFeatures()[k]) == placedElements.end()) 
					{
						placedElements.insert(tmpgrid.pixels[i][j]->getFeatures()[k]) ;
						sortedElements.push_back(tmpgrid.pixels[i][j]->getFeatures()[k]) ;
					}
				}
				
			}
		}
		
		for(std::vector<Geometry *>::iterator i = sortedElements.begin() ; i != sortedElements.end() ; ++i)
		{
			DelaunayTriangle * tri = dynamic_cast<DelaunayTriangle *>(*i) ;
			if(tri && tri->getBehaviour()->type != VOID_BEHAVIOUR)
			{
				for(size_t j = 0 ; j < tri->getBoundingPoints().size() ; j++)
				{
					if(tri->getBoundingPoint(j).id == -1)
						tri->getBoundingPoint(j).id = count++ ;
				}
			}
		}
		
		this->dtree->getLastNodeId() = count ;
		
		std::cerr << count*2 << " ...done " << std::endl ;

	}
	else if (is3D())
	{
		std::vector<DelaunayTetrahedron *> tets = this->dtree3D->getElements() ;
		size_t count = 0 ;
		std::cerr << " renumbering... " << std::flush ;

		for(std::vector<DelaunayTetrahedron *>::iterator i = tets.begin() ; i != tets.end() ; ++i)
		{
			for(size_t j = 0 ; j < (*i)->getBoundingPoints().size() ; j++)
			{
				(*i)->getBoundingPoint(j).id = -1 ;
			}
		}

		
		Grid3D tmpgrid =this->grid3d->getGrid(std::max(((size_t)round(pow(tets.size(), .333333)))/1024, (size_t)1)) ; //magic number such that the cache is full, but not too much
		for(std::vector<DelaunayTetrahedron *>::iterator i = tets.begin() ; i != tets.end() ; ++i)
		{
			tmpgrid.forceAdd((*i)->getPrimitive()) ;
		}
		
		std::vector<Geometry *> sortedElements ;
		std::set<Geometry *> placedElements ;
		for(size_t i = 0 ; i < tmpgrid.pixels.size() ; ++i)
		{
			
			for(size_t j = 0 ; j < tmpgrid.pixels[i].size() ; ++j)
			{
				for(size_t k = 0 ; k < tmpgrid.pixels[i][j].size() ; ++k)
				{
					for(size_t l = 0 ; l < tmpgrid.pixels[i][j][k]->getFeatures().size() ; ++l)
					{
						if(placedElements.find( tmpgrid.pixels[i][j][k]->getFeatures()[l]) == placedElements.end()) 
						{
							placedElements.insert(tmpgrid.pixels[i][j][k]->getFeatures()[l]) ;
							sortedElements.push_back(tmpgrid.pixels[i][j][k]->getFeatures()[l]) ;
						}
					}
				}
			}
		}
		
		
		for(std::vector<Geometry *>::iterator i = sortedElements.begin() ; i != sortedElements.end() ; ++i)
		{
			DelaunayTetrahedron * tet = dynamic_cast<DelaunayTetrahedron *>(*i) ;
			if(tet && tet->getBehaviour()->type != VOID_BEHAVIOUR)
			{
				for(size_t j = 0 ; j < tet->getBoundingPoints().size() ; j++)
				{
					if(tet->getBoundingPoint(j).id == -1)
						tet->getBoundingPoint(j).id = count++ ;
				}
			}
		}
		
		this->dtree3D->getLastNodeId() = count ;
		
		std::cerr << count*3 << " ...done " << std::endl ;


	}
	
	if(useMultigrid)
	{
		for(size_t g = 0 ; g < coarseTrees.size() ; g++)
		{
			std::vector<DelaunayTriangle *> triangles = coarseTrees[g]->getElements() ;
			size_t count = 0 ;
			std::cerr << " renumbering... " << std::flush ;

			for(std::vector<DelaunayTriangle *>::iterator i = triangles.begin() ; i != triangles.end() ; ++i)
			{
				for(size_t j = 0 ; j < (*i)->getBoundingPoints().size() ; j++)
				{
					(*i)->getBoundingPoint(j).id = -1 ;
				}
			}

			
			Grid tmpgrid =this->grid->getGrid(std::max((size_t)round(sqrt(triangles.size()))/1024, (size_t)1)) ; //magic number such that the cache is full, but not too much
			for(std::vector<DelaunayTriangle *>::iterator i = triangles.begin() ; i != triangles.end() ; ++i)
				tmpgrid.forceAdd((*i)->getPrimitive()) ;
			
			std::vector<Geometry *> sortedElements ;
			std::set<Geometry *> placedElements ;
			for(size_t i = 0 ; i < tmpgrid.pixels.size() ; ++i)
			{
				for(size_t j = 0 ; j < tmpgrid.pixels[i].size() ; ++j)
				{
					for(size_t k = 0 ; k < tmpgrid.pixels[i][j]->getFeatures().size() ; ++k)
					{
						if(placedElements.find( tmpgrid.pixels[i][j]->getFeatures()[k]) == placedElements.end()) 
						{
							placedElements.insert(tmpgrid.pixels[i][j]->getFeatures()[k]) ;
							sortedElements.push_back(tmpgrid.pixels[i][j]->getFeatures()[k]) ;
						}
					}
					
				}
			}
			
			for(std::vector<Geometry *>::iterator i = sortedElements.begin() ; i != sortedElements.end() ; ++i)
			{
				DelaunayTriangle * tri = dynamic_cast<DelaunayTriangle *>(*i) ;
				if(tri && tri->getBehaviour()->type != VOID_BEHAVIOUR)
				{
					for(size_t j = 0 ; j < tri->getBoundingPoints().size() ; j++)
					{
						if(tri->getBoundingPoint(j).id == -1)
							tri->getBoundingPoint(j).id = count++ ;
					}
				}
			}
			
			coarseTrees[g]->getLastNodeId() = count ;
			
			std::cerr << count*2 << " ...done " << std::endl ;
		}
		
		for(size_t g = 0 ; g < coarseTrees3D.size() ; g++)
		{
			std::vector<DelaunayTetrahedron *> tets = coarseTrees3D[g]->getElements() ;
			size_t count = 0 ;
			std::cerr << " renumbering... " << std::flush ;

			for(std::vector<DelaunayTetrahedron *>::iterator i = tets.begin() ; i != tets.end() ; ++i)
			{
				for(size_t j = 0 ; j < (*i)->getBoundingPoints().size() ; j++)
				{
					(*i)->getBoundingPoint(j).id = -1 ;
				}
			}

			
			Grid3D tmpgrid =this->grid3d->getGrid(std::max(((size_t)round(pow(tets.size(), .333333)))/1024, (size_t)1)) ; //magic number such that the cache is full, but not too much
			for(std::vector<DelaunayTetrahedron *>::iterator i = tets.begin() ; i != tets.end() ; ++i)
			{
				tmpgrid.forceAdd((*i)->getPrimitive()) ;
			}
			
			std::vector<Geometry *> sortedElements ;
			std::set<Geometry *> placedElements ;
			for(size_t i = 0 ; i < tmpgrid.pixels.size() ; ++i)
			{
				
				for(size_t j = 0 ; j < tmpgrid.pixels[i].size() ; ++j)
				{
					for(size_t k = 0 ; k < tmpgrid.pixels[i][j].size() ; ++k)
					{
						for(size_t l = 0 ; l < tmpgrid.pixels[i][j][k]->getFeatures().size() ; ++l)
						{
							if(placedElements.find( tmpgrid.pixels[i][j][k]->getFeatures()[l]) == placedElements.end()) 
							{
								placedElements.insert(tmpgrid.pixels[i][j][k]->getFeatures()[l]) ;
								sortedElements.push_back(tmpgrid.pixels[i][j][k]->getFeatures()[l]) ;
							}
						}
					}
				}
			}
			
			
			for(std::vector<Geometry *>::iterator i = sortedElements.begin() ; i != sortedElements.end() ; ++i)
			{
				DelaunayTetrahedron * tet = dynamic_cast<DelaunayTetrahedron *>(*i) ;
				if(tet && tet->getBehaviour()->type != VOID_BEHAVIOUR)
				{
					for(size_t j = 0 ; j < tet->getBoundingPoints().size() ; j++)
					{
						if(tet->getBoundingPoint(j).id == -1)
							tet->getBoundingPoint(j).id = count++ ;
					}
				}
			}
			
			coarseTrees3D[g]->getLastNodeId() = count ;
			
			std::cerr << count*3 << " ...done " << std::endl ;
		}
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
			dtree->setElementOrder(elemOrder) ;
			
			for(size_t j = 0 ; j < coarseTrees.size() ; j++)
				coarseTrees[j]->setElementOrder(elemOrder) ;
			
			stitched  = true ;
			
			Point a (0.2, 0.2) ;
			Point b (0.6, 0.2) ;
			Point c (0.2, 0.6) ;
			Point d (1./3., 1./3.) ;
			for(size_t j = 1 ; j < this->tree.size() ; j++)
			{
				if(!tree[j]->isEnrichmentFeature && tree[j]->getGeometryType() != TRIANGLE && tree[j]->getGeometryType() != RECTANGLE)
				{
					
					std::vector<DelaunayTriangle *> triangles = this->tree[j]->getElements(dtree) ;
	
					for(size_t i = 0 ; i < triangles.size() ; i++)
					{
						triangles[i]->refresh(father2D) ;
						if(triangles[i]->getPrimitive()->intersects(tree[j]))
						{
							
							Point proj_0(*triangles[i]->first) ;
							tree[j]->project(&proj_0) ;
							Point proj_1(*triangles[i]->second) ;
							tree[j]->project(&proj_1) ;
							Point proj_2(*triangles[i]->third) ;
							tree[j]->project(&proj_2) ;
							bool changed  = true;
							
							if(squareDist2D(&proj_0 , triangles[i]->first ) < POINT_TOLERANCE && 
								squareDist2D(&proj_1 , triangles[i]->second) < POINT_TOLERANCE && 
								squareDist2D(&proj_2 , triangles[i]->third) > 10.*POINT_TOLERANCE)
							{
								count+=changed ; 
								changed = false ;
								Point test = triangles[i]->getBoundingPoint(1) ;
								tree[j]->project(&test) ;
								if (inRoot(test))
								{
									Point orig(triangles[i]->getBoundingPoint(1)) ;
									tree[j]->project(&triangles[i]->getBoundingPoint(1)) ;
									if(triangles[i]->jacobianAtPoint(a) > 0 && 
									   triangles[i]->jacobianAtPoint(b) > 0 && 
									   triangles[i]->jacobianAtPoint(c) > 0 && 
									   triangles[i]->jacobianAtPoint(d) > 0
										)
									{
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
									else
									{
										triangles[i]->getBoundingPoint(1) = orig ;
									}
								}
		// 						std::cerr << "--> " << (*triangles)[i]->getBoundingPoint(1)->x << ", " << (*triangles)[i]->getBoundingPoint(1)->y << std::endl ;
							}
							if(squareDist2D(&proj_1 , triangles[i]->second) < POINT_TOLERANCE && 
								squareDist2D(&proj_2 , triangles[i]->third) < POINT_TOLERANCE && 
								squareDist2D(&proj_0 , triangles[i]->first) > 10.*POINT_TOLERANCE
								)
							{
								count+=changed ; 
								changed = false ;								
								Point test = triangles[i]->getBoundingPoint(3) ;
								tree[j]->project(&test) ;
								if (inRoot(test))
								{
									Point orig(triangles[i]->getBoundingPoint(3)) ;
									tree[j]->project(&triangles[i]->getBoundingPoint(3)) ;
									if(triangles[i]->jacobianAtPoint(a) > 0 && 
									   triangles[i]->jacobianAtPoint(b) > 0 && 
									   triangles[i]->jacobianAtPoint(c) > 0 && 
									   triangles[i]->jacobianAtPoint(d) > 0
										)
									{
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
									else
									{
										triangles[i]->getBoundingPoint(3) = orig ;
									}
								}
								
		// 						std::cerr << "--> " << (*triangles)[i]->getBoundingPoint(3)->x << ", " << (*triangles)[i]->getBoundingPoint(3)->y << std::endl ;
							}
							if(squareDist2D(&proj_2 , triangles[i]->third) < POINT_TOLERANCE && 
							   squareDist2D(&proj_0, triangles[i]->first) < POINT_TOLERANCE &&
							   squareDist2D(&proj_1, triangles[i]->second) > 10.*POINT_TOLERANCE
							   ) 
							{
								count+=changed ; 
								changed = false ;								
								Point test = triangles[i]->getBoundingPoint(5) ;
								tree[j]->project(&test) ;
								if (inRoot(test))
								{
									Point orig(triangles[i]->getBoundingPoint(5)) ;
									tree[j]->project(&triangles[i]->getBoundingPoint(5)) ;
									if(triangles[i]->jacobianAtPoint(a) > 0 && 
									   triangles[i]->jacobianAtPoint(b) > 0 && 
									   triangles[i]->jacobianAtPoint(c) > 0 && 
									   triangles[i]->jacobianAtPoint(d) > 0
										)
									{
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
									}
									else
									{
										triangles[i]->getBoundingPoint(5) = orig ;
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
			dtree3D->setElementOrder(elemOrder) ;
			
			stitched = true ;
			return ;
			std::vector<DelaunayTetrahedron *> tets = this->dtree3D->getElements() ;
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

		tree[0]->sample(n) ;

		for(size_t i  = 1 ; i < this->tree.size() ; i++)
		{
			double shape_factor =(sqrt(tree[0]->area())/(2.*M_PI*tree[0]->getRadius()))/(sqrt(tree[i]->area())/(2.*M_PI*tree[i]->getRadius()));
			if(shape_factor < POINT_TOLERANCE)
				continue ;
			size_t npoints = std::max((size_t)round(sqrt(tree[i]->area()/(total_area*shape_factor))*n), (size_t)8) ;

			if(npoints >= 8 && !tree[i]->isVirtualFeature && npoints < n)
				tree[i]->sample(npoints) ;
		}
	}
	else if (is3D())
	{
		std::cout << n << std::endl ;
		std::cerr << "\r 3D features... sampling feature 0/" << this->tree.size() << "          " << std::flush ;
		tree[0]->sample(n) ;

		double total_area = tree[0]->area()*tree[0]->area()/(4.*M_PI*tree[0]->getRadius()*tree[0]->getRadius()) ;
		int count = 0 ;
#pragma omp parallel for
		for(int i  = 1 ; i < (int)tree.size() ; i++)
		{
			std::cerr << "\r 3D features... sampling feature "<< count << "/" << this->tree.size() << "          " << std::flush ;
			
			double shape_factor = tree[i]->area()/(4.*M_PI*tree[i]->getRadius()*tree[i]->getRadius());
			size_t npoints = .15*(size_t)(((double)n*tree[i]->area()*shape_factor)/(total_area)) ;
			
			if(npoints > 14 && !tree[i]->isVirtualFeature)
				tree[i]->sample(npoints) ;

		}
		std::cerr << "\r 3D features... sampling feature "<< count << "/" << this->tree.size() << " ...done" << std::endl ;
	}
}

void FeatureTree::refine(size_t nit, SamplingCriterion *cri)
{
	for(size_t t = 0 ; t < 512 ; t++)
	{
		bool corrected = false ;
		std::vector <DelaunayTriangle *> triangles  =  dtree->getElements() ;
		
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
	double pointDensity = 0 ; 
	if(is2D())
		pointDensity = .5*sqrt(tree[0]->area())/meshPoints.size() ;
	else
		pointDensity = .5*pow(tree[0]->volume(), 1./3.)/meshPoints.size() ;
	
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
					std::vector<DelaunayTetrahedron *>  tet = this->dtree3D->getConflictingElements(zonesVec[i].first[j]) ;
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
								if(zonesVec[i].second->getFather()->getChild(l)->inBoundary(p0, pointDensity))
								{
									count_0++ ;
								}
								if(zonesVec[i].second->getFather()->getChild(l)->inBoundary(p1, pointDensity))
								{
									count_1++ ;
								}
								if(zonesVec[i].second->getFather()->getChild(l)->inBoundary(p1, pointDensity))
								{
									count_2++ ;
								}
								if(zonesVec[i].second->getFather()->getChild(l)->inBoundary(p2, pointDensity))
								{
									count_2++ ;
								}
								if(zonesVec[i].second->getFather()->getChild(l)->inBoundary(p3, pointDensity))
								{
									count_3++ ;
								}
								if(zonesVec[i].second->getFather()->getChild(l)->inBoundary(p4, pointDensity))
								{
									count_4++ ;
								}
							}
						}
						
						
	// 					for(size_t m = 0 ; m < enrichmentFeature.size() ; m++)
	// 					{
	// 						if(enrichmentFeature[m]->inBoundary(p0))
	// 							count_0++ ;
	// 						if(enrichmentFeature[m]->inBoundary(p1))
	// 							count_1++ ;
	// 						if(enrichmentFeature[m]->inBoundary(p2))
	// 							count_2++ ;
	// 						if(enrichmentFeature[m]->inBoundary(p3))
	// 							count_3++ ;
	// 						if(enrichmentFeature[m]->inBoundary(p4))
	// 							count_4++ ;
	// 						if(enrichmentFeature[m]->inBoundary(p5))
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
							if(enrichmentFeature[l]->inBoundary(*sample[k], pointDensity))
							{
								yes = false ;
								break ;
							}
						}
						for(size_t l = 0 ; l < zonesVec[i].second->getChildren().size() ; l++)
						{
							if(zonesVec[i].second->getChild(l)->inBoundary(*sample[k], pointDensity))
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

			std::vector<DelaunayTriangle *>  tri = this->dtree->getConflictingElements(zonesVec[i].first[j]) ;
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
							if(zonesVec[i].second->getFather()->getChild(l)->inBoundary(p0, pointDensity))
							{
								count_0++ ;
							}
							if(zonesVec[i].second->getFather()->getChild(l)->inBoundary(p1, pointDensity))
							{
								count_1++ ;
							}
							if(zonesVec[i].second->getFather()->getChild(l)->inBoundary(p1, pointDensity))
							{
								count_2++ ;
							}
						}
					}
				
			
			
				for(size_t m = 0 ; m < enrichmentFeature.size() ; m++)
				{
					if(enrichmentFeature[m]->inBoundary(p0, pointDensity))
						count_0++ ;
					if(enrichmentFeature[m]->inBoundary(p1, pointDensity))
						count_1++ ;
					if(enrichmentFeature[m]->inBoundary(p2, pointDensity))
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


Form * FeatureTree::getElementBehaviour(const DelaunayTriangle * t) const
{
	int root_box = 0 ;

	if(!inRoot(t->getCenter())) 
		return new VoidForm() ;

	if(t->getBoundingPoints().size()%3 != 0)
		return new VoidForm() ;
	
	for(size_t i = 0 ; i < t->getBoundingPoints().size() ; i++)
		if( t->getBoundingPoint(i).id == -1)
		{
			return new VoidForm() ;
		}
	
	std::vector<Geometry *> targetstmp = grid->coOccur(t->getPrimitive()) ;
	std::vector<Feature *> targets ;
	for(size_t i = 0 ; i < targetstmp.size() ; i++)
		targets.push_back(dynamic_cast<Feature *>(targetstmp[i]) ) ;
		
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

	
	if(tree[root_box]->getBehaviour(t->getCenter())->timeDependent())
	{
		if( !tree[root_box]->getBehaviour(t->getCenter())->spaceDependent())
			return tree[root_box]->getBehaviour(t->getCenter())->getCopy() ;
		else
		{
			Form * b = tree[root_box]->getBehaviour(t->getCenter())->getCopy() ;
			b->transform(t->getXTransform(), t->getYTransform()) ;
			return b ;
		}
	}
	else if(!tree[root_box]->getBehaviour(t->getCenter())->spaceDependent())
		return tree[root_box]->getBehaviour(t->getCenter())->getCopy() ;
	else
	{
		Form * b = tree[root_box]->getBehaviour(t->getCenter())->getCopy() ;
		b->transform(t->getXTransform(), t->getYTransform()) ;
		return b ;
	}
	
	return tree[root_box]->getBehaviour(t->getCenter())->getCopy() ;

}

Form * FeatureTree::getElementBehaviour(const DelaunayTetrahedron * t) const
{
	int root_box = 0 ;

	
	if(!inRoot(t->getCenter())) 
		return new VoidForm() ;
	
// 	for(size_t i = 0 ; i < t->getBoundingPoints().size() ; i++)
// 		if( t->getBoundingPoint(i).id == -1)
// 		{
// 			return new VoidForm() ;
// 		}
	
	std::vector<Geometry *> targetstmp = grid3d->coOccur(t->getCenter()) ;
	std::vector<Feature *> targets  ;
	for(size_t i = 0 ; i < targetstmp.size() ; ++i)
		targets.push_back(dynamic_cast<Feature *>(targetstmp[i])) ;
	
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
	
	if(tree[root_box]->getBehaviour(t->getCenter())->timeDependent())
	{
		if( !tree[root_box]->getBehaviour(t->getCenter())->spaceDependent())
			return tree[root_box]->getBehaviour(t->getCenter())->getCopy() ;
		else
		{
			Form * b = tree[root_box]->getBehaviour(t->getCenter())->getCopy() ;
			b->transform(t->getXTransform(), t->getYTransform(), t->getZTransform()) ;
			return b ;
		}
	}
	else if(!tree[root_box]->getBehaviour(t->getCenter())->spaceDependent())
		return tree[root_box]->getBehaviour(t->getCenter())->getCopy() ;
	else
	{
		Form * b = tree[root_box]->getBehaviour(t->getCenter())->getCopy() ;
		b->transform(t->getXTransform(), t->getYTransform(), t->getZTransform()) ;
		return b ;
	}
	
	return tree[root_box]->getBehaviour(t->getCenter())->getCopy() ;

}

Point * FeatureTree::checkElement( const DelaunayTetrahedron * t ) const
{
	double pointDensity = 0 ; 
	if(is2D())
		pointDensity = .6*sqrt(tree[0]->area())/meshPoints.size() ;
	else
		pointDensity = .6*pow(tree[0]->volume(), 1./3.)/meshPoints.size() ;
	
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
				
				count_in += tree[i]->inBoundary(*t->first, pointDensity) ;
				count_in += tree[i]->inBoundary(*t->second, pointDensity) ;
				count_in += tree[i]->inBoundary(*t->third, pointDensity) ;
				count_in += tree[i]->inBoundary(*t->fourth, pointDensity) ;
				
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
	double pointDensity = 0 ; 
	if(is2D())
		pointDensity = .6*sqrt(tree[0]->area())/meshPoints.size() ;
	else
		pointDensity = .6*pow(tree[0]->volume(), 1./3.)/meshPoints.size() ;
	
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
				
				count_in += tree[i]->inBoundary(*t->first, pointDensity) ;
				count_in += tree[i]->inBoundary(*t->second, pointDensity) ;
				count_in += tree[i]->inBoundary(*t->third, pointDensity) ;
				
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
	double n_void ;
	if(is2D())
	{
		std::vector<DelaunayTriangle *> triangles = this->dtree->getElements() ;
				
		std::cerr << " setting behaviours..." << std::flush ;
		for(size_t i = 0 ; i < triangles.size() ;i++)
		{
			if (i%1000 == 0)
				std::cerr << "\r setting behaviours... triangle " << i << "/" << triangles.size() << std::flush ;
			triangles[i]->refresh(father2D) ;
			if(!triangles[i]->getBehaviour())
				triangles[i]->setBehaviour(getElementBehaviour(triangles[i])) ;
			n_void++ ;
		}
		std::cerr << " ...done" << std::endl ;
		
		if(useMultigrid)
		{
			for(size_t i = 0 ; i < coarseTrees.size() ; i++)
			{
				
				triangles = coarseTrees[i]->getElements() ;
				for(size_t j = 0 ; j < triangles.size() ;j++)
				{
					if (j%1000 == 0)
						std::cerr << "\r setting behaviours... grid " << i << ", triangle " << j << "/" << triangles.size() << std::flush ;
					triangles[j]->refresh(father2D) ;
					std::vector<Geometry * > coocuring ;
					if(tree.size() > 1)
						coocuring = grid->coOccur(triangles[j]->getPrimitive()) ;
					if(coocuring.size() == 1 && !static_cast<Feature *>(coocuring[0])->getBehaviour(triangles[j]->getCenter())->spaceDependent())
					{
						if(coocuring[0]->in(*triangles[j]->first) && coocuring[0]->in(*triangles[j]->second) && coocuring[0]->in(*triangles[j]->third))
							triangles[j]->setBehaviour(static_cast<Feature *>(coocuring[0])->getBehaviour(triangles[j]->getCenter())->getCopy()) ;
						else
							triangles[j]->setBehaviour(new HomogeneisedBehaviour(dtree, triangles[j])) ;
					}
					else if (tree.size() == 1 && !tree[0]->getBehaviour(triangles[j]->getCenter())->spaceDependent())
						triangles[j]->setBehaviour(tree[0]->getBehaviour(triangles[j]->getCenter())->getCopy()) ;
					else
						triangles[j]->setBehaviour(new HomogeneisedBehaviour(dtree, triangles[j])) ;
				}
				std::cerr << " ...done" << std::endl ;
			}
		}
		
		setBehaviours = true ;
	}
	else
	{
		std::vector<DelaunayTetrahedron *> tetrahedrons = this->dtree3D->getElements() ;
		
		std::cerr << " setting behaviours..." << std::flush ;
		for(size_t i = 0 ; i < tetrahedrons.size() ;i++)
		{
			if (i%1000 == 0)
				std::cerr << "\r setting behaviours... tet " << i << "/" << tetrahedrons.size() << std::flush ;
			
			tetrahedrons[i]->refresh(father3D) ;
			if(!tetrahedrons[i]->getBehaviour())
				tetrahedrons[i]->setBehaviour(getElementBehaviour(tetrahedrons[i])) ;
			n_void++ ;
		}
		if(useMultigrid)
		{
			for(size_t i = 0 ; i < coarseTrees3D.size() ; i++)
			{
				
				tetrahedrons = coarseTrees3D[i]->getElements() ;
				for(size_t j = 0 ; j < tetrahedrons.size() ;j++)
				{
					tetrahedrons[j]->refresh(father3D) ;
					if (j%1000 == 0)
						std::cerr << "\r setting behaviours... grid " << i << ", triangle " << j << "/" << tetrahedrons.size() << std::flush ;
					tetrahedrons[j]->setBehaviour(new HomogeneisedBehaviour(dtree3D, tetrahedrons[j])) ;
				}
				std::cerr << " ...done" << std::endl ;
			}
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
	{
		this->generateElements() ;
	}	
	 
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
		initialized =true ;
	}
	
	
	std::cerr << " enriching..." << std::flush ;		
	
	if(enrichmentChanged())
	{
		if(this->dtree3D == NULL && this->dtree !=NULL)
		{
			initializeElements() ;
		}
		else
		{
			initializeElements() ;
		}
		
		for(size_t i = 1 ; i < this->tree.size() ; i++)
		{
			if(is3D())
			{
				if(this->tree[i]->isEnrichmentFeature && dynamic_cast<EnrichmentFeature *>(this->tree[i])->moved())
				{
					dynamic_cast<EnrichmentFeature *>(this->tree[i])->enrich(dtree3D->getLastNodeId(), dtree3D) ;
					if(useMultigrid)
					{
						for(size_t j =  0 ; j < coarseTrees.size() ; j++)
							dynamic_cast<EnrichmentFeature *>(this->tree[i])->enrich(coarseTrees[j]->getLastNodeId(), coarseTrees[j]) ;
					}
				}
			

			
				if(i%10 == 0)
					std::cerr << "\r enriching... feature " << i+1 <<"/" << this->tree.size() << std::flush ;
			}
			else
			{
				if(this->tree[i]->isEnrichmentFeature && dynamic_cast<EnrichmentFeature *>(this->tree[i])->moved())
				{
					dynamic_cast<EnrichmentFeature *>(this->tree[i])->enrich(dtree->getLastNodeId(), dtree) ;
					if(useMultigrid)
					{
						for(size_t j =  0 ; j < coarseTrees.size() ; j++)
							dynamic_cast<EnrichmentFeature *>(this->tree[i])->enrich(coarseTrees[j]->getLastNodeId(), coarseTrees[j]) ;
					}
				}
				if(i%10 == 0)
					std::cerr << "\r enriching... feature " << i+1 <<"/" << this->tree.size() << std::flush ;
			}
		}

		std::cerr << " ...done" << std::endl ;
	}
	
	if(dtree != NULL)
	{
		numdofs = dtree->getLastNodeId() ;
		triangles =dtree->getElements() ;
		
		for(size_t j = 0 ; j < triangles.size() ; j++)
		{
			if(	triangles[j]->getBehaviour()->type != VOID_BEHAVIOUR)
			{
				if(j%1000 == 0)
					std::cerr << "\r assembling stiffness matrix... triangle " << j+1 << "/" << triangles.size() << std::flush ;
				triangles[j]->refresh(father2D) ;
				K->add(triangles[j]) ;
			}
		}
		std::cerr << " ...done." << std::endl ;
		if(useMultigrid)
		{
			for(size_t i =  0 ; i < coarseTrees.size() ; i++)
			{
				triangles = coarseTrees[i]->getElements() ;
			
				for(size_t j = 0 ; j < triangles.size() ; j++)
				{
					if(	triangles[j]->getBehaviour()->type != VOID_BEHAVIOUR)
					{
						if(j%1000 == 0)
							std::cerr << "\r assembling stiffness matrix... grid " << i << " triangle " << j+1 << "/" << triangles.size() << std::flush ;
						triangles[j]->refresh(father2D) ;
						coarseAssemblies[i]->add(triangles[j]) ;
					}
				}
				std::cerr << " ...done." << std::endl ;
			}
		}
		
	}
	else
	{
		numdofs = dtree3D->getLastNodeId() ;
		std::vector<DelaunayTetrahedron *> tets = this->dtree3D->getElements() ;
		
		for(size_t j = 0 ; j < tets.size() ; j++)
		{
			if(	tets[j]->getBehaviour()->type != VOID_BEHAVIOUR)
			{
				
				if(j%1000 == 0)
					std::cerr << "\r assembling stiffness matrix... tetrahedron " << j+1 << "/" << tets.size() << std::flush ;
				
				tets[j]->refresh(father3D) ;
				K->add(tets[j]) ;
			}
		}
		
		std::cerr << " ...done." << std::endl ;
	}
}

std::vector<DelaunayTriangle> FeatureTree::getSnapshot2D() const
{
	std::vector<DelaunayTriangle> copy ;
	
	std::vector<DelaunayTriangle *> tris = dtree->getElements() ;
	
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
		std::vector<DelaunayTriangle *> elements = dtree->getElements() ;
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

const Vector & FeatureTree::getDisplacements(int g) const
{
	if(g == -1 || !useMultigrid)
		return K->getDisplacements() ;
	else if(coarseAssemblies.size() > g)
		return coarseAssemblies[g]->getDisplacements() ;
	else if(!coarseAssemblies.empty())
		return coarseAssemblies.back()->getDisplacements() ;
	else
		return  K->getDisplacements() ;
}


std::pair<Vector , Vector > FeatureTree::getStressAndStrain(int g)
{
	if(dtree != NULL)
	{
		std::vector<DelaunayTriangle *> elements = dtree->getElements() ;
		if(useMultigrid)
		{
			if(g != -1 && coarseTrees.size() > g)
				elements = coarseTrees[g]->getElements() ;
			else if(g != -1 && !coarseTrees.empty())
				elements = coarseTrees.back()->getElements() ;
		}
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
				if(i%1000 == 0)
					std::cerr << "\r computing strain+stress... element " << i+1 << "/" << elements.size() << std::flush ;
			}
		}
		std::cerr << " ...done." << std::endl ;
		return stress_strain ;
	}
	else
	{
		std::vector<DelaunayTetrahedron *> tets = dtree3D->getElements() ;
		if(g != -1)
			tets = coarseTrees3D[g]->getElements() ;
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
			if(i%1000 == 0)
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
			if(i%1000 == 0)
				std::cerr << "\r computing strain+stress... element " << i+1 << "/" << tets.size() << std::flush ;
		}
		std::cerr << " ...done." << std::endl ;
		return stress_strain ;
}

Vector FeatureTree::strainFromDisplacements() const
{
	if(dtree != NULL)
	{
		std::vector<DelaunayTriangle *> elements = dtree->getElements() ;
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
	double pointDensity = 0 ; 
	if(is2D())
		pointDensity = .6*sqrt(tree[0]->area())/meshPoints.size() ;
	else
		pointDensity = .6*pow(tree[0]->volume(), 1./3.)/meshPoints.size() ;
	
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
		if( mother->getChild(k)->inBoundary(*p, pointDensity))
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
		std::vector<DelaunayTriangle *> elements = dtree->getElements() ;
		
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
		std::vector<DelaunayTetrahedron *> elements = dtree3D->getElements() ;
		
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
		std::vector<DelaunayTriangle *> elements = dtree->getElements() ;
		
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
		
		std::vector<DelaunayTetrahedron *> elements = dtree3D->getElements() ;
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
	{
		this->K->clear() ;
		if(useMultigrid)
		{
			for(size_t j = 0 ; j < coarseAssemblies.size() ;j++)
			{
				coarseAssemblies[j]->clear() ;
			}
		}
	}

	assemble() ;
	
	if(dtree)
	{
		std::vector<DelaunayTriangle *> elements = dtree->getElements() ;
		for(size_t i = 0 ; i < elements.size() ; ++i)
		{
			elements[i]->applyBoundaryCondition(K) ;
		}
	}
	else
	{
		std::vector<DelaunayTetrahedron *> elements = dtree3D->getElements() ;
		
		for(size_t i = 0 ; i < elements.size() ; ++i)
		{
			elements[i]->applyBoundaryCondition(K) ;
		}
		
	}
	
	for(size_t i = 0 ; i < boundaryCondition.size() ; ++i)
	{
		if(dtree)
		{
			
			boundaryCondition[i]->apply(K, dtree) ;
			if(useMultigrid)
			{
				for(size_t j = 0 ; j < coarseTrees.size() ;j++)
				{
					boundaryCondition[i]->apply(coarseAssemblies[j], coarseTrees[j]) ;
				}
			}
		}
		if(dtree3D)
			boundaryCondition[i]->apply(K, dtree3D) ;
	}
	
	needAssembly = true ;
	meshChange = false ;

	if(solverConvergence || reuseDisplacements)
	{
		if(useMultigrid && !coarseAssemblies.empty())
		{
			if(is2D())
			{
				coarseTrees[0]->project(dtree, coarseAssemblies[0]->getDisplacements(), lastx) ;
				coarseAssemblies[0]->cgsolve(coarseAssemblies[0]->getDisplacements()) ;
			}
			else
			{
				coarseTrees3D[0]->project(dtree3D, coarseAssemblies[0]->getDisplacements(), lastx) ;
				coarseAssemblies[0]->cgsolve(coarseAssemblies[0]->getDisplacements()) ;
			}
			
			if(is2D())
			{
				std::cerr << " stepping through elements... grid " << 0 << std::flush ;
				std::vector<DelaunayTriangle *> elements = coarseTrees[0]->getElements() ;
				for(size_t i = 0 ; i < elements.size() ;i++)
				{	
					if(i%1000 == 0)
						std::cerr << "\r stepping through  grid " << 0 <<", elements... " << i << "/" << elements.size() << std::flush ;
					elements[i]->step(dt, &coarseAssemblies[0]->getDisplacements()) ;
				}
				std::cerr << " ...done" << std::endl ;

			}
			else if(is3D())
			{
				std::cerr << " stepping through elements... grid " << 0 << std::flush ;
				std::vector<DelaunayTetrahedron *> elements = coarseTrees3D[0]->getElements() ;
				for(size_t i = 0 ; i < elements.size() ;i++)
				{	
					if(i%1000 == 0)
						std::cerr << "\r stepping through  grid " << 0 <<", elements... " << i << "/" << elements.size() << std::flush ;
					elements[i]->step(dt, &coarseAssemblies[0]->getDisplacements()) ;
				}
				std::cerr << " ...done" << std::endl ;

			}
			for(size_t j = 1 ; j < coarseAssemblies.size() ;j++)
			{
				coarseTrees[j]->project(coarseTrees[j-1], coarseAssemblies[j]->getDisplacements(), coarseAssemblies[j-1]->getDisplacements());
				coarseAssemblies[j]->cgsolve(coarseAssemblies[j]->getDisplacements()) ;
				if(is2D())
				{
					std::cerr << " stepping through elements... grid " << j << std::flush ;
					std::vector<DelaunayTriangle *> elements = coarseTrees[j]->getElements() ;
					for(size_t i = 0 ; i < elements.size() ;i++)
					{	
						if(i%1000 == 0)
							std::cerr << "\r stepping through  grid " << j <<", elements... " << i << "/" << elements.size() << std::flush ;
						elements[i]->step(dt, &coarseAssemblies[j]->getDisplacements()) ;
					}
					std::cerr << " ...done" << std::endl ;

				}
				else if(is3D())
				{
					std::cerr << " stepping through elements... grid " << j << std::flush ;
					std::vector<DelaunayTetrahedron *> elements = coarseTrees3D[j]->getElements() ;
					for(size_t i = 0 ; i < elements.size() ;i++)
					{	
						if(i%1000 == 0)
							std::cerr << "\r stepping through  grid " << j <<", elements... " << i << "/" << elements.size() << std::flush ;
						elements[i]->step(dt, &coarseAssemblies[j]->getDisplacements()) ;
					}
					std::cerr << " ...done" << std::endl ;

				}
			}
		}


		if(useMultigrid && !coarseAssemblies.empty())
		{
			dtree->project(coarseTrees.back(), K->getDisplacements(), coarseAssemblies.back()->getDisplacements()) ;
			solverConvergence = K->cgsolve(K->getDisplacements()) ;
		}
		else
		{
			solverConvergence = K->cgsolve(lastx) ;
		}
		
		Vector r = K->getMatrix()*K->getDisplacements()-K->getForces() ;
		double perror = residualError ;
		residualError = sqrt(parallel_inner_product(&r[0], &r[0], r.size())) ;
		if(perror > residualError || solverConvergence)
			reuseDisplacements = true;
	}
	else
	{
		lastx = 0 ;
		if(useMultigrid && !coarseAssemblies.empty())
		{
			coarseAssemblies[0]->cgsolve() ;

			if(is2D())
			{
				std::cerr << " stepping through elements... grid " << 0 << std::flush ;
				std::vector<DelaunayTriangle *> elements = coarseTrees[0]->getElements() ;
				for(size_t i = 0 ; i < elements.size() ;i++)
				{	
					if(i%1000 == 0)
						std::cerr << "\r stepping through  grid " << 0 <<", elements... " << i << "/" << elements.size() << std::flush ;
					elements[i]->step(dt, &coarseAssemblies[0]->getDisplacements()) ;
				}
				std::cerr << " ...done" << std::endl ;
			}
			else if(is3D())
			{
				std::cerr << " stepping through elements... grid " << 0 << std::flush ;
				std::vector<DelaunayTetrahedron *> elements = coarseTrees3D[0]->getElements() ;
				for(size_t i = 0 ; i < elements.size() ;i++)
				{	
					if(i%1000 == 0)
						std::cerr << "\r stepping through  grid " << 0 <<", elements... " << i << "/" << elements.size() << std::flush ;
					elements[i]->step(dt, &coarseAssemblies[0]->getDisplacements()) ;
				}
				std::cerr << " ...done" << std::endl ;
			}
			
			for(size_t j = 1 ; j < coarseAssemblies.size() ;j++)
			{
				coarseTrees[j]->project(coarseTrees[j-1], coarseAssemblies[j]->getDisplacements(), coarseAssemblies[j-1]->getDisplacements()) ;
				coarseAssemblies[j]->cgsolve(coarseAssemblies[j]->getDisplacements()) ;
				if(is2D())
				{
					std::cerr << " stepping through elements... grid " << j << std::flush ;
					std::vector<DelaunayTriangle *> elements = coarseTrees[j]->getElements() ;
					for(size_t i = 0 ; i < elements.size() ;i++)
					{	
						if(i%1000 == 0)
							std::cerr << "\r stepping through  grid " << j <<", elements... " << i << "/" << elements.size() << std::flush ;
						elements[i]->step(dt, &coarseAssemblies[j]->getDisplacements()) ;
					}
					std::cerr << " ...done" << std::endl ;
				}
				else if(is3D())
				{
					std::cerr << " stepping through elements... grid " << j << std::flush ;
					std::vector<DelaunayTetrahedron *> elements = coarseTrees3D[j]->getElements() ;
					for(size_t i = 0 ; i < elements.size() ;i++)
					{	
						if(i%1000 == 0)
							std::cerr << "\r stepping through  grid " << j <<", elements... " << i << "/" << elements.size() << std::flush ;
						elements[i]->step(dt, &coarseAssemblies[j]->getDisplacements()) ;
					}
					std::cerr << " ...done" << std::endl ;
				}
			
			}
		}
	
		if(useMultigrid && !coarseAssemblies.empty())
		{
			
			K->mgprepare() ;
			std::vector<const CoordinateIndexedSparseMatrix *> coarseMatrices ;
			for(size_t j = 0 ; j < coarseAssemblies.size() ;j++)
				coarseMatrices.push_back(&coarseAssemblies[j]->getMatrix()) ;
 
			dtree->project(coarseTrees.back(), lastx, coarseAssemblies.back()->getDisplacements()) ;
// 			for(size_t j = 0 ; j < coarseAssemblies.size() ;j++)
// 			{
// 				ConjugateGradient cg(K->getMatrix(), K->getForces()) ;
// 				MultiGridStep<Mesh<DelaunayTriangle,DelaunayTreeItem>, DelaunayTriangle> mgs(dtree, 
// 																						 coarseTrees[j], 
// 																						 &K->getMatrix(), 
// 																						 coarseMatrices[j], coarseAssemblies[j]->getForces()) ;
// 				MultiGrid<Mesh<DelaunayTriangle,DelaunayTreeItem>, DelaunayTriangle> mg(K->getMatrix(), coarseMatrices, dtree, coarseTrees, K->getForces()) ;
//  				solverConvergence = K->mgsolve(&cg, lastx, &mgs) ;
// 				solverConvergence = K->mgsolve(&mg, lastx, NULL) ;
			solverConvergence = K->cgsolve(lastx) ;
// 			}
			
			
		}
		else
			solverConvergence = K->cgsolve() ;
		
// 		dtree->project(coarseTrees[3], K->getDisplacements(), coarseAssemblies[3]->getDisplacements(), false) ;
		Vector r = K->getMatrix()*K->getDisplacements()-K->getForces() ;
		double perror = residualError ;
		residualError = sqrt(parallel_inner_product(&r[0], &r[0], r.size())) ;
		if(perror > residualError || solverConvergence)
			reuseDisplacements = true;
	}
	
	enrichmentChange = false ;
	
	if(solverConvergence)
	{
		if(is2D())
		{
			std::vector<DelaunayTriangle *> elements = dtree->getElements() ;
			double volume = 0;	
			crackedVolume = 0 ;	
			damagedVolume = 0 ;	
			//this will update the state of all elements. This is necessary as 
			//the behaviour updates might depend on the global state of the 
			//simulation.
	// 		if(solverConverged())
	// 		{
				std::cerr << " stepping through elements... " << std::flush ;
				for(size_t i = 0 ; i < elements.size() ;i++)
				{	
					if(i%1000 == 0)
						std::cerr << "\r stepping through elements... " << i << "/" << elements.size() << std::flush ;
					elements[i]->step(dt, &K->getDisplacements()) ;
				}
				std::cerr << " ...done" << std::endl ;
				
				int fracturedCount = 0 ;
				int ccount = 0 ;
				for(size_t i = 0 ; i < elements.size() ;i++)
				{
					if(i%1000 == 0)
						std::cerr << "\r checking for fractures (1)... " << i << "/" << elements.size() << std::flush ;
					if(elements[i]->getBehaviour()->getFractureCriterion())
						elements[i]->getBehaviour()->getFractureCriterion()->step(elements[i]->getState()) ;
				}
				std::cerr << " ...done. " << std::endl ;
				for(size_t i = 0 ; i < elements.size() ;i++)
				{	
					if(i%10000 == 0)
						std::cerr << "\r checking for fractures (2)... " << i << "/" << elements.size() << std::flush ;
					if(elements[i]->getBehaviour()->type !=VOID_BEHAVIOUR )
					{
						volume += elements[i]->area() ;
						
						elements[i]->getBehaviour()->step(dt, elements[i]->getState()) ;
						if(elements[i]->getBehaviour()->changed())
						{
							needAssembly = true ;
							meshChange = true ;
							ccount++ ;
							ret = false ;
						}
						if(elements[i]->getBehaviour()->fractured())
						{
							fracturedCount++ ;
							crackedVolume +=  elements[i]->area() ;
						}
						else if(std::abs(elements[i]->getBehaviour()->getTensor(Point(1./3.,1./3.))[0][0] - elements[i]->getBehaviour()->param[0][0]) > 1e-12*elements[i]->getBehaviour()->getTensor(Point(1./3.,1./3.))[0][0] )
						{
							damagedVolume +=  elements[i]->area() ;
						}
					}
					else if (elements[i]->getBehaviour()->fractured())
						crackedVolume +=  elements[i]->area() ;
					
				}
				std::cerr << " ...done. " << ccount << " elements changed."<< std::endl ;
				for(size_t i = 0 ; i < elements.size() ;i++)
					elements[i]->clearVisited() ;
		// 		std::cout << " Fractured " << fracturedCount << " Elements" << std::endl ;
		// 		std::cout << " Fractured Fraction " <<  crackedVolume / volume << std::endl ;
				
				for(size_t i = 0 ; i< tree.size() ; i++)
				{
					if(tree[i]->isEnrichmentFeature)
					{
						dynamic_cast<EnrichmentFeature *>(tree[i])->step(dt, &K->getForces(), dtree) ;
						bool moved = dynamic_cast<EnrichmentFeature *>(tree[i])->moved() ;
						enrichmentChange = enrichmentChange || moved;
						if(moved)
						{
							if(useMultigrid)
							{
								for(size_t j = 0 ; j < coarseTrees.size() ; j++)
								{
									dynamic_cast<EnrichmentFeature *>(tree[i])->step(dt, &coarseAssemblies[j]->getForces(), coarseTrees[j]) ;
								}
							}
						}
						needAssembly = true ;
					}
				}
	// 		}

	// 		CrackInitiation().step(.01, dtree) ;
			
		}
		else if(is3D())
		{
			
			std::vector<DelaunayTetrahedron *> elements = dtree3D->getElements() ;
			
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
					bool moved = 
					enrichmentChange = enrichmentChange || dynamic_cast<EnrichmentFeature *>(tree[i])->moved() ;
					if(enrichmentChange)
						
					needAssembly = true ;
				}
			}
		}
	}
	else
	{
		Vector dummyx(0., K->getDisplacements().size()) ;
		if(is2D())
		{
			ret = false ;
			std::vector<DelaunayTriangle *> elements = dtree->getElements() ;
			double volume = 0;	
			crackedVolume = 0 ;	
			damagedVolume = 0 ;	
			//this will update the state of all elements. This is necessary as 
			//the behaviour updates might depend on the global state of the 
			//simulation.
	// 		if(solverConverged())
	// 		{
				std::cerr << " stepping through elements... " << std::flush ;
				for(size_t i = 0 ; i < elements.size() ;i++)
				{	
					if(i%1000 == 0)
						std::cerr << "\r stepping through elements... " << i << "/" << elements.size() << std::flush ;
					elements[i]->step(0., &dummyx) ;
				}
				std::cerr << " ...done" << std::endl ;
					
				for(size_t i = 0 ; i < elements.size() ;i++)
					elements[i]->clearVisited() ;

			
		}
		else if(is3D())
		{
			ret = false ;
			std::vector<DelaunayTetrahedron *> elements = dtree3D->getElements() ;
			
			//this will update the state of all elements. This is necessary as 
			//the behaviour updates might depend on the global state of the 
			//simulation.

			//this will update the state of all elements. This is necessary as 
			//the behaviour updates might depend on the global state of the 
			//simulation.
			std::cerr << " stepping through elements... " << std::flush ;
			for(size_t i = 0 ; i < elements.size() ;i++)
			{	
				if(i%1000 == 0)
					std::cerr << "\r stepping through elements... " << i << "/" << elements.size() << std::flush ;
				elements[i]->step(0., &dummyx) ;
			}
			std::cerr << " ...done" << std::endl ;

			// 		std::cout << " Fractured " << fracturedCount << " Elements" << std::endl ;
			// 		std::cout << " Fractured Fraction " <<  crackedVolume / volume << std::endl ;
			
		}
	}
	
	if(useMultigrid && meshChange)
	{
		if(is2D())
		{

			for(size_t j = 0 ; j < coarseTrees.size() ;j++)
			{
				std::cerr << " stepping through elements... grid " << j << std::flush ;
				std::vector<DelaunayTriangle *> elements = coarseTrees[j]->getElements() ;
				for(size_t i = 0 ; i < elements.size() ;i++)
				{	
					if(i%1000 == 0)
						std::cerr << "\r stepping through  grid " << j <<", elements... " << i << "/" << elements.size() << std::flush ;
					if(elements[i]->getBehaviour()->type !=VOID_BEHAVIOUR )
						elements[i]->getBehaviour()->step(dt, elements[i]->getState()) ;
				}
				std::cerr << " ...done" << std::endl ;
			}
			
			
		}
		else if(is3D())
		{
			
			for(size_t j = 0 ; j < coarseTrees3D.size() ;j++)
			{
				std::cerr << " stepping through elements... grid " << j << std::flush ;
				std::vector<DelaunayTetrahedron *> elements = coarseTrees3D[j]->getElements() ;
				for(size_t i = 0 ; i < elements.size() ;i++)
				{	
					if(i%1000 == 0)
						std::cerr << "\r stepping through  grid " << j <<", elements... " << i << "/" << elements.size() << std::flush ;
					if(elements[i]->getBehaviour()->type !=VOID_BEHAVIOUR )
						elements[i]->getBehaviour()->step(dt, elements[i]->getState()) ;
				}
				std::cerr << " ...done" << std::endl ;
			}
		}
	}
	
	return ret ;
	
}

bool FeatureTree::stable(double dt)
{
	bool needAssemblyinit = needAssembly ;
	bool meshChangeinit = meshChange ;
	bool enrichmentChangeinit = enrichmentChange ;
	double crackedVolumeinit = crackedVolume ;	
	double damagedVolumeinit = damagedVolume ;
	bool stab = step(dt) ;
	if(meshChanged() && solverConverged())
		stepBack() ;
	needAssembly = true ;
	meshChange = meshChangeinit ;
	enrichmentChange = enrichmentChangeinit ;
	crackedVolume = crackedVolumeinit ;	
	damagedVolume = damagedVolumeinit ;
	
	return stab ;
}

double FeatureTree::getMaximumDisplacement() const
{
	if(is2D())
	{
		std::vector<DelaunayTriangle *> tri = dtree->getElements() ;
		
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
		std::vector<DelaunayTetrahedron *> tets = dtree3D->getElements() ;
		
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
		std::vector<DelaunayTriangle *> tri = dtree->getElements() ;
		
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
		std::vector<DelaunayTetrahedron *> tets = dtree3D->getElements() ;
		
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
	if(initialized)
		return ;
	if(is2D())
	{
		std::vector<DelaunayTriangle *> triangles = this->dtree->getElements() ;
		std::cout << " initialising..." ;
		for(size_t i = 0 ; i < triangles.size() ;i++)
		{
			triangles[i]->refresh(father2D);
			if(i%1000 == 0)
				std::cout << "\r initialising... element " << i << "/" << triangles.size() << std::flush ;
			triangles[i]->getState().initialize() ;
		}
		std::cout << "\r initialising... element " << triangles.size() << "/" << triangles.size() << std::endl ;
		if(useMultigrid)
		{
			for(size_t i = 0 ; i < coarseTrees.size() ; i++)
			{
				triangles = coarseTrees[i]->getElements() ;
				
				for(size_t j = 0 ; j < triangles.size() ;j++)
				{
					triangles[j]->refresh(father2D);
					triangles[j]->getState().initialize() ;
				}
			}
		}
	}
	
	if(is3D())
	{
		std::vector<DelaunayTetrahedron *> tets = this->dtree3D->getElements() ;
		std::cout << " initialising..." ;
		for(size_t i = 0 ; i < tets.size() ;i++)
		{
			tets[i]->refresh(father3D);
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
	double pointDensity = 0 ; 
	if(is2D())
		pointDensity = .2*sqrt(tree[0]->area()/(tree[0]->getBoundingPoints().size()+tree[0]->getInPoints().size())) ;
	else
		pointDensity = .2*pow(tree[0]->volume()/(tree[0]->getBoundingPoints().size()+tree[0]->getInPoints().size()), .33333333333) ;
		
	std::cout << "space meshed with " << pointDensity << " points per unit length"<< std::endl ;

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
	std::vector<Feature *> nullFatherFeatures ;
	for(size_t i  = 1 ; i < tree.size() ; i++)
	{
		if(!tree[i]->isEnrichmentFeature && !tree[i]->isVirtualFeature)
		{
			if(tree[i]->getFather() == NULL)
				nullFatherFeatures.push_back(tree[i]);
		}
	}
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
				
				std::vector<Geometry *> potentialFeaturestmp  ;
				std::vector<Feature *> potentialFeatures ;
	
				if(is2D())
					potentialFeaturestmp = grid->coOccur(tree[i]->getBoundingPoint(j)) ;
				else
				{
					potentialFeaturestmp = grid3d->coOccur(tree[i]->getBoundingPoint(j)) ;
				}
				
				for(size_t l = 0 ; l < potentialFeaturestmp.size() ; l++)
					potentialFeatures.push_back(dynamic_cast<Feature *>(potentialFeaturestmp[l]) ) ;
				
				std::vector<Feature *> potentialChildren ;
				for(size_t l = 0 ; l < potentialFeatures.size() ; l++)
				{
					if(!potentialFeatures[l]->isEnrichmentFeature 
						&& std::binary_search(descendants.begin(), descendants.end(), potentialFeatures[l]))
					{
						potentialChildren.push_back(potentialFeatures[l]) ;
					}
				}
				
				for(size_t k  =  0 ; k <  potentialChildren.size() ; k++)
				{
					
					if((!potentialChildren[k]->isVirtualFeature 
					    && potentialChildren[k]->inBoundary(tree[i]->getBoundingPoint(j), pointDensity)) 
					   || (potentialChildren[k]->isVirtualFeature 
					       && tree[i]->isVirtualFeature 
					       && (dynamic_cast<VirtualFeature *>(potentialChildren[k])->getSource() 
					           != dynamic_cast<VirtualFeature *>(tree[i])->getSource())
					       && potentialChildren[k]->inBoundary(tree[i]->getBoundingPoint(j), pointDensity)
					      )
					   )
					{
						if(potentialChildren[k]->getBoundingPoints().size())
						{
							isIn = true ;
							break ;
						}
					}
				}
				
				if(i != 0 && !inRoot(tree[i]->getBoundingPoint(j)))
					isIn = true ;
				if(!isIn && tree[i]->isVirtualFeature && !tree[i]->in(tree[i]->getBoundingPoint(j)))
					isIn = true ;
				if(tree[i]->getFather() && tree[i]->getFather()->onBoundary(tree[i]->getBoundingPoint(j), pointDensity))
					isIn = true ;
				if(!isIn && i != 0 && tree[0]->onBoundary(tree[i]->getBoundingPoint(j), pointDensity))
				{
// 					Point proj(tree[i]->getBoundingPoint(j)) ;
// 					tree[0]->project(&proj) ;
// 					if(dist(proj, tree[i]->getBoundingPoint(j)) > 2.*POINT_TOLERANCE)
						isIn = true ;
// 					else
// 					{
// 						isIn = false ;
// 					}
				}
				
				if(tree[i]->getFather() == NULL && i != 0)
					isIn = false ;
				
				int nullFeatureIndex = -1 ;
				for(size_t k = 0 ; k < nullFatherFeatures.size() ; k++)
				{
					if(tree[i] == nullFatherFeatures[k])
					{
						nullFeatureIndex = k ;
						break ;
					}
				}
				for(size_t k = 0 ; k < nullFatherFeatures.size() ; k++)
				{
					if(tree[i] != nullFatherFeatures[k] && nullFeatureIndex < k)
					{
						Point proj(tree[i]->getBoundingPoint(j)) ;
						nullFatherFeatures[k]->project(&proj) ;
						if(dist(proj, tree[i]->getBoundingPoint(j)) < 2.*POINT_TOLERANCE)
						{
							isIn = true ;
							break ;
						}
						if(nullFatherFeatures[k]->in(tree[i]->getBoundingPoint(j)))
						{
							isIn = true ;
							break ;
						}
					}
				}
				
				if(!isIn && tree[i]->getFather() && tree[i]->getFather()->onBoundary(tree[i]->getBoundingPoint(j), pointDensity))
				{
					Point proj(tree[i]->getBoundingPoint(j)) ;
					tree[i]->getFather()->project(&proj) ;
					if(dist(proj, tree[i]->getBoundingPoint(j)) > 2.*POINT_TOLERANCE)
						isIn = true ;
					else
					{
						isIn = false ;
					}
				}
				if(!isIn)
				{
					if(tree[i]->getBoundingPoint(j) == Point(-250, 0))
						std::cout << i << " found" << std::endl ;
					bpcount++ ;
					meshPoints.push_back(std::pair<Point *, Feature *>(&tree[i]->getBoundingPoint(j), this->tree[i])) ;
					if(i == 0)
						basepoints++ ;
				}
			}

			for(size_t j  =  0 ; j <  tree[i]->getInPoints().size() ; j++)
			{
//                            if(i>0)
//                                tree[i]->getInPoint(j).print() ;
                            bool isIn = false ;
				std::vector<Geometry *> potentialFeaturestmp  ;
				if(is2D())
					potentialFeaturestmp = grid->coOccur(tree[i]->getInPoint(j)) ;
				else
					potentialFeaturestmp = grid3d->coOccur(tree[i]->getInPoint(j)) ;
				
				std::vector<Feature *> potentialFeatures ;
				for(size_t k = 0 ; k < potentialFeaturestmp.size() ; ++k)
					potentialFeatures.push_back(dynamic_cast<Feature *>(potentialFeaturestmp[k])) ;
				
				std::vector<Feature *> potentialChildren ;
				
				for(size_t l = 0 ; l < potentialFeatures.size() ; l++)
				{
					if(!potentialFeatures[l]->isVirtualFeature
					   &&!potentialFeatures[l]->isEnrichmentFeature 
					   && std::binary_search(descendants.begin(), descendants.end(), potentialFeatures[l] ))
						potentialChildren.push_back(potentialFeatures[l]) ;
				}
				for(size_t k  =  0 ; k <  potentialChildren.size() ; k++)
				{
					if(
					    (
					      !potentialChildren[k]->isVirtualFeature 
					      && potentialChildren[k]->inBoundary(tree[i]->getInPoint(j), pointDensity)
					    ) 
					   || 
					      (
					        potentialChildren[k]->isVirtualFeature 
					        && tree[i]->isVirtualFeature 
					        && 
					        (
					          dynamic_cast<VirtualFeature *>(potentialChildren[k])->getSource() 
					          != dynamic_cast<VirtualFeature *>(tree[i])->getSource()
					        )
					        && potentialChildren[k]->inBoundary(tree[i]->getInPoint(j), pointDensity)
					      )
					  )
					{
						if(potentialChildren[k]->getBoundingPoints().size())
						{
							isIn = true ;
							break ;
						}
					}
				}
				
				if(i != 0 && !inRoot(tree[i]->getInPoint(j)))
					isIn = true ;
				if(tree[i]->getFather() && tree[i]->getFather()->onBoundary(tree[i]->getInPoint(j), pointDensity))
					isIn = true ;
				if(tree[i]->isVirtualFeature && !tree[i]->in(tree[i]->getInPoint(j)))
					isIn = true ;
					
				if(i != 0 && tree[0]->onBoundary(tree[i]->getInPoint(j), pointDensity))
					isIn = true ;

				if(tree[i]->getFather() == NULL && i != 0)
					isIn = false ;
				
				int nullFeatureIndex = -1 ;
				for(size_t k = 0 ; k < nullFatherFeatures.size() ; k++)
				{
					if(tree[i] == nullFatherFeatures[k])
					{
						nullFeatureIndex = k ;
						break ;
					}
				}
				for(size_t k = 0 ; k < nullFatherFeatures.size() ; k++)
				{
					if(tree[i] != nullFatherFeatures[k] && nullFeatureIndex < k)
					{
						Point proj(tree[i]->getInPoint(j)) ;
						nullFatherFeatures[k]->project(&proj) ;
						if(dist(proj, tree[i]->getInPoint(j)) < 2.*POINT_TOLERANCE)
						{
							isIn = true ;
							break ;
						}
						if(nullFatherFeatures[k]->in(tree[i]->getInPoint(j)))
						{
							isIn = true ;
							break ;
						}
					}
				}
				
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
	
	if(is2D())
	{
		//this approach maximises the number of coincident points between the coarse grids.
		int ndivs = 4/*std::max(tree[0]->getBoundingPoints().size()/8)*/ ;
		while( ndivs*ndivs < meshPoints.size()/4)
		{
			coarseTrees.push_back(new StructuredMesh((max_x-min_x), (max_y-min_y), ndivs, Point((max_x+min_x)*.5, (max_y+min_y)*.5 ))) ;
			coarseAssemblies.push_back(new Assembly()) ;
			ndivs *= 4 ;
		}
	}
	else
	{
		
	}
	
	size_t count  = 0 ;

	if(computeIntersections)
	{
		for(size_t i = 1 ;  i < tree.size() ; i++)
		{
			if(!tree[i]->isEnrichmentFeature && !tree[i]->isVirtualFeature && tree[i]->getFather() != NULL)
			{
				std::vector<Geometry *> coOccuringFeaturestmp ;
				std::vector<Feature *> descendants = tree[i]->getDescendants() ;

				if(is3D())
					coOccuringFeaturestmp = grid3d->coOccur(tree[i]) ;
				else
					coOccuringFeaturestmp = grid->coOccur(tree[i]) ;
				
				std::vector<Feature *> coOccuringFeatures ;
				
				for(size_t k = 0 ; k < coOccuringFeaturestmp.size() ; k++)
					coOccuringFeatures.push_back(dynamic_cast<Feature *>(coOccuringFeaturestmp[k])) ;
				
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
								if(!descendants[l]->isVirtualFeature && descendants[l]->inBoundary(inter[k], pointDensity))
								{
									indescendants = true ;
									break ;
								}
							}
// 	
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
			if(!tree[i]->isEnrichmentFeature && tree[i]->getBoundingPoints().size() && !tree[i]->isVirtualFeature && tree[0]->intersects(tree[i])&& tree[i]->getFather() != NULL)
			{
				std::vector<Point> inter = tree[0]->intersection(tree[i]) ;
				std::vector<Feature *> descendants = tree[i]->getDescendants() ;
				std::vector<Feature *> fatherdescendants = tree[0]->getDescendants() ;
				for(size_t k = 0 ;  k < inter.size() ; k++)
				{
					bool indescendants = false ;
					for(size_t l = 0 ; l < descendants.size() ; l++)
					{
						if(descendants[l]->inBoundary(inter[k], pointDensity))
						{
							indescendants = true ;
							break ;
						}
					}
					for(size_t l = 0 ; l < fatherdescendants.size() ; l++)
					{
						if(fatherdescendants[l] != tree[i] && !fatherdescendants[l]->isVirtualFeature && fatherdescendants[l]->inBoundary(inter[k], pointDensity) && fatherdescendants[l]->getBoundingPoints().size() )
						{
							indescendants = true ;
							break ;
						}
					}
					
					if(is3D())
					{
					
						Point proj(inter[k]) ;
						tree[0]->project(&proj) ;
						Point proj0(inter[k]+Point(2.*POINT_TOLERANCE, 0, 0)) ;
						Point proj1(inter[k]+Point(-2.*POINT_TOLERANCE, 0, 0)) ;
						Point proj2(inter[k]+Point(0, 2.*POINT_TOLERANCE, 0)) ;
						Point proj3(inter[k]+Point(0, -2.*POINT_TOLERANCE, 0)) ;
						Point proj4(inter[k]+Point(0, 0, 2.*POINT_TOLERANCE)) ;
						Point proj5(inter[k]+Point(0, 0, -2.*POINT_TOLERANCE)) ;
						
						int position = tree[0]->in(proj0) 
						+ tree[0]->in(proj1)
						+ tree[0]->in(proj2)
						+ tree[0]->in(proj3)
						+ tree[0]->in(proj4)
						+ tree[0]->in(proj5) ;
						
						bool onSurface = (position == 5) ;
						bool onEdge = (position == 4) ;
						bool onVertex = (position == 3) ;
						proj0= (inter[k]+Point(pointDensity, 0, 0)) ;
						proj1= (inter[k]+Point(-pointDensity, 0, 0)) ;
						proj2= (inter[k]+Point(0, pointDensity, 0)) ;
						proj3= (inter[k]+Point(0, -pointDensity, 0)) ;
						proj4= (inter[k]+Point(0, 0, pointDensity)) ;
						proj5= (inter[k]+Point(0, 0, -pointDensity)) ;
						int tooClose =  tree[0]->in(proj0) 
						+ tree[0]->in(proj1)
						+ tree[0]->in(proj2)
						+ tree[0]->in(proj3)
						+ tree[0]->in(proj4)
						+ tree[0]->in(proj5) ;

						// no overlap with other features, intersection is indeed on the surface, and not too near another part of the surface
						if(!indescendants && squareDist3D(proj, inter[k]) < POINT_TOLERANCE*POINT_TOLERANCE && /*inRoot(inter[k]) && */((onSurface && tooClose == 5) || (onEdge && tooClose == 4) || onVertex))
						{
							Point *p = new Point(inter[k]) ;
							additionalPoints.push_back(p) ;
							++count ;
							meshPoints.push_back(std::make_pair(p, tree[i])) ;
						}
					}
					else
					{
						Point proj(inter[k]) ;
						tree[0]->project(&proj) ;
						Point proj0(inter[k]+Point(2.*POINT_TOLERANCE, 0, 0)) ;
						Point proj1(inter[k]+Point(-2.*POINT_TOLERANCE, 0, 0)) ;
						Point proj2(inter[k]+Point(0, 2.*POINT_TOLERANCE, 0)) ;
						Point proj3(inter[k]+Point(0, -2.*POINT_TOLERANCE, 0)) ;

						
						int position = tree[0]->in(proj0) 
						+ tree[0]->in(proj1)
						+ tree[0]->in(proj2)
						+ tree[0]->in(proj3);
						
						bool onEdge = (position == 3) ;
						bool onVertex = (position == 2) ;
						proj0= (inter[k]+Point(pointDensity, 0, 0)) ;
						proj1= (inter[k]+Point(-pointDensity, 0, 0)) ;
						proj2= (inter[k]+Point(0, pointDensity, 0)) ;
						proj3= (inter[k]+Point(0, -pointDensity, 0)) ;

						int tooClose =  tree[0]->in(proj0) 
						+ tree[0]->in(proj1)
						+ tree[0]->in(proj2)
						+ tree[0]->in(proj3);

						// no overlap with other features, intersection is indeed on the surface, and not too near another part of the surface
						if(!indescendants && squareDist3D(proj, inter[k]) < POINT_TOLERANCE*POINT_TOLERANCE && inRoot(inter[k]) && ( (onEdge && tooClose == 3) || onVertex))
						{
							Point *p = new Point(inter[k]) ;
							additionalPoints.push_back(p) ;
							++count ;
							meshPoints.push_back(std::make_pair(p, tree[i])) ;
						}
					}

				}
					
				if(count%100 == 0)
					std::cerr << "\r adding intersection points... " << count << std::flush ;
			}
	
		}
	}
	std::cerr << "\r adding intersection points... " << count << " ...done." << std::endl ;
	count = 0 ;

	//let us make sure we have no overlap
// 	std::stable_sort(meshPoints.begin(), meshPoints.end(), PairPointFeatureLess_Than_x()) ;
// 	std::stable_sort(meshPoints.begin(), meshPoints.end(), PairPointFeatureLess_Than_y()) ;
// 	if(is3D())
// 		std::stable_sort(meshPoints.begin(), meshPoints.end(), PairPointFeatureLess_Than_z()) ;
// 	std::deque<std::pair<Point *, Feature *> > ::iterator e = std::unique(meshPoints.begin(), meshPoints.end(), PairPointFeatureEqual());
// 	meshPoints.erase(e, meshPoints.end()) ;

// 	std::srand(1000) ;

	//shuffle for efficiency
	std::random_shuffle(meshPoints.begin(),meshPoints.end()) ;
	shuffleMeshPoints() ;
	
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
			if( (i - meshPoints.begin())%1000 == 0)
				std::cerr << "\r generating triangles... point " << count << "/" << meshPoints.size() << std::flush ;
			
			++count ;
			
			if(*i->first != bbox[0] &&
			   *i->first != bbox[2] &&
			   *i->first != bbox[4] &&
			   *i->first != bbox[6] && (inRoot(*i->first) || i->second->getFather() == NULL)
			  )
			{
				dtree->insert(i->first) ;
			}
		}
		
		std::cerr << "\r generating triangles.... point " << meshPoints.size()-3 << "/" << meshPoints.size()-4 << " ...done" << std::endl ;
		
		bool correct = false ;
		int tries = correctionSteps ;
		
		while(!correct && tries)
		{
			std::vector< DelaunayTriangle * > tets = dtree->getElements();
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

		std::pair<std::vector<int>,std::vector<int> > pb ;
		for(size_t i = 0 ; i < meshPoints.size() ; i++)
		{
			if(meshPoints[i].first == NULL)
				pb.first.push_back(i) ;
			if(meshPoints[i].second == NULL)
				pb.second.push_back(i) ;
		}
		std::vector<Point *> toInsert ;
		for( std::deque<std::pair<Point *, Feature *> >::iterator i = meshPoints.begin()+8 ; i != this->meshPoints.end(); ++i)
		{
			if((i - meshPoints.begin())%1000 == 0)
				std::cerr << "\r generating tetrahedrons... point " << ++count*1000 << "/" << meshPoints.size()-8 <<  std::flush ;
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
					toInsert.push_back(i->first) ;
				}
			}
		}
// 		int delta = 1 ;
// 		while(!toInsert.empty() && delta > 0)
// 		{
// 			std::random_shuffle(toInsert.begin(), toInsert.end()) ;
// 			std::vector<Point *> newToInsert ;
// 			for( std::vector<Point *>::iterator i = toInsert.begin() ; i != toInsert.end(); ++i)
// 			{
// 				if((i - toInsert.begin())%1000 == 0)
// 					std::cerr << "\r generating tetrahedrons... point " << ++count*1000 << "/" << toInsert.size() <<  std::flush ;
// 				if(**i != bbox[0] &&
// 				**i != bbox[1] &&
// 				**i != bbox[2] &&
// 				**i != bbox[3] &&
// 				**i != bbox[4] &&
// 				**i != bbox[5] &&
// 				**i != bbox[6] &&
// 				**i != bbox[7]
// 				)
// 				{
// 					dtree3D->insert(*i) ; 
// 					if((*i)->id == -1)
// 					{
// // 						std::cout << "insertion failed" << std::endl ;
// 						newToInsert.push_back(*i) ;
// 					}
// 				}
// 			}
// 			delta = toInsert.size()-newToInsert.size() ;
// 			toInsert = newToInsert ;
// 		}
		
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
			std::vector< DelaunayTetrahedron * > tets = dtree3D->getElements();
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
//		dtree3D->purge() ;
		
	}
}




std::vector<DelaunayTriangle *> FeatureTree::getTriangles(int g)
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
		if(g == -1)
			return dtree->getElements() ;
		else if(coarseTrees.size() > g)
			return coarseTrees[g]->getElements() ;
		else if(!coarseTrees.empty())
			return coarseTrees.back()->getElements() ;
		else
			return dtree->getElements() ;
	}
	else
		return std::vector<DelaunayTriangle *>(0) ;
}

std::vector<DelaunayTetrahedron *> FeatureTree::getTetrahedrons(int g)
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
		
		if(g == -1)
			return dtree3D->getElements() ;
		else if(coarseTrees3D.size() > g)
			return coarseTrees3D[g]->getElements() ;
		else if(!coarseTrees3D.empty())
			return coarseTrees3D.back()->getElements() ;
		else
			return dtree3D->getElements() ;
		
		return dtree3D->getElements() ;
	}
	else
		return std::vector<DelaunayTetrahedron *>(0) ;
}

void FeatureTree::shuffleMeshPoints()
{
	return ;
	std::cout << "shuffling mesh points... " ;
	std::deque<std::pair<Point *, Feature * > > shuffled ;
	for(size_t i = 0 ; i < meshPoints.size() ; i++)
		shuffled.push_back(meshPoints[i]) ;

	meshPoints.clear() ;
	
	std::random_shuffle(shuffled.begin(),shuffled.end()) ;

	std::vector<bool> visited ;
	for(size_t i = 0 ; i < shuffled.size() ; i++)
		visited.push_back(true) ;

	size_t ix = 0 ;
	size_t iy = 0 ;
	size_t iz = 0 ;

	size_t p = 0 ;

	size_t np = shuffled.size() / 2 ;
	if(is2D())
	{
		np =  std::pow(np, 0.5) + 1 ;
		Grid * shufflingGrid = new Grid(static_cast<Sample *>(tree[0])->width()*1.01,static_cast<Sample *>(tree[0])->height()*1.01,np,tree[0]->getCenter()) ;
		while(meshPoints.size() < shuffled.size())
		{
			Point ptest(shuffled[p].first->x,shuffled[p].first->y) ;
			if((visited[p]) && (shufflingGrid->pixels[ix][iy]->coOccur(ptest)))
			{
				visited[p] = false ;
				meshPoints.push_back(shuffled[p]) ;
				ix++ ;
				if(ix == shufflingGrid->getLengthX())
				{
					ix = 0 ;
					iy++ ;
					if(iy == shufflingGrid->getLengthY())
						iy = 0 ;
				}
				p = 0 ;
			} else {
				p++ ;
				if (p == shuffled.size())
				{
					p = 0 ;
					ix++ ;
					if(ix == shufflingGrid->getLengthX())
					{
						ix = 0 ;
						iy++ ;
						if(iy == shufflingGrid->getLengthY())
							iy = 0 ;
					}
//					shufflingGrid->pixels[ix][iy]->print() ;
				}
			}
		}
		delete shufflingGrid ;
	}	
	if(is3D())
	{
		np = std::pow(np, 0.3333333) +1 ;
		Grid3D * shufflingGrid = new Grid3D(static_cast<Sample3D *>(tree[0])->getXSize()*1.01,
						    static_cast<Sample3D *>(tree[0])->getYSize()*1.01,
						    static_cast<Sample3D *>(tree[0])->getZSize()*1.01,np,tree[0]->getCenter()) ;
		while(meshPoints.size() < shuffled.size())
		{
			Point ptest(shuffled[p].first->x,shuffled[p].first->y,shuffled[p].first->z) ;
			if((visited[p]) && (shufflingGrid->pixels[ix][iy][iz]->coOccur(ptest)))
			{
				visited[p] = false ;
				meshPoints.push_back(shuffled[p]) ;
				ix++ ;
				if(ix == shufflingGrid->getLengthX())
				{
					ix = 0 ;
					iy++ ;
					if(iy == shufflingGrid->getLengthY())
					{
						iy = 0 ;
						iz++ ;
						if(iz == shufflingGrid->getLengthY())
							iz = 0 ;
					}
				}
				p = 0 ;
			} else {
				p++ ;
				if (p == shuffled.size())
				{
					p = 0 ;
					ix++ ;
					if(ix == shufflingGrid->getLengthX())
					{
						ix = 0 ;
						iy++ ;
						if(iy == shufflingGrid->getLengthY())
						{
							iy = 0 ;
							iz++ ;
							if(iz == shufflingGrid->getLengthZ())
								iz = 0 ;
						}
					}
				}
			}
		}
	}
	std::cout << "done... " << std::endl ;
	
}


