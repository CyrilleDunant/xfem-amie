// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "features.h"
#include "crackinitiation.h"
#include "layeredinclusion.h"
#include "../physics/void_form.h"

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
	double nu = 0.32 ;
	double E = 2 ;
// 	double nu = 0.5 ;
// 	double E = 1 ;
	Matrix cg(3,3) ;
	
// 	double coef = E/(1.-nu*nu) ;
// 	
// 	cg[0][0] = 1. ; cg[0][1] = nu ; cg[0][2] = 0 ; 
// 	cg[1][0] = nu ; cg[1][1] = 1. ; cg[1][2] = 0 ; 
// 	cg[2][0] = 0 ; cg[2][1] = 0 ; cg[2][2] = (1.-nu)/2. ; 
// 	cg *=coef ;
	
	double coef = E/((1.+nu)*(1.-2.*nu)) ;
	
	cg[0][0] = 1.-nu ; cg[0][1] = nu ;    cg[0][2] = 0 ; 
	cg[1][0] = nu ;    cg[1][1] = 1.-nu ; cg[1][2] = 0 ; 
	cg[2][0] = 0 ;     cg[2][1] = 0 ;     cg[2][2] = (1.-2.*nu)/2. ; 
	cg *=coef ;
	
// 	cg[0][0] = E/(1.-nu*nu) ; cg[0][1] =E/(1.-nu*nu)*nu ; cg[0][2] = 0 ;
// 	cg[1][0] = E/(1.-nu*nu)*nu ; cg[1][1] = E/(1.-nu*nu) ; cg[1][2] = 0 ; 
// 	cg[2][0] = 0 ; cg[2][1] = 0 ; cg[2][2] = E/(1.-nu*nu)*(1.-nu)/2. ; 
// 	cg[0][0] = E*(1.-nu)/((1.-nu)*(1.+nu)) ; cg[0][1] =E/((1.-nu)*(1.+nu))*nu ; cg[0][2] = 0 ;
// 	cg[1][0] = E/((1.-nu)*(1.+nu))*nu ; cg[1][1] = E*(1.-nu)/((1.-nu)*(1.+nu)) ; cg[1][2] = 0 ; 
// 	cg[2][0] = 0 ; cg[2][1] = 0 ; cg[2][2] = E*(1.-2*nu)/((1.-nu)*(1.+nu)) ; 
		
	this->behaviour = new Stiffness/*AndFracture*/(cg/*, sigmaRupt*/) ;
	
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

void Feature::setBehaviour(Form * f)
{
	this->behaviour = f ;
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
	std::vector<Feature *> ret ;
	for(size_t i = 0 ; i < m_c.size() ; i++)
	{
		ret.push_back(m_c[i]) ;
		std::vector<Feature *> to_add = m_c[i]->getDescendants() ;
		ret.insert(ret.end(), to_add.begin(), to_add.end()) ;
	}
	std::stable_sort(ret.begin(), ret.end()) ;
	std::vector<Feature *>::iterator e = std::unique(ret.begin(), ret.end()) ;
	ret.erase(e, ret.end()) ;
	
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

// void Feature::sample(size_t n)
// {
// 	this->sampleSurface(n) ;
// }

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
	if(this->dtree == NULL)
		this->generateElements() ;
	return this->dtree ;
}

DelaunayTree_3D * FeatureTree::getDelaunayTree3D()
{
	if(this->dtree3D == NULL)
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
		grid = new Grid(first->getBoundingBox()[1].x-first->getBoundingBox()[0].x, first->getBoundingBox()[1].y-first->getBoundingBox()[2].y, 20) ;
	if(is3D())
		grid3d =new Grid3D(first->getBoundingBox()[7].x-first->getBoundingBox()[0].x, first->getBoundingBox()[7].y-first->getBoundingBox()[0].y, first->getBoundingBox()[7].z-first->getBoundingBox()[0].z, 20);
	this->father3D = NULL;
	this->father2D = NULL ;
	this->elemOrder = LINEAR ;
	this->stitched = false ;
	this->renumbered = false ;
	this->needAssembly = true ;
	this->initialized = false ;

	meshChange = true ;
	solverConvergence = false ;
	enrichmentChange = true ;

	K = new Assembly() ;
	
	if(is3D())
		K->set3D() ;
	else
		K->set2D() ;
}

void FeatureTree::twineFeature(CompositeFeature * father, CompositeFeature * f)
{
	std::vector<VirtualFeature *> fatherComponents = father->getComponents() ;
	std::vector<VirtualFeature *> childComponents = f->getComponents() ;	
	addFeature(father, f) ;
	
	std::sort(father->getChildren().begin(), father->getChildren().end()) ;
	std::vector<Feature *>::iterator child = std::find(father->getChildren().begin(), father->getChildren().end(),
	fatherComponents[0]) ;
	father->getChildren().erase(child) ;
	child = std::find(tree.begin(), tree.end(), fatherComponents[0]) ;
	tree.erase(child) ;
	for(size_t i = 1 ; i< fatherComponents.size() ; i++)
	{
		std::vector<Feature *>::iterator child = std::find(fatherComponents[i-1]->getChildren().begin(),
			fatherComponents[i-1]->getChildren().end(), 
			fatherComponents[i]) ;
		fatherComponents[i-1]->getChildren().erase(child) ;
		child = std::find(tree.begin(), tree.end(), fatherComponents[i]) ;
		tree.erase(child) ;
	}
	
	addFeature(f, childComponents[0]) ;
	for(size_t i = 0 ; i< fatherComponents.size() ; i++)
	{
		addFeature(childComponents[i], fatherComponents[i]) ;
		if(i < childComponents.size()-2)
			addFeature(fatherComponents[i], childComponents[i+1]) ;
	}
}


void FeatureTree::addFeature(Feature * father, Feature * f)
{
	if(!tree.empty() && f->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
		grid->forceAdd(f) ;
	else if(!tree.empty())
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
}

void FeatureTree::setOrder(Order ord)
{
	this->elemOrder = ord ;
	
	if(father3D)
		delete father3D ;
	
	father3D = new TetrahedralElement(this->elemOrder) ;
	
	for(size_t k = 0 ; k < father3D->getShapeFunctions().size() ; k++)
	{
		father3D->getShapeFunction(k).compile() ;
		std::vector<Variable> vars ;
		vars.push_back(XI) ;
		vars.push_back(ETA) ;
		vars.push_back(ZETA) ;
		father3D->getShapeFunction(k).preCalculate(father3D->getGaussPoints(), vars) ;
	}
	
	if(father2D)
		delete father2D ;
	
	father2D = new TriElement(this->elemOrder) ;
	
	for(size_t k = 0 ; k < father2D->getShapeFunctions().size() ; k++)
	{
		father2D->getShapeFunction(k).compile() ;
		std::vector<Variable> vars ;
		vars.push_back(XI) ;
		vars.push_back(ETA) ;
		father2D->getShapeFunction(k).preCalculate(father2D->getGaussPoints(), vars) ;
	}
}

void FeatureTree::renumber()
{
	if(dtree != NULL)
	{
		std::vector<DelaunayTriangle *> triangles = this->dtree->getTriangles() ;
		size_t count = 0 ;
		std::cerr << " renumbering... " << std::flush ;
	// 	for(size_t i = 0 ;  i < this->meshPoints.size() ; i++)
	// 	{
	// 		std::cerr << "\r reseting point IDs... point " << ++count << "/"<< this->meshPoints.size() <<  std::flush ;
	// 		this->meshPoints[i].first->id = -1 ;
	// 	}
	// 	std::cerr << " ...done."<< std::endl ;
		
		std::vector<Point *> points ;
		for(size_t i = 0 ; i < triangles.size() ; i++)
		{
			for(size_t j = 0 ; j < triangles[i]->getBoundingPoints().size() ; j++)
			{
				triangles[i]->getBoundingPoint(j).id = -1 ;
				points.push_back(&triangles[i]->getBoundingPoint(j)) ;
			}
		}
		
		std::stable_sort(points.begin(), points.end()) ;
		std::vector<Point *>::iterator e = std::unique(points.begin(), points.end()) ;
		points.erase(e, points.end()) ;
		
		for(size_t i = 0 ; i < points.size() ; i++)
		{
				points[i]->id = count++ ;
		}
		
		for(size_t i = 0 ; i < triangles.size() ; i++)
		{
			for(size_t j = 0 ; j < triangles[i]->getBoundingPoints().size() ; j++)
			{
				if(triangles[i]->getBoundingPoint(j).id == -1 && inRoot(triangles[i]->getBoundingPoint(j)))
					triangles[i]->getBoundingPoint(j).id = count++ ;
			}
		}
		
		
		this->dtree->global_counter = count ;
		
		std::cerr << count*2 << " ...done " << std::endl ;

	}
	else if (dtree3D !=NULL)
	{
		std::vector<DelaunayTetrahedron *> tets = this->dtree3D->getTetrahedrons() ;
		size_t count = 0 ;
		std::cerr << " renumbering... " << std::flush ;
// 	for(size_t i = 0 ;  i < this->meshPoints.size() ; i++)
// 	{
// 		std::cerr << "\r reseting point IDs... point " << ++count << "/"<< this->meshPoints.size() <<  std::flush ;
// 		this->meshPoints[i].first->id = -1 ;
// 	}
// 	std::cerr << " ...done."<< std::endl ;
		
		std::vector<Point *> points ;
		for(size_t i = 0 ; i < tets.size() ; i++)
		{
			for(size_t j = 0 ; j < tets[i]->getBoundingPoints().size() ; j++)
			{
				tets[i]->getBoundingPoint(j).id = -1 ;
				points.push_back(&tets[i]->getBoundingPoint(j)) ;
			}
		}
		std::stable_sort(points.begin(), points.end()) ;
		std::vector<Point *>::iterator e = std::unique(points.begin(), points.end()) ;
		points.erase(e, points.end()) ;
		
		for(size_t i = 0 ; i < points.size() ; i++)
		{
			bool inNothing =false;
			for(size_t k = 0 ; k < this->tree.size() ; k++)
			{
				if(tree[k]->isVoid(*points[i]))
				{
					inNothing = true ;
					break;
				}
			}
			if(!inNothing)
				points[i]->id = count++ ;
			
		}
		
		this->dtree3D->global_counter = count ;
		
		std::cerr << count*4 << " ...done " << std::endl ;
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
	if(dtree != NULL)
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
				
			}
			stitched  = true ;	
			return ;
			for(size_t j = 1 ; j < this->tree.size() ; j++)
			{
				if(!tree[j]->isEnrichmentFeature)
				{
					
					std::vector<DelaunayTriangle *> triangles = this->tree[j]->getTriangles(dtree) ;
	
					for(size_t i = 0 ; i < triangles.size() ; i++)
					{
						if(triangles[i]->Triangle::intersects(dynamic_cast<Geometry *>(tree[j])))
						{
							Point proj_0(*triangles[i]->first) ;
							tree[j]->project(&proj_0) ;
							Point proj_1(*triangles[i]->second) ;
							tree[j]->project(&proj_1) ;
							Point proj_2(*triangles[i]->third) ;
							tree[j]->project(&proj_2) ;
							bool changed  = true;
							
							if(squareDist(&proj_0 , triangles[i]->first ) < 0.0000001 && 
								squareDist(&proj_1 , triangles[i]->second) < 0.0000001 )
							{
								count+=changed ; 
								changed = false ;
								Point test = triangles[i]->getBoundingPoint(1) ;
								tree[j]->project(&test) ;
								if (inRoot(test))
								{
									tree[j]->project(&triangles[i]->getBoundingPoint(1)) ;
									if(elemOrder >= CONSTANT_TIME_LINEAR)
										tree[j]->project(&triangles[i]->getBoundingPoint(7)) ;
									triangles[i]->moved = true ;
								}
		// 						std::cerr << "--> " << (*triangles)[i]->getBoundingPoint(1)->x << ", " << (*triangles)[i]->getBoundingPoint(1)->y << std::endl ;
							}
							if(squareDist(&proj_1 , triangles[i]->second) < 0.0000001 && 
								squareDist(&proj_2 , triangles[i]->third) < 0.0000001 )
							{
								count+=changed ; 
								changed = false ;								
								Point test = triangles[i]->getBoundingPoint(3) ;
								tree[j]->project(&test) ;
								if (inRoot(test))
								{
									tree[j]->project(&triangles[i]->getBoundingPoint(3)) ;
									if(elemOrder >= CONSTANT_TIME_LINEAR)
										tree[j]->project(&triangles[i]->getBoundingPoint(9)) ;
									triangles[i]->moved = true ;
								}
								
		// 						std::cerr << "--> " << (*triangles)[i]->getBoundingPoint(3)->x << ", " << (*triangles)[i]->getBoundingPoint(3)->y << std::endl ;
							}
							if(squareDist(&proj_2 , triangles[i]->third) < 0.0000001 && 
							   squareDist(&proj_0, triangles[i]->first) < 0.0000001) 
							{
								count+=changed ; 
								changed = false ;								
								Point test = triangles[i]->getBoundingPoint(5) ;
								tree[j]->project(&test) ;
								if (inRoot(test))
								{
									tree[j]->project(&triangles[i]->getBoundingPoint(5)) ;
									if(elemOrder >= CONSTANT_TIME_LINEAR)
										tree[j]->project(&triangles[i]->getBoundingPoint(11)) ;
									triangles[i]->moved = true ;
								}
								
		// 						std::cerr << "--> " << (*triangles)[i]->getBoundingPoint(5)->x << ", " << (*triangles)[i]->getBoundingPoint(5)->y << std::endl ;
							}
							
						}
						if(count % 1000 == 0)
							std::cerr << "\r projecting points on boundaries... triangle " << count << "/" << triangles.size() << " feature " << i << std::flush ; 
						
					}
				}
			}
		}
	}
	else
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
				
			}
			stitched = true ;
#warning there is a problem when the points get projected
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
						    squareDist(&proj_0 , tets[i]->first ) < 1e-12 && 
						    squareDist(&proj_1 , tets[i]->second) < 1e-12 
						  )
						{
							count +=1; 
							tree[j]->project(&tets[i]->getBoundingPoint(1)) ;
							if(elemOrder >= CONSTANT_TIME_LINEAR)
								tree[j]->project(&tets[i]->getBoundingPoint(11)) ;
							tets[i]->moved = true ;
						}
						if(
						    squareDist(&proj_0 , tets[i]->first ) < 1e-12 && 
						    squareDist(&proj_2 , tets[i]->third) < 1e-12 
						  )
						{
							count +=1; 
							tree[j]->project(&tets[i]->getBoundingPoint(9)) ;
							if(elemOrder >= CONSTANT_TIME_LINEAR)
								tree[j]->project(&tets[i]->getBoundingPoint(19)) ;
							tets[i]->moved = true ;
						}
						if(
						    squareDist(&proj_0 , tets[i]->first ) < 1e-12 && 
						    squareDist(&proj_3 , tets[i]->fourth) < 1e-12 
						  )
						{
							count +=1; 
							tree[j]->project(&tets[i]->getBoundingPoint(7)) ;
							if(elemOrder >= CONSTANT_TIME_LINEAR)
								tree[j]->project(&tets[i]->getBoundingPoint(17)) ;
							tets[i]->moved = true ;
						}
						if(
						    squareDist(&proj_1 , tets[i]->second ) < 1e-12 && 
						    squareDist(&proj_3 , tets[i]->fourth) < 1e-12 
						  )
						{
							count +=1; 
							tree[j]->project(&tets[i]->getBoundingPoint(8)) ;
							if(elemOrder >= CONSTANT_TIME_LINEAR)
								tree[j]->project(&tets[i]->getBoundingPoint(18)) ;
							tets[i]->moved = true ;
						}
						if(
						    squareDist(&proj_1 , tets[i]->second ) < 1e-12 && 
						    squareDist(&proj_2 , tets[i]->third) < 1e-12
						  )
						{
							count +=1; 
							tree[j]->project(&tets[i]->getBoundingPoint(3)) ;
							if(elemOrder >= CONSTANT_TIME_LINEAR)
								tree[j]->project(&tets[i]->getBoundingPoint(13)) ;
							tets[i]->moved = true ;
						}
						if(
						    squareDist(&proj_3 , tets[i]->fourth ) < 1e-12 && 
						    squareDist(&proj_2 , tets[i]->third) < 1e-12 
						  )
						{
							count +=1; 
							tree[j]->project(&tets[i]->getBoundingPoint(5)) ;
							if(elemOrder >= CONSTANT_TIME_LINEAR)
								tree[j]->project(&tets[i]->getBoundingPoint(15)) ;
							tets[i]->moved = true ;
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
		for(size_t i  = 1 ; i < this->tree.size() ; i ++)
		{
// 			tree[i]->sample(std::max((size_t)(1.5*n*pow(tree[i]->area()/total_area, .6)),(size_t)10)) ;
// 			if((size_t)((double)n*sqrt(tree[i]->area()/(total_area))) >= 8)
// 				tree[i]->sample((size_t)((double)n*sqrt(tree[i]->area()/(total_area)))) ;
			tree[i]->sample(std::max((size_t)((double)n*sqrt(tree[i]->area()/(total_area))),(size_t)8)) ;
// 			std::cerr << std::max((size_t)((double)n*tree[i]->area()/(.9*total_area)),(size_t)12) << std::endl ;
	// 		tree[i]->sample(sqrt(n)*10) ;
	// 		tree[i]->sampleSurface(n) ;
		}
	}
	else if (is3D())
	{
		std::cerr << "3D features" << std::endl ;
		double total_area = tree[0]->area() ;
// 		total_volume = sqrt(total_volume) ;
		tree[0]->sample(n) ;
		for(size_t i  = 1 ; i < this->tree.size() ; i ++)
		{
			if(inRoot(this->tree[i]->getCenter()))
			{
				double a =tree[i]->area() ;

				tree[i]->sample((size_t)round(
					                               (2.*a/(total_area))*n
											)
								) ;
			}
			else
			{
				double a =tree[i]->area() ;

					tree[i]->sample((size_t)round(
					                               (2.*a/(total_area))*n
											)
									) ;
			}
// 			std::cerr << (size_t)((double)n*pow(v/total_volume, .6)/2.) << std::endl ;
	// 		tree[i]->sample(sqrt(n)*10) ;
	// 		tree[i]->sampleSurface(n) ;
		}
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
						//std::cerr << "inserting..." ; toAdd[i].print() ; std::cerr << std::endl ;
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
					double rand0 = 0 ;/*((2.*random()/(RAND_MAX+1.0))-1.)*0.1 ;*/
					double rand1 = 0 ;/*((2.*random()/(RAND_MAX+1.0))-1.)*0.1 ;*/
					double rand2 = 0 ;/*((2.*random()/(RAND_MAX+1.0))-1.)*0.1 ;*/
						
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

	for(size_t i = 0 ; i < t->getBoundingPoints().size() ; i++)
		if( t->getBoundingPoint(i).id == -1)
		{
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
				
				for(size_t j = 0 ; j < targets[i]->getChildren().size() ; j++)
				{
					if(!targets[i]->getChild(j)->isEnrichmentFeature && targets[i]->getChild(j)->in(t->getCenter()))
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
	
	for(size_t i = 0 ; i < t->getBoundingPoints().size() ; i++)
		if( t->getBoundingPoint(i).id == -1)
		{
			return new VoidForm() ;
		}
	
	std::vector<Feature *> targets = grid3d->coOccur(static_cast<const Tetrahedron *>(t)) ;
	if(!targets.empty())
	{
		for(int i = targets.size()-1 ; i >=0  ; i--)
		{
			if (!targets[i]->isEnrichmentFeature && targets[i]->in(t->getCenter()))
			{
				
				
				bool notInChildren  = true ;
				
				for(size_t j = 0 ; j < targets[i]->getChildren().size() ; j++)
				{
					if(!targets[i]->getChild(j)->isEnrichmentFeature && targets[i]->getChild(j)->in(t->getCenter()))
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

Point * FeatureTree::checkElement( const DelaunayTetrahedron * t ) const
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
					Point *p = new Point(*t->getCircumCenter()) ;
					tree[i]->project(p);
					return p;
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


void FeatureTree::assemble()
{	
	std::vector<DelaunayTriangle *> triangles ; 
	std::vector<DelaunayTetrahedron *> tetrahedrons ; 
	
	if(!father3D)
		father3D = new TetrahedralElement(this->elemOrder) ;

	for(size_t k = 0 ; k < father3D->getShapeFunctions().size() ; k++)
	{
		father3D->getShapeFunction(k).compile() ;
		std::vector<Variable> vars ;
		vars.push_back(XI) ;
		vars.push_back(ETA) ;
		vars.push_back(ZETA) ;
		father3D->getShapeFunction(k).preCalculate(father3D->getGaussPoints(), vars) ;
	}
	
	if(!father2D)
		father2D = new TriElement(this->elemOrder) ;
	
	for(size_t k = 0 ; k < father2D->getShapeFunctions().size() ; k++)
	{
		father2D->getShapeFunction(k).compile() ;
		std::vector<Variable> vars ;
		vars.push_back(XI) ;
		vars.push_back(ETA) ;
		father2D->getShapeFunction(k).preCalculate(father2D->getGaussPoints(), vars) ;
	}
	
	if( this->dtree == NULL && this->dtree3D == NULL)
		this->generateElements() ;
	
	size_t en_counter = 0 ;
	if(!initialized)
	{
		if(this->dtree3D == NULL && this->dtree !=NULL)
		{
			std::cerr << " stich..." << std::flush ;
			if(!stitched)
				stitch() ;

			if(!renumbered)
				renumber() ;
			
			std::cerr << " refreshing..." << std::flush ;
			this->dtree->refresh(father2D, false) ;
			std::cerr << " ...done" << std::endl ;
			
			en_counter = this->dtree->global_counter ;
			
			triangles = this->dtree->getTriangles() ;
			std::cerr << " setting behaviours..." << std::flush ;
			for(size_t i = 0 ; i < triangles.size() ;i++)
			{
				if (i%1000 == 0)
					std::cerr << "\r setting behaviours... triangle " << i << "/" << triangles.size() << std::flush ;
				
				if(!triangles[i]->getBehaviour())
					triangles[i]->setBehaviour(getElementBehaviour(triangles[i])) ;
				triangles[i]->getEnrichmentFunctions().clear() ;
// 				triangles[i]->getState().initialize() ;
			}
			std::cerr << " ...done" << std::endl ;
			
			initializeElements() ;
			
		}
		else if(this->dtree3D != NULL && this->dtree ==NULL)
		{
			
			std::cerr << " stich..." << std::flush ;
			if(!stitched)
				stitch() ;
			
			if(!renumbered)
				renumber() ;
			
			std::cerr << " refreshing..." << std::flush ;
			this->dtree3D->refresh(father3D) ;
			std::cerr << " ...done" << std::endl ;
			
			en_counter = this->dtree3D->global_counter ;
			
			tetrahedrons = this->dtree3D->getTetrahedrons() ;
			
			for(size_t i = 0 ; i < tetrahedrons.size() ;i++)
			{
				if(!tetrahedrons[i]->getBehaviour())
					tetrahedrons[i]->setBehaviour(getElementBehaviour(tetrahedrons[i])) ;
				tetrahedrons[i]->getEnrichmentFunctions().clear() ;
				
// 				triangles[i]->getState().initialize() ;
			}
			
			initializeElements() ;
		}
	}
	
	std::cerr << " enriching..." << std::flush ;		
	
	if(enrichmentChange)
	{

		if(this->dtree3D == NULL && this->dtree !=NULL)
		{
			en_counter = this->dtree->global_counter ;
			triangles = this->dtree->getTriangles() ;
			
			for(size_t i = 0 ; i < triangles.size() ;i++)
			{
				triangles[i]->getEnrichmentFunctions().clear() ;
				triangles[i]->getState().initialize() ;
			}
		}
		else
		{
			en_counter = this->dtree3D->global_counter ;
			tetrahedrons = this->dtree3D->getTetrahedrons() ;
			
			for(size_t i = 0 ; i < tetrahedrons.size() ;i++)
			{
				tetrahedrons[i]->getEnrichmentFunctions().clear() ;
				tetrahedrons[i]->getState().initialize() ;
			}
		}
		
		for(size_t i = 1 ; i < this->tree.size() ; i++)
		{

			if(this->tree[i]->isEnrichmentFeature)
			{
				dynamic_cast<EnrichmentFeature *>(this->tree[i])->enrich(en_counter, this->dtree) ;
			}
			
			if(i%10 == 0)
				std::cerr << "\r enriching... feature " << i+1 <<"/" << this->tree.size() << std::flush ;
		}

		std::cerr << " ...done" << std::endl ;

		for(size_t j = 0 ; j < triangles.size() ; j++)
		{
			if(	triangles[j]->getBehaviour()->type != VOID_BEHAVIOUR)
			{
				
				std::cerr << "\r compiling... triangle " << j+1 << "/" << triangles.size() << std::flush ;
				
				for(size_t k = 0 ; k < triangles[j]->getEnrichmentFunctions().size() ; k++)
				{
					triangles[j]->getEnrichmentFunction(k).compile() ;
				}
			}
		}

		std::cerr << " ...done" << std::endl ;
	}
	
	if(this->dtree3D == NULL && this->dtree !=NULL)
	{
		std::cerr << " refreshing..." << std::flush ;
		this->dtree->refresh(father2D, false) ;
		std::cerr << " ...done" << std::endl ;
	}
	else if(this->dtree3D != NULL && this->dtree ==NULL)
	{
		std::cerr << " refreshing..." << std::flush ;
		this->dtree3D->refresh(father3D) ;
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
				
				for(size_t k = 0 ; k < tets[j]->getEnrichmentFunctions().size() ; k++)
				{
					tets[j]->getEnrichmentFunction(k).compile() ;
				}
				
				K->add(tets[j]) ;
				elements3D.push_back(tets[j]) ;
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
		copy.rbegin()->setBehaviour(tris[i]->getBehaviour()->getCopy()) ;
		copy.rbegin()->getState().initialize() ;
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
		std::pair<Vector , Vector > stress_strain(Vector(0.f, 4*6*elements3D.size()), Vector(0.f, 4*6*elements3D.size())) ;
		
		for(size_t i  = 0 ; i < elements3D.size() ; i++)
		{
			std::valarray<Point *> pts(4) ;
			pts[0] =  elements3D[i]->first ;
			pts[1] =  elements3D[i]->second ;
			pts[2] =  elements3D[i]->third ;
			pts[4] =  elements3D[i]->fourth ;
			
			std::pair<Vector , Vector > str = elements3D[i]->getState().getStressAndStrain(pts) ;
			
			for(size_t j = 0 ; j < 24 ; j++)
			{
				stress_strain.first[i*4*6+j] = str.first[j] ;
				stress_strain.second[i*4*6+j] = str.second[j] ;
			}
			std::cerr << "\r computing strain+stress... element " << i+1 << "/" << elements3D.size() << std::flush ;
		}
		std::cerr << " ...done." << std::endl ;
		return stress_strain ;
	}
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

bool FeatureTree::step(double dt)
{
	bool ret = true ;
	if(true)
	{
		this->K->clear() ;
		assemble() ;
	}
	else
	{
		this->K->setBoundaryConditions() ;
	}
	
	needAssembly = true ;
// 	Vector displacements = this->K->solve(/**extforces*/Vector(0), 100000, true) ;
	
	meshChange = false ;
	solverConvergence = this->K->cgsolve() ;
	enrichmentChange = false ;

	if(is2D())
	{
		std::vector<DelaunayTriangle *> elements = dtree->getTriangles() ;
		
		double volume = 0;
		double crackedVolume = 0 ;
		
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
				if(elements[i]->getBehaviour()->changed() && elements[i]->getBehaviour()->fractured())
				{
					fracturedCount++ ;
					needAssembly = true ;
					meshChange = true ;
					ret = false ;
					crackedVolume +=  elements[i]->area() ;
				}
				else if(elements[i]->getBehaviour()->changed() )
				{
					meshChange = true ;
					needAssembly = true ;
					ret = false ;
				}
			}
			else if (elements[i]->getBehaviour()->type !=VOID_BEHAVIOUR && elements[i]->getBehaviour()->fractured())
				crackedVolume +=  elements[i]->area() ;
		}
		std::cerr << " ...done" << std::endl ;
		
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
		for(size_t i = 0 ; i < elements.size() ;i++)
		{	
			elements[i]->step(dt, &K->getDisplacements()) ;
		}
		
		int fracturedCount = 0 ;
		
		for(size_t i = 0 ; i < elements.size() ;i++)
		{	
			if(elements[i]->getBehaviour()->type !=VOID_BEHAVIOUR && !elements[i]->getBehaviour()->fractured())
			{
				elements[i]->getBehaviour()->step(dt, elements[i]->getState()) ;
				if(elements[i]->getBehaviour()->fractured())
				{
					fracturedCount++ ;
					needAssembly = true ;
					meshChange = true ;
					ret = false ;
				}
			}
		}
		
// 		for(size_t i = 0 ; i< tree.size() ; i++)
// 		{
// 			if(tree[i]->isEnrichmentFeature)
// 			{
// 				dynamic_cast<EnrichmentFeature *>(tree[i])->step(dt, &K->getForces(), dtree3D) ;
// 				needAssembly = dynamic_cast<EnrichmentFeature *>(tree[i])->moved() ;
// // 				ret = false ;
// 			}
// 		}
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
	if(dtree)
	{
		std::vector<DelaunayTriangle *> triangles = this->dtree->getTriangles() ;
		for(size_t i = 0 ; i < triangles.size() ;i++)
		{
			triangles[i]->getState().initialize() ;
		}
	}
	
	if(dtree3D)
	{
		std::vector<DelaunayTetrahedron *> triangles = this->dtree3D->getTetrahedrons() ;
		for(size_t i = 0 ; i < triangles.size() ;i++)
		{
			triangles[i]->getState().initialize() ;
		}
	}
	
	initialized = true ;
}

void FeatureTree::generateElements( size_t correctionSteps) 
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
	
	size_t basepoints = 0 ;
	std::cerr << " getting mesh points..." << std::flush ;
	for(size_t i  = 0 ; i < tree.size() ; i++)
	{
		std::cerr << "\r getting mesh points... feature " << i << "/"<< tree.size() << std::flush ;
		if(!tree[i]->isEnrichmentFeature)
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
						&& std::find(descendants.begin(), descendants.end(), potentialFeatures[l]) != descendants.end() )
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

					if( /*!potentialChildren[k]->isVirtualFeature
					    && */potentialChildren[k]->getBoundary()->in(tree[i]->getBoundingPoint(j)) )
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
					if(!potentialFeatures[l]->isEnrichmentFeature 
						&& std::find(descendants.begin(), descendants.end(), potentialFeatures[l]) != descendants.end() )
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
					if(/*!potentialChildren[k]->isVirtualFeature 
					   && */potentialChildren[k]->getBoundary()->in(tree[i]->getInPoint(j)) )
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

	
	for(size_t i = 0 ;  i < tree.size() ; i++)
	{
		if(!tree[i]->isEnrichmentFeature)
		{
			for(size_t j  = i+1 ; j < tree.size() ; j++)
			{
				if(!this->tree[j]->isEnrichmentFeature )
				{
					std::vector<Point> inter = tree[i]->intersection(tree[j]) ;
					for(size_t k = 0 ;  k < inter.size() ; k++)
					{
						Point * going_in = new Point(inter[k]) ;

						if(inRoot(*going_in))
						{
							++count ;
							meshPoints.push_back(std::make_pair(going_in, tree[i])) ;
						}
						else
							delete going_in ;

						if(count%100 == 0)
							std::cerr << "\r adding intersection points... " << count << std::flush ;
					}
				}
			}
		}
	}
	count = 0 ;
	std::cerr << " ...done." << std::endl ;
	
	//let us make sure we have no overlap
	std::deque<std::pair<Point *, Feature *> > ::iterator e = std::unique(meshPoints.begin(), meshPoints.end(), PairPointFeatureEqual());
	meshPoints.erase(e, meshPoints.end()) ;
	
	//shuffle for efficiency
	std::random_shuffle(meshPoints.begin(),meshPoints.end()) ;
	
	if(is2D())
	{		
		meshPoints.push_front(std::make_pair(new Point(bbox[0]), tree[0])) ;
		meshPoints.push_front(std::make_pair(new Point(bbox[2]), tree[0])) ;
		meshPoints.push_front(std::make_pair(new Point(bbox[4]), tree[0])) ;
		meshPoints.push_front(std::make_pair(new Point(bbox[6]), tree[0])) ;
		this->dtree = new DelaunayTree( meshPoints[0].first, meshPoints[1].first, meshPoints[2].first) ;
		this->dtree->insert(meshPoints[3].first) ;
		
		for( std::deque<std::pair<Point *, Feature *> >::iterator i = meshPoints.begin()+4 ; i != this->meshPoints.end(); ++i)
		{
			if( (i - meshPoints.begin())%1000 == 0)
				std::cerr << "\r generating triangles... point " << ++count*1000 << "/" << meshPoints.size() << std::flush ;
			
			
			if(*i->first != bbox[0] &&
			   *i->first != bbox[2] &&
			   *i->first != bbox[4] &&
			   *i->first != bbox[6] && inRoot(*i->first)
			  )
				dtree->insert(i->first) ;
		}
		
		std::cerr << "\r generating triangles.... point " << meshPoints.size()-3 << "/" << meshPoints.size()-4 << " ...done" << std::endl ;
		
	}
	else if (is3D())
	{
		meshPoints.push_front(std::make_pair(new Point(bbox[0]), tree[0])) ;
		meshPoints.push_front(std::make_pair(new Point(bbox[1]), tree[0])) ;
		meshPoints.push_front(std::make_pair(new Point(bbox[7]), tree[0])) ;
		meshPoints.push_front(std::make_pair(new Point(bbox[2]), tree[0])) ;
		meshPoints.push_front(std::make_pair(new Point(bbox[3]), tree[0])) ;
		meshPoints.push_front(std::make_pair(new Point(bbox[4]), tree[0])) ;
		meshPoints.push_front(std::make_pair(new Point(bbox[5]), tree[0])) ;
		meshPoints.push_front(std::make_pair(new Point(bbox[6]), tree[0])) ;
		
		
		this->dtree3D = new DelaunayTree_3D( meshPoints[0].first, meshPoints[1].first, meshPoints[2].first, meshPoints[3].first) ;

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
			if((i - meshPoints.begin())%100 == 0)
				std::cerr << "\r generating tetrahedrons... point " << ++count*100 << "/" << meshPoints.size()-8 << std::flush ;
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
					delete i->first ;
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
		
	}
}

std::vector<DelaunayTriangle *> FeatureTree::getTriangles()
{
	if(dtree == NULL)
		this->generateElements() ;
	
	if(dtree != NULL)
	{
		if(!stitched)
			stitch() ;
		
		if(!renumbered)
			renumber() ;

		return dtree->getTriangles() ;
	}
	else
		return std::vector<DelaunayTriangle *>(0) ;
}

std::vector<DelaunayTetrahedron *> FeatureTree::getTetrahedrons()
{
	if(dtree3D == NULL)
		this->generateElements() ;
	
	if(dtree3D != NULL)
	{
		
		if(!stitched)
			stitch() ;

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

std::vector<DelaunayTetrahedron *> Feature::getBoundingTetrahedrons( DelaunayTree_3D * dt)
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
		|| in(inc->getCenter());
}

void Voxel::remove(Feature * inc)
{
	features.erase(std::find(features.begin(), features.end(), inc)) ;
	filled = false ;
}

bool Voxel::add(Feature * inc)
{
	if(filled)
		return false;
	
	if(!features.empty())
	{
		for(size_t i = 0 ; i < this->features.size() ; i++)
		{
			if(
				this->features[i]->intersects(inc)
				)
				return false;
		}
		this->features.push_back(inc) ;
		
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
	this->features.push_back(inc) ;
}

void Voxel::print() const
{
	for(size_t i = 0 ; i < features.size() ; i++)
		features[i]->print() ;
}


Pixel::Pixel()
{
	filled = false ;
}

Pixel::Pixel(double x, double y, double s) : tl(x-s*.5, y+s*.5), tr(x+s*.5, y+s*.5), bl(x-s*.5, y-s*.5), br(x+s*.5, y-s*.5), filled(false) { } ;

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
// 	double size = (tr.x-bl.x)*.1;
// 		std::cout << (p.x >= tl.x-size)  << (p.x <= br.x+size) << (p.y >= br.y-size) << (p.y <= tl.y+size) << std::endl ;
	return (p.x >= tl.x/*-size*/)  && (p.x <= br.x/*+size*/) && (p.y >= br.y/*-size*/) && (p.y <= tl.y/*+size*/);
}

bool Pixel::coOccur(const Geometry * inc) const
{
	
	return inc->in(tl) 
		|| inc->in(tr) 
		|| inc->in(br) 
		|| inc->in(bl) 
		|| in(inc->getCenter()+Point(inc->getRadius(), 0))
		|| in(inc->getCenter())
		|| in(inc->getCenter()+Point(-inc->getRadius(), 0)) 
		|| in(inc->getCenter()+Point(0,inc->getRadius())) 
		|| in(inc->getCenter()+Point(0,-inc->getRadius())) ;
}

void Pixel::remove(Feature * inc)
{
	features.erase(std::find(features.begin(), features.end(), inc)) ;
	filled = false ;
}

bool Pixel::add(Feature * inc)
{
	if(filled)
		return false;

	if(!features.empty())
	{
		for(size_t i = 0 ; i < this->features.size() ; i++)
		{
			if(
			    this->features[i]->intersects(inc)
// 			    || inc->in(inclusions[i]->getCenter())
// 			    || inclusions[i]->in(inc->getCenter())
			  )
				return false;
		}
		this->features.push_back(inc) ;
		
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
}

void Pixel::print() const
{
	for(size_t i = 0 ; i < features.size() ; i++)
		features[i]->print() ;
}


Grid3D::Grid3D(double sizeX, double sizeY, double sizeZ, int div ): x(sizeX), y(sizeY) , z(sizeZ) 
{
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
				pixels[i][j][k] = new Voxel(x*(double)(i)/(double)lengthX-x*.5+psize*.5, y*(double)(j)/(double)lengthY-y*.5+psize*.5,z*(double)(k)/(double)lengthZ-z*.5+psize*.5,  psize) ;
			}
		}
	}
}

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

bool Grid3D::add(Feature * inc)
{
	bool ret = true ;
	std::vector<Voxel *> cleanup ;
	
	double startX = .5*x + inc->getCenter().x-inc->getRadius() ;
	int startI = std::max(0., startX/psize - 2) ;
	
	double endX =  startX+2.*inc->getRadius();
	int endI = std::min(endX/psize + 2, (double)lengthX);
	
	double startY = .5*y + inc->getCenter().y-inc->getRadius() ;
	int startJ = std::max(0., startY/psize - 2) ;
	
	double endY =  startY+2.*inc->getRadius();
	int endJ = std::min(endY/psize + 2, (double)lengthY);
	
	double startZ = .5*y + inc->getCenter().z-inc->getRadius() ;
	int startK = std::max(0., startZ/psize - 2) ;
	
	double endZ =  startZ+2.*inc->getRadius();
	int endK = std::min(endZ/psize + 2, (double)lengthZ);
	
	for(int i = startI ; i < endI ; i++)
	{
		for(int j = startJ ; j < endJ ; j++)
		{
			for(int k = startK ; k < endK ; k++)
			{

				if(pixels[i][j][k]->coOccur(inc))
				{
					if(pixels[i][j][k]->add(inc))
					{
						cleanup.push_back(pixels[i][j][k]) ;
					}
					else
					{
						for(size_t l = 0 ; l < cleanup.size() ; l++)
						{
							cleanup[l]->remove(inc) ;
						}
						return false ;
					}
					
				}
			}
		}
		
	}
	
	return ret ;
}

void Grid3D::forceAdd(Feature * inc)
{
	double startX = .5*x + inc->getCenter().x-inc->getRadius() ;
	int startI = std::max(0., startX/psize - 2) ;
	
	double endX =  startX+2.*inc->getRadius();
	int endI = std::min(endX/psize + 2, (double)lengthX);
	
	double startY = .5*y + inc->getCenter().y-inc->getRadius() ;
	int startJ = std::max(0., startY/psize - 2) ;
	
	double endY =  startY+2.*inc->getRadius();
	int endJ = std::min(endY/psize + 2, (double)lengthY);
	
	double startZ = .5*y + inc->getCenter().z-inc->getRadius() ;
	int startK = std::max(0., startZ/psize - 2) ;
	
	double endZ =  startZ+2.*inc->getRadius();
	int endK = std::min(endZ/psize + 2, (double)lengthZ);
	
	for(int i = startI ; i < endI ; i++)
	{
		for(int j = startJ ; j < endJ ; j++)
		{
			for(int k = startK ; k < endK ; k++)
			{
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
	double startX = .5*x + geo->getCenter().x-geo->getRadius() ;
	int startI = std::max(0., startX/psize - 2) ;
	
	double endX =  startX+2.*geo->getRadius();
	int endI = std::min(endX/psize + 2, (double)lengthX);
	
	double startY = .5*y + geo->getCenter().y-geo->getRadius() ;
	int startJ = std::max(0., startY/psize - 2) ;
	
	double endY =  startY+2.*geo->getRadius();
	int endJ = std::min(endY/psize + 2, (double)lengthY);
	
	double startZ = .5*z + geo->getCenter().z-geo->getRadius() ;
	int startK = std::max(0., startZ/psize - 2) ;
	
	double endZ =  startZ+2.*geo->getRadius();
	int endK = std::min(endZ/psize + 2, (double)lengthZ);
	
	for(int i = startI ; i < endI ; i++)
	{
		for(int j = startJ ; j < endJ ; j++)
		{
			for(int k = startK ; k < endK ; k++)
			{
				if(pixels[i][j][k]->coOccur(geo))
				{
					std::vector<Feature *> to_insert = pixels[i][j][k]->getFeatures() ;
					ret.insert(ret.end(), to_insert.begin(), to_insert.end()) ;
				}
			}
		}
	}
	
	std::stable_sort(ret.begin(), ret.end());
	std::vector<Feature *>::iterator e = std::unique(ret.begin(), ret.end()) ;
	ret.erase(e, ret.end()) ;
	return ret ;
}

const std::vector<Feature *> & Grid3D::coOccur(const Point & p) const
{
	double startX = .5*x + p.x ;
	int startI = std::max(0.,std::min(floor((startX/x)*lengthX), (double)lengthX-1)) ;
	
	double startY = .5*y + p.y ;
	int startJ = std::max(0.,std::min(floor((startY/y)*lengthY), (double)lengthY-1)) ;
	
	double startZ = .5*z + p.z ;
	int startK = std::max(0.,std::min(floor((startZ/z)*lengthZ), (double)lengthZ-1)) ;
	std::cout << "plouf" << std::endl ;
	return pixels[startI][startJ][startK]->getFeatures() ;
}

std::vector<Feature *> & Grid3D::coOccur(const Point & p)
{
	double startX = .5*x + p.x ;
	int startI = std::max(0.,std::min(floor((startX/x)*lengthX), (double)lengthX-1)) ;
	
	double startY = .5*y + p.y ;
	int startJ = std::max(0.,std::min(floor((startY/y)*lengthY), (double)lengthY-1)) ;
	
	double startZ = .5*z + p.z ;
	int startK = std::max(0.,std::min(floor((startZ/z)*lengthZ), (double)lengthZ-1)) ;
	
	return pixels[startI][startJ][startK]->getFeatures() ;
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

Grid::Grid(double sizeX, double sizeY, int div ) : x(sizeX), y(sizeY) 
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
	
	psize = x/lengthX;
	for(size_t i = 0 ; i < lengthX ; i++)
	{
		for(size_t j = 0 ; j < lengthY ; j++)
		{
			pixels[i][j] = new Pixel(x*(double)(i)/(double)lengthX-x*.5+psize*.5, y*(double)(j)/(double)lengthY-y*.5+psize*.5, psize) ;
		}
	}
}

std::vector<Feature *> Grid::coOccur(const Geometry * geo) const
{
	std::vector<Feature *> ret ;
	double startX = .5*x + geo->getCenter().x-geo->getRadius() ;
	int startI = std::max(0., startX/psize - 2) ;
	
	double endX =  startX+2.*geo->getRadius();
	int endI = std::min(endX/psize + 2, (double)lengthX);
	
	double startY = .5*y + geo->getCenter().y-geo->getRadius() ;
	int startJ = std::max(0., startY/psize - 2) ;
	
	double endY =  startY+2.*geo->getRadius();
	int endJ = std::min(endY/psize + 2, (double)lengthY);
	
	for(int i = startI ; i < endI ; i++)
	{
		for(int j = startJ ; j < endJ ; j++)
		{
			if(pixels[i][j]->coOccur(geo))
			{
				std::vector<Feature *> to_insert = pixels[i][j]->getFeatures() ;
				ret.insert(ret.end(), to_insert.begin(), to_insert.end()) ;
			}
		}
	}

	std::stable_sort(ret.begin(), ret.end());
	std::vector<Feature *>::iterator e = std::unique(ret.begin(), ret.end()) ;
	ret.erase(e, ret.end()) ;
	return ret ;
}

 std::vector<Feature *> & Grid::coOccur(const Point & p)
{
	double startX = .5*x + p.x ;
	int startI = std::max(0.,std::min(floor((startX/x)*lengthX), (double)lengthX-1)) ;

	double startY = .5*y + p.y ;
	int startJ = std::max(0.,std::min(floor((startY/y)*lengthY), (double)lengthY-1)) ;
	return pixels[startI][startJ]->getFeatures() ;
}

const std::vector<Feature *> & Grid::coOccur(const Point & p) const
{
	double startX = .5*x + p.x ;
	int startI = std::max(0.,floor((startX/x)*lengthX-1)) ;
	
	double startY = .5*y + p.y ;
	int startJ = std::max(0.,floor((startY/y)*lengthY-1)) ;
	
	return pixels[startI][startJ]->getFeatures() ;
}

void Grid::forceAdd(Feature * inc)
{
	double startX = .5*x + inc->getCenter().x-inc->getRadius() ;
	int startI = std::max(0., startX/psize - 2) ;
	
	double endX =  startX+2.*inc->getRadius();
	int endI = std::min(endX/psize + 2, (double)lengthX);
	
	double startY = .5*y + inc->getCenter().y-inc->getRadius() ;
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

bool Grid::add(Feature * inc)
{
	bool ret = true ;
	std::vector<Pixel *> cleanup ;
	
	double startX = .5*x + inc->getCenter().x-inc->getRadius() ;
	int startI = std::max(0., startX/psize - 2) ;
	
	double endX =  startX+2.*inc->getRadius();
	int endI = std::min(endX/psize + 2, (double)lengthX);
	
	double startY = .5*y + inc->getCenter().y-inc->getRadius() ;
	int startJ = std::max(0., startY/psize - 2) ;
	
	double endY =  startY+2.*inc->getRadius();
	int endJ = std::min(endY/psize + 2, (double)lengthY);
	
	for(int i = startI ; i < endI ; i++)
	{
		for(int j = startJ ; j < endJ ; j++)
		{
			if(pixels[i][j]->coOccur(inc))
			{
				if(pixels[i][j]->add(inc))
				{
					cleanup.push_back(pixels[i][j]) ;
				}
				else
				{
					for(size_t k = 0 ; k < cleanup.size() ; k++)
					{
						cleanup[k]->remove(inc) ;
					}
					return false ;
				}
				
			}
		}
		
	}
	
	return ret ;
}
