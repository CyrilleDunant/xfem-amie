

// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2009-2011
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
#include "../physics/kelvinvoight.h"
#ifdef HAVE_OPENMP
#include <omp.h>
#endif
#include "../physics/homogeneised_behaviour.h"
#include "../solvers/multigrid.h"
#include "../solvers/multigridstep.h"
#include <time.h>
#include <sys/time.h>



using namespace Amie ;



Mesh<DelaunayTriangle, DelaunayTreeItem> * FeatureTree::get2DMesh( int g )
{
	state.setStateTo( RENUMBERED, false ) ;

	if( g == -1 )
		return dtree ;
	else
		return layer2d[g] ;
}

Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * FeatureTree::get3DMesh( int g )
{
	state.setStateTo( RENUMBERED, false ) ;

	if( g == -1 )
		return dtree3D ;
	else
		return layer3d[g] ;
}

std::vector<DelaunayTriangle *> FeatureTree::getBoundingTriangles( const Feature *f )
{
	state.setStateTo( MESHED, false ) ;

	if( f == nullptr )
		return tree[0]->getBoundingElements2D( this ) ;
	else
		return f->getBoundingElements2D( this ) ;
}

std::vector<DelaunayTriangle *> FeatureTree::getElements2D( int g)
{
	state.setStateTo( MESHED, false ) ;

	if( is2D() )
	{
		if( g == -1 )
			return dtree->getElements() ;
		
		return layer2d[g]->getElements() ;
	}

	return std::vector<DelaunayTriangle *>() ;
}

std::vector<DelaunayTriangle *> FeatureTree::getElements2DInLayer( int l )
{
// 	if(l == -1)
// 		return getElements2D() ;
	
	state.setStateTo( MESHED, false ) ;


	if( is2D() && layer2d.find(l) != layer2d.end())
	{
		return layer2d[l]->getElements() ;
	}
	return std::vector<DelaunayTriangle *>() ;
}

std::vector<DelaunayTriangle *> FeatureTree::getActiveElements2D()
{
	std::vector<DelaunayTriangle *> elements  ;
	for(auto j = layer2d.begin() ; j != layer2d.end() ;j++)
	{
		std::vector<DelaunayTriangle *> elementstmp = j->second->getElements() ;
		for( size_t i = 0 ; i < elementstmp.size() ; i++ )
		{
			if(elementstmp[i]->getBehaviour())
				elements.push_back(elementstmp[i]);
		}
	}
	return elements ;
}

std::vector<DelaunayTetrahedron *> FeatureTree::getActiveElements3D()
{
	std::vector<DelaunayTetrahedron *> elements  ;
	for(auto j = layer3d.begin() ; j != layer3d.end() ;j++)
	{
		std::vector<DelaunayTetrahedron *> elementstmp = j->second->getElements() ;
		for( size_t i = 0 ; i < elementstmp.size() ; i++ )
		{
			if(elementstmp[i]->getBehaviour())
				elements.push_back(elementstmp[i]);
		}
	}
	return elements ;
}

std::vector<DelaunayTetrahedron *> FeatureTree::getElements3DInLayer( int l )
{
	if(l == -1)
		return getElements3D() ;
	
	state.setStateTo( MESHED, false ) ;

	if( is3D() && layer3d.find(l) != layer3d.end())
	{
			return layer3d[l]->getElements() ;
	}

	return std::vector<DelaunayTetrahedron *>() ;
}

std::vector<DelaunayTriangle *> FeatureTree::getElements2D( const Point *p, int g )
{
	state.setStateTo( MESHED, false ) ;

	if( is2D() )
	{
		if( g == -1 )
			return dtree->getConflictingElements( p ) ;

		return layer2d[g]->getConflictingElements( p ) ;

	}

	return std::vector<DelaunayTriangle *>() ;
}

std::vector<DelaunayTriangle *> FeatureTree::getElements2D( const Geometry *p, int g )
{
	state.setStateTo( MESHED, false ) ;

	if( is2D() )
	{
		if( g == -1 )
			return dtree->getConflictingElements( p ) ;
		return layer2d[g]->getConflictingElements( p ) ;

	}

	return std::vector<DelaunayTriangle *>() ;
}

std::vector<DelaunayTetrahedron *> FeatureTree::getElements3D( int g )
{
	state.setStateTo( MESHED, false ) ;

	if( is3D() )
	{
		if( g == -1 )
			return dtree3D->getElements() ;
		
		return layer3d[g]->getElements() ;
	}

	return std::vector<DelaunayTetrahedron *>() ;
}

std::vector<DelaunayTetrahedron *> FeatureTree::getElements3D( const Point *p, int g )
{
	state.setStateTo( MESHED, false ) ;

	if( is3D() )
	{
		if( g == -1 )
			return dtree3D->getConflictingElements( p ) ;
		
		return layer3d[g]->getConflictingElements( p ) ;

	}

	return std::vector<DelaunayTetrahedron *>() ;
}

std::vector<DelaunayTetrahedron *> FeatureTree::getElements3D( const Geometry *p, int g )
{
	state.setStateTo( MESHED, false ) ;

	if( is3D() )
	{
		if( g == -1 )
			return dtree3D->getConflictingElements( p ) ;

		return layer3d[g]->getConflictingElements( p ) ;

	}

	return std::vector<DelaunayTetrahedron *>() ;
}

FeatureTree::FeatureTree( Feature *first, int layer, double fraction, size_t gridsize ) : grid( nullptr ), grid3d( nullptr ), state( this ), nodes(0) 
{
	initialValue = 0 ;
	deltaTime = 0 ;
	minDeltaTime = 0.001 ; 
	reuseDisplacements = false ;
	foundCheckPoint = true ;
	averageDamage = 0 ;
	damageConverged = false ;
	stateConverged = false ;
	dtree = nullptr ;
	dtree3D = nullptr ;
	
	samplingRestriction = SAMPLE_NO_RESTRICTION ;

	if( first )
	{
		addFeature( nullptr, first, layer, fraction ) ;
	}

	if( is2D() )
		grid = new Grid( ( first->getBoundingBox()[1].getX() - first->getBoundingBox()[0].getX() ) * 1.1,
		                 ( first->getBoundingBox()[1].getY() - first->getBoundingBox()[2].getY() ) * 1.1, gridsize,
		                 Point( ( first->getBoundingBox()[1].getX() + first->getBoundingBox()[0].getX() )*.5,
		                        ( first->getBoundingBox()[1].getY() + first->getBoundingBox()[2].getY() )*.5
		                      ) ) ;

	if( is3D() )
		grid3d = new Grid3D( ( first->getBoundingBox()[7].getX() - first->getBoundingBox()[0].getX() ) * 1.1,
		                     ( first->getBoundingBox()[7].getY() - first->getBoundingBox()[0].getY() ) * 1.1,
		                     ( first->getBoundingBox()[7].getZ() - first->getBoundingBox()[0].getZ() ) * 1.1, gridsize / 5, ( first->getBoundingBox()[7] + first->getBoundingBox()[0] )*.5 );

	father3D = nullptr;
	father2D = nullptr ;
	elemOrder = LINEAR ;
	renumbered = false ;
	needAssembly = true ;
	setBehaviours = false ;
	behaviourChange = true ;
	solverConvergence = false ;
	enrichmentChange = true ;
	needMeshing = true ;

	elastic = false ;
	projectOnBoundaries = true ;

	K = new Assembly() ;

	if( is2D() )
		K->setSpaceDimension( SPACE_TWO_DIMENSIONAL ) ;
	else
		K->setSpaceDimension( SPACE_THREE_DIMENSIONAL ) ;

	crackedVolume = 0 ;
	damagedVolume = 0 ;
	residualError = 1e9 ;
	samplingNumber = 0 ;
	previousSamplingNumber = 0 ;


	lastNodeId = 0;
	lastEnrichmentId = 0;
	maxitPerStep = 200 ;
	deltaTime = .1 ;
	realDeltaTime = deltaTime ;
	now = 0 ;

	setElementGenerationMethod() ;

}

void FeatureTree::twineFeature( CompositeFeature *father, CompositeFeature *f )
{
	std::vector<Feature *> chain = father->getDescendants() ;
	std::vector<VirtualFeature *> fatherComponents = father->getComponents() ;
	std::vector<VirtualFeature *> childComponents = f->getComponents() ;

	//we look for the father of the second descendant of the "father" feature
	Feature *headerEnd = ( *std::find( chain.begin(), chain.end(), fatherComponents[0] ) ) ;
	Feature *headerEndChild = ( *headerEnd->getChildren().rbegin() ) ;

	//we attach the first component of the twinee to the twinee
	addFeature( f, childComponents[0] ) ;

	//the header end loses its child and gets a new one
	headerEnd->removeChild( *headerEnd->getChildren().rbegin() );
	addFeature( headerEnd, f ) ;

	//we re-attach the new header to the rest
	childComponents[0]->addChild( headerEndChild ) ;
	headerEndChild->setFather( childComponents[0] ) ;

	//now that the header is complete, we find each father component in the chain, and insert at that
	//point a new link: the corresponding child component

	for( size_t i = 1 ; i < fatherComponents.size() ; i++ )
	{
		Feature *attachPoint = ( *std::find( chain.begin(), chain.end(), fatherComponents[i] ) ) ;
		Feature *attachPointChild = nullptr;

		if( !attachPoint->getChildren().empty() )
			attachPointChild = ( *attachPoint->getChildren().rbegin() ) ;

		//we detach the attach point from its father
		attachPoint->removeChild( attachPointChild ) ;

		//we attach the corresponding child component to the loose end
		addFeature( attachPoint, childComponents[i] ) ;

		//we re-attach the attach point to the new feature
		if( attachPointChild )
		{
			childComponents[i]->addChild( attachPointChild ) ;
			attachPointChild->setFather( childComponents[i] ) ;

		}

	}
}

void FeatureTree::addPoint( Point *p )
{
	extraPoints.push_back(p);
}

void FeatureTree::addFeature( Feature *father, Feature *f, int layer, double fraction )
{
	f->setLayer(layer) ;
	f->setFraction(fraction) ;
	scalingFactors[layer] = fraction ;
	if( !f->isEnrichmentFeature )
		needMeshing = true ;

	if( !tree.empty() && f->spaceDimensions() == SPACE_TWO_DIMENSIONAL && !f->isEnrichmentFeature )
		grid->forceAdd( f ) ;
	else if( !tree.empty() && !f->isEnrichmentFeature )
		grid3d->forceAdd( f ) ;

	if( f->isCompositeFeature && father && !father->isCompositeFeature )
	{
		std::vector<VirtualFeature *> pile = dynamic_cast<CompositeFeature *>( f )->getComponents();
		f->setFather( father ) ;

		if( father != nullptr )
			father->addChild( f ) ;

		this->tree.push_back( f ) ;
		addFeature( f, pile[0] ) ;

		for( size_t i = 0 ; i < pile.size() - 1 ; i++ )
		{
			addFeature( pile[i], pile[i + 1] ) ;
		}

		return ;

	}

	f->setFather( father ) ;

	
	if( father != nullptr )
		father->addChild( f ) ;

	this->tree.push_back( f ) ;
	

}

FeatureTree::~FeatureTree()
{
	delete father3D ;
	delete father2D ;
	delete grid ;
	delete grid3d ;
	for(auto j = layer2d.begin() ; j!=layer2d.end() ; ++j)
	{
		delete j->second ;
	}
	for(auto j = layer3d.begin() ; j!=layer3d.end() ; ++j)
	{
		delete j->second ;
	}
	delete this->dtree3D ;
	delete this->K ;

	for( size_t i = 0 ; i < additionalPoints.size() ; i++ )
		delete additionalPoints[i] ;

	for( size_t i = 0 ; i < boundaryCondition.size() ; ++i )
		delete boundaryCondition[i] ;
	for(size_t i = 0 ; i < extraPoints.size() ; ++i)
		delete extraPoints[i] ;
	
}

void FeatureTree::scaleBoundaryConditions(double scale)
{
	for(size_t i = 0 ; i < boundaryCondition.size() ; i++)
	{
		boundaryCondition[i]->setScale(scale);
	}
}

void FeatureTree::addBoundaryCondition( BoundaryCondition *bc )
{
	boundaryCondition.push_back( bc ) ;
}

void FeatureTree::removeBoundaryCondition( BoundaryCondition *bc )
{
	std::vector<BoundaryCondition *>::iterator toDelete = std::find( boundaryCondition.begin(), boundaryCondition.end(), bc ) ;
	boundaryCondition.erase( toDelete ) ;
}

void FeatureTree::setOrder( Order ord )
{
	state.stitched = false ;
	state.renumbered = false ;
	state.initialised = false ;

	elemOrder = ord ;

	if( father3D )
		delete father3D ;

	father3D = new TetrahedralElement( elemOrder ) ;
	father3D->compileAndPrecalculate() ;


	if( father2D )
		delete father2D ;

	father2D = new TriElement( elemOrder ) ;
	father2D->compileAndPrecalculate() ;
	
	if(ord >= CONSTANT_TIME_LINEAR)
	{
		addBoundaryCondition( new TimeContinuityBoundaryCondition(initialValue) ) ;
	}
	
}

void FeatureTree::renumber()
{
	if( is2D() )
	{
			std::vector<DelaunayTriangle *> triangles = dtree->getElements() ;
			size_t count = 0 ;
			std::cerr << " renumbering... " << std::flush ;

			for( auto i = triangles.begin() ; i != triangles.end() ; ++i )
			{
				for( size_t j = 0 ; j < ( *i )->getBoundingPoints().size() ; j++ )
				{
					( *i )->getBoundingPoint( j ).setId( -1 ) ;
				}
			}

			std::vector<DelaunayTriangle *> sortedElements ;
			sortedElements = triangles ;

			for( auto i = sortedElements.begin() ; i != sortedElements.end() ; ++i )
			{
				for( size_t j = 0 ; j < (*i)->getBoundingPoints().size()/(*i)->timePlanes() ; j++ )
				{
					if( (*i)->getBoundingPoint( j ).getId() == -1 )
						(*i)->getBoundingPoint( j ).setId(count++) ;
				}
			}
			

			lastNodeId = count ;

			for( auto i = sortedElements.begin() ; i != sortedElements.end() ; ++i )
			{
				for(size_t k = 1 ; k < (*i)->timePlanes() ; k++)
				{
					for( size_t j = 0 ; j < (*i)->getBoundingPoints().size()/(*i)->timePlanes() ; j++ )
					{
						if( (*i)->getBoundingPoint( j + k*(*i)->getBoundingPoints().size()/(*i)->timePlanes() ).getId() == -1 )
						{
							(*i)->getBoundingPoint( j + k*(*i)->getBoundingPoints().size()/(*i)->timePlanes()).setId( (*i)->getBoundingPoint( j).getId() + lastNodeId*k) ;
							count++ ;
						}
					}
				}
			}
			
			lastNodeId = count ;
			
			nodes.resize(count) ;
			for( auto i = sortedElements.begin() ; i != sortedElements.end() ; ++i )
			{
				
				if(*i)
				{
					for(size_t j = 0 ; j < (*i)->getBoundingPoints().size() ; j++)
					{
						nodes[ (*i)->getBoundingPoint( j ).getId() ] = &(*i)->getBoundingPoint( j ) ;
					}
				}
				else
				{
					std::cerr << "nullTri" << std::endl ;
				}

			}

			std::cerr << count * 2 << " ...done " << std::endl ;

	}
	else if( is3D() )
	{
		std::vector<DelaunayTetrahedron *> tets = dtree3D->getElements() ;
		size_t count = 0 ;
		std::cerr << " renumbering... " << std::flush ;

		for( auto i = tets.begin() ; i != tets.end() ; ++i )
		{
			for( size_t j = 0 ; j < ( *i )->getBoundingPoints().size() ; j++ )
			{
				( *i )->getBoundingPoint( j ).setId( -1 ) ;
			}
		}

		std::vector< DelaunayTetrahedron *> sortedElements ;
		std::set<const Geometry *> placedElements ;

		sortedElements = tets ;

		for( auto i = sortedElements.begin() ; i != sortedElements.end() ; ++i )
		{
			DelaunayTetrahedron *tet = const_cast<DelaunayTetrahedron *>( *i ) ;

			for( size_t j = 0 ; j < tet->getBoundingPoints().size()/tet->timePlanes() ; j++ )
			{
				if( tet->getBoundingPoint( j ).getId() == -1 )
					tet->getBoundingPoint( j ).getId() = count++ ;
			}
		}

		lastNodeId = count ;

		for( auto i = sortedElements.begin() ; i != sortedElements.end() ; ++i )
		{
			for(size_t k = 1 ; k < (*i)->timePlanes() ; k++)
			{
				for( size_t j = 0 ; j < (*i)->getBoundingPoints().size()/(*i)->timePlanes() ; j++ )
				{
					if( (*i)->getBoundingPoint( j + k*(*i)->getBoundingPoints().size()/(*i)->timePlanes() ).getId() == -1 )
					{
						(*i)->getBoundingPoint( j + k*(*i)->getBoundingPoints().size()/(*i)->timePlanes()).setId( (*i)->getBoundingPoint( j).getId() + lastNodeId*k );
						count++ ;
					}
				}
			}
		}
		
		lastNodeId = count ;
		
		nodes.resize(count) ;
		for( auto i = sortedElements.begin() ; i != sortedElements.end() ; ++i )
		{
			
			if(*i)
			{
				for(size_t j = 0 ; j < (*i)->getBoundingPoints().size() ; j++)
				{
					nodes[ (*i)->getBoundingPoint( j ).getId() ] = &(*i)->getBoundingPoint( j ) ;
				}
			}
			else
			{
				std::cerr << "nullTri" << std::endl ;
			}
		}

		std::cerr << count * 3 << " ...done " << std::endl ;


	}

	renumbered = true ;

}

bool FeatureTree::inRoot( const Point &p ) const
{
	if( is2D() )
	{
		Point p0( p.getX(), p.getY() + POINT_TOLERANCE_2D ) ;
		Point p1( p.getX(), p.getY() - POINT_TOLERANCE_2D ) ;
		Point p2( p.getX() + POINT_TOLERANCE_2D, p.getY() ) ;
		Point p3( p.getX() - POINT_TOLERANCE_2D, p.getY() ) ;
		return ( tree[0]->in( p ) || tree[0]->in( p0 ) || tree[0]->in( p1 ) || tree[0]->in( p2 ) || tree[0]->in( p3 ) ) ;
	}
	else
	{
		Point p0( p.getX(), p.getY() + POINT_TOLERANCE_3D, p.getZ() ) ;
		Point p1( p.getX(), p.getY() - POINT_TOLERANCE_3D, p.getZ() ) ;
		Point p2( p.getX() + POINT_TOLERANCE_3D, p.getY(), p.getZ() ) ;
		Point p3( p.getX() - POINT_TOLERANCE_3D, p.getY(), p.getZ() ) ;
		Point p4( p.getX(), p.getY(), p.getZ() + POINT_TOLERANCE_3D ) ;
		Point p5( p.getX(), p.getY(), p.getZ() - POINT_TOLERANCE_3D ) ;
		return ( tree[0]->in( p ) || tree[0]->in( p0 ) || tree[0]->in( p1 ) || tree[0]->in( p2 ) || tree[0]->in( p3 ) || tree[0]->in( p4 ) || tree[0]->in( p5 ) ) ;
	}
}

void FeatureTree::projectTetrahedronsOnBoundaries( size_t edge, size_t time )
{
	if( edge + time == 0 )
		return ;

	size_t first = 0 ;
	size_t second = ( edge + 1 ) ;
	size_t third = ( edge + 1 ) * 2  ;
	size_t fourth = ( edge + 1 ) * 3  ;

	std::vector<size_t> indexes( edge * ( time + 1 ) * 6 ) ;



	size_t count = 0 ;

	
	Point a( 0.25, 0.25, 0.25 ) ;
	Point b( 0.166666666666667, 0.166666666666667, 0.166666666666667 ) ;
	Point c( 0.5, 0.166666666666667, 0.166666666666667 ) ;
	Point d( 0.166666666666667, 0.5, 0.166666666666667 ) ;
	Point e( 0.166666666666667, 0.166666666666667, 0.5 ) ;

	for( size_t j = 1 ; j < this->tree.size() ; j++ )
	{
		if( !tree[j]->isEnrichmentFeature )
		{
			//In two pass
		  std::vector<DelaunayTetrahedron *> tets = this->tree[j]->getElements3D( this ) ;
		  std::valarray<Point> originalPoints( tets[0]->getBoundingPoints().size() ) ;

		 for( size_t i = 0 ; i < tets.size() ; i++ )
		 {

			tets[i]->refresh( father3D ) ;
			Point proj_0( *tets[i]->first ) ;
			tree[j]->project( &proj_0 ) ;
			Point proj_1( *tets[i]->second ) ;
			tree[j]->project( &proj_1 ) ;
			Point proj_2( *tets[i]->third ) ;
			tree[j]->project( &proj_2 ) ;
			Point proj_3( *tets[i]->fourth ) ;
			tree[j]->project( &proj_3 ) ;

			if(
			    squareDist3D( proj_0 , *tets[i]->first  ) < POINT_TOLERANCE_3D * POINT_TOLERANCE_3D &&
			    squareDist3D( proj_1 , *tets[i]->second ) < POINT_TOLERANCE_3D * POINT_TOLERANCE_3D
			)
			{
				  count++;
				  indexes.clear() ;
				  for(size_t k = 0 ; k < tets[i]->getBoundingPoints().size() ; k++)
				  {
					  if(tets[i]->getBoundingPoint(k) == Point(0.5*(tets[i]->first->getX()+tets[i]->second->getX()),0.5*(tets[i]->first->getY()+tets[i]->second->getY()),0.5*(tets[i]->first->getZ()+tets[i]->second->getZ()),tets[i]->getBoundingPoint(k).getT()))
					  {
					    indexes.push_back(k) ;
					  }
				  }
				  Point test = tets[i]->getBoundingPoint( indexes[0] ) ;
				  tree[j]->project( &test ) ;

				  if( inRoot( test ) )
				  {
					  for( size_t ni = 0 ; ni < indexes.size() ; ni++ )
					  {
						  size_t k = indexes[ni] ;
						  originalPoints[ni] = tets[i]->getBoundingPoint( k ) ;
						  tree[j]->project( &tets[i]->getBoundingPoint( k ) ) ;
					  }

					  if( tets[i]->jacobianAtPoint( a ) > 0 &&
						  tets[i]->jacobianAtPoint( b ) > 0 &&
						  tets[i]->jacobianAtPoint( c ) > 0 &&
						  tets[i]->jacobianAtPoint( d ) > 0 &&
						  tets[i]->jacobianAtPoint( e ) > 0
					    )
					  {
						  tets[i]->moved = true ;

						  for( size_t j = 0 ; j < 4 ; j++ )
						  {
							  if( tets[i]->getNeighbour( j )->isTetrahedron() )
							  {
								  dynamic_cast<DelaunayTetrahedron *>( tets[i]->getNeighbour( j ) )->moved = true ;
							  }
						  }
					  }
					  else
					  {
						  for( size_t ni = 0 ; ni < indexes.size() ; ni++ )
						  {
							  size_t k = indexes[ni] ;
							  tets[i]->getBoundingPoint( k ) = originalPoints[ni] ;
						  }
					  }
				  }
			  }

			  if(
			      squareDist3D( proj_1 , *tets[i]->second ) < POINT_TOLERANCE_3D * POINT_TOLERANCE_3D &&
			      squareDist3D( proj_2 , *tets[i]->third  ) < POINT_TOLERANCE_3D * POINT_TOLERANCE_3D
			  )
			  {
				count++;
				indexes.clear() ;
				for(size_t k = 0 ; k < tets[i]->getBoundingPoints().size() ; k++)
				{
					  if(tets[i]->getBoundingPoint(k) == Point(0.5*(tets[i]->third->getX()+tets[i]->second->getX()),0.5*(tets[i]->third->getY()+tets[i]->second->getY()),0.5*(tets[i]->third->getZ()+tets[i]->second->getZ()),tets[i]->getBoundingPoint(k).getT()))
					  {
					    indexes.push_back(k) ;
					  }
				}
				Point test = tets[i]->getBoundingPoint( indexes[0] ) ;
				tree[j]->project( &test ) ;

				if( inRoot( test ) )
				{
					  for( size_t ni = 0 ; ni < indexes.size() ; ni++ )
					  {
						  size_t k = indexes[ni] ;
						  originalPoints[ni] = tets[i]->getBoundingPoint( k ) ;
						  tree[j]->project( &tets[i]->getBoundingPoint( k ) ) ;
					  }

					  if( tets[i]->jacobianAtPoint( a ) > 0 &&
						  tets[i]->jacobianAtPoint( b ) > 0 &&
						  tets[i]->jacobianAtPoint( c ) > 0 &&
						  tets[i]->jacobianAtPoint( d ) > 0 &&
						  tets[i]->jacobianAtPoint( e ) > 0
					    )
					  {
						  tets[i]->moved = true ;

						  for( size_t j = 0 ; j < 4 ; j++ )
						  {
							  if( tets[i]->getNeighbour( j )->isTetrahedron() )
							  {
								  dynamic_cast<DelaunayTetrahedron *>( tets[i]->getNeighbour( j ) )->moved = true ;
							  }
						  }
					  }
					  else
					  {
						  for( size_t ni = 0 ; ni < indexes.size() ; ni++ )
						  {
							  size_t k = indexes[ni] ;
							  tets[i]->getBoundingPoint( k ) = originalPoints[ni] ;
						  }
					  }
				  }
			  }

			  if(
			      squareDist3D( proj_3 , *tets[i]->fourth ) < POINT_TOLERANCE_3D * POINT_TOLERANCE_3D &&
			      squareDist3D( proj_2 , *tets[i]->third  ) < POINT_TOLERANCE_3D * POINT_TOLERANCE_3D
			  )
			  {
				count++;
				indexes.clear() ;
				for(size_t k = 0 ; k < tets[i]->getBoundingPoints().size() ; k++)
				{
					  if(tets[i]->getBoundingPoint(k) == Point(0.5*(tets[i]->third->getX()+tets[i]->fourth->getX()),0.5*(tets[i]->third->getY()+tets[i]->fourth->getY()),0.5*(tets[i]->third->getZ()+tets[i]->fourth->getZ()),tets[i]->getBoundingPoint(k).getT()))
					  {
					    indexes.push_back(k) ;
					  }
				}
				Point test = tets[i]->getBoundingPoint( indexes[0] ) ;
				tree[j]->project( &test ) ;

				if( inRoot( test ) )
				{
					  for( size_t ni = 0 ; ni < indexes.size() ; ni++ )
					  {
						  size_t k = indexes[ni] ;
						  originalPoints[ni] = tets[i]->getBoundingPoint( k ) ;
						  tree[j]->project( &tets[i]->getBoundingPoint( k ) ) ;
					  }

					  if( tets[i]->jacobianAtPoint( a ) > 0 &&
						  tets[i]->jacobianAtPoint( b ) > 0 &&
						  tets[i]->jacobianAtPoint( c ) > 0 &&
						  tets[i]->jacobianAtPoint( d ) > 0 &&
						  tets[i]->jacobianAtPoint( e ) > 0
					    )
					  {
						  tets[i]->moved = true ;

						for( size_t j = 0 ; j < 4 ; j++ )
						{
							if( tets[i]->getNeighbour( j )->isTetrahedron() )
							{
							  dynamic_cast<DelaunayTetrahedron *>( tets[i]->getNeighbour( j ) )->moved = true ;
							}
						}
					  }
					  else
					  {
						for( size_t ni = 0 ; ni < indexes.size() ; ni++ )
						{
							  size_t k = indexes[ni] ;
							  tets[i]->getBoundingPoint( k ) = originalPoints[ni] ;
						}
					}
				}
			  }

			  if(
			      squareDist3D( proj_0 , *tets[i]->first ) < POINT_TOLERANCE_3D * POINT_TOLERANCE_3D &&
			      squareDist3D( proj_3 , *tets[i]->fourth ) < POINT_TOLERANCE_3D * POINT_TOLERANCE_3D
			  )
			  {
				  count++;
				  indexes.clear() ;
				  for(size_t k = 0 ; k < tets[i]->getBoundingPoints().size() ; k++)
				  {
					  if(tets[i]->getBoundingPoint(k) == Point(0.5*(tets[i]->first->getX()+tets[i]->fourth->getX()),0.5*(tets[i]->first->getY()+tets[i]->fourth->getY()),0.5*(tets[i]->first->getZ()+tets[i]->fourth->getZ()),tets[i]->getBoundingPoint(k).getT()))
					  {
					    indexes.push_back(k) ;
					  }
				  }
				  Point test = tets[i]->getBoundingPoint( indexes[0] ) ;
				  tree[j]->project( &test ) ;

				  if( inRoot( test ) )
				  {
					  for( size_t ni = 0 ; ni < indexes.size() ; ni++ )
					  {
						  size_t k = indexes[ni] ;
						  originalPoints[ni] = tets[i]->getBoundingPoint( k ) ;
						  tree[j]->project( &tets[i]->getBoundingPoint( k ) ) ;
					  }

					  if( tets[i]->jacobianAtPoint( a ) > 0 &&
						  tets[i]->jacobianAtPoint( b ) > 0 &&
						  tets[i]->jacobianAtPoint( c ) > 0 &&
						  tets[i]->jacobianAtPoint( d ) > 0 &&
						  tets[i]->jacobianAtPoint( e ) > 0
					    )
					  {
						  tets[i]->moved = true ;

						  for( size_t j = 0 ; j < 4 ; j++ )
						  {
							  if( tets[i]->getNeighbour( j )->isTetrahedron() )
							  {
								  dynamic_cast<DelaunayTetrahedron *>( tets[i]->getNeighbour( j ) )->moved = true ;
							  }
						  }
					  }
					  else
					  {
						  for( size_t ni = 0 ; ni < indexes.size() ; ni++ )
						  {
							  size_t k = indexes[ni] ;
							  tets[i]->getBoundingPoint( k ) = originalPoints[ni] ;
						  }
					  }
				  }
			  }

			  if(
			      squareDist3D( proj_1 , *tets[i]->second ) < POINT_TOLERANCE_3D * POINT_TOLERANCE_3D &&
			      squareDist3D( proj_3 , *tets[i]->fourth ) < POINT_TOLERANCE_3D * POINT_TOLERANCE_3D
			  )
			  {
				  count++;
				  indexes.clear() ;
				  for(size_t k = 0 ; k < tets[i]->getBoundingPoints().size() ; k++)
				  {
					  if(tets[i]->getBoundingPoint(k) == Point(0.5*(tets[i]->second->getX()+tets[i]->fourth->getX()),0.5*(tets[i]->second->getY()+tets[i]->fourth->getY()),0.5*(tets[i]->second->getZ()+tets[i]->fourth->getZ()),tets[i]->getBoundingPoint(k).getT()))
					  {
					    indexes.push_back(k) ;
					  }
				  }
				  Point test = tets[i]->getBoundingPoint( indexes[0] ) ;
				  tree[j]->project( &test ) ;

				  if( inRoot( test ) )
				  {
					  for( size_t ni = 0 ; ni < indexes.size() ; ni++ )
					  {
						  size_t k = indexes[ni] ;
						  originalPoints[ni] = tets[i]->getBoundingPoint( k ) ;
						  tree[j]->project( &tets[i]->getBoundingPoint( k ) ) ;
					  }

					  if( tets[i]->jacobianAtPoint( a ) > 0 &&
						  tets[i]->jacobianAtPoint( b ) > 0 &&
						  tets[i]->jacobianAtPoint( c ) > 0 &&
						  tets[i]->jacobianAtPoint( d ) > 0 &&
						  tets[i]->jacobianAtPoint( e ) > 0
					    )
					  {
						  tets[i]->moved = true ;

						  for( size_t j = 0 ; j < 4 ; j++ )
						  {
							  if( tets[i]->getNeighbour( j )->isTetrahedron() )
							  {
								  dynamic_cast<DelaunayTetrahedron *>( tets[i]->getNeighbour( j ) )->moved = true ;
							  }
						  }
					  }
					  else
					  {
						  for( size_t ni = 0 ; ni < indexes.size() ; ni++ )
						  {
							  size_t k = indexes[ni] ;
							  tets[i]->getBoundingPoint( k ) = originalPoints[ni] ;
						  }
					  }
				  }
			  }

			  if(
			      squareDist3D( proj_0 , *tets[i]->first ) < POINT_TOLERANCE_3D * POINT_TOLERANCE_3D &&
			      squareDist3D( proj_2 , *tets[i]->third ) < POINT_TOLERANCE_3D * POINT_TOLERANCE_3D
			  )
			  {
				  count++;
				  indexes.clear() ;
				  for(size_t k = 0 ; k < tets[i]->getBoundingPoints().size() ; k++)
				  {
					  if(tets[i]->getBoundingPoint(k) == Point(0.5*(tets[i]->first->getX()+tets[i]->third->getX()),0.5*(tets[i]->first->getY()+tets[i]->third->getY()),0.5*(tets[i]->first->getZ()+tets[i]->third->getZ()),tets[i]->getBoundingPoint(k).getT()))
					  {
					    indexes.push_back(k) ;
					  }
				  }
				  Point test = tets[i]->getBoundingPoint( indexes[0] ) ;
				  tree[j]->project( &test ) ;

				  if( inRoot( test ) )
				  {
					  for( size_t ni = 0 ; ni < indexes.size() ; ni++ )
					  {
						  size_t k = indexes[ni] ;
						  originalPoints[ni] = tets[i]->getBoundingPoint( k ) ;
						  tree[j]->project( &tets[i]->getBoundingPoint( k ) ) ;
					  }

					  if( tets[i]->jacobianAtPoint( a ) > 0 &&
						  tets[i]->jacobianAtPoint( b ) > 0 &&
						  tets[i]->jacobianAtPoint( c ) > 0 &&
						  tets[i]->jacobianAtPoint( d ) > 0 &&
						  tets[i]->jacobianAtPoint( e ) > 0
					    )
					  {
						  tets[i]->moved = true ;

						  for( size_t j = 0 ; j < 4 ; j++ )
						  {
							  if( tets[i]->getNeighbour( j )->isTetrahedron() )
							  {
								  dynamic_cast<DelaunayTetrahedron *>( tets[i]->getNeighbour( j ) )->moved = true ;
							  }
						  }
					  }
					  else
					  {
						  for( size_t ni = 0 ; ni < indexes.size() ; ni++ )
						  {
							  size_t k = indexes[ni] ;
							  tets[i]->getBoundingPoint( k ) = originalPoints[ni] ;
						  }
					  }
				  }
			  }

			  if( count % 1000 == 0 )
				  std::cerr << "\r projecting points on boundaries... point " << count << "/xx" << " feature " << j << std::flush ;
		  }
	    
		  
		}
	}

	std::cerr << "\r projecting points on boundaries... point " << count << "/xx" << " ...done." << std::endl ;

}

void FeatureTree::projectTrianglesOnBoundaries( size_t edge, size_t time )
{
	if( edge + time == 0 )
		return ;

	std::valarray<size_t> indexes( edge * ( time + 1 ) * 3 ) ;

	for( size_t j = 0 ; j < 3 ; j++ )
	{
		for( size_t e = 0 ; e < edge ; e++ )
		{
			indexes[3 * e + j] = j * ( edge + 1 ) + e + 1 ;

			for( size_t t = 0 ; t < time ; t++ )
			{
				indexes[3 * ( edge * ( t + 1 ) + e ) + j] = j * ( edge + 1 ) + e + 1 + 3 * ( edge + 1 ) * ( t + 1 ) ;
			}
		}
	}

	size_t count = 0 ;
	size_t pd = 0 ;
	size_t k = 0 ;
	size_t n = indexes.size() / 3 ;

	std::valarray<Point> originalPoints( n ) ;

	Point a( 0.2, 0.2 ) ;
	Point b( 0.6, 0.2 ) ;
	Point c( 0.2, 0.6 ) ;
	Point d( 1. / 3., 1. / 3. ) ;

	for( size_t j = 1 ; j < this->tree.size() ; j++ )
	{
		if( !tree[j]->isEnrichmentFeature && tree[j]->getGeometryType() != TRIANGLE && tree[j]->getGeometryType() != RECTANGLE )
		{

			std::vector<DelaunayTriangle *> triangles = this->tree[j]->getElements2D( this ) ;

			for( size_t i = 0 ; i < triangles.size() ; i++ )
			{
				triangles[i]->refresh( father2D ) ;

				if( triangles[i]->getPrimitive()->intersects( tree[j] ) )
				{

					Point proj_0( *triangles[i]->first ) ;
					tree[j]->project( &proj_0 ) ;
					Point proj_1( *triangles[i]->second ) ;
					tree[j]->project( &proj_1 ) ;
					Point proj_2( *triangles[i]->third ) ;
					tree[j]->project( &proj_2 ) ;
					bool changed  = true;

					if( squareDist2D( &proj_0 , triangles[i]->first ) < POINT_TOLERANCE_2D &&
					        squareDist2D( &proj_1 , triangles[i]->second ) < POINT_TOLERANCE_2D &&
					        squareDist2D( &proj_2 , triangles[i]->third ) > 10.*POINT_TOLERANCE_2D )
					{
						count += changed ;
						changed = false ;
						Point test = triangles[i]->getBoundingPoint( indexes[0] ) ;
						tree[j]->project( &test ) ;

						if( inRoot( test ) )
						{
							for( size_t ni = 0 ; ni < n ; ni++ )
							{
								k = indexes[3 * ni + 0] ;
								originalPoints[ni] = triangles[i]->getBoundingPoint( k ) ;
								tree[j]->project( &triangles[i]->getBoundingPoint( k ) ) ;
							}

							if( triangles[i]->jacobianAtPoint( a ) > 0 &&
							        triangles[i]->jacobianAtPoint( b ) > 0 &&
							        triangles[i]->jacobianAtPoint( c ) > 0 &&
							        triangles[i]->jacobianAtPoint( d ) > 0
							  )
							{
								triangles[i]->moved = true ;

								for( size_t j = 0 ; j < 3 ; j++ )
								{
									if( triangles[i]->getNeighbour( j )->isTriangle )
									{
										dynamic_cast<DelaunayTriangle *>( triangles[i]->getNeighbour( j ) )->moved = true ;
									}
								}
							}
							else
							{
								for( size_t ni = 0 ; ni < n ; ni++ )
								{
									k = indexes[3 * ni + 0] ;
									triangles[i]->getBoundingPoint( k ) = originalPoints[ni] ;
								}
							}
						}

// 						std::cerr << "--> " << (*triangles)[i]->getBoundingPoint(1)->getX() << ", " << (*triangles)[i]->getBoundingPoint(1)->getY() << std::endl ;
					}

					if( squareDist2D( &proj_1 , triangles[i]->second ) < POINT_TOLERANCE_2D &&
					        squareDist2D( &proj_2 , triangles[i]->third ) < POINT_TOLERANCE_2D &&
					        squareDist2D( &proj_0 , triangles[i]->first ) > 10.*POINT_TOLERANCE_2D
					  )
					{
						count += changed ;
						changed = false ;
						Point test = triangles[i]->getBoundingPoint( indexes[1] ) ;
						tree[j]->project( &test ) ;

						if( inRoot( test ) )
						{
							for( size_t ni = 0 ; ni < n ; ni++ )
							{
								k = indexes[3 * ni + 1] ;
								originalPoints[ni] = triangles[i]->getBoundingPoint( k ) ;
								tree[j]->project( &triangles[i]->getBoundingPoint( k ) ) ;
							}

							if( triangles[i]->jacobianAtPoint( a ) > 0 &&
							        triangles[i]->jacobianAtPoint( b ) > 0 &&
							        triangles[i]->jacobianAtPoint( c ) > 0 &&
							        triangles[i]->jacobianAtPoint( d ) > 0
							  )
							{
								triangles[i]->moved = true ;

								for( size_t j = 0 ; j < 3 ; j++ )
								{
									if( triangles[i]->getNeighbour( j )->isTriangle )
									{
										dynamic_cast<DelaunayTriangle *>( triangles[i]->getNeighbour( j ) )->moved = true ;
									}
								}
							}
							else
							{
								for( size_t ni = 0 ; ni < n ; ni++ )
								{
									k = indexes[3 * ni + 1] ;
									triangles[i]->getBoundingPoint( k ) = originalPoints[ni] ;
								}
							}
						}

// 						std::cerr << "--> " << (*triangles)[i]->getBoundingPoint(3)->getX() << ", " << (*triangles)[i]->getBoundingPoint(3)->getY() << std::endl ;
					}

					if( squareDist2D( &proj_2 , triangles[i]->third ) < POINT_TOLERANCE_2D &&
					        squareDist2D( &proj_0, triangles[i]->first ) < POINT_TOLERANCE_2D &&
					        squareDist2D( &proj_1, triangles[i]->second ) > 10.*POINT_TOLERANCE_2D
					  )
					{
						count += changed ;
						changed = false ;
						Point test = triangles[i]->getBoundingPoint( indexes[2] ) ;
						tree[j]->project( &test ) ;

						if( inRoot( test ) )
						{
							for( size_t ni = 0 ; ni < n ; ni++ )
							{
								k = indexes[3 * ni + 2] ;
								originalPoints[ni] = triangles[i]->getBoundingPoint( k ) ;
								tree[j]->project( &triangles[i]->getBoundingPoint( k ) ) ;
							}

							if( triangles[i]->jacobianAtPoint( a ) > 0 &&
							        triangles[i]->jacobianAtPoint( b ) > 0 &&
							        triangles[i]->jacobianAtPoint( c ) > 0 &&
							        triangles[i]->jacobianAtPoint( d ) > 0
							  )
							{
								triangles[i]->moved = true ;

								for( size_t j = 0 ; j < 3 ; j++ )
								{
									if( triangles[i]->getNeighbour( j )->isTriangle )
									{
										dynamic_cast<DelaunayTriangle *>( triangles[i]->getNeighbour( j ) )->moved = true ;
									}
								}
							}
							else
							{
								for( size_t ni = 0 ; ni < n ; ni++ )
								{
									k = indexes[3 * ni + 2] ;
									triangles[i]->getBoundingPoint( k ) = originalPoints[ni] ;
								}
							}
						}
					}

				}

				if( count % 1000 == 0 )
					std::cerr << "\r projecting points on boundaries... triangle " << count << "/" << triangles.size() << " feature " << i << std::flush ;

			}
		}
	}

	std::cerr << "\r projecting points on boundaries... point " << count << "/" << pd << " ...done." << std::endl ;
}

void FeatureTree::stitch()
{
	size_t count = 0 ;
	size_t pd = 0 ;

	if( is2D() )
	{
		if( elemOrder >= QUADRATIC )
		{

			layer2d.begin()->second->setElementOrder( elemOrder, realDeltaTime ) ;
			for(auto i = ++layer2d.begin() ; i != layer2d.end() ; i++)
				dynamic_cast<DelaunayTree *>(i->second)->addSharedNodes( dynamic_cast<DelaunayTree *>(layer2d.begin()->second) ) ;
			
			if( projectOnBoundaries)
			{
				switch( elemOrder )
				{
					case QUADRATIC:
						projectTrianglesOnBoundaries( 1, 0 ) ;
						break ;
					case CUBIC:
						projectTrianglesOnBoundaries( 2, 0 ) ;
						break ;
					case QUADRIC:
					case QUINTIC:
						projectTrianglesOnBoundaries( 3, 0 ) ;
						break ;
					case QUADRATIC_TIME_LINEAR:
						projectTrianglesOnBoundaries( 1, 1 ) ;
						break ;
					case QUADRATIC_TIME_QUADRATIC:
						projectTrianglesOnBoundaries( 1, 2 ) ;
						break ;
					case CUBIC_TIME_LINEAR:
						projectTrianglesOnBoundaries( 2, 1 ) ;
						break ;
					case CUBIC_TIME_QUADRATIC:
						projectTrianglesOnBoundaries( 2, 2 ) ;
						break ;
					case QUADRIC_TIME_LINEAR:
					case QUINTIC_TIME_LINEAR:
						projectTrianglesOnBoundaries( 3, 1 ) ;
						break ;
					case QUADRIC_TIME_QUADRATIC:
					case QUINTIC_TIME_QUADRATIC:
						projectTrianglesOnBoundaries( 3, 2 ) ;
						break ;
				}
			}
		}
	}
	else if( is3D() )
	{
		if( elemOrder >= QUADRATIC )
		{
			layer3d.begin()->second->setElementOrder( elemOrder, realDeltaTime ) ;
			for(auto i = ++layer3d.begin() ; i != layer3d.end() ; i++)
				dynamic_cast<DelaunayTree3D *>(i->second)->addSharedNodes( dynamic_cast<DelaunayTree3D *>(layer3d.begin()->second) ) ;
			
			if( projectOnBoundaries )
			{
				switch( elemOrder )
				{
					case QUADRATIC:
						projectTetrahedronsOnBoundaries( 1, 0 ) ;
						break ;
					case CUBIC:
						projectTetrahedronsOnBoundaries( 2, 0 ) ;
						break ;
					case QUADRIC:
					case QUINTIC:
						projectTetrahedronsOnBoundaries( 3, 0 ) ;
						break ;
					case QUADRATIC_TIME_LINEAR:
						projectTetrahedronsOnBoundaries( 1, 1 ) ;
						break ;
					case QUADRATIC_TIME_QUADRATIC:
						projectTetrahedronsOnBoundaries( 1, 2 ) ;
						break ;
					case CUBIC_TIME_LINEAR:
						projectTetrahedronsOnBoundaries( 2, 1 ) ;
						break ;
					case CUBIC_TIME_QUADRATIC:
						projectTetrahedronsOnBoundaries( 2, 2 ) ;
						break ;
					case QUADRIC_TIME_LINEAR:
					case QUINTIC_TIME_LINEAR:
						projectTetrahedronsOnBoundaries( 3, 1 ) ;
						break ;
					case QUADRIC_TIME_QUADRATIC:
					case QUINTIC_TIME_QUADRATIC:
						projectTetrahedronsOnBoundaries( 3, 2 ) ;
						break ;
				}
			}

		}
	}
	
	if(instants.size() > 2 && elemOrder >= CONSTANT_TIME_LINEAR && is2D())
		for(auto i = layer2d.begin() ; i != layer2d.end() ; i++)
			i->second->extrude(instants) ;
	else if(instants.size() > 2 && elemOrder >= CONSTANT_TIME_LINEAR && is3D())
		for(auto i = layer3d.begin() ; i != layer3d.end() ; i++)
			i->second->extrude(instants) ;
}

void FeatureTree::setSamplingNumber( size_t news )
{
	samplingNumber = news ;
	needMeshing = true ;
	state.enriched = false ;
}

void FeatureTree::quadTreeRefine(const Geometry * location)
{
	if(location->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
	{
		std::cerr << "quadtree refine... " << std::flush ;
		std::vector<DelaunayTriangle *> conflictingElements = location? dtree->getConflictingElements(location):dtree->getElements();
		std::cerr << conflictingElements.size() << " elements... " << std::flush ;
		std::vector<Point> pointsToAdd ;
		std::vector<Point> illegalPoints ;
		for(size_t i = 0 ; i < conflictingElements.size() ; i++)
		{
			if(!MinimumAngle(M_PI/7.).meetsCriterion(conflictingElements[i]))
			{
				bool inrefinedFeature= false ;
				for(size_t j = 0 ; j < refinedFeatures.size() ; j++)
				{
					if(!refinedFeatures[j]->isVirtualFeature && refinedFeatures[j]->in(conflictingElements[i]->getCenter() ))
					{
						inrefinedFeature = true ;
						break ;
					}
				}
				if(inrefinedFeature)
					continue ;
				
	// 			conflictingElements[i]->print() ;
				Point a = *conflictingElements[i]->first*.5+*conflictingElements[i]->second*.5 ;
				Point b = *conflictingElements[i]->first*.5+*conflictingElements[i]->third*.5 ;
				Point c = (a+b+*conflictingElements[i]->second+*conflictingElements[i]->third)*.25 ;
				double d0 = dist(conflictingElements[i]->first, conflictingElements[i]->second) ;
				double d1 = dist(conflictingElements[i]->first, conflictingElements[i]->third) ;
				double d2 = dist(conflictingElements[i]->third, conflictingElements[i]->second) ;
				if(d0 < d1 && d0 < d2)
				{
					illegalPoints.push_back(a);
					a = *conflictingElements[i]->second*.5 + *conflictingElements[i]->third*.5 ;
					c = (a+b+*conflictingElements[i]->second+*conflictingElements[i]->first)*.25 ;
				}
				if(d1 < d0 && d1 < d2)
				{
					illegalPoints.push_back(b);
					b = *conflictingElements[i]->third*.5 + *conflictingElements[i]->second*.5 ;
					c = (a+b+*conflictingElements[i]->third+*conflictingElements[i]->first)*.25 ;
				}
				
				bool uniquea = true ;
				bool uniqueb = true ;
				bool uniquec = true ;
				
				for(size_t j = 0 ; j < pointsToAdd.size() ; j++)
				{
					if(uniquea && dist(a, pointsToAdd[j]) < POINT_TOLERANCE_2D)
					{
						uniquea = false ;
					}
					if(uniqueb && dist(b, pointsToAdd[j]) < POINT_TOLERANCE_2D)
					{
						uniqueb = false ;
					}
					if(uniquec && dist(c, pointsToAdd[j]) < POINT_TOLERANCE_2D)
					{
						uniquec = false ;
					}
					if(!uniquea && !uniqueb && !uniquec)
						break ;
				}
				
				for(size_t j = 0 ; j < illegalPoints.size() ; j++)
				{
					if(uniquea && dist(a, illegalPoints[j]) < POINT_TOLERANCE_2D)
					{
						uniquea = false ;
					}
					if(uniqueb && dist(b, illegalPoints[j]) < POINT_TOLERANCE_2D)
					{
						uniqueb = false ;
					}
					if(uniquec && dist(c, illegalPoints[j]) < POINT_TOLERANCE_2D)
					{
						uniquec = false ;
					}
					if(!uniquea && !uniqueb && !uniquec)
						break ;
				}
				
				if(uniquea)
					pointsToAdd.push_back(a);
				if(uniqueb)
					pointsToAdd.push_back(b);
			}
			
			bool inrefinedFeature= false ;
			for(size_t j = 0 ; j < refinedFeatures.size() ; j++)
			{
				if(refinedFeatures[j]->in(conflictingElements[i]->getCenter() ))
				{
					Point proj = conflictingElements[i]->getCenter() ;
					tree[0]->project(&proj) ;
					if(dist(proj,conflictingElements[i]->getCenter() ) > 2.*conflictingElements[i]->getRadius() )
					{
						inrefinedFeature = true ;
						break ;
					}
				}
			}
			if(inrefinedFeature)
				continue ;
			
// 			conflictingElements[i]->print() ;
			Point a = *conflictingElements[i]->first*.5+*conflictingElements[i]->second*.5 ;
			Point b = *conflictingElements[i]->first*.5+*conflictingElements[i]->third*.5 ;
			Point c = *conflictingElements[i]->third*.5+*conflictingElements[i]->second*.5 ;
			
			bool uniquea = true ;
			bool uniqueb = true ;
			bool uniquec = true ;
			
			for(size_t j = 0 ; j < pointsToAdd.size() ; j++)
			{
				if(uniquea && dist(a, pointsToAdd[j]) < POINT_TOLERANCE_2D)
				{
					uniquea = false ;
				}
				if(uniqueb && dist(b, pointsToAdd[j]) < POINT_TOLERANCE_2D)
				{
					uniqueb = false ;
				}
				if(uniquec && dist(c, pointsToAdd[j]) < POINT_TOLERANCE_2D)
				{
					uniquec = false ;
				}
				if(!uniquea && !uniqueb && !uniquec)
					break ;
			}
			
			for(size_t j = 0 ; j < illegalPoints.size() ; j++)
			{
				if(uniquea && dist(a, illegalPoints[j]) < POINT_TOLERANCE_2D)
				{
					uniquea = false ;
				}
				if(uniqueb && dist(b, illegalPoints[j]) < POINT_TOLERANCE_2D)
				{
					uniqueb = false ;
				}
				if(uniquec && dist(c, illegalPoints[j]) < POINT_TOLERANCE_2D)
				{
					uniquec = false ;
				}
				if(!uniquea && !uniqueb && !uniquec)
					break ;
			}
			
			if(uniquea)
				pointsToAdd.push_back(a);
			if(uniqueb)
				pointsToAdd.push_back(b);
			if(uniquec)
				pointsToAdd.push_back(c);
		}
		
		for(size_t i = 0 ; i < pointsToAdd.size() ; i++)
		{
			additionalPoints.push_back(new Point(pointsToAdd[i]));
			for(auto j = layer2d.begin() ; j != layer2d.end() ; j++)
			{
				j->second->insert(additionalPoints.back()) ;
			}
		}
		std::cerr <<  " ...done. " << std::endl ;
	}
}

void FeatureTree::sample()
{


	if( samplingNumber != previousSamplingNumber )
	{
		meshPoints.clear();
		previousSamplingNumber = samplingNumber ;

		if( is2D() )
		{
			std::cerr << "2D features " << tree.size() << std::endl ;
			double total_area = tree[0]->area() ;
			
			double correctionfactor = 1. ;
			if(samplingFactors.find(tree[0]) != samplingFactors.end())
				correctionfactor =  samplingFactors[tree[0]] ;
			tree[0]->sample( correctionfactor * samplingNumber * 4) ;
			int count = 0 ; 
			
// 			#pragma omp parallel for reduction(+:count) schedule(auto)
			for( size_t i  = 1 ; i < tree.size() ; i++ )
			{
				double shape_factor = ( sqrt( tree[0]->area() ) / ( 2.*M_PI * tree[0]->getRadius() ) ) / ( sqrt( tree[i]->area() ) / ( 2.*M_PI * tree[i]->getRadius() ) );

				if(!tree[i]->isVirtualFeature)
					tree[i]->isUpdated = false ;
				
				if( shape_factor < POINT_TOLERANCE_2D )
				{
					continue ;
				}

				size_t npoints = ( size_t )round( sqrt( tree[i]->area() / ( total_area * shape_factor ) ) * samplingNumber ) ;
				if(npoints < 8 && npoints >= 5)
					npoints = 8 ;
// 				size_t npoints = std::max( ( size_t )round( sqrt( tree[i]->area() / ( total_area * shape_factor ) ) * samplingNumber ), (size_t) 8 ) ;
				correctionfactor = 1. ;
				if(samplingFactors.find(tree[i]) != samplingFactors.end())
				{
					correctionfactor = samplingFactors[tree[i]] ;
					npoints = ( size_t )round(correctionfactor*npoints) ;
				}
				
				if(!tree[i]->isVirtualFeature)
				{
					for(size_t j = 0 ;  j < refinementZones.size() ; j++)
					{
						if(tree[i]->intersects(refinementZones[j]) || refinementZones[j]->in(tree[i]->getCenter()))
						{
							npoints *= 2 ;

							refinedFeatures.push_back(tree[i]);

						}
					}
				}
				
				if(samplingRestriction == SAMPLE_RESTRICT_8)
				{
					if( npoints >= 8 && !tree[i]->isVirtualFeature /*&& npoints < correctionfactor*samplingNumber */)
					{
						count++ ;
						tree[i]->sample( npoints ) ;
					}
				}
				else if(samplingRestriction == SAMPLE_RESTRICT_4)
				{
					if( npoints >= 4 && !tree[i]->isVirtualFeature /*&& npoints < correctionfactor*samplingNumber */)
					{
						count++ ;
						tree[i]->sample( npoints ) ;
					}
				}
				else if(samplingRestriction == SAMPLE_RESTRICT_16)
				{
					if( npoints >= 16 && !tree[i]->isVirtualFeature /*&& npoints < correctionfactor*samplingNumber */)
					{
						count++ ;
						tree[i]->sample( npoints ) ;
					}
				}
				else
				{
					if( !tree[i]->isVirtualFeature )
					{
						count++ ;
						tree[i]->sample( npoints ) ;
					}
				}

//				tree[i]->addMeshPointsInFather() ;
			}
			std::cout << count << " particles meshed" << std::endl ;
		}
		else if( is3D() )
		{
//			std::cout << samplingNumber << std::endl ;
			std::cerr << "\r 3D features... sampling feature 0/" << this->tree.size() << "          " << std::flush ;
			tree[0]->sample( samplingNumber ) ;

			double total_area = tree[0]->area() * tree[0]->area() / ( 4.*M_PI * tree[0]->getRadius() * tree[0]->getRadius() ) * ( tree[0]->area() / ( 4.*M_PI * tree[0]->getRadius() * tree[0]->getRadius() ) ) ;
			int count = 0 ;
			
			#pragma omp parallel for reduction(+:count) schedule(runtime)
			for( int i  = 1 ; i < ( int )tree.size() ; i++ )
			{
				std::cerr << "\r 3D features... sampling feature " << count << "/" << this->tree.size() << "          " << std::flush ;

				double shape_factor = tree[i]->area() / ( 4.*M_PI * tree[i]->getRadius() * tree[i]->getRadius() );
				size_t npoints = ( size_t )round( ( 1.5 * samplingNumber * tree[i]->area() * shape_factor ) / ( total_area ) ) ;
				if(samplingFactors.find(tree[i]) != samplingFactors.end())
				{
					npoints = ( size_t )round(samplingFactors[tree[i]]*npoints) ;
				}
				
				if(samplingRestriction == SAMPLE_RESTRICT_8)
				{
					if( npoints >= 8 && !tree[i]->isVirtualFeature /*&& npoints < correctionfactor*samplingNumber */)
					{
						count++ ;
						tree[i]->sample( npoints ) ;
					}
				}
				else if(samplingRestriction == SAMPLE_RESTRICT_4)
				{
					if( npoints >= 4 && !tree[i]->isVirtualFeature /*&& npoints < correctionfactor*samplingNumber */)
					{
						count++ ;
						tree[i]->sample( npoints ) ;
					}
				}
				else if(samplingRestriction == SAMPLE_RESTRICT_16)
				{
					if( npoints >= 16 && !tree[i]->isVirtualFeature /*&& npoints < correctionfactor*samplingNumber */)
					{
						count++ ;
						tree[i]->sample( npoints ) ;
					}
				}
				else
				{
					if( !tree[i]->isVirtualFeature )
					{
						count++ ;
						tree[i]->sample( npoints ) ;
					}
				}
				if(!tree[i]->isVirtualFeature)
					tree[i]->isUpdated = false ;

			}

			std::cerr << "\r 3D features... sampling feature " << count << "/" << this->tree.size() << " ...done" << std::endl ;
		}
	}
	else
	{

		for( size_t i = 0 ; i < additionalPoints.size() ; i++ )
		{
			delete additionalPoints[i] ;
		}

		additionalPoints.clear();

		meshPoints.clear();
		previousSamplingNumber = samplingNumber ;

		if( is2D() )
		{
			std::cerr << "2D features (updating sampling)" << std::endl ;
			double total_area = tree[0]->area() ;

			if( tree[0]->isUpdated )
				tree[0]->sample( samplingNumber ) ;

			#pragma omp parallel for schedule(runtime) 
			for( size_t i  = 1 ; i < tree.size() ; i++ )
			{
				if( tree[i]->isUpdated )
				{
					grid->remove( tree[i] ) ;
					grid->forceAdd( tree[i] );

					double shape_factor = ( sqrt( tree[0]->area() ) / ( 2.*M_PI * tree[0]->getRadius() ) ) / ( sqrt( tree[i]->area() ) / ( 2.*M_PI * tree[i]->getRadius() ) );

					if( shape_factor < POINT_TOLERANCE_2D )
						continue ;

					size_t npoints = std::max( ( size_t )round( sqrt( tree[i]->area() / ( total_area * shape_factor ) ) * samplingNumber ), ( size_t )8 ) ;
					if(samplingFactors.find(tree[i]) != samplingFactors.end())
					{
						npoints = ( size_t )round(samplingFactors[tree[i]]*npoints) ;
					}
					if( npoints >= 8 && !tree[i]->isVirtualFeature && npoints < samplingNumber )
						tree[i]->sample( npoints ) ;
				}
			}
		}
		else if( is3D() )
		{
//			std::cout << samplingNumber << std::endl ;
			std::cerr << "\r 3D features... sampling feature 0/" << this->tree.size() << "          " << std::flush ;

			if( tree[0]->isUpdated )
				tree[0]->sample( 2.5 * samplingNumber ) ;

			double total_area = tree[0]->area() * tree[0]->area() / ( 4.*M_PI * tree[0]->getRadius() * tree[0]->getRadius() ) * ( tree[0]->area() / ( 4.*M_PI * tree[0]->getRadius() * tree[0]->getRadius() ) ) ;
			int count = 0 ;
			
			#pragma omp parallel for schedule(runtime)
			for( int i  = 1 ; i < ( int )tree.size() ; i++ )
			{
				if( tree[i]->isUpdated )
				{
					grid3d->remove( tree[i] ) ;
					grid3d->forceAdd( tree[i] );
					std::cerr << "\r 3D features... sampling feature " << count << "/" << this->tree.size() << "          " << std::flush ;

					double shape_factor = tree[i]->area() / ( 4.*M_PI * tree[i]->getRadius() * tree[i]->getRadius() );
					size_t npoints = ( size_t )round( ( 1.5 * samplingNumber * tree[i]->area() * shape_factor ) / ( total_area ) ) ;
					if(samplingFactors.find(tree[i]) != samplingFactors.end())
					{
						npoints = ( size_t )round(samplingFactors[tree[i]]*npoints) ;
					}
					if( npoints > 4 && !tree[i]->isVirtualFeature )
						tree[i]->sample( npoints ) ;

					count++ ;
				}
			}

			std::cerr << "\r 3D features... sampling feature " << count << "/" << this->tree.size() << " ...done" << std::endl ;
		}
	}
}

void FeatureTree::refine( size_t nit, SamplingCriterion *cri )
{
	for( size_t t = 0 ; t < 512 ; t++ )
	{
		bool corrected = false ;
		std::vector <DelaunayTriangle *> triangles  =  dtree->getElements() ;

		int count = 0 ;

		for( size_t j = 0;  j < triangles.size() ; j++ )
		{
			if( !cri->meetsCriterion( triangles[j] ) )
			{
				count++ ;
			}
		}

		std::cout << count << " non-conformant triangles " << std::endl ;

		for( size_t j = 0;  j < triangles.size() ; j++ )
		{
			if( !cri->meetsCriterion( triangles[j] ) )
			{
				std::vector<Point> temp = cri->suggest( triangles[j] ) ;

				if( !temp.empty() )
				{
					std::random_shuffle( temp.begin(), temp.end() ) ;
					std::cout << "inserting " << temp.size() << " points" << std::endl ;

					for( size_t i = 0 ; i < temp.size() ; i++ )
					{
						dtree->insert( new Point( temp[i] ) ) ;
						corrected = true ;
					}

					break ;
				}
			}
		}

		if( !corrected )
			break ;
	}
}

void FeatureTree::refine( size_t level )
{
	state.setStateTo( RENUMBERED, false ) ;
	double pointDensity = 0 ;

	if( is2D() )
		pointDensity = .5 * sqrt( tree[0]->area() ) / meshPoints.size() ;
	else
		pointDensity = .5 * pow( tree[0]->volume(), 1. / 3. ) / meshPoints.size() ;

	if( level < 1 )
		return ;

	if( this->dtree == nullptr && this->dtree3D != nullptr )
	{
		std::vector<std::pair<std::vector<Geometry *>, Feature *> >zonesVec ;

		for( size_t j = 1;  j < this->tree.size() ; j++ )
		{
			zonesVec.push_back( std::pair<std::vector<Geometry *>, Feature *>( this->tree[j]->getRefinementZones( level ),  this->tree[j] ) ) ;
		}

		std::vector<Feature *> enrichmentFeature ;

		for( size_t i  = 0 ; i < this->tree.size() ; i++ )
		{
			if( tree[i]->isEnrichmentFeature )
			{
				enrichmentFeature.push_back( tree[i] ) ;
			}
		}

		size_t points_added = 0 ;

		for( size_t i = 0 ; i < zonesVec.size() ; i++ )
		{
			for( size_t j = 0 ; j < zonesVec[i].first.size() ; j++ )
			{
				if( !zonesVec[i].second->isEnrichmentFeature )
				{


					std::vector<Point> toAdd ;
					std::vector<DelaunayTetrahedron *>  tet = this->dtree3D->getConflictingElements( zonesVec[i].first[j] ) ;
					// 			std::vector<DelaunayTriangle *>  * tri_in = this->getBoundingTriangles(zonesVec[i].second) ;

					for( size_t k = 0 ; k < tet.size() ; k++ )
					{
						size_t count_0 = 0 ;
						size_t count_1 = 0 ;
						size_t count_2 = 0 ;
						size_t count_3 = 0 ;
						size_t count_4 = 0 ;
						size_t count_5 = 0 ;

						Point p0  = *tet[k]->first * ( 0.5 ) + *tet[k]->second * ( 0.5 ) ;
						p0.setId( -1 ) ;
						Point p1  = *tet[k]->first * ( 0.5 ) + *tet[k]->third * ( 0.5 ) ;
						p1.setId( -1 ) ;
						Point p2  = *tet[k]->first * ( 0.5 ) + *tet[k]->fourth * ( 0.5 ) ;
						p2.setId( -1 ) ;
						Point p3  = *tet[k]->third * ( 0.5 ) + *tet[k]->second * ( 0.5 ) ;
						p3.setId( -1 ) ;
						Point p4  = *tet[k]->fourth * ( 0.5 ) + *tet[k]->third * ( 0.5 ) ;
						p4.setId( -1 ) ;
						Point p5  = *tet[k]->second * ( 0.5 ) + *tet[k]->fourth * ( 0.5 ) ;
						p5.setId( -1 ) ;


						for( size_t l = 0 ;  l < zonesVec[i].second->getFather()->getChildren().size() ; l++ )
						{
							if( !zonesVec[i].second->getFather()->getChild( l )->isEnrichmentFeature )
							{
								if( zonesVec[i].second->getFather()->getChild( l )->inBoundary( p0, pointDensity ) )
								{
									count_0++ ;
								}

								if( zonesVec[i].second->getFather()->getChild( l )->inBoundary( p1, pointDensity ) )
								{
									count_1++ ;
								}

								if( zonesVec[i].second->getFather()->getChild( l )->inBoundary( p1, pointDensity ) )
								{
									count_2++ ;
								}

								if( zonesVec[i].second->getFather()->getChild( l )->inBoundary( p2, pointDensity ) )
								{
									count_2++ ;
								}

								if( zonesVec[i].second->getFather()->getChild( l )->inBoundary( p3, pointDensity ) )
								{
									count_3++ ;
								}

								if( zonesVec[i].second->getFather()->getChild( l )->inBoundary( p4, pointDensity ) )
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

						if( count_0 == 0 )
						{
							toAdd.push_back( p0 ) ;
						}

						if( count_1 == 0 )
						{
							toAdd.push_back( p1 ) ;
						}

						if( count_2 == 0 )
						{
							toAdd.push_back( p2 ) ;
						}

						if( count_3 == 0 )
						{
							toAdd.push_back( p3 ) ;
						}

						if( count_4 == 0 )
						{
							toAdd.push_back( p4 ) ;
						}

						if( count_5 == 0 )
						{
							toAdd.push_back( p5 ) ;
						}

					}


					std::sort( toAdd.begin(), toAdd.end() ) ;
					auto e = std::unique( toAdd.begin(), toAdd.end() ) ;
					toAdd.erase( e, toAdd.end() ) ;

					// 			std::cerr << "we have " << toAdd.size() << " points for refinement" << std::endl ;

					std::random_shuffle( toAdd.begin(), toAdd.end() ) ;

					for( size_t k = 0 ; k < toAdd.size() ; k++ )
					{
						std::cerr << "\r refining feature " << i + 1 << "/" << zonesVec.size() << "...added " << points_added << " points" << std::flush ;
						Point *p = new Point( toAdd[k] ) ;
						points_added++ ;
						meshPoints.push_back( std::pair<Point *, const Feature *>( p, zonesVec[i].second->getFather() ) ) ;
						dtree3D->insert( p ) ;
					}
				}
			}
		}

		std::cerr << "...done" << std::endl ;
	}
	else
	{
		std::vector<std::pair<std::vector<Geometry *>, Feature *> >zonesVec ;

		for( size_t j = 1;  j < this->tree.size() ; j++ )
		{
			zonesVec.push_back( std::pair<std::vector<Geometry *>, Feature *>( this->tree[j]->getRefinementZones( level ),  this->tree[j] ) ) ;
		}

		std::vector<Feature *> enrichmentFeature ;

		for( size_t i  = 0 ; i < this->tree.size() ; i++ )
		{
			if( tree[i]->isEnrichmentFeature )
			{
				enrichmentFeature.push_back( tree[i] ) ;
			}
		}

		size_t points_added = 0 ;

		for( size_t i = 0 ; i < zonesVec.size() ; i++ )
		{

			for( size_t j = 0 ; j < zonesVec[i].first.size() ; j++ )
			{

				if( !zonesVec[i].second->isEnrichmentFeature )
				{
					std::vector<Point *> sample = zonesVec[i].second->doubleSurfaceSampling() ;

					for( size_t k = 0 ; k < sample.size() ; k++ )
					{
						if( tree[0]->in( *sample[k] ) )
						{
							bool yes = true ;

							for( size_t l = 0 ; l < enrichmentFeature.size() ; l++ )
							{
								if( enrichmentFeature[l]->inBoundary( *sample[k], pointDensity ) )
								{
									yes = false ;
									break ;
								}
							}

							for( size_t l = 0 ; l < zonesVec[i].second->getChildren().size() ; l++ )
							{
								if( zonesVec[i].second->getChild( l )->inBoundary( *sample[k], pointDensity ) )
								{
									yes = false ;
									break ;
								}
							}

							if( yes )
							{
								dtree->insert( sample[k] ) ;
								this->meshPoints.push_back( std::make_pair( sample[k], zonesVec[i].second ) ) ;
							}
						}
					}
				}

				std::vector<Point> toAdd ;

				std::vector<DelaunayTriangle *>  tri = this->dtree->getConflictingElements( zonesVec[i].first[j] ) ;
// 			std::vector<DelaunayTriangle *>  * tri_in = this->getBoundingTriangles(zonesVec[i].second) ;

				for( size_t k = 0 ; k < tri.size() ; k++ )
				{
// 				if((*tri)[k]->area() > 4e-4)
// 				{
					size_t count_0 = 0 ;
					size_t count_1 = 0 ;
					size_t count_2 = 0 ;
					double rand0 = 0 ;/*((2.*rand()/(RAND_MAX+1.0))-1.)*0.1 ;*/
					double rand1 = 0 ;/*((2.*rand()/(RAND_MAX+1.0))-1.)*0.1 ;*/
					double rand2 = 0 ;/*((2.*rand()/(RAND_MAX+1.0))-1.)*0.1 ;*/

					Point p0  = *tri[k]->first * ( 0.5 + rand0 ) + *tri[k]->second * ( 0.5 - rand0 ) ;
					p0.setId( -1 ) ;
					Point p1  = *tri[k]->first * ( 0.5 + rand1 ) + *tri[k]->third * ( 0.5 - rand1 ) ;
					p1.setId( -1 ) ;
					Point p2  = *tri[k]->second * ( 0.5 + rand2 ) + *tri[k]->third * ( 0.5 - rand2 ) ;
					p2.setId( -1 ) ;
//

					for( size_t l = 0 ;  l < zonesVec[i].second->getFather()->getChildren().size() ; l++ )
					{
						if( !zonesVec[i].second->getFather()->getChild( l )->isEnrichmentFeature )
						{
							if( zonesVec[i].second->getFather()->getChild( l )->inBoundary( p0, pointDensity ) )
							{
								count_0++ ;
							}

							if( zonesVec[i].second->getFather()->getChild( l )->inBoundary( p1, pointDensity ) )
							{
								count_1++ ;
							}

							if( zonesVec[i].second->getFather()->getChild( l )->inBoundary( p1, pointDensity ) )
							{
								count_2++ ;
							}
						}
					}



					for( size_t m = 0 ; m < enrichmentFeature.size() ; m++ )
					{
						if( enrichmentFeature[m]->inBoundary( p0, pointDensity ) )
							count_0++ ;

						if( enrichmentFeature[m]->inBoundary( p1, pointDensity ) )
							count_1++ ;

						if( enrichmentFeature[m]->inBoundary( p2, pointDensity ) )
							count_2++ ;
					}

					if( count_0 == 0 && zonesVec[i].first[j]->in( *tri[k]->first ) && zonesVec[i].first[j]->in( *tri[k]->second ) )
					{
						toAdd.push_back( p0 ) ;
					}

					if( count_1 == 0 && zonesVec[i].first[j]->in( *tri[k]->first ) && zonesVec[i].first[j]->in( *tri[k]->third ) )
					{
						toAdd.push_back( p1 ) ;
					}

					if( count_2 == 0 && zonesVec[i].first[j]->in( *tri[k]->second ) && zonesVec[i].first[j]->in( *tri[k]->third ) )
					{
						toAdd.push_back( p2 ) ;
					}

// 				}

				}


				std::sort( toAdd.begin(), toAdd.end() ) ;
				auto e = std::unique( toAdd.begin(), toAdd.end() );
				toAdd.erase( e, toAdd.end() ) ;

// 			std::cerr << "we have " << toAdd.size() << " points for refinement" << std::endl ;
				std::random_shuffle( toAdd.begin(), toAdd.end() ) ;

				for( size_t k = 0 ; k < toAdd.size() ; k++ )
				{
					std::cerr << "\r refining feature " << i + 1 << "/" << zonesVec.size() << "...added " << points_added << " points" << std::flush ;
					Point *p = new Point( toAdd[k] ) ;
					points_added++ ;
					meshPoints.push_back( std::pair<Point *, const Feature *>( p, zonesVec[i].second->getFather() ) ) ;
					dtree->insert( p ) ;
				}
			}
		}

		std::cerr << "...done" << std::endl ;
	}
}

Form * FeatureTree::getElementBehaviour( const DelaunayTriangle *t, int layer,  bool onlyUpdate ) const
{
	int root_box = 0 ;

	if( !inRoot( t->getCenter() ))
	{
		return new VoidForm();
	}

	if( t->getBoundingPoints().size() % 3 != 0 )
		return new VoidForm() ;

	std::vector<Feature *> targets ;

	if( tree.size() > 32 )
	{
		std::vector<const Geometry *> targetstmp = grid->coOccur( t->getPrimitive() ) ;

		for( size_t i = 0 ; i < targetstmp.size() ; i++ )
		{
			const Feature * tmp = dynamic_cast<const Feature *>( targetstmp[i] ) ;
			if(tmp->getLayer() == layer && (tmp->getBoundingPoints().size() || tmp->isVirtualFeature || samplingRestriction == SAMPLE_NO_RESTRICTION))
				targets.push_back( const_cast<Feature *>(tmp) ) ;
		}
	}
	else
	{
		for( size_t i = 0 ; i < tree.size() ; i++ )
		{
			if(tree[i]->getLayer() == layer && (tree[i]->getBoundingPoints().size() || tree[i]->isVirtualFeature || samplingRestriction == SAMPLE_NO_RESTRICTION))
				targets.push_back( tree[i] ) ;
		}
	}

	Form * found = nullptr ;
	
	
	if( !targets.empty() )
	{

		for( int i = targets.size() - 1 ; i >= 0  ; i-- )
		{
			
			if( !targets[i]->isEnrichmentFeature && targets[i]->in( t->getCenter() ) && ( !onlyUpdate || onlyUpdate && targets[i]->isUpdated ) )
			{
				bool notInChildren  = true ;

				std::vector<Feature *> descendants = targets[i]->getDescendants() ;

				for( size_t j = 0 ; j < descendants.size() ; j++ )
				{
					if(descendants[j]->getLayer() == layer && !descendants[j]->isEnrichmentFeature && descendants[j]->in( t->getCenter() ) )
					{
						notInChildren = false ;
						break ;
					}
				}

				if( notInChildren)
				{
					if( targets[i]->getBehaviour( t->getCenter() )->timeDependent() )
					{
						if( !targets[i]->getBehaviour( t->getCenter() )->spaceDependent() )
						{
							 Form *b = targets[i]->getBehaviour( t->getCenter() )->getCopy() ;
							 if(targets[i]->getBehaviourSource())
								b->setSource(targets[i]->getBehaviourSource());
							 else
								b->setSource(targets[i]);
							 found = b ;
						}
						else
						{
							Form *b = targets[i]->getBehaviour( t->getCenter() )->getCopy() ;
							b->transform( t ) ;
							if(targets[i]->getBehaviourSource())
								b->setSource(targets[i]->getBehaviourSource());
							else
								b->setSource(targets[i]);
							found = b ;
						}
					}
					else if( !targets[i]->getBehaviour( t->getCenter() )->spaceDependent() )
					{
						Form *b = targets[i]->getBehaviour( t->getCenter() )->getCopy() ;
						if(targets[i]->getBehaviourSource())
							b->setSource(targets[i]->getBehaviourSource());
						else
							b->setSource(targets[i]);
						 found = b ;
					}
					else
					{
						Form *b = targets[i]->getBehaviour( t->getCenter() )->getCopy() ;
						b->transform( t) ;
						if(targets[i]->getBehaviourSource())
							b->setSource(targets[i]->getBehaviourSource());
						else
							b->setSource(targets[i]);
						found = b ;
					}
				}
			}
		}
	}
	
	if(found)
		return found ;

	if( !onlyUpdate && tree[root_box]->getLayer() == layer)
	{
		if( tree[root_box]->getBehaviour( t->getCenter() )->timeDependent() )
		{
			if( !tree[root_box]->getBehaviour( t->getCenter() )->spaceDependent() )
			{
				Form *b = tree[root_box]->getBehaviour( t->getCenter() )->getCopy() ;
				if(tree[root_box]->getBehaviourSource())
					b->setSource(tree[root_box]->getBehaviourSource());
				else
					b->setSource(tree[root_box]);
				return b ;
			}
			else
			{
				Form *b = tree[root_box]->getBehaviour( t->getCenter() )->getCopy() ;
				if(tree[root_box]->getBehaviourSource())
					b->setSource(tree[root_box]->getBehaviourSource());
				else
					b->setSource(tree[root_box]);
				b->transform( t ) ;
				return b ;
			}
		}
		else if( !tree[root_box]->getBehaviour( t->getCenter() )->spaceDependent() )
		{
			Form *b = tree[root_box]->getBehaviour( t->getCenter() )->getCopy() ;
				if(tree[root_box]->getBehaviourSource())
					b->setSource(tree[root_box]->getBehaviourSource());
				else
					b->setSource(tree[root_box]);
			return b ;
		}
		else
		{
			Form *b = tree[root_box]->getBehaviour( t->getCenter() )->getCopy() ;
				if(tree[root_box]->getBehaviourSource())
					b->setSource(tree[root_box]->getBehaviourSource());
				else
					b->setSource(tree[root_box]);
			b->transform( t ) ;
			return b ;
		}
		Form *b = tree[root_box]->getBehaviour( t->getCenter() )->getCopy() ;
				if(tree[root_box]->getBehaviourSource())
					b->setSource(tree[root_box]->getBehaviourSource());
				else
					b->setSource(tree[root_box]);
		return b ;
	}
	else if(!onlyUpdate && tree[root_box]->getLayer() == layer)
	{
		return new VoidForm() ;
	}

	return new VoidForm() ;

}

Form * FeatureTree::getElementBehaviour( const Amie::DelaunayTetrahedron *t, int layer,  bool onlyUpdate ) const
{
	int root_box = 0 ;

	if( !inRoot( t->getCenter() ) )
	{
		return new VoidForm() ;
	}

// 	std::vector<Geometry *> targetstmp = grid3d->coOccur(t->getPrimitive()) ;
	std::vector<Feature *> targets ;

	if( tree.size() > 32 )
	{
		std::vector<const Geometry *> targetstmp = grid3d->coOccur( t->getPrimitive() ) ;

		for( size_t i = 0 ; i < targetstmp.size() ; i++ )
		{
			const Feature * tmp = dynamic_cast<const Feature *>( targetstmp[i] ) ;
			if(tmp->getLayer() == layer)
				targets.push_back( const_cast<Feature *>(tmp) ) ;
		}
	}
	else
	{
		std::vector<Feature *> targetstmp = tree ;
		for( size_t i = 0 ; i < targetstmp.size() ; i++ )
		{
			if(targetstmp[i]->getLayer() == layer)
				targets.push_back( targetstmp[i] ) ;
		}
	}

	Form * found = nullptr ;
	
	
	if( !targets.empty() )
	{

		for( int i = targets.size() - 1 ; i >= 0  ; i-- )
		{
			
			if( !targets[i]->isEnrichmentFeature && targets[i]->in( t->getCenter() ) && ( !onlyUpdate || onlyUpdate && targets[i]->isUpdated ) )
			{
				bool notInChildren  = true ;

				std::vector<Feature *> descendants = targets[i]->getDescendants() ;

				for( size_t j = 0 ; j < descendants.size() ; j++ )
				{
					if(descendants[j]->getLayer() == layer && !descendants[j]->isEnrichmentFeature && descendants[j]->in( t->getCenter() ) )
					{
						notInChildren = false ;
						break ;
					}
				}

				if( notInChildren)
				{
					if( targets[i]->getBehaviour( t->getCenter() )->timeDependent() )
					{
						if( !targets[i]->getBehaviour( t->getCenter() )->spaceDependent() )
						{
							 Form *b = targets[i]->getBehaviour( t->getCenter() )->getCopy() ;
							 if(targets[i]->getBehaviourSource())
								b->setSource(targets[i]->getBehaviourSource());
							 else
								b->setSource(targets[i]);
							 found = b ;
						}
						else
						{
							Form *b = targets[i]->getBehaviour( t->getCenter() )->getCopy() ;
							b->transform( t ) ;
							if(targets[i]->getBehaviourSource())
								b->setSource(targets[i]->getBehaviourSource());
							else
								b->setSource(targets[i]);
							found = b ;
						}
					}
					else if( !targets[i]->getBehaviour( t->getCenter() )->spaceDependent() )
					{
						Form *b = targets[i]->getBehaviour( t->getCenter() )->getCopy() ;
						if(targets[i]->getBehaviourSource())
							b->setSource(targets[i]->getBehaviourSource());
						else
							b->setSource(targets[i]);
						 found = b ;
					}
					else
					{
						Form *b = targets[i]->getBehaviour( t->getCenter() )->getCopy() ;
						b->transform( t) ;
						if(targets[i]->getBehaviourSource())
							b->setSource(targets[i]->getBehaviourSource());
						else
							b->setSource(targets[i]);
						found = b ;
					}
				}
			}
		}
	}
	
	if(found)
		return found ;

	if( !onlyUpdate && tree[root_box]->getLayer() == layer)
	{
		if( tree[root_box]->getBehaviour( t->getCenter() )->timeDependent() )
		{
			if( !tree[root_box]->getBehaviour( t->getCenter() )->spaceDependent() )
			{
				Form *b = tree[root_box]->getBehaviour( t->getCenter() )->getCopy() ;
				if(tree[root_box]->getBehaviourSource())
					b->setSource(tree[root_box]->getBehaviourSource());
				else
					b->setSource(tree[root_box]);
				return b ;
			}
			else
			{
				Form *b = tree[root_box]->getBehaviour( t->getCenter() )->getCopy() ;
				if(tree[root_box]->getBehaviourSource())
					b->setSource(tree[root_box]->getBehaviourSource());
				else
					b->setSource(tree[root_box]);
				b->transform( t ) ;
				return b ;
			}
		}
		else if( !tree[root_box]->getBehaviour( t->getCenter() )->spaceDependent() )
		{
			Form *b = tree[root_box]->getBehaviour( t->getCenter() )->getCopy() ;
				if(tree[root_box]->getBehaviourSource())
					b->setSource(tree[root_box]->getBehaviourSource());
				else
					b->setSource(tree[root_box]);
			return b ;
		}
		else
		{
			Form *b = tree[root_box]->getBehaviour( t->getCenter() )->getCopy() ;
				if(tree[root_box]->getBehaviourSource())
					b->setSource(tree[root_box]->getBehaviourSource());
				else
					b->setSource(tree[root_box]);
			b->transform( t ) ;
			return b ;
		}
		Form *b = tree[root_box]->getBehaviour( t->getCenter() )->getCopy() ;
				if(tree[root_box]->getBehaviourSource())
					b->setSource(tree[root_box]->getBehaviourSource());
				else
					b->setSource(tree[root_box]);
		return b ;
	}
	else if(!onlyUpdate && tree[root_box]->getLayer() == layer)
	{
		return new VoidForm() ;
	}

	return new VoidForm() ;

}

Point *FeatureTree::checkElement( const DelaunayTetrahedron *t ) const
{
	double pointDensity = 0 ;

	if( is2D() )
		pointDensity = .6 * sqrt( tree[0]->area() ) / meshPoints.size() ;
	else
		pointDensity = .6 * pow( tree[0]->volume(), 1. / 3. ) / meshPoints.size() ;

	if( !inRoot( t->getCenter() ) )
		return nullptr;

	for( int i = tree.size() - 1 ; i >= 0 ; i-- )
	{
		int inCount = tree[i]->in( *t->first )
		              + tree[i]->in( *t->second )
		              + tree[i]->in( *t->third )
		              + tree[i]->in( *t->fourth )
		              + tree[i]->in( t->getCenter() );

		if( inCount > 2 && inRoot( t->getCenter() ) )
		{
			bool inChild = false ;

			std::vector<Feature *> tocheck = tree[i]->getDescendants();

			for( size_t j = 0 ; j < tocheck.size() ; j++ )
			{
				if( tocheck[j]->in( t->getCenter() ) )
				{
					inChild = true ;
					break ;
				}
			}


			if( !inChild )
			{

				size_t count_in = 0 ;

				count_in += tree[i]->inBoundary( *t->first, pointDensity ) ;
				count_in += tree[i]->inBoundary( *t->second, pointDensity ) ;
				count_in += tree[i]->inBoundary( *t->third, pointDensity ) ;
				count_in += tree[i]->inBoundary( *t->fourth, pointDensity ) ;

				if( count_in == 4 && tree[i]->in( t->getCenter() ) )
				{
					return nullptr;
				}

				else
				{
					Point p1( t->getCenter() ) ;
					Point p0( t->getCircumCenter() ) ;
					tree[i]->project( &p1 );
					tree[i]->project( &p0 );
					double d0 = std::min( dist( &p0, t->first ), std::min( dist( &p0, t->second ), dist( &p0, t->third ) ) ) ;
					double d1 = std::min( dist( &p1, t->first ), std::min( dist( &p1, t->second ), dist( &p1, t->third ) ) ) ;

					if( t->inCircumSphere( p0 )
					        && d0 > 1e-8
					        && d0 > d1
					  )
						return new Point( p0 );
					else if( t->inCircumSphere( p1 )
					         && d1 > 1e-8
					         && d1 > d0
					       )
						return new Point( p1 );
					else
						return nullptr ;
				}

			}
		}
	}

	return nullptr;
}

Point *FeatureTree::checkElement( const DelaunayTriangle *t ) const
{
	double pointDensity = 0 ;

	if( is2D() )
		pointDensity = .6 * sqrt( tree[0]->area() ) / meshPoints.size() ;
	else
		pointDensity = .6 * pow( tree[0]->volume(), 1. / 3. ) / meshPoints.size() ;

	if( !inRoot( t->getCenter() ) )
		return nullptr;

	for( int i = tree.size() - 1 ; i >= 0 ; i-- )
	{
		int inCount = tree[i]->in( *t->first )
		              + tree[i]->in( *t->second )
		              + tree[i]->in( *t->third ) + tree[i]->in( t->getCenter() );

		if( inCount > 1 && inRoot( t->getCenter() ) )
		{
			bool inChild = false ;

			std::vector<Feature *> tocheck = tree[i]->getDescendants();

			for( size_t j = 0 ; j < tocheck.size() ; j++ )
			{
				if( tocheck[j]->in( t->getCenter() ) )
				{
					inChild = true ;
					break ;
				}
			}


			if( !inChild )
			{

				size_t count_in = 0 ;

				count_in += tree[i]->inBoundary( *t->first, pointDensity ) ;
				count_in += tree[i]->inBoundary( *t->second, pointDensity ) ;
				count_in += tree[i]->inBoundary( *t->third, pointDensity ) ;

				if( count_in == 3 && tree[i]->in( t->getCenter() ) )
				{
					return nullptr;
				}

				else
				{
					Point p1( t->getCenter() ) ;
					Point p0( t->getCircumCenter() ) ;
					tree[i]->project( &p1 );
					tree[i]->project( &p0 );
					double d0 = std::min( dist( &p0, t->first ), std::min( dist( &p0, t->second ), dist( &p0, t->third ) ) ) ;
					double d1 = std::min( dist( &p1, t->first ), std::min( dist( &p1, t->second ), dist( &p1, t->third ) ) ) ;

					if( t->inCircumCircle( p0 )
					        && d0 > 1e-4
					        && d0 > d1
					  )
						return new Point( p0 );
					else if( t->inCircumCircle( p1 )
					         && d1 > 1e-4
					         && d1 > d0
					       )
						return new Point( p1 );
					else
						return nullptr ;
				}

			}
		}
	}

	return nullptr;
}

Feature *FeatureTree::getFeatForTetra( const DelaunayTetrahedron *t ) const
{

	if( !inRoot( t->getCenter() ) )
		return nullptr;

	for( int i = tree.size() - 1 ; i >= 0 ; i-- )
	{
		if( tree[i]->in( t->getCenter() ) && ( inRoot( t->getCenter() ) ) )
		{
			bool inChild = false ;

			std::vector<Feature *> tocheck = tree[i]->getChildren();
			std::vector<Feature *> tocheckNew =  tocheck;

			while( !tocheckNew.empty() )
			{
				std::vector<Feature *> tocheckTemp ;

				for( size_t k = 0 ; k < tocheckNew.size() ; k++ )
				{
					tocheckTemp.insert(
					    tocheckTemp.end(),
					    tocheckNew[k]->getChildren().begin(),
					    tocheckNew[k]->getChildren().end()
					) ;
				}

				tocheck.insert( tocheck.end(), tocheckTemp.begin(), tocheckTemp.end() ) ;
				tocheckNew = tocheckTemp ;
			}

			for( size_t j = 0 ; j < tocheck.size() ; j++ )
			{
				if( tocheck[j]->in( t->getCenter() ) )
				{
					inChild = true ;
					break ;
				}
			}

			return tree[i];
		}
	}

	return nullptr;
}

void FeatureTree::setElementBehaviours()
{
	if( !father3D )
		father3D = new TetrahedralElement( elemOrder ) ;

	father3D->compileAndPrecalculate() ;

	if( !father2D )
		father2D = new TriElement( elemOrder ) ;

	father2D->compileAndPrecalculate() ;

	if( is2D() )
	{
		double remainder = 0 ;
		for(auto i = layer2d.begin() ; i != layer2d.end() ; i++)
		{
			if(i->first != -1)
				remainder += scalingFactors[i->first] ;
		}
		scalingFactors[-1] = 1.-remainder ;
		layer2d[-1] = dtree ;

		

		for(auto i = layer2d.begin() ; i != layer2d.end() ; i++)
		{

			
			int setcount = 0 ;
			std::vector<DelaunayTriangle *> tris = i->second->getElements() ;
			std::cerr << "\r setting behaviours... triangle : layer (" << scalingFactors[i->first] << ") "<<i->first << "  " << setcount++ << "/" << tris.size() << "    " << std::flush ;
			for( size_t j = 0 ; j < tris.size() ; j++ )
			{
				if( setcount++ % 1000 == 0 )
					std::cerr << "\r setting behaviours... triangle : layer (" << scalingFactors[i->first] << ") "<<i->first << "  " << setcount << "/" << tris.size() << "    " << std::flush ;
				
				tris[j]->refresh( father2D ) ;
				Form * bf =  getElementBehaviour( tris[j], i->first );
				
				tris[j]->setBehaviour(bf) ;

				
				if(!tris[j]->getBehaviour())
				{
					tris[j]->setBehaviour( new VoidForm()) ;
				}
				
			}
			std::cerr << " ...done" << std::endl ;
			
		}


	}
	else
	{
		std::vector<DelaunayTetrahedron *> tetrahedrons = dtree3D->getElements() ;

		std::cerr << " setting behaviours..." << std::flush ;
		int setcount = 0 ;

// 		for( size_t i = 0 ; i < tetrahedrons.size() ; i++ )
// 			tetrahedrons[i]->refresh( father3D ) ;

		for(auto i = layer3d.begin() ; i != layer3d.end() ; i++)
		{
			std::vector<DelaunayTetrahedron *> tets = i->second->getElements() ;
// 			for( size_t j = 0 ; j < tets.size() ; j++ )
// 			{
// 				tets[j]->refresh( father3D ) ;
// 			}
		}
		
		
		for( size_t i = 0 ; i < tetrahedrons.size() ; i++ )
		{
			if( setcount % 1000 == 0 )
				std::cerr << "\r setting behaviours : base layer : tet " << setcount << "/" << tetrahedrons.size() << std::flush ;

			tetrahedrons[i]->refresh( father3D ) ;
			Form * bf =  getElementBehaviour( tetrahedrons[i]);
			
			tetrahedrons[i]->setBehaviour( bf) ;
			
			if(!tetrahedrons[i]->getBehaviour())
			{
				tetrahedrons[i]->setBehaviour( new VoidForm()) ;
			}
			setcount++ ;
		}


		std::cerr << " ...done" << std::endl ;

		
	}

}

void FeatureTree::updateElementBehaviours()
{
	double n_void ;

	if( is2D() )
	{
		double remainder = 0 ;
		for(auto i = layer2d.begin() ; i != layer2d.end() ; i++)
		{
			if(i->first != -1)
				remainder += scalingFactors[i->first] ;
		}
		scalingFactors[-1] = 1.-remainder ;
		layer2d[-1] = dtree ;

		
		for(auto i = layer2d.begin() ; i != layer2d.end() ; i++)
		{
			std::vector<double> scales ;
			if(layer2d.size() > 1)
			{
				for(auto j = ++layer2d.begin() ; j != layer2d.end() ; j++)
				{
					if(i != j)
					{
						scales.push_back(scalingFactors[j->first]/scalingFactors[i->first]);
					}
				}
			}
			
			int setcount = 0 ;
			std::vector<DelaunayTriangle *> tris = i->second->getElements() ;
			std::cerr << "\r updating behaviours... triangle : layer (" << scalingFactors[i->first] << ") "<<i->first << "  " << setcount++ << "/" << tris.size() << "    " << std::flush ;
			for( size_t j = 0 ; j < tris.size() ; j++ )
			{
				
				if( setcount++ % 1000 == 0 )
					std::cerr << "\r updating behaviours... triangle : layer (" << scalingFactors[i->first] << ") "<<i->first << "  " << setcount << "/" << tris.size() << "    " << std::flush ;
				
				tris[j]->refresh( father2D ) ;
				Form * bf =  getElementBehaviour( tris[j], i->first, true );
				
				if(bf )
				{
					tris[j]->setBehaviour( bf ) ;
				}
			}
			std::cerr << " ...done" << std::endl ;
			
		}

		setBehaviours = true ;
	}
	else
	{
		std::vector<DelaunayTetrahedron *> tetrahedrons = this->dtree3D->getElements() ;

		std::cerr << " setting behaviours..." << std::flush ;
		int setcount = 0 ;

		for( size_t i = 0 ; i < tetrahedrons.size() ; i++ )
			tetrahedrons[i]->refresh( father3D ) ;


		
		#pragma omp parallel for schedule(runtime)
		for( size_t i = 0 ; i < tetrahedrons.size() ; i++ )
		{
			if( setcount % 1000 == 0 )
				std::cerr << "\r setting behaviours... tet " << setcount << "/" << tetrahedrons.size() << std::flush ;

			Form *b = getElementBehaviour( tetrahedrons[i], true ) ;

			if( b )
				tetrahedrons[i]->setBehaviour( b) ;

			if( !tetrahedrons[i]->getBehaviour() || tetrahedrons[i]->getBehaviour()->type == VOID_BEHAVIOUR )
				tetrahedrons[i]->setBehaviour( getElementBehaviour( tetrahedrons[i] )) ;

			n_void++ ;
			setcount++ ;
		}

		std::cerr << " ...done" << std::endl ;

		setBehaviours = true ;
	}


	for( size_t i = 0 ; i < tree.size() ; i++ )
	{
		if( !tree[i]->isEnrichmentFeature )
			tree[i]->isUpdated = false ;
	}
}

void FeatureTree::enrich()
{
	enrichmentChange = false ;
	lastEnrichmentId = getNodes().size() ;

	std::cerr << "\r enriching... feature " << 0 << "/" << this->tree.size() << std::flush ;

	for( size_t i = 1 ; i < this->tree.size() ; i++ )
	{
		if( is3D() )
		{
			std::vector<Mesh <DelaunayTetrahedron, DelaunayTreeItem3D > *> extra3dMeshes ;
			if(layer3d.size() > 1)
			{
			  for(auto i = ++layer3d.begin() ; i != layer3d.end() ; i++)
			    extra3dMeshes.push_back(i->second) ;
			}
			if( tree[i]->isEnrichmentFeature && ( dynamic_cast<EnrichmentFeature *>( tree[i] )->moved() || !state.enriched ) )
			{
				if( !state.enriched )
					dynamic_cast<EnrichmentFeature *>( tree[i] )->update( dtree3D ) ;

				dynamic_cast<EnrichmentFeature *>( tree[i] )->enrich( lastEnrichmentId, dtree3D) ;

				enrichmentChange = true ;
				reuseDisplacements = false ;

			}

			if( i % 10 == 0 )
				std::cerr << "\r enriching... feature " << i + 1 << "/" << this->tree.size() << std::flush ;
		}
		else
		{
			std::vector<Mesh <DelaunayTriangle, DelaunayTreeItem > *> extra2dMeshes ;
			if(layer2d.size() > 1)
			{
			  for(auto i = ++layer2d.begin() ; i != layer2d.end() ; i++)
			    extra2dMeshes.push_back(i->second) ;
			}
			
			if( tree[i]->isEnrichmentFeature && ( dynamic_cast<EnrichmentFeature *>( tree[i] )->moved() || !state.enriched ) )
			{
				if( !state.enriched )
					dynamic_cast<EnrichmentFeature *>( tree[i] )->update( dtree ) ;
				dynamic_cast<EnrichmentFeature *>( tree[i] )->enrich( lastEnrichmentId, dtree) ;

				enrichmentChange = true ;
				reuseDisplacements = false ;

			}

			if( i % 10 == 0 )
				std::cerr << "\r enriching... feature " << i + 1 << "/" << this->tree.size() << std::flush ;
		}
	}

	std::cerr << " ...done" << std::endl ;
}

void FeatureTree::assemble()
{
	K->getElements2d().clear();
	K->getElements3d().clear();
	K->getScales().clear();

	if( is2D() )
	{
		numdofs = dtree->getLastNodeId() ;
//		std::cout << deltaTime << std::endl ;
		for(auto i = layer2d.begin() ; i != layer2d.end() ; i++)
		{
			std::vector<DelaunayTriangle *> tris = i->second->getElements() ;
			for( size_t j = 0 ; j < tris.size() ; j++ )
			{
				if( j % 1000 == 0 )
					std::cerr << "\r assembling stiffness matrix, layer "<< i->first << " ... triangle " << j + 1 << "/" << tris.size() << std::flush ;
						
				if(tris[j] && tris[j]->getBehaviour() && tris[j]->getBehaviour()->type != VOID_BEHAVIOUR  )
				{
					tris[j]->refresh( father2D ) ;
					tris[j]->getBehaviour()->preProcess( deltaTime, tris[j]->getState() ) ;
// 					if(!tris[j]->getBehaviour()->fractured())
						K->add( tris[j], scalingFactors[i->first] ) ;
				}
				
			}
			K->setMaxDof( std::max(getNodes().size(),lastEnrichmentId) ) ;
			std::cerr << " ...done." << std::endl ;
		}

	}
	else
	{
		numdofs = dtree3D->getLastNodeId() ;

		for(auto i = layer3d.begin() ; i!= layer3d.end() ; i++)
		{
			std::vector<DelaunayTetrahedron *> tets = i->second->getElements() ;
			for( size_t j = 0 ; j < tets.size() ; j++ )
			{
				if( tets[j]->getBehaviour() && tets[j]->getBehaviour()->type != VOID_BEHAVIOUR )
				{
					if( j % 1000 == 0 )
						std::cerr << "\r assembling stiffness matrix... tetrahedron " << j + 1 << "/" << tets.size() << std::flush ;

					tets[j]->refresh( father3D ) ;
					tets[j]->getBehaviour()->preProcess( deltaTime, tets[j]->getState() ) ;
					K->add( tets[j] ) ;
				}
			}
			K->setMaxDof( std::max(getNodes().size(),lastEnrichmentId) ) ;
			std::cerr << " ...done." << std::endl ;
		}

	}
}

std::vector<DelaunayTriangle> FeatureTree::getSnapshot2D() const
{
	std::vector<DelaunayTriangle> copy ;

	std::vector<DelaunayTriangle *> tris = dtree->getElements() ;
	
	for( size_t i = 0 ; i < tris.size() ; i++ )
	{
		copy.push_back( *tris[i] ) ;
		copy.back().setBehaviour( tris[i]->getBehaviour()->getCopy()) ;
		copy.back().getState().initialize(false) ;
	}

	return copy ;
}

Vector FeatureTree::stressFromDisplacements()
{
	state.setStateTo( XFEM_STEPPED, false ) ;

	VirtualMachine vm ;
	if( dtree)
	{
		std::vector<DelaunayTriangle *> elements = dtree->getElements() ;
		Vector stress( 0.f, 3 * 3 * elements.size() ) ;

		for( size_t i  = 0 ; i < elements.size() ; i++ )
		{
			if( elements[i]->getBehaviour() && elements[i]->getBehaviour()->type != VOID_BEHAVIOUR )
			{
				std::valarray<Point *> pts( 3 ) ;
				pts[0] =  elements[i]->first ;
				pts[1] =  elements[i]->second ;
				pts[2] =  elements[i]->third ;
				
				Vector str(0., 3*3) ;
				elements[i]->getState().getField(REAL_STRESS_FIELD, pts, str, false, &vm) ;

				for( size_t j = 0 ; j < 9 ; j++ )
					stress[i * 3 * 3 + j] = str[j] ;

				std::cerr << "\r computing stress... element " << i + 1 << "/" << elements.size() << std::flush ;
			}

		}

		std::cerr << " ...done." << std::endl ;
		return stress ;
	}
	else
	{
		std::vector<DelaunayTetrahedron *> elements3D = dtree3D->getElements() ;
		Vector stress( 0., 4 * 6 * elements3D.size() ) ;

		for( size_t i  = 0 ; i < elements3D.size() ; i++ )
		{
			std::valarray<Point *> pts( 4 ) ;
			pts[0] =  elements3D[i]->first ;
			pts[1] =  elements3D[i]->second ;
			pts[2] =  elements3D[i]->third ;
			pts[2] =  elements3D[i]->fourth ;

			Vector str(0., 24) ;
			elements3D[i]->getState().getField(REAL_STRESS_FIELD, pts, str, false, &vm) ;

			for( size_t j = 0 ; j < 24 ; j++ )
				stress[i * 4 * 6 + j] = str[j] ;

			std::cerr << "\r computing stress... element " << i + 1 << "/" << elements3D.size() << std::flush ;

		}

		std::cerr << " ...done." << std::endl ;
		return stress ;
	}
}

const Vector &FeatureTree::getDisplacements( int g, bool stepTree )
{
	if(stepTree)
		state.setStateTo( XFEM_STEPPED, false ) ;

	return  K->getDisplacements() ;
}

Vector FeatureTree::getDisplacements(Point * pt, int g , bool stepTree )
{
	if(stepTree)
		state.setStateTo( XFEM_STEPPED, false ) ;
	
	if(is2D())
	{
		Vector ret ;
		std::vector<DelaunayTriangle *> elements = dtree->getConflictingElements(pt) ;
		for(size_t i = 0 ; i < elements.size() ; i++)
		{
			
			if(elements[i]->in(*pt))
			{
				ret.resize( elements[i]->getBehaviour()->getNumberOfDegreesOfFreedom() , 0.) ;
				elements[i]->getState().getField(DISPLACEMENT_FIELD, *pt, ret, false) ;
				return ret ;
			}
		}
		ret.resize( 2 , 0.) ;
		return ret ;
	}
	else
	{
		Vector ret ;
		std::vector<DelaunayTetrahedron *> elements = dtree3D->getConflictingElements(pt) ;
		for(size_t i = 0 ; i < elements.size() ; i++)
		{
			if(elements[i]->in(*pt))
			{
				ret.resize( elements[i]->getBehaviour()->getNumberOfDegreesOfFreedom() , 0.) ;
				elements[i]->getState().getField(DISPLACEMENT_FIELD, *pt, ret, false) ;
				return ret ;
			}
		}
		ret.resize( 2 , 0.) ;
		return ret ;
	}
}

std::pair<Vector , Vector > FeatureTree::getStressAndStrain( int g, bool stepTree )
{
	if(stepTree)
		state.setStateTo( XFEM_STEPPED, false ) ;

	VirtualMachine vm ;
	if( dtree)
	{
		std::vector<DelaunayTriangle *> elements = dtree->getElements() ;

		std::pair<Vector , Vector > stress_strain( Vector( 0., elements[0]->getBoundingPoints().size() * 3 * elements.size() ), Vector( 0., elements[0]->getBoundingPoints().size() * 3 * elements.size() ) ) ;
		int donecomputed = 0 ;
		
		#pragma omp parallel for shared(donecomputed) schedule(runtime)
		for( size_t i  = 0 ; i < elements.size() ; i++ )
		{
			if( elements[i]->getBehaviour() && elements[i]->getBehaviour()->type != VOID_BEHAVIOUR )
			{
// 				std::valarray<Point *> pts(3) ;
// 				pts[0] =  elements[i]->first ;
// 				pts[1] =  elements[i]->second ;
// 				pts[2] =  elements[i]->third ;

				Vector strain(0., 3*elements[i]->getBoundingPoints().size()) ;
				Vector stress(0., 3*elements[i]->getBoundingPoints().size()) ;
				elements[i]->getState().getField( STRAIN_FIELD, REAL_STRESS_FIELD, elements[i]->getBoundingPoints(), strain, stress, false, &vm ) ;

				for( size_t j = 0 ; j < elements[0]->getBoundingPoints().size() * 3 ; j++ )
				{
					stress_strain.first[i * elements[0]->getBoundingPoints().size() * 3 + j] = stress[j] ;
					stress_strain.second[i * elements[0]->getBoundingPoints().size() * 3 + j] = strain[j] ;
				}

				if( donecomputed % 10000 == 0 )
					std::cerr << "\r computing strain+stress... element " << donecomputed + 1 << "/" << elements.size() << std::flush ;
			}

			donecomputed++ ;
		}

		std::cerr << " ...done." << std::endl ;
		return stress_strain ;
	}
	else
	{
		std::vector<DelaunayTetrahedron *> tets = dtree3D->getElements() ;

		std::pair<Vector , Vector > stress_strain( Vector( 0.f, 4 * 6 * tets.size() ), Vector( 0.f, 4 * 6 * tets.size() ) ) ;
		int donecomputed = 0 ;

		#pragma omp parallel for shared(donecomputed) schedule(runtime)
		for( size_t i  = 0 ; i < tets.size() ; i++ )
		{
			if( tets[i]->getBehaviour() && tets[i]->getBehaviour()->type != VOID_BEHAVIOUR )
			{
				std::valarray<Point *> pts( 4 ) ;
				pts[0] =  tets[i]->first ;
				pts[1] =  tets[i]->second ;
				pts[2] =  tets[i]->third ;
				pts[3] =  tets[i]->fourth ;

				Vector strain(0., 24) ;
				Vector stress(0., 24) ;
				tets[i]->getState().getField( STRAIN_FIELD, REAL_STRESS_FIELD, pts, strain, stress, false, &vm ) ;

				for( size_t j = 0 ; j < 24 ; j++ )
				{
					stress_strain.first[i * 4 * 6 + j] = strain[j] ;
					stress_strain.second[i * 4 * 6 + j] = stress[j] ;
				}
			}

			if( donecomputed % 1000 == 0 )
				std::cerr << "\r computing strain+stress... element " << donecomputed + 1 << "/" << tets.size() << std::flush ;

			donecomputed++ ;
		}

		std::cerr << " ...done." << std::endl ;
		return stress_strain ;
	}
}

std::pair<Vector , Vector > FeatureTree::getStressAndStrainInLayer( int g, bool stepTree )
{
	if(stepTree)
		state.setStateTo( XFEM_STEPPED, false ) ;
	
	VirtualMachine vm ;

	if( dtree )
	{
		std::vector<DelaunayTriangle *> elements = dtree->getElements() ;
		if(g != -1 && layer2d.find(g) != layer2d.end())
			elements = layer2d[g]->getElements() ;

		std::pair<Vector , Vector > stress_strain( Vector( 0., elements[0]->getBoundingPoints().size() * 3 * elements.size() ), Vector( 0., elements[0]->getBoundingPoints().size() * 3 * elements.size() ) ) ;
		int donecomputed = 0 ;
		
		#pragma omp parallel for shared(donecomputed) schedule(runtime)
		for( size_t i  = 0 ; i < elements.size() ; i++ )
		{
			if( elements[i]->getBehaviour() && elements[i]->getBehaviour()->type != VOID_BEHAVIOUR )
			{
// 				std::valarray<Point *> pts(3) ;
// 				pts[0] =  elements[i]->first ;
// 				pts[1] =  elements[i]->second ;
// 				pts[2] =  elements[i]->third ;

				Vector strain(0., 3*elements[i]->getBoundingPoints().size()) ;
				Vector stress(0., 3*elements[i]->getBoundingPoints().size()) ;
				elements[i]->getState().getField( STRAIN_FIELD, REAL_STRESS_FIELD, elements[i]->getBoundingPoints(), strain, stress, false, &vm ) ;

				for( size_t j = 0 ; j < elements[0]->getBoundingPoints().size() * 3 ; j++ )
				{
					stress_strain.first[i * elements[0]->getBoundingPoints().size() * 3 + j] = stress[j] ;
					stress_strain.second[i * elements[0]->getBoundingPoints().size() * 3 + j] = strain[j] ;
				}

				if( donecomputed % 10000 == 0 )
					std::cerr << "\r computing strain+stress... element " << donecomputed + 1 << "/" << elements.size() << std::flush ;
			}

			donecomputed++ ;
		}

		std::cerr << " ...done." << std::endl ;
		return stress_strain ;
	}
	else
	{
		std::vector<DelaunayTetrahedron *> tets = dtree3D->getElements() ;

		if(g != -1 && layer3d.find(g) != layer3d.end())
			tets = layer3d[g]->getElements() ;

		std::pair<Vector , Vector > stress_strain( Vector( 0.f, 4 * 6 * tets.size() ), Vector( 0.f, 4 * 6 * tets.size() ) ) ;
		int donecomputed = 0 ;

		#pragma omp parallel for shared(donecomputed) schedule(runtime)
		for( size_t i  = 0 ; i < tets.size() ; i++ )
		{
			if( tets[i]->getBehaviour() && tets[i]->getBehaviour()->type != VOID_BEHAVIOUR )
			{
				std::valarray<Point *> pts( 4 ) ;
				pts[0] =  tets[i]->first ;
				pts[1] =  tets[i]->second ;
				pts[2] =  tets[i]->third ;
				pts[3] =  tets[i]->fourth ;

				Vector strain(0., 24) ;
				Vector stress(0., 24) ;
				tets[i]->getState().getField( STRAIN_FIELD, REAL_STRESS_FIELD, pts, strain, stress, false, &vm ) ;

				for( size_t j = 0 ; j < 24 ; j++ )
				{
					stress_strain.first[i * 4 * 6 + j] = stress[j] ;
					stress_strain.second[i * 4 * 6 + j] = strain[j] ;
				}
			}

			if( donecomputed % 1000 == 0 )
				std::cerr << "\r computing strain+stress... element " << donecomputed + 1 << "/" << tets.size() << std::flush ;

			donecomputed++ ;
		}

		std::cerr << " ...done." << std::endl ;
		return stress_strain ;
	}
}

std::pair<Vector , Vector > FeatureTree::getStressAndStrainInAllLayers( bool stepTree)
{
	if(stepTree)
		state.setStateTo( XFEM_STEPPED, false ) ;

	VirtualMachine vm ;
	if( dtree != nullptr )
	{
		std::vector<DelaunayTriangle *> elements = getActiveElements2D() ;


		std::pair<Vector , Vector > stress_strain( Vector( 0., elements[0]->getBoundingPoints().size() * 3 * elements.size() ), 
																							 Vector( 0., elements[0]->getBoundingPoints().size() * 3 * elements.size() ) ) ;
		int donecomputed = 0 ;
		
		#pragma omp parallel for shared(donecomputed) schedule(runtime) 
		for( size_t i  = 0 ; i < elements.size() ; i++ )
		{
			if( elements[i]->getBehaviour() && elements[i]->getBehaviour()->type != VOID_BEHAVIOUR )
			{
// 				std::valarray<Point *> pts(3) ;
// 				pts[0] =  elements[i]->first ;
// 				pts[1] =  elements[i]->second ;
// 				pts[2] =  elements[i]->third ;

				Vector strain(0., 3*elements[i]->getBoundingPoints().size()) ;
				Vector stress(0., 3*elements[i]->getBoundingPoints().size()) ;
				elements[i]->getState().getField( STRAIN_FIELD, REAL_STRESS_FIELD, elements[i]->getBoundingPoints(), strain, stress, false, &vm) ;
				for( size_t j = 0 ; j < elements[0]->getBoundingPoints().size() * 3 ; j++ )
				{
					stress_strain.first[i * elements[0]->getBoundingPoints().size() * 3 + j] = stress[j] ;
					stress_strain.second[i * elements[0]->getBoundingPoints().size() * 3 + j] = strain[j] ;
				}

				if( donecomputed % 10000 == 0 )
					std::cerr << "\r computing strain+stress... element " << donecomputed + 1 << "/" << elements.size() << std::flush ;
			}

			donecomputed++ ;
		}

		std::cerr << " ...done." << std::endl ;
		return stress_strain ;
	}
	else
	{
		std::vector<DelaunayTetrahedron *> tets = dtree3D->getElements() ;

		std::pair<Vector , Vector > stress_strain( Vector( 0.f, 4 * 6 * tets.size() ), Vector( 0.f, 4 * 6 * tets.size() ) ) ;
		int donecomputed = 0 ;

		#pragma omp parallel for shared(donecomputed) schedule(runtime)
		for( size_t i  = 0 ; i < tets.size() ; i++ )
		{
			if( tets[i]->getBehaviour() && tets[i]->getBehaviour()->type != VOID_BEHAVIOUR )
			{
				std::valarray<Point *> pts( 4 ) ;
				pts[0] =  tets[i]->first ;
				pts[1] =  tets[i]->second ;
				pts[2] =  tets[i]->third ;
				pts[3] =  tets[i]->fourth ;

				Vector strain(0., 24) ;
				Vector stress(0., 24) ;
				tets[i]->getState().getField( STRAIN_FIELD, REAL_STRESS_FIELD, tets[i]->getBoundingPoints(), strain, stress, false) ;

				for( size_t j = 0 ; j < 24 ; j++ )
				{
					stress_strain.first[i * 4 * 6 + j] = stress[j] ;
					stress_strain.second[i * 4 * 6 + j] = strain[j] ;
				}
			}

			if( donecomputed % 1000 == 0 )
				std::cerr << "\r computing strain+stress... element " << donecomputed + 1 << "/" << tets.size() << std::flush ;

			donecomputed++ ;
		}

		std::cerr << " ...done." << std::endl ;
		return stress_strain ;
	}
}

std::pair<Vector , Vector > FeatureTree::getGradientAndFlux( int g , bool stepTree)
{
	if(stepTree)
		state.setStateTo( XFEM_STEPPED, false ) ;

	if( dtree != nullptr )
	{
		std::vector<DelaunayTriangle *> elements = dtree->getElements() ;

		std::pair<Vector , Vector > grad_flux( Vector( 0., elements[0]->getBoundingPoints().size() * 2 * elements.size() ), Vector( 0., elements[0]->getBoundingPoints().size() * 2 * elements.size() ) ) ;

		for( size_t i  = 0 ; i < elements.size() ; i++ )
		{
			if( elements[i]->getBehaviour() && elements[i]->getBehaviour()->type != VOID_BEHAVIOUR )
			{
// 				std::valarray<Point *> pts(3) ;
// 				pts[0] =  elements[i]->first ;
// 				pts[1] =  elements[i]->second ;
// 				pts[2] =  elements[i]->third ;

				Vector gradient(0., 2*elements[i]->getBoundingPoints().size()) ;
				Vector flux(0., 2*elements[i]->getBoundingPoints().size()) ;
				elements[i]->getState().getField( GRADIENT_FIELD, FLUX_FIELD, elements[i]->getBoundingPoints(), gradient, flux, false) ;

				for( size_t j = 0 ; j < elements[0]->getBoundingPoints().size() * 2 ; j++ )
				{
					grad_flux.first[i * elements[0]->getBoundingPoints().size() * 2 + j] = gradient[j] ;
					grad_flux.second[i * elements[0]->getBoundingPoints().size() * 2 + j] = flux[j] ;
				}

				if( i % 1000 == 0 )
					std::cerr << "\r computing gradient+flux... element " << i + 1 << "/" << elements.size() << std::flush ;
			}
		}

		std::cerr << " ...done." << std::endl ;
		return grad_flux ;
	}
	else
	{
		std::vector<DelaunayTetrahedron *> tets = dtree3D->getElements() ;

		size_t npoints = tets[0]->getBoundingPoints().size() ;
		std::pair<Vector , Vector > grad_flux( Vector( 0.f, npoints * 3 * tets.size() ), Vector( 0.f, npoints * 3 * tets.size() ) ) ;

		for( size_t i  = 0 ; i < tets.size() ; i++ )
		{

			Vector gradient(0., 3*tets[i]->getBoundingPoints().size()) ;
			Vector flux(0., 3*tets[i]->getBoundingPoints().size()) ;
			tets[i]->getState().getField( GRADIENT_FIELD, FLUX_FIELD, tets[i]->getBoundingPoints(), gradient, flux, false) ;

			for( size_t j = 0 ; j < npoints * 3 ; j++ )
			{
				grad_flux.first[i * npoints * 3 + j] = gradient[j] ;
				grad_flux.second[i * npoints * 3 + j] = flux[j] ;
			}

			if( i % 1000 == 0 )
				std::cerr << "\r computing gradient+flux... element " << i + 1 << "/" << tets.size() << std::flush ;

//				std::cout << grflx.first.size() << std::endl ;
		}

		std::cerr << " ...done." << std::endl ;
		return grad_flux ;
	}
}

std::vector<int>FeatureTree:: listLayers() const
{
	std::vector<int> ret ;
	if(is2D())
	{
		for(auto i = layer2d.begin() ; i!= layer2d.end() ; ++i)
			ret.push_back(i->first);
	}
	else
	{
		for(auto i = layer3d.begin() ; i!= layer3d.end() ; ++i)
			ret.push_back(i->first);
	}
	
	return ret ;
}

std::pair<Vector , Vector > FeatureTree::getGradientAndFluxInLayer( int g, bool stepTree )
{
	if(stepTree)
		state.setStateTo( XFEM_STEPPED, false ) ;

	if( dtree != nullptr )
	{
		std::vector<DelaunayTriangle *> elements = layer2d[g]->getElements() ;

		std::pair<Vector , Vector > grad_flux( Vector( 0., elements[0]->getBoundingPoints().size() * 2 * elements.size() ), Vector( 0., elements[0]->getBoundingPoints().size() * 2 * elements.size() ) ) ;

		for( size_t i  = 0 ; i < elements.size() ; i++ )
		{
			if( elements[i]->getBehaviour() && elements[i]->getBehaviour()->type != VOID_BEHAVIOUR )
			{
// 				std::valarray<Point *> pts(3) ;
// 				pts[0] =  elements[i]->first ;
// 				pts[1] =  elements[i]->second ;
// 				pts[2] =  elements[i]->third ;

				Vector gradient(0., 2*elements[i]->getBoundingPoints().size()) ;
				Vector flux(0., 2*elements[i]->getBoundingPoints().size()) ;
				elements[i]->getState().getField( GRADIENT_FIELD, FLUX_FIELD, elements[i]->getBoundingPoints(), gradient, flux, false) ;

				for( size_t j = 0 ; j < elements[0]->getBoundingPoints().size() * 2 ; j++ )
				{
					grad_flux.first[i * elements[0]->getBoundingPoints().size() * 2 + j] = gradient[j] ;
					grad_flux.second[i * elements[0]->getBoundingPoints().size() * 2 + j] = flux[j] ;
				}

				if( i % 1000 == 0 )
					std::cerr << "\r computing gradient+flux... element " << i + 1 << "/" << elements.size() << std::flush ;
			}
		}

		std::cerr << " ...done." << std::endl ;
		return grad_flux ;
	}
	else
	{
		std::vector<DelaunayTetrahedron *> tets = dtree3D->getElements() ;

		if( g != -1 && layer3d.find(g) != layer3d.end() )
			tets = layer3d[g]->getElements() ;

		size_t npoints = tets[0]->getBoundingPoints().size() ;
		std::pair<Vector , Vector > grad_flux( Vector( 0.f, npoints * 3 * tets.size() ), Vector( 0.f, npoints * 3 * tets.size() ) ) ;

		for( size_t i  = 0 ; i < tets.size() ; i++ )
		{

			Vector gradient(0., 3*tets[i]->getBoundingPoints().size()) ;
			Vector flux(0., 3*tets[i]->getBoundingPoints().size()) ;
			tets[i]->getState().getField( GRADIENT_FIELD, FLUX_FIELD, tets[i]->getBoundingPoints(), gradient, flux, false) ;

			for( size_t j = 0 ; j < npoints * 3 ; j++ )
			{
				grad_flux.first[i * npoints * 3 + j] = gradient[j] ;
				grad_flux.second[i * npoints * 3 + j] = flux[j] ;
			}

			if( i % 1000 == 0 )
				std::cerr << "\r computing gradient+flux... element " << i + 1 << "/" << tets.size() << std::flush ;

//				std::cout << grflx.first.size() << std::endl ;
		}

		std::cerr << " ...done." << std::endl ;
		return grad_flux ;
	}
}

std::pair<Vector , Vector > FeatureTree::getGradientAndFlux( const std::vector<DelaunayTetrahedron *> & tets , bool stepTree)
{
	if(stepTree)
		state.setStateTo( XFEM_STEPPED, false ) ;
	std::pair<Vector , Vector > stress_strain( Vector( 4 * 3 * tets.size() ), Vector( 4 * 3 * tets.size() ) ) ;

	for( size_t i  = 0 ; i < tets.size() ; i++ )
	{
		std::valarray<Point *> pts( 4 ) ;
		pts[0] =  tets[i]->first ;
		pts[1] =  tets[i]->second ;
		pts[2] =  tets[i]->third ;
		pts[3] =  tets[i]->fourth ;

		Vector gradient(0., 12) ;
		Vector flux(0., 12) ;
		tets[i]->getState().getField( GRADIENT_FIELD, FLUX_FIELD, tets[i]->getBoundingPoints(), gradient, flux, false) ;
		
		for( size_t j = 0 ; j < 4 ; j++ )
		{
			for( size_t k = 0 ; k < 3 ; k++ )
			{
				stress_strain.first[i * 4 * 3 + j * 3 + k] = gradient[j * 3 + k] ;
				stress_strain.second[i * 4 * 3 + j * 3 + k] = flux[j * 3 + k] ;
			}
		}

//			if(i%1000 == 0)
//				std::cerr << "\r computing gradient+flux... element " << i+1 << "/" << tets.size() << std::flush ;
	}

	std::cerr << " ...done." << std::endl ;
	return stress_strain ;
}

std::pair<Vector , Vector > FeatureTree::getStressAndStrain( const std::vector<DelaunayTetrahedron *> & tets, bool stepTree )
{
	if(stepTree)
		state.setStateTo( XFEM_STEPPED, false ) ;
	std::pair<Vector , Vector > stress_strain( Vector( tets[0]->getBoundingPoints().size() * 6 * tets.size() ), Vector( tets[0]->getBoundingPoints().size() * 6 * tets.size() ) ) ;

	for( size_t i  = 0 ; i < tets.size() ; i++ )
	{
		std::valarray<Point *> pts( 4 ) ;
		pts[0] =  tets[i]->first ;
		pts[1] =  tets[i]->second ;
		pts[2] =  tets[i]->third ;
		pts[3] =  tets[i]->fourth ;

		Vector strain(0., tets[0]->getBoundingPoints().size() * 6) ;
		Vector stress(0., tets[0]->getBoundingPoints().size() * 6) ;
		tets[i]->getState().getField( STRAIN_FIELD, REAL_STRESS_FIELD, tets[i]->getBoundingPoints(), strain, stress, false) ;

		for( size_t j = 0 ; j < tets[i]->getBoundingPoints().size() ; j++ )
		{
			for( size_t k = 0 ; k < 6 ; k++ )
			{
				stress_strain.first[i * tets[0]->getBoundingPoints().size() * 6 + j * 6 + k] = stress[j * 6 + k] ;
				stress_strain.second[i * tets[0]->getBoundingPoints().size() * 6 + j * 6 + k] = strain[j * 6 + k] ;
			}
		}

		if( i % 1000 == 0 )
			std::cerr << "\r computing strain+stress... element " << i + 1 << "/" << tets.size() << std::flush ;
	}

	std::cerr << " ...done." << std::endl ;
	return stress_strain ;
}

Vector FeatureTree::strainFromDisplacements()
{
	state.setStateTo( XFEM_STEPPED, false ) ;

	if( dtree != nullptr )
	{
		std::vector<DelaunayTriangle *> elements = dtree->getElements() ;
		Vector strain( 0.f, 3 * 3 * elements.size() ) ;

		for( size_t i  = 0 ; i < elements.size() ; i++ )
		{
			if( elements[i]->getBehaviour() && elements[i]->getBehaviour()->type != VOID_BEHAVIOUR )
			{
				std::valarray<Point *> pts( 3 ) ;
				pts[0] =  elements[i]->first ;
				pts[1] =  elements[i]->second ;
				pts[2] =  elements[i]->third ;

				Vector str(0., 9) ;
				elements[i]->getState().getField( STRAIN_FIELD, pts, str, false) ;

				for( size_t j = 0 ; j < 9 ; j++ )
					strain[i * 3 * 3 + j] = str[j] ;

				std::cerr << "\r computing strain... element " << i + 1 << "/" << elements.size() << std::flush ;
			}
		}

		std::cerr << " ...done." << std::endl ;
		return strain ;
	}
	else
	{
		std::vector<DelaunayTetrahedron *> elements3D = dtree3D->getElements() ;
		Vector strain( 0., 4 * 6 * elements3D.size() ) ;

		for( size_t i  = 0 ; i < elements3D.size() ; i++ )
		{
			std::valarray<Point *>  pts( 4 ) ;
			pts[0] =  elements3D[i]->first ;
			pts[1] =  elements3D[i]->second ;
			pts[2] =  elements3D[i]->third ;
			pts[3] =  elements3D[i]->fourth ;

			Vector str(0., 24) ;
			elements3D[i]->getState().getField( STRAIN_FIELD, pts, str, false) ;

			for( size_t j = 0 ; j < 24 ; j++ )
				strain[i * 4 * 6 + j] = str[j] ;

			std::cerr << "\r computing strain... element " << i + 1 << "/" << elements3D.size() << std::flush ;
		}

		std::cerr << " ...done." << std::endl ;
		return strain ;
	}

}

Assembly *FeatureTree::getAssembly(bool forceReassembly)
{
	if(forceReassembly)
	{
		state.setStateTo( ASSEMBLED, false ) ;
	}
	return K ;
}

void FeatureTree::insert( Point *p )
{
	double pointDensity = 0 ;

	if( is2D() )
		pointDensity = .6 * sqrt( tree[0]->area() ) / meshPoints.size() ;
	else
		pointDensity = .6 * pow( tree[0]->volume(), 1. / 3. ) / meshPoints.size() ;

	Feature *mother = nullptr;

	for( size_t i  = 0 ; i < this->tree.size() ; i ++ )
	{
		if( this->tree[i]->in( ( *p ) ) )
		{
			bool yes = true ;

			for( size_t k  =  0 ; k < this->tree[i]->getChildren().size() && yes; k++ )
			{
				if( this->tree[i]->getChild( k )->in( ( *p ) ) )
					yes = false ;
			}

			if( yes )
			{
				mother = this->tree[i] ;
				break ;
			}
		}
	}

	bool yes = true ;

	for( size_t k  =  0 ; k <  mother->getChildren().size() && yes; k++ )
	{
		if( mother->getChild( k )->inBoundary( *p, pointDensity ) )
			yes = false ;
	}

	if( yes )
	{
		this->meshPoints.push_back( std::pair<Point *, Feature *>( p, mother ) ) ;

		if( dtree != nullptr )
			this->dtree->insert( p ) ;
		else
			this->dtree3D->insert( p ) ;
	}
}

bool FeatureTree::solverConverged() const
{
	return solverConvergence ;
}

bool FeatureTree::behaviourChanged() const
{
	return behaviourChange ;
}

bool FeatureTree::enrichmentChanged() const
{
	return enrichmentChange ;
}

void FeatureTree::forceEnrichmentChange()
{
	enrichmentChange = true ;
}

void FeatureTree::elasticStep()
{
	Vector lastx( K->getDisplacements() ) ;
	this->K->clear() ;
	assemble() ;
	solve() ;
	bool prevElastic = elastic ;
	elastic = true ;
	stepElements() ;
	elastic = prevElastic ;
}

void FeatureTree::solve()
{
	Vector lastx( K->getDisplacements() ) ;

	if( enrichmentChange || needMeshing)
	{
		K->clear() ;
	}
	
	timeval time0, time1 ;
	gettimeofday( &time0, nullptr );

	if( dtree )
	{
		K->initialiseElementaryMatrices(father2D);
		std::vector<DelaunayTriangle *> elements ;

		std::cerr << "finding nodes for boundary conditions... " << std::flush ;
		for(auto j = layer2d.begin() ; j != layer2d.end() ;j++)
		{
			std::vector<DelaunayTriangle *> elementstmp = j->second->getElements() ;
			elements.insert(elements.end(), elementstmp.begin(), elementstmp.end()) ;
		}
		for( size_t i = 0 ; i < elements.size() ; ++i )
		{
			elements[i]->applyBoundaryCondition( K ) ;
		}
	}
	else
	{
		K->initialiseElementaryMatrices(father3D);
		std::vector<DelaunayTetrahedron *> elements = dtree3D->getElements() ;

		std::cerr << "finding nodes for boundary conditions... " << std::flush ;
		for( size_t i = 0 ; i < elements.size() ; ++i )
		{
			elements[i]->applyBoundaryCondition( K ) ;
		}
	}

	gettimeofday( &time1, nullptr );
	double delta = time1.tv_sec * 1000000 - time0.tv_sec * 1000000 + time1.tv_usec - time0.tv_usec ;
	std::cerr << "...done. Time (s) " << delta / 1e6 << std::endl ;

	gettimeofday( &time0, nullptr );
	
	std::cerr << "Applying coundary conditions... " << std::flush ;
	for( size_t i = 0 ; i < boundaryCondition.size() ; ++i )
	{
		if(i%20 == 0)
			std::cerr << "\rApplying coundary conditions... " << i+1 << "/" << boundaryCondition.size() << std::flush ;
		
		if( dtree )
		{
			boundaryCondition[i]->apply( K, dtree ) ;
		}

		if( dtree3D )
		{
			boundaryCondition[i]->apply( K, dtree3D ) ;
		}
	}
	gettimeofday( &time1, nullptr );
	delta = time1.tv_sec * 1000000 - time0.tv_sec * 1000000 + time1.tv_usec - time0.tv_usec ;
	std::cerr << "\rApplying coundary conditions... " << boundaryCondition.size() << "/" << boundaryCondition.size() << "...done  Time (s) " << delta / 1e6 << std::endl  ;

	needAssembly = true ;

	if( solverConvergence || reuseDisplacements )
	{

		solverConvergence = K->cgsolve( lastx ) ;

		Vector r = K->getMatrix() * K->getDisplacements() - K->getForces() ;
		double perror = residualError ;
		residualError = sqrt( parallel_inner_product( &r[0], &r[0], r.size() ) ) ;

		if( perror > residualError || solverConvergence )
			reuseDisplacements = true;
	}
	else
	{
		lastx = 0 ;


		solverConvergence = K->cgsolve() ;

// 		dtree->project(coarseTrees[3], K->getDisplacements(), coarseAssemblies[3]->getDisplacements(), false) ;
		Vector r = K->getMatrix() * K->getDisplacements() - K->getForces() ;
		double perror = residualError ;
		residualError = sqrt( parallel_inner_product( &r[0], &r[0], r.size() ) ) ;

		if( perror > residualError || solverConvergence )
			reuseDisplacements = true;
	}

}

void FeatureTree::stepXfem()
{
	enrichmentChange = false ;

	if( solverConvergence )
	{
		if( is2D() )
		{
			std::vector<DelaunayTriangle *> elements = dtree->getElements() ;
			
// 			#pragma omp parallel for schedule(runtime)
			for( size_t i = 0 ; i < tree.size() ; i++ )
			{
				if( tree[i]->isEnrichmentFeature )
				{
					
					dynamic_cast<EnrichmentFeature *>( tree[i] )->step( deltaTime, &K->getForces(), dtree ) ;
					bool moved = dynamic_cast<EnrichmentFeature *>( tree[i] )->moved() ;
					enrichmentChange = enrichmentChange || moved;

					if( moved )
					{
						reuseDisplacements = false ;
						needAssembly = true ;
					}

					
				}
				else if( tree[i]->isUpdated )
				{
// 					tree[i]->print() ;
// 					std::cout << "update ! " << std::endl ;
					needAssembly = true ;
					needMeshing = true ;
					reuseDisplacements = false ;
				}
			}

		}
		else if( is3D() )
		{

// 			#pragma omp parallel for schedule(runtime)
			for( size_t i = 0 ; i < tree.size() ; i++ )
			{
				if( tree[i]->isEnrichmentFeature )
				{
					dynamic_cast<EnrichmentFeature *>( tree[i] )->step( deltaTime, &K->getForces(), dtree ) ;
					bool moved =
					    enrichmentChange = enrichmentChange || dynamic_cast<EnrichmentFeature *>( tree[i] )->moved() ;

					if( enrichmentChange )
						needAssembly = true ;

					if( moved )
						reuseDisplacements = false ;
				}
				else if( tree[i]->isUpdated )
				{
					needAssembly = true ;
					needMeshing = true ;
					reuseDisplacements = false ;
				}
			}
		}
	}

	if(enrichmentChange)
		residualError = 1e9 ;
}

bool sortByScore( DelaunayTriangle * tri1, DelaunayTriangle * tri2)
{
	if(tri1->getBehaviour()->getFractureCriterion() && tri2->getBehaviour()->getFractureCriterion())
		return tri1->getBehaviour()->getFractureCriterion()->getNonLocalScoreAtState() > tri2->getBehaviour()->getFractureCriterion()->getNonLocalScoreAtState();
	return false ;
}

bool FeatureTree::stepElements()
{
	behaviourChange = false ;
	needAssembly = false ;
	stateConverged = false ;
	double maxScore = -1 ;
	double maxTolerance = 0 ;
	if( solverConvergence )
	{
		if( is2D() )
		{

			std::vector<DelaunayTriangle *> elements  ;
			for(auto j = layer2d.begin() ; j != layer2d.end() ;j++)
			{
				std::vector<DelaunayTriangle *> elementstmp = j->second->getElements() ;
				for( size_t i = 0 ; i < elementstmp.size() ; i++ )
				{
					if(elementstmp[i]->getBehaviour())
						elements.push_back(elementstmp[i]);
				}
			}
			if(cachedVolumes.empty())
			{
				for( size_t i = 0 ; i < elements.size() ; i++ )
				{
					if(elements[i]->getBehaviour() && elements[i]->getBehaviour()->type != VOID_BEHAVIOUR)
						cachedVolumes.push_back(elements[i]->area()) ;
					else
						cachedVolumes.push_back(0.) ;
				}
			}
			double volume = std::accumulate(cachedVolumes.begin(), cachedVolumes.end(), double(0)) ;
			double previousAverageDamage = averageDamage ;
			double adamage = 0 ;
			if(!elastic)
			{
				crackedVolume = 0 ;
				damagedVolume = 0 ;
				averageDamage = 0. ;
			}
			//this will update the state of all elements. This is necessary as
			//the behaviour updates might depend on the global state of the
			//simulation.
 			std::cerr << " stepping through elements... " << std::flush ;
//#pragma omp parallel for schedule(runtime)
			for( size_t i = 0 ; i < elements.size() ; i++ )
			{
				if( i % 1000 == 0 )
					std::cerr << "\r stepping through elements... " << i << "/" << elements.size() << std::flush ;
				
				elements[i]->step( deltaTime, &K->getDisplacements() ) ;
			}

			std::cerr << " ...done" << std::endl ;

			int fracturedCount = 0 ;
			int ccount = 0 ;
			size_t changecount = 0 ;				

			if( !elastic )
			{
				double maxScoreInit = -1;
				for( size_t i = 0 ; i < elements.size() ; i++ )
				{
					if( i % 500 == 0 )
						std::cerr << "\r checking for fractures (1)... " << i << "/" << elements.size() << std::flush ;

					if( elements[i]->getBehaviour()->getFractureCriterion() )
					{
						elements[i]->getBehaviour()->getFractureCriterion()->step( elements[i]->getState() ) ;
						elements[i]->getBehaviour()->getFractureCriterion()->computeNonLocalState( elements[i]->getState(), NULL_SMOOTH ) ;
						maxScoreInit = std::max(elements[i]->getBehaviour()->getFractureCriterion()->getNonLocalScoreAtState(), maxScoreInit) ;
						
					}
				}
				
				std::cerr << " ...done. " << std::endl ;

// 				std::stable_sort(elements.begin(), elements.end(), sortByScore) ;
				
// #pragma omp parallel for reduction(+:volume,adamage) 
				for( size_t i = 0 ; i < elements.size() ; i++ )
				{

					double are = cachedVolumes[i] ;

					if( i % 10000 == 0 )
						std::cerr << "\r checking for fractures (2)... " << i << "/" << elements.size() << std::flush ;

					if( elements[i]->getBehaviour()->type != VOID_BEHAVIOUR )
					{
						DamageModel * dmodel = elements[i]->getBehaviour()->getDamageModel() ;
						bool wasFractured = elements[i]->getBehaviour()->fractured() ;
						
						elements[i]->getBehaviour()->step( deltaTime, elements[i]->getState(), maxScoreInit ) ;
						if( dmodel )
						{
							if( !elements[i]->getBehaviour()->fractured() )
							{
								adamage += are  * (dmodel->getState().max() > 0.) ;
// 								std::cout << dmodel->getState()[0] << " " << dmodel->getState()[1]  << " " << dmodel->getState()[2] << " " << dmodel->getState()[3] << std::endl ;
// 								std::cout << are << " * " << dmodel->getState().max() << std::endl ;
							}
							else
								adamage += are ;// * dmodel->getState().max() ;

						}
						if( elements[i]->getBehaviour()->changed() )
						{
							needAssembly = true ;
							behaviourChange = true ;
							ccount++ ;
						}
						if( elements[i]->getBehaviour()->fractured() )
						{
							fracturedCount++ ;
							crackedVolume += are ;
							
							if(!wasFractured)
							{
								needAssembly = true ;
								behaviourChange = true ;
							}
						}
						else if( dmodel && dmodel->getState().max() > POINT_TOLERANCE_3D )
						{
							damagedVolume += are ;
						}
					}
				}
				averageDamage = adamage/volume ;

				std::cerr << " ...done. " << ccount << " elements changed." << std::endl ;

				for( size_t i = 0 ; i < elements.size() ; i++ )
				{
					if( i % 1000 == 0 )
						std::cerr << "\r checking for fractures (3)... " << i << "/" << elements.size() << std::flush ;

					if( elements[i]->getBehaviour()->getDamageModel() )
					{
						elements[i]->getBehaviour()->getDamageModel()->postProcess() ;
						if(elements[i]->getBehaviour()->changed())
						{
							needAssembly = true ;
							behaviourChange = true ;
						}
					}
				}
				foundCheckPoint = true ;
				for( size_t i = 0 ; i < elements.size() ; i++ )
				{
					if( elements[i]->getBehaviour()->getDamageModel() && !elements[i]->getBehaviour()->getDamageModel()->converged )
					{
						foundCheckPoint = false ;
						maxScore = maxScoreInit ;
						break ;
					}
				}
				
				#pragma omp parallel for
				for( size_t i = 0 ; i < elements.size() ; i++ )
				{
					if( elements[i]->getBehaviour()->getDamageModel())
						elements[i]->getBehaviour()->getFractureCriterion()->setCheckpoint(foundCheckPoint) ;
				}
				
				
				if(!behaviourChange)
				{
					for( size_t i = 0 ; i < elements.size() ; i++ )
					{
						if( elements[i]->getBehaviour()->getFractureCriterion() && elements[i]->getBehaviour()->getFractureCriterion()->met() )
						{
							behaviourChange = true ;
							break ;
						}
					}
				}
			}
			
			if( !elastic && foundCheckPoint )
			{
				needAssembly = true ;
				std::cout << "[" << averageDamage << " ; " << ccount << " ; " <<  std::flush ;
				maxScore = -1. ;
				maxTolerance = 1 ;

				for( size_t i = 0 ; i < elements.size() ; i++ )
				{
					if( elements[i]->getBehaviour()->getFractureCriterion() )
					{
						elements[i]->getBehaviour()->getFractureCriterion()->setCheckpoint( true ) ;
						maxScore = std::max(elements[i]->getBehaviour()->getFractureCriterion()->getNonLocalScoreAtState(), maxScore) ;
						maxTolerance = std::max(elements[i]->getBehaviour()->getFractureCriterion()->getScoreTolerance(), maxTolerance) ;

					}
				}
				
				std::cout << maxScore << "]" << std::flush ;
				if(elements[0]->getOrder() >= LINEAR_TIME_LINEAR && maxScore > 0 && maxScore < 1.-POINT_TOLERANCE_2D)
				{
					std::cerr << "adjusting time step..." << std::endl ;
					double begin = elements[0]->getBoundingPoint(0).getT() ;
					double end = elements[0]->getBoundingPoint( elements[0]->getBoundingPoints().size() -1).getT() ;
					if(maxScore*(end-begin) > minDeltaTime) 
						moveFirstTimePlanes( (1.-maxScore)*(end-begin) , elements) ;
					else if(end - begin > minDeltaTime)
					{
						moveFirstTimePlanes( end-begin-minDeltaTime , elements) ;
					}
					else
					{
						std::cout << "negative time step: setting to 0..." << std::endl ;
						moveFirstTimePlanes( 0., elements) ;
					}	
				}

				if(maxScore > 1.-POINT_TOLERANCE_2D)
				{
					moveFirstTimePlanes( 0., elements) ;
				}
				
				
			}
			else if(!elastic)
			{
					
// 				if(elements[0]->getOrder() >= LINEAR_TIME_LINEAR && maxScore > 0)
// 				{
// 					moveFirstTimePlanes( 0. , elements) ;
// 				}

			    
				#pragma omp parallel for
				for( size_t i = 0 ; i < elements.size() ; i++ )
				{
					if( elements[i]->getBehaviour()->getFractureCriterion() )
						elements[i]->getBehaviour()->getFractureCriterion()->setCheckpoint( false ) ;
				}


			  
			}

/*			Vector inter(12) ;
			inter = 0 ;
			inter[0] = nodes[0]->getT() ;
			Vector stress = getAverageField(REAL_STRESS_FIELD, -1, -1) ;
			Vector strain = getAverageField(STRAIN_FIELD, -1, -1) ;
			inter[3] = strain[0] ; inter[4] = strain[1] ; inter[5] = strain[2] ; 
			inter[6] = stress[0] ; inter[7] = stress[1] ; inter[8] = stress[2] ; 
			inter[9] = damageAreaInAggregates( elements) ;
			inter[10] = damageAreaInPaste( elements ) ;
			inter[11] = averageDamage ;
			intermediateStates.push_back(inter) ;*/

			std::cerr << " ...done. " << std::endl ;
				
			

		}
		else if( is3D() )
		{
			std::vector<DelaunayTetrahedron *> elements = dtree3D->getElements() ;

			if(cachedVolumes.empty())
			{
				for( size_t i = 0 ; i < elements.size() ; i++ )
				{
					if(elements[i]->getBehaviour() && elements[i]->getBehaviour()->type != VOID_BEHAVIOUR)
						cachedVolumes.push_back(elements[i]->area()) ;
					else
						cachedVolumes.push_back(0.) ;
				}
			}
			double volume = std::accumulate(cachedVolumes.begin(), cachedVolumes.end(), double(0)) ;
			double previousAverageDamage = averageDamage ;
			double adamage = 0 ;
			if(!elastic)
			{
				crackedVolume = 0 ;
				damagedVolume = 0 ;
				averageDamage = 0. ;
			}
			//this will update the state of all elements. This is necessary as
			//the behaviour updates might depend on the global state of the
			//simulation.
 			std::cerr << " stepping through elements... " << std::flush ;
			#pragma omp parallel for
			for( size_t i = 0 ; i < elements.size() ; i++ )
			{
				if( i % 1000 == 0 )
					std::cerr << "\r stepping through elements... " << i << "/" << elements.size() << std::flush ;
				
				elements[i]->step( deltaTime, &K->getDisplacements() ) ;
			}

			std::cerr << " ...done" << std::endl ;

			int fracturedCount = 0 ;
			int ccount = 0 ;

			if( !elastic )
			{
				double maxScoreInit = -1;
				for( size_t i = 0 ; i < elements.size() ; i++ )
				{
					if( i % 500 == 0 )
						std::cerr << "\r checking for fractures (1)... " << i << "/" << elements.size() << std::flush ;

					if( elements[i]->getBehaviour()->getFractureCriterion() )
					{
						elements[i]->getBehaviour()->getFractureCriterion()->step( elements[i]->getState() ) ;
						elements[i]->getBehaviour()->getFractureCriterion()->computeNonLocalState( elements[i]->getState(), NULL_SMOOTH ) ;
						maxScoreInit = std::max(elements[i]->getBehaviour()->getFractureCriterion()->getNonLocalScoreAtState(), maxScoreInit) ;
						
					}
				}
				
				std::cerr << " ...done. " << std::endl ;

// 				std::stable_sort(elements.begin(), elements.end(), sortByScore) ;
				
// #pragma omp parallel for reduction(+:volume,adamage) 
				for( size_t i = 0 ; i < elements.size() ; i++ )
				{

					double are = cachedVolumes[i] ;

					if( i % 10000 == 0 )
						std::cerr << "\r checking for fractures (2)... " << i << "/" << elements.size() << std::flush ;

					if( elements[i]->getBehaviour()->type != VOID_BEHAVIOUR )
					{
						DamageModel * dmodel = elements[i]->getBehaviour()->getDamageModel() ;
						bool wasFractured = elements[i]->getBehaviour()->fractured() ;
						
						elements[i]->getBehaviour()->step( deltaTime, elements[i]->getState(), maxScoreInit ) ;
						if( dmodel )
						{
							if( !elements[i]->getBehaviour()->fractured() )
							{
								adamage += are * dmodel->getState().max() ;
// 								std::cout << dmodel->getState()[0] << " " << dmodel->getState()[1]  << " " << dmodel->getState()[2] << " " << dmodel->getState()[3] << std::endl ;
// 								std::cout << are << " * " << dmodel->getState().max() << std::endl ;
							}
						}
						if( elements[i]->getBehaviour()->changed() )
						{
							needAssembly = true ;
							behaviourChange = true ;
							ccount++ ;
						}
						if( elements[i]->getBehaviour()->fractured() )
						{
							fracturedCount++ ;
							crackedVolume += are ;
							
							if(!wasFractured)
							{
								needAssembly = true ;
								behaviourChange = true ;
							}
						}
						else if( dmodel && dmodel->getState().max() > POINT_TOLERANCE_3D )
						{
							damagedVolume += are ;
						}
					}
				}
				averageDamage = adamage/volume ;

				std::cerr << " ...done. " << ccount << " elements changed." << std::endl ;
				
				#pragma omp parallel for
				for( size_t i = 0 ; i < elements.size() ; i++ )
				{
					if( i % 1000 == 0 )
						std::cerr << "\r checking for fractures (3)... " << i << "/" << elements.size() << std::flush ;

					if( elements[i]->getBehaviour()->getDamageModel() )
					{
						elements[i]->getBehaviour()->getDamageModel()->postProcess() ;
						if(elements[i]->getBehaviour()->changed())
						{
							#pragma omp critical
							{
								needAssembly = true ;
								behaviourChange = true ;
							}
						}
					}
				}
				foundCheckPoint = true ;
				for( size_t i = 0 ; i < elements.size() ; i++ )
				{
					if( elements[i]->getBehaviour()->getDamageModel() && !elements[i]->getBehaviour()->getDamageModel()->converged )
					{
						foundCheckPoint = false ;
						maxScore = maxScoreInit ;
						break ;
					}
				}
				
				if(!behaviourChange)
				{
					for( size_t i = 0 ; i < elements.size() ; i++ )
					{
						if( elements[i]->getBehaviour()->getFractureCriterion() && elements[i]->getBehaviour()->getFractureCriterion()->met() )
						{
							behaviourChange = true ;
							break ;
						}
					}
				}
			}
			
			if( !elastic && foundCheckPoint )
			{
				std::cout << "[" << averageDamage << " ; " << std::flush ;
				maxScore = -1. ;
				maxTolerance = 1 ;
// 				double maxs = -1 ;
// 				double maxtol = 1 ;
// 				#pragma omp parallel lastshared(maxs,maxtol)
// 				{
// 
// 				  #pragma omp for nowait 
				  for( size_t i = 0 ; i < elements.size() ; i++ )
				  {
					  if( elements[i]->getBehaviour()->getFractureCriterion() )
					  {
						  //std::cout << "." << std::flush ;
						  elements[i]->getBehaviour()->getFractureCriterion()->setCheckpoint( true ) ;
						  maxScore = std::max(elements[i]->getBehaviour()->getFractureCriterion()->getNonLocalScoreAtState(), maxScore) ;
						  maxTolerance = std::max(elements[i]->getBehaviour()->getFractureCriterion()->getScoreTolerance(), maxTolerance) ;
	
					  }
				  }
// 				  #pragma omp critical
// 				  {
// 				    maxScore = std::max(maxScore, maxs) ;
// 				    maxTolerance = std::max(maxTolerance, maxtol) ;
// 				  }
// 				}
				
				std::cout << maxScore << "]" << std::flush ;
				if(elements[0]->getOrder() >= LINEAR_TIME_LINEAR && maxScore > 0. && maxScore < 1.)
				{
					std::cerr << "adjusting time step..." << std::endl ;
					double begin = elements[0]->getBoundingPoint(0).getT() ;
					double end = elements[0]->getBoundingPoint( elements[0]->getBoundingPoints().size() -1).getT() ;
					if(maxScore*(end-begin) > minDeltaTime) 
						moveFirstTimePlanes( (1.-maxScore)*(end-begin) , elements) ;
					else if(end - begin > minDeltaTime)
					{
						moveFirstTimePlanes( end-begin-minDeltaTime , elements) ;
					}
					else
					{
						std::cout << "negative time step: setting to 0..." << std::endl ;
						this->moveFirstTimePlanes( 0., elements) ;
					}	
				}
				
				
			}
			else if(!elastic)
			{
					
				if(elements[0]->getOrder() >= LINEAR_TIME_LINEAR && maxScore > 0)
				{
					moveFirstTimePlanes( 0. , elements) ;
				}

			    
				#pragma omp parallel for schedule(runtime)
				for( size_t i = 0 ; i < elements.size() ; i++ )
				{
					if( elements[i]->getBehaviour()->getFractureCriterion() )
						elements[i]->getBehaviour()->getFractureCriterion()->setCheckpoint( false ) ;
				}


			  
			}

			std::cerr << " ...done. " << std::endl ;
				
		}
	}
	else
	{
		Vector dummyx( 0., K->getDisplacements().size() ) ;

		if( is2D() )
		{
			std::vector<DelaunayTriangle *> elements ;
			for(auto j = layer2d.begin() ; j != layer2d.end() ;j++)
			{
				std::vector<DelaunayTriangle *> elementstmp = j->second->getElements() ;
				elements.insert(elements.end(), elementstmp.begin(), elementstmp.end()) ;
			}
			double volume = 0;
			crackedVolume = 0 ;
			damagedVolume = 0 ;
			//this will update the state of all elements. This is necessary as
			//the behaviour updates might depend on the global state of the
			//simulation.
			std::cerr << " stepping through elements... " << std::flush ;

			for( size_t i = 0 ; i < elements.size() ; i++ )
			{
				if( i % 1000 == 0 )
					std::cerr << "\r stepping through elements... " << i << "/" << elements.size() << std::flush ;

				elements[i]->step( 0., &dummyx ) ;
			}

			std::cerr << " ...done" << std::endl ;

			for( size_t i = 0 ; i < elements.size() ; i++ )
				elements[i]->clearVisited() ;


		}
		else if( is3D() )
		{
			std::vector<DelaunayTetrahedron *> elements = dtree3D->getElements() ;
			for(auto j = layer3d.begin() ; j != layer3d.end() ;j++)
			{
				std::vector<DelaunayTetrahedron *> elementstmp = j->second->getElements() ;
				elements.insert(elements.end(), elementstmp.begin(), elementstmp.end()) ;
			}
			//this will update the state of all elements. This is necessary as
			//the behaviour updates might depend on the global state of the
			//simulation.

			//this will update the state of all elements. This is necessary as
			//the behaviour updates might depend on the global state of the
			//simulation.
			std::cerr << " stepping through elements... " << std::flush ;

			for( size_t i = 0 ; i < elements.size() ; i++ )
			{
				if( i % 1000 == 0 )
					std::cerr << "\r stepping through elements... " << i << "/" << elements.size() << std::flush ;

				elements[i]->step( 0., &dummyx ) ;
			}

			std::cerr << " ...done" << std::endl ;

			// 		std::cout << " Fractured " << fracturedCount << " Elements" << std::endl ;
			// 		std::cout << " Fractured Fraction " <<  crackedVolume / volume << std::endl ;

		}
	}

	stateConverged = foundCheckPoint && maxScore < maxTolerance ;
	if(behaviourChange)
		residualError = 1e9 ;
	return foundCheckPoint && maxScore < maxTolerance;
}


void FeatureTree::State::setStateTo( StateType s, bool stepChanged )
{
	bool behaviourChanged = ft->behaviourChanged() ;
	bool xfemChanged = ft->enrichmentChanged() ;
	bool samplingChanged = ft->needMeshing ;
	bool initialiseFractureCache = false ;

	if( samplingChanged )
	{
		initialiseFractureCache = true ;
		sampled = false ;
		meshed = false ;
		behaviourSet = behaviourSet ;
		behaviourUpdated = false;
		stitched = false ;
		renumbered = false ;
		initialised = false ;
		enriched = false ;
		assembled = false ;
		solved = false ;
		behaviourStepped = false;
		xfemStepped = false ;
		featureStepped = false;
	}

	if( stepChanged )
	{
		if( ft->deltaTime > POINT_TOLERANCE_2D )
			enriched = false ;

		assembled = false ;
		solved = false ;
		behaviourStepped = false;
		xfemStepped = false ;
		featureStepped = false;
	}

	if( xfemChanged )
	{
		initialised = false ;
		enriched = false ;
		assembled = false ;
		solved = false ;
		behaviourStepped = false;
		xfemStepped = false ;
		featureStepped = false;
	}

	if( behaviourChanged )
	{
		assembled = false ;
		solved = false ;
		behaviourStepped = false;
		xfemStepped = false ;
		featureStepped = false;
	}

// 			std::cout << 				sampled <<
// 				meshed <<
// 				behaviourSet <<
// 				behaviourUpdated <<
// 				stitched <<
// 				renumbered <<
// 				initialised <<
// 				enriched <<
// 				assembled <<
// 				solved <<
// 				behaviourStepped <<
// 				xfemStepped <<
// 				featureStepped <<  std::endl ;

#ifdef HAVE_OMP
	if( !sampled )
	{
		double t0 = omp_get_wtime() ;
		ft->sample();
		sampled = true ;
		std::cerr << "... time to sample (s) " << omp_get_wtime()-t0 << std::endl ;
	}
#else
	if( !sampled )
	{
		ft->sample();
		sampled = true ;
	}
#endif
	
	if( s == SAMPLED )
		return ;

#ifdef HAVE_OMP
	if( !meshed )
	{
		double t0 = omp_get_wtime() ;
		ft->generateElements();
		meshed = true ;
		std::cerr << "... time to mesh (s) " << omp_get_wtime()-t0 << std::endl ;
	}
#else
	if( !meshed )
	{
		ft->generateElements();
		meshed = true ;
	}
#endif
	
	if( s == MESHED )
		return ;

	
	if( !stitched )
	{
		ft->stitch();
		stitched = true ;
	}
	if( s == STITCHED )
		return ;
	
	if( !behaviourSet )
	{
		ft->setElementBehaviours() ;
		behaviourSet = true;
		behaviourUpdated = true;
		initialised = false ;
	}
	else if( !behaviourUpdated && behaviourSet )
	{
		ft->updateElementBehaviours();
		behaviourUpdated = true;
	}
	if( s == BEHAVIOUR_SET )
		return ;
	
	if( !renumbered )
	{
		ft->renumber();
		renumbered = true ;
	}
	if( s == RENUMBERED )
		return ;

	if( !initialised )
	{
		ft->initializeElements( initialiseFractureCache );
		initialised = true ;
	}
	if( s == INITIALISED )
		return ;

	if( !enriched )
	{
		ft->enrich() ;
		enriched = true ;
	}
	if( s == ENRICHED )
		return ;

	if( !assembled )
	{
		ft->assemble();
		assembled = true ;
	}
	if( s == ASSEMBLED )
		return ;

	if( !solved )
	{
		ft->solve() ;
		solved = true ;
	}
	if( s == SOLVED )
		return ;

	if( !behaviourStepped )
	{
		ft->stepElements();
		behaviourStepped = true ;
	}
	if( s == BEHAVIOUR_STEPPED )
		return ;
	
	if( !xfemStepped )
	{
		ft->stepXfem();
		xfemStepped = true ;
	}
	if( s == XFEM_STEPPED )
		return ;


	
}

void FeatureTree::resetBoundaryConditions() 
{
	boundaryCondition.clear() ; 
	needAssembly = true ; 
	if(K)
		K->clear() ;
} ;

bool FeatureTree::step()
{
	double realdt = deltaTime ;
	
	if( stateConverged && state.meshed && solverConverged() && !behaviourChanged())
	{
		now += deltaTime ;
 		for(size_t i = 0 ; i < nodes.size() ; i++)
 			nodes[i]->getT() += deltaTime ;
	}
	else
	{
		deltaTime = 0 ;
	}

	for(size_t i = 0 ; i < boundaryCondition.size() ; i++)
	{
		TimeContinuityBoundaryCondition * timec = dynamic_cast<TimeContinuityBoundaryCondition *>(boundaryCondition[i]) ;
		if(timec != nullptr)
			timec->goToNext = stateConverged ;
	}

	bool ret = true ;
	size_t it = 0 ;
	int notConvergedCounts = 0 ;
	int betweenCheckpointCount = 0 ;
	int maxBetweenCheckPoints = 2 ;
	
//	std::cout << it << "/" << maxitPerStep << "." << std::flush ;
	bool needexit = false ;
	do
	{
		if(it == 1)
		{
			for(size_t k = 0 ; k < boundaryCondition.size() ; k++)
			{
				TimeContinuityBoundaryCondition * timec = dynamic_cast<TimeContinuityBoundaryCondition *>(boundaryCondition[k]) ;
				if(timec != nullptr)
					timec->goToNext = stateConverged ;
			}
		}
		state.setStateTo( XFEM_STEPPED, true ) ;

 		deltaTime = 0 ;
		if( solverConverged() )
		{
			std::cout << "." << std::flush ;
			notConvergedCounts = 0 ;
			if(!foundCheckPoint)
				betweenCheckpointCount++ ;
			maxBetweenCheckPoints = std::max(betweenCheckpointCount, maxBetweenCheckPoints) ;
		}
		else
		{
			notConvergedCounts++ ;
			std::cout << "+" << std::flush ;
		}

		if( enrichmentChange || needMeshing )
		{
			K->clear() ;
		}
		
		if(needexit && foundCheckPoint && it%(maxBetweenCheckPoints-1) == 0)
			break ;
		
		if( it > maxitPerStep && foundCheckPoint)
		{
			ret = false ;
			needexit = true ;
		}
		++it ;
		
	} while ( ( behaviourChanged() || !solverConverged() || enrichmentChange) && 
	         !( !solverConverged() && !reuseDisplacements ) && 
	          ( notConvergedCounts < 20 ) 
					) ;

	if(notConvergedCounts >= 20)
		ret = false ;
	std::cout << std::endl ;
	if(ret)
		setDeltaTime(realDeltaTime) ;
	std::cout << it << "/" << maxitPerStep << "." << std::flush ;
	damageConverged = solverConverged() && !behaviourChanged() /*stateConverged*/ && ret && (it < maxitPerStep) ;	

	return solverConverged() && !behaviourChanged() /*stateConverged*/ && ret ;
}

bool FeatureTree::stepToCheckPoint()
{
	scaleBoundaryConditions(1);
	double realdt = deltaTime ;

	
	if( solverConverged() && !behaviourChanged() )
		now += deltaTime ;
	else
		deltaTime = 0 ;

	if( enrichmentChange || needMeshing )
	{
		K->clear() ;
	}

	state.setStateTo( XFEM_STEPPED, true ) ;
	int notConvergedCounts = 0 ;
	
	do
	{
		deltaTime = 0 ;
		if( solverConverged() )
		{
			std::cout << "." << std::flush ;
			notConvergedCounts = 0 ;
		}
		else
		{
			notConvergedCounts++ ;
			std::cout << "+" << std::flush ;
		}

		if( enrichmentChange || needMeshing )
		{
			K->clear() ;
		}

		state.setStateTo( XFEM_STEPPED, true ) ;

	}while ( !foundCheckPoint && ( behaviourChanged() || !solverConverged() )  && !( !solverConverged() && !reuseDisplacements ) && notConvergedCounts < 4 ) ;
	
	if(behaviourChanged())
	{
		double upmultiplier = 1 ;
		double currentmultiplier = 0.5 ;
		double downmultiplier = 0 ;
		scaleBoundaryConditions(currentmultiplier);
		while(std::abs(upmultiplier-downmultiplier) > 1./pow(2, 16) )
		{
			if(!isStable())
			{
				upmultiplier = currentmultiplier ;
				currentmultiplier = (upmultiplier+downmultiplier)*.5 ;
			}
			else
			{
				downmultiplier = currentmultiplier ;
				currentmultiplier = (upmultiplier+downmultiplier)*.5 ;
			}
			scaleBoundaryConditions(currentmultiplier);
//			std::cout << currentmultiplier << std::endl ;
		}
		
		scaleBoundaryConditions(downmultiplier);
		deltaTime = realdt ;
		elasticStep();
// 		state.setStateTo( XFEM_STEPPED, true ) ;
		if( solverConverged() )
		{
			std::cout << ":" << std::flush ;
			notConvergedCounts = 0 ;
		}
		else
		{
			notConvergedCounts++ ;
			std::cout << ";" << std::flush ;
		}
		scaleBoundaryConditions(1);
	}

	std::cout  << std::endl ;
	setDeltaTime(realdt) ;
	return solverConverged();
}

bool orderPointsByID( Point * p1, Point * p2)
{
	return p1->getId() < p2->getId() ;
}

std::vector<Point *> FeatureTree::getNodes(int grid) 
{
	if(nodes.size() > 0)
		return nodes ;
  
	std::vector<Point *> pts ;
	if(is2D())
	{
		std::vector<DelaunayTriangle *> elements = this->getElements2D( grid ) ;
		std::valarray<bool> done(elements.size()*elements[0]->getBoundingPoints().size()) ;
		done = false ;
		for(size_t i = 0 ; i < elements.size() ; i++)
		{
			for(size_t j = 0 ; j < elements[i]->getBoundingPoints().size() ; j++)
			{
				if(!done[ elements[i]->getBoundingPoint(j).getId()])
				{
					done[ elements[i]->getBoundingPoint(j).getId()] = true ;
					pts.push_back(&elements[i]->getBoundingPoint(j)) ;
				}
			} 
		}
	}
	if(is3D())
	{
		std::vector<DelaunayTetrahedron *> elements = this->getElements3D( grid ) ;
		std::valarray<bool> done(elements.size()*elements[0]->getBoundingPoints().size()) ;
		done = false ;
		for(size_t i = 0 ; i < elements.size() ; i++)
		{
			for(size_t j = 0 ; j < elements[i]->getBoundingPoints().size() ; j++)
			{
				if(!done[ elements[i]->getBoundingPoint(j).getId()])
				{
					done[ elements[i]->getBoundingPoint(j).getId()] = true ;
					pts.push_back(&elements[i]->getBoundingPoint(j)) ;
				}
			} 
		}
	}
	std::sort( pts.begin(), pts.end(), orderPointsByID) ;
	return pts ;
}

Vector FeatureTree::getAverageField( FieldType f, int grid , double t) 
{
	Vector avg ;
	Vector buffer ;
	double volume = 0 ;
	VirtualMachine vm ;
	if(is2D())
	{
		std::vector<DelaunayTriangle *> elements = getElements2D( grid ) ;
		size_t blocks = elements[0]->getBehaviour()->getNumberOfDegreesOfFreedom()/2 ;
		avg.resize(fieldTypeElementarySize(f, SPACE_TWO_DIMENSIONAL, blocks)) ; 
		buffer.resize(fieldTypeElementarySize(f, SPACE_TWO_DIMENSIONAL, blocks)) ; 
		avg = 0 ; buffer = 0 ;
		for(size_t i = 0 ; i < elements.size() ; i++)
		{
			if(elements[i]->getBehaviour()->type != VOID_BEHAVIOUR)
			{
				elements[i]->getState().getAverageField( f, buffer,&vm, -1, t) ;

				avg += buffer * elements[i]->area() ;
				volume += elements[i]->area() ;
			}
		}
	}
	else
	{
		std::vector<DelaunayTetrahedron *> elements = getElements3D( grid ) ;
		size_t blocks = elements[0]->getBehaviour()->getNumberOfDegreesOfFreedom()/3 ;
		avg.resize(fieldTypeElementarySize(f, SPACE_THREE_DIMENSIONAL, blocks)) ; 
		buffer.resize(fieldTypeElementarySize(f, SPACE_THREE_DIMENSIONAL, blocks)) ; 
		avg = 0 ; buffer = 0 ;
		for(size_t i = 0 ; i < elements.size() ; i++)
		{
			elements[i]->getState().getAverageField( f, buffer,&vm, -1, t ) ;
			avg += buffer * elements[i]->volume() ;
			volume += elements[i]->volume() ;
		}
	  
	}
	return avg/volume ;
}

std::vector<double>  FeatureTree::getMacroscopicStrain(const Geometry * base, double tol) 
{
	Vector disps = getDisplacements() ;
	std::vector<double> ret ;
	if(!disps.size())
		return ret ;
	if(is2D())
	{
		double size_x = dynamic_cast<const Rectangle *>(base)->width() ;
		double size_y = dynamic_cast<const Rectangle *>(base)->height() ;
		std::vector<DelaunayTriangle *> elements = dtree->getElements() ;
		double dxp = 0;
		double dxpc = 0 ;
		double dyp = 0;
		double dypc = 0 ;
		double dxm = 0;
		double dxmc = 0 ;
		double dym = 0;
		double dymc = 0 ;
		std::set<int> doneIds ;
		for(size_t i = 0 ; i < elements.size() ; i++)
		{
			if(elements[i]->getBehaviour() && elements[i]->getBehaviour()->type != VOID_BEHAVIOUR)
			{
				for(size_t j = 0 ; j < elements[i]->getBoundingPoints().size() ; j++)
				{
					
					Point test(elements[i]->getBoundingPoint(j)) ;
					
					base->project(&test);
					if(dist(test, elements[i]->getBoundingPoint(j)) < tol && std::abs(elements[i]->getBoundingPoint(j).getX() - (base->getCenter().getX()+size_x*.5)) < tol && doneIds.find(2*elements[i]->getBoundingPoint(j).getId()) == doneIds.end())
					{
						dxp += disps[2*elements[i]->getBoundingPoint(j).getId()] ;
						dxpc++ ;
						doneIds.insert(2*elements[i]->getBoundingPoint(j).getId());
					}
					if(dist(test, elements[i]->getBoundingPoint(j)) < tol && std::abs(elements[i]->getBoundingPoint(j).getY() - (base->getCenter().getY()+size_y*.5)) < tol&& doneIds.find(2*elements[i]->getBoundingPoint(j).getId()+1) == doneIds.end())
					{
						dyp += disps[2*elements[i]->getBoundingPoint(j).getId()+1] ;
						doneIds.insert(2*elements[i]->getBoundingPoint(j).getId()+1);
						dypc++ ;
					}
					if(dist(test, elements[i]->getBoundingPoint(j)) < tol && std::abs(elements[i]->getBoundingPoint(j).getX() - (base->getCenter().getX()-size_x*.5)) < tol&& doneIds.find(2*elements[i]->getBoundingPoint(j).getId()) == doneIds.end())
					{
						dxm += disps[2*elements[i]->getBoundingPoint(j).getId()] ;
						doneIds.insert(2*elements[i]->getBoundingPoint(j).getId());
						dxmc++ ;
					}
					if(dist(test, elements[i]->getBoundingPoint(j)) < tol && std::abs(elements[i]->getBoundingPoint(j).getY() - (base->getCenter().getY()-size_y*.5)) < tol&& doneIds.find(2*elements[i]->getBoundingPoint(j).getId()+1) == doneIds.end())
					{
						dym += disps[2*elements[i]->getBoundingPoint(j).getId()+1] ;
						doneIds.insert(2*elements[i]->getBoundingPoint(j).getId()+1);
						dymc++ ;
					}
				}
			}
		}
		
		ret.push_back((dxp/dxpc-dxm/dxmc)/size_x);
		ret.push_back((dyp/dypc-dym/dymc)/size_y);
	}
	else
	{
		double size_x = dynamic_cast<const Hexahedron *>(base)->getXSize() ;
		double size_y = dynamic_cast<const Hexahedron *>(base)->getYSize() ;
		double size_z = dynamic_cast<const Hexahedron *>(base)->getZSize() ;
		std::vector<DelaunayTetrahedron *> elements = dtree3D->getElements() ;
		double dxp = 0;
		double dxpc = 0 ;
		double dyp = 0;
		double dypc = 0 ;
		double dzp = 0;
		double dzpc = 0 ;
		double dxm = 0;
		double dxmc = 0 ;
		double dym = 0;
		double dymc = 0 ;
		double dzm = 0;
		double dzmc = 0 ;
		std::set<int> doneIds ;
		for(size_t i = 0 ; i < elements.size() ; i++)
		{
			if(elements[i]->getBehaviour()->type != VOID_BEHAVIOUR)
			{
				for(size_t j = 0 ; j < elements[i]->getBoundingPoints().size() ; j++)
				{
					Point test(elements[i]->getBoundingPoint(j)) ;
					base->project(&test);
					if(dist(test, elements[i]->getBoundingPoint(j)) < tol && std::abs(elements[i]->getBoundingPoint(j).getX() - (base->getCenter().getX()-size_x*.5)) < tol && doneIds.find(3*elements[i]->getBoundingPoint(j).getId()) == doneIds.end())
					{
						dxp = disps[3*elements[i]->getBoundingPoint(j).getId()] ;
						dxpc++ ;
						doneIds.insert(3*elements[i]->getBoundingPoint(j).getId());
					}
					if(dist(test, elements[i]->getBoundingPoint(j)) < tol && std::abs(elements[i]->getBoundingPoint(j).getY() - (base->getCenter().getY()-size_y*.5)) < tol&& doneIds.find(3*elements[i]->getBoundingPoint(j).getId()+1) == doneIds.end())
					{
						dyp = disps[3*elements[i]->getBoundingPoint(j).getId()+1] ;
						doneIds.insert(3*elements[i]->getBoundingPoint(j).getId()+1);
						dypc++ ;
					}
					if(dist(test, elements[i]->getBoundingPoint(j)) < tol && std::abs(elements[i]->getBoundingPoint(j).getZ() - (base->getCenter().getZ()-size_z*.5)) < tol&& doneIds.find(3*elements[i]->getBoundingPoint(j).getId()+2) == doneIds.end())
					{
						dzp = disps[3*elements[i]->getBoundingPoint(j).getId()+2] ;
						doneIds.insert(3*elements[i]->getBoundingPoint(j).getId()+2);
						dzpc++ ;
					}
					if(dist(test, elements[i]->getBoundingPoint(j)) < tol && std::abs(elements[i]->getBoundingPoint(j).getX() - (base->getCenter().getX()-size_x*.5)) < tol&& doneIds.find(3*elements[i]->getBoundingPoint(j).getId()) == doneIds.end())
					{
						dxm = disps[3*elements[i]->getBoundingPoint(j).getId()] ;
						dxmc++ ;
						doneIds.insert(3*elements[i]->getBoundingPoint(j).getId());
					}
					if(dist(test, elements[i]->getBoundingPoint(j)) < tol && std::abs(elements[i]->getBoundingPoint(j).getY() - (base->getCenter().getY()-size_y*.5)) < tol&& doneIds.find(3*elements[i]->getBoundingPoint(j).getId()+1) == doneIds.end())
					{
						dym = disps[3*elements[i]->getBoundingPoint(j).getId()+1] ;
						doneIds.insert(3*elements[i]->getBoundingPoint(j).getId()+1);
						dymc++ ;
					}
					if(dist(test, elements[i]->getBoundingPoint(j)) < tol && std::abs(elements[i]->getBoundingPoint(j).getY() - (base->getCenter().getZ()-size_z*.5)) < tol&& doneIds.find(3*elements[i]->getBoundingPoint(j).getId()+2) == doneIds.end())
					{
						dzm = disps[3*elements[i]->getBoundingPoint(j).getId()+2] ;
						doneIds.insert(3*elements[i]->getBoundingPoint(j).getId()+2);
						dzmc++ ;
					}
				}
			}
		}
		
		
		ret.push_back((dxp/dxpc-dxm/dxmc)/size_x);
		ret.push_back((dyp/dypc-dym/dymc)/size_y);
		ret.push_back((dzp/dzpc-dzm/dzmc)/size_z);
	}
	
	return ret ;
}

std::vector<double> FeatureTree::getCutMacroscopicStrain(const Geometry * base, double tol, double cut) 
{
	Vector disps = getDisplacements() ;
	std::vector<double> ret ;
	if(!disps.size())
		return ret ;
	if(is2D())
	{
		double size_x = dynamic_cast<const Rectangle *>(base)->width() ;
		double size_y = dynamic_cast<const Rectangle *>(base)->height() ;
		std::vector<DelaunayTriangle *> elements = dtree->getElements() ;
		tol = .05*std::min(size_x, size_y) ;
		std::vector<double> dxp ;
		std::vector<double> dyp ;
		std::vector<double> dxm ;
		std::vector<double> dym ;
		std::set<int> doneIds ;
		for(size_t i = 0 ; i < elements.size() ; i++)
		{
			if(elements[i]->getBehaviour() && elements[i]->getBehaviour()->type != VOID_BEHAVIOUR)
			{
				for(size_t j = 0 ; j < elements[i]->getBoundingPoints().size() ; j++)
				{
					
					Point test(elements[i]->getBoundingPoint(j)) ;
					
					base->project(&test);
					if(dist(test, elements[i]->getBoundingPoint(j)) < tol && std::abs(elements[i]->getBoundingPoint(j).getX() - (base->getCenter().getX()+size_x*.5)) < tol && doneIds.find(2*elements[i]->getBoundingPoint(j).getId()) == doneIds.end())
					{
						dxp.push_back(disps[2*elements[i]->getBoundingPoint(j).getId()]) ;
						doneIds.insert(2*elements[i]->getBoundingPoint(j).getId());
					}
					if(dist(test, elements[i]->getBoundingPoint(j)) < tol && std::abs(elements[i]->getBoundingPoint(j).getY() - (base->getCenter().getY()+size_y*.5)) < tol&& doneIds.find(2*elements[i]->getBoundingPoint(j).getId()+1) == doneIds.end())
					{
						dyp.push_back(disps[2*elements[i]->getBoundingPoint(j).getId()+1]) ;
						doneIds.insert(2*elements[i]->getBoundingPoint(j).getId()+1);
					}
					if(dist(test, elements[i]->getBoundingPoint(j)) < tol && std::abs(elements[i]->getBoundingPoint(j).getX() - (base->getCenter().getX()-size_x*.5)) < tol&& doneIds.find(2*elements[i]->getBoundingPoint(j).getId()) == doneIds.end())
					{
						dxm.push_back(disps[2*elements[i]->getBoundingPoint(j).getId()]) ;
						doneIds.insert(2*elements[i]->getBoundingPoint(j).getId());
					}
					if(dist(test, elements[i]->getBoundingPoint(j)) < tol && std::abs(elements[i]->getBoundingPoint(j).getY() - (base->getCenter().getY()-size_y*.5)) < tol&& doneIds.find(2*elements[i]->getBoundingPoint(j).getId()+1) == doneIds.end())
					{
						dym.push_back(disps[2*elements[i]->getBoundingPoint(j).getId()+1]) ;
						doneIds.insert(2*elements[i]->getBoundingPoint(j).getId()+1);
					}
				}
			}
		}
		std::sort(dxp.begin(), dxp.end()) ;
		std::sort(dxm.begin(), dxm.end()) ;
		std::sort(dyp.begin(), dyp.end()) ;
		std::sort(dym.begin(), dym.end()) ;
		double xpmax = dxp[dxp.size()*cut]; 
		double xpmin = dxp[dxp.size()*(1.-cut)];
		double xmmax = dxm[dxm.size()*cut]; 
		double xmmin = dxm[dxm.size()*(1.-cut)];
		double ypmax = dyp[dyp.size()*cut]; 
		double ypmin = dyp[dyp.size()*(1.-cut)];
		double ymmax = dym[dym.size()*cut]; 
		double ymmin = dym[dym.size()*(1.-cut)];
		
		double dxpv = 0;
		double dxpc = 0;
		double dxmv = 0;
		double dxmc = 0;
		double dypv = 0;
		double dypc = 0;
		double dymv = 0;
		double dymc = 0;
		
		for(size_t i = 0 ; i < dxp.size() ; i++)
		{
			if(dxp[i] > xpmin && dxp[i] < xpmax )
			{
				dxpv += dxp[i] ;
				dxpc++ ;
			}
		}
		for(size_t i = 0 ; i < dxm.size() ; i++)
		{
			if(dxm[i] > xmmin && dxm[i] < xmmax )
			{
				dxmv += dxm[i] ;
				dxmc++ ;
			}
		}

		for(size_t i = 0 ; i < dyp.size() ; i++)
		{
			if(dyp[i] > ypmin && dyp[i] < ypmax )
			{
				dypv += dyp[i] ;
				dypc++ ;
			}
		}
		for(size_t i = 0 ; i < dym.size() ; i++)
		{
			if(dym[i] > ymmin && dym[i] < ymmax )
			{
				dymv += dym[i] ;
				dymc++ ;
			}
		}
		
		ret.push_back((dxpv/dxpc-dxmv/dxmc)/size_x);
		ret.push_back((dypv/dypc-dymv/dymc)/size_y);
	}
	else
	{
		double size_x = dynamic_cast<const Hexahedron *>(base)->getXSize() ;
		double size_y = dynamic_cast<const Hexahedron *>(base)->getYSize() ;
		double size_z = dynamic_cast<const Hexahedron *>(base)->getZSize() ;
		std::vector<DelaunayTetrahedron *> elements = dtree3D->getElements() ;
		std::vector<double> dxp ;
		std::vector<double> dyp ;
		std::vector<double> dzp ;
		std::vector<double> dxm ;
		std::vector<double> dym ;
		std::vector<double> dzm ;
		std::set<int> doneIds ;
		for(size_t i = 0 ; i < elements.size() ; i++)
		{
			if(elements[i]->getBehaviour()->type != VOID_BEHAVIOUR)
			{
				for(size_t j = 0 ; j < elements[i]->getBoundingPoints().size() ; j++)
				{
					Point test(elements[i]->getBoundingPoint(j)) ;
					base->project(&test);
					if(dist(test, elements[i]->getBoundingPoint(j)) < tol && std::abs(elements[i]->getBoundingPoint(j).getX() - (base->getCenter().getX()-size_x*.5)) < tol && doneIds.find(3*elements[i]->getBoundingPoint(j).getId()) == doneIds.end())
					{
						dxp.push_back(disps[3*elements[i]->getBoundingPoint(j).getId()]) ;
					}
					if(dist(test, elements[i]->getBoundingPoint(j)) < tol && std::abs(elements[i]->getBoundingPoint(j).getY() - (base->getCenter().getY()-size_y*.5)) < tol && doneIds.find(3*elements[i]->getBoundingPoint(j).getId()) == doneIds.end())
					{
						dyp.push_back(disps[3*elements[i]->getBoundingPoint(j).getId()+1]) ;
					}
					if(dist(test, elements[i]->getBoundingPoint(j)) < tol && std::abs(elements[i]->getBoundingPoint(j).getZ() - (base->getCenter().getZ()-size_z*.5)) < tol && doneIds.find(3*elements[i]->getBoundingPoint(j).getId()) == doneIds.end())
					{
						dzp.push_back(disps[3*elements[i]->getBoundingPoint(j).getId()+2]) ;
					}
					if(dist(test, elements[i]->getBoundingPoint(j)) < tol && std::abs(elements[i]->getBoundingPoint(j).getX() - (base->getCenter().getX()-size_x*.5)) < tol && doneIds.find(3*elements[i]->getBoundingPoint(j).getId()) == doneIds.end())
					{
						dxm.push_back(disps[3*elements[i]->getBoundingPoint(j).getId()]) ;
					}
					if(dist(test, elements[i]->getBoundingPoint(j)) < tol && std::abs(elements[i]->getBoundingPoint(j).getY() - (base->getCenter().getY()-size_y*.5)) < tol && doneIds.find(3*elements[i]->getBoundingPoint(j).getId()) == doneIds.end())
					{
						dym.push_back(disps[3*elements[i]->getBoundingPoint(j).getId()+1]) ;
					}
					if(dist(test, elements[i]->getBoundingPoint(j)) < tol && std::abs(elements[i]->getBoundingPoint(j).getY() - (base->getCenter().getZ()-size_z*.5)) < tol && doneIds.find(3*elements[i]->getBoundingPoint(j).getId()) == doneIds.end())
					{
						dzm.push_back(disps[3*elements[i]->getBoundingPoint(j).getId()+2]) ;
					}
				}
			}
		}
		
		std::sort(dxp.begin(), dxp.end()) ;
		std::sort(dxm.begin(), dxm.end()) ;
		std::sort(dyp.begin(), dyp.end()) ;
		std::sort(dym.begin(), dym.end()) ;
		std::sort(dzp.begin(), dzp.end()) ;
		std::sort(dzm.begin(), dzm.end()) ;
		
		double xpmax = dxp[dxp.size()*cut]; 
		double xpmin = dxp[dxp.size()*(1.-cut)];
		double xmmax = dxm[dxm.size()*cut]; 
		double xmmin = dxm[dxm.size()*(1.-cut)];
		double ypmax = dyp[dyp.size()*cut]; 
		double ypmin = dyp[dyp.size()*(1.-cut)];
		double ymmax = dym[dym.size()*cut]; 
		double ymmin = dym[dym.size()*(1.-cut)];
		double zpmax = dzp[dzp.size()*cut]; 
		double zpmin = dzp[dzp.size()*(1.-cut)];
		double zmmax = dzm[dzm.size()*cut]; 
		double zmmin = dzm[dzm.size()*(1.-cut)];
		
		double dxpv = 0;
		double dxpc = 0;
		double dxmv = 0;
		double dxmc = 0;
		double dypv = 0;
		double dypc = 0;
		double dymv = 0;
		double dymc = 0;
		double dzpv = 0;
		double dzpc = 0;
		double dzmv = 0;
		double dzmc = 0;
		
		for(size_t i = 0 ; i < dxp.size() ; i++)
		{
			if(dxp[i] > xpmin && dxp[i] < xpmax )
			{
				dxpv += dxp[i] ;
				dxpc++ ;
			}
		}
		for(size_t i = 0 ; i < dxm.size() ; i++)
		{
			if(dxm[i] > xmmin && dxm[i] < xmmax )
			{
				dxmv += dxm[i] ;
				dxmc++ ;
			}
		}

		for(size_t i = 0 ; i < dyp.size() ; i++)
		{
			if(dyp[i] > ypmin && dyp[i] < ypmax )
			{
				dypv += dyp[i] ;
				dypc++ ;
			}
		}
		for(size_t i = 0 ; i < dym.size() ; i++)
		{
			if(dym[i] > ymmin && dym[i] < ymmax )
			{
				dymv += dym[i] ;
				dymc++ ;
			}
		}
		
		for(size_t i = 0 ; i < dzp.size() ; i++)
		{
			if(dzp[i] > zpmin && dzp[i] < zpmax )
			{
				dzpv += dzp[i] ;
				dzpc++ ;
			}
		}
		for(size_t i = 0 ; i < dym.size() ; i++)
		{
			if(dzm[i] > zmmin && dym[i] < zmmax )
			{
				dzmv += dzm[i] ;
				dzmc++ ;
			}
		}
		
		
		ret.push_back((dxpv/dxpc-dxmv/dxmc)/size_x);
		ret.push_back((dypv/dypc-dymv/dymc)/size_y);
		ret.push_back((dzpv/dzpc-dzmv/dzmc)/size_z);
	}
	
	return ret ;
}

std::vector<double> FeatureTree::getMedianMacroscopicStrain(const Geometry * base, double tol) 
{
	Vector disps = getDisplacements() ;
	std::vector<double> ret ;
	if(!disps.size())
		return ret ;
	if(is2D())
	{
		double size_x = dynamic_cast<const Rectangle *>(base)->width() ;
		double size_y = dynamic_cast<const Rectangle *>(base)->height() ;
		std::vector<DelaunayTriangle *> elements = dtree->getElements() ;
		tol = .05*std::min(size_x, size_y) ;
		std::vector<double> dxp ;
		std::vector<double> dyp ;
		std::vector<double> dxm ;
		std::vector<double> dym ;
		std::set<int> doneIds ;
		for(size_t i = 0 ; i < elements.size() ; i++)
		{
			if(elements[i]->getBehaviour() && elements[i]->getBehaviour()->type != VOID_BEHAVIOUR)
			{
				for(size_t j = 0 ; j < elements[i]->getBoundingPoints().size() ; j++)
				{
					
					Point test(elements[i]->getBoundingPoint(j)) ;
					
					base->project(&test);
					if(dist(test, elements[i]->getBoundingPoint(j)) < tol && std::abs(elements[i]->getBoundingPoint(j).getX() - (base->getCenter().getX()+size_x*.5)) < tol && doneIds.find(2*elements[i]->getBoundingPoint(j).getId()) == doneIds.end())
					{
						dxp.push_back(disps[2*elements[i]->getBoundingPoint(j).getId()]) ;
						doneIds.insert(2*elements[i]->getBoundingPoint(j).getId());
					}
					if(dist(test, elements[i]->getBoundingPoint(j)) < tol && std::abs(elements[i]->getBoundingPoint(j).getY() - (base->getCenter().getY()+size_y*.5)) < tol&& doneIds.find(2*elements[i]->getBoundingPoint(j).getId()+1) == doneIds.end())
					{
						dyp.push_back(disps[2*elements[i]->getBoundingPoint(j).getId()+1]) ;
						doneIds.insert(2*elements[i]->getBoundingPoint(j).getId()+1);
					}
					if(dist(test, elements[i]->getBoundingPoint(j)) < tol && std::abs(elements[i]->getBoundingPoint(j).getX() - (base->getCenter().getX()-size_x*.5)) < tol&& doneIds.find(2*elements[i]->getBoundingPoint(j).getId()) == doneIds.end())
					{
						dxm.push_back(disps[2*elements[i]->getBoundingPoint(j).getId()]) ;
						doneIds.insert(2*elements[i]->getBoundingPoint(j).getId());
					}
					if(dist(test, elements[i]->getBoundingPoint(j)) < tol && std::abs(elements[i]->getBoundingPoint(j).getY() - (base->getCenter().getY()-size_y*.5)) < tol&& doneIds.find(2*elements[i]->getBoundingPoint(j).getId()+1) == doneIds.end())
					{
						dym.push_back(disps[2*elements[i]->getBoundingPoint(j).getId()+1]) ;
						doneIds.insert(2*elements[i]->getBoundingPoint(j).getId()+1);
					}
				}
			}
		}
		std::sort(dxp.begin(), dxp.end()) ;
		std::sort(dxm.begin(), dxm.end()) ;
		std::sort(dyp.begin(), dyp.end()) ;
		std::sort(dym.begin(), dym.end()) ;
		
		double dxpv = dxp[dxp.size()/2];
		double dxpc = 1;
		double dxmv = dxm[dxm.size()/2];
		double dxmc = 1;
		double dypv = dyp[dyp.size()/2];
		double dypc = 1;
		double dymv = dym[dym.size()/2] ;
		double dymc = 1;
		
		
		ret.push_back((dxpv/dxpc-dxmv/dxmc)/size_x);
		ret.push_back((dypv/dypc-dymv/dymc)/size_y);
	}
	else
	{
		double size_x = dynamic_cast<const Hexahedron *>(base)->getXSize() ;
		double size_y = dynamic_cast<const Hexahedron *>(base)->getYSize() ;
		double size_z = dynamic_cast<const Hexahedron *>(base)->getZSize() ;
		std::vector<DelaunayTetrahedron *> elements = dtree3D->getElements() ;
		std::vector<double> dxp ;
		std::vector<double> dyp ;
		std::vector<double> dzp ;
		std::vector<double> dxm ;
		std::vector<double> dym ;
		std::vector<double> dzm ;
		std::set<int> doneIds ;
		for(size_t i = 0 ; i < elements.size() ; i++)
		{
			if(elements[i]->getBehaviour()->type != VOID_BEHAVIOUR)
			{
				for(size_t j = 0 ; j < elements[i]->getBoundingPoints().size() ; j++)
				{
					Point test(elements[i]->getBoundingPoint(j)) ;
					base->project(&test);
					if(dist(test, elements[i]->getBoundingPoint(j)) < tol && std::abs(elements[i]->getBoundingPoint(j).getX() - (base->getCenter().getX()-size_x*.5)) < tol && doneIds.find(3*elements[i]->getBoundingPoint(j).getId()) == doneIds.end())
					{
						dxp.push_back(disps[3*elements[i]->getBoundingPoint(j).getId()]) ;
					}
					if(dist(test, elements[i]->getBoundingPoint(j)) < tol && std::abs(elements[i]->getBoundingPoint(j).getY() - (base->getCenter().getY()-size_y*.5)) < tol && doneIds.find(3*elements[i]->getBoundingPoint(j).getId()) == doneIds.end())
					{
						dyp.push_back(disps[3*elements[i]->getBoundingPoint(j).getId()+1]) ;
					}
					if(dist(test, elements[i]->getBoundingPoint(j)) < tol && std::abs(elements[i]->getBoundingPoint(j).getZ() - (base->getCenter().getZ()-size_z*.5)) < tol && doneIds.find(3*elements[i]->getBoundingPoint(j).getId()) == doneIds.end())
					{
						dzp.push_back(disps[3*elements[i]->getBoundingPoint(j).getId()+2]) ;
					}
					if(dist(test, elements[i]->getBoundingPoint(j)) < tol && std::abs(elements[i]->getBoundingPoint(j).getX() - (base->getCenter().getX()-size_x*.5)) < tol && doneIds.find(3*elements[i]->getBoundingPoint(j).getId()) == doneIds.end())
					{
						dxm.push_back(disps[3*elements[i]->getBoundingPoint(j).getId()]) ;
					}
					if(dist(test, elements[i]->getBoundingPoint(j)) < tol && std::abs(elements[i]->getBoundingPoint(j).getY() - (base->getCenter().getY()-size_y*.5)) < tol && doneIds.find(3*elements[i]->getBoundingPoint(j).getId()) == doneIds.end())
					{
						dym.push_back(disps[3*elements[i]->getBoundingPoint(j).getId()+1]) ;
					}
					if(dist(test, elements[i]->getBoundingPoint(j)) < tol && std::abs(elements[i]->getBoundingPoint(j).getY() - (base->getCenter().getZ()-size_z*.5)) < tol && doneIds.find(3*elements[i]->getBoundingPoint(j).getId()) == doneIds.end())
					{
						dzm.push_back(disps[3*elements[i]->getBoundingPoint(j).getId()+2]) ;
					}
				}
			}
		}
		
		std::sort(dxp.begin(), dxp.end()) ;
		std::sort(dxm.begin(), dxm.end()) ;
		std::sort(dyp.begin(), dyp.end()) ;
		std::sort(dym.begin(), dym.end()) ;
		std::sort(dzp.begin(), dzp.end()) ;
		std::sort(dzm.begin(), dzm.end()) ;
		
		
		double dxpv = dxp[dxp.size()/2];
		double dxpc = 1;
		double dxmv = dxm[dxm.size()/2];
		double dxmc = 1;
		double dypv = dyp[dyp.size()/2];
		double dypc = 1;
		double dymv = dym[dym.size()/2] ;
		double dymc = 1;
		double dzpv = dzp[dzp.size()/2];
		double dzpc = 1;
		double dzmv = dzm[dzm.size()/2] ;
		double dzmc = 1;
		
		ret.push_back((dxpv/dxpc-dxmv/dxmc)/size_x);
		ret.push_back((dypv/dypc-dymv/dymc)/size_y);
		ret.push_back((dzpv/dzpc-dzmv/dzmc)/size_z);
	}
	
	return ret ;
}


Vector FeatureTree::getAverageField( FieldType f, const std::vector<DelaunayTriangle *> & tri ) 
{
	Vector avg ;
	Vector buffer ;
	double volume = 0 ;
	avg.resize(fieldTypeElementarySize(f, SPACE_TWO_DIMENSIONAL)) ; buffer.resize(fieldTypeElementarySize(f, SPACE_TWO_DIMENSIONAL)) ; 
	avg = 0 ; buffer = 0 ;
	for(size_t i = 0 ; i < tri.size() ; i++)
	{
		tri[i]->getState().getAverageField( f, buffer ) ;
		avg += buffer * tri[i]->area() ;
		volume += tri[i]->area() ;
	}
	return avg/volume ;
}

Vector FeatureTree::getAverageField( FieldType f, const std::vector<DelaunayTetrahedron *> & tet ) 
{
	Vector avg ;
	Vector buffer ;
	double volume = 0 ;
	avg.resize(fieldTypeElementarySize(f, SPACE_THREE_DIMENSIONAL)) ; buffer.resize(fieldTypeElementarySize(f, SPACE_THREE_DIMENSIONAL)) ; 
	avg = 0 ; buffer = 0 ;
	for(size_t i = 0 ; i < tet.size() ; i++)
	{
		tet[i]->getState().getAverageField( f, buffer ) ;
		avg += buffer * tet[i]->volume() ;
		volume += tet[i]->volume() ;
	}
	return avg/volume ;
}

std::pair<Vector, Vector> FeatureTree::getFieldMinMax( FieldType f, int grid, double t)
{
	Vector min ;
	Vector max ;
	Vector buffer ;
	VirtualMachine vm ;
	if(is2D())
	{
		std::vector<DelaunayTriangle *> elements = getElements2D( grid ) ;
		size_t blocks = elements[0]->getBehaviour()->getNumberOfDegreesOfFreedom()/2 ;
		min.resize(fieldTypeElementarySize(f, SPACE_TWO_DIMENSIONAL, blocks)) ; 
		max.resize(fieldTypeElementarySize(f, SPACE_TWO_DIMENSIONAL, blocks)) ; 
		buffer.resize(fieldTypeElementarySize(f, SPACE_TWO_DIMENSIONAL, blocks)) ; 
		buffer = 0 ;
		size_t start = 0 ;
		while(elements[start]->getBehaviour()->type == VOID_BEHAVIOUR && start < elements.size())
			start++ ;
		elements[start]->getState().getAverageField( f, buffer ) ;
		min = buffer ;
		max = buffer ;
		for(size_t i = start+1 ; i < elements.size() ; i++)
		{
			if(elements[i]->getBehaviour()->type != VOID_BEHAVIOUR)
			{
				elements[i]->getState().getAverageField( f, buffer,&vm, -1, t) ;
				for(size_t j = 0 ; j < min.size() ; j++)
				{
					min[j] = std::min(min[j], buffer[j]) ;
					max[j] = std::max(max[j], buffer[j]) ;
				}
			}
		}
	}
	else
	{
		std::vector<DelaunayTetrahedron *> elements = getElements3D( grid ) ;
		size_t blocks = elements[0]->getBehaviour()->getNumberOfDegreesOfFreedom()/3 ;
		min.resize(fieldTypeElementarySize(f, SPACE_THREE_DIMENSIONAL, blocks)) ;
		max.resize(fieldTypeElementarySize(f, SPACE_THREE_DIMENSIONAL, blocks)) ;
		buffer.resize(fieldTypeElementarySize(f, SPACE_THREE_DIMENSIONAL, blocks)) ; 
		buffer = 0 ;
		size_t start = 0 ;
		while(elements[start]->getBehaviour()->type == VOID_BEHAVIOUR && start < elements.size())
			start++ ;
		elements[start]->getState().getAverageField( f, buffer ) ;
		min = buffer ;
		max = buffer ;
		for(size_t i = start+1 ; i < elements.size() ; i++)
		{
			if(elements[i]->getBehaviour()->type != VOID_BEHAVIOUR)
			{
				elements[i]->getState().getAverageField( f, buffer,&vm, -1, t) ;
				for(size_t j = 0 ; j < min.size() ; j++)
				{
					min[j] = std::min(min[j], buffer[j]) ;
					max[j] = std::max(max[j], buffer[j]) ;
				}
			}
		}
	}
	return std::make_pair(min, max) ;
}


std::pair<Vector, Vector> FeatureTree::getFieldMinMax( FieldType f, const std::vector<DelaunayTriangle *> & tri)
{
	Vector min ;
	Vector max ;
	Vector buffer ;
	min.resize(fieldTypeElementarySize(f, SPACE_TWO_DIMENSIONAL)) ; 
	max.resize(fieldTypeElementarySize(f, SPACE_TWO_DIMENSIONAL)) ; 
	buffer.resize(fieldTypeElementarySize(f, SPACE_TWO_DIMENSIONAL)) ; 
	tri[0]->getState().getAverageField( f, buffer ) ;
	min = buffer ;
	max = buffer ;
	for(size_t i = 1 ; i < tri.size() ; i++)
	{
		tri[i]->getState().getAverageField( f, buffer ) ;
		for(size_t j = 0 ; j < min.size() ; j++)
		{
			min[j] = std::min(min[j], buffer[j]) ;
			max[j] = std::max(max[j], buffer[j]) ;
		}
	}
	return std::make_pair(min, max) ;
}

std::pair<Vector, Vector> FeatureTree::getFieldMinMax( FieldType f, const std::vector<DelaunayTetrahedron *> & tet)
{
	Vector min ;
	Vector max ;
	Vector buffer ;
	min.resize(fieldTypeElementarySize(f, SPACE_THREE_DIMENSIONAL)) ; 
	max.resize(fieldTypeElementarySize(f, SPACE_THREE_DIMENSIONAL)) ; 
	buffer.resize(fieldTypeElementarySize(f, SPACE_THREE_DIMENSIONAL)) ; 
	tet[0]->getState().getAverageField( f, buffer ) ;
	min = buffer ;
	max = buffer ;
	for(size_t i = 1 ; i < tet.size() ; i++)
	{
		tet[i]->getState().getAverageField( f, buffer ) ;
		for(size_t j = 0 ; j < min.size() ; j++)
		{
			min[j] = std::min(min[j], buffer[j]) ;
			max[j] = std::max(max[j], buffer[j]) ;
		}
	}
	return std::make_pair(min, max) ;
}



bool FeatureTree::isStable()
{
	bool needAssemblyinit = needAssembly ;
	bool meshChangeinit = behaviourChange ;
	bool enrichmentChangeinit = enrichmentChange ;
	double crackedVolumeinit = crackedVolume ;
	double damagedVolumeinit = damagedVolume ;
	size_t maxits = maxitPerStep ;
	setMaxIterationsPerStep( 0 ) ;
	elasticStep();
	bool stable = true ;
	
	
	if(is2D())
	{
		std::vector<DelaunayTriangle *> elements ;
		for(auto j = layer2d.begin() ; j!=layer2d.end() ; ++j)
		{
			auto etmp = j->second->getElements() ;
			elements.insert(elements.end(), etmp.begin(), etmp.end()) ;
		}
		
		for(size_t i = 0 ; i < elements.size() ; i++)
		{
			if(elements[i]->getBehaviour() && elements[i]->getBehaviour()->getFractureCriterion() && !elements[i]->getBehaviour()->fractured())
			{
				if(elements[i]->getBehaviour()->getFractureCriterion()->grade(elements[i]->getState()) > 0)
				{
					stable = false ;
					break ;
				}
			}
		}
	}
	else
	{
		std::vector<DelaunayTetrahedron *> elements = dtree3D->getElements() ;
		for(size_t i = 0 ; i < elements.size() ; i++)
		{
			if(elements[i]->getBehaviour() && elements[i]->getBehaviour()->getFractureCriterion())
			{
				if(elements[i]->getBehaviour()->getFractureCriterion()->grade(elements[i]->getState()) > 0)
				{
					stable = false ;
					break ;
				}
			}
		}
	}
	
	
	needAssembly = true ;
	setMaxIterationsPerStep( maxitPerStep ) ;
	behaviourChange = meshChangeinit ;
	enrichmentChange = enrichmentChangeinit ;
	crackedVolume = crackedVolumeinit ;
	damagedVolume = damagedVolumeinit ;

	return stable ;
}

double FeatureTree::getMaximumDisplacement()
{
	state.setStateTo( RENUMBERED, false ) ;

	if( is2D() )
	{
		std::vector<DelaunayTriangle *> tri = dtree->getElements() ;

		double max = 0 ;

		for( std::vector<DelaunayTriangle *>::const_iterator i = tri.begin() ; i != tri.end() ; ++i )
		{
			if( !( *i )->getBehaviour())
			continue ;
			
			if( ( *i )->getBehaviour()->type != VOID_BEHAVIOUR )
				max = std::max( max, ( *i )->getState().getDisplacements().max() ) ;
		}

		return max ;
	}
	else if( is3D() )
	{
		std::vector<DelaunayTetrahedron *> tets = dtree3D->getElements() ;

		double max = 0 ;

		for( std::vector<DelaunayTetrahedron *>::const_iterator i = tets.begin() ; i != tets.end() ; ++i )
		{
			if( !( *i )->getBehaviour())
				continue ;
			
			if( ( *i )->getBehaviour()->type != VOID_BEHAVIOUR )
				max = std::max( max, ( *i )->getState().getDisplacements().max() ) ;
		}

		return max ;
	}

	return 0 ;
}

double FeatureTree::getMinimumDisplacement()
{
	state.setStateTo( RENUMBERED, false ) ;

	if( is2D() )
	{
		std::vector<DelaunayTriangle *> tri = dtree->getElements() ;

		double max = 0 ;

		for( std::vector<DelaunayTriangle *>::const_iterator i = tri.begin() ; i != tri.end() ; ++i )
		{
			if( !( *i )->getBehaviour())
				continue ;
			if( ( *i )->getBehaviour()->type != VOID_BEHAVIOUR )
				max = std::min( max, ( *i )->getState().getDisplacements().min() ) ;
		}

		return max ;
	}
	else if( is3D() )
	{
		std::vector<DelaunayTetrahedron *> tets = dtree3D->getElements() ;

		double max = 0 ;

		for( std::vector<DelaunayTetrahedron *>::const_iterator i = tets.begin() ; i != tets.end() ; ++i )
		{
			if( !( *i )->getBehaviour())
				continue ;
			if( ( *i )->getBehaviour()->type != VOID_BEHAVIOUR )
				max = std::min( max, ( *i )->getState().getDisplacements().min() ) ;
		}

		return max ;
	}

	return 0 ;
}

size_t FeatureTree::numPoints() const
{
	return this->numdofs ;
}

void FeatureTree::print() const
{
	printForFeature( tree[0] );

}

void FeatureTree::printForFeature( const Feature *f ) const
{
	f->print();
	std::vector<Feature *> children = f->getChildren();

	for( size_t i = 0; i != children.size(); ++i )
	{
// 		if ( !(*children)[i]->getChildren().empty())
		printForFeature( children[i] );
	}

}

void FeatureTree::reMesh()
{
	needMeshing = true ;
}

bool FeatureTree::is3D() const
{
	return tree[0]->spaceDimensions() == SPACE_THREE_DIMENSIONAL ;
}

bool FeatureTree::is2D() const
{
	return tree[0]->spaceDimensions() == SPACE_TWO_DIMENSIONAL ;
}

void FeatureTree::initializeElements( bool initialiseFractureCache )
{

	if( !father3D )
		father3D = new TetrahedralElement( elemOrder ) ;

	father3D->compileAndPrecalculate() ;

	if( !father2D )
		father2D = new TriElement( elemOrder ) ;

	father2D->compileAndPrecalculate() ;

	timeval time0, time1 ;
	gettimeofday( &time0, nullptr );


	if( is2D() )
	{

		std::cerr << " initialising..." << std::flush;
		
		for(auto j = layer2d.begin() ; j != layer2d.end() ;j++)
		{
			std::vector<DelaunayTriangle *> tris = j->second->getElements() ;

			int ecounter = 0 ;
			for( size_t i = 0 ; i < tris.size() ; i++ )
			{
				if(!tris[i]->getBehaviour())
				{
					std::cout << "ouch" << std::endl ;
				}
				else
				{
					tris[i]->refresh( father2D );
					tris[i]->getState().initialize( initialiseFractureCache) ;
					#pragma omp critical
					{
						ecounter++ ;
						if(ecounter % 100 == 0)
						std::cerr << "\r initialising... element " << ecounter << "/" << tris.size() << std::flush ;
					}
				}
			}
		}

		gettimeofday( &time1, nullptr );
		int numtris = layer2d.begin()->second->getElements().size() ;
		double delta = time1.tv_sec * 1000000 - time0.tv_sec * 1000000 + time1.tv_usec - time0.tv_usec ;
		std::cerr << "\r initialising... element " << numtris << "/" << numtris << ". Time to initialise (s) " << delta / 1e6 << std::endl ;

		
	}

	if( is3D() )
	{
	  
		std::vector<Mesh <DelaunayTetrahedron, DelaunayTreeItem3D > *> extra3dMeshes ;
		if(layer3d.size() > 1)
		{
		  for(auto i = ++layer3d.begin() ; i != layer3d.end() ; i++)
		    extra3dMeshes.push_back(i->second) ;
		}
		std::vector<DelaunayTetrahedron *> tets = dtree3D->getElements() ;
		std::cout << " initialising..." ;

		#pragma omp parallel for schedule(runtime)
		for( size_t i = 0 ; i < tets.size() ; i++ )
		{
			if(!tets[i]->getBehaviour())
				continue ;
			tets[i]->refresh( father3D );
			tets[i]->getState().initialize( initialiseFractureCache) ;
		}

		gettimeofday( &time1, nullptr );
		double delta = time1.tv_sec * 1000000 - time0.tv_sec * 1000000 + time1.tv_usec - time0.tv_usec ;
		std::cout << "\r initialising... element " << tets.size() << "/" << tets.size() << ". Time to initialise (s) " << delta / 1e6 << std::endl ;

		for(auto j = layer3d.begin() ; j != layer3d.end() ;j++)
		{
			int ecounter = 0 ;
			std::vector<DelaunayTetrahedron *> tetras = j->second->getElements() ;

			for( size_t i = 0 ; i < tetras.size() ; i++ )
			{
				if(!tetras[i]->getBehaviour())
				{
					std::cout << "ouch" << std::endl ;
				}
				tetras[i]->refresh( father3D );
				tetras[i]->getState().initialize( initialiseFractureCache) ;
				#pragma omp critical
				{
					ecounter++ ;
					if(ecounter % 100 == 0)
					std::cerr << "\r initialising... element " << ecounter << "/" << tetras.size() << std::flush ;
				}
			}
		}
	}

}

void FeatureTree::setDeltaTime(double d) 
{
	double prev = deltaTime ;
	deltaTime = d ; realDeltaTime = d ;
	if(dtree)
	{
		for(auto j = layer2d.begin() ; j != layer2d.end() ;j++)
		{
			std::vector<DelaunayTriangle *> triangles = j->second->getElements() ;
			prev = triangles[0]->getBoundingPoint( triangles[0]->getBoundingPoints().size() -1 ).getT() - triangles[0]->getBoundingPoint(0).getT() ;
			double end = triangles[0]->getBoundingPoint( triangles[0]->getBoundingPoints().size() -1 ).getT() ;
			double begin = triangles[0]->getBoundingPoint( 0 ).getT() ;
			if(triangles.size() && triangles[0]->timePlanes() > 1)
			{
				for(size_t i = 0 ; i < triangles.size() ; i++)
				{
					size_t k0 = triangles[i]->getBoundingPoints().size()/triangles[i]->timePlanes() ;
					for(size_t t = 0 ; t < triangles[i]->timePlanes() -1 ; t++)
					{
						for(size_t k = 0 ; k < k0 ; k++)
						{
							triangles[i]->getBoundingPoint(k+k0*t).getT() = end - d + d*t/(triangles[i]->timePlanes()-1) ;
						}
					}
					
					if(triangles[i]->getBehaviour() && triangles[i]->getBehaviour()->type != VOID_BEHAVIOUR)
							triangles[i]->adjustElementaryMatrix( prev, d ) ;
				}
			}
		}
	}
	else
	{
		for(auto j = layer3d.begin() ; j != layer3d.end() ;j++)
		{
			std::vector<DelaunayTetrahedron *> triangles = j->second->getElements() ;
			prev = triangles[0]->getBoundingPoint( triangles[0]->getBoundingPoints().size() -1 ).getT() - triangles[0]->getBoundingPoint(0).getT() ;
			double end = triangles[0]->getBoundingPoint( triangles[0]->getBoundingPoints().size() -1 ).getT() ;
			double begin = triangles[0]->getBoundingPoint( 0 ).getT() ;
			if(triangles.size() && triangles[0]->timePlanes() > 1)
			{
				for(size_t i = 0 ; i < triangles.size() ; i++)
				{
					size_t k0 = triangles[i]->getBoundingPoints().size()/triangles[i]->timePlanes() ;
					for(size_t t = 0 ; t < triangles[i]->timePlanes() -1 ; t++)
					{
						for(size_t k = 0 ; k < k0 ; k++)
						{
							triangles[i]->getBoundingPoint(k+k0*t).getT() = end - d + d*t/(triangles[i]->timePlanes()-1) ;
						}
					}
					
					if(triangles[i]->getBehaviour() && triangles[i]->getBehaviour()->type != VOID_BEHAVIOUR)
							triangles[i]->adjustElementaryMatrix( prev, d ) ;
				}
			}
		}
	}
}


void FeatureTree::moveFirstTimePlanes(double d, std::vector<DelaunayTriangle *> & triangles ) 
{
	double prev = 0.; 
	size_t i = 0 ;
	size_t ndof = 0 ;
	while(ndof == 0 && i < triangles.size())
	{
		ndof = triangles[i]->getBehaviour()->getNumberOfDegreesOfFreedom() ;
		i++ ;
	}
	Vector buff(0.,ndof) ;

	if(dtree)
	{
		prev = triangles[0]->getBoundingPoint( triangles[0]->getBoundingPoints().size() -1 ).getT() - triangles[0]->getBoundingPoint(0).getT() ;
		VirtualMachine vm ;
		if(triangles.size() && triangles[0]->timePlanes() > 1)
		{

			for(size_t i = 0 ; i < triangles.size() ; i++)
			{
				if(triangles[i]->getBehaviour() && triangles[i]->getBehaviour()->type != VOID_BEHAVIOUR)
				{
					size_t k0 = triangles[i]->getBoundingPoints().size()/triangles[i]->timePlanes() ;
					for(size_t t = 0 ; t < triangles[i]->timePlanes() -1 ; t++)
					{
						for(size_t k = 0 ; k < k0 ; k++)
						{
//							std::cout << i << ";" << k << std::endl ;
							Point p(triangles[i]->getBoundingPoint(k+k0*t).getX(),
											triangles[i]->getBoundingPoint(k+k0*t).getY(),
											0.,
							        triangles[i]->getBoundingPoint(k+k0*t).getT() + d*( triangles[i]->timePlanes()-t )/triangles[i]->timePlanes()) ;
							triangles[i]->getStatePointer()->getField( GENERALIZED_VISCOELASTIC_DISPLACEMENT_FIELD, p, buff, false, &vm) ;

							for(size_t n = 0 ; n < ndof ; n++)
							{
								double z = buff[n] ;
								K->setDisplacementByDof( triangles[i]->getBoundingPoint((t+1)*k0+k).getId() * ndof + n, z  );
							}
//							std::cout << i << ";" << k << std::endl ;
						}
					}
				}
			}
		}
		
		if(std::abs(d) > POINT_TOLERANCE_2D)
			  setDeltaTime( prev - d ) ;
		
	}
}

void FeatureTree::moveFirstTimePlanes(double d, std::vector<DelaunayTetrahedron *> & tets ) 
{
	double prev = 0.; 
	size_t i = 0 ;
	size_t ndof = 0 ;
	while(ndof == 0 && i < tets.size())
	{
		ndof = tets[i]->getBehaviour()->getNumberOfDegreesOfFreedom() ;
		i++ ;
	}
	Vector buff(0.,ndof) ;

	if(dtree)
	{
		prev = tets[0]->getBoundingPoint( tets[0]->getBoundingPoints().size() -1 ).getT() - tets[0]->getBoundingPoint(0).getT() ;
		VirtualMachine vm ;
		if(tets.size() && tets[0]->timePlanes() > 1)
		{

			for(size_t i = 0 ; i < tets.size() ; i++)
			{
				if(tets[i]->getBehaviour() && tets[i]->getBehaviour()->type != VOID_BEHAVIOUR)
				{
					size_t k0 = tets[i]->getBoundingPoints().size()/tets[i]->timePlanes() ;
					for(size_t t = 0 ; t < tets[i]->timePlanes() -1 ; t++)
					{
						for(size_t k = 0 ; k < k0 ; k++)
						{
//							std::cout << i << ";" << k << std::endl ;
							Point p(tets[i]->getBoundingPoint(k+k0*t).getX(),
											tets[i]->getBoundingPoint(k+k0*t).getY(),
											tets[i]->getBoundingPoint(k+k0*t).getZ(),
							        tets[i]->getBoundingPoint(k+k0*t).getT() + d*( tets[i]->timePlanes()-t )/tets[i]->timePlanes()) ;
							tets[i]->getStatePointer()->getField( GENERALIZED_VISCOELASTIC_DISPLACEMENT_FIELD, p, buff, false, &vm) ;

							for(size_t n = 0 ; n < ndof ; n++)
							{
								double z = buff[n] ;
								K->setDisplacementByDof( tets[i]->getBoundingPoint((t+1)*k0+k).getId() * ndof + n, z  );
							}
//							std::cout << i << ";" << k << std::endl ;
						}
					}
				}
			}
		}
		
		if(std::abs(d) > POINT_TOLERANCE_2D)
			  setDeltaTime( prev - d ) ;
		
	}
}

void FeatureTree::generateElements()
{
	for( size_t i = 0 ; i < boundaryCondition.size() ; i++ )
		boundaryCondition[i]->clearCache() ;

	if( dtree || dtree3D )
	{
		if( K )
			K->clear() ;

	}

	needMeshing = false ;

	double pointDensity = 0 ;

	if( is2D() )
		pointDensity = sqrt( tree[0]->area() / ( tree[0]->getBoundingPoints().size() + tree[0]->getInPoints().size() ) ) ;
	else
		pointDensity = pow( tree[0]->volume() / ( tree[0]->getBoundingPoints().size() + tree[0]->getInPoints().size() ), .33333333333 ) ;

	std::cout << "space meshed with " << pointDensity << " points per unit length" << std::endl ;

	std::valarray<Point> bbox( 8 ) ;
	double min_x = 0, min_y = 0, max_x = 0, max_y = 0, max_z = 0, min_z = 0;

	for( size_t j  =  0 ; j <  this->tree[0]->getBoundingPoints().size() ; j++ )
	{
		if( this->tree[0]->getBoundingPoint( j ).getY() < min_y )
			min_y = this->tree[0]->getBoundingPoint( j ).getY() ;

		if( this->tree[0]->getBoundingPoint( j ).getY() > max_y )
			max_y = this->tree[0]->getBoundingPoint( j ).getY() ;

		if( this->tree[0]->getBoundingPoint( j ).getX() < min_x )
			min_x = this->tree[0]->getBoundingPoint( j ).getX() ;

		if( this->tree[0]->getBoundingPoint( j ).getX() > max_x )
			max_x = this->tree[0]->getBoundingPoint( j ).getX() ;

		if( this->tree[0]->getBoundingPoint( j ).getZ() < min_z )
			min_z = this->tree[0]->getBoundingPoint( j ).getZ() ;

		if( this->tree[0]->getBoundingPoint( j ).getZ() > max_z )
			max_z = this->tree[0]->getBoundingPoint( j ).getZ() ;
	}

	bbox[0] = Point( min_x, min_y, min_z ) ;

	bbox[1] = Point( min_x, min_y, max_z ) ;

	bbox[2] = Point( min_x, max_y, min_z ) ;

	bbox[3] = Point( min_x, max_y, max_z ) ;

	bbox[4] = Point( max_x, min_y, min_z ) ;

	bbox[5] = Point( max_x, min_y, max_z ) ;

	bbox[6] = Point( max_x, max_y, min_z ) ;

	bbox[7] = Point( max_x, max_y, max_z ) ;

	int bpcount = 0 ;
	size_t basepoints = 0 ;
	std::cerr << " getting mesh points..." << std::flush ;
	std::vector<Feature *> nullFatherFeatures ;

	for( size_t i  = 1 ; i < tree.size() ; i++ )
	{
		if( !tree[i]->isEnrichmentFeature && !tree[i]->isVirtualFeature )
		{
			if( !tree[i]->getFather() )
				nullFatherFeatures.push_back( tree[i] );
		}
	}


	for( size_t i  = 0 ; i < tree.size() ; i++ )
	{
		std::cerr << "\r getting mesh points... feature " << i << "/" << tree.size() << std::flush ;

		if( !tree[i]->isEnrichmentFeature && !tree[i]->isVirtualFeature )
		{
			std::vector<Feature *> descendants = tree[i]->getDescendants() ;
			std::stable_sort( descendants.begin(), descendants.end() ) ;

			for( size_t j  =  0 ; j <  tree[i]->getBoundingPoints().size() ; j++ )
			{
				bool isIn = false ;

				std::vector<const Geometry *> potentialFeaturestmp  ;
				std::vector<Feature *> potentialFeatures ;

				if( is2D() )
					potentialFeaturestmp = grid->coOccur( tree[i]->getBoundingPoint( j ) ) ;
				else
					potentialFeaturestmp = grid3d->coOccur( tree[i]->getBoundingPoint( j ) ) ;

				for( size_t l = 0 ; l < potentialFeaturestmp.size() ; l++ )
					potentialFeatures.push_back( const_cast<Feature *>(dynamic_cast<const Feature *>( potentialFeaturestmp[l] )) ) ;

				std::vector<Feature *> potentialChildren ;

				for( size_t l = 0 ; l < potentialFeatures.size() ; l++ )
				{
					if( !potentialFeatures[l]->isEnrichmentFeature
					        && std::binary_search( descendants.begin(), descendants.end(), potentialFeatures[l] ) )
					{
						potentialChildren.push_back( potentialFeatures[l] ) ;
					}
				}
				
				if(tree.size() < 128)
					potentialChildren = descendants ; 

				for( size_t k  =  0 ; k <  potentialChildren.size() ; k++ )
				{
					if( ( !potentialChildren[k]->isVirtualFeature
					        && potentialChildren[k]->inBoundary( tree[i]->getBoundingPoint( j ), pointDensity * .25 ) )
					        || ( potentialChildren[k]->isVirtualFeature
					             && tree[i]->isVirtualFeature
					             && ( dynamic_cast<VirtualFeature *>( potentialChildren[k] )->getSource()
					                  != dynamic_cast<VirtualFeature *>( tree[i] )->getSource() )
					             && potentialChildren[k]->inBoundary( tree[i]->getBoundingPoint( j ), pointDensity * .25 )
					           )
					  )
					{
						if( potentialChildren[k]->getBoundingPoints().size() )
						{
							isIn = true ;
							break ;
						}
					}
				}

				if( i && tree[i]->getFather() && !inRoot( tree[i]->getBoundingPoint( j ) ) )
					isIn = true ;

				if( !isIn && tree[i]->isVirtualFeature && !tree[i]->in( tree[i]->getBoundingPoint( j ) ) )
					isIn = true ;

				if( tree[i]->getFather() && tree[i]->getFather()->onBoundary( tree[i]->getBoundingPoint( j ), pointDensity * .25 ) )
					isIn = true ;

				if( tree[i]->getFather() && !isIn && i && tree[0]->onBoundary( tree[i]->getBoundingPoint( j ), pointDensity * .25 ) )
					isIn = true ;

				if( !tree[i]->getFather() && i)
				{
					isIn = false ;
					for( size_t k = 0 ; k < nullFatherFeatures.size() ; k++ )
					{
						if( nullFatherFeatures[k] == tree[i] )
							break ;

						if( nullFatherFeatures[k]->inBoundary(tree[i]->getBoundingPoint( j ), 2.*POINT_TOLERANCE_2D) )
						{
							isIn = true ;
							break ;
						}
					}
					
					Point proj( tree[i]->getBoundingPoint( j ) ) ;
					tree[0]->project( &proj ) ;
					if( dist( proj, tree[i]->getBoundingPoint( j ) ) < 2.*POINT_TOLERANCE_2D )
						isIn = false ;
				}
				
				//the border is always defined by non-root features of nullptr father
				if(!i && !isIn)
				{
					for( size_t k = 0 ; k < nullFatherFeatures.size() ; k++ )
					{
						Point proj( tree[i]->getBoundingPoint( j ) ) ;
						nullFatherFeatures[k]->project( &proj ) ;
						if( dist( proj, tree[i]->getBoundingPoint( j ) ) < 2.*POINT_TOLERANCE_2D )
						{
							isIn = true ;
							break ;
						}
					}
				}
				
				
				if( !isIn )
				{
					meshPoints.push_back( std::pair<Point *, Feature *>( &tree[i]->getBoundingPoint( j ), this->tree[i] ) ) ;

					if( i == 0 )
						basepoints++ ;
				}
			}

			for( size_t j  =  0 ; j <  tree[i]->getInPoints().size() ; j++ )
			{
				bool isIn = false ;
				std::vector<const Geometry *> potentialFeaturestmp  ;

				if( is2D() )
					potentialFeaturestmp = grid->coOccur( tree[i]->getInPoint( j ) ) ;
				else
					potentialFeaturestmp = grid3d->coOccur( tree[i]->getInPoint( j ) ) ;

				std::vector<Feature *> potentialFeatures ;

				for( size_t k = 0 ; k < potentialFeaturestmp.size() ; ++k )
					potentialFeatures.push_back( const_cast<Feature *>(dynamic_cast<const Feature *>( potentialFeaturestmp[k] )) ) ;

				std::vector<Feature *> potentialChildren ;

				for( size_t l = 0 ; l < potentialFeatures.size() ; l++ )
				{
					if( !potentialFeatures[l]->isVirtualFeature
					        && !potentialFeatures[l]->isEnrichmentFeature
					        && std::binary_search( descendants.begin(), descendants.end(), potentialFeatures[l] ) )
						potentialChildren.push_back( potentialFeatures[l] ) ;
				}
				
				if(tree.size() < 128)
					potentialChildren = descendants ; 

				for( size_t k  =  0 ; k <  potentialChildren.size() ; k++ )
				{
					if(
					    (
					        !potentialChildren[k]->isVirtualFeature
					        && potentialChildren[k]->inBoundary( tree[i]->getInPoint( j ), pointDensity )
					    )
					    ||
					    (
					        potentialChildren[k]->isVirtualFeature
					        && tree[i]->isVirtualFeature
					        &&(
					            dynamic_cast<VirtualFeature *>( potentialChildren[k] )->getSource()
					            != dynamic_cast<VirtualFeature *>( tree[i] )->getSource()
					        )
					        && potentialChildren[k]->inBoundary( tree[i]->getInPoint( j ), pointDensity )
					    )
					)
					{
						if( potentialChildren[k]->getBoundingPoints().size() )
						{
							isIn = true ;
							break ;
						}
					}
				}

				if( i && !inRoot( tree[i]->getInPoint( j ) ) )
					isIn = true ;

				if( tree[i]->getFather() && tree[i]->getFather()->onBoundary( tree[i]->getInPoint( j ), pointDensity ) )
					isIn = true ;

				if( tree[i]->isVirtualFeature && !tree[i]->in( tree[i]->getInPoint( j ) ) )
					isIn = true ;

				if( tree[i]->getFather() && i && tree[0]->onBoundary( tree[i]->getInPoint( j ), pointDensity ) )
					isIn = true ;

				if( !tree[i]->getFather()&& i )
					isIn = false ;


				for( size_t k = 0 ; k < nullFatherFeatures.size() ; k++ )
				{
					if( tree[i] == nullFatherFeatures[k] )
						break ;

					if( nullFatherFeatures[k]->inBoundary(tree[i]->getInPoint( j ), 2.*POINT_TOLERANCE_2D) )
					{
						isIn = true ;
						break ;
					}
				}

				if( !isIn )
				{
					meshPoints.push_back( std::pair<Point *, Feature *>( &tree[i]->getInPoint( j ), tree[i] ) ) ;

					if( i == 0 )
						basepoints++ ;
				}
			}

		}
	}

	std::cerr << "...done" << std::endl ;

	size_t count  = 0 ;

	if( computeIntersections )
	{
		for( size_t i = 1 ;  i < tree.size() ; i++ )
		{
			if( !tree[i]->isEnrichmentFeature && !tree[i]->isVirtualFeature && tree[i]->getFather() != nullptr )
			{
				std::vector<const Geometry *> coOccuringFeaturestmp ;
				std::vector<Feature *> descendants = tree[i]->getDescendants() ;

				if( is3D() )
					coOccuringFeaturestmp = grid3d->coOccur( tree[i] ) ;
				else
					coOccuringFeaturestmp = grid->coOccur( tree[i] ) ;

				std::vector<const Feature *> coOccuringFeatures ;

				for( size_t k = 0 ; k < coOccuringFeaturestmp.size() ; k++ )
					coOccuringFeatures.push_back( dynamic_cast<const Feature *>( coOccuringFeaturestmp[k] ) ) ;

				for( size_t j  = 0 ; j < coOccuringFeatures.size() ; j++ )
				{
					if( !coOccuringFeatures[j]->isEnrichmentFeature
					        && !coOccuringFeatures[j]->isVirtualFeature
					        && tree[i] != coOccuringFeatures[j]
					        && tree[i]->intersects( coOccuringFeatures[j] ) )
					{
						std::vector<Point> inter = tree[i]->intersection( coOccuringFeatures[j] ) ;

						for( size_t k = 0 ;  k < inter.size() ; k++ )
						{

							bool indescendants = false ;

							for( size_t l = 0 ; l < descendants.size() ; l++ )
							{
								if( !descendants[l]->isVirtualFeature && descendants[l]->inBoundary( inter[k], pointDensity ) )
								{
									indescendants = true ;
									break ;
								}
							}

//
							if( !indescendants )
							{
								if( inRoot( inter[k] ) )
								{
									Point *p = new Point( inter[k] ) ;
									additionalPoints.push_back( p ) ;
									++count ;
									meshPoints.push_back( std::make_pair( p, tree[i] ) ) ;
								}
							}

							if( count % 100 == 0 )
								std::cerr << "\r adding intersection points... " << count << std::flush ;
						}
					}
				}
			}
		}

		for( size_t i = 1 ;  i < tree.size() ; i++ )
		{
			if( !tree[i]->isEnrichmentFeature && tree[i]->getBoundingPoints().size() && !tree[i]->isVirtualFeature && tree[0]->intersects( tree[i] ) && tree[i]->getFather() != nullptr )
			{
				std::vector<Point> inter = tree[0]->intersection( tree[i] ) ;
				std::vector<Feature *> descendants = tree[i]->getDescendants() ;
				std::vector<Feature *> fatherdescendants = tree[0]->getDescendants() ;

				for( size_t k = 0 ;  k < inter.size() ; k++ )
				{
					bool indescendants = false ;

					for( size_t l = 0 ; l < descendants.size() ; l++ )
					{
						if( descendants[l]->inBoundary( inter[k], pointDensity ) )
						{
							indescendants = true ;
							break ;
						}
					}

					for( size_t l = 0 ; l < fatherdescendants.size() ; l++ )
					{
						if( fatherdescendants[l] != tree[i] && !fatherdescendants[l]->isVirtualFeature && fatherdescendants[l]->inBoundary( inter[k], pointDensity ) && fatherdescendants[l]->getBoundingPoints().size() )
						{
							indescendants = true ;
							break ;
						}
					}

					if( is3D() )
					{

						Point proj( inter[k] ) ;
						tree[0]->project( &proj ) ;
						Point proj0( inter[k] + Point( 2.*POINT_TOLERANCE_3D, 0, 0 ) ) ;
						Point proj1( inter[k] + Point( -2.*POINT_TOLERANCE_3D, 0, 0 ) ) ;
						Point proj2( inter[k] + Point( 0, 2.*POINT_TOLERANCE_3D, 0 ) ) ;
						Point proj3( inter[k] + Point( 0, -2.*POINT_TOLERANCE_3D, 0 ) ) ;
						Point proj4( inter[k] + Point( 0, 0, 2.*POINT_TOLERANCE_3D ) ) ;
						Point proj5( inter[k] + Point( 0, 0, -2.*POINT_TOLERANCE_3D ) ) ;

						int position = tree[0]->in( proj0 )
						               + tree[0]->in( proj1 )
						               + tree[0]->in( proj2 )
						               + tree[0]->in( proj3 )
						               + tree[0]->in( proj4 )
						               + tree[0]->in( proj5 ) ;

						bool onSurface = ( position == 5 ) ;
						bool onEdge = ( position == 4 ) ;
						bool onVertex = ( position == 3 ) ;
						proj0 = ( inter[k] + Point( pointDensity, 0, 0 ) ) ;
						proj1 = ( inter[k] + Point( -pointDensity, 0, 0 ) ) ;
						proj2 = ( inter[k] + Point( 0, pointDensity, 0 ) ) ;
						proj3 = ( inter[k] + Point( 0, -pointDensity, 0 ) ) ;
						proj4 = ( inter[k] + Point( 0, 0, pointDensity ) ) ;
						proj5 = ( inter[k] + Point( 0, 0, -pointDensity ) ) ;
						int tooClose =  tree[0]->in( proj0 )
						                + tree[0]->in( proj1 )
						                + tree[0]->in( proj2 )
						                + tree[0]->in( proj3 )
						                + tree[0]->in( proj4 )
						                + tree[0]->in( proj5 ) ;

						// no overlap with other features, intersection is indeed on the surface, and not too near another part of the surface
						if( !indescendants && squareDist3D( proj, inter[k] ) < POINT_TOLERANCE_3D * POINT_TOLERANCE_3D && /*inRoot(inter[k]) && */( ( onSurface && tooClose == 5 ) || ( onEdge && tooClose == 4 ) || onVertex ) )
						{
							Point *p = new Point( inter[k] ) ;
							additionalPoints.push_back( p ) ;
							++count ;
							meshPoints.push_back( std::make_pair( p, tree[i] ) ) ;
						}
					}
					else
					{
						Point proj( inter[k] ) ;
						tree[0]->project( &proj ) ;
						Point proj0( inter[k] + Point( 2.*POINT_TOLERANCE_3D, 0, 0 ) ) ;
						Point proj1( inter[k] + Point( -2.*POINT_TOLERANCE_3D, 0, 0 ) ) ;
						Point proj2( inter[k] + Point( 0, 2.*POINT_TOLERANCE_3D, 0 ) ) ;
						Point proj3( inter[k] + Point( 0, -2.*POINT_TOLERANCE_3D, 0 ) ) ;


						int position = tree[0]->in( proj0 )
						               + tree[0]->in( proj1 )
						               + tree[0]->in( proj2 )
						               + tree[0]->in( proj3 );

						bool onEdge = ( position == 3 ) ;
						bool onVertex = ( position == 2 ) ;
						proj0 = ( inter[k] + Point( pointDensity, 0, 0 ) ) ;
						proj1 = ( inter[k] + Point( -pointDensity, 0, 0 ) ) ;
						proj2 = ( inter[k] + Point( 0, pointDensity, 0 ) ) ;
						proj3 = ( inter[k] + Point( 0, -pointDensity, 0 ) ) ;

						int tooClose =  tree[0]->in( proj0 )
						                + tree[0]->in( proj1 )
						                + tree[0]->in( proj2 )
						                + tree[0]->in( proj3 );

						// no overlap with other features, intersection is indeed on the surface, and not too near another part of the surface
						if( !indescendants && squareDist3D( proj, inter[k] ) < POINT_TOLERANCE_3D * POINT_TOLERANCE_3D && inRoot( inter[k] ) && ( ( onEdge && tooClose == 3 ) || onVertex ) )
						{
							Point *p = new Point( inter[k] ) ;
							additionalPoints.push_back( p ) ;
							++count ;
							meshPoints.push_back( std::make_pair( p, tree[i] ) ) ;
						}
					}

				}

				if( count % 100 == 0 )
					std::cerr << "\r adding intersection points... " << count << std::flush ;
			}

		}
	}

	std::cerr << "\r adding intersection points... " << count << " ...done." << std::endl ;
	count = 0 ;

	//shuffle for efficiency
	shuffleMeshPoints() ;


	if( is2D() )
	{
		additionalPoints.push_back( new Point( bbox[0] ) ) ;
		additionalPoints.push_back( new Point( bbox[2] ) ) ;
		additionalPoints.push_back( new Point( bbox[4] ) ) ;
		additionalPoints.push_back( new Point( bbox[6] ) ) ;
		for(auto i = extraPoints.begin() ; i!= extraPoints.end() ; ++i)
		{
			meshPoints.push_front( std::make_pair( *i, tree[0] ) ) ;
		}
		meshPoints.push_front( std::make_pair( additionalPoints[additionalPoints.size() - 4], tree[0] ) ) ;
		meshPoints.push_front( std::make_pair( additionalPoints[additionalPoints.size() - 3], tree[0] ) ) ;
		meshPoints.push_front( std::make_pair( additionalPoints[additionalPoints.size() - 2], tree[0] ) ) ;
		meshPoints.push_front( std::make_pair( additionalPoints[additionalPoints.size() - 1], tree[0] ) ) ;


		Mesh<DelaunayTriangle, DelaunayTreeItem> * oldDtree = dtree ;

		dtree = new DelaunayTree( meshPoints[0].first, meshPoints[1].first, meshPoints[2].first ) ;
		dtree->insert( meshPoints[3].first ) ;
		layer2d[-1] = dtree ;

		for( size_t i  = 0 ; i < tree.size() ; i++ )
		{
			if( !tree[i]->isEnrichmentFeature && !tree[i]->isVirtualFeature && tree[i]->getLayer() != -1)
			{
				if(layer2d.find(tree[i]->getLayer()) == layer2d.end())
				{
					layer2d[tree[i]->getLayer()] = new DelaunayTree( meshPoints[0].first, meshPoints[1].first, meshPoints[2].first) ;
					layer2d[tree[i]->getLayer()]->insert( meshPoints[3].first ) ;
				}
			}
		}
		
		std::vector<size_t> iterators(meshPoints.size()-4) ;
		for(size_t i = 0 ;i < iterators.size() ; i++)
		  iterators[i] = i+4 ;

		std::random_shuffle(iterators.begin(), iterators.end());
		
		for( size_t i = 0 ; i < iterators.size() ; i++ )
		{
			if( ( i ) % 1000 == 0 )
				std::cerr << "\r generating triangles... point " << count << "/" << meshPoints.size() << std::flush ;

			++count ;

			if( *meshPoints[iterators[i]].first != bbox[0] &&
			        *meshPoints[iterators[i]].first != bbox[2] &&
			        *meshPoints[iterators[i]].first != bbox[4] &&
			        *meshPoints[iterators[i]].first != bbox[6] && ( inRoot( *meshPoints[iterators[i]].first ) || meshPoints[iterators[i]].second->getFather() == nullptr )
			  )
			{
				for(auto j = layer2d.begin() ; j != layer2d.end() ; j++)
				{
					j->second->insert(meshPoints[iterators[i]].first) ;
				}
			}
		}

		std::cerr << "\r generating triangles.... point " << meshPoints.size() - 3 << "/" << meshPoints.size() - 4 << " ...done" << std::endl ;
		

		bool correct = false ;
		int tries = correctionSteps ;

		while( !correct && tries )
		{
			std::vector< DelaunayTriangle * > tets = dtree->getElements();
			std::vector< Point *> to_insert ;

			for( size_t i = 0 ; i < tets.size() ; i++ )
			{
				Point *test = checkElement( tets[i] );

				if( test )
				{
					to_insert.push_back( test );
				}
			}

			if( to_insert.empty() )
				correct = true ;

			for( size_t i = 0 ; i < to_insert.size() ; i++ )
			{
				std::cerr << "\r generating triangles.... point " << ++count << "/" << meshPoints.size() - 3 << std::flush ;

				if( *to_insert[i] != bbox[0] &&
				        *to_insert[i] != bbox[2] &&
				        *to_insert[i] != bbox[4] &&
				        *to_insert[i] != bbox[6] &&
				        inRoot( *to_insert[i] )
				  )
				{
					for(auto j = layer2d.begin() ; j != layer2d.end() ; j++)
					{
						j->second->insert(to_insert[i] ) ;
					}
				}

				if( to_insert[i]->getId() == -1 )
					delete to_insert[i] ;
			}

			tries-- ;

		}

		std::cerr << " ...done." << std::endl ;

		for(size_t i = 0 ; i < refinementZones.size() ; i++)
			quadTreeRefine(refinementZones[i]) ;
		
		if( oldDtree )
		{
			setElementBehavioursFromMesh<Mesh<DelaunayTriangle, DelaunayTreeItem>, DelaunayTriangle>( oldDtree, dtree ) ;
			delete oldDtree ;
		}

	}
	else if( is3D() )
	{
		additionalPoints.push_back( new Point( bbox[0] ) ) ;
		additionalPoints.push_back( new Point( bbox[1] ) ) ;
		additionalPoints.push_back( new Point( bbox[7] ) ) ;
		additionalPoints.push_back( new Point( bbox[2] ) ) ;
		additionalPoints.push_back( new Point( bbox[3] ) ) ;
		additionalPoints.push_back( new Point( bbox[4] ) ) ;
		additionalPoints.push_back( new Point( bbox[5] ) ) ;
		additionalPoints.push_back( new Point( bbox[6] ) ) ;
		for(auto i = extraPoints.begin() ; i!= extraPoints.end() ; ++i)
		{
			meshPoints.push_front( std::make_pair( *i, tree[0] ) ) ;
		}
		meshPoints.push_front( std::make_pair( additionalPoints[additionalPoints.size() - 8], tree[0] ) ) ;
		meshPoints.push_front( std::make_pair( additionalPoints[additionalPoints.size() - 7], tree[0] ) ) ;
		meshPoints.push_front( std::make_pair( additionalPoints[additionalPoints.size() - 6], tree[0] ) ) ;
		meshPoints.push_front( std::make_pair( additionalPoints[additionalPoints.size() - 5], tree[0] ) ) ;
		meshPoints.push_front( std::make_pair( additionalPoints[additionalPoints.size() - 4], tree[0] ) ) ;
		meshPoints.push_front( std::make_pair( additionalPoints[additionalPoints.size() - 3], tree[0] ) ) ;
		meshPoints.push_front( std::make_pair( additionalPoints[additionalPoints.size() - 2], tree[0] ) ) ;
		meshPoints.push_front( std::make_pair( additionalPoints[additionalPoints.size() - 1], tree[0] ) ) ;


		Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * oldDtree = dtree3D ;
//         std::vector<Geometry *> domains ;
//         domains.push_back(new Hexahedron(200, 200, 200, 100, 100, 100));
//         domains.push_back(new Hexahedron(200, 200, 200, 100, 100, 300));
//         domains.push_back(new Hexahedron(200, 200, 200, 100, 300, 100));
//         domains.push_back(new Hexahedron(200, 200, 200, 300, 100, 100));
//         domains.push_back(new Hexahedron(200, 200, 200, 300, 300, 100));
//         domains.push_back(new Hexahedron(200, 200, 200, 100, 300, 300));
//         domains.push_back(new Hexahedron(200, 200, 200, 300, 100, 300));
//         domains.push_back(new Hexahedron(200, 200, 200, 300, 300, 300));
		
		dtree3D = new /*Parallel*/DelaunayTree3D( meshPoints[0].first, meshPoints[1].first, meshPoints[2].first, meshPoints[3].first/*, domains*/) ;
		dtree3D->insert( meshPoints[4].first ) ;
		dtree3D->insert( meshPoints[5].first ) ;
		dtree3D->insert( meshPoints[6].first ) ;
		dtree3D->insert( meshPoints[7].first ) ;
		layer3d[-1] = dtree3D ;
		for( size_t i  = 0 ; i < tree.size() ; i++ )
		{
			if( !tree[i]->isEnrichmentFeature && !tree[i]->isVirtualFeature && tree[i]->getLayer() != -1)
			{
				if(layer3d.find(tree[i]->getLayer()) == layer3d.end())
				{
					layer3d[tree[i]->getLayer()] = new /*Parallel*/DelaunayTree3D( meshPoints[0].first, meshPoints[1].first, meshPoints[2].first, meshPoints[3].first/*, domains*/) ;
					layer3d[tree[i]->getLayer()]->insert( meshPoints[4].first ) ;
					layer3d[tree[i]->getLayer()]->insert( meshPoints[5].first ) ;
					layer3d[tree[i]->getLayer()]->insert( meshPoints[6].first ) ;
					layer3d[tree[i]->getLayer()]->insert( meshPoints[7].first ) ;
				}
			}
		}

		std::pair<std::vector<int>, std::vector<int> > pb ;

		for( size_t i = 0 ; i < meshPoints.size() ; i++ )
		{
			if( meshPoints[i].first == nullptr )
				pb.first.push_back( i ) ;

			if( meshPoints[i].second == nullptr )
				pb.second.push_back( i ) ;
		}

		std::vector<Point *> toInsert ;

		for( auto i = meshPoints.begin() + 8 ; i != this->meshPoints.end(); ++i )
		{
			if( ( i - meshPoints.begin() ) % 1000 == 0 )
				std::cerr << "\r generating tetrahedrons... point " << ++count * 1000 << "/" << meshPoints.size() - 8 <<  std::flush ;

			if( *i->first != bbox[0] &&
			        *i->first != bbox[1] &&
			        *i->first != bbox[2] &&
			        *i->first != bbox[3] &&
			        *i->first != bbox[4] &&
			        *i->first != bbox[5] &&
			        *i->first != bbox[6] &&
			        *i->first != bbox[7]
			  )
			{
				for(auto j = layer3d.begin() ; j != layer3d.end() ; j++)
				{
					j->second->insert(i->first) ;
				}

				if( i->first->getId() == -1 )
				{
					std::cout << "insertion failed" << std::endl ;
					toInsert.push_back( i->first ) ;
				}
			}
		}

		bool correct = false ;
		int tries = correctionSteps ;

		while( !correct && tries )
		{
			std::vector< DelaunayTetrahedron * > tets = dtree3D->getElements();
			std::vector< Point *> to_insert ;

			for( size_t i = 0 ; i < tets.size() ; i++ )
			{
				Point *test = checkElement( tets[i] );

				if( test )
				{
					to_insert.push_back( test );
				}
			}

			if( to_insert.empty() )
				correct = true ;

			for( size_t i = 0 ; i < to_insert.size() ; i++ )
			{
				std::cerr << "\r generating tetrahedrons.... point " << ++count << "/" << meshPoints.size() - 3 << std::flush ;

				if( *to_insert[i] != bbox[0] &&
				        *to_insert[i] != bbox[1] &&
				        *to_insert[i] != bbox[2] &&
				        *to_insert[i] != bbox[3] &&
				        *to_insert[i] != bbox[4] &&
				        *to_insert[i] != bbox[5] &&
				        *to_insert[i] != bbox[6] &&
				        *to_insert[i] != bbox[7] &&
				        inRoot( *to_insert[i] )
				  )
					dtree3D->insert( to_insert[i] ) ;

				if( to_insert[i]->getId() == -1 )
					delete to_insert[i] ;
			}

			tries-- ;

		}

		std::cerr << " ...done." << std::endl ;
//		dtree3D->purge() ;

		if( oldDtree )
		{
			setElementBehavioursFromMesh<Mesh<DelaunayTetrahedron, DelaunayTreeItem3D>, DelaunayTetrahedron>( oldDtree, dtree3D ) ;
			delete oldDtree ;
		}

	}
}

void FeatureTree::shuffleMeshPoints()
{
	std::random_shuffle( meshPoints.begin(), meshPoints.end() ) ;
	return ;
	std::cout << "shuffling mesh points... " ;

	std::deque<std::pair<Point *, const Feature * > > shuffled ;

	for( size_t i = 0 ; i < meshPoints.size() ; i++ )
		shuffled.push_back( meshPoints[i] ) ;

	meshPoints.clear() ;

	std::random_shuffle( shuffled.begin(), shuffled.end() ) ;

	std::vector<bool> visited ;

	for( size_t i = 0 ; i < shuffled.size() ; i++ )
		visited.push_back( true ) ;

	size_t ix = 0 ;
	size_t iy = 0 ;
	size_t iz = 0 ;

	size_t p = 0 ;

	size_t np = shuffled.size() / 2 ;

	if( is2D() )
	{
		np =  std::pow( np, 0.5 ) + 1 ;
		Grid *shufflingGrid = new Grid( static_cast<Sample *>( tree[0] )->width() * 1.01, static_cast<Sample *>( tree[0] )->height() * 1.01, np, tree[0]->getCenter() ) ;

		while( meshPoints.size() < shuffled.size() )
		{
			Point ptest( shuffled[p].first->getX(), shuffled[p].first->getY() ) ;

			if( ( visited[p] ) && ( shufflingGrid->pixels[ix][iy]->coOccur( ptest ) ) )
			{
				visited[p] = false ;
				meshPoints.push_back( shuffled[p] ) ;
				ix++ ;

				if( ix == shufflingGrid->getLengthX() )
				{
					ix = 0 ;
					iy++ ;

					if( iy == shufflingGrid->getLengthY() )
						iy = 0 ;
				}

				p = 0 ;
			}
			else
			{
				p++ ;

				if( p == shuffled.size() )
				{
					p = 0 ;
					ix++ ;

					if( ix == shufflingGrid->getLengthX() )
					{
						ix = 0 ;
						iy++ ;

						if( iy == shufflingGrid->getLengthY() )
							iy = 0 ;
					}

//					shufflingGrid->pixels[ix][iy]->print() ;
				}
			}
		}

		delete shufflingGrid ;
	}

	if( is3D() )
	{
		np = std::pow( np, 0.3333333 ) + 1 ;
		Grid3D *shufflingGrid = new Grid3D( static_cast<Sample3D *>( tree[0] )->getXSize() * 1.01,
		                                    static_cast<Sample3D *>( tree[0] )->getYSize() * 1.01,
		                                    static_cast<Sample3D *>( tree[0] )->getZSize() * 1.01, np, tree[0]->getCenter() ) ;

		while( meshPoints.size() < shuffled.size() )
		{
			Point ptest( shuffled[p].first->getX(), shuffled[p].first->getY(), shuffled[p].first->getZ() ) ;

			if( ( visited[p] ) && ( shufflingGrid->pixels[ix][iy][iz]->coOccur( ptest ) ) )
			{
				visited[p] = false ;
				meshPoints.push_back( shuffled[p] ) ;
				ix++ ;

				if( ix == shufflingGrid->getLengthX() )
				{
					ix = 0 ;
					iy++ ;

					if( iy == shufflingGrid->getLengthY() )
					{
						iy = 0 ;
						iz++ ;

						if( iz == shufflingGrid->getLengthY() )
							iz = 0 ;
					}
				}

				p = 0 ;
			}
			else
			{
				p++ ;

				if( p == shuffled.size() )
				{
					p = 0 ;
					ix++ ;

					if( ix == shufflingGrid->getLengthX() )
					{
						ix = 0 ;
						iy++ ;

						if( iy == shufflingGrid->getLengthY() )
						{
							iy = 0 ;
							iz++ ;

							if( iz == shufflingGrid->getLengthZ() )
								iz = 0 ;
						}
					}
				}
			}
		}
	}

	std::random_shuffle( meshPoints.begin(), meshPoints.end() );
	std::cout << "done... " << std::endl ;
}

void FeatureTree::homothety(double before, double now, double after)
{
	std::valarray<bool> nodes(getDisplacements().size()/2) ;
	if(nodes.size() == 0)
		return ;
	nodes = false ;
	std::vector<DelaunayTriangle *> tri = getElements2D() ;
	if(tri[0]->timePlanes() != 2)
		return ;
	double stepa = after/now ;
	double stepb = now/before ;
	for(size_t i = 0 ; i < tri.size() ; i++)
	{
		for(size_t j = 0 ; j < tri[i]->getBoundingPoints().size()/2 ; j++)
		{
			size_t id = tri[i]->getBoundingPoint(j).getId() ;
			if(!nodes[id])
			{
				nodes[id] = true ;
				tri[i]->getBoundingPoint(j).getX() = tri[i]->getBoundingPoint(tri[i]->getBoundingPoints().size()/2+j).getX() ;
				tri[i]->getBoundingPoint(j).getY() = tri[i]->getBoundingPoint(tri[i]->getBoundingPoints().size()/2+j).getY() ;
			}
		}
		for(size_t j = tri[i]->getBoundingPoints().size()/2 ; j < tri[i]->getBoundingPoints().size() ; j++)
		{
			size_t id = tri[i]->getBoundingPoint(j).getId() ;
			if(!nodes[id])
			{
				nodes[id] = true ;
				if(tri[i]->getBoundingPoint(j).getX() != 0.5 && tri[i]->getBoundingPoint(j).getX() != -0.5 && tri[i]->getBoundingPoint(j).getY() != 0.5 && tri[i]->getBoundingPoint(j).getY() != -0.5)
				{
					double r = sqrt(tri[i]->getBoundingPoint(j).getX()*tri[i]->getBoundingPoint(j).getX() + tri[i]->getBoundingPoint(j).getY()*tri[i]->getBoundingPoint(j).getY()) ;
					if(r < now)
					{
						tri[i]->getBoundingPoint(j).getX() /= now ;
						tri[i]->getBoundingPoint(j).getY() /= now ;
						tri[i]->getBoundingPoint(j).getX() *= after ;
						tri[i]->getBoundingPoint(j).getY() *= after ;
					}
					else
					{
						Point inter(tri[i]->getBoundingPoint(j).getX(), tri[i]->getBoundingPoint(j).getY()) ;
						dynamic_cast<Rectangle *>(tree[0])->project(&inter) ;
						double rl = sqrt(inter.getX()*inter.getX() + inter.getY()*inter.getY()) ;
						tri[i]->getBoundingPoint(j).getX() *= 1.-(1.-stepa)*(rl-r)/(rl-now) ;
						tri[i]->getBoundingPoint(j).getY() *= 1.-(1.-stepa)*(rl-r)/(rl-now) ;
					}
				}
			}
		}
		tri[i]->clearElementaryMatrix() ;
	}
//	tri[0]->getBoundingPoint(0).print() ;
}
