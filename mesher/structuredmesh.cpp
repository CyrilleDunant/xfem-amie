// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010-2011

#include "structuredmesh.h"

using namespace Amie ;

StructuredMesh::StructuredMesh(double sizeX, double sizeY, int div, const Point & center ): Mesh< Amie::DelaunayTriangle, Amie::DelaunayTreeItem >(SPACE_TWO_DIMENSIONAL), grid( sizeX,  sizeY,  div, center)
{
	global_counter = 0 ;
	for(size_t i = 0 ; i < grid.pixels.size()+1 ; i++)
	{
		for(size_t j = 0 ; j < grid.pixels[0].size()+1 ; j++)
		{
			points.push_back(new Point(center.getX()-sizeX*.5 +(sizeX)*i/(grid.pixels.size()), (center.getY()-sizeY*.5) + (sizeY)*j/(grid.pixels[0].size()))) ;
			if(i != 0 && i!= grid.pixels.size() && j != 0 && j != grid.pixels[0].size())
			{
				double dx = .2*(2.*((double)rand()/(double)RAND_MAX)-1.)*(sizeX)/(grid.pixels.size()) ;
				double dy = .2*(2.*((double)rand()/(double)RAND_MAX)-1.)*(sizeY)/(grid.pixels[0].size()) ;
				*points[global_counter] += Point(dx, dy) ;
			}
			points.back()->getId() = global_counter++ ;
		}
	}
	
	for(size_t i = 0 ; i < grid.pixels.size() ; i++)
	{
		for(size_t j = 0 ; j < grid.pixels[0].size() ; j++)
		{
			new DelaunayTriangle(this,nullptr, points[(grid.pixels[0].size()+1)*i+j], points[(grid.pixels[0].size()+1)*i+j+1], points[(grid.pixels[0].size()+1)*(i+1)+j+1], nullptr) ;
			new DelaunayTriangle(this,nullptr, points[(grid.pixels[0].size()+1)*i+j], points[(grid.pixels[0].size()+1)*(i+1)+j], points[(grid.pixels[0].size()+1)*(i+1)+j+1], nullptr) ;
		}
	}
	
	
	for(size_t i = 1 ; i < grid.pixels.size()-1 ; i++)
	{
		for(size_t j = 1 ; j < grid.pixels[i].size()-1 ; j++)
		{
			std::vector<int> idx ;
			idx.push_back(grid.pixels[i].size()*(i-1)+j-1) ;
			idx.push_back(grid.pixels[i].size()*(i-1)+j) ;
			idx.push_back(grid.pixels[i].size()*(i-1)+j+1) ;
			idx.push_back(grid.pixels[i].size()*(i)+j-1) ;
			idx.push_back(grid.pixels[i].size()*(i)+j) ;
			idx.push_back(grid.pixels[i].size()*(i)+j+1) ;
			idx.push_back(grid.pixels[i].size()*(i+1)+j-1) ;
			idx.push_back(grid.pixels[i].size()*(i+1)+j) ;
			idx.push_back(grid.pixels[i].size()*(i+1)+j+1) ;
			
			std::vector<DelaunayTriangle *> tris ;
			for(size_t k = 0 ; k < idx.size() ; k++)
			{
				tris.push_back(static_cast<DelaunayTriangle *>(tree[idx[k]*2])) ;
				tris.push_back(static_cast<DelaunayTriangle *>(tree[idx[k]*2+1])) ;
			}
			
			for(size_t k = 0 ; k < tris.size() ; k++)
			{
				for(size_t l = k+1 ; l < tris.size() ; l++)
				{
					if(tris[k]->isNeighbour(tris[l]))
						makeNeighbours(tris[k], tris[l]);
					if(tris[k]->isInNeighbourhood(tris[l]))
					{
						tris[k]->addNeighbourhood(tris[l]) ;
						tris[l]->addNeighbourhood(tris[k]) ;
					}
				}
			}
		}
	}
	
}

size_t StructuredMesh::getLastNodeId() const
{
	return global_counter ;
}

StructuredMesh::~StructuredMesh() 
{
	std::valarray<Point *> nularray(0) ;
	for(size_t i = 0 ;  i < this->tree.size() ; i++)
	{

		DelaunayTriangle * t = static_cast<DelaunayTriangle *>(tree[i]) ;
			
		t->setBoundingPoints(nularray) ;
	}
	
	for(size_t i = 0 ;  i < this->tree.size() ; i++)
	{
		delete this->tree[i] ;
	}

	for(size_t i = 0 ; i < points.size() ; i++)
		delete points[i] ;
	
}

std::vector<DelaunayTriangle *> StructuredMesh::getElements() 
{
	std::vector<DelaunayTriangle *> triangles ;
	for(size_t i = 0 ; i < tree.size() ; i++)
		triangles.push_back( static_cast<DelaunayTriangle *>(tree[i])) ;
	return triangles ;
}

std::vector<DelaunayTriangle *> StructuredMesh::getConflictingElements(const Point  * p) 
{
	std::vector<DelaunayTriangle *> ret ;
	double startX = grid.getX()*.5-grid.getCenter().getX() + p->getX() ;
	int startI = std::max(0., startX/grid.getPixelSize() - 2) ;
	
	double endX =  startX+.05*grid.getX();
	int endI = std::min(endX/grid.getPixelSize() + 2, (double)grid.pixels.size());
	
	double startY = grid.getY()*.5-grid.getCenter().getY() + p->getY() ;
	int startJ = std::max(0., startY/grid.getPixelSize() - 2) ;
	
	double endY =  startY+.05*grid.getY();
	int endJ = std::min(endY/grid.getPixelSize() + 2, (double)grid.pixels[0].size());
	
	for(int i = startI ; i < endI ; i++)
	{
		for(int j = startJ ; j < endJ ; j++)
		{
			Point proj0(*p) ;
			Point proj1(*p) ;
			static_cast<DelaunayTriangle *>(tree[i*2*grid.pixels[0].size()+j*2])->project(&proj0) ;
			static_cast<DelaunayTriangle *>(tree[i*2*grid.pixels[0].size()+j*2+1])->project(&proj1) ;
			if(tree[i*2*grid.pixels[0].size()+j*2]->isAlive())
			{
				if(static_cast<DelaunayTriangle *>(tree[i*2*grid.pixels[0].size()+j*2])->in(*p) || dist(proj0, *p) < 128.*POINT_TOLERANCE_2D)
					ret.push_back(static_cast<DelaunayTriangle *>(tree[i*2*grid.pixels[0].size()+j*2])) ;
			}
			
			if(tree[i*2*grid.pixels[0].size()+j*2+1]->isAlive())
			{
				if(static_cast<DelaunayTriangle *>(tree[i*2*grid.pixels[0].size()+j*2+1])->in(*p)|| dist(proj1, *p) < 128.*POINT_TOLERANCE_2D)
					ret.push_back(static_cast<DelaunayTriangle *>(tree[i*2*grid.pixels[0].size()+j*2+1])) ;
			}
		}
	}
	
	std::stable_sort(ret.begin(), ret.end());
	auto e = std::unique(ret.begin(), ret.end()) ;
	ret.erase(e, ret.end()) ;
	return ret ;
}
std::vector<DelaunayTriangle *> StructuredMesh::getConflictingElements(const Geometry * geo) 
{
	std::vector<DelaunayTriangle *> ret ;
	double startX = grid.getX()*.5-grid.getCenter().getX() + geo->getCenter().getX()-2.*geo->getRadius() ;
	int startI = std::max(0., startX/grid.getPixelSize() - 2) ;
	
	double endX =  startX+2.*geo->getRadius();
	int endI = std::min(endX/grid.getPixelSize() + 2, (double)grid.pixels.size());
	
	double startY = grid.getY()*.5-grid.getCenter().getY() + geo->getCenter().getY()-2.*geo->getRadius() ;
	int startJ = std::max(0., startY/grid.getPixelSize() - 2) ;
		
	double endY =  startY+2.*geo->getRadius();
	int endJ = std::min(endY/grid.getPixelSize() + 2, (double)grid.pixels[0].size());
	for(int i = startI ; i < endI ; i++)
	{
		for(int j = startJ ; j < endJ ; j++)
		{
// 	for(int i = 0 ; i < grid.pixels.size() ; i++)
// 	{
// 		for(int j = 0 ; j < grid.pixels[0].size() ; j++)
// 		{
			if(tree[i*2*grid.pixels[0].size()+j*2]->isAlive())
			{
				if(static_cast<DelaunayTriangle *>(tree[i*2*grid.pixels[0].size()+j*2])->intersects(geo) 
					|| geo->in(*tree[i*2*grid.pixels[0].size()+j*2]->first)
					|| geo->in(*tree[i*2*grid.pixels[0].size()+j*2]->second)
					|| geo->in(*tree[i*2*grid.pixels[0].size()+j*2]->third)
					|| geo->in(static_cast<DelaunayTriangle *>(tree[i*2*grid.pixels[0].size()+j*2])->getCenter())
					|| static_cast<DelaunayTriangle *>(tree[i*2*grid.pixels[0].size()+j*2])->in(geo->getCenter())
					)
					ret.push_back(static_cast<DelaunayTriangle *>(tree[i*2*grid.pixels[0].size()+j*2])) ;
			}
			if(tree[i*2*grid.pixels[0].size()+j*2+1]->isAlive())
			{
				if(static_cast<DelaunayTriangle *>(tree[i*2*grid.pixels[0].size()+j*2+1])->intersects(geo) 
					|| geo->in(*tree[i*2*grid.pixels[0].size()+j*2+1]->first)
					|| geo->in(*tree[i*2*grid.pixels[0].size()+j*2+1]->second)
					|| geo->in(*tree[i*2*grid.pixels[0].size()+j*2+1]->third)
					|| geo->in(static_cast<DelaunayTriangle *>(tree[i*2*grid.pixels[0].size()+j*2+1])->getCenter())
					|| static_cast<DelaunayTriangle *>(tree[i*2*grid.pixels[0].size()+j*2+1])->in(geo->getCenter())
					)
					ret.push_back(static_cast<DelaunayTriangle *>(tree[i*2*grid.pixels[0].size()+j*2+1])) ;
			}
		}
	}
	return ret ;
}
void StructuredMesh::setElementOrder(Order o, double dt)
{
    if(allElementsCacheID != -1)
    {
        caches[allElementsCacheID].clear() ;
        coefs[allElementsCacheID].clear() ;
        allElementsCacheID = -1 ;
    }
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
void StructuredMesh::insert(Point *) 
{
}

void StructuredMesh::addSharedNodes(size_t nodes_per_side, size_t time_planes, double timestep)
{

    std::valarray<bool> visited(false, size()) ;
	for(auto & i : tree)
	{
            delete i->cachedGps ;
            i->cachedGps = nullptr ;
            i->behaviourUpdated = true ;
		visited[i->index] = true ;
			
		size_t nodes_per_plane = nodes_per_side*3+3 ;
		
		std::valarray<Point *> newPoints(nodes_per_plane*time_planes) ;
		std::valarray<bool> done(false, nodes_per_plane*time_planes) ;
		
		for(size_t plane = 0 ; plane < time_planes ; plane++)
		{
			for(size_t side = 0 ; side < 3 ; side++)
			{
				Point a(i->getBoundingPoint(side)) ;
				Point b(i->getBoundingPoint((side+1)%3)) ;
				
				if(time_planes> 1)
				{
					a.getT() = (double)plane*(timestep/(double)(time_planes-1))-timestep/2.;
					b.getT() = (double)plane*(timestep/(double)(time_planes-1))-timestep/2.;
				}
				for(size_t node = 0 ; node < nodes_per_side+1 ; node++)
				{
					double fraction = (double)(node)/((double)nodes_per_side+1) ;
					Point proto = a*(1.-fraction) + b*fraction ;
					Point * foundPoint = nullptr ;
					
					for(size_t j = 0 ; j< i->getBoundingPoints().size() ; j++)
					{
						if(i->getBoundingPoint(j) == proto)
						{
							foundPoint = &i->getBoundingPoint(j) ;
							break ;
						}
					}
					
					if(!foundPoint)
					{
                        std::vector<DelaunayTriangle *> neighbourhood = getNeighbourhood(i) ;
						for(auto & n : neighbourhood)
						{
							if(visited[n->index])
							{
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
							newPoints[nodes_per_plane*plane+side*(nodes_per_side+1)+node]->getId() = global_counter++ ;
						}
						
						done[nodes_per_plane*plane+side*(nodes_per_side+1)+node] = true ;
					}
				}
			}
		}
		
		i->setBoundingPoints(newPoints) ;
	}


}

MicDerivedMesh::MicDerivedMesh(const char * voxelSource, std::map<unsigned char,Form *> behaviourMap) : Mesh<DelaunayTetrahedron, DelaunayTreeItem3D>(SPACE_THREE_DIMENSIONAL)
{
    VoxelFilter vx ;
    vx.behaviourMap = behaviourMap ;
    vx.read(voxelSource, this);
    tree  = vx.getElements() ;
//     additionalPoints = vx.getPoints() ;
    global_counter = vx.getPoints().size() ;
    generateCache() ;
    
    for(size_t i = 0 ; i < tree.size() ; i++)
    {
        tree[i]->tree = this ;
    }
}

void MicDerivedMesh::extrude(double dt)
{
    std::map<Point *, Point *> points ;

    std::vector<DelaunayTetrahedron *> tri = getElements() ;
    double beginning = tri[0]->getBoundingPoint(0).getT() ;
    double end = tri[0]->getBoundingPoint(0).getT() ;
    for(size_t i = 1 ; i < tri[0]->getBoundingPoints().size() ; i++)
    {
        if(tri[0]->getBoundingPoint(i).getT() < beginning)
            beginning = tri[0]->getBoundingPoint(i).getT() ;
        if(tri[0]->getBoundingPoint(i).getT() > end)
            end = tri[0]->getBoundingPoint(i).getT() ;
    }

    int indexOfLastTimePlane = (tri[0]->timePlanes()-1)*tri[0]->getBoundingPoints().size()/tri[0]->timePlanes() ;
    int pointsPerTimePlane = tri[0]->getBoundingPoints().size()/tri[0]->timePlanes() ;

    for(size_t i = 0 ; i < tri.size() ; i++)
    {
        for(size_t j = 0 ; j < tri[i]->getBoundingPoints().size() ; j++)
        {
            Point * next = new Point(tri[i]->getBoundingPoint(j).getX(), tri[i]->getBoundingPoint(j).getY(), tri[i]->getBoundingPoint(j).getZ()) ;
            additionalPoints.push_back(next);
            next->getT() = tri[i]->getBoundingPoint(j).getT() ;
            next->getT() = end + dt * (next->getT() - beginning) / (end - beginning) ;
            bool increment = true ;
            if(std::abs(next->getT() - end) < POINT_TOLERANCE_2D)
            {
                next = &tri[i]->getBoundingPoint(j+indexOfLastTimePlane) ;
                increment = false ;
//              next->print() ;
            }
            if(increment && !points.find(&tri[i]->getBoundingPoint(j))->second)
            {
                next->setId(global_counter++) ;
            }
            points.insert(std::pair<Point *, Point *>(&tri[i]->getBoundingPoint(j), next)) ;
        }
    }

    std::map<Point *, Point *>::iterator finder ;
    for(size_t i = 0 ; i < tri.size() ; i++)
    {
        Point * a = points.find(&tri[i]->getBoundingPoint(0))->second ;
        Point * b = points.find(&tri[i]->getBoundingPoint(pointsPerTimePlane/4))->second ;
        Point * c = points.find(&tri[i]->getBoundingPoint(2*pointsPerTimePlane/4))->second ;
        Point * d = points.find(&tri[i]->getBoundingPoint(3*pointsPerTimePlane/4))->second ;

        std::valarray<Point *> newPoints(tri[i]->getBoundingPoints().size()) ;
        for(size_t j = 0 ; j < newPoints.size() ; j++)
        {
            newPoints[j] = points.find(&tri[i]->getBoundingPoint(j))->second ;
//          newPoints[j]->print() ;
        }

        DelaunayTetrahedron * toInsert = new DelaunayTetrahedron(this, nullptr, a,b,c,d, a) ;
        toInsert->setOrder(tri[i]->getOrder()) ;
        toInsert->setBoundingPoints(newPoints) ;
        toInsert->setBehaviour(this, tri[i]->getBehaviour()->getCopy()) ;
    }
}

int MicDerivedMesh::addToTree ( DelaunayTreeItem3D * toAdd ) 
{ 
    tree.push_back(static_cast<DelaunayTetrahedron *>(toAdd)); 
    return tree.size()-1 ;
    
} ;

void MicDerivedMesh::addSharedNodes( size_t nodes_per_side, size_t time_planes, double timestep, const TetrahedralElement *father )
{

    std::vector<DelaunayTetrahedron *> tet = getElements() ;
    std::valarray<bool> visited(false, tree.size()) ;

    if( nodes_per_side > 1 )
        nodes_per_side = 1 ;

    std::vector<size_t> positions ;

    if( nodes_per_side )
    {
        positions.push_back( 0 ) ;
        positions.push_back( 1 ) ;
        positions.push_back( 2 ) ;

        positions.push_back( 2 ) ;
        positions.push_back( 3 ) ;
        positions.push_back( 4 ) ;

        positions.push_back( 4 ) ;
        positions.push_back( 5 ) ;
        positions.push_back( 6 ) ;

        positions.push_back( 0 ) ;
        positions.push_back( 7 ) ;
        positions.push_back( 6 ) ;

        positions.push_back( 6 ) ;
        positions.push_back( 8 ) ;
        positions.push_back( 2 ) ;

        positions.push_back( 0 ) ;
        positions.push_back( 9 ) ;
        positions.push_back( 4 ) ;
    }
    else
    {
        positions.push_back( 0 ) ;
        positions.push_back( 1 ) ;

        positions.push_back( 1 ) ;
        positions.push_back( 2 ) ;

        positions.push_back( 2 ) ;
        positions.push_back( 3 ) ;

        positions.push_back( 3 ) ;
        positions.push_back( 0 ) ;

        positions.push_back( 3 ) ;
        positions.push_back( 1 ) ;

        positions.push_back( 0 ) ;
        positions.push_back( 2 ) ;
    }

    for( size_t i = 0 ; i < tet.size() ; i++ )
    {
        delete tet[i]->cachedGps ;
        tet[i]->cachedGps = nullptr ;
        tet[i]->behaviourUpdated = true ;
//      if( tet[i]->volume() < 0 )
//      {
//          for( int j = 0 ; j < tet[i]->getBoundingPoints().size() ; j++ )
//              tet[i]->getBoundingPoint( j ).print() ;
//
//          exit( 0 ) ;
//      }

        std::vector<std::pair<Point, Point> > sides ;
        sides.push_back( std::make_pair( tet[i]->getBoundingPoint( 0 ), tet[i]->getBoundingPoint( 1 ) ) ) ;
        sides.push_back( std::make_pair( tet[i]->getBoundingPoint( 1 ), tet[i]->getBoundingPoint( 2 ) ) ) ;
        sides.push_back( std::make_pair( tet[i]->getBoundingPoint( 2 ), tet[i]->getBoundingPoint( 3 ) ) ) ;
        sides.push_back( std::make_pair( tet[i]->getBoundingPoint( 3 ), tet[i]->getBoundingPoint( 0 ) ) ) ;
        sides.push_back( std::make_pair( tet[i]->getBoundingPoint( 3 ), tet[i]->getBoundingPoint( 1 ) ) ) ;
        sides.push_back( std::make_pair( tet[i]->getBoundingPoint( 0 ), tet[i]->getBoundingPoint( 2 ) ) ) ;


        visited[tet[i]->index] = true ;

        size_t nodes_per_plane = nodes_per_side * 6 + 4 ;

        std::valarray<Point *> newPoints( ( Point * )nullptr, nodes_per_plane * time_planes ) ;
        std::valarray<bool> done( false, nodes_per_plane * time_planes ) ;

        for( size_t plane = 0 ; plane < time_planes ; plane++ )
        {
            size_t current = 0 ;

            for( size_t side = 0 ; side < 6 ; side++ )
            {
                Point a( sides[side].first ) ;
                Point b( sides[side].second ) ;

                if( time_planes > 1 )
                {
                    a.getT() = ( double )plane * ( timestep / ( double )( time_planes - 1 ) ) - timestep / 2.;
                    b.getT() = ( double )plane * ( timestep / ( double )( time_planes - 1 ) ) - timestep / 2.;
                }

                for( size_t node = 0 ; node < nodes_per_side + 2 ; node++ )
                {
                    double fraction = ( double )( node ) / ( ( double )nodes_per_side + 1 ) ;
                    Point proto = a * ( 1. - fraction ) + b * fraction ;
                    Point *foundPoint = nullptr ;

                    for( size_t j = 0 ; j < tet[i]->getBoundingPoints().size() ; j++ )
                    {
                        if( tet[i]->getBoundingPoint( j ) == proto )
                        {
                            foundPoint = &tet[i]->getBoundingPoint( j ) ;
                            break ;
                        }
                    }

                    if( !foundPoint )
                    {
                        for( size_t j = 0 ; j < tet[i]->neighbourhood.size() ; j++ )
                        {
                            if( visited[tet[i]->neighbourhood[j]] )
                            {
                                for( size_t k = 0 ; k < tet[i]->getNeighbourhood( j )->getBoundingPoints().size() ; k++ )
                                {
                                    if( tet[i]->getNeighbourhood( j )->getBoundingPoint( k ) == proto )
                                    {
                                        foundPoint = &tet[i]->getNeighbourhood( j )->getBoundingPoint( k ) ;
                                        break ;
                                    }
                                }

                                if( foundPoint )
                                {
                                    break ;
                                }
                            }
                        }
                    }

                    if( !done[positions[current] + plane * nodes_per_plane] )
                    {
                        if( foundPoint )
                        {
                            newPoints[positions[current] + plane * nodes_per_plane]  = foundPoint ;
                        }
                        else
                        {
                            additionalPoints.push_back( new Point( proto ) ) ;
                            newPoints[positions[current] + plane * nodes_per_plane]  = additionalPoints.back() ;
                            newPoints[positions[current] + plane * nodes_per_plane]->getId() = global_counter++ ;
                        }

                        done[positions[current] + plane * nodes_per_plane] = true ;
                    }

                    current++ ;
                }
            }
        }

        for( size_t k = 0 ; k < newPoints.size() ; k++ )
            if( !newPoints[k] )
            {
                std::cout << "ouch !" << std::endl ;

                for( size_t k = 0 ; k < newPoints.size() ; k++ )
                    if( newPoints[k] )
                        newPoints[k]->print() ;

                exit( 0 ) ;
            }

        tet[i]->setBoundingPoints( newPoints ) ;

    }

}


void MicDerivedMesh::extrude(const Vector & dt)
{
    std::map<Point *, std::vector<Point *> > points ;
    std::map<Point *, Point *> pointsInTriangle ;
    std::map<Point *, std::vector<Point *> >::iterator finder ;

    std::vector<DelaunayTetrahedron *> tri = getElements() ;
    double beginning = tri[0]->getBoundingPoint(0).getT() ;
    double end = tri[0]->getBoundingPoint(0).getT() ;
    for(size_t i = 1 ; i < tri[0]->getBoundingPoints().size() ; i++)
    {
        if(tri[0]->getBoundingPoint(i).getT() < beginning)
            beginning = tri[0]->getBoundingPoint(i).getT() ;
        if(tri[0]->getBoundingPoint(i).getT() > end)
            end = tri[0]->getBoundingPoint(i).getT() ;
    }

    int indexOfLastTimePlane = (tri[0]->timePlanes()-1)*tri[0]->getBoundingPoints().size()/tri[0]->timePlanes() ;
    int pointsPerTimePlane = tri[0]->getBoundingPoints().size()/tri[0]->timePlanes() ;

    for(size_t i = 0 ; i < tri.size() ; i++)
    {
        for(size_t j = 0 ; j < tri[i]->getBoundingPoints().size() ; j++)
        {
            Point * current = &tri[i]->getBoundingPoint(j) ;
            current->getT() = dt[0]+(dt[1]-dt[0])*(current->getT() - beginning)/(end-beginning) ;
//          tri[i]->setBoundingPoint(j, current) ;
        }
    }


    for(size_t i = 0 ; i < tri.size() ; i++)
    {
        for(size_t j = 0 ; j < tri[i]->getBoundingPoints().size() ; j++)
        {
            Point * current = &tri[i]->getBoundingPoint(j) ;
            finder = points.find(current) ;
            if(finder == points.end())
            {
                Point * current = &tri[i]->getBoundingPoint(j) ;
                current->getT() = dt[0]+(dt[1]-dt[0])*(current->getT() - beginning)/(end-beginning) ;

                std::vector<Point *> newPoints ;
                if( (int)j < indexOfLastTimePlane)
                {
                    if( (int)j < pointsPerTimePlane )
                    {
                        newPoints.push_back(&tri[i]->getBoundingPoint(indexOfLastTimePlane+j)) ;
                    }
                    size_t ifirst = newPoints.size() ;
                    for(size_t k = ifirst+1 ; k < dt.size()-1 ; k++)
                    {
                        Point * next = new Point(current->getX(), current->getY(), current->getZ()) ;
                        additionalPoints.push_back(next);
                        next->getT() = dt[k]+(dt[k+1]-dt[k])*(current->getT()-beginning)/(end-beginning) ;
                        next->setId(global_counter++) ;
                        newPoints.push_back(next) ;
                    }
                }
                else
                {
                    size_t jp = j - indexOfLastTimePlane ;
                    Point * previous = &tri[i]->getBoundingPoint(jp) ;
                    std::vector<Point *> previousPoints = points.find(previous)->second ;
                    for(size_t k = 1 ; k < previousPoints.size() ; k++)
                        newPoints.push_back(previousPoints[k]) ;
                    Point * next = new Point(current->getX(), current->getY(), current->getZ()) ;
                    size_t k = dt.size()-2 ;
                    next->getT() = dt[k]+(dt[k+1]-dt[k])*(current->getT()-beginning)/(end-beginning) ;
                    next->setId(global_counter++) ;
                    newPoints.push_back(next) ;
                }

                points.insert(std::pair<Point *, std::vector<Point *> >(current, newPoints)) ;
            }
        }
    }

    for(size_t i = 0 ; i < tri.size() ; i++)
    {
        std::valarray<Point *> newPoints(tri[i]->getBoundingPoints().size()) ;
        for(size_t k = 1 ; k < dt.size()-1 ; k++)
        {
            for(size_t j = 0 ; j < newPoints.size() ; j++)
            {
                std::vector<Point *> found = points.find(&tri[i]->getBoundingPoint(j))->second ;
                newPoints[j] = found[k-1] ;
            }

            DelaunayTetrahedron * toInsert = new DelaunayTetrahedron(this, nullptr, newPoints[0],newPoints[pointsPerTimePlane/4],newPoints[pointsPerTimePlane*2/4],newPoints[pointsPerTimePlane*3/4],newPoints[0]) ;
            toInsert->setOrder(tri[i]->getOrder()) ;
            toInsert->setBoundingPoints(newPoints) ;
            toInsert->setBehaviour(tri[i]->getState().getMesh3D(), tri[i]->getBehaviour()->getCopy()) ;
        }
    }

    std::cout << tri.size() << "\t" << tri.size() << std::endl ;

}

    
std::vector<DelaunayTetrahedron *> MicDerivedMesh::getNeighbourhood(DelaunayTetrahedron * element) const
{
    std::vector<DelaunayTetrahedron *> ret ;
    for(const auto & idx : element->neighbourhood)
    {
        ret.push_back((DelaunayTetrahedron *)tree[idx]);
    }
    return ret ;
};

std::vector<DelaunayTetrahedron *> MicDerivedMesh::getConflictingElements ( const Point  * p )  
{
    std::vector<DelaunayTetrahedron *> ret ;
    for(auto & e : tree)
        if(e->inCircumSphere(*p))
            ret.push_back(e);
    return ret ;
};

std::vector<DelaunayTetrahedron *> MicDerivedMesh::getConflictingElements ( const Geometry * g )
{
    std::vector<DelaunayTetrahedron *> ret ;
    for(auto & e : tree)
        if(g->in(e->getCenter()) || e->in(g->getCenter()))
            ret.push_back(e);
    return ret ;
} ;
    

void MicDerivedMesh::setElementOrder ( Order elemOrder, double dt  ){
    if(allElementsCacheID != -1)
    {
        caches[allElementsCacheID].clear() ;
        coefs[allElementsCacheID].clear() ;
        allElementsCacheID = -1 ;
    }
    switch( elemOrder )
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
        addSharedNodes( 1, 1, 0 ) ;
        break ;
    }
    case CUBIC:
    {
        addSharedNodes( 2, 1, 0 ) ;
        break ;
    }
    case QUADRIC:
    {
        addSharedNodes( 3, 1, 0 ) ;
        break ;
    }
    case QUINTIC:
    {
        addSharedNodes( 3, 1, 0 ) ;
        break ;
    }
    case CONSTANT_TIME_LINEAR:
    {
        addSharedNodes( 0, 2, dt ) ;
        break ;
    }
    case CONSTANT_TIME_QUADRATIC:
    {
        addSharedNodes( 0, 3, dt ) ;
        break ;
    }
    case LINEAR_TIME_LINEAR:
    {
        addSharedNodes( 0, 2, dt ) ;
        break ;
    }
    case LINEAR_TIME_QUADRATIC:
    {
        addSharedNodes( 0, 3, dt ) ;
        break ;
    }
    case QUADRATIC_TIME_LINEAR:
    {
        addSharedNodes( 1, 2, dt ) ;
        break ;
    }
    case QUADRATIC_TIME_QUADRATIC:
    {
        addSharedNodes( 1, 3, dt ) ;
        break ;
    }
    case CUBIC_TIME_LINEAR:
    {
        addSharedNodes( 2, 2, dt ) ;
        break ;
    }
    case CUBIC_TIME_QUADRATIC:
    {
        addSharedNodes( 2, 3, dt ) ;
        break ;
    }
    case QUADRIC_TIME_LINEAR:
    {
        addSharedNodes( 3, 2, dt ) ;
        break ;
    }
    case QUADRIC_TIME_QUADRATIC:
    {
        addSharedNodes( 3, 3, dt ) ;
        break ;
    }
    case QUINTIC_TIME_LINEAR:
    {
        addSharedNodes( 3, 2, dt ) ;
        break ;
    }
    case QUINTIC_TIME_QUADRATIC:
    {
        addSharedNodes( 3, 3, dt ) ;
        break ;
    }
    default:
        break ;

    }
}
