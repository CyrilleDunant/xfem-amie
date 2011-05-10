// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010-2011

#include "structuredmesh.h"

using namespace Mu ;

StructuredMesh::StructuredMesh(double sizeX, double sizeY, int div, const Point & center ): grid( sizeX,  sizeY,  div, center)
{
	global_counter = 0 ;
	for(size_t i = 0 ; i < grid.pixels.size()+1 ; i++)
	{
		for(size_t j = 0 ; j < grid.pixels[0].size()+1 ; j++)
		{
			points.push_back(new Point(center.x-sizeX*.5 +(sizeX)*i/(grid.pixels.size()), (center.y-sizeY*.5) + (sizeY)*j/(grid.pixels[0].size()))) ;
			if(i != 0 && i!= grid.pixels.size() && j != 0 && j != grid.pixels[0].size())
			{
				double dx = .2*(2.*((double)rand()/(double)RAND_MAX)-1.)*(sizeX)/(grid.pixels.size()) ;
				double dy = .2*(2.*((double)rand()/(double)RAND_MAX)-1.)*(sizeY)/(grid.pixels[0].size()) ;
				*points[global_counter] += Point(dx, dy) ;
			}
			points.back()->id = global_counter++ ;
		}
	}
	
	for(size_t i = 0 ; i < grid.pixels.size() ; i++)
	{
		for(size_t j = 0 ; j < grid.pixels[0].size() ; j++)
		{
			new DelaunayTriangle(this,NULL, points[(grid.pixels[0].size()+1)*i+j], points[(grid.pixels[0].size()+1)*i+j+1], points[(grid.pixels[0].size()+1)*(i+1)+j+1], NULL) ;
			new DelaunayTriangle(this,NULL, points[(grid.pixels[0].size()+1)*i+j], points[(grid.pixels[0].size()+1)*(i+1)+j], points[(grid.pixels[0].size()+1)*(i+1)+j+1], NULL) ;
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

const size_t & StructuredMesh::getLastNodeId() const
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
std::vector<DelaunayTriangle *> StructuredMesh::getConflictingElements(const Point  * p) const
{
	std::vector<DelaunayTriangle *> ret ;
	double startX = grid.getX()*.5-grid.getCenter().x + p->x ;
	int startI = std::max(0., startX/grid.getPixelSize() - 2) ;
	
	double endX =  startX+.05*grid.getX();
	int endI = std::min(endX/grid.getPixelSize() + 2, (double)grid.pixels.size());
	
	double startY = grid.getY()*.5-grid.getCenter().y + p->y ;
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
std::vector<DelaunayTriangle *> StructuredMesh::getConflictingElements(const Geometry * geo) const
{
	std::vector<DelaunayTriangle *> ret ;
	double startX = grid.getX()*.5-grid.getCenter().x + geo->getCenter().x-2.*geo->getRadius() ;
	int startI = std::max(0., startX/grid.getPixelSize() - 2) ;
	
	double endX =  startX+2.*geo->getRadius();
	int endI = std::min(endX/grid.getPixelSize() + 2, (double)grid.pixels.size());
	
	double startY = grid.getY()*.5-grid.getCenter().y + geo->getCenter().y-2.*geo->getRadius() ;
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

	for(auto i = tree.begin() ; i != tree.end() ; ++i)
	{
		
		(*i)->visited = true ;
			
		size_t nodes_per_plane = nodes_per_side*3+3 ;
		
		std::valarray<Point *> newPoints(nodes_per_plane*time_planes) ;
		std::valarray<bool> done(false, nodes_per_plane*time_planes) ;
		
		for(size_t plane = 0 ; plane < time_planes ; plane++)
		{
			for(size_t side = 0 ; side < 3 ; side++)
			{
				Point a(static_cast<DelaunayTriangle *>(*i)->getBoundingPoint(side)) ;
				Point b(static_cast<DelaunayTriangle *>(*i)->getBoundingPoint((side+1)%3)) ;
				
				if(time_planes> 1)
				{
					a.t = (double)plane*(timestep/(double)(time_planes-1))-timestep/2.;
					b.t = (double)plane*(timestep/(double)(time_planes-1))-timestep/2.;
				}
				for(size_t node = 0 ; node < nodes_per_side+1 ; node++)
				{
					double fraction = (double)(node)/((double)nodes_per_side+1) ;
					Point proto = a*(1.-fraction) + b*fraction ;
					Point * foundPoint = NULL ;
					
					for(size_t j = 0 ; j< static_cast<DelaunayTriangle *>(*i)->getBoundingPoints().size() ; j++)
					{
						if(static_cast<DelaunayTriangle *>(*i)->getBoundingPoint(j) == proto)
						{
							foundPoint = &static_cast<DelaunayTriangle *>(*i)->getBoundingPoint(j) ;
							break ;
						}
					}
					
					if(!foundPoint)
					{
						for(size_t j = 0 ; j < static_cast<DelaunayTriangle *>(*i)->neighbourhood.size() ; j++)
						{
							if(static_cast<DelaunayTriangle *>(*i)->getNeighbourhood(j)->visited)
							{
								DelaunayTriangle * n = static_cast<DelaunayTriangle *>(*i)->getNeighbourhood(j) ;
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
		
		static_cast<DelaunayTriangle *>(*i)->setBoundingPoints(newPoints) ;
	}
			
	
	for(auto i = tree.begin() ; i != tree.end() ; ++i)
	{
		(*i)->clearVisited() ;
	}

}
