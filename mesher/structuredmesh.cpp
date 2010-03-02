#include "structuredmesh.h"

using namespace Mu ;

StructuredMesh::StructuredMesh(double sizeX, double sizeY, int div, const Point & center ): grid( sizeX,  sizeY,  div, center)
{
	global_counter = 0 ;
	for(size_t i = 0 ; i < grid.pixels.size()+1 ; i++)
	{
		for(size_t j = 0 ; j < grid.pixels[i].size()+1 ; j++)
		{
			points.push_back(new Point((center.x-sizeX*.5)*i/(grid.pixels.size()+1), (center.x-sizeX*.5)*j/(grid.pixels[i].size()+1))) ;
			points.back()->id = global_counter++ ;
		}
	}
	
	for(size_t i = 0 ; i < grid.pixels.size() ; i++)
	{
		for(size_t j = 0 ; j < grid.pixels[i].size() ; j++)
		{
			triangles.push_back(new DelaunayTriangle(NULL,NULL, points[(grid.pixels[i].size()+1)*i+j], points[(grid.pixels[i].size()+1)*i+j+1], points[(grid.pixels[i].size()+1)*(i+1)+j+1], NULL)) ;
			triangles.push_back(new DelaunayTriangle(NULL,NULL, points[(grid.pixels[i].size()+1)*i+j], points[(grid.pixels[i].size()+1)*(i+1)+j], points[(grid.pixels[i].size()+1)*(i+1)+j+1], NULL)) ;
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
				tris.push_back(triangles[idx[k]*2]) ;
				tris.push_back(triangles[idx[k]*2+1]) ;
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
StructuredMesh::~StructuredMesh() 
{
	for(size_t i = 0 ; i < points.size() ; i++)
		delete points[i] ;
	
	for(size_t i = 0 ; i < triangles.size() ; i++)
		delete triangles[i] ;
}
std::vector<DelaunayTriangle *> StructuredMesh::getElements()
{
	return triangles ;
}
std::vector<DelaunayTriangle *> StructuredMesh::getConflictingElements(const Point  * p) 
{
	std::vector<DelaunayTriangle *> ret ;
	return ret ;
}
std::vector<DelaunayTriangle *> StructuredMesh::getConflictingElements(const Geometry * g) 
{
	std::vector<DelaunayTriangle *> ret ;
	return ret ;
}
void StructuredMesh::setElementOrder(Order o)
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
			addSharedNodes(0,2,2) ;
			break ;
		}
	case CONSTANT_TIME_QUADRATIC:
		{
			addSharedNodes(0,3,2) ;
			break ;
		}
	case LINEAR_TIME_LINEAR:
		{
			addSharedNodes(0,2,2) ;
			break ;
		}
	case LINEAR_TIME_QUADRATIC:
		{
			addSharedNodes(0,3,2) ;
			break ;
		}
	case QUADRATIC_TIME_LINEAR:
		{
			addSharedNodes(1,2,2) ;
			break ;
		}
	case QUADRATIC_TIME_QUADRATIC:
		{
			addSharedNodes(1,3,2) ;
			break ;
		}
	case CUBIC_TIME_LINEAR:
		{
			addSharedNodes(2,2,2) ;
			break ;
		}
	case CUBIC_TIME_QUADRATIC:
		{
			addSharedNodes(2,3,2) ;
			break ;
		}
	case QUADRIC_TIME_LINEAR:
		{
			addSharedNodes(3,2,2) ;
			break ;
		}
	case QUADRIC_TIME_QUADRATIC:
		{
			addSharedNodes(3,3,2) ;
			break ;
		}
	case QUINTIC_TIME_LINEAR:
		{
			addSharedNodes(3,2,2) ;
			break ;
		}
	case QUINTIC_TIME_QUADRATIC:
		{
			addSharedNodes(3,3,2) ;
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

	for(std::vector<DelaunayTriangle *>::iterator i = triangles.begin() ; i != triangles.end() ; ++i)
	{
		
		(*i)->visited = true ;
			
		size_t nodes_per_plane = nodes_per_side*3+3 ;
		
		std::valarray<Point *> newPoints(nodes_per_plane*time_planes) ;
		std::valarray<bool> done(false, nodes_per_plane*time_planes) ;
		
		for(size_t plane = 0 ; plane < time_planes ; plane++)
		{
			for(size_t side = 0 ; side < 3 ; side++)
			{
				Point a((*i)->getBoundingPoint(side)) ;
				Point b((*i)->getBoundingPoint((side+1)%3)) ;
				
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
					
					for(size_t j = 0 ; j< (*i)->getBoundingPoints().size() ; j++)
					{
						if((*i)->getBoundingPoint(j) == proto)
						{
							foundPoint = &(*i)->getBoundingPoint(j) ;
							break ;
						}
					}
					
					if(!foundPoint)
					{
						for(size_t j = 0 ; j < (*i)->neighbourhood.size() ; j++)
						{
							if((*i)->getNeighbourhood(j)->visited)
							{
								DelaunayTriangle * n = (*i)->getNeighbourhood(j) ;
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
		
		(*i)->setBoundingPoints(newPoints) ;
	}
			
	
	for(std::vector<DelaunayTriangle *>::iterator i = triangles.begin() ; i != triangles.end() ; ++i)
	{
		(*i)->clearVisited() ;
	}

}
