// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010-2011

#ifndef STRUCTURED_MESH_H
#define STRUCTURED_MESH_H

#include "mesh.h"
#include "delaunay.h"
#include "../utilities/grid.h"

namespace Amie
{
	class StructuredMesh : public Mesh<DelaunayTriangle, DelaunayTreeItem>
	{
	protected:
		std::vector< Point *> points ;
		Grid grid ;
		void addSharedNodes(size_t nodes_per_side, size_t time_planes, double timestep) ;
		size_t global_counter ;
		std::vector <DelaunayTriangle *> tree ;
        std::vector< std::vector<DelaunayTriangle *> > caches ;
        std::vector<Vector> coefs ;
	public:
		StructuredMesh(double sizeX, double sizeY, int div, const Point & center ) ;
        virtual size_t size() const { return tree.size() ; } ;
		virtual ~StructuredMesh() ;
		virtual std::vector<DelaunayTriangle *> getElements();
		virtual std::vector<DelaunayTriangle *> getConflictingElements(const Point  * p) ;
		virtual std::vector<DelaunayTriangle *> getConflictingElements(const Geometry * g) ;
		virtual std::vector<Point * > & getAdditionalPoints(){ return points ;} ;
		virtual const std::vector<Point * > & getAdditionalPoints() const {return points ;} ;
		virtual void setElementOrder(Order o, double dt = 0);
		virtual void insert(Point *) ;
		virtual size_t getLastNodeId() const;
		
        virtual unsigned int generateCache(const std::vector<DelaunayTriangle *> original) 
        { 
            caches.push_back(original);
            return caches.size()-1 ;
        } ;
        
        
    virtual unsigned int generateCache(const Geometry * locus, const Geometry * source = nullptr, Function smoothing = Function("1"))
    {
        VirtualMachine vm ;
        std::vector<double> co ;
        for(auto & element : tree)
        {
            if(source && element->getBehaviour()->getSource() != source)
                continue ;
            
            if(locus->in(element->getCenter()))
            {
                       
                Function x = element->getXTransform() ;
                Function y = element->getYTransform() ;
                Function z = element->getZTransform() ;
                Function t = element->getTTransform() ;
                for(size_t i = 0 ; i < co.size() ; i++)
                {
                    double xx= vm.eval(x, element->getGaussPoints().gaussPoints[i].first) ;
                    double xy = vm.eval(y, element->getGaussPoints().gaussPoints[i].first) ;
                    double xz = vm.eval(z, element->getGaussPoints().gaussPoints[i].first) ;
                    double xt = vm.eval(t, element->getGaussPoints().gaussPoints[i].first) ;
                    co.push_back(vm.eval(smoothing, xx, xy, xz, xt));
                }
            }
        }
        
        Vector cf(co.size()) ;
        std::copy(co.begin(), co.end(), &cf[0]) ;
        coefs.push_back(cf);
        return coefs.size()-1 ;
        
    } ;
        
        virtual std::vector<DelaunayTriangle *> getCache(unsigned int cacheID) 
        {
            return caches[cacheID] ; 
        } ;
	} ;
}

#endif // STRUCTURED_MESH_H
