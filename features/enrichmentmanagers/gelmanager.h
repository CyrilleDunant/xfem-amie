// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2014
// Copyright: See COPYING file that comes with this distribution
//
//

#include "../features.h"
#include "../enrichmentInclusion.h"
#include "../expansiveZone.h"

#ifndef GEL_MANAGER
#define GEL_MANAGER
namespace Amie
{
    
    class GelManager : public EnrichmentManager
    {
    protected:
        std::vector<std::pair<ExpansiveZone *, Feature *> > zones ;
        double deltaRadius ;
        double reactedArea ;
        double aggregateArea ;
        double reactiveFraction ;
        FeatureTree * ftree ;
    public:
            GelManager(FeatureTree * ftree, double zonedensity, const std::vector<Feature *> & aggregates,double reactiveFraction = 0.1, double deltaRadius = -1, double initalRadius = 0.00025) ;
              
            virtual bool step(double dt, Vector * v, Mesh< DelaunayTriangle, DelaunayTreeItem >* dtree) ;
            virtual bool step(double dt, Vector * v, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * dtree) ;
            void setDeltaRadius(double dr) ;
            double getReactedFraction() const ;
            double getAggregateArea() const ;
            double getReactiveFraction() const ;
            void setReactiveFraction(double r) ;
            bool converged() {return true ;} ;
    } ;
    
} 


#endif // GEL_MANAGER