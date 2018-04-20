// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2014
// Copyright: See COPYING file that comes with this distribution
//
//

#include "../features.h"
#include "../enrichmentInclusion.h"
#include "../expansiveZone.h"
#include "../../utilities/inclusion_family.h"

#ifndef GEL_MANAGER
#define GEL_MANAGER
namespace Amie
{
    
/*PARSE Gel EnrichmentManager
        @object<FeatureTree>[feature_tree] // main feature tree of the simulation
        @object<InclusionFamily>[zones] // contains both the gel pockets already placed and already associated with a father inclusion
        @value[reactive_fraction] 0.1 // fraction of the aggregate covered at the end of the reaction
        @value[radius_increment] -1 // delta_r applied at each time step of the simulation (warning: independent of actual time step!)
    */
    class GelManager : public EnrichmentManager
    {
    protected:
        double steps = 0 ;
        double totalSteps = 100 ;
        std::vector<std::pair<ExpansiveZone *, Feature *> > zones ;
        double deltaRadius ;
        double reactedArea ;
        double aggregateArea ;
        double reactiveFraction ;
        FeatureTree * ftree ;
    public:
            GelManager(FeatureTree * ftree, double zonedensity, const std::vector<Feature *> & aggregates, double reactiveFraction = 0.03) ;
            GelManager( FeatureTree * f, InclusionFamily * zones, double reactiveFraction = 0.03) ;
            GelManager(FeatureTree * f) ;
              
            virtual bool step(double dt, Vector * v, Mesh< DelaunayTriangle, DelaunayTreeItem >* dtree) ;
            virtual bool step(double dt, Vector * v, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * dtree) ;
            virtual void setDeltaRadius(double dr) ;
            virtual double getReactedFraction() const ;
            virtual double getAggregateArea() const ;
            virtual double getReactiveFraction() const ;
            virtual void setReactiveFraction(double r) ;
            virtual bool converged() {return true ;} ;
    } ;


/*PARSE FunctionBasedGel EnrichmentManager
        @object<FeatureTree>[feature_tree] // main feature tree of the simulation
        @object<InclusionFamily>[zones] // contains both the gel pockets already placed and already associated with a father inclusion
        @string[radius] // radius of the zones as a function of time
        @value[reactive_fraction] 0.1 // fraction of the aggregate covered at the end of the reaction
    */
    class FunctionBasedGelManager : public GelManager
    {
    protected:
        Function radius ;

    public:
        FunctionBasedGelManager( FeatureTree * f, InclusionFamily * zones, std::string function, double reactiveFraction = 0.1) ;

        virtual bool step(double dt, Vector * v, Mesh< DelaunayTriangle, DelaunayTreeItem >* dtree) ;

    } ;

/*PARSE SpaceTimeGel EnrichmentManager
        @object<FeatureTree>[feature_tree] // main feature tree of the simulation
        @object<InclusionFamily>[zones] // contains both the gel pockets already placed and already associated with a father inclusion
        @string[radius] // radius of the zones as a function of time
        @value[reactive_fraction] 0.1 // fraction of the aggregate covered at the end of the reaction
    */
    class SpaceTimeGelManager : public GelManager
    {
        std::vector<std::pair<GrowingExpansiveZone *, Feature *> > stzones ;

    public:
        SpaceTimeGelManager( FeatureTree * f, InclusionFamily * zones, std::string function, double reactiveFraction = 0.1)  ; 

        virtual bool step(double dt, Vector * v, Mesh< DelaunayTriangle, DelaunayTreeItem >* dtree) ;

    } ;


    
} 


#endif // GEL_MANAGER
