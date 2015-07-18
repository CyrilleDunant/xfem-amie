// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2014
// Copyright: See COPYING file that comes with this distribution
//
//

#include "../features.h"
#include "../crack.h"

#ifndef CRACK_MANAGER
#define CRACK_MANAGER
namespace Amie
{
    
    class CrackManager : public EnrichmentManager
    {
        friend class BranchedCrack ;
    protected:
        std::vector<std::pair<std::pair<Point *, double>, std::pair<Point *, double>>> movingTips ;
        std::vector<std::pair<Point *, BranchedCrack *>> movingSet ;
        std::vector<BranchedCrack *> featureSet ;
        std::vector<double> startAngles ;
        std::vector<Point> centers ;
        double criticalEnergy;
        double maxExpansion ;
        double minCRadius ;
        bool findRadius ;
        bool findExtension ;
        double currentRadius ;
        double upRadius ;
        double downRadius ;
        double currentExtension ;
        double upExtension ;
        double downExtension ;
        int iteration ;
    public:
            CrackManager(FeatureTree * ft, Feature * father, BranchedCrack * first, double criticalEnergy = 1e6, double minCRadius = 1e-2, double maxExpansion = 1e-2) ;
              
            virtual bool step(double dt, Vector * v, Mesh< DelaunayTriangle, DelaunayTreeItem >* dtree) ;
            virtual bool step(double dt, Vector * v, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * dtree) ;
            void setCriticalEnergy(double) ;
            void setMaximumRadius(double) ;
            void setMinimumRadius(double) ;
            void setMaximumExpansion(double) ;

            bool converged()  ;
    } ;
    
} 


#endif // CRACK_MANAGER