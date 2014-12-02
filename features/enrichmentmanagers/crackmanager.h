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
    protected:
        std::vector<BranchedCrack *> movingSet ;
        std::vector<std::pair<Point *, Point *>> movingTips ;
        double criticalEnergy;
    public:
            CrackManager(BranchedCrack * first, double criticalEnergy = 1e6) : EnrichmentManager(first), criticalEnergy(criticalEnergy) {} ;
              
            virtual bool step(double dt, Vector * v, Mesh< DelaunayTriangle, DelaunayTreeItem >* dtree) 
            { 
                if(stable) //checkpoint: start iteration
                {
                    stable = false ;
                    
                    std::vector<double> scores ;
                    double maxScore
                }
                
                bool moved = false ;
                for(auto i : featureSet)
                {
                    i->step(dt, v, dtree) ;
                    moved = moved || i->moved() ;
                }
                return moved ;
                
            };

            virtual bool step(double dt, Vector * v, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * dtree) 
            { 

            };

            bool converged() { return stable ; } ;
    } ;
    
} ;


#endif // CRACK_MANAGER