//
// C++ Interface: fracturecriterion
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2018-
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef MUCOLLISIONDETECTOR_H
#define MUCOLLISIONDETECTOR_H

#include "../../utilities/matrixops.h"
#include "../../elements/integrable_entity.h"
#include "../../mesher/mesh.h"
#include "../../features/features.h"
#include "../fracturecriteria/fracturecriterion.h"

namespace Amie {

class DelaunayTriangle ;
class DelaunayTetrahedron ;

// typedef enum {
//     NULL_SMOOTH,
//     MAX_PROXIMITY_SMOOTH,
//     GAUSSIAN_SMOOTH
// } NonLocalSmoothingType ;



/**
Abstract definition of a fracture criterion

	@author Cyrille Dunant <cyrille.dunant@gmail.com>
*/
class CollisionDetector : virtual public FractureCriterion
{
    
public:

    CollisionDetector() ;

    virtual ~CollisionDetector();

    virtual void step(ElementState& s) ;

    virtual std::pair<double, double> setChange( Amie::ElementState& s, double thresholdScore )  ;
    
    virtual double getMaxScoreInNeighbourhood(ElementState & s) ;
};

} 

#endif
