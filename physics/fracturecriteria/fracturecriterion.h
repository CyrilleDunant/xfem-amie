//
// C++ Interface: fracturecriterion
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef MUFRACTURECRITERION_H
#define MUFRACTURECRITERION_H

#include "../../utilities/matrixops.h"
#include "../../elements/integrable_entity.h"
#include "../../mesher/mesh.h"
#include "../../features/features.h"

namespace Amie {

class DelaunayTriangle ;
class DelaunayTetrahedron ;

typedef enum {
    NO_MIRROR,
    MIRROR_X,
    MIRROR_Y,
    MIRROR_Z,
    MIRROR_XY,
    MIRROR_XZ,
    MIRROR_YZ
} MirrorState ;

// typedef enum {
//     NULL_SMOOTH,
//     MAX_PROXIMITY_SMOOTH,
//     GAUSSIAN_SMOOTH
// } NonLocalSmoothingType ;


typedef enum {
    QUARTIC_COMPACT,
    GAUSSIAN_NONCOMPACT
} SmoothingFunctionType ;

/**
Abstract definition of a fracture criterion

	@author Cyrille Dunant <cyrille.dunant@epfl.ch>
*/
class FractureCriterion
{
    friend class DamageModel ;
protected:

    std::vector<unsigned int> damagingSet ;
    std::vector<unsigned int> proximitySet ;

    double initialScore ;
    double physicalCharacteristicRadius ;
    double scoreAtState ;
    double deltaScoreAtState ;
    double scoreAtTimeStepEnd = -2 ;

    double criterionDamageDifferential ;
    MirrorState mirroring ;
    double delta_x ;
    double delta_y ;
    double delta_z ;

    bool metAtStep ;
    bool stable ;

    double minDeltaInNeighbourhood ;
    int maxModeInNeighbourhood ;
    double maxScoreInNeighbourhood ;
    double maxAngleShiftInNeighbourhood ;

    double scoreTolerance ;
    bool checkpoint ;
    bool inset ;
    SmoothingFunctionType smoothingType ;

    double cachedInfluenceRatio ;
    int cacheID ;
    int cachecoreID ;

public:
    bool inIteration ;

    virtual bool directionInTension(size_t direction) {
        return true ;
    }
    virtual bool directionInCompression(size_t direction) {
        return true ;
    }
    virtual bool directionMet(size_t direction) {
        return metAtStep ;
    }

    std::pair<Vector, Vector> getSmoothedFields(FieldType f0, FieldType f1,  ElementState &s ,double t = 0) ;
    Vector getSmoothedField( FieldType f0,ElementState &s , double t = 0) ;

public:

    Mesh<DelaunayTriangle, DelaunayTreeItem>  *mesh2d ;
    Mesh<DelaunayTetrahedron,DelaunayTreeItem3D>  *mesh3d ;

    double getCriterionDamageDifferential()  const {
        return criterionDamageDifferential ;
    }

    double getScoreTolerance() const {
        return scoreTolerance ;
    }
    double getMinDeltaInNeighbourhood() const {
        return minDeltaInNeighbourhood ;
    }
    double getMaxScoreInNeighbourhood( Amie::ElementState& s ) ;

    FractureCriterion(MirrorState mirroring = NO_MIRROR, double delta_x = 0, double delta_y = 0, double delta_z = 0) ;


    virtual void initialiseCache( Amie::ElementState& s ) ;
    virtual ~FractureCriterion();

    void step(Amie::ElementState& s) ;
    
    bool isAtCheckpoint() const {
        return checkpoint ;
    }
    bool isInDamagingSet() const {
        return inset ;
    }
    void setCheckpoint( bool c) {
        checkpoint = c ;
    }

    void setScoreTolerance(double f) {
        scoreTolerance = f ;
    } ;

    /** \brief Return true if the fracture criterion is met
     *
     * @param s ElementState ton consider
     * @return true if the fracture criterion is met
     */
    virtual bool met() const ;

    /** \brief Return a normalised distance to the fracture surface,
     *
     * The returned value lies between -1 and 1
     * @param  ElementState to consider
     * @return a value between -1 and 1
     */
    virtual double grade(ElementState & s) = 0 ;

    virtual std::pair<double, double> setChange( Amie::ElementState& s, double thresholdScore )  ;

    /** \brief Produce a copy of the fracture criterion
     *
     * @return a new FractureCriterion
     */
    virtual FractureCriterion * getCopy() const = 0;


    virtual void setMaterialCharacteristicRadius(double r) ;
    double getMaterialCharacteristicRadius() const {
        return physicalCharacteristicRadius ;
    } ;

    double getScoreAtState() const ;

    double getScoreAtTimeStepEnd() const { return scoreAtTimeStepEnd < -1 ? getScoreAtState() : scoreAtTimeStepEnd ; }

    virtual double getTensileLimit(const ElementState & s) const = 0 ;
};

} 

#endif
