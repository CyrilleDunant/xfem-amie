//
// C++ Interface: mohrcoulomb
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef MCFT_H
#define MCFT_H

#include "fracturecriterion.h"

namespace Amie {

typedef enum
{
    UPPER_BOUND,
    LOWER_BOUND,
    AVERAGE
} RedistributionType ;

/*PARSE NonLocalMCFT FractureCriterion
    @value[compressive_strength] // compressive strength of the material
    @value[young_modulus] // Young modulus of the material
    @value[material_characteristic_radius] // characteristic radius of the non-local damage model
    @string<RedistributionType>[redistribution_type] // indicates how the strength is redistributed towards the rebars
*/
class NonLocalMCFT : public FractureCriterion
{
protected:
    double getBareConcreteTensileCriterion(const ElementState & s, double pseudoYoung, double tstrain, double tstress) ;
    double getRebarConcreteTensileCriterion(const ElementState & s, double pseudoYoung, double tstrain, double tstress, double value = 342, double deltaCriterion = 0) ; //Bentz model 342
    double getConcreteTensileCriterion(const ElementState & s, double pseudoYoung, double tstrain, double tstress) ;
    double getConcreteCompressiveCriterion(const ElementState & s, double pseudoYoung, double cstrain, double tstress, double cstress) ;
    void initialise(ElementState &s) ;
    RedistributionType rtype ;
    bool inRebarInfluence  ;
    double distanceToRebar  ;
    double effectiveInfluenceDistance ;
    bool crackInitiated ;
public:
    bool initialised ;
    bool radiusInitialised ;
    bool firstTension ;
    bool firstCompression ;
    bool secondTension ;
    bool secondCompression ;
    bool firstMet ;
    bool secondMet ;
    double upVal ;
    double downVal ;
    double tensionCritStrain ;
    double critStrain ;
    double youngModulus ;
    double k ;
    double strain_ch ;
    double strain_te ;


    virtual bool directionInTension(size_t direction, double t = 0)
    {
        if(direction == 0)
            return firstTension ;
        if(direction == 1)
            return secondTension ;
        return false ;
    }
    virtual bool directionInCompression(size_t direction, double t = 0)
    {
        if(direction == 0)
            return firstCompression ;
        if(direction == 1)
            return secondCompression ;
        return false ;
    }
    virtual bool directionMet(size_t direction, double t = 0)
    {
        if(direction == 0)
            return firstMet ;
        if(direction == 1)
            return secondMet ;
        return false ;
    }
    std::vector<std::pair<double, double> > rebarLocationsAndDiameters ;
    /** \brief Constructor, set the maximum and minimum strain
     * @param up Maximum stress (tension)
     * @param down Minimum stress (compression)
    */
    NonLocalMCFT(double down, double youngModulus, double charDistance, RedistributionType r = UPPER_BOUND) ;

    virtual ~NonLocalMCFT();

    /** \brief Return a copy of this fracture criterion*/
    virtual FractureCriterion * getCopy() const;

    /** \brief return the normalised distance to the fracture surface
     *
     * The distance is computed as: \f$ 1.-|\frac{max\; principal\; strain\; in\; element}{Limit\; strain}|  \f$
     * @param s ElementState to consider
    */
    virtual double grade(ElementState &s)  ;
    virtual double gradeAtTime(ElementState &s, double t) ;

    virtual double getTensileLimit(const ElementState & s) const ;
};

/*PARSE NonLocalSpaceTimeMCFT FractureCriterion
    @value[compressive_strength] // compressive strength of the material
    @value[young_modulus] // Young modulus of the material
    @value[material_characteristic_radius] // characteristic radius of the non-local damage model
    @string<RedistributionType>[redistribution_type] // indicates how the strength is redistributed towards the rebars
*/
class NonLocalSpaceTimeMCFT : public NonLocalMCFT
{
public:
    NonLocalSpaceTimeMCFT(double down, double youngModulus, double charDistance, RedistributionType r = UPPER_BOUND) ;
    virtual ~NonLocalSpaceTimeMCFT();

    virtual double grade(ElementState &s)  ;

    virtual FractureCriterion * getCopy() const;
} ;

}

#endif
