//
// C++ Interface: isotropiclineardamage
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2008-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef MU_ROTATING_CRACK_H
#define MU_ROTATING_CRACK_H

#include "damagemodel.h"

extern Amie::Matrix E;
namespace Amie {

/** \brief Rotating crack damage model. The stifness of an affected element is scaled by a factor between 1 and 0
	@author Cyrille Dunant <cyrille.dunant@epfl.ch>
*/
struct OrthotropicStiffness ;
class RotatingCrack final: public DamageModel
{
protected:

    double E ;
    double nu ;
    double factor ;

// 	std::vector< std::pair<double, double> > compressionAngles ;
// 	std::vector< std::pair<double, double> > tensionAngles ;
//   std::vector<double> compressionweights ;
// 	std::vector<double> tensionweights ;
    bool firstTension ;
    bool secondTension ;
    bool firstMet ;
    bool secondMet ;
    bool firstTensionFailure ;
    bool secondTensionFailure ;
    bool firstCompressionFailure ;
    bool secondCompressionFailure ;
    bool broken ;
    double currentAngle ;
    double initialAngle ;
    

    
    bool newtonIteration ;
    std::vector<std::pair<double, double>> angles_scores ;
    ElementState * es ;

    OrthotropicStiffness * stiff ;
    bool postprocheck ;

    bool roughsampling ;
    int iterationcount ;

public:
  
      double damage0 ;
    double damage1 ;
    double damage2 ;
    double damage3 ;


    /** \brief Constructor. Set the number of degrees of freedom
     *
     * @param numDof number of degrees of freedom
     */
    RotatingCrack(double E, double nu) ;

    virtual ~RotatingCrack();
    virtual void scale(double s) {
        factor = s ;
    } ;
    virtual void computeDelta( ElementState &s) ;

    /** \brief Increment the damage
     *
     * @param s ElementState passed as a parameter
     */
    virtual std::pair<Vector, Vector> computeDamageIncrement(ElementState & s) /*override*/;

    /** \brief compute the new stifness matrix after damage
     *
     * \f$ K' = K(1-d) \f$
     * @param m Matrix to modify
     * @return the new Matrix
     */
    virtual Matrix apply(const Matrix & m, const Point & p = Point(), const IntegrableEntity * e = nullptr, int g = -1) const;

    /** \brief return true is the element concerned is fractured
    	*/
    virtual bool fractured(int direction = -1) const  ;

    virtual void postProcess() ;

    virtual void step(ElementState & s, double maxscore)  ;

    virtual DamageModel * getCopy() const ; 

    virtual int getMode() const ;

    virtual double getAngleShift() const ;
};

/** \brief fixed crack damage model. The stiffness of an affected element is scaled by a factor between 1 and 0
	@author Cyrille Dunant <cyrille.dunant@epfl.ch>
*/
class FixedCrack : public DamageModel
{
protected:

    double E ;
    double nu ;
    double factor ;

// 	std::vector< std::pair<double, double> > compressionAngles ;
// 	std::vector< std::pair<double, double> > tensionAngles ;
//   std::vector<double> compressionweights ;
// 	std::vector<double> tensionweights ;
    bool firstTension ;
    bool secondTension ;
    bool firstTensionFailure ;
    bool secondTensionFailure ;
    bool firstCompressionFailure ;
    bool secondCompressionFailure ;
    ElementState * es ;
    bool angleset ;
public:

    double currentAngle ;

    /** \brief Constructor. Set the number of degrees of freedom
     *
     * @param numDof number of degrees of freedom
     */
    FixedCrack(double E, double nu) ;

    virtual ~FixedCrack();
    virtual void scale(double s) {
        factor = s ;
    } ;
    virtual void computeDelta(ElementState &s) ;

    /** \brief Increment the damage
     *
     * @param s ElementState passed as a parameter
     */
    virtual std::pair<Vector, Vector> computeDamageIncrement(ElementState & s) /*override*/;

    /** \brief compute the new stifness matrix after damage
     *
     * \f$ K' = K(1-d) \f$
     * @param m Matrix to modify
     * @return the new Matrix
     */

    virtual Matrix apply(const Matrix & m, const Point & p = Point(), const IntegrableEntity * e = nullptr, int g = -1) const;
    /** \brief return true is the element concerned is fractured
    	*/
    virtual bool fractured(int direction = -1) const  ;

    virtual void postProcess() ;

    virtual DamageModel * getCopy() const ; 

    virtual int getMode() const ;

    virtual double getAngleShift() const ;
};

}

#endif
