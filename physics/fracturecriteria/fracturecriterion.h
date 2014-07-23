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

#include "../../elements/integrable_entity.h"
#include "../../mesher/mesh.h"

namespace Mu {

class DelaunayTriangle ;
class DelaunayTetrahedron ;

typedef enum{
  NO_MIRROR,
  MIRROR_X,
  MIRROR_Y,
  MIRROR_Z,
  MIRROR_XY,
  MIRROR_XZ,
  MIRROR_YZ
} MirrorState ;

typedef enum{
	NULL_SMOOTH,
	MAX_PROXIMITY_SMOOTH,
	GAUSSIAN_SMOOTH
} NonLocalSmoothingType ;


typedef enum{
	FROM_STRESS_STRAIN,
	FROM_PRINCIPAL_STRESS_STRAIN
} SmoothingSourceType ;

typedef enum{
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

		std::valarray<unsigned int> physicalcache ;
		void initialiseFactors(const Mu::ElementState& s) ;
		std::valarray<double> factors ;
		
		std::vector<unsigned int> damagingSet ;
		std::vector<unsigned int> proximitySet ;

		double initialScore ;
		double neighbourhoodvolume ;
		double physicalCharacteristicRadius ;
		double scoreAtState ;
		double nonLocalScoreAtState ;
		double deltaScoreAtState ;
		double deltaEnergyAtState ;
		double energyDamageDifferential ;
		double criterionDamageDifferential ;
		double currentEnergy ;
		double previousEnergy ;
		MirrorState mirroring ;
		double delta_x ;
		double delta_y ;
		double delta_z ;
		
		bool energyIndexed ;
		bool noEnergyUpdate ;
		bool metAtStep ;
		bool stable ;
		
		double currentAngle ;
		double minDeltaInNeighbourhood ;
		int maxModeInNeighbourhood ;
		double maxScoreInNeighbourhood ;
		double maxAngleShiftInNeighbourhood ;
		
		double scoreTolerance ;
		bool checkpoint ;
		bool inset ;
		SmoothingFunctionType smoothingType ;
		
		double getDeltaEnergy(const ElementState & s, double delta_d) ;
		
		// as described in C. Giry et al. / International Journal of Solids and Structures 48 (2011) 3431–3443
		double getSquareInfluenceRatio( ElementState & s, const Point & direction) ;
		
		double cachedInfluenceRatio ;
		
	public:
		bool inIteration ;
		
		virtual bool directionInTension(size_t direction) {return true ;}
		virtual bool directionInCompression(size_t direction) {return true ;}
		virtual bool directionMet(size_t direction) {return metAtStep ;}
		
		Vector smoothedPrincipalStress( ElementState &s, StressCalculationMethod m = REAL_STRESS) ;
		double smoothedScore(ElementState& s) ;
		Vector smoothedPrincipalStrain( ElementState &s) ;
		std::pair<Vector, Vector> smoothedStressAndStrain( ElementState &s , StressCalculationMethod m = REAL_STRESS, double t = 0) ;
		std::pair<Vector, Vector> smoothedPrincipalStressAndStrain( ElementState &s, SmoothingSourceType ss = FROM_STRESS_STRAIN , StressCalculationMethod m = REAL_STRESS,  double t = 0) ;
		double smoothedPrincipalStressAngle( ElementState &s, StressCalculationMethod m = REAL_STRESS) ;
		double smoothedCrackAngle( ElementState &s) const ;
		double getCurrentAngle() const {return currentAngle ; }
		
		std::pair<double, double> getCrackOpeningAndSlip(const ElementState & s) ;
		
	public:
		std::vector<unsigned int> cache ;
		Mesh<DelaunayTriangle, DelaunayTreeItem>  *mesh2d ;
		Mesh<DelaunayTetrahedron,DelaunayTreeItem3D>  *mesh3d ;
		
		double getEnergyDamageDifferential()  const {return energyDamageDifferential ;}
		double getCriterionDamageDifferential()  const {return criterionDamageDifferential ;}
		double getDeltaEnergyAtState() const {return deltaEnergyAtState ;}
		double getScoreTolerance() const { return scoreTolerance ; }
		double getMinDeltaInNeighbourhood() const { return minDeltaInNeighbourhood ;}
		double getMaxScoreInNeighbourhood() const ;
		
		FractureCriterion(MirrorState mirroring = NO_MIRROR, double delta_x = 0, double delta_y = 0, double delta_z = 0) ;
		virtual void initialiseCache(const Mu::ElementState& s)
 ;
	
		virtual ~FractureCriterion();
		
		void step(Mu::ElementState& s) ;
		void computeNonLocalState(ElementState &s, NonLocalSmoothingType st = MAX_PROXIMITY_SMOOTH) ;
		bool isAtCheckpoint() const {return checkpoint ;}
		bool isInDamagingSet() const {return inset ;}
		void setCheckpoint( bool c) {checkpoint = c ;} 
		
		void setScoreTolerance(double f) { scoreTolerance = f ;} ;
		
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
		
		virtual std::pair<double, double> setChange(const ElementState &s, double maxscore)  ;
		
		/** \brief Produce a copy of the fracture criterion
		 * 
		 * @return a new FractureCriterion
		 */
		virtual FractureCriterion * getCopy() const = 0;

		/** \brief set the neighbourhood in which to consider other elements to check for failure.
		 * 
		 * @param r new radius
		 */
		double getNeighbourhoodVolume() const { return neighbourhoodvolume ;} ;

		virtual void setMaterialCharacteristicRadius(double r) ;
		double getMaterialCharacteristicRadius() const { return physicalCharacteristicRadius ;} ;
		
		virtual Material toMaterial() ;
		
		const std::vector<unsigned int> & getCache() const { return cache ; } ;
		
		double getScoreAtState() const { return scoreAtState ; }
		double getNonLocalScoreAtState() const { return nonLocalScoreAtState ;}
		
		std::pair<double, double> getDeltaEnergyDeltaCriterion(const ElementState & s, double delta_d) const ;
		
		/** \brief Set true to compute the energy state at each time step.
		 * Set true to compute the energy state at each time step. This operation is expensive
		 * so the default is set to false. 
		 */
		void setEnergyIndexed(bool t) {energyIndexed = t ;};
		
		const std::valarray<double> & getFactors() const {return factors ;}
		
		virtual double getTensileLimit(const ElementState & s) const = 0 ;
	};

} ;

#endif
