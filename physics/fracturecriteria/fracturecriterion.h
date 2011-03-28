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

	/**
	Abstract definition of a fracture criterion
	
		@author Cyrille Dunant <cyrille.dunant@epfl.ch>
	*/
	class FractureCriterion
	{
	protected:
		std::vector<unsigned int> cache ;
		std::vector<unsigned int> physicalcache ;
		std::vector<double> area ;

		double neighbourhoodradius ;
		double neighbourhoodvolume ;
		double physicalCharacteristicRadius ;
		double scoreAtState ;
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
		
	double getDeltaEnergy(const ElementState & s, double delta_d) ;
		
	public:
		std::vector<DelaunayTreeItem *>  *mesh2d ;
		std::vector<DelaunayTreeItem3D *>  *mesh3d ;
		
		double getEnergyDamageDifferential()  const {return energyDamageDifferential ;}
		double getCriterionDamageDifferential()  const {return criterionDamageDifferential ;}
		double getDeltaEnergyAtState() const {return deltaEnergyAtState ;}
		bool metInTension ;
		bool metInCompression ;
		
		FractureCriterion(MirrorState mirroring = NO_MIRROR, double delta_x = 0, double delta_y = 0, double delta_z = 0) ;
		virtual void initialiseCache(const ElementState & s) ;
	
		virtual ~FractureCriterion();
		
		void step(const Mu::ElementState& s) ;
		
		/** \brief Return true if the fracture criterion is met
		 * 
		 * @param s ElementState ton consider
		 * @return true if the fracture criterion is met
		 */
		virtual bool met(const ElementState & s) ;

		/** \brief Return a normalised distance to the fracture surface, 
		 * 
		 * The returned value lies between -1 and 1
		 * @param  ElementState to consider
		 * @return a value between -1 and 1
		 */
		virtual double grade(const ElementState & s) = 0 ;
		
		virtual int getRank(int fractiles, const ElementState &s) const ;
		
		virtual double getSteppedScore() const {return scoreAtState ;} ;
		
		/** \brief Produce a copy of the fracture criterion
		 * 
		 * @return a new FractureCriterion
		 */
		virtual FractureCriterion * getCopy() const = 0;

		/** \brief set the neighbourhood in which to consider other elements to check for failure.
		 * 
		 * @param r new radius
		 */
		virtual void setNeighbourhoodRadius(double r) ;
		double getNeighbourhoodRadius() const { return neighbourhoodradius ;} ;
		double getNeighbourhoodVolume() const { return neighbourhoodvolume ;} ;

		virtual void setMaterialCharacteristicRadius(double r) ;
		double getMaterialCharacteristicRadius() const { return physicalCharacteristicRadius ;} ;
		
		virtual Material toMaterial() ;
		
		const std::vector<unsigned int> & getCache() const { return cache ; } ;
		
		double getScoreAtState() const { return scoreAtState ; }
		
		std::pair<double, double> getDeltaEnergyDeltaCriterion(const ElementState & s, double delta_d) const ;
		
		/** \brief Set true to compute the energy state at each time step.
		 * Set true to compute the energy state at each time step. This operation is expensive
		 * so the default is set to false. 
		 */
		void setEnergyIndexed(bool t) {energyIndexed = t ;};
		
		virtual void multiply(double f) { } ;
	
	};

} ;

#endif
