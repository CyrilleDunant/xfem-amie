//
// C++ Interface: granulo
//
// Description: 
//
//
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2010-2011
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include <iostream>
#include <vector>
#include "../features/inclusion.h"
#include "../features/expansiveZone.h"
#include "../features/inclusion3d.h"
#include "../physics/stiffness_with_imposed_deformation.h"
#include "../geometry/geometry_base.h"


namespace Mu
{

/**
 * \brief Type of particle size distribution. If you add a new PSDType, please complete ParticleSizeDistribution::getPSD() static method and create the appropriate class.
 */
typedef enum
{
	BOLOME_A,
	BOLOME_B,
	BOLOME_C,
	BOLOME_D
} PSDType ;

typedef enum
{
	CIRCLE_INCLUSION = 0,
	SPHERE_INCLUSION,
	ELLIPSE_INCLUSION,
} TypeInclusion ;


/**
 * \brief Criterion for ending generation of particles.
 */
struct PSDEndCriteria
{
      double rmin ; // end the PSD when a certain radius is met (all radius are taken if negative)
      double placedFraction ; // end the PSD when enough aggregates have been placed in term of volume/surface fraction (between 0 and 1, all aggregates are taken if greater than one)
      size_t nmax ; // end the PSD when enough aggregates have been placed
      
      PSDEndCriteria(double r, double f, size_t n) 
      {
	    rmin = r ;
	    placedFraction = f ;
	    nmax = n ;
      }
      
      /**
       * \brief Checks if the minimum radius, the minimum fraction or the desired number of aggregates have been reached. 
       * @param radius the current radius in the PSD
       * @param fraction the remaining fraction of aggregates
       * @param n the number of aggregates generated
       * @return true if at least one of these conditions have been reached, false otherwise
       */ 
      bool meets(double radius, double fraction, size_t n) const 
      {
	    return radius < rmin || fraction < placedFraction || n > nmax ;
      }
      
      /**
       * \brief Prints which criteria is met (if any)
       * @param radius the current radius in the PSD
       * @param fraction the remaining fraction of aggregates
       * @param n the number of aggregates generated
       */
      void print(double radius, double fraction, size_t n) const
      {
	    if(radius < rmin)
		std::cout << "minimum radius reached!" << std::endl ;
	    if(fraction < placedFraction)
		std::cout << "fraction of aggregates reached!" << std::endl ;
	    if(n > nmax)
		std::cout << "max number of aggregates reached!" << std::endl ;
      }
} ;

/**
 * \brief Basic class for generation of particle size distribution
 */
class ParticleSizeDistribution
{
public:
	ParticleSizeDistribution() ;
	
	/**
	 * \brief Generic method to get 2D inclusions
	 * @param rmax maximum radius in the particle size distribution
	 * @param surface target surface covered by the aggregates
	 * @param type type of PSD to use
	 * @param crit criteria to stop the generation
	 * @return vector of Inclusion*
	 */
	static std::vector<Inclusion *> get2DInclusions(double rmax, double surface, PSDType type, PSDEndCriteria crit) ;

	/**
	 * \brief Generic method to get 3D inclusions
	 * @param rmax maximum radius in the particle size distribution
	 * @param area target area filled by the aggregates
	 * @param type type of PSD to use
	 * @param crit criteria to stop the generation
	 * @return vector of Inclusion3D*
	 */
	static std::vector<Inclusion3D *> get3DInclusions(double rmax, double area, PSDType type, PSDEndCriteria crit) ;
	
	/**
	 * \brief Creates PSD for mortar square samples
	 * @param width width of the sample
	 * @param rmax maximum radius in the particle size distribution
	 * @param n maximum number of aggregates
	 * @param type type of PSD to use
	 * @return vector of Inclusion*
	 */
	static std::vector<Inclusion *> get2DMortar(double rmax = 0.0025, double width = 0.04, size_t n = 4000, PSDType type = BOLOME_D) ;
	
	/**
	 * \brief Creates PSD for concrete square samples
	 * @param width width of the sample
	 * @param rmax maximum radius in the particle size distribution
	 * @param n maximum number of aggregates
	 * @param type type of PSD to use
	 * @return vector of Inclusion*
	 */
	static std::vector<Inclusion *> get2DConcrete(double rmax = 0.008, double width = 0.07, size_t n = 6000, PSDType type = BOLOME_A) ;

	/**
	 * \brief Creates concrete PSD, set the behaviour and place the inclusions in the sample
	 * @param F Feature tree
	 * @param behaviour aggregate behaviour
	 * @param rmax maximum radius in the particle size distribution
	 * @param n maximum number of aggregates
	 * @param type type of PSD to use
	 * @param tries number of tries for the placement
	 * @param seed seed for random generator
	 * @return vector of Inclusion*
	 */
	static std::vector<Inclusion *> get2DConcrete(FeatureTree * F, Form * behaviour, double rmax = 0.008, size_t n = 6000, PSDType type = BOLOME_A, size_t tries = 10000, size_t seed = 0) ;

		/**
	 * \brief Creates mortar PSD, set the behaviour and place the inclusions in the sample
	 * @param F Feature tree
	 * @param behaviour aggregate behaviour
	 * @param rmax maximum radius in the particle size distribution
	 * @param n maximum number of aggregates
	 * @param type type of PSD to use
	 * @param tries number of tries for the placement
	 * @param seed seed for random generator
	 * @return vector of Inclusion*
	 */
	static std::vector<Inclusion *> get2DMortar(FeatureTree * F, Form * behaviour, double rmax = 0.0025, size_t n = 4000, PSDType type = BOLOME_D, size_t tries = 10000, size_t seed = 0) ;

	/** \brief Creates expansive zones randomly distribued in the aggregates
	 * @param F feature tree
	 * @param aggregates reactive aggregates
	 * @param behaviour expansive zones behaviour
	 * @param radius radius of the expansive zones
	 * @param n number of generated zones
	 * @param max maximum number of kept zones
	 */
	static std::vector<std::pair<ExpansiveZone *, Inclusion *> > get2DExpansiveZonesInAggregates(FeatureTree * F, std::vector<Inclusion *> aggregates, StiffnessWithImposedDeformation * behaviour, double radius, size_t n, size_t max, int maxPerAgg = -1) ;
	
	/**
	 * \brief Returns appropriate ParticleSizeDistribution* object corresponding to the given type
	 * @param type PSD type to use
	 * @return pointer to the corresponding ParticleSizeDistribution object
	 */
	static ParticleSizeDistribution * getPSD(PSDType type) ;
	
	/**
	 * \brief Gives the next diameter for 2D inclusions. This method is to be overloaded by inherited classes
	 * @param diameter previous diameter found in the distribution
	 * @param fraction remaining fraction of aggregates to place
	 * @param dmax maximum diameter of the distribution
	 * @return next diameter
	 */
	virtual double getNext2DDiameter(double diameter, double fraction, double dmax) ;

	/**
	 * \brief Gives the next diameter for 3D inclusions. This method is to be overloaded by inherited classes
	 * @param diameter previous diameter found in the distribution
	 * @param fraction remaining fraction of aggregates to place
	 * @param dmax maximum diameter of the distribution
	 * @return next diameter
	 */
	virtual double getNext3DDiameter(double diameter, double fraction, double dmax) ;
} ;

class PSDBolomeA : public ParticleSizeDistribution
{
public:
	virtual double getNext2DDiameter(double diameter, double fraction, double dmax) ;
	virtual double getNext3DDiameter(double diameter, double fraction, double dmax) ;
} ;

class PSDBolomeB : public ParticleSizeDistribution
{
public:
	virtual double getNext2DDiameter(double diameter, double fraction, double dmax) ;
	virtual double getNext3DDiameter(double diameter, double fraction, double dmax) ;
} ;

class PSDBolomeC : public ParticleSizeDistribution
{
public:
	virtual double getNext2DDiameter(double diameter, double fraction, double dmax) ;
	virtual double getNext3DDiameter(double diameter, double fraction, double dmax) ;
} ;

class PSDBolomeD : public ParticleSizeDistribution
{
public:
	virtual double getNext2DDiameter(double diameter, double fraction, double dmax) ;
	virtual double getNext3DDiameter(double diameter, double fraction, double dmax) ;
} ;



class GranuloFromFile
{
private:
	std::string filename;
	std::vector<std::string> fields;
	std::vector<double> values ;

public:
	GranuloFromFile(std::string fname, std::vector<std::string> columns) ;
	int numberOfField() {return this->fields.size() ; } ;
	bool verifyField(std::vector<std::string> columns) ;
	int getFieldNumber(std::string column) ;
	std::vector<double> getFieldValues(std::string column) ;
	std::vector<double> getFieldValues(int) ;
	std::vector<Feature *> getFeatures(TypeInclusion type, int ninc) ;
	std::vector<Inclusion3D *> getInclusion3D(int ninc, double scale = 1) ;
} ;

}

