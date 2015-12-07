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

#ifndef GRANULO_H
#define GRANULO_H

#include <iostream>
#include <vector>
#include "../features/microstructuregenerator.h"
#include "../features/inclusion.h"
#include "../features/expansiveZone.h"
#include "../features/growingExpansiveZone.h"
#include "../features/timeDependentEnrichmentInclusion.h"
#include "../features/inclusion3d.h"
#include "../physics/stiffness_with_imposed_deformation.h"
#include "../geometry/geometry_base.h"


namespace Amie
{

/**
 * \brief Type of particle size distribution. If you add a new PSDType, please complete ParticleSizeDistribution::getPSD() static method and create the appropriate class.
 */
// typedef enum
// {
// 	BOLOME_A,
// 	BOLOME_B,
// 	BOLOME_C,
// 	BOLOME_D,
// 	PSD_UNIFORM,
// } PSDType ;

typedef enum
{
	CIRCLE_INCLUSION = 0,
	SPHERE_INCLUSION,
	ELLIPSE_INCLUSION,
} TypeInclusion ;

typedef enum { 
	CUMULATIVE_PERCENT,
	CUMULATIVE_FRACTION,
	CUMULATIVE_ABSOLUTE,
	CUMULATIVE_PERCENT_REVERSE,
	CUMULATIVE_FRACTION_REVERSE,
	CUMULATIVE_ABSOLUTE_REVERSE
} PSDSpecificationType;



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

class ParticleSizeDistribution ;

struct VoronoiGrain
{
    Form * behaviour ;
    double radius ;
    double fraction ;
    double correctionFactor ;

    VoronoiGrain(Form * b, double r, double f, double c = 1.) : behaviour(b), radius(r), fraction(f), correctionFactor(c) { } ;

} ;

/**
 * \brief Basic class for generation of particle size distribution
 */
class PSDGenerator
{
public:
	PSDGenerator() ;

	/**
	 * \brief Generic method to get 2D inclusions
	 * @param rmax maximum radius in the particle size distribution
	 * @param surface target surface covered by the aggregates
	 * @param type type of PSD to use
	 * @param crit criteria to stop the generation
	 * @return vector of Inclusion*
	 */
	static std::vector<Inclusion *> get2DInclusions(double rmax, double surface, ParticleSizeDistribution * PSDType, PSDEndCriteria crit) ;

	/**
	 * \brief Generic method to get 3D inclusions
	 * @param rmax maximum radius in the particle size distribution
	 * @param area target area filled by the aggregates
	 * @param type type of PSD to use
	 * @param crit criteria to stop the generation
	 * @return vector of Inclusion3D*
	 */
	static std::vector<Inclusion3D *> get3DInclusions(double rmax, double area, ParticleSizeDistribution * type, PSDEndCriteria crit) ;
		
        static std::vector<std::vector<PolygonalSample *> > get2DSourceVoronoiPolygons(Rectangle * box, std::vector<VoronoiGrain> & grains, size_t n, double minDist, double delta = 0, size_t seed = 1) ;

        static std::vector<std::vector<Feature *> > get2DVoronoiPolygons(FeatureTree * F, std::vector<VoronoiGrain> & grains, size_t n, double minDist, double border = 0, size_t nmax = 16, bool copy = false, double delta = 0, size_t seed = 1) ;

        static std::vector<std::vector<Feature *> > get2DVoronoiPolygons(Rectangle * box, std::vector<VoronoiGrain> & grains, size_t n, double minDist, double border = 0, size_t nmax = 16, bool copy = false, double delta = 0, size_t seed = 1) ;

        static std::vector<std::vector<Feature *> > get2DVoronoiPolygons(Feature * feat, std::vector<VoronoiGrain> & grains, size_t n, double minDist, double border = 0, size_t nmax = 16, bool copy = false, double delta = 0, size_t seed = 1) ;

        static std::vector<std::vector<Feature *> > get2DVoronoiPolygons(FeatureTree * F, std::vector<VoronoiGrain> & grains, std::vector<Feature *> feats, size_t n, double minDist, double border = 0, size_t nmax = 16, bool copy = false, double delta = 0, size_t seed = 1) ;

        static std::vector<std::vector<Feature *> > get2DVoronoiPolygons(Rectangle * box, std::vector<VoronoiGrain> & grains, std::vector<Feature *> feats, size_t n, double minDist, double border = 0, size_t nmax = 16, bool copy = false, double delta = 0, size_t seed = 1) ;

	/**
	 * \brief Creates PSD for mortar square samples
	 * @param width width of the sample
	 * @param rmax maximum radius in the particle size distribution
	 * @param n maximum number of aggregates
	 * @param type type of PSD to use
	 * @return vector of Inclusion*
	 */
	static std::vector<Inclusion *> get2DMortar(double rmax = 0.0025, double width = 0.04, size_t n = 4000, ParticleSizeDistribution * type = nullptr) ; //D
	
	/**
	 * \brief Creates PSD for concrete square samples
	 * @param width width of the sample
	 * @param rmax maximum radius in the particle size distribution
	 * @param n maximum number of aggregates
	 * @param type type of PSD to use
	 * @return vector of Inclusion*
	 */
	static std::vector<Inclusion *> get2DConcrete(double rmax = 0.008, double width = 0.07, size_t n = 6000, ParticleSizeDistribution * type = nullptr, double percent = 0.8) ; //A

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
	static std::vector<Feature *> get2DConcrete(FeatureTree * F, Form * behaviour,  size_t n = 6000, double rmax = 0.008, double itz = 0, ParticleSizeDistribution * type = nullptr, InclusionGenerator * geometry = nullptr, size_t tries = 100000, double percent=0.8, Geometry * placement = nullptr,  std::vector<Geometry *> exclusionZones = std::vector<Geometry *>(), size_t seed = 0) ;

	static std::vector<Feature *> get2DEmbeddedInclusions(FeatureTree * F, Form * behaviour,  std::vector<Feature *> base, size_t n = 6000, double rmax = 0.008, double itz = 0, ParticleSizeDistribution * type = nullptr, InclusionGenerator * geometry = nullptr, size_t tries = 100000, double percent=0.1, Geometry * placement = nullptr,  std::vector<Geometry *> exclusionZones = std::vector<Geometry *>(), size_t seed = 0) ;

	static std::vector<Feature *> get2DMaskedInclusions(FeatureTree * F, Form * behaviour,  std::vector<Feature *> mask, size_t n = 6000, double rmax = 0.008, double itz = 0, ParticleSizeDistribution * type = nullptr, InclusionGenerator * geometry = nullptr, size_t tries = 100000, double percent=0.1, Geometry * placement = nullptr,  std::vector<Geometry *> exclusionZones = std::vector<Geometry *>(), size_t seed = 0) ;

	static std::vector<Feature *> get2DInclusionsOnEdge(FeatureTree * F, Form * behaviour,  std::vector<Feature *> base, bool checkMask = true, bool onVertex = false, size_t n = 6000, double rmax = 0.008, double itz = 0, ParticleSizeDistribution * type = nullptr, InclusionGenerator * geometry = nullptr, size_t tries = 100000, double percent=0.1, Geometry * placement = nullptr,  std::vector<Geometry *> exclusionZones = std::vector<Geometry *>(), size_t seed = 0) ;

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
	static std::vector<Inclusion *> get2DMortar(FeatureTree * F, Form * behaviour, double rmax = 0.0025, size_t n = 4000, ParticleSizeDistribution * type = nullptr, size_t tries = 10000, size_t seed = 0) ; //D

	/** \brief Creates expansive zones randomly distribued in the aggregates
	 * @param F feature tree
	 * @param aggregates reactive aggregates
	 * @param behaviour expansive zones behaviour
	 * @param radius radius of the expansive zones
	 * @param n number of generated zones
	 * @param max maximum number of kept zones
	 */
	static std::vector<std::pair<ExpansiveZone *, Inclusion *> > get2DExpansiveZonesInAggregates(FeatureTree * F, std::vector<Inclusion *> aggregates, StiffnessWithImposedDeformation * behaviour, double radius, size_t n, size_t max, int maxPerAgg = -1) ;

	static std::vector<std::pair<TimeDependentHomogenisingInclusion *, Inclusion *> > get2DGrowingExpansiveZonesInAggregates(FeatureTree * F, std::vector<Inclusion *> aggregates, ViscoelasticityAndImposedDeformation * behaviour, Function radius, double rmax, size_t n, size_t max, int maxPerAgg = -1) ;

} ;

class ParticleSizeDistribution
{
public:
	virtual ~ParticleSizeDistribution() {} ;
		/**
	 * \brief Gives the next diameter for 2D inclusions. This method is to be overloaded by inherited classes
	 * @param diameter previous diameter found in the distribution
	 * @param fraction remaining fraction of aggregates to place
	 * @param dmax maximum diameter of the distribution
	 * @return next diameter
	 */
		virtual double getNext3DDiameter(double diameter, double fraction, double dmax) = 0;
		
			/**
	 * \brief Gives the next diameter for 3D inclusions. This method is to be overloaded by inherited classes
	 * @param diameter previous diameter found in the distribution
	 * @param fraction remaining fraction of aggregates to place
	 * @param dmax maximum diameter of the distribution
	 * @return next diameter
	 */
		virtual double getNext2DDiameter(double diameter, double fraction, double dmax) = 0;

} ;

/*PARSE PSDBolomeA ParticleSizeDistribution */
class PSDBolomeA : public ParticleSizeDistribution
{
public:
	virtual double getNext2DDiameter(double diameter, double fraction, double dmax) ;
	virtual double getNext3DDiameter(double diameter, double fraction, double dmax) ;
} ;

/*PARSE PSDBolomeB ParticleSizeDistribution */
class PSDBolomeB : public ParticleSizeDistribution
{
public:
	virtual double getNext2DDiameter(double diameter, double fraction, double dmax) ;
	virtual double getNext3DDiameter(double diameter, double fraction, double dmax) ;
} ;

/*PARSE PSDBolomeC ParticleSizeDistribution */
class PSDBolomeC : public ParticleSizeDistribution
{
public:
	virtual double getNext2DDiameter(double diameter, double fraction, double dmax) ;
	virtual double getNext3DDiameter(double diameter, double fraction, double dmax) ;
} ;

/*PARSE PSDBolomeD ParticleSizeDistribution */
class PSDBolomeD : public ParticleSizeDistribution
{
public:
	virtual double getNext2DDiameter(double diameter, double fraction, double dmax) ;
	virtual double getNext3DDiameter(double diameter, double fraction, double dmax) ;
} ;

/*PARSE PSDFuller ParticleSizeDistribution 
	@value[radius_minimum] 0 // minimum radius of the Fuller distribution
	@value[exponent] 0.5 // exponent of the Fuller distribution
*/
class PSDFuller : public ParticleSizeDistribution
{
public:
	double dmin ;
	double exponent ;

	PSDFuller(double min = 0., double exp = 0.5) : ParticleSizeDistribution(), dmin(min), exponent(exp) { } ;
	virtual ~PSDFuller() {} ;

public:
	virtual double getNext2DDiameter(double diameter, double fraction, double dmax) ;
	virtual double getNext3DDiameter(double diameter, double fraction, double dmax) ;
} ;


/*PARSE ConstantSizeDistribution ParticleSizeDistribution */
class ConstantSizeDistribution : public ParticleSizeDistribution
{
public:
	virtual double getNext2DDiameter(double diameter, double fraction, double dmax) 
	{
		return diameter ;
	}
	virtual double getNext3DDiameter(double diameter, double fraction, double dmax) 
	{
		return diameter ;
	}
} ;

/*PARSE GranuloFromCumulativePSD ParticleSizeDistribution 
	@string[file_name]  // path to the file containing the particle size distribution
	@string<PSDSpecificationType>[specification] // how to read the file
	@value[factor] 1 // coefficient by which all radii will be multiplied
	@value[radius_maximum] -1 // cuts off the distribution above the specified radius (if positive)
	@value[radius_minimum] -1 // cuts off the distribution below the specified radius (if positive)
*/
class GranuloFromCumulativePSD : public ParticleSizeDistribution
{
private:
	std::vector<double> fraction ;
	std::vector<double> radius ;
public:
	GranuloFromCumulativePSD(const std::string & filename, PSDSpecificationType t, double factor = 1., double cutOffUp = -1, double cutOffDown = -1) ;
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
	GranuloFromFile(const std::string & fname, std::vector<std::string> columns) ;
	int numberOfField() {return this->fields.size() ; } ;
	bool verifyField(std::vector<std::string> columns) ;
	int getFieldNumber(std::string column) ;
	std::vector<double> getFieldValues(std::string column) ;
	std::vector<double> getFieldValues(int) ;
	std::vector<Feature *> getFeatures(TypeInclusion type, int ninc) ;
	std::vector<Inclusion3D *> getInclusion3D(int ninc, double scale = 1) ;


} ;




}

#endif
