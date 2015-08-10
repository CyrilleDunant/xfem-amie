#ifndef INCLUSION_FAMILY_H
#define INCLUSION_FAMILY_H

#include "granulo.h"
#include "../features/microstructuregenerator.h"

namespace Amie
{

/*PARSE . InclusionFamily 
    @value[number] // number of inclusions in the family
    @value[radius_maximum] // maximum radius in the family
    @value[surface] // surface occupied by the family
    @object<ParticleSizeDistribution>[particle_size_distribution]  // particle size distribution of the family
    @object<InclusionGenerator>[geometry] // rule to generate the shape and orientation of the inclusions
*/
struct InclusionFamily
{
	bool keepNoFatherFeatures ;
	std::vector<std::vector<Feature *> > features ;
	std::vector<double> rotations ;
	size_t fatherIndex ;
	InclusionFamily * father = nullptr ;

	InclusionFamily() { keepNoFatherFeatures = true ; } 
	InclusionFamily( size_t n, double rmax, double surface, ParticleSizeDistribution * type = nullptr, InclusionGenerator * geometry = nullptr) ;

	void setFather( InclusionFamily * f, size_t i ) { father = f ; fatherIndex = i ; }

	std::vector<Feature *> getFeatures(size_t i) ;

	std::vector<Geometry *> getFeaturesAsGeometry(size_t i) ;

	virtual void setBehaviour( Form * behaviour, int i = -1, bool copy = false) ;

	virtual void concatenate( InclusionFamily * brother ) ;

	virtual void place( Rectangle * box, double spacing, size_t tries, size_t seed ) ;

	virtual void addToFeatureTree( FeatureTree * f) ;
} ;

/*PARSE Embedded InclusionFamily 
    @value[number] // number of inclusions in the family
    @value[radius_maximum] // maximum radius in the family
    @value[surface] // surface occupied by the family
    @object<ParticleSizeDistribution>[particle_size_distribution]  // particle size distribution of the family
    @object<InclusionGenerator>[geometry] // rule to generate the shape and orientation of the inclusions
*/
struct EmbeddedInclusionFamily : public InclusionFamily
{
	EmbeddedInclusionFamily( size_t n, double rmax, double surface, ParticleSizeDistribution * type = nullptr, InclusionGenerator * geometry = nullptr ) ;

	virtual void place( Rectangle * box, double spacing, size_t tries, size_t seed ) ;
} ;

/*PARSE Masked InclusionFamily 
    @value[number] // number of inclusions in the family
    @value[radius_maximum] // maximum radius in the family
    @value[surface] // surface occupied by the family
    @object<ParticleSizeDistribution>[particle_size_distribution]  // particle size distribution of the family
    @object<InclusionGenerator>[geometry] // rule to generate the shape and orientation of the inclusions
*/
struct MaskedInclusionFamily : public InclusionFamily
{
	MaskedInclusionFamily( size_t n, double rmax, double surface, ParticleSizeDistribution * type = nullptr, InclusionGenerator * geometry = nullptr ) ;

	virtual void place( Rectangle * box, double spacing, size_t tries, size_t seed ) ;
} ;

/*PARSE Concentric InclusionFamily 
    @value[layer_width] // width of the layer between the different inclusions
*/
struct ConcentricInclusionFamily : public InclusionFamily
{
	size_t width ;

	ConcentricInclusionFamily( double width ) ;

	virtual void place( Rectangle * box, double spacing, size_t tries, size_t seed ) ;
} ;

/*PARSE Voronoi InclusionFamily 
    @value[radius] // target radius of the current phase
    @value[surface_fraction] // target fraction of the current phase
    @value[correction_factor] // correction_factor of the current phase
    @value[number_of_grains] // number of grains in the Voronoi diagram
    @value[maximum_vertex] 20 // maximum number of vertexes of the polygons
*/
struct VoronoiInclusionFamily : public InclusionFamily
{
	std::vector<VoronoiGrain> grains ;
	size_t n ;
	size_t nmax ;
	bool copy ;
	bool generated ;

	VoronoiInclusionFamily( double radius, double fraction, double correction, size_t n, size_t maxvertex = 20) ;

	virtual void concatenate( InclusionFamily * brothers ) ;

	virtual void place( Rectangle * box, double spacing, size_t tries, size_t seed ) ;

	virtual void setBehaviour( Form * behaviour, int i = 0, bool c = false) { grains[i].behaviour = behaviour ; copy = c ; } 

} ;



}

#endif // INCLUSION_FAMILY

