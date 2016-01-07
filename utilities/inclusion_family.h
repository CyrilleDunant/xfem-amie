#ifndef INCLUSION_FAMILY_H
#define INCLUSION_FAMILY_H

#include "granulo.h"
#include "../features/microstructuregenerator.h"

namespace Amie
{

struct InclusionFamily ;
struct InclusionFamilyTree
{
	std::vector<std::vector< Feature *> > features ;
	std::vector<InclusionFamily *> brothers ;
	InclusionFamily * father ;

	InclusionFamilyTree( InclusionFamily * f = nullptr ) : father(f) { } ;

	void insert( InclusionFamily * family ) ; 
	std::vector<size_t> getAllVoronoiFamilies() ;
	void place( Rectangle * box, double spacing, size_t tries, size_t seed ) ;
	void concatenate() ;
	std::vector<Geometry *> getFeaturesAsGeometry(size_t index) ;
} ;


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
	bool placed ;
	std::vector<Feature *> features ;
	std::vector<bool> added ;
	double rotations ;
	double factors ;
	InclusionFamilyTree * brothers ;
	InclusionFamilyTree * sons ;

	InclusionFamily() ; 
	InclusionFamily( size_t n, double rmax, double surface, ParticleSizeDistribution * type = nullptr, InclusionGenerator * geometry = nullptr) ;
        InclusionFamily(Feature * matrix) ;

	InclusionFamily * getFather() ;
	std::vector<Geometry *> getFeaturesAsGeometry() ;
	std::vector<Feature *> getAddedFeatures() ;
	size_t size() const { return features.size() ; }
	virtual void setBehaviour( Form * behaviour, bool copy = false) ;
	virtual void place( Rectangle * box, double spacing, size_t tries, size_t seed, std::vector<Geometry *> & placedFeatures ) ;
	virtual void addToFeatureTree( FeatureTree * f) ;
	void setSamplingFactor(double f) { factors = f ; }
	virtual bool inNeighbours( Feature * f ) ;
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

	virtual void place( Rectangle * box, double spacing, size_t tries, size_t seed, std::vector<Geometry *> & placed ) ;
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

	virtual void place( Rectangle * box, double spacing, size_t tries, size_t seed, std::vector<Geometry *> & placed ) ;
} ;

/*PARSE Concentric InclusionFamily 
    @value[layer_width] // width of the layer between the different inclusions
*/
struct ConcentricInclusionFamily : public InclusionFamily
{
	size_t width ;

	ConcentricInclusionFamily( double width ) ;

	virtual void place( Rectangle * box, double spacing, size_t tries, size_t seed, std::vector<Geometry *> & placed ) ;
} ;

/*PARSE Voronoi InclusionFamily 
    @value[radius] // target radius of the current phase
    @value[surface_fraction] // target fraction of the current phase
    @value[correction_factor] // correction_factor of the current phase
    @value[number_of_grains] // number of grains in the Voronoi diagram
    @value[outside_layer] -1 // width of the layer outside of the box
    @value[interface] 0 // width of the interface between the grains
    @value[maximum_vertex] 20 // maximum number of vertexes of the polygons
*/
struct VoronoiInclusionFamily : public InclusionFamily
{
	VoronoiGrain grains ;
	size_t nmax ;
	double outside ;
	double interface ;
	bool copy ;

	VoronoiInclusionFamily( double radius, double fraction, double correction, size_t n, double out = -1, double width = 0, size_t maxvertex = 20) ;

	virtual void place( Rectangle * box, double spacing, size_t tries, size_t seed, std::vector<Geometry *> & placed ) ;
	virtual void setBehaviour( Form * behaviour, bool c = false) { grains.behaviour = behaviour ; copy = c ; } 
	virtual bool inNeighbours( Feature * f ) ;

} ;

/*PARSE FileDefinedCircle InclusionFamily
    @value[number] // maximum number of inclusions of the family
    @string[file_name] // file in which the circles information are stored
    @string[column_1] radius // description of the first column in the file
    @string[column_2] center_x // description of the second column in the file
    @string[column_3] center_y // description of the third column in the file
*/
struct FileDefinedCircleInclusionFamily : public InclusionFamily
{
	FileDefinedCircleInclusionFamily( int n, std::string file, std::string column_1 = "radius", std::string column_2 = "center_x", std::string column_3 = "center_y" ) ;

	virtual void place( Rectangle * box, double spacing, size_t tries, size_t seed, std::vector<Geometry *> & placed ) ;

} ;

/*PARSE FileDefinedPolygon InclusionFamily
    @value[number] // maximum number of inclusions of the family
    @string[file_name] // file in which the circles information are stored
*/
struct FileDefinedPolygonInclusionFamily : public InclusionFamily
{
	FileDefinedPolygonInclusionFamily( int n, std::string file) ;

	virtual void place( Rectangle * box, double spacing, size_t tries, size_t seed, std::vector<Geometry *> & placed ) ;

} ;



}

#endif // INCLUSION_FAMILY

