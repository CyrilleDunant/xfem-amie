#include "inclusion_family.h"
#include "placement.h"

namespace Amie
{

InclusionFamily::InclusionFamily( size_t n, double rmax, double fraction, ParticleSizeDistribution * type, InclusionGenerator * geometry) 
{
	keepNoFatherFeatures = true ;

	if(type == nullptr)
		type = new ConstantSizeDistribution() ;
	if(geometry == nullptr)
		geometry = new InclusionGenerator() ;

	std::vector<Inclusion *> incs = PSDGenerator::get2DConcrete( rmax, std::sqrt(fraction), n, type, 1. ) ;
	std::vector<Feature *> feats = geometry->convert( incs ) ;
	features.push_back(feats) ;
	rotations.push_back( geometry->authorizeRotationsDuringPlacement ) ;
}

std::vector<Feature *> InclusionFamily::getFeatures(size_t i)
{
	if(i < features.size()) { return features[i] ; }
	if(features.size() == 1) { return features[0] ; }
	std::vector<Feature *> ret ;
	for(size_t j = 0 ; j < features.size() ; j++)
	{
		for(size_t k = 0 ; k < features[j].size() ; k++)
			ret.push_back(features[j][k]) ;
	}
	return ret ;
}

void InclusionFamily::setBehaviour( Form * behaviour, int i, bool copy)
{
	size_t j = (i >= 0 ? i : features.size()-1) ;
	for(size_t k = 0 ; k < features[j].size() ; k++)
		features[j][k]->setBehaviour( copy ? behaviour->getCopy() : behaviour ) ;
}

void InclusionFamily::concatenate( InclusionFamily * brother ) 
{
	for(size_t i = 0 ; i < brother->features.size() ; i++)
	{
		features.push_back( brother->features[i] ) ;
		rotations.push_back( brother->rotations[i] ) ;
	}
}

void InclusionFamily::addToFeatureTree(FeatureTree * f)
{
	for(size_t i = 0 ; i < features.size() ; i++)
	{
		for(size_t j = 0 ; j < features[i].size() ; j++)
		{
			if( features[i][j]->getFather() != nullptr)
				f->addFeature( features[i][j]->getFather(), features[i][j] ) ;
			else if( keepNoFatherFeatures )
				f->addFeature( f->getFeature(0), features[i][j] ) ;
		}
	}
}

std::vector<Geometry *> InclusionFamily::getFeaturesAsGeometry(size_t index) 
{
	std::vector<Geometry *> ret ;
	if(index < features.size())
	{
		for(size_t j = 0 ; j < features[index].size() ; j++)
			ret.push_back( dynamic_cast<Geometry *>( features[index][j] ) ) ;
		return ret ;
	}

	for(size_t i = 0 ; i < features.size() ; i++)
	{
		for(size_t j = 0 ; j < features[i].size() ; j++)
			ret.push_back( dynamic_cast<Geometry *>( features[i][j] ) ) ;
	}
	return ret ;
}

void InclusionFamily::place( Rectangle * box, double spacing, size_t tries, size_t seed )
{
	srand(seed) ;
	std::vector<Geometry *> brothers ;
	for(size_t i = 0 ; i < features.size() ; i++)
	{
		std::vector<Feature *> placed = placement2D( box, features[i], spacing, 0, tries, rotations[i], brothers ) ;
		features[i] = placed ;
		for(size_t j = 0 ; j < placed.size() ; j++)
			brothers.push_back(dynamic_cast<Geometry *>(placed[j])) ;
	}
}

EmbeddedInclusionFamily::EmbeddedInclusionFamily( size_t n, double rmax, double fraction, ParticleSizeDistribution * type, InclusionGenerator * geometry) : InclusionFamily(n,rmax,fraction,type,geometry)
{
	keepNoFatherFeatures = false ;

}

void EmbeddedInclusionFamily::place( Rectangle * box, double spacing, size_t tries, size_t seed )
{
	srand(seed) ;
	std::vector<Geometry *> brothers ;
	std::vector<Geometry *> fathers = father->getFeaturesAsGeometry(fatherIndex) ;
	for(size_t i = 0 ; i < features.size() ; i++)
	{
		std::vector<Feature *> placed = placement2DInInclusions( box, fathers, features[i], spacing, 0, tries, rotations[i], brothers ) ;
		features[i] = placed ;
		for(size_t j = 0 ; j < placed.size() ; j++)
		{
			brothers.push_back(dynamic_cast<Geometry *>(placed[j])) ;
			bool found = true ;
			for(size_t f = 0 ; f < father->features.size() && !found; f++)
			{
				for(size_t k = 0 ; k < father->features[f].size() && !found ; k++)
				{
					if( father->features[f][k]->in( placed[j]->getCenter() ) )
					{
						placed[j]->setFather( father->features[f][k] ) ;
						found = false ;
					}
				}
			}
		}
	}
}

MaskedInclusionFamily::MaskedInclusionFamily( size_t n, double rmax, double fraction, ParticleSizeDistribution * type, InclusionGenerator * geometry) : InclusionFamily(n,rmax,fraction,type,geometry)
{
	keepNoFatherFeatures = false ;

}

void MaskedInclusionFamily::place( Rectangle * box, double spacing, size_t tries, size_t seed )
{
	srand(seed) ;
	std::vector<Geometry *> brothers ;
	std::vector<Geometry *> fathers = father->getFeaturesAsGeometry( fatherIndex ) ;
	for(size_t i = 0 ; i < features.size() ; i++)
	{
		std::vector<Feature *> placed = placement2D( box, features[i], spacing, 0, tries, rotations[i], brothers ) ;
		features[i] = placed ;
		for(size_t j = 0 ; j < placed.size() ; j++)
		{
			brothers.push_back(dynamic_cast<Geometry *>(placed[j])) ;
			bool found = true ;
			for(size_t f = 0 ; f < father->features.size() && !found; f++)
			{
				for(size_t k = 0 ; k < father->features[f].size() && !found ; k++)
				{
					if( father->features[f][k]->in( placed[j]->getCenter() ) || father->features[f][k]->intersects(placed[j]) )
					{
						placed[j]->addToMask( father->features[f][k] ) ;
						placed[j]->setFather( father->features[f][k] ) ;
						found = false ;
					}
				}
			}
		}
	}
}

ConcentricInclusionFamily::ConcentricInclusionFamily( double r) : InclusionFamily(), width(r)
{

}

void ConcentricInclusionFamily::place( Rectangle * box, double spacing, size_t tries, size_t seed )
{
	std::vector<Feature *> feats ;
	std::vector<Feature *> fatherFeats = father->getFeatures(fatherIndex) ;
	for(size_t j = 0 ; j < fatherFeats.size() ; j++)
	{
		if(fatherFeats[j]->getRadius() > width)
			feats.push_back( new Inclusion( fatherFeats[j], fatherFeats[j]->getRadius() - width, fatherFeats[j]->getCenter() ) ) ;
	}
	features.push_back(feats) ;
}

VoronoiInclusionFamily::VoronoiInclusionFamily( double radius, double fraction, double correction, size_t n_, size_t max)
{
	n = n_ ;
	nmax = max ;
	grains.push_back( VoronoiGrain( nullptr, radius, fraction, correction ) ) ;
	copy = false ;
	generated = false ;
}

void VoronoiInclusionFamily::concatenate( InclusionFamily * brother )   
{
	if(generated)
	{
		InclusionFamily::concatenate( brother ) ;
		return ;
	}
	
	VoronoiInclusionFamily * voronoi = dynamic_cast<VoronoiInclusionFamily *>(brother) ;
	if(voronoi != nullptr)
	{
		n += voronoi->n ;
		nmax = std::max(nmax, voronoi->nmax) ;
		copy |= voronoi->copy ;
		for(size_t i = 0 ; i < voronoi->grains.size() ; i++)
			grains.push_back( voronoi->grains[i] ) ;
	}
}



void VoronoiInclusionFamily::place( Rectangle * box, double spacing, size_t tries, size_t seed )
{
	srand(seed) ;
	
	if(father == nullptr)
		features = PSDGenerator::get2DVoronoiPolygons( box, grains, n, -1, spacing, nmax, copy, spacing ) ;
	else
		features = PSDGenerator::get2DVoronoiPolygons(box, grains, father->getFeatures(fatherIndex), n, -1, spacing, nmax, copy, spacing) ;
	generated = true ;
}











}
