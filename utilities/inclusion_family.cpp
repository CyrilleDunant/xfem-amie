#include "inclusion_family.h"
#include "placement.h"
#include "configuration.h"
#include "enumeration_translator.h"

namespace Amie
{

void InclusionFamilyTree::insert( InclusionFamily * family ) 
{
	brothers.push_back( family ) ;
	family->brothers = this ;
}

std::vector<size_t> InclusionFamilyTree::getAllVoronoiFamilies() 
{
	std::vector<size_t> index ;
	for(size_t i = 0 ; i < brothers.size() ; i++)
	{
		if( dynamic_cast<VoronoiInclusionFamily *>( brothers[i] ) != nullptr )
			index.push_back(i) ;
	}
	return index ;
}

void InclusionFamilyTree::place( Rectangle * box, double spacing, size_t tries, size_t seed ) 
{
	std::vector<Geometry *> placedFeatures ;
	for(size_t i = 0 ; i < brothers.size() ; i++)
	{
		if(!brothers[i]->placed)
			brothers[i]->place( box, spacing, tries, seed, placedFeatures ) ;
		if(brothers[i]->placed)
		{
			std::vector<Geometry *> current = brothers[i]->getFeaturesAsGeometry() ;
			for(size_t j = 0 ; j < current.size() ; j++)
				placedFeatures.push_back( current[j] ) ;
		}
	}
}

void InclusionFamilyTree::concatenate()
{
	features.clear() ;
	for(size_t i = 0 ; i < brothers.size() ; i++)
	{
		if(brothers[i]->placed)
		{
			features.push_back( brothers[i]->getAddedFeatures() ) ;
			if(brothers[i]->sons != nullptr)
			{
				brothers[i]->sons->concatenate() ;
				for(size_t j = 0 ; j < brothers[i]->sons->features.size() ; j++)
					features.push_back( brothers[i]->sons->features[j] ) ;
			}
		}
	}
}

std::vector<Feature *> InclusionFamily::getAddedFeatures() 
{
	std::vector<Feature *> ret ;
	if(added.size() < features.size() || !placed)
		return ret ;

	for(size_t i = 0 ; i < features.size() ; i++)
	{
		if(added[i])
			ret.push_back(features[i]) ;
	}

	return ret ;
}

std::vector<Geometry *> InclusionFamilyTree::getFeaturesAsGeometry(size_t index) 
{
	std::vector<Geometry *> ret ;
	if(index >= features.size())
		return ret ;
	for(size_t i = 0 ; i < features[index].size() ; i++)
		ret.push_back( dynamic_cast<Geometry *>( features[index][i] ) ) ;
	return ret ;
}

InclusionFamily::InclusionFamily() 
{
	keepNoFatherFeatures = true ;
	placed = false ;
	brothers = nullptr ;
	sons = nullptr ;
	sampler = nullptr ;
	rotations = 0 ;	
	factors = -1 ;
}

InclusionFamily::InclusionFamily( size_t n, double rmax, double fraction, ParticleSizeDistribution * type, InclusionGenerator * geometry) 
{
	keepNoFatherFeatures = true ;
	placed = false ;

	if(type == nullptr)
		type = new ConstantSizeDistribution() ;
	if(geometry == nullptr)
		geometry = new InclusionGenerator() ;

	std::vector<Inclusion *> incs = PSDGenerator::get2DConcrete( rmax, std::sqrt(fraction), n, type, 1. ) ;
	features = geometry->convert( incs ) ;
	rotations = geometry->authorizeRotationsDuringPlacement ;
	factors = -1 ;

	brothers = nullptr ;
	sons = nullptr ;
	sampler = nullptr ;
}

InclusionFamily::InclusionFamily(Feature * matrix) 
{
	keepNoFatherFeatures = true ;
	placed = false ;

	features.push_back(matrix) ;
	rotations = 0 ;
	factors = -1 ;

	brothers = nullptr ;
	sons = nullptr ;
	sampler = nullptr ;
}

std::vector<Geometry *> InclusionFamily::getFeaturesAsGeometry() 
{
	std::vector<Geometry *> ret ;
	for(size_t i = 0 ; i < features.size() ; i++)
		ret.push_back( dynamic_cast<Geometry *>( features[i] ) ) ;
	return ret ;
}

void InclusionFamily::setBehaviour( Form * behaviour, bool copy)
{
	for(size_t k = 0 ; k < features.size() ; k++)
		features[k]->setBehaviour( copy ? behaviour->getCopy() : behaviour ) ;
}

void InclusionFamily::place( Rectangle * box, double spacing, size_t tries, size_t seed, std::vector<Geometry *> & placedFeatures )
{
	if(placed)
		return ;

	features = placement2D( box, features, spacing, 0, tries, rotations, placedFeatures, seed ) ;
	placed = true ;
}

InclusionFamily * InclusionFamily::getFather() 
{
	if( brothers )
		return brothers->father ;
	return nullptr ;
}

void InclusionFamily::addToFeatureTree(FeatureTree * f)
{

	if(!placed)
		return ;

	added.resize( features.size() ) ;

	for(size_t i = 0 ; i < features.size() ; i++)
	{
		added[i] = false ;

		std::vector<Point> bounds = features[i]->getBoundingBox() ;
		bool in = false ;
		for(size_t i = 0 ; i < bounds.size() ; i++)
			in |= f->inRoot( bounds[i] ) ;

		if(!in)
			continue ;

		if( features[i]->getFather() != nullptr || keepNoFatherFeatures )
		{
			std::vector<Feature *> feats = f->getCoOccuringFeatures( features[i] ) ;
			for(size_t j = 0 ; j < feats.size() ; j++)
			{
				if( this->inNeighbours( feats[j] ) )
					continue ;

				if( feats[j] != f->getFeature(0) && (feats[j]->intersects( features[i] ) || features[i]->in(feats[j]->getCenter()) || feats[j]->in(features[i]->getCenter()) ) && feats[j]->getFather() != features[i] )
				{
					features[i]->addChild( feats[j] ) ;
					if(feats[j]->getFather() == f->getFeature(0) || feats[j]->getFather() == nullptr)
					{
						feats[j]->setFather( features[i] ) ;
						f->getFeature(0)->removeChild( features[i] ) ;
					}
				}
			}
			added[i] = true ;
		}
		if( features[i]->getFather() != nullptr)
			f->addFeature( features[i]->getFather(), features[i] ) ;
		else if( keepNoFatherFeatures )
			f->addFeature( f->getFeature(0), features[i] ) ;
		if(factors > 0)
			f->setSamplingFactor( features[i], factors ) ;
		if(sampler != nullptr)
			f->setSampler( features[i], sampler ) ;
	}
}

bool InclusionFamily::inNeighbours( Feature * f ) 
{
	return std::find( features.begin(), features.end(), f ) != features.end() ;
}


EmbeddedInclusionFamily::EmbeddedInclusionFamily( size_t n, double rmax, double fraction, ParticleSizeDistribution * type, InclusionGenerator * geometry) : InclusionFamily(n,rmax,fraction,type,geometry)
{
	keepNoFatherFeatures = false ;
}

void EmbeddedInclusionFamily::place( Rectangle * box, double spacing, size_t tries, size_t seed, std::vector<Geometry *> & placedFeatures )
{
	InclusionFamily * father = getFather() ;
	if(placed || father == nullptr)
		return ;

	features = placement2DInInclusions( box, father->getFeaturesAsGeometry(), features, spacing, 0, tries, rotations, placedFeatures, seed ) ;
	for(size_t j = 0 ; j < features.size() ; j++)
	{
		bool found = true ;
		for(size_t k = 0 ; k < father->size() && found; k++)
		{
			if( father->features[k]->in( features[j]->getCenter() ) )
			{
				features[j]->setFather( father->features[k] ) ;
				found = false ;
			}
		}
	}
	placed = true ;
}

MaskedInclusionFamily::MaskedInclusionFamily( size_t n, double rmax, double fraction, ParticleSizeDistribution * type, InclusionGenerator * geometry) : InclusionFamily(n,rmax,fraction,type,geometry)
{
	keepNoFatherFeatures = false ;

}

void MaskedInclusionFamily::place( Rectangle * box, double spacing, size_t tries, size_t seed, std::vector<Geometry *> & placedFeatures )
{
	InclusionFamily * father = getFather() ;
	if(placed || father == nullptr)
		return ;

	features = placement2D( box, features, spacing, 0, tries, rotations, placedFeatures, seed ) ;
	for(size_t j = 0 ; j < features.size() ; j++)
	{
		bool found = true ;
		for(size_t k = 0 ; k < father->size() && found; k++)
		{
			if( father->features[k]->in( features[j]->getCenter() ) )
			{
				features[j]->addToMask( father->features[k] ) ;
				features[j]->setFather( father->features[k] ) ;
				found = false ;
			}
		}
	}
	placed = true ;
}

ConcentricInclusionFamily::ConcentricInclusionFamily( double r) : InclusionFamily(), width(r)
{

}

void ConcentricInclusionFamily::place( Rectangle * box, double spacing, size_t tries, size_t seed, std::vector<Geometry *> & placedFeatures )
{
	InclusionFamily * father = getFather() ;
	if(placed || father == nullptr)
		return ;

	for(size_t j = 0 ; j < father->size() ; j++)
	{
		if(father->features[j]->getRadius() > width && dynamic_cast<Inclusion *>(father->features[j]) )
			features.push_back( new Inclusion( father->features[j], father->features[j]->getRadius() - width, father->features[j]->getCenter() ) ) ;
	}
	placed = true ;
}

VoronoiInclusionFamily::VoronoiInclusionFamily( double radius, double fraction, double correction, size_t n_, double out, double width, size_t max) : InclusionFamily(), grains(nullptr, radius, fraction, n_, correction)
{
	nmax = max ;
	copy = false ;
	outside = out ;
	if(out < 0) { outside = radius*1.5 ; }
	interface = width ;
	rotations = 0 ;
	factors = -1 ;
}


void VoronoiInclusionFamily::place( Rectangle * box, double spacing, size_t tries, size_t seed, std::vector<Geometry *> & placedFeatures )
{
	if(placed || brothers == nullptr)
		return ;

	std::vector<size_t> index = brothers->getAllVoronoiFamilies() ;
	std::vector<VoronoiGrain> phases ;
	double rmin = grains.radius*2 ;
	int idx = -1 ;
	for(size_t i = 0 ; i < index.size() ; i++)
	{
		VoronoiInclusionFamily * phase = dynamic_cast<VoronoiInclusionFamily *>( brothers->brothers[index[i]] ) ;
		phases.push_back( phase->grains ) ;
		outside = std::max( outside, phase->outside ) ;
		interface = std::max( interface, phase->interface ) ;
		if( phase->grains.radius < rmin )
		{
			idx = i ;
			rmin = phase->grains.radius ;
		}
	}

	double realn = phases[idx].n / phases[idx].fraction ;

	std::vector<std::vector<Feature *> > polys ;

	if(getFather() == nullptr)
		polys = PSDGenerator::get2DVoronoiPolygons( box, phases, realn, spacing > 0 ? spacing : -1, outside, nmax, copy, interface, seed ) ;
	else
		polys = PSDGenerator::get2DVoronoiPolygons( box, phases, getFather()->features, realn, spacing > 0 ? spacing : -1, outside, nmax, copy, interface, seed) ;

	for(size_t i = 0 ; i < index.size() ; i++)
	{
		VoronoiInclusionFamily * phase = dynamic_cast<VoronoiInclusionFamily *>( brothers->brothers[index[i]] ) ;
		phase->features = polys[i] ;
		phase->placed = true ;
	}

}

bool VoronoiInclusionFamily::inNeighbours( Feature * f ) 
{
	std::vector<size_t> index = brothers->getAllVoronoiFamilies() ;
	for(size_t i = 0 ; i < index.size() ; i++)
		if( brothers->brothers[index[i]]->InclusionFamily::inNeighbours( f ) ) { return true ; }
	return InclusionFamily::inNeighbours( f ) ;
}

FileDefinedCircleInclusionFamily::FileDefinedCircleInclusionFamily( int n, std::string file, std::string column_1, std::string column_2, std::string column_3 ) : InclusionFamily()
{
	std::vector<std::string> columns ;
	columns.push_back(column_1) ;
	columns.push_back(column_2) ;
	columns.push_back(column_3) ;
	GranuloFromFile reader( file, columns ) ;

	int realn = reader.getFieldValues("radius").size() ;
	if(n >= 0 && n < realn ) { realn = n ; }
	features = reader.getFeatures( CIRCLE, realn ) ;
	rotations = 0 ;
	factors = -1 ;
	placed = true ;
}

void FileDefinedCircleInclusionFamily::place( Rectangle * box, double spacing, size_t tries, size_t seed, std::vector<Geometry *> & placedFeatures ) { } 

FileDefinedPolygonInclusionFamily::FileDefinedPolygonInclusionFamily( int n, std::string file )  : InclusionFamily()
{
	PolygonGranuloFromFile reader( file ) ;

	features = reader.getFeatures( SPACE_TWO_DIMENSIONAL, nullptr ) ;
	while((int) features.size() > n && n >= 0)
		features.pop_back() ;
	rotations = 0 ;
	factors = -1 ;
	placed = true ;
}

void FileDefinedPolygonInclusionFamily::place( Rectangle * box, double spacing, size_t tries, size_t seed, std::vector<Geometry *> & placedFeatures ) { } 


}
