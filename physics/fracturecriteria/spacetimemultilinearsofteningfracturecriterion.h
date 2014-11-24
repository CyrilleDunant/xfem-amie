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
#ifndef SPACE_TIME_MULTILINEAR_SOFTENING_FRACTURE_CRITERION_H__
#define SPACE_TIME_MULTILINEAR_SOFTENING_FRACTURE_CRITERION_H__

#include "fracturecriterion.h"
#include "maxstrain.h"
#include "../../mesher/delaunay_3d.h"

namespace Amie {


class SpaceTimeNonLocalMultiLinearSofteningFractureCriterion : public MaximumStrain
{
public:
	SegmentedLine * stressStrainCurve ;
	Line * asymptote ;
	double E ;
	double renormStrain ;
	double renormStress ;

	SpaceTimeNonLocalMultiLinearSofteningFractureCriterion( std::string file, double E_, double strain = 1e4, double stress = 1e-6, MirrorState mirroring = NO_MIRROR, double delta_x = 0, double delta_y = 0, double delta_z = 0) ;
	SpaceTimeNonLocalMultiLinearSofteningFractureCriterion( const std::vector<Point> & p, double E_, double strain = 1e4, double stress = 1e-6, MirrorState mirroring = NO_MIRROR, double delta_x = 0, double delta_y = 0, double delta_z = 0) ;
	SpaceTimeNonLocalMultiLinearSofteningFractureCriterion( SegmentedLine * s, Line * l, double E_, double strain = 1e4, double stress = 1e-6, MirrorState mirroring = NO_MIRROR, double delta_x = 0, double delta_y = 0, double delta_z = 0) ;

	virtual ~SpaceTimeNonLocalMultiLinearSofteningFractureCriterion() { delete stressStrainCurve ; delete asymptote ; } ;

	virtual FractureCriterion * getCopy() const ;

	virtual double grade(ElementState &s)  ;

	void setMaximumStress(double s = -1, double smax = -1) ;
	void setMaximumStrain(double e = -1, double emax = -1) ;

} ;

class AsymmetricSpaceTimeNonLocalMultiLinearSofteningFractureCriterion : public MaximumStrain
{
public:
	SegmentedLine * tensileStressStrainCurve ;
	Line * tensileAsymptote ;
	SegmentedLine * compressiveStressStrainCurve ;
	Line * compressiveAsymptote ;
	double E ;
	double renormStrain ;
	double renormStress ;

	AsymmetricSpaceTimeNonLocalMultiLinearSofteningFractureCriterion( std::string ftension, std::string fcompression, double E_, double strain = 1e4, double stress = 1e-6, MirrorState mirroring = NO_MIRROR, double delta_x = 0, double delta_y = 0, double delta_z = 0) ;
	AsymmetricSpaceTimeNonLocalMultiLinearSofteningFractureCriterion( const std::vector<Point> & ptension, const std::vector<Point> & pcompression, double E_, double strain = 1e4, double stress = 1e-6, MirrorState mirroring = NO_MIRROR, double delta_x = 0, double delta_y = 0, double delta_z = 0) ;

	virtual ~AsymmetricSpaceTimeNonLocalMultiLinearSofteningFractureCriterion() { delete tensileStressStrainCurve ; delete tensileAsymptote ; delete compressiveStressStrainCurve ; delete compressiveAsymptote ; } ;

	virtual FractureCriterion * getCopy() const ;

	virtual double grade(ElementState &s)  ;

	void setMaximumTensileStress(double s = -1, double smax = -1) ;
	void setMaximumTensileStrain(double e = -1, double emax = -1) ;
	void setMaximumCompressiveStress(double s = 1, double smax = 1) ;
	void setMaximumCompressiveStrain(double e = 1, double emax = 1) ;
} ;


}

#endif
