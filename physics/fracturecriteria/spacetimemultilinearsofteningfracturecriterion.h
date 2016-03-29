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


/*PARSE SpaceTimeNonLocalMultiLinearSofteningFractureCriterion FractureCriterion
	@string[tension_file_name] // location of the file containing the linear interpolation of stress vs strain
	@value[young_modulus] // Young modulus of the material
	@value[strain_renormalization_factor] 1e4 // renormalization factor to facilitate the geometry interactions
	@value[stress_renormalization_factor] 1e-6 // renormalization factor to facilitate the geometry interactions
*/
class SpaceTimeNonLocalMultiLinearSofteningFractureCriterion : public MaximumStrain
{
public:
	double E ;
	double renormStrain ;
	double renormStress ;
	SegmentedLine * stressStrainCurve ;
	Line * asymptote ;

	SpaceTimeNonLocalMultiLinearSofteningFractureCriterion( std::string file, double E_, double strain = 1e4, double stress = 1e-6) ;
	SpaceTimeNonLocalMultiLinearSofteningFractureCriterion( const std::vector<Point> & p, double E_, double strain = 1e4, double stress = 1e-6) ;
	SpaceTimeNonLocalMultiLinearSofteningFractureCriterion( SegmentedLine * s, Line * l, double E_, double strain = 1e4, double stress = 1e-6) ;

	virtual ~SpaceTimeNonLocalMultiLinearSofteningFractureCriterion() { delete stressStrainCurve ; delete asymptote ; } ;

	virtual FractureCriterion * getCopy() const ;

	virtual double grade(ElementState &s)  ;
    virtual double gradeAtTime(ElementState &s, double t)  ;

	void setMaximumStress(double s = -1, double smax = -1) ;
	void setMaximumStrain(double e = -1, double emax = -1) ;

} ;

/*PARSE AsymmetricSpaceTimeNonLocalMultiLinearSofteningFractureCriterion FractureCriterion
	@string[tension_file_name] // location of the file containing the linear interpolation of stress vs strain in tension
	@string[compression_file_name] // location of the file containing the linear interpolation of stress vs strain in compression
	@value[young_modulus] // Young modulus of the material
	@value[strain_renormalization_factor] 1e4 // renormalization factor to facilitate the geometry interactions
	@value[stress_renormalization_factor] 1e-6 // renormalization factor to facilitate the geometry interactions
*/
class AsymmetricSpaceTimeNonLocalMultiLinearSofteningFractureCriterion : public MaximumStrain
{
public:
	double E ;
	double renormStrain ;
	double renormStress ;
	SegmentedLine * tensileStressStrainCurve ;
	Line * tensileAsymptote ;
	SegmentedLine * compressiveStressStrainCurve ;
	Line * compressiveAsymptote ;

	AsymmetricSpaceTimeNonLocalMultiLinearSofteningFractureCriterion( std::string ftension, std::string fcompression, double E_, double strain = 1e4, double stress = 1e-6) ;
	AsymmetricSpaceTimeNonLocalMultiLinearSofteningFractureCriterion( const std::vector<Point> & ptension, const std::vector<Point> & pcompression, double E_, double strain = 1e4, double stress = 1e-6) ;

	virtual ~AsymmetricSpaceTimeNonLocalMultiLinearSofteningFractureCriterion() ; //{ delete tensileStressStrainCurve ; delete tensileAsymptote ; delete compressiveStressStrainCurve ; delete compressiveAsymptote ; } ;

	virtual FractureCriterion * getCopy() const ;

	virtual double grade(ElementState &s)  ;
	virtual double gradeAtTime(ElementState &s, double t)  ;

	void setMaximumTensileStress(double s = -1, double smax = -1) ;
	void setMaximumTensileStrain(double e = -1, double emax = -1) ;
	void setMaximumCompressiveStress(double s = 1, double smax = 1) ;
	void setMaximumCompressiveStrain(double e = 1, double emax = 1) ;

	double getTensileStrength() const ;
	double getTensileStrain() const ;
	double getCompressiveStrength() const ;
	double getCompressiveStrain() const ;
} ;


}

#endif
