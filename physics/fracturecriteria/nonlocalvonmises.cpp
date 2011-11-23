
//
// C++ Implementation: vonmises
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "nonlocalvonmises.h"
#include "../../mesher/delaunay.h"
#include "../damagemodels/damagemodel.h"
namespace Mu {

NonLocalVonMises::NonLocalVonMises(double thresh, double E, double radius, MirrorState mirroring, double delta_x, double delta_y, double delta_z) : FractureCriterion(mirroring, delta_x, delta_y, delta_z)
	, threshold(thresh), E(E)
{
	setMaterialCharacteristicRadius(radius);
}


NonLocalVonMises::~NonLocalVonMises()
{
}

double NonLocalVonMises::grade(ElementState &s)
{
	metInCompression = true ;
	metInTension = true ;
	std::pair<Vector, Vector> str( smoothedPrincipalStressAndStrain(s) ) ;
	double effectiveStiffness = E ;
	if(s.getParent()->getBehaviour() && s.getParent()->getBehaviour()->getDamageModel())
		effectiveStiffness = E*(1.-s.getParent()->getBehaviour()->getDamageModel()->getState().max()) ;
	
	double maxStress = 0 ;
	double maxStrain = 0 ;
	if( s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
	{           
		maxStress = sqrt( ( ( str.first[0] - str.first[1] ) * ( str.first[0] - str.first[1] ) + str.first[0] * str.first[0] + str.first[1] * str.first[1] ) * .5 ) ;
		maxStrain = sqrt( ( ( str.second[0] - str.second[1] ) * ( str.second[0] - str.second[1] ) + str.second[0] * str.second[0] + str.second[1] * str.second[1] ) / 6.) ;
	}
	else if( s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL )
	{
		maxStress = sqrt( ( str.first[0] - str.first[1] ) * ( str.first[0] - str.first[1] ) + ( str.first[0] - str.first[2] ) * ( str.first[0] - str.first[2] ) + ( str.first[1] - str.first[2] ) * ( str.first[1] - str.first[2] ) ) / 6 ;
		maxStrain =sqrt( ( str.second[0] - str.second[1] ) * ( str.second[0] - str.second[1] ) + ( str.second[0] - str.second[2] ) * ( str.second[0] - str.second[2] ) + ( str.second[1] - str.second[2] ) * ( str.second[1] - str.second[2] ) ) *2./ 3. ;
	}
	
	std::vector<double> scores ;
	scores.push_back(-1);
	if( maxStress >= threshold  /*|| maxStrain > threshold/effectiveStiffness && maxStrain > 0*/)
	{
		metInTension = true;
		scores.push_back(std::min(1. - std::abs( threshold / maxStress ), 1. - std::abs( threshold / maxStress ) ));
	}
	else 
			scores.push_back(-1. + std::abs( maxStress / threshold ));

	std::sort(scores.begin(), scores.end()) ;
	return scores.back() ;

}

FractureCriterion * NonLocalVonMises::getCopy() const
{
	return new NonLocalVonMises(threshold,E,  getMaterialCharacteristicRadius()) ;
}

Material NonLocalVonMises::toMaterial()
{
	Material mat ;
	return mat ;
}

}
