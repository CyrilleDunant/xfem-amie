//
// C++ Interface: mechanical homogenization
//
// Description: 
//
//
// Author:  Alain Giorla, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef PROPERTIES_BASE_H
#define PROPERTIES_BASE_H

#include "../../utilities/matrixops.h"

namespace Mu
{


typedef enum
{
	TAG_NULL, 
	TAG_UNIVERSAL, 
	TAG_VOLUME,	
	TAG_VOLUME_FRACTION, 
	TAG_VOLUME_TOTAL,
	TAG_MASS,
	TAG_MASS_FRACTION,
	TAG_MASS_TOTAL,
	TAG_DENSITY,
	TAG_YOUNG_MODULUS,
	TAG_POISSON_RATIO,
	TAG_BULK_MODULUS,	
	TAG_SHEAR_MODULUS,
	TAG_LAME_COEFFICIENT,
	TAG_STRAIN,
	TAG_STRESS,
	TAG_EXPANSION_COEFFICIENT, 
	TAG_CRACK_DENSITY, 
	TAG_CRACK_LENGTH,
	TAG_ELLIPSE_A,
	TAG_ELLIPSE_B,
	TAG_AREA,
	TAG_PERIMETER,
	TAG_ELLIPSE_FIRST_COMPLETE_INTEGRAL,
	TAG_ELLIPSE_SECOND_COMPLETE_INTEGRAL,
	TAG_CIRCLE_RADIUS,
} Tag ;

typedef enum
{
	MAT_DUMMY,
	MAT_AGGREGATE,
	MAT_CEMENT,
} PredefinedMaterial ;

class Properties
{
protected:
	Tag ptag ;
	double p ;

public:
	/* \brief default constructor for a VOID_PROP properties item*/
	Properties() ;
	/* \brief constructor for an ABSTRACT properties item */
	Properties(double v) ;
	/* \brief constructor for an specific properties item */
	Properties(Tag p, double v) ;
	/* \brief copy constructor */
	Properties(const Properties & p) ;

	double val() const {return p ; } ;
	Tag tag() const {return ptag ; } ;
	bool is(Tag t) const {return ptag == t ; } ;
	void kill() {ptag = TAG_NULL ; } ;
	bool isNull() const {return ptag == TAG_NULL ; } ;

	void set(Tag t) {ptag = t ; } ;
	void set(double v) {p = v ; } ;

	void print() ;
} ;


/* \brief vector of Properties */
class Material : public std::vector<Properties>
{
protected:
	std::vector<Tag> tagset ;

public:
	Material() ;
	Material(PredefinedMaterial mat) ;
	Material(const Properties & p) ;
	Material(const std::vector<Properties> & p) ;
	Material(const Matrix & cauchy) ;

	std::vector<int> getIndex(Tag t) const ;
	int getIndex(Tag t, int i) const ;
	double val(Tag t, int i) const ;

	bool replace(Properties p) ;
	bool replaceForce(Properties p) ;
	bool isSet(Tag t) const ;
	bool set(Tag t) ;
	bool set(Tag t, int i) ;

	bool combine(Material m, std::vector<Tag> compare, Tag combine) ;
	bool merge(Material m) ;

	bool findMissing(std::vector<Tag> t) ;
	bool findMissing(Tag t) ;

	void print() ;
} ;


} ;


#endif

