#ifndef HOMOGENIZATION_BASE_H
#define HOMOGENIZATION_BASE_H

#include "../../utilities/matrixops.h"
#include "../../geometry/geometry_base.h" 

namespace Mu
{

typedef enum
{
	P_BAD_INDEX,
	P_NULL_PROPERTIES,

	P_BULK_MODULUS,
	P_POISSON_RATIO,
	P_SHEAR_MODULUS,
	P_YOUNG_MODULUS,
	
	P_EXPANSION_COEFFICIENT,
	
	P_VOLUME,
	P_VOLUME_FRACTION,
	P_MASS,
	P_MASS_FRACTION,
	P_DENSITY,

} PType ;

typedef enum
{
	SEARCH_VOLUME_FROM_MASS_DENSITY,
	SEARCH_VOLUME_FROM_FATHER_VOLUME,
	SEARCH_VOLUME_FROM_CHILDREN_VOLUME,
	
	SEARCH_MASS_FROM_VOLUME_DENSITY,
	SEARCH_MASS_FROM_FATHER_MASS,
	SEARCH_MASS_FROM_CHILDREN_MASS,
	
	SEARCH_TRY_ANYTHING,
	SEARCH_TRY_NOTHING,

} SearchPattern ;

typedef enum
{
	ELASTICITY_DILUTED,
	ELASTICITY_GENERALIZED_DILUTED,
	ELASTICITY_INCREMENTAL,
	ELASTICITY_MORI_TANAKA,
	ELASTICITY_GENERALIZED_MORI_TANAKA,
	ELASTICITY_SELF_CONSISTENT,
	ELASTICITY_GENERALIZED_SELF_CONSISTENT,
	
	EXPANSION_HOBBS,
	EXPANSION_KERNER,
	EXPANSION_TURNER,
	
	SCHEME_DO_NOTHING,
} HomogenizationScheme ;


/* A Properties is a typed double.
 * The type of the Properties item is set one and for all at the creation.
 * However, its value may change.
 */
class Properties
{
	PType tag ;
	double value ;

public:
// CONSTRUCTORS

	/* Simple constructor, creates a useless P_NULL_PROPERTIES item */
	Properties() ;
	
	/* Simple constructor, creates a Properties item 
	 * @param t the type of the item
	 * @param v the initial value of the item
	 */
	Properties(PType t, double v) ;
	
	/* Copy constructor */
	Properties(const Properties & p) ;

// OPERATIONS ON THE TYPE
	
	/* Returns the Properties item type
	 * @return the type
	 */
	const PType type() const ;
	
	/* Checks if the item is of a specific type 
	 * @param t the type to check
	 * @return <b>true</b> if the item is of type <i>t</i>
	 */
	const bool is(PType t) const ;
	
	/* Checks if the item is of an extensible type. Two Properties of extensible types max be added together.
	 * Example extensible types are <i>P_VOLUME</i> and <i>P_MASS</i>
	 * @return <b>true</b> if the Properties is extensible
	 */
	const bool isExtensible() const ;

// OPERATIONS ON THE VALUE
	
	/* Returns the Properties value. This is not an allocator!
	 * @return the value
	 */
	double val() const ;

	/* Returns the Properties value if the item is of type <i>t</i>. This is not an allocator!
	 * @param t the type to check
	 * @return the value if the item is of type <i>t</i>, 0 otherwise
	 */
	double val(PType t) const ;

	/* Sets the item value
	 * @param v the new value
	 */
	void set(double v) ;

// OPERATIONS ON OTHER PROPERTIES

	/* Checks if two Properties item are equal. Set the <i>extensible</i> flag to <b>true</b> to check if two Properties have the same type
	 * @param p another Properties item
	 * @param extensible if true, then the equality of the values will not be checked, but only the equality of the the types
	 * @return true if the two Properties are of same type and their values are equals
	 */
	bool equals(Properties p, bool extensible = false) const ;	

} ;


/* A Material is a collection of Properties, using the standard c++ vector type. 
 * Public methods allow the direct access to a specific Properties type. 
 * <br>
 * It is possible to have in the same Material different Properties item with the same type.
 * However, this is highly not recommanded.
 * In general, only the first Properties of a specific type will be detected and used.
 * The other items of the same type will be ignored. 
 * <br>
 * Materials are organized as a tree, which helps to represents matrix-inclusions or similar morphoplogies.
 * They may have an unlimited number of children (stored in a standard vector).
 * However, information to the upper levels (father) is not stored.
 */
class Material
{
	std::vector<Properties> prop ;
	std::vector<Material> phases ;

public:

// CONSTRUCTORS

	/* Simple constructor, creates an empty Material with no Properties, no phases */
	Material() ;

	/* Simple constructor, creates a Material from an existing set of Properties. 
	 * This constructor does not check is a Properties type is defined multiple time.
	 * @param p the Properties to use
	 */
	Material(std::vector<Properties> p) ;

	/* Simple constructor, creates a Material from an existing set of Materials. 
	 * @param p the Material to use
	 */
	Material(std::vector<Material> p) ;
	
	Material(const Matrix & cauchy) ;
	
// OPERATIONS ON PROPERTIES
	
	/* Adds a Properties to the collection. This method does not check the unicity of the type. 
	 * @param p a Properties item
	 */
	void addProperties(Properties p) ;

	/* Creates a new Properties and add it to the collection. This method does not check the unicity of the type. 
	 * @param t the type of the new Properties item
	 * @param v the value of the new Properties item
	 */
	void addProperties(PType t, double v) ;

	/* Makes sure that all Properties type defined in the Material are only defined once. 
	 * This method will keep the first instance of each Properties type. 
	 * If the <i>last</i> flag is set to <b>false</b>, this method will keep the last instance of each Properties type.
	 * @param last indicates which instance to keep
	 */
	void aloneProperties(bool last = true) ;

	/* Makes sure that the chosen Properties type is only defined once in the Material. 
	 * This method will delete the last instance of this type until unicity is reached. 
	 * If the <i>last</i> flag is set to <b>false</b>, this method will delete the first instance of each Properties type.
	 * @param t the Properties type to check
	 * @param last indicates which instance to delete
	 */
	void aloneProperties(PType t, bool last = true) ;

	/* Clears all Properties */
	void clearProperties() ;
	
	/* Access the <i>ith</i> properties in the collection. If <i>i</i> is out of bound, it returns a P_BAD_INDEX properties.
	 * @param i index of the Properties to reach
	 * @return the <i>ith</i> properties
	 */
	Properties getProperties(int i) const ;

	/* Access to the first Properties of given type in the collection, 
	 * or a P_BAD_INDEX properties if there is no Properties of the given type.
	 * @param t the Properties type to get
	 * @return a Properties of type <i>t</i> or P_BAD_INDEX
	 */
	Properties getProperties(PType t) const ;

	/* Checks if a Properties of type <i>t</i> is contained in the collection. 
	 * This method uses the result of indexOfProperties(t).
	 * Use hasProperties(t) when you don't need to access the specific Properties, but just need to check the existence.
	 * @param t the Properties type to check
	 * @return <b>true</b> if there is at least one Properties item of type <i>t</i> in the collection, <b>false</b> otherwise.
	 */
	bool hasProperties(PType t) const ;

	/* Search for the first or the last Properties of type <i>t</i> contained in the collection.
	 * @param t the Properties type to check
	 * @param first if <b>true</b> (default), get the first Properties, get the last otherwise
	 * @return the index of the first (or last) Properties item of type <i>t</i> met in the collection, or -1 if not found
	 */
	int indexOfProperties(PType t, bool first = true) const ;
	
	/* Removes the <i>ith</i> Properties in the collection, do nothing if <i>i</i> is out of bound.
	 * @param i the index to remove
	 */
	void removeProperties(int i) ;

	/* Removes the first Properties of type <i>t</i> in the collection, do nothing if no Properties of type <i>t</i> is to be found
	 * @param t the Properties type to remove
	 */
	void removeProperties(PType t) ;

	/* Search for a Properties of type <i>t</i>, and tries to deduce it from other informations when possible.
	 * Various search strategies are possible, depending on the Properties type to find.
	 * <ul>	<li>If p = SEARCH_TRY_ANYTHING , all known search strategies for the type <i>t</i> will be applied.</li>
	 *	<li>If p = SEARCH_TRY_NOTHING , this method just checks if a Properties of type <i>t</i> exists within the collection.</li></ul>
	 * <br>
	 * Different search strategies are possible depending on the givent type <i>t</i>
	 * <ul>	<li>P_BULK_MODULUS and P_SHEAR_MODULUS are built from P_YOUNG_MODULUS and P_POISSON_RATIO (and vice-versa)</li>
	 *	<li>P_EXPANSION_COEFFICIENT is created and set to 0 if not found.</li>
	 *	<li>P_DENSITY is built from P_VOLUME and P_MASS</li>
	 *	<li>P_VOLUME_FRACTION (or P_MASS_FRACTION) are built from the P_VOLUME (or P_MASS) of the father and itself</li>
	 *	<li>P_VOLUME may be built from one of the following strategies:
	 *	  <ul>	<li>SEARCH_VOLUME_FROM_MASS_DENSITY uses the P_MASS and the P_DENSITY of the current Material</li>
	 *		<li>SEARCH_VOLUME_FROM_FATHER_VOLUME uses the P_VOLUME_FRACTION of the current Material and the P_VOLUME of the father</li>
	 *		<li>SEARCH_VOLUME_FROM_CHILDREN_VOLUME makes the sum of the P_VOLUME of all children of the current Material</li></ul>
	 *	<li>Similar strategies exist for P_MASS</li></ul>
	 * <br>
	 * 
	 *
	 */
	Properties searchProperties(PType t, SearchPattern p = SEARCH_TRY_ANYTHING, bool next = false, Material father = Material()) ;
	void setProperties(PType t, double v) ;
	void setProperties(Properties p) ;
	void setProperties(Material m) ;
	int sizeProperties() const ;
	double valProperties(int i) const ;
	double valProperties(PType t) const ;
	
// OPERATIONS ON CHILDREN	
	
	void addPhase(Material m) ;
	void clearPhase() ;
	Material getPhase(int i) const ;
	void mergePhase() ;
	void removePhase(int i) ;
	int sizePhase() const ;
	double valPhase(int i, PType t) const ;
	std::vector<double> valPhase(PType t) const ;

// OPERATIONS ON OTHER MATERIAL

	bool equals(Material m, bool extensible = false) const ;

// HOMOGENIZATION

	bool homogenize(HomogenizationScheme s) ;

// CONVENIENCE STATIC METHOD

	static Matrix cauchyGreen(std::pair<double,double> prop, bool hooke, SpaceDimensionality dim) ;



} ;

} ;


#endif
