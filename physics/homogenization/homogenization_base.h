#ifndef HOMOGENIZATION_BASE_H
#define HOMOGENIZATION_BASE_H

#include "../../utilities/matrixops.h"
#include "../../geometry/geometry_base.h" 

namespace Mu
{

typedef enum
{
	P_BAD_INDEX,
	P_NO_FATHER,
	P_NULL_PROPERTIES,

	P_BULK_MODULUS,
	P_POISSON_RATIO,
	P_SHEAR_MODULUS,
	P_YOUNG_MODULUS,
	
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
 * They may have an unlimited number of children (stored in a standard vector), and one (or no) father.
 */
class Material
{
	std::vector<Properties> prop ;
	std::vector<Material> phases ;
	Material * father ;

public:

// CONSTRUCTORS

	/* Simple constructor, creates an empty Material with no Properties, no phases, no father */
	Material() ;

	/* Simple constructor, creates a Material from an existing set of Properties. 
	 * This constructor does not check is a Properties type is defined multiple time.
	 * @param p the Properties to use
	 */
	Material(std::vector<Properties> p) ;

	/* Simple constructor, creates a Material from an existing set of Materials. 
	 * The father of each children Material will be set to <b>this</b>
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
	
	Properties getProperties(int i) const ;
	Properties getProperties(PType t) const ;
	bool hasProperties(PType t) const ;
	int indexOfProperties(PType t, bool first = true) const ;
	void removeProperties(int i) ;
	void removeProperties(PType t) ;
	Properties searchProperties(PType t, SearchPattern p = SEARCH_TRY_ANYTHING, bool next = false) ;
	void setProperties(PType t, double v) ;
	void setProperties(Properties p) ;
	void setProperties(Material m) ;
	int sizeProperties() const ;
	double valProperties(int i) const ;
	double valProperties(PType t) const ;
	
	void addPhase(Material m) ;
	void clearPhase() ;
	Material getPhase(int i) const ;
	void mergePhase() ;
	void removePhase(int i) ;
	int sizePhase() const ;
	double valPhase(int i, PType t) const ;
	std::vector<double> valPhase(PType t) const ;

	void setFather(Material * m = NULL) ;
	Properties getFatherProperties(PType t) const ;
	bool hasFatherProperties(PType t) const ;
	double valFatherProperties(PType t) const ;

	bool equals(Material m, bool extensible = false) const ;

	bool homogenize(HomogenizationScheme s) ;

	static Matrix cauchyGreen(std::pair<double,double> prop, bool hooke, SpaceDimensionality dim) ;



} ;

} ;


#endif
