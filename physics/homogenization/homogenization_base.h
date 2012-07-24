#ifndef HOMOGENIZATION_BASE_H
#define HOMOGENIZATION_BASE_H

#include "../../utilities/matrixops.h"
#include "../../geometry/geometry_base.h" 

namespace Mu
{

typedef enum
{
	P_BAD_INDEX,
	P_nullptr_PROPERTIES,

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
	PLANE_STRESS,
	PLANE_STRAIN,
	PLANE_STRESS_FREE_G,
} planeType ; 

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


/** \brief Typed double for the Homogenization library.
 *
 * A Properties is a typed double.
 * The type of the Properties item is set one and for all at the creation.
 * However, its value may change.
 */
class Properties
{
	PType tag ;
	double value ;

public:
// CONSTRUCTORS

	/** \brief Simple constructor, creates a useless P_nullptr_PROPERTIES item */
	Properties() ;
	
	/** \brief Simple constructor, creates a Properties item 
	 * @param t the type of the item
	 * @param v the initial value of the item
	 */
	Properties(PType t, double v) ;
	
	/** \brief Copy constructor */
	Properties(const Properties & p) ;

// OPERATIONS ON THE TYPE
	
	/** \brief Returns the Properties item type
	 * @return the type
	 */
	PType type() const ;
	
	/** \brief Checks if the item is of a specific type 
	 * @param t the type to check
	 * @return <b>true</b> if the item is of type <i>t</i>
	 */
	bool is(PType t) const ;
	
	/** \brief Checks if the item is of an extensible type. 
	 *
	 * Two Properties of extensible types max be added together.
	 * Example extensible types are <i>P_VOLUME</i> and <i>P_MASS</i>
	 * @return <b>true</b> if the Properties is extensible
	 */
	bool isExtensible() const ;

// OPERATIONS ON THE VALUE
	
	/** \brief Returns the Properties value. This is not an allocator!
	 * @return the value
	 */
	double val() const ;

	/** \brief Returns the Properties value if the item is of type <i>t</i>. This is not an allocator!
	 * @param t the type to check
	 * @return the value if the item is of type <i>t</i>, 0 otherwise
	 */
	double val(PType t) const ;

	/** \brief Sets the item value
	 * @param v the new value
	 */
	void set(double v) ;

// OPERATIONS ON OTHER PROPERTIES

	/** \brief Checks if two Properties item are equal. 
	 *
	 * Set the <i>extensible</i> flag to <b>true</b> to check if two Properties have the same type
	 * @param p another Properties item
	 * @param extensible if true, then the equality of the values will not be checked, but only the equality of the the types
	 * @return true if the two Properties are of same type and their values are equals
	 */
	bool equals(Properties p, bool extensible = false) const ;	

} ;


/** \brief A Material is a collection of Properties
 *
 * It uses the standard c++ vector type. 
 * Public methods allow the direct access to a specific Properties type. 
 * <br>
 * It is possible to have in the same Material different Properties item with the same type.
 * However, this is highly not recommanded.
 * In general, only the first Properties of a specific type will be detected and used.
 * The other items of the same type will be ignored. 
 * <br>
 * Material items are able to describe multi-phase (composite) materials.
 * A Material can have any number of phases, stored in a vector.
 * Each phase in itself is a Material item, thus creating a tree-like architecture.
 */
class Material
{
	std::vector<Properties> prop ;
	std::vector<Material> phases ;

public:

// CONSTRUCTORS

	/** \brief Simple constructor, creates an empty Material with no Properties, no phases */
	Material() ;

	/** \brief Simple constructor, creates a Material from an existing set of Properties. 
	 * This constructor does not check is a Properties type is defined multiple time.
	 * @param p the Properties to use
	 */
	Material(std::vector<Properties> p) ;

	/** \brief Simple constructor, creates a Material from an existing set of Materials, using them as phases.
	 * @param p the Material to use
	 */
	Material(std::vector<Material> p) ;
	
	/** \brief Simple constructor, creates a Material from a Cauchy-Green tensor.
	 * The created Material will have mechanical Properties (Young Modulus, Poisson ration, Bulk and Shear Moduli)
	 * corresponding to the given matrix.
	 * @param cauchy a Cauchy-Green elastic tensor
	 */
	Material(const Matrix & cauchy) ;
	
// OPERATIONS ON PROPERTIES
	
	/** \brief Adds a Properties to the collection. This method does not check the unicity of the type. 
	 * @param p a Properties item
	 */
	void addProperties(Properties p) ;

	/** \brief Creates a new Properties and add it to the collection. This method does not check the unicity of the type. 
	 * @param t the type of the new Properties item
	 * @param v the value of the new Properties item
	 */
	void addProperties(PType t, double v) ;

	/** \brief Makes sure that all Properties type defined in the Material are only defined once. 
	 * This method will keep the first instance of each Properties type. 
	 * If the <i>last</i> flag is set to <b>false</b>, this method will keep the last instance of each Properties type.
	 * @param last indicates which instance to keep
	 */
	void aloneProperties(bool last = true) ;

	/** \brief Makes sure that the chosen Properties type is only defined once in the Material. 
	 * This method will delete the last instance of this type until unicity is reached. 
	 * If the <i>last</i> flag is set to <b>false</b>, this method will delete the first instance of each Properties type.
	 * @param t the Properties type to check
	 * @param last indicates which instance to delete
	 */
	void aloneProperties(PType t, bool last = true) ;

	/** \brief Clears all Properties */
	void clearProperties() ;
	
	/** \brief Access the <i>ith</i> properties in the collection. If <i>i</i> is out of bound, it returns a P_BAD_INDEX properties.
	 * @param i index of the Properties to reach
	 * @return the <i>ith</i> properties
	 */
	Properties getProperties(int i) const ;

	/** \brief Access to the first Properties of given type in the collection, 
	 * or a P_BAD_INDEX properties if there is no Properties of the given type.
	 * @param t the Properties type to get
	 * @return a Properties of type <i>t</i> or P_BAD_INDEX
	 */
	Properties getProperties(PType t) const ;

	/** \brief Checks if a Properties of type <i>t</i> is contained in the collection. 
	 * This method uses the result of indexOfProperties(t).
	 * Use hasProperties(t) when you don't need to access the specific Properties, but just need to check the existence.
	 * @param t the Properties type to check
	 * @return <b>true</b> if there is at least one Properties item of type <i>t</i> in the collection, <b>false</b> otherwise.
	 */
	bool hasProperties(PType t) const ;

	/** \brief Search for the first or the last Properties of type <i>t</i> contained in the collection.
	 * @param t the Properties type to check
	 * @param first if <b>true</b> (default), get the first Properties, get the last otherwise
	 * @return the index of the first (or last) Properties item of type <i>t</i> met in the collection, or -1 if not found
	 */
	int indexOfProperties(PType t, bool first = true) const ;
	
	/** \brief Removes the <i>ith</i> Properties in the collection, do nothing if <i>i</i> is out of bound.
	 * @param i the index to remove
	 */
	void removeProperties(int i) ;

	/** \brief Removes the first Properties of type <i>t</i> in the collection, do nothing if no Properties of type <i>t</i> is to be found
	 * @param t the Properties type to remove
	 */
	void removeProperties(PType t) ;

	/** \brief Search for a Properties of type <i>t</i>, and tries to deduce it from other informations when possible.
	 * If a Properties of type <i>t</i> already exists, then no further search is performed.
	 * <br>
	 * Otherwise, the method search for other Properties that are required to deduce <i>t</i>.
	 * The behaviour regarding this search depends on the search strategy <i>p</i>.
	 * <ul>	<li>SEARCH_TRY_ANYTHING	=> searchProperties() is called to found the required Properties.</li>
	 *	<li>SEARCH_TRY_NOTHING	=> no further search is performed if the required Properties is not found.</li></ul>
	 *
	 * Other search strategies exists for specific Properties type.
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
	 * In order to avoid an infinite search recursion, the flag <i>next</i> is set to <b>false</b> after the first recursion, 
	 * except in some specific cases, where a search strategy different from SEARCH_TRY_ANYTHING is specified.
	 * <br>
	 * If the search was successful, a new Properties of type <i>t</i> is added to the collection and returned. 
	 * If not, a P_BAD_INDEX item is returned.
	 *
	 * @param t the Properties type to build
	 * @param p the search strategy to use
	 * @param next if <b>true</b>, then the method is allowed to try to continue the search and build the Properties required to build <i>t</i>
	 * @param father the upper level Material (not necesseraly needed)
	 * 
	 * @return a Properties item of type <i>t</i> (if successful) or type P_BAD_INDEX (if not)
	 */
	Properties searchProperties(PType t, SearchPattern p = SEARCH_TRY_ANYTHING, bool next = false, Material father = Material()) ;
	
	/** \brief Set the first Properties item of type <i>t</i> in the collection to the new value <i>v</i>. 
	 * If no item of type <i>t</i> exists in the collection, a new Properties item will be created, with type <i>t</i> and value <i>v</i>.
	 * @param t the Properties type to set
	 * @param v the new value
	 */
	void setProperties(PType t, double v) ;
	
	/** \brief Set the first Properties item of the same type than <i>p</i> to the value of <i>p</i>, or adds <i>p</i> to the collection
	 * if no item of the same Properties type are to be found
	 * @param p a Properties item
	 */
	void setProperties(Properties p) ;

	/** \brief Set all the Properties of the current Material to the value of the Properties found in <i>m</i>.
	 * If <i>m</i> contains a Properties type which is not found in the current Material, then this Properties is added to the current Material.
	 * @param m a Material to copy
	 */
	void setProperties(Material m) ;

	/** \brief Convenience method to access the size of the Properties collection
	 * @return the number of Properties in the collection
	 */
	int sizeProperties() const ;

	/** \brief Convenience method to get the value of the <i>ith</i> Properties in the collection
	 * @param i the index to reach
	 * @return the value of the Properties, or 0 if <i>i</i> is out of bound
	 */
	double valProperties(int i) const ;

	/** \brief Convenience method to get the value of the first Properties of type <i>t</i> in the collection
	 * @param t the Properties type to reach
	 * @return the value of the Properties, or 0 if no Properties of type <i>t</i> is found
	 */
	double valProperties(PType t) const ;
	
// OPERATIONS ON CHILDREN	
	
	/** \brief Add a phase in the Material. 
	 * If the current Material already has a phase which has the same intrinsec Properties of the phase to add
	 * and if <i>merge</i> is set to <b>true</b>, then the phase <i>p</i> is not added to the collection,
	 * but is taken account in the equivalent phase.
	 * @param p the Material to add as a new phase
	 * @param merge if <b>true</b>, checks if a phase with equal properties already exists, otherwise add without any check
	 */
	void addPhase(Material p, bool merge = true) ;
	
	/** \brief Clears the phase collection */
	void clearPhase() ;
	
	/** \brief Acces to the <i>ith</i> phase, or returns an empty Material if <i>i</i> is out of bound
	 * @param i the index to reach
	 * @return the <i>ith</i> phase if possible
	 */ 
	Material getPhase(int i) const ;

	/** \brief Ensures that every phase in the current Material is defined only once.
	 * If two phases are detected as equals (ignoring extensible Properties like volume or mass), then the
	 * two phases are merged, and any extensible Properties are added together
	 */ 
	void mergePhase() ;
	
	/** \brief Removes the <i>ith</i> phase in the collection, if <i>i</i> is a valid index
	 * @param i the index to remove
	 */
	void removePhase(int i) ;

	/** \brief Convenience method to reach the number of phases in the current Material
	 * @return the number of phases
	 */ 
	int sizePhase() const ;

	/** \brief Convenience method to access the value of the Properties of type <i>t</i> in the <i>ith</i> phase.
	 * This method returns 0 if <i>i</i> is out of bounds or if the <i>ith</i> phase does not have a Properties of type <i>t</i>
	 * @param i the index of the phase to reach
	 * @param t the Properties type to get
	 * @return the value of the first Properties of type <i>t</i> in the <i>ith</i> phase, or 0 if not possible
	 */
	double valPhase(int i, PType t) const ;
	
	/** \brief Returns the values of the first Properties item of type <i>t</i> in all phases contained in the Material.
	 * The index of the result vector corresponds to the index in the current Properties collection.
	 * If the desired Properties is not found in a phase, then the corresponding returned value is 0
	 * @param t the Properties type to get
	 * @return all values of <i>t</i> in the children phases
	 */
	std::vector<double> valPhase(PType t) const ;

// OPERATIONS ON OTHER MATERIAL

	/** \brief Checks if two materials are equals. Two materials are equals if all their Properties are equals.
	 * @param m the Material to be compared with
	 * @param extensible if <b>true</b>, then equalities of extensible Properties (like Volume or Mass) will not be checked
	 * @return <b>false</b> is at least one Properties in the material are not equals, <b>true</b> otherwise
	 */
	bool equals(Material m, bool extensible = false) const ;

// HOMOGENIZATION

	/** \brief Applies an homogenization scheme to the material.
	 * To be applied, a scheme first needs to extract Properties from the phases of the current Material.
	 * <br>
	 * If the Material does not have enough phases, then nothing is done. Otherwise, the scheme will find all
	 * required Properties in the different phase. If a Properties is missing within a phase, then it will try to build
	 * it from existing informations.
	 * <br>
	 * If all required Properties are found or built, then the scheme is applied. Otherwise nothing is done.
	 * After the scheme calculations, new Properties corresponding to the scheme results are built within the current Material.
	 * <br>
	 * The scheme does not check if the resulting Properties already exist within the current Material.
	 * Any existing Properties will be replaced by the result of the scheme.
	 * 
	 * @param s the scheme to apply
	 *
	 * @return <b>true</b> if the scheme was applied successfully, <b>false</b> otherwise
	 */
	bool homogenize(HomogenizationScheme s) ;

// STATIC METHODS

	/** \brief Convenience method to build a Cauchy-Green matrix.
	 * @param prop a pair of values representing either the Young Modulus and the Poisson Ratio, or the Bulk Modulus and the Shear Modulus
	 * @param hooke indicates whether to consider <i>prop</i> as the young Modulus and Poisson Ratio (if <b>true</b>), 
	 *		or the Bulk and Shear Moduli (if <b>false</b>)
	 * @param dim indicates whether to build the matrix for a 2D or a 3D problem
	 * @return a 3x3 (for 2D) or 6x6 (for 3D) matrix
	 */
	static Matrix cauchyGreen(std::pair<double,double> prop, bool hooke, SpaceDimensionality dim, planeType pt = PLANE_STRESS) ;
	static Matrix orthothropicCauchyGreen(double E_1, double E_2, double G,  double nu, planeType pt= PLANE_STRESS) ;
	static Matrix orthothropicCauchyGreen(double E_1, double E_2, double E_3, double G_1, double G_2, double G_3,  double nu) ;

} ;

} ;


#endif
