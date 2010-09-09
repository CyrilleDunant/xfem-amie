#ifndef SCHEME_TEMPLATE_H
#define SCHEME_TEMPLATE_H

#include "homogenization_base.h"

namespace Mu
{

/** \brief A PhaseTemplate is a convenience object which allows formated extraction of data from one or more Materials.
 * 
 * The main feature of a PhaseTemplate is a collection of specific Properties.
 * As opposed to Materials, in which the order of Properties is not relevant, the order of Properties in a PhaseTemplate collection
 * is particularly important.
 * <br>
 * When a PhaseTemplate will analyse a Material, it will returns a vector of values will correspond to the values of the Properties
 * of the Material, in the same order as they appear in the PhaseTemplate.
 * <br>
 * A PhaseTemplate acts as a buffer. When a Material is processed (using the <i>cast()</i> method), the values of the Properties of
 * the PhaseTemplate are set to the values found in the Material.  If these values are not extracted (with the method <i>val()</i>) before
 * the next Material is read, then they are lost, since the values of the later one replace the values of the previous.
 * <br>
 * A PhaseTemplate can process a limited number of Materials. Upon reaching this limit, the PhaseTemplate is "filled" and cannot
 * format Materials anymore. It is possible to set an indefinite number of Materials that can be processed, by setting the <i>max</i> 
 * flag to -1.
 * <br>
 * Static methods are given to create automatically usual PhaseTemplates.
 */
class PhaseTemplate
{
	std::vector<Properties> prop ;
	int max ;
	int done ;
	
public:

// CONSTRUCTORS

	/** \brief Simple constructor, creates an empty PhaseTemplate.
	 * @param n the number of Materials that the PhaseTemplate can process before being filled (-1 for infinite Materials).
	 */
	PhaseTemplate(int n = 1) ;

	/** \brief Simple constructor, creates a PhaseTemplate with pre-selected Properties
	 * @param n the number of Materials that the PhaseTemplate can process before being filled (-1 for infinite Materials).
	 * @param t the Properties type to use in this PhaseTemplate
	 */
	PhaseTemplate(int n, std::vector<PType> t) ;

// OPERATIONS ON PHASE COUNT
	
	/** \brief Checks is a PhaseTemplate is filled or not. A filled PhaseTemplate cannot process Materials anymore.
	 * If <i>max</i> is set to -1, then the PhaseTemplate is never filled.
	 * @return <b>true</b> if the PhaseTemplate is filled, <b>false</b> otherwise
	 */
	bool filled() const ;

	/** \brief Returns the maximum number of phases that the PhaseTemplate can process before being filled
	 *
	 * This method takes into account the number of phases that the PhaseTemplate has already processed.
	 * If the PhaseTemplate can process an infinite number of Materials, then the result is equal to the input <i>count</i>
	 * @param count a number of Materials
	 * @return the number of Materials that the PhaseTemplate can process before being filled.
	 */
	int phases(int count = 0) const ;

// OPERATIONS ON PROPERTIES
	
	/** \brief Returns the <i>ith</i> Properties in the collection, or a P_BAD_INDEX Properties if <i>i</i> is out of bounds,
	 * @param i the index to reach
	 * @return the <i>ith</i> Properties, when possible
	 */
	Properties get(int i) const ;

	/** \brief Checks if the Properties type <i>t</i> exists within the collection
	 * @param t the Properties type to check
	 * @return <b>true</b> if <i>t</i> exists within the collection, <b>false</b> otherwise
	 */
	bool hasProperties(PType t) const ;

	/** \brief Search for the index of the Properties of type <i>t</i> contained in the collection.
	 * @param t the Properties type to check
	 * @return the index of the first Properties item of type <i>t</i> met in the collection, or -1 if not found
	 */
	int indexOfProperties(PType t) const ;	
	
	/** \brief Convenience method to access the size of the Properties collection
	 * @return the number of Properties in the collection
	 */	
	int size() const ;

	/** \brief Returns all the values buffered in the PhaseTemplate.
	 *
	 * The values are ordered in the same order of the Properties in the PhaseTemplate
	 * @return the vector of values
	 */ 
	std::vector<double> val() const ;

// BUFFERING OPERATIONS

	/** \brief Forces the values of the Properties contained in the PhaseTemplate to the values contained in <i>data</i>
	 *
	 * (or 0 is the input vector is too small. This method does not check if the PhaseTemplate is filled or not, and will not 
	 * increment the number of extractions done so far.
	 * @param data the value to set
	 */
	void buffer(std::vector<double> data) ;
	
	/** \brief Extract data contained into a Material and buffers them.
	 *
	 * The PhaseTemplate will try to find the required Properties in the Material.
	 * If a specific Properties is not found, then the PhaseTemplate will try to build it
	 * using the existing information, and the higher-level Material (if possible).
	 * If a Properties can either be found or built, then the buffering operation fails.
	 * <br>
	 * This operation increments the number of buffering done so far.
	 * @param m the Material to buffer
	 * @param father the higher-level Material, used to build unknown Properties if needed
	 * @return <b>true</b> if the buffering operation succeeded, <b>false<b> otherwise
	 */
	bool buffer(Material m, Material father = Material()) ;

	/** \brief Set the <i>ith</i> Properties in the PhaseTemplate to the value of <i>p</i>, if <i>p</i> is a suitable Properties.
	 * @param i the index of the Properties in the PhaseTemplate to reach
	 * @param p the Properties to buffer
	 * @return <b>true</b> if the operation was successful, <b>false</b> otherwise (usually because <i>i</i> is out of bounds
	 * or <i>p</i> is not of the appropriate Properties type.
	 */
	bool buffer(int i, Properties p) ;

// STATIC METHODS

	/** \brief Convenience method to build a PhaseTemplate containing the following Properties (in the following order):
	 * 	P_BULK_MODULUS
	 *	P_SHEAR_MODULUS
	 * @param n the maximum number of phases in the buffer
	 * @return the corresponding PhaseTemplate
	 */
	static PhaseTemplate makeBulkShearTemplate(int n = -1) ;

	/** \brief Convenience method to build a PhaseTemplate containing the following Properties (in the following order):
	 *	P_VOLUME_FRACTION
	 * 	P_BULK_MODULUS
	 *	P_SHEAR_MODULUS
	 * @param n the maximum number of phases in the buffer
	 * @return the corresponding PhaseTemplate
	 */
	static PhaseTemplate makeVolumeBulkShearTemplate(int n = -1) ;

	/** \brief Convenience method to build a PhaseTemplate containing the following Properties (in the following order):
	 * 	P_EXPANSION_COEFFICIENT
	 * @param n the maximum number of phases in the buffer
	 * @return the corresponding PhaseTemplate
	 */
	static PhaseTemplate makeExpansionTemplate(int n = -1) ;

	/** \brief Convenience method to build a PhaseTemplate containing the following Properties (in the following order):
	 * 	P_BULK_MODULUS
	 *	P_EXPANSION_COEFFICIENT
	 * @param n the maximum number of phases in the buffer
	 * @return the corresponding PhaseTemplate
	 */
	static PhaseTemplate makeBulkExpansionTemplate(int n = -1) ;

	/** \brief Convenience method to build a PhaseTemplate containing the following Properties (in the following order):
	 *	P_VOLUME_FRACTION
	 * 	P_BULK_MODULUS
	 *	P_EXPANSION_COEFFICIENT
	 * @param n the maximum number of phases in the buffer
	 * @return the corresponding PhaseTemplate
	 */
	static PhaseTemplate makeVolumeBulkExpansionTemplate(int n = -1) ;

	/** \brief Convenience method to build a PhaseTemplate containing the following Properties (in the following order):
	 *	P_VOLUME_FRACTION
	 * 	P_BULK_MODULUS
	 *	P_SHEAR_MODULUS
	 *	P_EXPANSION_COEFFICIENT
	 * @param n the maximum number of phases in the buffer
	 * @return the corresponding PhaseTemplate
	 */
	static PhaseTemplate makeVolumeBulkShearExpansionTemplate(int n = -1) ;

} ;

/** \brief Base class for homogenization scheme application.
 *
 * A SchemeTemplate role is to extract the data from a Material phases,
 * to make the computation, and to restore the results into the composite material.
 * <br>
 * A SchemeTemplate contains the following information : <ul>
 *	<li>The id of the scheme to be applied</li>
 *	<li>A list of PhaseTemplate which are the input data of the scheme</li>
 *	<li>A single PhaseTemplate which represents the results of the scheme</li></ul>
 * <br>
 */
class SchemeTemplate
{
	HomogenizationScheme scheme ;
	std::vector<PhaseTemplate> phases ;
	PhaseTemplate result ;
	std::vector<std::vector<double> > data ;
	
public:

// CONSTRUCTOR

	/** \brief Simple constructor.
	 *
	 * @param r the PhaseTemplate to use as results
	 * @param p the list of PhaseTemplate to use as input
	 * @param s the scheme to apply
	 */
	SchemeTemplate(PhaseTemplate r, std::vector<PhaseTemplate> p, HomogenizationScheme s = SCHEME_DO_NOTHING) ;

	/** \brief Simple constructor.
	 *
	 * This constructor will automatically build the appropriate input and result PhaseTemplate
	 * @param s the scheme to apply
	 */	
	SchemeTemplate(HomogenizationScheme s) ;

// OPERATIONS ON MATERIALS
	
	/** \brief Applies the scheme to a multi-phase Material
	 *
	 * The Material will first be analyzed in order to extract all relevant data.
	 * The scheme will fail if there is not enough phases, or if there is at least one missing Properties in the Material phases.
	 * Missing Properties will be built when possible during the process.
	 * After extracting the data, they will be sent to the homogenization scheme, which will give back the results as a vector.
	 * This vector will then be formated into Properties and buffered in the result PhaseTemplate
	 *
	 * @param m the composite Material on which to apply the scheme
	 * @return <b>true</b> if the scheme computation runs successfully, <b>false</b> otherwise
	 */
	bool cast(Material m) ;

// OPERATIONS ON RESULTS
	
	/** \brief Returns the number of output Properties of the scheme
	 * @return the number of Properties in the result PhaseTemplate
	 */
	int sizeResult() const ;

	/** \brief Returns a formatted data output.
	 *  
	 * Note that if no homogenization scheme were run successfully before calling this method, it will only returns 0.
	 * @param i the index to reach
	 * @return the <i>ith</i> Properties in the result PhaseTemplate, or a P_BAD_INDEX Properties if <i>i</i> is out of bound.
	 */
	Properties getResult(int i) const ;

// STATIC METHODS

	/** \brief Calls the appropriate homogenization scheme
	 *
	 * Data are formatted using the following rule :<ul>
	 *	<li>The <i>ith</i> row corresponds to the <i>ith</i> phase in the Material</li>
	 *	<li>The <i>jth</i> column on line <i> correspond to the <i>jth</i> Properties
	 *		found in the PhaseTemplate used to build the <i>ith</i> row.</li></ul>
	 * @param d a vector of vector of formatted extracted data.
	 * @return a vector of output data
	 */
	static std::vector<double> applyScheme(HomogenizationScheme s, std::vector<std::vector<double> > d) ;	

// ELASTIC HOMOGENIZATION SCHEMES

	/** \brief Diluted bi-phasic elastic homognization.
	 * First phase is the matrix, second phase are the inclusions.
	 */
	static std::vector<double> elasticityDilutedScheme(std::vector<std::vector<double> > d) ;

	/** \brief Diluted multi-phasic elastic homognization.
	 * First phase is the matrix, subsequent phases are the inclusions
	 */
	static std::vector<double> elasticityGeneralizedDilutedScheme(std::vector<std::vector<double> > d) ;

	/** \brief Incremental bi-phasic elastic homognization.
	 * First phase is the matrix, second phase are the inclusions.
	 */
	static std::vector<double> elasticityIncrementalScheme(std::vector<std::vector<double> > d, double dalpha = 1e-5) ;

	/** \brief Bi-phasic Mori-Tanaka elastic homognization.
	 * First phase is the matrix, second phase are the inclusions.
	 */
	static std::vector<double> elasticityMoriTanakaScheme(std::vector<std::vector<double> > d) ;

	/** \brief Multi-phasic elastic homognization.
	 * First phase is the matrix, subsequent phases are the inclusions
	 */
	static std::vector<double> elasticityGeneralizedMoriTanakaScheme(std::vector<std::vector<double> > d) ;

	/** \brief Bi-phasic self-consistent elastic homogenization 
	 * First and second phases play the same role.
	 */
	static std::vector<double> elasticitySelfConsistentScheme(std::vector<std::vector<double> > d) ;

	/** \brief Multi-phasic self-consistent elastic homogenization 
	 * All phases play the same role.
	 */
	static std::vector<double> elasticityGeneralizedSelfConsistentScheme(std::vector<std::vector<double> > d) ;

// THERMAL EXPANSION HOMOGENIZATION SCHEMES
	
	/** \brief Bi-phasic Hobbs homognization of expansion coefficient.
	 * First phase is the matrix, second phase are the inclusions.
	 */
	static std::vector<double> expansionHobbsScheme(std::vector<std::vector<double> > d) ;

	/** \brief Bi-phasic Kerner homognization of expansion coefficient.
	 * First phase is the matrix, second phase are the inclusions.
	 */
	static std::vector<double> expansionKernerScheme(std::vector<std::vector<double> > d) ;

	/** \brief Bi-phasic Turner homognization of expansion coefficient.
	 * First phase is the matrix, second phase are the inclusions.
	 */
	static std::vector<double> expansionTurnerScheme(std::vector<std::vector<double> > d) ;

} ;
	
	
} ;

#endif
