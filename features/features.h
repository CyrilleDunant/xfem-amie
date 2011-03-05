// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2010
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __FEATURES_H_
#define __FEATURES_H_

#include "../geometry/geometry_base.h"
#include "../elements/elements.h" 
#include "../mesher/delaunay.h"
#include "../mesher/delaunay_3d.h"
#include "../mesher/structuredmesh.h"
#include "../utilities/samplingcriterion.h"
#include "../utilities/grid.h"
#include "../physics/void_form.h"
#include "../physics/damagemodels/damagemodel.h"
#include "../solvers/assembly.h"
#include "boundarycondition.h"
#include "feature_base.h"

#include <valarray>
#include <deque>
#include <iostream>
#include <algorithm>



namespace Mu
{



/** \brief Container for the features defining the setup.
 * 
 * The feature tree is responsible for all global operations: meshing, matrix assembly, 
 * applying boundary conditions.
 * 
 * Postprocessing the results is also done here.
 * 
 */
class FeatureTree 
{
protected:
	typedef enum
	{
		SAMPLED,
		MESHED,
		RENUMBERED,
		BEHAVIOUR_SET,
		STITCHED,
		INITIALISED,
		ENRICHED,
		ASSEMBLED,
		SOLVED,
		BEHAVIOUR_STEPPED,
		XFEM_STEPPED,
		FEATURE_STEPPED
	} StateType ;
	
	struct State
	{
		bool sampled ;
		bool meshed ;
		bool renumbered ;
		bool behaviourSet ;
		bool behaviourUpdated ;
		bool stitched ;
		bool initialised ;
		bool enriched ;
		bool assembled ;
		bool solved ;
		bool behaviourStepped ;
		bool xfemStepped ; 
		bool featureStepped ; 
		
		FeatureTree * ft ;
		
		
		State(FeatureTree * ft) : ft(ft) 
		{
			sampled = false ;
			meshed = false;
			renumbered = false;
			behaviourSet = false;
			behaviourUpdated = false;
			stitched = false;
			initialised = false;
			enriched = false;
			assembled = false;
			solved = false;
			behaviourStepped = false;
			xfemStepped = false; 
			featureStepped = false; 
		} ;
		void setStateTo(StateType s,bool stepChanged );
	} ;
// 

	State state ;
protected:
	
	void duplicate2DMeshPoints() ;
	void duplicate3DMeshPoints() ;
	
	std::vector< BoundaryCondition * > boundaryCondition ;
	/** \brief Contains all the features. */
	std::vector<Feature *> tree ;

	/** \brief For fast Access*/
	Grid * grid ;
	Grid3D * grid3d ;
	
	std::vector<Mesh<DelaunayTriangle, DelaunayTreeItem> *> coarseTrees ;
	std::vector<Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> *> coarseTrees3D ;
	
	/** \brief  Contains the mesh in the form of a delaunay tree. 
	 * The mesh is generated with linear triangles, and when it is final, midpoints are added and 
	 * projected. No operations should add midpoints before meshing is complete.
	 */
	Mesh<DelaunayTriangle, DelaunayTreeItem> * dtree ;
	Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * dtree3D ;
	
	TetrahedralElement *father3D  ;
	TriElement *father2D  ;


	double now ;
	double deltaTime ; 
	
	size_t numdofs ;
	size_t samplingNumber ;
	size_t previousSamplingNumber ;
	size_t maxitPerStep ;
	size_t lastNodeId ;
	std::vector<size_t> coarseLastNodeId ;
	size_t lastEnrichmentId ;
	std::vector<size_t> coarseLastEnrichmentId ;
	
	bool renumbered ;
	bool needAssembly ;
	bool needMeshing ;
	bool reuseDisplacements ;
	bool behaviourChange ;
	bool solverConvergence ;
	bool setBehaviours ;
	bool enrichmentChange ;
	
	size_t correctionSteps ;
	bool computeIntersections ;
	
	/** \brief  List of points used for the mesh.
	 * Each point is associated with the feature from whose discretiation it was generated.
	 */
	std::deque<std::pair<Point *, Feature *> > meshPoints;
	std::vector<Point *> additionalPoints ;
	
	/** \brief  Assembly used for the generation of the stiffness matrix and the solving of the problem.
	 */
	Assembly * K ;
	std::vector<Assembly *> coarseAssemblies ;
	
	//Assembly2D
	/** \brief  Order to which Elements should be brought once all operations are accomplished.
	 */
	Order elemOrder ;
	
	/** \brief  Insert triangle midpoints. */
	void makeToOrder() ;
	
	/** \brief  Project all points on their respective boundaries.*/
	void stitch() ;
	
	void renumber() ;
	void enrich() ;
	/** \brief  Generate the sample points for all the features. The features are passed a sampling 
	 * argument proportionnal to their area compared with the area of the root feature. 
	 * If the number is lower than 10, than the argument passed is 10.
	 * 
	 * @param npoints number of sampling points on the boundary of the sample. This number must be 
	 * divisable by four if the sample is a rectangle.
	 */
	void sample() ;

	void setElementBehaviours() ;
	void updateElementBehaviours() ;
	void solve() ;
	void stepElements() ;
	void stepXfem() ;
	
	/** \brief  Generate the triangulation.
	 * Once the sampling is done, the sampling points are fed into a Delaunay Triangulation algorithm, 
	 * which generates the triangles. The mesh is still composed of 3 points triangles at this point.
	 * 
	 * @param correctionSteps additional steps where points are inserted in incorrect tetrahedrons.
	 */
	void generateElements() ;
	
	/** \brief  Perform the assembly.
	 * 
	 * The mesh is completed with eventual intermediate points (if higher order elements were asked for), and elements 
	 * are assembled in K.
	 * 
	 */
	void assemble() ;
	
	template<class MTYPE, class ETYPE>
	void setElementBehavioursFromMesh( const MTYPE * source, MTYPE * destination)
	{
		std::cout << "setting element behaviours from previous mesh... " << std::flush ;
		std::vector<ETYPE * > elems = destination->getElements() ;
#pragma omp parallel for
		for(size_t i = 0 ; i < elems.size() ; i++)
		{
// 			std::cout << "\r element " << i << "/" << elems.size() << std::flush ;s
// 			Circle c(elems[i]->getRadius(), elems[i]->getCircumCenter()) ;
// 			Sphere s(elems[i]->getRadius(), elems[i]->getCircumCenter()) ;
			std::vector<ETYPE * > conflicts =  source->getConflictingElements(&(elems[i]->getCircumCenter()) );
// 			if(elems[i]->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
// 				conflicts = source->getConflictingElements(&c) ;
// 			else
// 				conflicts = source->getConflictingElements(&s) ;
			
			ETYPE * main = NULL ;
			if(!conflicts.empty())
				main = conflicts.front() ;
			else
				conflicts = source->getConflictingElements(elems[i]->getPrimitive()) ;
			
			if(conflicts.size() == 1)
			{
				for(size_t j = 0 ; j < conflicts[0]->neighbourhood.size() ; j++)
					conflicts.push_back(conflicts[0]->getNeighbourhood(j)) ;
			}
			
			std::vector<double> fractions ; 
			for(size_t j = 0 ; j < conflicts.size() ; j++)
				fractions.push_back(conflicts[j]->overlapFraction(elems[i]->getPrimitive())) ;
			
			if(!conflicts.empty())
			{
				main = conflicts.front() ;
				
				double overlap = fractions[0] ;
				for(size_t j = 1 ; j < conflicts.size() ; j++)
				{
					if(fractions[j] > overlap)
					{
						main = conflicts[j] ;
						overlap = fractions[j] ;
					}
				}
				
				elems[i]->setBehaviour(main->getBehaviour()->getCopy()) ;
				
				if(elems[i]->getBehaviour()->getDamageModel())
				{
					elems[i]->getBehaviour()->getDamageModel()->state = 0 ;
					elems[i]->getBehaviour()->getDamageModel()->previousstate = 0 ;
					double renorm(0) ;
					for(size_t j = 0 ; j < conflicts.size() ; j++)
					{
						if(conflicts[j]->getBehaviour() &&
							conflicts[j]->getBehaviour()->getDamageModel() && 
							conflicts[j]->getBehaviour()->getDamageModel()->state.size() == elems[i]->getBehaviour()->getDamageModel()->state.size())
						{
							elems[i]->getBehaviour()->getDamageModel()->state += conflicts[j]->getBehaviour()->getDamageModel()->state*fractions[j] ;
							elems[i]->getBehaviour()->getDamageModel()->previousstate += conflicts[j]->getBehaviour()->getDamageModel()->previousstate*fractions[j] ;
							renorm += fractions[j] ;
						}
					}
					if(std::abs(renorm) > POINT_TOLERANCE)
					{
						elems[i]->getBehaviour()->getDamageModel()->state /= renorm ;
						elems[i]->getBehaviour()->getDamageModel()->previousstate /= renorm ;
					}
				}
			}
			else
				elems[i]->setBehaviour(NULL) ;
		}
		
		std::cout << " ...done."<< std::endl ;
	} ;
	
/** \brief return the 2D mesh*/
	Mesh<DelaunayTriangle, DelaunayTreeItem> * get2DMesh(int grid = -1) ;

/** \brief return the 3D mesh*/
	Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * get3DMesh(int grid = -1) ;
	
	void shuffleMeshPoints() ;

public:
	

	State & getState() {return state ;}
	const State & getState() const {return state ;}
	const std::vector<Feature *> & getFeatures() const {return tree ;}
	
	void setElementGenerationMethod(size_t correctionSteps_ = 0, bool computeIntersections_ = true) { correctionSteps = correctionSteps_ ; computeIntersections = computeIntersections_ ; }
	
	void print() const;
	void printForFeature(const Feature *f) const;
	void reMesh() ;
	
	Point * checkElement( const DelaunayTetrahedron * t ) const;
	Point * checkElement( const DelaunayTriangle * t ) const ;
	Feature * getFeatForTetra( const DelaunayTetrahedron * t ) const;

	
	std::vector<DelaunayTriangle *> getElements2D(int g = -1) ;
	std::vector<DelaunayTriangle *> getElements2D(const Point *p, int g = -1) ;
	std::vector<DelaunayTriangle *> getElements2D(const Geometry *p, int g = -1) ;
	
	std::vector<DelaunayTetrahedron *> getElements3D(int g = -1) ;
	std::vector<DelaunayTetrahedron *> getElements3D(const Point *p, int g = -1) ;
	std::vector<DelaunayTetrahedron *> getElements3D(const Geometry *p, int g = -1) ;
	
	
	/** \brief The solver converged at the last step.
	* This condition checks whether a numerical solution was found. If false, this means the solver diverged, or 
	* could not find a solution with the prescribed precision in the given maximum number of iteration
	*/
	bool solverConverged() const ;
	
	/** \brief At least an element changed behaviour.
	* This typically occurs when an element is damaged, or a virtual feature changed its geometry.
	*/
	bool behaviourChanged() const ;
	
	/** \brief At least an enrichment feature changed its geometry/behaviour
	* This typically occurs when a crack grows
	*/
	bool enrichmentChanged() const ;

	void forceEnrichmentChange() ;

	double residualError ;
	double crackedVolume ;
	double damagedVolume ;
	bool useMultigrid;
	
public:
	
	void addBoundaryCondition(BoundaryCondition * bc) ;
	void removeBoundaryCondition(BoundaryCondition * bc) ;
	void resetBoundaryConditions() {  boundaryCondition.clear() ; } ;
	
	
public:
	
	/** \brief  Initialise the feature tree with the parent feature.
	 * 
	 * @param first Parent of all features. This shoud typically be the sample itself.
	 * @return 
	 */
	FeatureTree(Feature *first, size_t gridsize = 100) ;
	virtual ~FeatureTree() ;
	
	void dumpFeatures(std::string filename) ;
	void importFeatures(std::string filename) ;
	
	/** \brief  Add feature as the daughter of another.
	 * 
	 * A feature being the daughter of another typically implies that it is fully contained therein. 
	 * This is however not necessarilly the case.
	 * 
	 * @param father Parent feature.
	 * @param t daughter feature.
	 */
	void addFeature(Feature *father, Feature * f) ;
	
	void addPoint(Point * p) ;

	void twineFeature(CompositeFeature * father, CompositeFeature * f) ;
	
	/** \brief  Attempt to enhance the mesh, based on a sampling citerion.
	 * 
	 * Criterions are typically made to check the adequacy of the geometry of the triangles generated 
	 * during the triangulation phase. They also provide a hint as to the placement of new points which could 
	 * enhance the mesh quality.
	 * 
	 * @param nit Number of refinement iterations to perform.
	 * @param cri Citerion to use.
	 */
	void refine(size_t nit, SamplingCriterion *cri) ;
	
	
	const Vector & getDisplacements(int g = -1)  ;
	
	
	/** \brief  Refine the mesh around the features.
	 * 
	 * Features provide a set of geometries which are targets for successive refinement. Refinement is 
	 * done by inserting a point in the center of each triangle in the zone.
	 * 
	 */
	void refine(size_t level) ;

	double getResidualError() const {return residualError ;}
	
	/** \brief  Return true if the physics is modified by the BC
	 * 
	 * For example, if damage will occur, cracks will grow, etc.
	 * 
	 */
	bool stable(double dt) ;

	/** \brief  Set the constitutive law of the given feature.
	 * 
	 * <b>Beware!</b> this function <i>deletes</i> the previous tensor. Make sure is is not in use by another feature. 
	 * 
	 * @param m pointer to the tensor. 
	 * @param f Feature to affect.
	 */
	void setStrainTensor(Matrix * m, Feature * f) ;
	
	/** \brief  set the target order of Elements
	 * 
	 * @param ord order of elements to use.
	 * 
	 * \todo Make it work for values other than 1 and 2.
	 */
	void setOrder(Order ord) ;
	
	void setDeltaTime(double d) {deltaTime = d ;}
	
	void setMaxIterationsPerStep(size_t its) {maxitPerStep = its ;}
	
		/** \brief  set Sampling parameter
	 */
	void setSamplingNumber(size_t news);
	
	/**  \brief  Postprocess the result.
	 * 
	 * Given a vector containing the displacements at each point (containing n times as many elements as 
	 * there are points, n being the number of degrees of liberty) it returns an array containing the 
	 * strain values at the mesh points.
	 * 
	 * @return strain.
	 */
	Vector strainFromDisplacements()  ;
	
	/**  \brief  Postprocess the result.
	 * 
	 * Given a vector containing the displacements at each point (containing n times as many elements as 
	 * there are points, n being the number of degrees of liberty) it returns an array containing the 
	 * stress values at the mesh points.
	 * @return stress.
	 */
	Vector stressFromDisplacements()  ;
	
/** \brief Return the stress and strain of a vector of Tetrahedrons*/
	std::pair<Vector , Vector > getStressAndStrain(const std::vector<DelaunayTetrahedron *> &) ;

/** \brief Return the stress and strain of a vector of Tetrahedrons*/
	std::pair<Vector , Vector > getGradientAndFlux(const std::vector<DelaunayTetrahedron *> &) ;

/** \brief Return the temperature/concentration gradient and the flux of a vector of Tetrahedrons*/
	std::pair<Vector , Vector > getGradientAndFlux(int grid = -1) ;

/** \brief Return the stress and strain of the elements of the current mesh*/
	std::pair<Vector , Vector > getStressAndStrain(int grid = -1) ;
	
	size_t numPoints() const ;
	
/** \brief Step in time
*/
	bool step() ;

/** \brief annul the last timestep*/
	void stepBack() ;

/** \brief Perform a time step, but do not update the features*/
	void elasticStep() ;
	
/** \brief Return the Assembly*/
	Assembly * getAssembly() ;

/** \brief return the triangles of the mesh*/
// 	std::vector<DelaunayTriangle *> getTriangles(int grid = -1);

/** \brief return the tetrahedrons of the mesh*/
// 	std::vector<DelaunayTetrahedron *> getTetrahedrons(int grid = -1) ;
		
/** \brief return the triangles lying next to a mesh border*/
	std::vector<DelaunayTriangle *> getBoundingTriangles(Feature * f = NULL) ;	
	
/** \brief return the Behaviour of the argument, deduced from the Feature s*/
	Form * getElementBehaviour(const DelaunayTriangle *t, bool onlyUpdate = false) const ;

/** \brief return the Behaviour of the argument, deduced from the Feature s*/
	Form * getElementBehaviour(const DelaunayTetrahedron *t, bool onlyUpdate = false) const ;
	
/** \brief insert a point in the mesh*/
	void insert(Point * p ) ;
	
/** \brief Return true id the argument lies in the root feature*/
	bool inRoot(const Point &p) const ;
	
	Feature * getFeature(size_t i)
	{
		return tree[i] ;
	}
	
/** \brief initialise the element states*/
	void initializeElements() ;
	
	double getMaximumDisplacement()  ;
	double getMinimumDisplacement()  ;
	
/** \brief return true if the currently defined problem is 3D*/
	bool is3D() const ;

/** \brief return true if the currently defined problem is 3D*/
	bool is2D() const ;

	std::vector<DelaunayTriangle> getSnapshot2D() const ;

} ;


/** \brief functor for usage with STL containers. order pairs of point, features by point location*/
struct PairPointFeatureEqual
{
	bool operator()(std::pair<Point *, Feature *> p1, std::pair<Point *, Feature *> p2)
	{
		return *p1.first == *p2.first ;
	}
} ;

/** \brief functor for usage with STL containers. order pairs of point, features by point location*/
struct PairPointFeatureLess_Than
{
	bool operator()(std::pair<Point *, Feature *> p1, std::pair<Point *, Feature *> p2)
	{
		return *p1.first < *p2.first ;
	}
} ;

/** \brief functor for usage with STL containers. order pairs of point, features by point location*/
struct PairPointFeatureLess_Than_x
{
	bool operator()(std::pair<Point *, Feature *> p1, std::pair<Point *, Feature *> p2)
	{
		return p1.first->x < p2.first->x ;
	}
} ;

/** \brief functor for usage with STL containers. order pairs of point, features by point location*/
struct PairPointFeatureLess_Than_y
{
	bool operator()(std::pair<Point *, Feature *> p1, std::pair<Point *, Feature *> p2)
	{
		return p1.first->y < p2.first->y ;
	}
} ;

/** \brief functor for usage with STL containers. order pairs of point, features by point location*/
struct PairPointFeatureLess_Than_z
{
	bool operator()(std::pair<Point *, Feature *> p1, std::pair<Point *, Feature *> p2)
	{
		return p1.first->z < p2.first->z ;
	}
} ;

} ;

#endif // __FEATURES_H_
