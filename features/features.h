// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2014
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2009-2014
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
#include "../mesher/parallel_delaunay_3d.h"
#include "../mesher/structuredmesh.h"
#include "../utilities/samplingcriterion.h"
#include "../utilities/grid.h"
#include "../physics/void_form.h"
#include "../physics/viscoelasticity.h"
#include "../physics/damagemodels/damagemodel.h"
#include "../solvers/assembly.h"
#include "boundarycondition.h"
#include "feature_base.h"

#include <valarray>
#include <deque>
#include <iostream>
#include <algorithm>



namespace Amie
{

typedef enum {
    SAMPLE_RESTRICT_4,
    SAMPLE_RESTRICT_8,
    SAMPLE_RESTRICT_16,
    SAMPLE_NO_RESTRICTION
} SamplingRestrictionType ;

class ConfigTreeItem ;
class EnrichmentManager ;

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
    typedef enum {
        SPACETIME_SWITCH,
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

    struct State {
        bool spaceTimeSwitched ;
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


        State ( FeatureTree * ft ) : ft ( ft ) {
            spaceTimeSwitched = false ;
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
        void setStateTo ( StateType s,bool stepChanged );
    } ;
//

    State state ;
protected:

    std::vector<Vector> reportValues ; 
    SamplingRestrictionType samplingRestriction ;

    std::vector<std::vector<double>> cachedVolumes ;
    std::vector<Point *> extraPoints ;
    std::vector<Point *> nodes ;

    std::vector< BoundaryCondition * > boundaryCondition ;
    /** \brief Contains all the features. */
    std::vector<Feature *> tree ;
    std::vector<EnrichmentManager *> enrichmentManagers ;
    std::vector<Feature *> refinedFeatures ;

    /** \brief For fast Access*/
    Grid * grid ;
    Grid3D * grid3d ;
    double initialValue ;

    /** \brief  Contains the mesh in the form of a delaunay tree.
     * The mesh is generated with linear triangles, and when it is final, midpoints are added and
     * projected. No operations should add midpoints before meshing is complete.
     */
    Mesh<DelaunayTriangle, DelaunayTreeItem> * dtree ;
    Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * dtree3D ;

    std::vector<const Geometry * > refinementZones ;
    std::vector<const Geometry *>  domains ;
    void quadTreeRefine ( const Geometry * location = nullptr ) ;

    std::map<int, Mesh<DelaunayTriangle, DelaunayTreeItem> *> layer2d ;

    TetrahedralElement *father3D  ;
    TriElement *father2D  ;


    double now ;
    double deltaTime ;
    double previousDeltaTime ;
    double realDeltaTime ;
    double minDeltaTime ;

    size_t numdofs ;
    size_t samplingNumber ;
    size_t previousSamplingNumber ;
    size_t maxitPerStep ;
    size_t lastNodeId ;
    size_t lastEnrichmentId ;

    bool renumbered ;
    bool needMeshing ;
    bool reuseDisplacements ;
    bool behaviourChange ;
    bool behaviourSet ;
    bool solverConvergence ;
    bool setBehaviours ;
    bool enrichmentChange ;
    bool stateConverged ;
    bool damageConverged ;

    bool elastic ;
    bool projectOnBoundaries ;

    size_t correctionSteps ;
    bool computeIntersections ;

    /** \brief  List of points used for the mesh.
     * Each point is associated with the feature from whose discretiation it was generated.
     */
    std::deque<std::pair<Point *, const Feature *> > meshPoints;
    std::vector<Point *> additionalPoints ;
    std::map<Feature *, double> samplingFactors ;
    std::map<int, double > scalingFactors ;

    /** \brief  Assembly used for the generation of the stiffness matrix and the solving of the problem.
     */
    Assembly * K ;

    //Assembly2D
    /** \brief  Order to which Elements should be brought once all operations are accomplished.
     */
    Order elemOrder ;
    bool structuredMesh;

    /** \brief  Insert triangle midpoints. */
    void makeToOrder() ;

    /** \brief  Project all points on their respective boundaries.*/
    void stitch() ;

    /** \brief  Check element order and behaviour for space-time calculations.*/
    void checkSpaceTimeConsistency() ;

    /** \brief Project all points on their respective boundaries (triangles only)
     * @param edge number of additionnal points on the edge of the triangles
     * @param time number of additionnal time planes in the elements*/
    void projectTrianglesOnBoundaries ( size_t edge, size_t time ) ;

    /** \brief Project all points on their respective boundaries (tetrahedrons only)
     * @param edge number of additionnal points on the edge of the triangles
     * @param time number of additionnal time planes in the elements*/
    void projectTetrahedronsOnBoundaries ( size_t edge, size_t time ) ;

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
    bool stepElements() ;
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
    void setElementBehavioursFromMesh ( MTYPE * source, MTYPE * destination ) {
        std::cout << "setting element behaviours from previous mesh... " << std::flush ;
// // 		std::vector<ETYPE * > elems = destination->getElements() ;
// #pragma omp parallel for schedule(runtime)
        for ( auto i = source->begin() ; i != source->end()  ; i++ ) {
// 			std::cout << "\r element " << i << "/" << elems.size() << std::flush ;s
// 			Circle c(elems[i]->getRadius(), elems[i]->getCircumCenter()) ;
// 			Sphere s(elems[i]->getRadius(), elems[i]->getCircumCenter()) ;
            std::vector<ETYPE * > conflicts =  source->getConflictingElements ( & ( i->getCircumCenter() ) );
// 			if(elems[i]->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
// 				conflicts = source->getConflictingElements(&c) ;
// 			else
// 				conflicts = source->getConflictingElements(&s) ;

            ETYPE * main = nullptr ;
            if ( !conflicts.empty() ) {
                main = conflicts.front() ;
            } else {
                conflicts = source->getConflictingElements ( i->getPrimitive() ) ;
            }

            if ( conflicts.size() == 1 ) {
                std::vector<ETYPE *> neighbourhood = source->getNeighbourhood ( conflicts[0] ) ;
                for ( auto & n : neighbourhood ) {
                    conflicts.push_back ( n ) ;
                }
            }

            std::vector<double> fractions ;
            for ( size_t j = 0 ; j < conflicts.size() ; j++ ) {
                fractions.push_back ( conflicts[j]->overlapFraction ( i->getPrimitive() ) ) ;
            }

            if ( !conflicts.empty() ) {
                main = conflicts.front() ;

                double overlap = fractions[0] ;
                for ( size_t j = 1 ; j < conflicts.size() ; j++ ) {
                    if ( fractions[j] > overlap ) {
                        main = conflicts[j] ;
                        overlap = fractions[j] ;
                    }
                }

                i->setBehaviour ( destination,main->getBehaviour()->getCopy() ) ;

                if ( i->getBehaviour()->getDamageModel() ) {
                    i->getBehaviour()->getDamageModel()->getState ( true ) = 0 ;
                    double renorm ( 0 ) ;
                    for ( size_t j = 0 ; j < conflicts.size() ; j++ ) {
                        if ( conflicts[j]->getBehaviour() &&
                                conflicts[j]->getBehaviour()->getDamageModel() &&
                                conflicts[j]->getBehaviour()->getDamageModel()->getState().size() == i->getBehaviour()->getDamageModel()->getState().size() ) {
                            i->getBehaviour()->getDamageModel()->getState ( true ) += conflicts[j]->getBehaviour()->getDamageModel()->getState() *fractions[j] ;
                            renorm += fractions[j] ;
                        }
                    }
                    if ( std::abs ( renorm ) > POINT_TOLERANCE ) {
                        i->getBehaviour()->getDamageModel()->getState ( true ) /= renorm ;
                    }
                }
            } else {
                i->setBehaviour ( nullptr,nullptr ) ;
            }
        }

        std::cout << " ...done."<< std::endl ;
    } 

public:
    /** \brief return the 2D mesh*/
    Mesh<DelaunayTriangle, DelaunayTreeItem> * get2DMesh ( int grid = -1 ) ;

    /** \brief return the 3D mesh*/
    Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * get3DMesh ( ) ;

    std::map< BoundingBoxPosition, std::pair< Segment *, unsigned int> > boundaryCache ;

    void shuffleMeshPoints() ;

    void setInitialValue ( double v ) {
        initialValue = v ;
    } ;
    double getInitialValue() const {
        return initialValue ;
    } ;

    void setPartition ( size_t partitionNumber ) ;
    const std::vector<const Geometry *> getPartition() {
        return domains ;
    } ;
public:
    Vector instants ;

    void setDiscretizationParameters ( ConfigTreeItem * config, ConfigTreeItem * def = nullptr )  ;
    Vector setSteppingParameters ( ConfigTreeItem * config, ConfigTreeItem * def = nullptr )  ;

    State & getState() {
        return state ;
    }
    const State & getState() const {
        return state ;
    }

    void setElastic ( bool e ) {
        elastic = e ;
    }
    bool getElastic() const {
        return elastic ;
    }

    void setProjectionOnBoundaries ( bool p ) {
        projectOnBoundaries = p ;
    }
    bool getProjectionOnBoundaries() const {
        return projectOnBoundaries ;
    }

    const std::vector<Feature *> & getFeatures() const {
        return tree ;
    }

    void setElementGenerationMethod ( size_t correctionSteps_ = 0, bool computeIntersections_ = true ) {
        correctionSteps = correctionSteps_ ;
        computeIntersections = computeIntersections_ ;
    }

    void print() const;
    void printForFeature ( const Feature *f ) const;
    void reMesh() ;

    void addRefinementZone ( const Geometry * geo ) {
        refinementZones.push_back ( geo );
    }

    bool hasLayers() const {
        return layer2d.size() ;
    }

    void setSamplingRestriction ( SamplingRestrictionType sr ) {
        samplingRestriction = sr ;
    }

    std::vector<int> listLayers() const ;

    Point * checkElement ( const DelaunayTetrahedron * t ) const;
    Point * checkElement ( const DelaunayTriangle * t ) const ;
    Feature * getFeatForTetra ( const DelaunayTetrahedron * t ) const;

    void setSamplingFactor ( Feature * f, double a ) {
        samplingFactors[f] = a ;
    }


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
    double averageDamage;
    bool foundCheckPoint ;

public:

    void addBoundaryCondition ( BoundaryCondition * bc ) ;
    void removeBoundaryCondition ( BoundaryCondition * bc ) ;
    void resetBoundaryConditions() ;
    void scaleBoundaryConditions ( double scale ) ;

public:

    std::vector<Vector> intermediateStates ;

    double damageAreaInAggregates ( Mesh<DelaunayTriangle, DelaunayTreeItem>::iterator begin, Mesh<DelaunayTriangle, DelaunayTreeItem>::iterator end ) {
        double total = 0 ;
        double dam = 0 ;
        for ( auto i = begin ; i != end ; i++ ) {
            if ( ( dynamic_cast<Viscoelasticity *> ( i->getBehaviour() ) )->model == PURE_ELASTICITY ) {
                total += i->area() ;
                if ( i->getBehaviour()->getDamageModel() != nullptr ) {
                    dam += i->area() * i->getBehaviour()->getDamageModel()->getState().max() ;
                }
            }
        }
        if ( total < POINT_TOLERANCE ) {
            return 0 ;
        }
        return dam/total ;
    }

    double damageAreaInPaste ( Mesh<DelaunayTriangle, DelaunayTreeItem>::iterator begin, Mesh<DelaunayTriangle, DelaunayTreeItem>::iterator end ) {
        double total = 0 ;
        double dam = 0 ;
        for ( auto i = begin ; i != end ; i++ ) {
            if ( ( dynamic_cast<Viscoelasticity *> ( i->getBehaviour() ) )->model == GENERALIZED_KELVIN_VOIGT ) {
                total += i->area() ;
                if ( i->getBehaviour()->getDamageModel() != nullptr ) {
                    dam += i->area() * i->getBehaviour()->getDamageModel()->getState().max() ;
                }
            }
        }
        if ( total < POINT_TOLERANCE ) {
            return 0 ;
        }
        return dam/total ;
    }

    /** \brief  Initialise the feature tree with the parent feature.
     *
     * @param first Parent of all features. This shoud typically be the sample itself.
     * @return
     */
    FeatureTree ( Feature *first, int layer = -1, double fraction = 1,  size_t gridsize = 100 ) ;
    
    /** \brief construct a featuretree from a µic-output
     * There are 2 modes: if times is empty, a single file is assumed. Otherwise, as time moves forward, successive meshes are loaded as the hydration advances.
     * In the first case, the voxelSource is a voxel file, in the second it is the directory where the µic pixel files are located.
     * 
     */
    FeatureTree ( const char * voxelSource, std::map<unsigned char, Form *> behaviourMap, const std::vector<double> & times = std::vector<double>()) ;
    virtual ~FeatureTree() ;

    void dumpFeatures ( std::string filename ) ;
    void importFeatures ( std::string filename ) ;

    /** \brief  Add feature as the daughter of another.
     *
     * A feature being the daughter of another typically implies that it is fully contained therein.
     * This is however not necessarilly the case.
     *
     * @param father Parent feature.
     * @param t daughter feature.
     */
    void addFeature ( Feature * const father, Feature * const f, int layer = -1, double fraction = 1 ) ;
    void addFeature ( Feature * const father, EnrichmentManager * fm, int layer = -1, double fraction = 1 ) ;

    void removeFeature ( Feature * f ) ;

    void addPoint ( Point * p ) ;

    void twineFeature ( CompositeFeature * const father, CompositeFeature * const f ) ;

    void homothety ( double before, double now, double after ) ;

    /** \brief  Attempt to enhance the mesh, based on a sampling citerion.
     *
     * Criterions are typically made to check the adequacy of the geometry of the triangles generated
     * during the triangulation phase. They also provide a hint as to the placement of new points which could
     * enhance the mesh quality.
     *
     * @param nit Number of refinement iterations to perform.
     * @param cri Citerion to use.
     */
    void refine ( size_t nit, SamplingCriterion *cri ) ;


    double getCurrentTime() const {
        return now ;
    }

    const Vector & getDisplacements ( int g = -1, bool stepTree = false )  ;

    Vector getDisplacements ( Point * pt,  int g = -1, bool stepTree = false )  ;

    /** \brief  Refine the mesh around the features.
     *
     * Features provide a set of geometries which are targets for successive refinement. Refinement is
     * done by inserting a point in the center of each triangle in the zone.
     *
     */
    void refine ( size_t level ) ;

    double getResidualError() const {
        return residualError ;
    }

    /** \brief  Return true if the physics is modified by the BC
     *
     * For example, if damage will occur, cracks will grow, etc.
     *
     */
    bool isStable() ;

    /** \brief  Set the constitutive law of the given feature.
     *
     * <b>Beware!</b> this function <i>deletes</i> the previous tensor. Make sure is is not in use by another feature.
     *
     * @param m pointer to the tensor.
     * @param f Feature to affect.
     */
    void setStrainTensor ( Matrix * m, Feature * const f ) ;

    /** \brief  set the target order of Elements
     *
     * @param ord order of elements to use.
     *
     * \todo Make it work for values other than 1 and 2.
     */
    void setOrder ( Order ord ) ;

    Order getOrder() const {
        return elemOrder  ;
    }

    void setDeltaTime ( double d ) ;
    void setMinDeltaTime ( double d ) {
        minDeltaTime = d ;
    }

    void moveFirstTimePlanes ( double d, std::vector<DelaunayTriangle *> & triangles ) ;
    void moveFirstTimePlanes ( double d, std::vector<DelaunayTetrahedron *> & tets )  ;

    void moveFirstTimePlanes ( double d, const Mesh<DelaunayTetrahedron, DelaunayTreeItem3D >::iterator & begin,  const Mesh<DelaunayTetrahedron, DelaunayTreeItem3D>::iterator & end ) ;
    void moveFirstTimePlanes ( double d, const Mesh<DelaunayTriangle, DelaunayTreeItem>::iterator & begin,  const Mesh<DelaunayTriangle, DelaunayTreeItem>::iterator & end ) ;

    void setMaxIterationsPerStep ( size_t its ) {
        maxitPerStep = its ;
    }

    /** \brief  set Sampling parameter
    */
    void setSamplingNumber ( size_t news );

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
    std::pair<Vector , Vector > getStressAndStrain ( Amie::Mesh< Amie::DelaunayTetrahedron, Amie::DelaunayTreeItem3D >::iterator begin, Amie::Mesh< Amie::DelaunayTetrahedron, Amie::DelaunayTreeItem3D >::iterator end, bool stepTree = true ) ;

    /** \brief Return the stress and strain of a vector of Tetrahedrons*/
    std::pair<Vector , Vector > getGradientAndFlux ( Amie::Mesh< Amie::DelaunayTetrahedron, Amie::DelaunayTreeItem3D >::iterator begin, Amie::Mesh< Amie::DelaunayTetrahedron, Amie::DelaunayTreeItem3D >::iterator end, bool stepTree = true ) ;

    /** \brief Return the temperature/concentration gradient and the flux of a vector of Tetrahedrons*/
    std::pair<Vector , Vector > getGradientAndFlux ( bool stepTree = true ) ;
    std::pair<Vector , Vector > getGradientAndFluxInLayer ( int layer, bool stepTree = true ) ;

    /** \brief Return the stress and strain of the elements of the current mesh*/
    std::pair<Vector , Vector > getStressAndStrain ( bool stepTree = true ) ;
    std::pair<Vector , Vector > getStressAndStrainInLayer ( int layer, bool stepTree = true ) ;
// 	std::pair<Vector , Vector > getStressAndStrainInAllLayers( bool stepTree = true) ;

    Vector getAverageField ( FieldType f, int grid = -1, double t = 0 ) ;
// 	Vector getAverageField( FieldType f, const std::vector<DelaunayTriangle *> & tri) ;
// 	Vector getAverageField( FieldType f, const std::vector<DelaunayTetrahedron *> & tet) ;

    Vector getAverageFieldOnBoundary( BoundingBoxPosition position, FieldType f, int dummy = 0, double t = 0) ;

    std::pair<Vector, Vector> getFieldMinMax ( FieldType f, int grid = -1, double t = 0 ) ;
    std::pair<Vector, Vector> getFieldMinMax ( FieldType f, const std::vector<DelaunayTriangle *> & tri ) ;
    std::pair<Vector, Vector> getFieldMinMax ( FieldType f, const std::vector<DelaunayTetrahedron *> & tet ) ;

    /** \brief Assuming the base sample is a rectangle, this computed the apparent macro strain based on the displacement along the border */
    std::vector<double> getMacroscopicStrain ( const Geometry * base, double tol = 0.001 )  ;
    std::vector<double> getCutMacroscopicStrain ( const Geometry * base, double tol = 0.001, double cut = .99 )  ;
    std::vector<double> getMedianMacroscopicStrain ( const Geometry * base, double tol = 0.001 ) ;

    std::vector<Point *> getNodes() ;

    size_t numPoints() const ;

    /** \brief Step in time
    */
    bool step() ;

    /** \brief Step to next checkpoint.
     * This finds the load such that the sample is at equilibrium after incrementing the damage
    */
    bool stepToCheckPoint() ;

    /** \brief Perform a time step, but do not update the features*/
    void elasticStep() ;

    /** \brief Return the Assembly*/
    Assembly * getAssembly ( bool forceReassembly = true ) ;

    /** \brief return the triangles lying next to a mesh border*/
    std::vector<DelaunayTriangle *> getBoundingTriangles ( const Feature * f = nullptr ) ;

    Form * getElementBehaviour ( Mesh<DelaunayTriangle, DelaunayTreeItem>::iterator & t, int layer = -1,  bool onlyUpdate = false ) const ;

    Form * getElementBehaviour ( Mesh<DelaunayTetrahedron, DelaunayTreeItem3D>::iterator & t, int layer = -1,  bool onlyUpdate = false ) const ;
    /** \brief insert a point in the mesh*/
    void insert ( Point * p ) ;

    /** \brief Return true id the argument lies in the root feature*/
    bool inRoot ( const Point &p ) const ;

    Feature * getFeature ( size_t i ) 
    {
        if(tree.size() > i)
            return tree[i] ;
        
        return nullptr ;

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
    void stepMesh();
    
    void printReport(bool printHeader = true, bool vertical = true) ;
    void printReport(const std::vector<FieldType> & fields, bool vertical = true) ;


} ;

class EnrichmentManager
{
protected :
    std::vector<EnrichmentFeature *> featureSet ;
    bool stable ;
public:
    EnrichmentManager(EnrichmentFeature * first) : stable(true) {featureSet.push_back(first);} ;
    void addFeature(EnrichmentFeature * f) {featureSet.push_back(f);};
    void removeFeature(EnrichmentFeature * f) 
    {
        for(auto fs = featureSet.begin() ; fs != featureSet.end() ; fs++)
        {
            if(*fs == f)
            {
                featureSet.erase(fs) ;
                return ;
            }
        }
    };
    
        /** \brief update enrichment geometry*/
    virtual bool step(double dt, Vector * v, Mesh< DelaunayTriangle, DelaunayTreeItem >* dtree) 
    { 
        bool moved = false ;
        for(auto i : featureSet)
        {
            i->step(dt, v, dtree) ;
            moved = moved || i->moved() ;
        }
        return moved ;
        
    };

        /** \brief update enrichment geometry*/
    virtual bool step(double dt, Vector * v, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * dtree) 
    { 
        
        bool moved = false ;
        for(auto i : featureSet)
        {
            i->step(dt, v, dtree) ;
            moved = moved || i->moved() ;
        }
        return moved ;
    };
        
    virtual std::vector<EnrichmentFeature *> & getFeatures() { return featureSet ;}
    virtual const std::vector<EnrichmentFeature *> & getFeatures() const { return featureSet ;}
    
    bool converged() { return stable ; } ;
} ;


/** \brief functor for usage with STL containers. order pairs of point, features by point location*/
struct PairPointFeatureEqual {
    bool operator() ( std::pair<Point *, Feature *> p1, std::pair<Point *, Feature *> p2 ) {
        return *p1.first == *p2.first ;
    }
} ;

/** \brief functor for usage with STL containers. order pairs of point, features by point location*/
struct PairPointFeatureLess_Than {
    bool operator() ( std::pair<Point *, Feature *> p1, std::pair<Point *, Feature *> p2 ) {
        return *p1.first < *p2.first ;
    }
} ;

/** \brief functor for usage with STL containers. order pairs of point, features by point location*/
struct PairPointFeatureLess_Than_x {
    bool operator() ( std::pair<Point *, Feature *> p1, std::pair<Point *, Feature *> p2 ) {
        return p1.first->getX() < p2.first->getX() ;
    }
} ;

/** \brief functor for usage with STL containers. order pairs of point, features by point location*/
struct PairPointFeatureLess_Than_y {
    bool operator() ( std::pair<Point *, Feature *> p1, std::pair<Point *, Feature *> p2 ) {
        return p1.first->getY() < p2.first->getY() ;
    }
} ;

/** \brief functor for usage with STL containers. order pairs of point, features by point location*/
struct PairPointFeatureLess_Than_z {
    bool operator() ( std::pair<Point *, Feature *> p1, std::pair<Point *, Feature *> p2 ) {
        return p1.first->getZ() < p2.first->getZ() ;
    }
} ;

} 

#endif // __FEATURES_H_
