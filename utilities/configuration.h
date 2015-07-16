
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//

#ifndef __CONFIGURATION_H_
#define __CONFIGURATION_H_

#include "../features/features.h"
#include "../features/sample.h"
#include "../features/boundarycondition.h"
#include "../features/microstructuregenerator.h"
#include "../geometry/geometry_base.h"
#include "../elements/integrable_entity.h"
#include "../solvers/assembly.h"
#include "../physics/material_laws/material_laws.h"
#include "../physics/damagemodels/damagemodel.h"
#include "../physics/fracturecriteria/fracturecriterion.h"
#include "../physics/finite_difference_viscoelasticity.h"
#include "matrixops.h"
#include "granulo.h"
#include "writer/triangle_writer.h"

#include <string>
#include <vector>


namespace Amie
{

/** Configuration atom for the configuration tree.
 *
 */
class ConfigTreeItem
{
protected:
    /**  Parent configuration atom. */
    ConfigTreeItem * father ;

    /**  Children coniguration atom. */
    std::vector< ConfigTreeItem * > children ;

    /** Label for the configuration entry. Serves also as an identifyier.*/
    std::string label ;

    /** Numeral data container. Should be reserved for end-of-tree item.*/
    double data ;

    /** Character string data container. Should be reserved for end-of-tree item.*/
    std::string str ;

public:
    /** Constructor for trunk item.*/
    ConfigTreeItem() ;
    /** Constructor for intermediate branch item.*/
    ConfigTreeItem(ConfigTreeItem * father, std::string label) ;
    /** Constructor for leaf item with double value.*/
    ConfigTreeItem(ConfigTreeItem * father, std::string label, double d) ;
    /** Constructor for leaf item with string value.*/
    ConfigTreeItem(ConfigTreeItem * father, std::string label, std::string s) ;

    /** Accessor, returns the father of the current item */
    ConfigTreeItem * getFather() const ;

    /** Accessor, returns the root of the current tree */
    ConfigTreeItem * getRoot() const ;

    /** Accessor, returns the first direct child of the current item with the corresponding label */
    ConfigTreeItem * getChild(std::string childLabel) const ;

    /** Accessor, returns all direct children of the current item with the corresponding label */
    std::vector<ConfigTreeItem *> getAllChildren(std::string childLabel) const ;

    /** Accessor, returns all direct children of the current item */
    std::vector<ConfigTreeItem *> getAllChildren() const ;

    /** Accessor, returns the tree item below the current item with the corresponding label.
     * This allows access to the direct children as well as grand-children of any level.
     */
    ConfigTreeItem * getChildFromFullLabel(std::string childLabel) const ;

    /** Accessor, returns the tree item below the current item with the corresponding label.
     * This allows access to the direct children as well as grand-children of any level.
     */
    ConfigTreeItem * getChild(std::vector<std::string> childLabelDecomposed) const ;

    /** Accessor, returns the last direct child.*/
    ConfigTreeItem * getLastChild() const ;

    void removeAllChildren() ;
    void removeChild(std::string childLabel) ;
    void removeChild(std::vector<std::string> childLabelDecomposed) ;
    void removeChildFromFullLabel(std::string childLabel) ;

    /** Adds a child.*/
    void addChild(ConfigTreeItem * c) ;

    void addChildren(std::vector<ConfigTreeItem *> c) ;

    /** Sets the father.*/
    void setFather(ConfigTreeItem * f) ;

    /** Sets the father.*/
    void forceSetFather(ConfigTreeItem * f) ;

    /** Returns true if the current item has the input label.*/
    bool is(std::string label) const ;

    /** Returns true if the item has no father.*/
    bool isTrunk() const ;

    /** Returns true if the item has no child.*/
    bool isLeaf() const ;

    /** Returns true if the current item has at least one direct child with the corresponding label.*/
    bool hasChild( std::string label ) const ;

    /** Returns true if the current item has at least one child in its descendents with the corresponding label.
     * This allows access to the direct children as well as grand-children of any level.
     */
    bool hasChildFromFullLabel( std::string childLabel ) const ;

    /** Returns true if the current item has at least one child in its descendents with the corresponding label.
     * This allows access to the direct children as well as grand-children of any level.
     */
    bool hasChild( std::vector<std::string> childLabelDecomposed ) const ;

    /** Return the label of this entity */
    std::string getLabel() const {
        return label ;
    }

    /** Return the label of this entity with all its ascendency.
      * It is written in the form grandparent.parent.self .
      */
    std::string getFullLabel() const ;

    /** Set the label for this atom. */
    void setLabel(std::string label) ;

    void setData(double d) {
        data = d ;
    }

    void setStringData(std::string s) {
        str = s ;
    }

    /** Returns the numeral data of the current item */
    double getData() const {
        return data ;
    }

    /** Returns the numeral data of the descendent indexed by path, or the default value if no descendent with the appropriate path is found */
    double getData(std::string path, double defaultValue = 0.) const ;

    /** Returns the string data of the current item */
    std::string getStringData() const {
        return str ;
    }

    /** Returns the string data of the descendent indexed by path, or the default value if no descendent with the appropriate path is found */
    std::string getStringData(std::string path, std::string defaultValue = std::string()) const ;

    /** Print the current item.*/
    void print() const ;

    /** Print the current item and all its descendents.*/
    void printTree() const ;

    /** Read a vector of double from a file*/
    Vector readVectorFromFile() const ;

    std::map<std::string, double> makeDefinition() const ;

    void define(ConfigTreeItem * definition, bool over) ;

    void define(ConfigTreeItem * definition) ;

    void definePath(std::string baseDirectory) ;

    /** Translates the current item in a function*/
    Function getFunction() const ;

    /** Translates the current item in a 2D sample*/
    Sample * getSample() ;

    /** Translates the current item in a mechanical behaviour*/
    Form * getBehaviour(SpaceDimensionality dim, bool spaceTime = false) ;

    /** Translates the current item in a fracture criterion*/
    FractureCriterion * getFractureCriterion(bool spaceTime = false) ;

    /** Translates the current item in a damage model*/
    DamageModel * getDamageModel(bool spaceTime = false) ;

    /** Translates the current item in a Cauchy-Green stiffness matrix*/
    Matrix getStiffnessMatrix(SpaceDimensionality dim, planeType pt = PLANE_STRESS) const ;

    /** Translates the current item in a vector of imposed strain*/
    Vector getImposedDeformation(SpaceDimensionality dim) const ;

    /** Translates the current item in a particle size distribution*/
    ParticleSizeDistribution * getParticleSizeDistribution() const ;

    /** Translates the current item in a material law for log-creep behaviour*/
    ExternalMaterialLaw * getExternalMaterialLaw() const ;

    /** Translates the current item in a vector of inclusions, and places them into a FeatureTree object*/
    std::vector<std::vector<Feature *> > getInclusions(FeatureTree * F, std::vector<Feature *> base, std::vector<Geometry *> brothers ) ;

    /** Translates the current item in a boundary condition*/
    BoundaryCondition * getBoundaryCondition(FeatureTree * f) const ;

    /** Checks if data must be extracted and output at the i^th time step*/
    bool isAtTimeStep(int i, int nmax) const ;

    /** Writes average field values taken from f in a determined output file*/
    void writeOutput(FeatureTree * f, int n, int nmax, std::vector<unsigned int> cacheIndex) ;

    /** Writes field values taken from f in a determined triangle file*/
    void exportTriangles(FeatureTree * f, int n, int nmax) ;

    /** Writes field values taken from f in a determined svg triangle file*/
    void exportSvgTriangles(MultiTriangleWriter * trg, FeatureTree * f, int n, int nmax) ;

    /** Reads a template file and modifies it according to a set of rules*/
    ConfigTreeItem * makeTemplate() ;

    InclusionGenerator * getInclusionGenerator() const ;

    /** Translate a string in an element order*/
    static Order translateOrder(std::string order) ;

    /** Translate a string in a field type*/
    static FieldType translateFieldType(std::string field, bool & ok) ;

    static planeType translatePlaneType(std::string type) ;

    static SmoothingFunctionType translateSmoothingFunctionType(std::string type) ;

    /** Translate a string in a field type*/
    static TWFieldType translateTriangleWriterFieldType(std::string field, bool & ok) ;

    /** Translate a string in an inclusion type*/
    static TypeInclusion translateInclusionType(std::string type) ;

    /** Translate a string in a geometry type*/
    static GeometryType translateGeometryType(std::string type) ;

    /** Translate a string in a Lagrange multiplier type*/
    static LagrangeMultiplierType translateLagrangeMultiplierType(std::string type) ;

    static EMLOperation translateEMLOperation(std::string op) ;

    static ViscoelasticFiniteDifferenceIntegration translateViscoelasticFiniteDifferenceIntegration(std::string op) ;

    /** Translate a string in a PSD specification type*/
    static PSDSpecificationType translatePSDSpecificationType( std::string specification ) ;

    /** Translate a string in a bounding box position*/
    static BoundingBoxPosition translateBoundingBoxPosition( std::string position ) ;

    /** Decomposes a path as father.current.child into a vector of strings {father, current, child} */
    static std::vector<std::string> decompose(std::string path) ;

    static Vector readLineAsVector(std::string line, char sep = ',') ;
    static std::vector<double> readLineAsStdVector(std::string line, char sep = ',') ;

#ifdef __WIN32
    void makeWindowsPath() ;
#endif

} ;

}



#endif // __CONFIGURATION_H_
