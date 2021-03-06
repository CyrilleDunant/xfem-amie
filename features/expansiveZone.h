
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef __EXPANSIVE_ZONE_H__
#define __EXPANSIVE_ZONE_H__

#include "enrichmentInclusion.h"

namespace Amie
{

struct StiffnessWithImposedStrain ;

/** \brief Expansive zone.
 *
 * This enrichement feature will introduce a
 * soft discontinuity in the mesh, as well as
 * attribute an imposed-strain elastic behaviour
 * to the material contained within
*/
class ExpansiveZone :  public EnrichmentInclusion
{
    std::set<DelaunayTriangle *> bimateralInterfaced ;
    std::set<DelaunayTriangle *> expansive ;
    bool homogeneized ;

public:

    /** \brief Constructor. construct the zone
    *
    * @param father supporting feature. This features behaviour is used for the bimaterial interfaces
    * @param radius External radius
    * @param x center x
    * @param y center y
    * @param cgTensor Stifness of the expensive material of the ring, This is the complete experssion of the CG stress tensor
    * @param deformation Vector of the imposed strain
    */
    ExpansiveZone(Feature *father, double radius, double x, double y, const Matrix & cgTensor, Vector deformation) ;
    ExpansiveZone(Feature *father, double radius, double x, double y, Form * gel) ;
    ExpansiveZone(Feature *father, double radius, double x, double y) : EnrichmentInclusion( father, radius, x, y ),  homogeneized(false) { isVirtualFeature = true ; } ;
    virtual ~ExpansiveZone() ;

    /** \brief enrich elements and change their Behaviour if required*/
    virtual void enrich(size_t & , Mesh<DelaunayTriangle, DelaunayTreeItem> * dtree);

    /** \brief return empty vector*/
    virtual std::vector<Amie::DelaunayTetrahedron*> getElements3D(FeatureTree*) {
        return std::vector<Amie::DelaunayTetrahedron*>() ;
    }
    virtual std::vector<Amie::DelaunayTriangle*> getElements2D(FeatureTree*) {
        return std::vector<Amie::DelaunayTriangle*>() ;
    }

    virtual void setBehaviour(Form * const b) ;

    virtual void print() const
    {
        std::cout << "I am an expansive zone" << std::endl ;
    }

    void reset() ;

    virtual void addMeshPointsInFather() ;

    bool isHomogeneized() const {
        return homogeneized ;
    }

    void setExpansion(Vector a) ;

public:
    GEO_DERIVED_OBJECT(Circle) ;

} ;

/** \brief Material inclusion
 *
 * XFEM define inclusion with arbitrary behaviour
*/
class MaterialInclusion :  public EnrichmentInclusion
{
    std::set<DelaunayTriangle *> bimateralInterfaced ;
    std::set<DelaunayTriangle *> internal ;
    Form * inclusionBehaviour ;
public:

    /** \brief Constructor. construct the zone
    *
    * @param father supporting feature. This features behaviour is used for the bimaterial interfaces
    * @param radius External radius
    * @param x center x
    * @param y center y
    * @param inclusionBehaviour Inclusion behaviour
    */
    MaterialInclusion(Feature *father, double radius, double x, double y, Form * inclusionBehaviour) ;

    MaterialInclusion(Feature *father, double radius, double x, double y) : EnrichmentInclusion( father, radius, x, y ) { } ;

    virtual ~MaterialInclusion() ;

    /** \brief enrich elements and change their Behaviour if required*/
    virtual void enrich(size_t & , Mesh<DelaunayTriangle, DelaunayTreeItem> * dtree);

    /** \brief return empty vector*/
    virtual std::vector<Amie::DelaunayTetrahedron*> getElements(Mesh<DelaunayTetrahedron, DelaunayTreeItem3D>*) {
        return std::vector<Amie::DelaunayTetrahedron*>() ;
    }
    virtual std::vector<Amie::DelaunayTriangle*> getElements(Mesh<DelaunayTriangle, DelaunayTreeItem>*) {
        return std::vector<Amie::DelaunayTriangle*>() ;
    }

    virtual void setBehaviour(Form * const b) { inclusionBehaviour = b ; Feature::setBehaviour(b) ; }

    virtual void print() const
    {
        std::cout << "I am a material inclusion" << std::endl ;
    }

    void reset() ;

public:
    GEO_DERIVED_OBJECT(Circle) ;

} ;

}

#endif
