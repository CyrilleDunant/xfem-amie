// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2009-2011
//
// Copyright: See COPYING file that comes with this distribution
//


#ifndef GENERALIZED_SPACETIME_VISCOELASTIC_ELEMENT_STATE
#define GENERALIZED_SPACETIME_VISCOELASTIC_ELEMENT_STATE

#include "elements.h"
#include "integrable_entity.h"

namespace Amie
{


class GeneralizedSpaceTimeViscoElasticElementState : public ElementState
{
    Vector averagestressbefore ;
    Vector averagestressafter ;

    Vector averagestrainbefore ;
    Vector averagestrainafter ;

    Vector averagestrainratebefore ;
    Vector averagestrainrateafter ;

    Vector genStrainAtGaussPointBefore ;
    Vector genStrainAtGaussPointAfter ;
    Vector genStrainRateAtGaussPointBefore ;
    Vector genStrainRateAtGaussPointAfter ;

    Vector strainBuffer ;
    Vector principalBuffer ;
    Vector generalizedBuffer ;
    Vector generalizedBufferSecond ;

private:
   void getCache( FieldType f, GaussPointArray & gp, Vector & ret, VirtualMachine * vm ) ;

public:

    void getEssentialAverageFields(FieldType f, Vector & stress, Vector & strain, Vector & strain_rate, VirtualMachine * vm, double t) ;
    GeneralizedSpaceTimeViscoElasticElementState(IntegrableEntity * e, int blocks) ;
    GeneralizedSpaceTimeViscoElasticElementState( Amie::GeneralizedSpaceTimeViscoElasticElementState& s ) ;
    GeneralizedSpaceTimeViscoElasticElementState & operator =( GeneralizedSpaceTimeViscoElasticElementState & s) ;

    virtual void getField( FieldType f, const Point & p, Vector & ret, bool local, VirtualMachine * vm = nullptr, int i = 0)  final;

    virtual void getFieldAtNodes( FieldType f, Vector & ret, VirtualMachine * vm = nullptr, int i = 0) final;

    virtual void getField( FieldType f1, FieldType f2, const Point & p, Vector & ret1, Vector & ret2, bool local, VirtualMachine * vm = nullptr, int i = 0, int j = 0)  final ;

    virtual void getFieldAtNodes( FieldType f1, FieldType f2, Vector & ret1, Vector & ret2, VirtualMachine * vm = nullptr, int i = 0, int j = 0) final;

    virtual double getAverageField( FieldType f, Vector & ret, VirtualMachine * vm = nullptr, double t = 1, const std::vector< double > & weights = std::vector<double>(), int i = 0) final ;

    virtual double getAverageField( FieldType f, FieldType f_, Vector & ret, Vector & ret_, VirtualMachine * vm = nullptr, double t = 1, const std::vector< double > & weights = std::vector<double>(), int i = 0) final ;

    virtual void step( double dt, const Vector *d ) final;

    Vector getCachedFieldAtGaussPointBefore( FieldType f, GaussPointArray & gp, size_t i, VirtualMachine * vm) ;
    Vector getCachedFieldAtGaussPointAfter( FieldType f, GaussPointArray & gp, size_t g, VirtualMachine * vm) ;


    static GaussPointArray genEquivalentGaussPointArray2D ( TriElement * trg, double time ) ;

    static GaussPointArray genEquivalentGaussPointArray3D ( TetrahedralElement * tet, double time ) ;

    static GaussPointArray genEquivalentGaussPointArray ( IntegrableEntity * e, double time )
    {
        if(e->spaceDimensions() == SPACE_TWO_DIMENSIONAL && dynamic_cast<TriElement *>(e))
            return GeneralizedSpaceTimeViscoElasticElementState::genEquivalentGaussPointArray2D( dynamic_cast<TriElement *>(e), time) ;
        if(e->spaceDimensions() == SPACE_THREE_DIMENSIONAL && dynamic_cast<TetrahedralElement *>(e))
            return GeneralizedSpaceTimeViscoElasticElementState::genEquivalentGaussPointArray3D( dynamic_cast<TetrahedralElement *>(e), time) ;
        return e->getGaussPoints() ;
    }

} ;

class GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables final: public GeneralizedSpaceTimeViscoElasticElementState
{
    std::map<std::string, double> variables ;

public:
    GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables(IntegrableEntity * e, std::map<std::string, double> & external, int blocks) ;
    GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables( Amie::GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables& s ) ;
    GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & operator =( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s) ;

    bool has(std::string v) const ;
    double get(std::string v, std::map<std::string, double> & defaultValues) ;
    void synchronize(std::map<std::string, double> & defaultValues) ;
    void set(std::string v, double d) ;
    void add(std::string v, double d) ;
    void multiply(std::string v, double d) ;
    std::map<std::string, double> getVariables() const {
        return variables ;
    }


} ;


}


#endif // GENERALIZED_SPACETIME_VISCOELASTIC_ELEMENT_STATE
