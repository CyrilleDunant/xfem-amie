/* this is an auto-generated file created on 10/3/2017 at 10:24  */

#ifndef __ENUMERATION_TRANSLATOR_H__
#define __ENUMERATION_TRANSLATOR_H__

#include "../elements/integrable_entity.h"
#include "../features/boundarycondition.h"
#include "../geometry/geometry_base.h"
#include "../physics/finite_difference_viscoelasticity.h"
#include "../physics/materials/csh_behaviour.h"
#include "../physics/fracturecriteria/fracturecriterion.h"
#include "../physics/fracturecriteria/spacetimelimitsurfacefracturecriterion.h"
#include "../physics/fracturecriteria/mcft.h"
#include "../physics/material_laws/material_laws.h"
#include "../physics/viscoelasticity.h"
#include "../physics/damagemodels/damagemodel.h"
#include "../polynomial/vm_token.h"
#include "../polynomial/variable.h"
#include "../polynomial/typeref.h"
#include "../polynomial/vm_function_base.h"
#include "../solvers/assembly.h"
#include "../utilities/configuration.h"
#include "../utilities/writer/voxel_writer.h"
#include "../utilities/writer/triangle_writer.h"
#include "../utilities/parser/config_parser.h"
#include "../utilities/granulo.h"
#include "../utilities/tensor.h"

namespace Amie
{

struct Enum
{

    // parsed from header file: ../elements/integrable_entity.h
    static Order getOrder(std::string type, bool * ok = 0)
    {
        if(ok) { *ok = true ; }
        if( type == "CONSTANT") { return CONSTANT ; }
        if( type == "LINEAR") { return LINEAR ; }
        if( type == "QUADRATIC") { return QUADRATIC ; }
        if( type == "CUBIC") { return CUBIC ; }
        if( type == "QUADRIC") { return QUADRIC ; }
        if( type == "QUINTIC") { return QUINTIC ; }
        if( type == "CONSTANT_TIME_LINEAR") { return CONSTANT_TIME_LINEAR ; }
        if( type == "CONSTANT_TIME_QUADRATIC") { return CONSTANT_TIME_QUADRATIC ; }
        if( type == "LINEAR_TIME_LINEAR") { return LINEAR_TIME_LINEAR ; }
        if( type == "LINEAR_TIME_QUADRATIC") { return LINEAR_TIME_QUADRATIC ; }
        if( type == "QUADRATIC_TIME_LINEAR") { return QUADRATIC_TIME_LINEAR ; }
        if( type == "QUADRATIC_TIME_QUADRATIC") { return QUADRATIC_TIME_QUADRATIC ; }
        if( type == "CUBIC_TIME_LINEAR") { return CUBIC_TIME_LINEAR ; }
        if( type == "CUBIC_TIME_QUADRATIC") { return CUBIC_TIME_QUADRATIC ; }
        if( type == "QUADRIC_TIME_LINEAR") { return QUADRIC_TIME_LINEAR ; }
        if( type == "QUADRIC_TIME_QUADRATIC") { return QUADRIC_TIME_QUADRATIC ; }
        if( type == "QUINTIC_TIME_LINEAR") { return QUINTIC_TIME_LINEAR ; }
        if( type == "QUINTIC_TIME_QUADRATIC") { return QUINTIC_TIME_QUADRATIC ; }
        if( type == "QUADTREE_REFINED") { return QUADTREE_REFINED ; }
        if( type == "REGULAR_GRID") { return REGULAR_GRID ; }
        if(ok) { *ok = false ; }
        return CONSTANT ;
    }
    static std::string fromOrder(Order value)
    {
        switch(value)
        {
            case CONSTANT: return "CONSTANT" ;
            case LINEAR: return "LINEAR" ;
            case QUADRATIC: return "QUADRATIC" ;
            case CUBIC: return "CUBIC" ;
            case QUADRIC: return "QUADRIC" ;
            case QUINTIC: return "QUINTIC" ;
            case CONSTANT_TIME_LINEAR: return "CONSTANT_TIME_LINEAR" ;
            case CONSTANT_TIME_QUADRATIC: return "CONSTANT_TIME_QUADRATIC" ;
            case LINEAR_TIME_LINEAR: return "LINEAR_TIME_LINEAR" ;
            case LINEAR_TIME_QUADRATIC: return "LINEAR_TIME_QUADRATIC" ;
            case QUADRATIC_TIME_LINEAR: return "QUADRATIC_TIME_LINEAR" ;
            case QUADRATIC_TIME_QUADRATIC: return "QUADRATIC_TIME_QUADRATIC" ;
            case CUBIC_TIME_LINEAR: return "CUBIC_TIME_LINEAR" ;
            case CUBIC_TIME_QUADRATIC: return "CUBIC_TIME_QUADRATIC" ;
            case QUADRIC_TIME_LINEAR: return "QUADRIC_TIME_LINEAR" ;
            case QUADRIC_TIME_QUADRATIC: return "QUADRIC_TIME_QUADRATIC" ;
            case QUINTIC_TIME_LINEAR: return "QUINTIC_TIME_LINEAR" ;
            case QUINTIC_TIME_QUADRATIC: return "QUINTIC_TIME_QUADRATIC" ;
            case QUADTREE_REFINED: return "QUADTREE_REFINED" ;
            case REGULAR_GRID: return "REGULAR_GRID" ;
        }
        return "CONSTANT" ;
    }
   
    // parsed from header file: ../elements/integrable_entity.h
    static ParametersType getParametersType(std::string type, bool * ok = 0)
    {
        if(ok) { *ok = true ; }
        if( type == "PURE_LINEAR") { return PURE_LINEAR ; }
        if( type == "LINEAR_AND_CONSTANT") { return LINEAR_AND_CONSTANT ; }
        if( type == "NON_LINEAR") { return NON_LINEAR ; }
        if( type == "VOID_BEHAVIOUR") { return VOID_BEHAVIOUR ; }
        if(ok) { *ok = false ; }
        return PURE_LINEAR ;
    }
    static std::string fromParametersType(ParametersType value)
    {
        switch(value)
        {
            case PURE_LINEAR: return "PURE_LINEAR" ;
            case LINEAR_AND_CONSTANT: return "LINEAR_AND_CONSTANT" ;
            case NON_LINEAR: return "NON_LINEAR" ;
            case VOID_BEHAVIOUR: return "VOID_BEHAVIOUR" ;
        }
        return "PURE_LINEAR" ;
    }
   
    // parsed from header file: ../elements/integrable_entity.h
    static StressCalculationMethod getStressCalculationMethod(std::string type, bool * ok = 0)
    {
        if(ok) { *ok = true ; }
        if( type == "REAL_STRESS") { return REAL_STRESS ; }
        if( type == "EFFECTIVE_STRESS") { return EFFECTIVE_STRESS ; }
        if(ok) { *ok = false ; }
        return REAL_STRESS ;
    }
    static std::string fromStressCalculationMethod(StressCalculationMethod value)
    {
        switch(value)
        {
            case REAL_STRESS: return "REAL_STRESS" ;
            case EFFECTIVE_STRESS: return "EFFECTIVE_STRESS" ;
        }
        return "REAL_STRESS" ;
    }
   
    // parsed from header file: ../elements/integrable_entity.h
    static FieldType getFieldType(std::string type, bool * ok = 0)
    {
        if(ok) { *ok = true ; }
        if( type == "DISPLACEMENT_FIELD") { return DISPLACEMENT_FIELD ; }
        if( type == "ENRICHED_DISPLACEMENT_FIELD") { return ENRICHED_DISPLACEMENT_FIELD ; }
        if( type == "SPEED_FIELD") { return SPEED_FIELD ; }
        if( type == "FLUX_FIELD") { return FLUX_FIELD ; }
        if( type == "GRADIENT_FIELD") { return GRADIENT_FIELD ; }
        if( type == "TOTAL_STRAIN_FIELD") { return TOTAL_STRAIN_FIELD ; }
        if( type == "STRAIN_RATE_FIELD") { return STRAIN_RATE_FIELD ; }
        if( type == "MECHANICAL_STRAIN_FIELD") { return MECHANICAL_STRAIN_FIELD ; }
        if( type == "EFFECTIVE_STRESS_FIELD") { return EFFECTIVE_STRESS_FIELD ; }
        if( type == "REAL_STRESS_FIELD") { return REAL_STRESS_FIELD ; }
        if( type == "PRINCIPAL_TOTAL_STRAIN_FIELD") { return PRINCIPAL_TOTAL_STRAIN_FIELD ; }
        if( type == "PRINCIPAL_MECHANICAL_STRAIN_FIELD") { return PRINCIPAL_MECHANICAL_STRAIN_FIELD ; }
        if( type == "PRINCIPAL_EFFECTIVE_STRESS_FIELD") { return PRINCIPAL_EFFECTIVE_STRESS_FIELD ; }
        if( type == "PRINCIPAL_REAL_STRESS_FIELD") { return PRINCIPAL_REAL_STRESS_FIELD ; }
        if( type == "NON_ENRICHED_STRAIN_FIELD") { return NON_ENRICHED_STRAIN_FIELD ; }
        if( type == "NON_ENRICHED_STRAIN_RATE_FIELD") { return NON_ENRICHED_STRAIN_RATE_FIELD ; }
        if( type == "NON_ENRICHED_EFFECTIVE_STRESS_FIELD") { return NON_ENRICHED_EFFECTIVE_STRESS_FIELD ; }
        if( type == "NON_ENRICHED_REAL_STRESS_FIELD") { return NON_ENRICHED_REAL_STRESS_FIELD ; }
        if( type == "VON_MISES_STRAIN_FIELD") { return VON_MISES_STRAIN_FIELD ; }
        if( type == "VON_MISES_REAL_STRESS_FIELD") { return VON_MISES_REAL_STRESS_FIELD ; }
        if( type == "VON_MISES_EFFECTIVE_STRESS_FIELD") { return VON_MISES_EFFECTIVE_STRESS_FIELD ; }
        if( type == "PRINCIPAL_STRESS_ANGLE_FIELD") { return PRINCIPAL_STRESS_ANGLE_FIELD ; }
        if( type == "PRINCIPAL_STRAIN_ANGLE_FIELD") { return PRINCIPAL_STRAIN_ANGLE_FIELD ; }
        if( type == "INTERNAL_VARIABLE_FIELD") { return INTERNAL_VARIABLE_FIELD ; }
        if( type == "IMPOSED_STRESS_FIELD") { return IMPOSED_STRESS_FIELD ; }
        if( type == "IMPOSED_STRAIN_FIELD") { return IMPOSED_STRAIN_FIELD ; }
        if( type == "PRINCIPAL_IMPOSED_STRAIN_FIELD") { return PRINCIPAL_IMPOSED_STRAIN_FIELD ; }
        if( type == "PRINCIPAL_IMPOSED_STRESS_FIELD") { return PRINCIPAL_IMPOSED_STRESS_FIELD ; }
        if( type == "SCALAR_DAMAGE_FIELD") { return SCALAR_DAMAGE_FIELD ; }
        if( type == "TENSOR_DAMAGE_FIELD") { return TENSOR_DAMAGE_FIELD ; }
        if( type == "GENERALIZED_VISCOELASTIC_DISPLACEMENT_FIELD") { return GENERALIZED_VISCOELASTIC_DISPLACEMENT_FIELD ; }
        if( type == "GENERALIZED_VISCOELASTIC_ENRICHED_DISPLACEMENT_FIELD") { return GENERALIZED_VISCOELASTIC_ENRICHED_DISPLACEMENT_FIELD ; }
        if( type == "GENERALIZED_VISCOELASTIC_SPEED_FIELD") { return GENERALIZED_VISCOELASTIC_SPEED_FIELD ; }
        if( type == "GENERALIZED_VISCOELASTIC_STRAIN_FIELD") { return GENERALIZED_VISCOELASTIC_STRAIN_FIELD ; }
        if( type == "GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD") { return GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD ; }
        if( type == "GENERALIZED_VISCOELASTIC_EFFECTIVE_STRESS_FIELD") { return GENERALIZED_VISCOELASTIC_EFFECTIVE_STRESS_FIELD ; }
        if( type == "GENERALIZED_VISCOELASTIC_REAL_STRESS_FIELD") { return GENERALIZED_VISCOELASTIC_REAL_STRESS_FIELD ; }
        if( type == "GENERALIZED_VISCOELASTIC_PRINCIPAL_STRAIN_FIELD") { return GENERALIZED_VISCOELASTIC_PRINCIPAL_STRAIN_FIELD ; }
        if( type == "GENERALIZED_VISCOELASTIC_PRINCIPAL_EFFECTIVE_STRESS_FIELD") { return GENERALIZED_VISCOELASTIC_PRINCIPAL_EFFECTIVE_STRESS_FIELD ; }
        if( type == "GENERALIZED_VISCOELASTIC_PRINCIPAL_REAL_STRESS_FIELD") { return GENERALIZED_VISCOELASTIC_PRINCIPAL_REAL_STRESS_FIELD ; }
        if( type == "GENERALIZED_VISCOELASTIC_NON_ENRICHED_STRAIN_FIELD") { return GENERALIZED_VISCOELASTIC_NON_ENRICHED_STRAIN_FIELD ; }
        if( type == "GENERALIZED_VISCOELASTIC_NON_ENRICHED_STRAIN_RATE_FIELD") { return GENERALIZED_VISCOELASTIC_NON_ENRICHED_STRAIN_RATE_FIELD ; }
        if( type == "GENERALIZED_VISCOELASTIC_NON_ENRICHED_EFFECTIVE_STRESS_FIELD") { return GENERALIZED_VISCOELASTIC_NON_ENRICHED_EFFECTIVE_STRESS_FIELD ; }
        if( type == "GENERALIZED_VISCOELASTIC_NON_ENRICHED_REAL_STRESS_FIELD") { return GENERALIZED_VISCOELASTIC_NON_ENRICHED_REAL_STRESS_FIELD ; }
        if(ok) { *ok = false ; }
        return DISPLACEMENT_FIELD ;
    }
    static std::string fromFieldType(FieldType value)
    {
        switch(value)
        {
            case DISPLACEMENT_FIELD: return "DISPLACEMENT_FIELD" ;
            case ENRICHED_DISPLACEMENT_FIELD: return "ENRICHED_DISPLACEMENT_FIELD" ;
            case SPEED_FIELD: return "SPEED_FIELD" ;
            case FLUX_FIELD: return "FLUX_FIELD" ;
            case GRADIENT_FIELD: return "GRADIENT_FIELD" ;
            case TOTAL_STRAIN_FIELD: return "TOTAL_STRAIN_FIELD" ;
            case STRAIN_RATE_FIELD: return "STRAIN_RATE_FIELD" ;
            case MECHANICAL_STRAIN_FIELD: return "MECHANICAL_STRAIN_FIELD" ;
            case EFFECTIVE_STRESS_FIELD: return "EFFECTIVE_STRESS_FIELD" ;
            case REAL_STRESS_FIELD: return "REAL_STRESS_FIELD" ;
            case PRINCIPAL_TOTAL_STRAIN_FIELD: return "PRINCIPAL_TOTAL_STRAIN_FIELD" ;
            case PRINCIPAL_MECHANICAL_STRAIN_FIELD: return "PRINCIPAL_MECHANICAL_STRAIN_FIELD" ;
            case PRINCIPAL_EFFECTIVE_STRESS_FIELD: return "PRINCIPAL_EFFECTIVE_STRESS_FIELD" ;
            case PRINCIPAL_REAL_STRESS_FIELD: return "PRINCIPAL_REAL_STRESS_FIELD" ;
            case NON_ENRICHED_STRAIN_FIELD: return "NON_ENRICHED_STRAIN_FIELD" ;
            case NON_ENRICHED_STRAIN_RATE_FIELD: return "NON_ENRICHED_STRAIN_RATE_FIELD" ;
            case NON_ENRICHED_EFFECTIVE_STRESS_FIELD: return "NON_ENRICHED_EFFECTIVE_STRESS_FIELD" ;
            case NON_ENRICHED_REAL_STRESS_FIELD: return "NON_ENRICHED_REAL_STRESS_FIELD" ;
            case VON_MISES_STRAIN_FIELD: return "VON_MISES_STRAIN_FIELD" ;
            case VON_MISES_REAL_STRESS_FIELD: return "VON_MISES_REAL_STRESS_FIELD" ;
            case VON_MISES_EFFECTIVE_STRESS_FIELD: return "VON_MISES_EFFECTIVE_STRESS_FIELD" ;
            case PRINCIPAL_STRESS_ANGLE_FIELD: return "PRINCIPAL_STRESS_ANGLE_FIELD" ;
            case PRINCIPAL_STRAIN_ANGLE_FIELD: return "PRINCIPAL_STRAIN_ANGLE_FIELD" ;
            case INTERNAL_VARIABLE_FIELD: return "INTERNAL_VARIABLE_FIELD" ;
            case IMPOSED_STRESS_FIELD: return "IMPOSED_STRESS_FIELD" ;
            case IMPOSED_STRAIN_FIELD: return "IMPOSED_STRAIN_FIELD" ;
            case PRINCIPAL_IMPOSED_STRAIN_FIELD: return "PRINCIPAL_IMPOSED_STRAIN_FIELD" ;
            case PRINCIPAL_IMPOSED_STRESS_FIELD: return "PRINCIPAL_IMPOSED_STRESS_FIELD" ;
            case SCALAR_DAMAGE_FIELD: return "SCALAR_DAMAGE_FIELD" ;
            case TENSOR_DAMAGE_FIELD: return "TENSOR_DAMAGE_FIELD" ;
            case GENERALIZED_VISCOELASTIC_DISPLACEMENT_FIELD: return "GENERALIZED_VISCOELASTIC_DISPLACEMENT_FIELD" ;
            case GENERALIZED_VISCOELASTIC_ENRICHED_DISPLACEMENT_FIELD: return "GENERALIZED_VISCOELASTIC_ENRICHED_DISPLACEMENT_FIELD" ;
            case GENERALIZED_VISCOELASTIC_SPEED_FIELD: return "GENERALIZED_VISCOELASTIC_SPEED_FIELD" ;
            case GENERALIZED_VISCOELASTIC_STRAIN_FIELD: return "GENERALIZED_VISCOELASTIC_STRAIN_FIELD" ;
            case GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD: return "GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD" ;
            case GENERALIZED_VISCOELASTIC_EFFECTIVE_STRESS_FIELD: return "GENERALIZED_VISCOELASTIC_EFFECTIVE_STRESS_FIELD" ;
            case GENERALIZED_VISCOELASTIC_REAL_STRESS_FIELD: return "GENERALIZED_VISCOELASTIC_REAL_STRESS_FIELD" ;
            case GENERALIZED_VISCOELASTIC_PRINCIPAL_STRAIN_FIELD: return "GENERALIZED_VISCOELASTIC_PRINCIPAL_STRAIN_FIELD" ;
            case GENERALIZED_VISCOELASTIC_PRINCIPAL_EFFECTIVE_STRESS_FIELD: return "GENERALIZED_VISCOELASTIC_PRINCIPAL_EFFECTIVE_STRESS_FIELD" ;
            case GENERALIZED_VISCOELASTIC_PRINCIPAL_REAL_STRESS_FIELD: return "GENERALIZED_VISCOELASTIC_PRINCIPAL_REAL_STRESS_FIELD" ;
            case GENERALIZED_VISCOELASTIC_NON_ENRICHED_STRAIN_FIELD: return "GENERALIZED_VISCOELASTIC_NON_ENRICHED_STRAIN_FIELD" ;
            case GENERALIZED_VISCOELASTIC_NON_ENRICHED_STRAIN_RATE_FIELD: return "GENERALIZED_VISCOELASTIC_NON_ENRICHED_STRAIN_RATE_FIELD" ;
            case GENERALIZED_VISCOELASTIC_NON_ENRICHED_EFFECTIVE_STRESS_FIELD: return "GENERALIZED_VISCOELASTIC_NON_ENRICHED_EFFECTIVE_STRESS_FIELD" ;
            case GENERALIZED_VISCOELASTIC_NON_ENRICHED_REAL_STRESS_FIELD: return "GENERALIZED_VISCOELASTIC_NON_ENRICHED_REAL_STRESS_FIELD" ;
        }
        return "DISPLACEMENT_FIELD" ;
    }
   
    // parsed from header file: ../features/boundarycondition.h
    static LoadingState getLoadingState(std::string type, bool * ok = 0)
    {
        if(ok) { *ok = true ; }
        if( type == "LOADING") { return LOADING ; }
        if( type == "UNLOADING") { return UNLOADING ; }
        if( type == "ULTIMATE_STRAIN") { return ULTIMATE_STRAIN ; }
        if( type == "ULTIMATE_STRESS") { return ULTIMATE_STRESS ; }
        if(ok) { *ok = false ; }
        return LOADING ;
    }
    static std::string fromLoadingState(LoadingState value)
    {
        switch(value)
        {
            case LOADING: return "LOADING" ;
            case UNLOADING: return "UNLOADING" ;
            case ULTIMATE_STRAIN: return "ULTIMATE_STRAIN" ;
            case ULTIMATE_STRESS: return "ULTIMATE_STRESS" ;
        }
        return "LOADING" ;
    }
   
    // parsed from header file: ../geometry/geometry_base.h
    static GeometryType getGeometryType(std::string type, bool * ok = 0)
    {
        if(ok) { *ok = true ; }
        if( type == "nullptr_GEOMETRY") { return nullptr_GEOMETRY ; }
        if( type == "CIRCLE") { return CIRCLE ; }
        if( type == "LAYERED_CIRCLE") { return LAYERED_CIRCLE ; }
        if( type == "TRIANGLE") { return TRIANGLE ; }
        if( type == "RECTANGLE") { return RECTANGLE ; }
        if( type == "PARALLELOGRAMME") { return PARALLELOGRAMME ; }
        if( type == "CONVEX_POLYGON") { return CONVEX_POLYGON ; }
        if( type == "SEGMENTED_LINE") { return SEGMENTED_LINE ; }
        if( type == "POLYGON") { return POLYGON ; }
        if( type == "ORIENTABLE_CIRCLE") { return ORIENTABLE_CIRCLE ; }
        if( type == "CLOSED_NURB") { return CLOSED_NURB ; }
        if( type == "TETRAHEDRON") { return TETRAHEDRON ; }
        if( type == "HEXAHEDRON") { return HEXAHEDRON ; }
        if( type == "SPHERE") { return SPHERE ; }
        if( type == "LAYERED_SPHERE") { return LAYERED_SPHERE ; }
        if( type == "REGULAR_OCTAHEDRON") { return REGULAR_OCTAHEDRON ; }
        if( type == "POLYGON_PRISM") { return POLYGON_PRISM ; }
        if( type == "LOFTED_POLYGON") { return LOFTED_POLYGON ; }
        if( type == "ELLIPSE") { return ELLIPSE ; }
        if( type == "LEVEL_SET") { return LEVEL_SET ; }
        if( type == "TIME_DEPENDENT_CIRCLE") { return TIME_DEPENDENT_CIRCLE ; }
        if(ok) { *ok = false ; }
        return nullptr_GEOMETRY ;
    }
    static std::string fromGeometryType(GeometryType value)
    {
        switch(value)
        {
            case nullptr_GEOMETRY: return "nullptr_GEOMETRY" ;
            case CIRCLE: return "CIRCLE" ;
            case LAYERED_CIRCLE: return "LAYERED_CIRCLE" ;
            case TRIANGLE: return "TRIANGLE" ;
            case RECTANGLE: return "RECTANGLE" ;
            case PARALLELOGRAMME: return "PARALLELOGRAMME" ;
            case CONVEX_POLYGON: return "CONVEX_POLYGON" ;
            case SEGMENTED_LINE: return "SEGMENTED_LINE" ;
            case POLYGON: return "POLYGON" ;
            case ORIENTABLE_CIRCLE: return "ORIENTABLE_CIRCLE" ;
            case CLOSED_NURB: return "CLOSED_NURB" ;
            case TETRAHEDRON: return "TETRAHEDRON" ;
            case HEXAHEDRON: return "HEXAHEDRON" ;
            case SPHERE: return "SPHERE" ;
            case LAYERED_SPHERE: return "LAYERED_SPHERE" ;
            case REGULAR_OCTAHEDRON: return "REGULAR_OCTAHEDRON" ;
            case POLYGON_PRISM: return "POLYGON_PRISM" ;
            case LOFTED_POLYGON: return "LOFTED_POLYGON" ;
            case ELLIPSE: return "ELLIPSE" ;
            case LEVEL_SET: return "LEVEL_SET" ;
            case TIME_DEPENDENT_CIRCLE: return "TIME_DEPENDENT_CIRCLE" ;
        }
        return "nullptr_GEOMETRY" ;
    }
   
    // parsed from header file: ../geometry/geometry_base.h
    static GeometricTransformationType getGeometricTransformationType(std::string type, bool * ok = 0)
    {
        if(ok) { *ok = true ; }
        if( type == "ROTATE") { return ROTATE ; }
        if( type == "SCALE") { return SCALE ; }
        if( type == "TRANSLATE") { return TRANSLATE ; }
        if(ok) { *ok = false ; }
        return ROTATE ;
    }
    static std::string fromGeometricTransformationType(GeometricTransformationType value)
    {
        switch(value)
        {
            case ROTATE: return "ROTATE" ;
            case SCALE: return "SCALE" ;
            case TRANSLATE: return "TRANSLATE" ;
        }
        return "ROTATE" ;
    }
   
    // parsed from header file: ../geometry/geometry_base.h
    static BoundingBoxPosition getBoundingBoxPosition(std::string type, bool * ok = 0)
    {
        if(ok) { *ok = true ; }
        if( type == "TOP") { return TOP ; }
        if( type == "LEFT") { return LEFT ; }
        if( type == "BOTTOM") { return BOTTOM ; }
        if( type == "RIGHT") { return RIGHT ; }
        if( type == "FRONT") { return FRONT ; }
        if( type == "BACK") { return BACK ; }
        if( type == "BEFORE") { return BEFORE ; }
        if( type == "NOW") { return NOW ; }
        if( type == "AFTER") { return AFTER ; }
        if( type == "TOP_LEFT") { return TOP_LEFT ; }
        if( type == "TOP_RIGHT") { return TOP_RIGHT ; }
        if( type == "BOTTOM_LEFT") { return BOTTOM_LEFT ; }
        if( type == "BOTTOM_RIGHT") { return BOTTOM_RIGHT ; }
        if( type == "FRONT_LEFT") { return FRONT_LEFT ; }
        if( type == "FRONT_RIGHT") { return FRONT_RIGHT ; }
        if( type == "BACK_LEFT") { return BACK_LEFT ; }
        if( type == "BACK_RIGHT") { return BACK_RIGHT ; }
        if( type == "FRONT_TOP") { return FRONT_TOP ; }
        if( type == "FRONT_BOTTOM") { return FRONT_BOTTOM ; }
        if( type == "BOTTOM_BACK") { return BOTTOM_BACK ; }
        if( type == "TOP_BACK") { return TOP_BACK ; }
        if( type == "TOP_LEFT_FRONT") { return TOP_LEFT_FRONT ; }
        if( type == "TOP_LEFT_BACK") { return TOP_LEFT_BACK ; }
        if( type == "BOTTOM_LEFT_FRONT") { return BOTTOM_LEFT_FRONT ; }
        if( type == "BOTTOM_LEFT_BACK") { return BOTTOM_LEFT_BACK ; }
        if( type == "TOP_RIGHT_FRONT") { return TOP_RIGHT_FRONT ; }
        if( type == "TOP_RIGHT_BACK") { return TOP_RIGHT_BACK ; }
        if( type == "BOTTOM_RIGHT_FRONT") { return BOTTOM_RIGHT_FRONT ; }
        if( type == "BOTTOM_RIGHT_BACK") { return BOTTOM_RIGHT_BACK ; }
        if( type == "TOP_BEFORE") { return TOP_BEFORE ; }
        if( type == "LEFT_BEFORE") { return LEFT_BEFORE ; }
        if( type == "BOTTOM_BEFORE") { return BOTTOM_BEFORE ; }
        if( type == "RIGHT_BEFORE") { return RIGHT_BEFORE ; }
        if( type == "FRONT_BEFORE") { return FRONT_BEFORE ; }
        if( type == "BACK_BEFORE") { return BACK_BEFORE ; }
        if( type == "TOP_LEFT_BEFORE") { return TOP_LEFT_BEFORE ; }
        if( type == "TOP_RIGHT_BEFORE") { return TOP_RIGHT_BEFORE ; }
        if( type == "BOTTOM_LEFT_BEFORE") { return BOTTOM_LEFT_BEFORE ; }
        if( type == "BOTTOM_RIGHT_BEFORE") { return BOTTOM_RIGHT_BEFORE ; }
        if( type == "FRONT_LEFT_BEFORE") { return FRONT_LEFT_BEFORE ; }
        if( type == "FRONT_RIGHT_BEFORE") { return FRONT_RIGHT_BEFORE ; }
        if( type == "BACK_LEFT_BEFORE") { return BACK_LEFT_BEFORE ; }
        if( type == "BACK_RIGHT_BEFORE") { return BACK_RIGHT_BEFORE ; }
        if( type == "FRONT_TOP_BEFORE") { return FRONT_TOP_BEFORE ; }
        if( type == "FRONT_BOTTOM_BEFORE") { return FRONT_BOTTOM_BEFORE ; }
        if( type == "TOP_LEFT_FRONT_BEFORE") { return TOP_LEFT_FRONT_BEFORE ; }
        if( type == "TOP_LEFT_BACK_BEFORE") { return TOP_LEFT_BACK_BEFORE ; }
        if( type == "BOTTOM_LEFT_FRONT_BEFORE") { return BOTTOM_LEFT_FRONT_BEFORE ; }
        if( type == "BOTTOM_LEFT_BACK_BEFORE") { return BOTTOM_LEFT_BACK_BEFORE ; }
        if( type == "TOP_RIGHT_FRONT_BEFORE") { return TOP_RIGHT_FRONT_BEFORE ; }
        if( type == "TOP_RIGHT_BACK_BEFORE") { return TOP_RIGHT_BACK_BEFORE ; }
        if( type == "BOTTOM_RIGHT_FRONT_BEFORE") { return BOTTOM_RIGHT_FRONT_BEFORE ; }
        if( type == "BOTTOM_RIGHT_BACK_BEFORE") { return BOTTOM_RIGHT_BACK_BEFORE ; }
        if( type == "BOTTOM_BACK_BEFORE") { return BOTTOM_BACK_BEFORE ; }
        if( type == "TOP_BACK_BEFORE") { return TOP_BACK_BEFORE ; }
        if( type == "TOP_NOW") { return TOP_NOW ; }
        if( type == "LEFT_NOW") { return LEFT_NOW ; }
        if( type == "BOTTOM_NOW") { return BOTTOM_NOW ; }
        if( type == "RIGHT_NOW") { return RIGHT_NOW ; }
        if( type == "FRONT_NOW") { return FRONT_NOW ; }
        if( type == "BACK_NOW") { return BACK_NOW ; }
        if( type == "TOP_LEFT_NOW") { return TOP_LEFT_NOW ; }
        if( type == "TOP_RIGHT_NOW") { return TOP_RIGHT_NOW ; }
        if( type == "BOTTOM_LEFT_NOW") { return BOTTOM_LEFT_NOW ; }
        if( type == "BOTTOM_RIGHT_NOW") { return BOTTOM_RIGHT_NOW ; }
        if( type == "FRONT_LEFT_NOW") { return FRONT_LEFT_NOW ; }
        if( type == "FRONT_RIGHT_NOW") { return FRONT_RIGHT_NOW ; }
        if( type == "BACK_LEFT_NOW") { return BACK_LEFT_NOW ; }
        if( type == "BACK_RIGHT_NOW") { return BACK_RIGHT_NOW ; }
        if( type == "FRONT_TOP_NOW") { return FRONT_TOP_NOW ; }
        if( type == "BOTTOM_BACK_NOW") { return BOTTOM_BACK_NOW ; }
        if( type == "TOP_BACK_NOW") { return TOP_BACK_NOW ; }
        if( type == "FRONT_BOTTOM_NOW") { return FRONT_BOTTOM_NOW ; }
        if( type == "TOP_LEFT_FRONT_NOW") { return TOP_LEFT_FRONT_NOW ; }
        if( type == "TOP_LEFT_BACK_NOW") { return TOP_LEFT_BACK_NOW ; }
        if( type == "BOTTOM_LEFT_FRONT_NOW") { return BOTTOM_LEFT_FRONT_NOW ; }
        if( type == "BOTTOM_LEFT_BACK_NOW") { return BOTTOM_LEFT_BACK_NOW ; }
        if( type == "TOP_RIGHT_FRONT_NOW") { return TOP_RIGHT_FRONT_NOW ; }
        if( type == "TOP_RIGHT_BACK_NOW") { return TOP_RIGHT_BACK_NOW ; }
        if( type == "BOTTOM_RIGHT_FRONT_NOW") { return BOTTOM_RIGHT_FRONT_NOW ; }
        if( type == "BOTTOM_RIGHT_BACK_NOW") { return BOTTOM_RIGHT_BACK_NOW ; }
        if( type == "TOP_AFTER") { return TOP_AFTER ; }
        if( type == "LEFT_AFTER") { return LEFT_AFTER ; }
        if( type == "BOTTOM_AFTER") { return BOTTOM_AFTER ; }
        if( type == "RIGHT_AFTER") { return RIGHT_AFTER ; }
        if( type == "FRONT_AFTER") { return FRONT_AFTER ; }
        if( type == "BACK_AFTER") { return BACK_AFTER ; }
        if( type == "TOP_LEFT_AFTER") { return TOP_LEFT_AFTER ; }
        if( type == "TOP_RIGHT_AFTER") { return TOP_RIGHT_AFTER ; }
        if( type == "BOTTOM_LEFT_AFTER") { return BOTTOM_LEFT_AFTER ; }
        if( type == "BOTTOM_RIGHT_AFTER") { return BOTTOM_RIGHT_AFTER ; }
        if( type == "FRONT_LEFT_AFTER") { return FRONT_LEFT_AFTER ; }
        if( type == "FRONT_RIGHT_AFTER") { return FRONT_RIGHT_AFTER ; }
        if( type == "BACK_LEFT_AFTER") { return BACK_LEFT_AFTER ; }
        if( type == "BACK_RIGHT_AFTER") { return BACK_RIGHT_AFTER ; }
        if( type == "FRONT_TOP_AFTER") { return FRONT_TOP_AFTER ; }
        if( type == "FRONT_BOTTOM_AFTER") { return FRONT_BOTTOM_AFTER ; }
        if( type == "TOP_LEFT_FRONT_AFTER") { return TOP_LEFT_FRONT_AFTER ; }
        if( type == "TOP_LEFT_BACK_AFTER") { return TOP_LEFT_BACK_AFTER ; }
        if( type == "BOTTOM_LEFT_FRONT_AFTER") { return BOTTOM_LEFT_FRONT_AFTER ; }
        if( type == "BOTTOM_LEFT_BACK_AFTER") { return BOTTOM_LEFT_BACK_AFTER ; }
        if( type == "TOP_RIGHT_FRONT_AFTER") { return TOP_RIGHT_FRONT_AFTER ; }
        if( type == "TOP_RIGHT_BACK_AFTER") { return TOP_RIGHT_BACK_AFTER ; }
        if( type == "BOTTOM_RIGHT_FRONT_AFTER") { return BOTTOM_RIGHT_FRONT_AFTER ; }
        if( type == "BOTTOM_RIGHT_BACK_AFTER") { return BOTTOM_RIGHT_BACK_AFTER ; }
        if( type == "BOTTOM_BACK_AFTER") { return BOTTOM_BACK_AFTER ; }
        if( type == "TOP_BACK_AFTER") { return TOP_BACK_AFTER ; }
        if(ok) { *ok = false ; }
        return TOP ;
    }
    static std::string fromBoundingBoxPosition(BoundingBoxPosition value)
    {
        switch(value)
        {
            case TOP: return "TOP" ;
            case LEFT: return "LEFT" ;
            case BOTTOM: return "BOTTOM" ;
            case RIGHT: return "RIGHT" ;
            case FRONT: return "FRONT" ;
            case BACK: return "BACK" ;
            case BEFORE: return "BEFORE" ;
            case NOW: return "NOW" ;
            case AFTER: return "AFTER" ;
            case TOP_LEFT: return "TOP_LEFT" ;
            case TOP_RIGHT: return "TOP_RIGHT" ;
            case BOTTOM_LEFT: return "BOTTOM_LEFT" ;
            case BOTTOM_RIGHT: return "BOTTOM_RIGHT" ;
            case FRONT_LEFT: return "FRONT_LEFT" ;
            case FRONT_RIGHT: return "FRONT_RIGHT" ;
            case BACK_LEFT: return "BACK_LEFT" ;
            case BACK_RIGHT: return "BACK_RIGHT" ;
            case FRONT_TOP: return "FRONT_TOP" ;
            case FRONT_BOTTOM: return "FRONT_BOTTOM" ;
            case BOTTOM_BACK: return "BOTTOM_BACK" ;
            case TOP_BACK: return "TOP_BACK" ;
            case TOP_LEFT_FRONT: return "TOP_LEFT_FRONT" ;
            case TOP_LEFT_BACK: return "TOP_LEFT_BACK" ;
            case BOTTOM_LEFT_FRONT: return "BOTTOM_LEFT_FRONT" ;
            case BOTTOM_LEFT_BACK: return "BOTTOM_LEFT_BACK" ;
            case TOP_RIGHT_FRONT: return "TOP_RIGHT_FRONT" ;
            case TOP_RIGHT_BACK: return "TOP_RIGHT_BACK" ;
            case BOTTOM_RIGHT_FRONT: return "BOTTOM_RIGHT_FRONT" ;
            case BOTTOM_RIGHT_BACK: return "BOTTOM_RIGHT_BACK" ;
            case TOP_BEFORE: return "TOP_BEFORE" ;
            case LEFT_BEFORE: return "LEFT_BEFORE" ;
            case BOTTOM_BEFORE: return "BOTTOM_BEFORE" ;
            case RIGHT_BEFORE: return "RIGHT_BEFORE" ;
            case FRONT_BEFORE: return "FRONT_BEFORE" ;
            case BACK_BEFORE: return "BACK_BEFORE" ;
            case TOP_LEFT_BEFORE: return "TOP_LEFT_BEFORE" ;
            case TOP_RIGHT_BEFORE: return "TOP_RIGHT_BEFORE" ;
            case BOTTOM_LEFT_BEFORE: return "BOTTOM_LEFT_BEFORE" ;
            case BOTTOM_RIGHT_BEFORE: return "BOTTOM_RIGHT_BEFORE" ;
            case FRONT_LEFT_BEFORE: return "FRONT_LEFT_BEFORE" ;
            case FRONT_RIGHT_BEFORE: return "FRONT_RIGHT_BEFORE" ;
            case BACK_LEFT_BEFORE: return "BACK_LEFT_BEFORE" ;
            case BACK_RIGHT_BEFORE: return "BACK_RIGHT_BEFORE" ;
            case FRONT_TOP_BEFORE: return "FRONT_TOP_BEFORE" ;
            case FRONT_BOTTOM_BEFORE: return "FRONT_BOTTOM_BEFORE" ;
            case TOP_LEFT_FRONT_BEFORE: return "TOP_LEFT_FRONT_BEFORE" ;
            case TOP_LEFT_BACK_BEFORE: return "TOP_LEFT_BACK_BEFORE" ;
            case BOTTOM_LEFT_FRONT_BEFORE: return "BOTTOM_LEFT_FRONT_BEFORE" ;
            case BOTTOM_LEFT_BACK_BEFORE: return "BOTTOM_LEFT_BACK_BEFORE" ;
            case TOP_RIGHT_FRONT_BEFORE: return "TOP_RIGHT_FRONT_BEFORE" ;
            case TOP_RIGHT_BACK_BEFORE: return "TOP_RIGHT_BACK_BEFORE" ;
            case BOTTOM_RIGHT_FRONT_BEFORE: return "BOTTOM_RIGHT_FRONT_BEFORE" ;
            case BOTTOM_RIGHT_BACK_BEFORE: return "BOTTOM_RIGHT_BACK_BEFORE" ;
            case BOTTOM_BACK_BEFORE: return "BOTTOM_BACK_BEFORE" ;
            case TOP_BACK_BEFORE: return "TOP_BACK_BEFORE" ;
            case TOP_NOW: return "TOP_NOW" ;
            case LEFT_NOW: return "LEFT_NOW" ;
            case BOTTOM_NOW: return "BOTTOM_NOW" ;
            case RIGHT_NOW: return "RIGHT_NOW" ;
            case FRONT_NOW: return "FRONT_NOW" ;
            case BACK_NOW: return "BACK_NOW" ;
            case TOP_LEFT_NOW: return "TOP_LEFT_NOW" ;
            case TOP_RIGHT_NOW: return "TOP_RIGHT_NOW" ;
            case BOTTOM_LEFT_NOW: return "BOTTOM_LEFT_NOW" ;
            case BOTTOM_RIGHT_NOW: return "BOTTOM_RIGHT_NOW" ;
            case FRONT_LEFT_NOW: return "FRONT_LEFT_NOW" ;
            case FRONT_RIGHT_NOW: return "FRONT_RIGHT_NOW" ;
            case BACK_LEFT_NOW: return "BACK_LEFT_NOW" ;
            case BACK_RIGHT_NOW: return "BACK_RIGHT_NOW" ;
            case FRONT_TOP_NOW: return "FRONT_TOP_NOW" ;
            case BOTTOM_BACK_NOW: return "BOTTOM_BACK_NOW" ;
            case TOP_BACK_NOW: return "TOP_BACK_NOW" ;
            case FRONT_BOTTOM_NOW: return "FRONT_BOTTOM_NOW" ;
            case TOP_LEFT_FRONT_NOW: return "TOP_LEFT_FRONT_NOW" ;
            case TOP_LEFT_BACK_NOW: return "TOP_LEFT_BACK_NOW" ;
            case BOTTOM_LEFT_FRONT_NOW: return "BOTTOM_LEFT_FRONT_NOW" ;
            case BOTTOM_LEFT_BACK_NOW: return "BOTTOM_LEFT_BACK_NOW" ;
            case TOP_RIGHT_FRONT_NOW: return "TOP_RIGHT_FRONT_NOW" ;
            case TOP_RIGHT_BACK_NOW: return "TOP_RIGHT_BACK_NOW" ;
            case BOTTOM_RIGHT_FRONT_NOW: return "BOTTOM_RIGHT_FRONT_NOW" ;
            case BOTTOM_RIGHT_BACK_NOW: return "BOTTOM_RIGHT_BACK_NOW" ;
            case TOP_AFTER: return "TOP_AFTER" ;
            case LEFT_AFTER: return "LEFT_AFTER" ;
            case BOTTOM_AFTER: return "BOTTOM_AFTER" ;
            case RIGHT_AFTER: return "RIGHT_AFTER" ;
            case FRONT_AFTER: return "FRONT_AFTER" ;
            case BACK_AFTER: return "BACK_AFTER" ;
            case TOP_LEFT_AFTER: return "TOP_LEFT_AFTER" ;
            case TOP_RIGHT_AFTER: return "TOP_RIGHT_AFTER" ;
            case BOTTOM_LEFT_AFTER: return "BOTTOM_LEFT_AFTER" ;
            case BOTTOM_RIGHT_AFTER: return "BOTTOM_RIGHT_AFTER" ;
            case FRONT_LEFT_AFTER: return "FRONT_LEFT_AFTER" ;
            case FRONT_RIGHT_AFTER: return "FRONT_RIGHT_AFTER" ;
            case BACK_LEFT_AFTER: return "BACK_LEFT_AFTER" ;
            case BACK_RIGHT_AFTER: return "BACK_RIGHT_AFTER" ;
            case FRONT_TOP_AFTER: return "FRONT_TOP_AFTER" ;
            case FRONT_BOTTOM_AFTER: return "FRONT_BOTTOM_AFTER" ;
            case TOP_LEFT_FRONT_AFTER: return "TOP_LEFT_FRONT_AFTER" ;
            case TOP_LEFT_BACK_AFTER: return "TOP_LEFT_BACK_AFTER" ;
            case BOTTOM_LEFT_FRONT_AFTER: return "BOTTOM_LEFT_FRONT_AFTER" ;
            case BOTTOM_LEFT_BACK_AFTER: return "BOTTOM_LEFT_BACK_AFTER" ;
            case TOP_RIGHT_FRONT_AFTER: return "TOP_RIGHT_FRONT_AFTER" ;
            case TOP_RIGHT_BACK_AFTER: return "TOP_RIGHT_BACK_AFTER" ;
            case BOTTOM_RIGHT_FRONT_AFTER: return "BOTTOM_RIGHT_FRONT_AFTER" ;
            case BOTTOM_RIGHT_BACK_AFTER: return "BOTTOM_RIGHT_BACK_AFTER" ;
            case BOTTOM_BACK_AFTER: return "BOTTOM_BACK_AFTER" ;
            case TOP_BACK_AFTER: return "TOP_BACK_AFTER" ;
        }
        return "TOP" ;
    }
   
    // parsed from header file: ../geometry/geometry_base.h
    static SpaceDimensionality getSpaceDimensionality(std::string type, bool * ok = 0)
    {
        if(ok) { *ok = true ; }
        if( type == "SPACE_ONE_DIMENSIONAL") { return SPACE_ONE_DIMENSIONAL ; }
        if( type == "SPACE_TWO_DIMENSIONAL") { return SPACE_TWO_DIMENSIONAL ; }
        if( type == "SPACE_THREE_DIMENSIONAL") { return SPACE_THREE_DIMENSIONAL ; }
        if(ok) { *ok = false ; }
        return SPACE_ONE_DIMENSIONAL ;
    }
    static std::string fromSpaceDimensionality(SpaceDimensionality value)
    {
        switch(value)
        {
            case SPACE_ONE_DIMENSIONAL: return "SPACE_ONE_DIMENSIONAL" ;
            case SPACE_TWO_DIMENSIONAL: return "SPACE_TWO_DIMENSIONAL" ;
            case SPACE_THREE_DIMENSIONAL: return "SPACE_THREE_DIMENSIONAL" ;
        }
        return "SPACE_ONE_DIMENSIONAL" ;
    }
   
    // parsed from header file: ../physics/finite_difference_viscoelasticity.h
    static ViscoelasticFiniteDifferenceIntegration getViscoelasticFiniteDifferenceIntegration(std::string type, bool * ok = 0)
    {
        if(ok) { *ok = true ; }
        if( type == "FORWARD_EULER") { return FORWARD_EULER ; }
        if( type == "BACKWARD_EULER") { return BACKWARD_EULER ; }
        if( type == "CENTRAL_DIFFERENCE") { return CENTRAL_DIFFERENCE ; }
        if( type == "NEWMARK") { return NEWMARK ; }
        if( type == "ZIENKIEWICZ") { return ZIENKIEWICZ ; }
        if(ok) { *ok = false ; }
        return FORWARD_EULER ;
    }
    static std::string fromViscoelasticFiniteDifferenceIntegration(ViscoelasticFiniteDifferenceIntegration value)
    {
        switch(value)
        {
            case FORWARD_EULER: return "FORWARD_EULER" ;
            case BACKWARD_EULER: return "BACKWARD_EULER" ;
            case CENTRAL_DIFFERENCE: return "CENTRAL_DIFFERENCE" ;
            case NEWMARK: return "NEWMARK" ;
            case ZIENKIEWICZ: return "ZIENKIEWICZ" ;
        }
        return "FORWARD_EULER" ;
    }
   
    // parsed from header file: ../physics/materials/csh_behaviour.h
    static CSHType getCSHType(std::string type, bool * ok = 0)
    {
        if(ok) { *ok = true ; }
        if( type == "INNER_CSH") { return INNER_CSH ; }
        if( type == "OUTER_CSH") { return OUTER_CSH ; }
        if(ok) { *ok = false ; }
        return INNER_CSH ;
    }
    static std::string fromCSHType(CSHType value)
    {
        switch(value)
        {
            case INNER_CSH: return "INNER_CSH" ;
            case OUTER_CSH: return "OUTER_CSH" ;
        }
        return "INNER_CSH" ;
    }
   
    // parsed from header file: ../physics/fracturecriteria/fracturecriterion.h
    static SmoothingFunctionType getSmoothingFunctionType(std::string type, bool * ok = 0)
    {
        if(ok) { *ok = true ; }
        if( type == "QUARTIC_COMPACT") { return QUARTIC_COMPACT ; }
        if( type == "GAUSSIAN_NONCOMPACT") { return GAUSSIAN_NONCOMPACT ; }
        if( type == "LINEAR_COMPACT") { return LINEAR_COMPACT ; }
        if(ok) { *ok = false ; }
        return QUARTIC_COMPACT ;
    }
    static std::string fromSmoothingFunctionType(SmoothingFunctionType value)
    {
        switch(value)
        {
            case QUARTIC_COMPACT: return "QUARTIC_COMPACT" ;
            case GAUSSIAN_NONCOMPACT: return "GAUSSIAN_NONCOMPACT" ;
            case LINEAR_COMPACT: return "LINEAR_COMPACT" ;
        }
        return "QUARTIC_COMPACT" ;
    }
   
    // parsed from header file: ../physics/fracturecriteria/spacetimelimitsurfacefracturecriterion.h
    static ReferenceFrame getReferenceFrame(std::string type, bool * ok = 0)
    {
        if(ok) { *ok = true ; }
        if( type == "FRAME_CARTESIAN") { return FRAME_CARTESIAN ; }
        if( type == "FRAME_PRINCIPAL") { return FRAME_PRINCIPAL ; }
        if( type == "FRAME_MODEL") { return FRAME_MODEL ; }
        if( type == "FRAME_MATERIAL") { return FRAME_MATERIAL ; }
        if(ok) { *ok = false ; }
        return FRAME_CARTESIAN ;
    }
    static std::string fromReferenceFrame(ReferenceFrame value)
    {
        switch(value)
        {
            case FRAME_CARTESIAN: return "FRAME_CARTESIAN" ;
            case FRAME_PRINCIPAL: return "FRAME_PRINCIPAL" ;
            case FRAME_MODEL: return "FRAME_MODEL" ;
            case FRAME_MATERIAL: return "FRAME_MATERIAL" ;
        }
        return "FRAME_CARTESIAN" ;
    }
   
    // parsed from header file: ../physics/fracturecriteria/mcft.h
    static RedistributionType getRedistributionType(std::string type, bool * ok = 0)
    {
        if(ok) { *ok = true ; }
        if( type == "UPPER_BOUND") { return UPPER_BOUND ; }
        if( type == "LOWER_BOUND") { return LOWER_BOUND ; }
        if( type == "AVERAGE") { return AVERAGE ; }
        if(ok) { *ok = false ; }
        return UPPER_BOUND ;
    }
    static std::string fromRedistributionType(RedistributionType value)
    {
        switch(value)
        {
            case UPPER_BOUND: return "UPPER_BOUND" ;
            case LOWER_BOUND: return "LOWER_BOUND" ;
            case AVERAGE: return "AVERAGE" ;
        }
        return "UPPER_BOUND" ;
    }
   
    // parsed from header file: ../physics/material_laws/material_laws.h
    static EMLOperation getEMLOperation(std::string type, bool * ok = 0)
    {
        if(ok) { *ok = true ; }
        if( type == "SET") { return SET ; }
        if( type == "ADD") { return ADD ; }
        if( type == "MULTIPLY") { return MULTIPLY ; }
        if( type == "SUBSTRACT") { return SUBSTRACT ; }
        if( type == "DIVIDE") { return DIVIDE ; }
        if(ok) { *ok = false ; }
        return SET ;
    }
    static std::string fromEMLOperation(EMLOperation value)
    {
        switch(value)
        {
            case SET: return "SET" ;
            case ADD: return "ADD" ;
            case MULTIPLY: return "MULTIPLY" ;
            case SUBSTRACT: return "SUBSTRACT" ;
            case DIVIDE: return "DIVIDE" ;
        }
        return "SET" ;
    }
   
    // parsed from header file: ../physics/viscoelasticity.h
    static ViscoelasticModel getViscoelasticModel(std::string type, bool * ok = 0)
    {
        if(ok) { *ok = true ; }
        if( type == "PURE_ELASTICITY") { return PURE_ELASTICITY ; }
        if( type == "PURE_VISCOSITY") { return PURE_VISCOSITY ; }
        if( type == "KELVIN_VOIGT") { return KELVIN_VOIGT ; }
        if( type == "MAXWELL") { return MAXWELL ; }
        if( type == "BURGER") { return BURGER ; }
        if( type == "GENERALIZED_KELVIN_VOIGT") { return GENERALIZED_KELVIN_VOIGT ; }
        if( type == "GENERALIZED_MAXWELL") { return GENERALIZED_MAXWELL ; }
        if( type == "GENERAL_VISCOELASTICITY") { return GENERAL_VISCOELASTICITY ; }
        if(ok) { *ok = false ; }
        return PURE_ELASTICITY ;
    }
    static std::string fromViscoelasticModel(ViscoelasticModel value)
    {
        switch(value)
        {
            case PURE_ELASTICITY: return "PURE_ELASTICITY" ;
            case PURE_VISCOSITY: return "PURE_VISCOSITY" ;
            case KELVIN_VOIGT: return "KELVIN_VOIGT" ;
            case MAXWELL: return "MAXWELL" ;
            case BURGER: return "BURGER" ;
            case GENERALIZED_KELVIN_VOIGT: return "GENERALIZED_KELVIN_VOIGT" ;
            case GENERALIZED_MAXWELL: return "GENERALIZED_MAXWELL" ;
            case GENERAL_VISCOELASTICITY: return "GENERAL_VISCOELASTICITY" ;
        }
        return "PURE_ELASTICITY" ;
    }
   
    // parsed from header file: ../physics/viscoelasticity.h
    static CreepComplianceModel getCreepComplianceModel(std::string type, bool * ok = 0)
    {
        if(ok) { *ok = true ; }
        if( type == "LOGPOWER_CREEP") { return LOGPOWER_CREEP ; }
        if( type == "ACI_CREEP") { return ACI_CREEP ; }
        if( type == "CEB_CREEP") { return CEB_CREEP ; }
        if( type == "B3_DRYING_CREEP") { return B3_DRYING_CREEP ; }
        if( type == "JSCE_CREEP") { return JSCE_CREEP ; }
        if( type == "FIB_CREEP") { return FIB_CREEP ; }
        if(ok) { *ok = false ; }
        return LOGPOWER_CREEP ;
    }
    static std::string fromCreepComplianceModel(CreepComplianceModel value)
    {
        switch(value)
        {
            case LOGPOWER_CREEP: return "LOGPOWER_CREEP" ;
            case ACI_CREEP: return "ACI_CREEP" ;
            case CEB_CREEP: return "CEB_CREEP" ;
            case B3_DRYING_CREEP: return "B3_DRYING_CREEP" ;
            case JSCE_CREEP: return "JSCE_CREEP" ;
            case FIB_CREEP: return "FIB_CREEP" ;
        }
        return "LOGPOWER_CREEP" ;
    }
   
    // parsed from header file: ../physics/damagemodels/damagemodel.h
    static ConvergenceType getConvergenceType(std::string type, bool * ok = 0)
    {
        if(ok) { *ok = true ; }
        if( type == "DISSIPATIVE") { return DISSIPATIVE ; }
        if( type == "CONSERVATIVE") { return CONSERVATIVE ; }
        if(ok) { *ok = false ; }
        return DISSIPATIVE ;
    }
    static std::string fromConvergenceType(ConvergenceType value)
    {
        switch(value)
        {
            case DISSIPATIVE: return "DISSIPATIVE" ;
            case CONSERVATIVE: return "CONSERVATIVE" ;
        }
        return "DISSIPATIVE" ;
    }
   
    // parsed from header file: ../polynomial/vm_token.h
    static TokenOperationType getTokenOperationType(std::string type, bool * ok = 0)
    {
        if(ok) { *ok = true ; }
        if( type == "TOKEN_OPERATION_CONSTANT") { return TOKEN_OPERATION_CONSTANT ; }
        if( type == "TOKEN_OPERATION_X") { return TOKEN_OPERATION_X ; }
        if( type == "TOKEN_OPERATION_Y") { return TOKEN_OPERATION_Y ; }
        if( type == "TOKEN_OPERATION_Z") { return TOKEN_OPERATION_Z ; }
        if( type == "TOKEN_OPERATION_T") { return TOKEN_OPERATION_T ; }
        if( type == "TOKEN_OPERATION_U") { return TOKEN_OPERATION_U ; }
        if( type == "TOKEN_OPERATION_V") { return TOKEN_OPERATION_V ; }
        if( type == "TOKEN_OPERATION_W") { return TOKEN_OPERATION_W ; }
        if( type == "TOKEN_OPERATION_PLUS") { return TOKEN_OPERATION_PLUS ; }
        if( type == "TOKEN_OPERATION_INPLACE_PLUS") { return TOKEN_OPERATION_INPLACE_PLUS ; }
        if( type == "TOKEN_OPERATION_MINUS") { return TOKEN_OPERATION_MINUS ; }
        if( type == "TOKEN_OPERATION_INPLACE_MINUS") { return TOKEN_OPERATION_INPLACE_MINUS ; }
        if( type == "TOKEN_OPERATION_TIMES") { return TOKEN_OPERATION_TIMES ; }
        if( type == "TOKEN_OPERATION_INPLACE_TIMES") { return TOKEN_OPERATION_INPLACE_TIMES ; }
        if( type == "TOKEN_OPERATION_DIVIDES") { return TOKEN_OPERATION_DIVIDES ; }
        if( type == "TOKEN_OPERATION_INPLACE_DIVIDES") { return TOKEN_OPERATION_INPLACE_DIVIDES ; }
        if( type == "TOKEN_OPERATION_POWER") { return TOKEN_OPERATION_POWER ; }
        if( type == "TOKEN_OPERATION_INPLACE_POWER") { return TOKEN_OPERATION_INPLACE_POWER ; }
        if( type == "TOKEN_OPERATION_ABS") { return TOKEN_OPERATION_ABS ; }
        if( type == "TOKEN_OPERATION_INPLACE_ABS") { return TOKEN_OPERATION_INPLACE_ABS ; }
        if( type == "TOKEN_OPERATION_COS") { return TOKEN_OPERATION_COS ; }
        if( type == "TOKEN_OPERATION_INPLACE_COS") { return TOKEN_OPERATION_INPLACE_COS ; }
        if( type == "TOKEN_OPERATION_SIN") { return TOKEN_OPERATION_SIN ; }
        if( type == "TOKEN_OPERATION_INPLACE_SIN") { return TOKEN_OPERATION_INPLACE_SIN ; }
        if( type == "TOKEN_OPERATION_TAN") { return TOKEN_OPERATION_TAN ; }
        if( type == "TOKEN_OPERATION_INPLACE_TAN") { return TOKEN_OPERATION_INPLACE_TAN ; }
        if( type == "TOKEN_OPERATION_COSH") { return TOKEN_OPERATION_COSH ; }
        if( type == "TOKEN_OPERATION_INPLACE_COSH") { return TOKEN_OPERATION_INPLACE_COSH ; }
        if( type == "TOKEN_OPERATION_SINH") { return TOKEN_OPERATION_SINH ; }
        if( type == "TOKEN_OPERATION_INPLACE_SINH") { return TOKEN_OPERATION_INPLACE_SINH ; }
        if( type == "TOKEN_OPERATION_TANH") { return TOKEN_OPERATION_TANH ; }
        if( type == "TOKEN_OPERATION_INPLACE_TANH") { return TOKEN_OPERATION_INPLACE_TANH ; }
        if( type == "TOKEN_OPERATION_EXP") { return TOKEN_OPERATION_EXP ; }
        if( type == "TOKEN_OPERATION_INPLACE_EXP") { return TOKEN_OPERATION_INPLACE_EXP ; }
        if( type == "TOKEN_OPERATION_SIGN") { return TOKEN_OPERATION_SIGN ; }
        if( type == "TOKEN_OPERATION_INPLACE_SIGN") { return TOKEN_OPERATION_INPLACE_SIGN ; }
        if( type == "TOKEN_OPERATION_POSITIVITY") { return TOKEN_OPERATION_POSITIVITY ; }
        if( type == "TOKEN_OPERATION_INPLACE_POSITIVITY") { return TOKEN_OPERATION_INPLACE_POSITIVITY ; }
        if( type == "TOKEN_OPERATION_NEGATIVITY") { return TOKEN_OPERATION_NEGATIVITY ; }
        if( type == "TOKEN_OPERATION_INPLACE_NEGATIVITY") { return TOKEN_OPERATION_INPLACE_NEGATIVITY ; }
        if( type == "TOKEN_OPERATION_LOG") { return TOKEN_OPERATION_LOG ; }
        if( type == "TOKEN_OPERATION_INPLACE_LOG") { return TOKEN_OPERATION_INPLACE_LOG ; }
        if( type == "TOKEN_OPERATION_SQRT") { return TOKEN_OPERATION_SQRT ; }
        if( type == "TOKEN_OPERATION_INPLACE_SQRT") { return TOKEN_OPERATION_INPLACE_SQRT ; }
        if( type == "TOKEN_OPERATION_BESSEL") { return TOKEN_OPERATION_BESSEL ; }
        if( type == "TOKEN_OPERATION_INPLACE_BESSEL") { return TOKEN_OPERATION_INPLACE_BESSEL ; }
        if( type == "TOKEN_OPERATION_MIN") { return TOKEN_OPERATION_MIN ; }
        if( type == "TOKEN_OPERATION_INPLACE_MIN") { return TOKEN_OPERATION_INPLACE_MIN ; }
        if( type == "TOKEN_OPERATION_MAX") { return TOKEN_OPERATION_MAX ; }
        if( type == "TOKEN_OPERATION_INPLACE_MAX") { return TOKEN_OPERATION_INPLACE_MAX ; }
        if( type == "TOKEN_OPERATION_ATAN2") { return TOKEN_OPERATION_ATAN2 ; }
        if( type == "TOKEN_OPERATION_INPLACE_ATAN2") { return TOKEN_OPERATION_INPLACE_ATAN2 ; }
        if( type == "TOKEN_OPERATION_INTERPOLATE") { return TOKEN_OPERATION_INTERPOLATE ; }
        if( type == "TOKEN_OPERATION_INPLACE_INTERPOLATE") { return TOKEN_OPERATION_INPLACE_INTERPOLATE ; }
        if( type == "TOKEN_OPERATION_GEO_OPERATION") { return TOKEN_OPERATION_GEO_OPERATION ; }
        if(ok) { *ok = false ; }
        return TOKEN_OPERATION_CONSTANT ;
    }
    static std::string fromTokenOperationType(TokenOperationType value)
    {
        switch(value)
        {
            case TOKEN_OPERATION_CONSTANT: return "TOKEN_OPERATION_CONSTANT" ;
            case TOKEN_OPERATION_X: return "TOKEN_OPERATION_X" ;
            case TOKEN_OPERATION_Y: return "TOKEN_OPERATION_Y" ;
            case TOKEN_OPERATION_Z: return "TOKEN_OPERATION_Z" ;
            case TOKEN_OPERATION_T: return "TOKEN_OPERATION_T" ;
            case TOKEN_OPERATION_U: return "TOKEN_OPERATION_U" ;
            case TOKEN_OPERATION_V: return "TOKEN_OPERATION_V" ;
            case TOKEN_OPERATION_W: return "TOKEN_OPERATION_W" ;
            case TOKEN_OPERATION_PLUS: return "TOKEN_OPERATION_PLUS" ;
            case TOKEN_OPERATION_INPLACE_PLUS: return "TOKEN_OPERATION_INPLACE_PLUS" ;
            case TOKEN_OPERATION_MINUS: return "TOKEN_OPERATION_MINUS" ;
            case TOKEN_OPERATION_INPLACE_MINUS: return "TOKEN_OPERATION_INPLACE_MINUS" ;
            case TOKEN_OPERATION_TIMES: return "TOKEN_OPERATION_TIMES" ;
            case TOKEN_OPERATION_INPLACE_TIMES: return "TOKEN_OPERATION_INPLACE_TIMES" ;
            case TOKEN_OPERATION_DIVIDES: return "TOKEN_OPERATION_DIVIDES" ;
            case TOKEN_OPERATION_INPLACE_DIVIDES: return "TOKEN_OPERATION_INPLACE_DIVIDES" ;
            case TOKEN_OPERATION_POWER: return "TOKEN_OPERATION_POWER" ;
            case TOKEN_OPERATION_INPLACE_POWER: return "TOKEN_OPERATION_INPLACE_POWER" ;
            case TOKEN_OPERATION_ABS: return "TOKEN_OPERATION_ABS" ;
            case TOKEN_OPERATION_INPLACE_ABS: return "TOKEN_OPERATION_INPLACE_ABS" ;
            case TOKEN_OPERATION_COS: return "TOKEN_OPERATION_COS" ;
            case TOKEN_OPERATION_INPLACE_COS: return "TOKEN_OPERATION_INPLACE_COS" ;
            case TOKEN_OPERATION_SIN: return "TOKEN_OPERATION_SIN" ;
            case TOKEN_OPERATION_INPLACE_SIN: return "TOKEN_OPERATION_INPLACE_SIN" ;
            case TOKEN_OPERATION_TAN: return "TOKEN_OPERATION_TAN" ;
            case TOKEN_OPERATION_INPLACE_TAN: return "TOKEN_OPERATION_INPLACE_TAN" ;
            case TOKEN_OPERATION_COSH: return "TOKEN_OPERATION_COSH" ;
            case TOKEN_OPERATION_INPLACE_COSH: return "TOKEN_OPERATION_INPLACE_COSH" ;
            case TOKEN_OPERATION_SINH: return "TOKEN_OPERATION_SINH" ;
            case TOKEN_OPERATION_INPLACE_SINH: return "TOKEN_OPERATION_INPLACE_SINH" ;
            case TOKEN_OPERATION_TANH: return "TOKEN_OPERATION_TANH" ;
            case TOKEN_OPERATION_INPLACE_TANH: return "TOKEN_OPERATION_INPLACE_TANH" ;
            case TOKEN_OPERATION_EXP: return "TOKEN_OPERATION_EXP" ;
            case TOKEN_OPERATION_INPLACE_EXP: return "TOKEN_OPERATION_INPLACE_EXP" ;
            case TOKEN_OPERATION_SIGN: return "TOKEN_OPERATION_SIGN" ;
            case TOKEN_OPERATION_INPLACE_SIGN: return "TOKEN_OPERATION_INPLACE_SIGN" ;
            case TOKEN_OPERATION_POSITIVITY: return "TOKEN_OPERATION_POSITIVITY" ;
            case TOKEN_OPERATION_INPLACE_POSITIVITY: return "TOKEN_OPERATION_INPLACE_POSITIVITY" ;
            case TOKEN_OPERATION_NEGATIVITY: return "TOKEN_OPERATION_NEGATIVITY" ;
            case TOKEN_OPERATION_INPLACE_NEGATIVITY: return "TOKEN_OPERATION_INPLACE_NEGATIVITY" ;
            case TOKEN_OPERATION_LOG: return "TOKEN_OPERATION_LOG" ;
            case TOKEN_OPERATION_INPLACE_LOG: return "TOKEN_OPERATION_INPLACE_LOG" ;
            case TOKEN_OPERATION_SQRT: return "TOKEN_OPERATION_SQRT" ;
            case TOKEN_OPERATION_INPLACE_SQRT: return "TOKEN_OPERATION_INPLACE_SQRT" ;
            case TOKEN_OPERATION_BESSEL: return "TOKEN_OPERATION_BESSEL" ;
            case TOKEN_OPERATION_INPLACE_BESSEL: return "TOKEN_OPERATION_INPLACE_BESSEL" ;
            case TOKEN_OPERATION_MIN: return "TOKEN_OPERATION_MIN" ;
            case TOKEN_OPERATION_INPLACE_MIN: return "TOKEN_OPERATION_INPLACE_MIN" ;
            case TOKEN_OPERATION_MAX: return "TOKEN_OPERATION_MAX" ;
            case TOKEN_OPERATION_INPLACE_MAX: return "TOKEN_OPERATION_INPLACE_MAX" ;
            case TOKEN_OPERATION_ATAN2: return "TOKEN_OPERATION_ATAN2" ;
            case TOKEN_OPERATION_INPLACE_ATAN2: return "TOKEN_OPERATION_INPLACE_ATAN2" ;
            case TOKEN_OPERATION_INTERPOLATE: return "TOKEN_OPERATION_INTERPOLATE" ;
            case TOKEN_OPERATION_INPLACE_INTERPOLATE: return "TOKEN_OPERATION_INPLACE_INTERPOLATE" ;
            case TOKEN_OPERATION_GEO_OPERATION: return "TOKEN_OPERATION_GEO_OPERATION" ;
        }
        return "TOKEN_OPERATION_CONSTANT" ;
    }
   
    // parsed from header file: ../polynomial/variable.h
    static Variable getVariable(std::string type, bool * ok = 0)
    {
        if(ok) { *ok = true ; }
        if( type == "ONE") { return ONE ; }
        if( type == "XI") { return XI ; }
        if( type == "ETA") { return ETA ; }
        if( type == "ZETA") { return ZETA ; }
        if( type == "TIME_VARIABLE") { return TIME_VARIABLE ; }
        if( type == "U_VARIABLE") { return U_VARIABLE ; }
        if( type == "V_VARIABLE") { return V_VARIABLE ; }
        if( type == "W_VARIABLE") { return W_VARIABLE ; }
        if(ok) { *ok = false ; }
        return ONE ;
    }
    static std::string fromVariable(Variable value)
    {
        switch(value)
        {
            case ONE: return "ONE" ;
            case XI: return "XI" ;
            case ETA: return "ETA" ;
            case ZETA: return "ZETA" ;
            case TIME_VARIABLE: return "TIME_VARIABLE" ;
            case U_VARIABLE: return "U_VARIABLE" ;
            case V_VARIABLE: return "V_VARIABLE" ;
            case W_VARIABLE: return "W_VARIABLE" ;
        }
        return "ONE" ;
    }
   
    // parsed from header file: ../polynomial/typeref.h
    static TypeRef getTypeRef(std::string type, bool * ok = 0)
    {
        if(ok) { *ok = true ; }
        if( type == "ZERO_TO_ONE") { return ZERO_TO_ONE ; }
        if( type == "MIN_ONE_TO_ONE") { return MIN_ONE_TO_ONE ; }
        if(ok) { *ok = false ; }
        return ZERO_TO_ONE ;
    }
    static std::string fromTypeRef(TypeRef value)
    {
        switch(value)
        {
            case ZERO_TO_ONE: return "ZERO_TO_ONE" ;
            case MIN_ONE_TO_ONE: return "MIN_ONE_TO_ONE" ;
        }
        return "ZERO_TO_ONE" ;
    }
   
    // parsed from header file: ../polynomial/vm_function_base.h
    static PositionTokenType getPositionTokenType(std::string type, bool * ok = 0)
    {
        if(ok) { *ok = true ; }
        if( type == "POSITION_TOKEN") { return POSITION_TOKEN ; }
        if( type == "PROJECTION_TOKEN") { return PROJECTION_TOKEN ; }
        if(ok) { *ok = false ; }
        return POSITION_TOKEN ;
    }
    static std::string fromPositionTokenType(PositionTokenType value)
    {
        switch(value)
        {
            case POSITION_TOKEN: return "POSITION_TOKEN" ;
            case PROJECTION_TOKEN: return "PROJECTION_TOKEN" ;
        }
        return "POSITION_TOKEN" ;
    }
   
    // parsed from header file: ../polynomial/vm_function_base.h
    static TemporayUsageType getTemporayUsageType(std::string type, bool * ok = 0)
    {
        if(ok) { *ok = true ; }
        if( type == "NO_TEMPORARY") { return NO_TEMPORARY ; }
        if( type == "SET_TEMPORARY") { return SET_TEMPORARY ; }
        if( type == "SET_GET_TEMPORARY_A") { return SET_GET_TEMPORARY_A ; }
        if( type == "SET_GET_TEMPORARY_B") { return SET_GET_TEMPORARY_B ; }
        if( type == "GET_TEMPORARY_A") { return GET_TEMPORARY_A ; }
        if( type == "GET_TEMPORARY_B") { return GET_TEMPORARY_B ; }
        if(ok) { *ok = false ; }
        return NO_TEMPORARY ;
    }
    static std::string fromTemporayUsageType(TemporayUsageType value)
    {
        switch(value)
        {
            case NO_TEMPORARY: return "NO_TEMPORARY" ;
            case SET_TEMPORARY: return "SET_TEMPORARY" ;
            case SET_GET_TEMPORARY_A: return "SET_GET_TEMPORARY_A" ;
            case SET_GET_TEMPORARY_B: return "SET_GET_TEMPORARY_B" ;
            case GET_TEMPORARY_A: return "GET_TEMPORARY_A" ;
            case GET_TEMPORARY_B: return "GET_TEMPORARY_B" ;
        }
        return "NO_TEMPORARY" ;
    }
   
    // parsed from header file: ../solvers/assembly.h
    static LagrangeMultiplierType getLagrangeMultiplierType(std::string type, bool * ok = 0)
    {
        if(ok) { *ok = true ; }
        if( type == "GENERAL") { return GENERAL ; }
        if( type == "FIX_ALONG_ALL") { return FIX_ALONG_ALL ; }
        if( type == "FIX_ALONG_XI") { return FIX_ALONG_XI ; }
        if( type == "SET_ALONG_XI") { return SET_ALONG_XI ; }
        if( type == "INCREMENT_ALONG_XI") { return INCREMENT_ALONG_XI ; }
        if( type == "FIX_ALONG_ETA") { return FIX_ALONG_ETA ; }
        if( type == "SET_ALONG_ETA") { return SET_ALONG_ETA ; }
        if( type == "INCREMENT_ALONG_ETA") { return INCREMENT_ALONG_ETA ; }
        if( type == "FIX_ALONG_ZETA") { return FIX_ALONG_ZETA ; }
        if( type == "SET_ALONG_ZETA") { return SET_ALONG_ZETA ; }
        if( type == "INCREMENT_ALONG_ZETA") { return INCREMENT_ALONG_ZETA ; }
        if( type == "FIX_ALONG_XI_ETA") { return FIX_ALONG_XI_ETA ; }
        if( type == "SET_ALONG_XI_ETA") { return SET_ALONG_XI_ETA ; }
        if( type == "INCREMENT_ALONG_XI_ETA") { return INCREMENT_ALONG_XI_ETA ; }
        if( type == "FIX_ALONG_XI_ZETA") { return FIX_ALONG_XI_ZETA ; }
        if( type == "SET_ALONG_XI_ZETA") { return SET_ALONG_XI_ZETA ; }
        if( type == "INCREMENT_ALONG_XI_ZETA") { return INCREMENT_ALONG_XI_ZETA ; }
        if( type == "FIX_ALONG_ETA_ZETA") { return FIX_ALONG_ETA_ZETA ; }
        if( type == "SET_ALONG_ETA_ZETA") { return SET_ALONG_ETA_ZETA ; }
        if( type == "INCREMENT_ALONG_ETA_ZETA") { return INCREMENT_ALONG_ETA_ZETA ; }
        if( type == "FIX_ALONG_INDEXED_AXIS") { return FIX_ALONG_INDEXED_AXIS ; }
        if( type == "SET_ALONG_INDEXED_AXIS") { return SET_ALONG_INDEXED_AXIS ; }
        if( type == "INCREMENT_ALONG_INDEXED_AXIS") { return INCREMENT_ALONG_INDEXED_AXIS ; }
        if( type == "SET_PROPORTIONAL_DISPLACEMENT") { return SET_PROPORTIONAL_DISPLACEMENT ; }
        if( type == "SET_PROPORTIONAL_DISPLACEMENT_XI_ETA") { return SET_PROPORTIONAL_DISPLACEMENT_XI_ETA ; }
        if( type == "SET_PROPORTIONAL_DISPLACEMENT_XI_ZETA") { return SET_PROPORTIONAL_DISPLACEMENT_XI_ZETA ; }
        if( type == "SET_PROPORTIONAL_DISPLACEMENT_ETA_XI") { return SET_PROPORTIONAL_DISPLACEMENT_ETA_XI ; }
        if( type == "SET_PROPORTIONAL_DISPLACEMENT_ETA_ZETA") { return SET_PROPORTIONAL_DISPLACEMENT_ETA_ZETA ; }
        if( type == "SET_PROPORTIONAL_DISPLACEMENT_ZETA_XI") { return SET_PROPORTIONAL_DISPLACEMENT_ZETA_XI ; }
        if( type == "SET_PROPORTIONAL_DISPLACEMENT_ZETA_ETA") { return SET_PROPORTIONAL_DISPLACEMENT_ZETA_ETA ; }
        if( type == "FIX_NORMAL_DISPLACEMENT") { return FIX_NORMAL_DISPLACEMENT ; }
        if( type == "FIX_TANGENT_DISPLACEMENT") { return FIX_TANGENT_DISPLACEMENT ; }
        if( type == "SET_NORMAL_DISPLACEMENT") { return SET_NORMAL_DISPLACEMENT ; }
        if( type == "SET_TANGENT_DISPLACEMENT") { return SET_TANGENT_DISPLACEMENT ; }
        if( type == "SET_FORCE_XI") { return SET_FORCE_XI ; }
        if( type == "SET_FORCE_ETA") { return SET_FORCE_ETA ; }
        if( type == "SET_FORCE_ZETA") { return SET_FORCE_ZETA ; }
        if( type == "SET_FORCE_INDEXED_AXIS") { return SET_FORCE_INDEXED_AXIS ; }
        if( type == "INCREMENT_FORCE_INDEXED_AXIS") { return INCREMENT_FORCE_INDEXED_AXIS ; }
        if( type == "SET_FLUX_XI") { return SET_FLUX_XI ; }
        if( type == "SET_FLUX_ETA") { return SET_FLUX_ETA ; }
        if( type == "SET_FLUX_ZETA") { return SET_FLUX_ZETA ; }
        if( type == "SET_VOLUMIC_STRESS_XI") { return SET_VOLUMIC_STRESS_XI ; }
        if( type == "SET_VOLUMIC_STRESS_XI_ETA") { return SET_VOLUMIC_STRESS_XI_ETA ; }
        if( type == "SET_VOLUMIC_STRESS_ETA") { return SET_VOLUMIC_STRESS_ETA ; }
        if( type == "SET_VOLUMIC_STRESS_ZETA") { return SET_VOLUMIC_STRESS_ZETA ; }
        if( type == "SET_VOLUMIC_STRESS_XI_ZETA") { return SET_VOLUMIC_STRESS_XI_ZETA ; }
        if( type == "SET_VOLUMIC_STRESS_ETA_ZETA") { return SET_VOLUMIC_STRESS_ETA_ZETA ; }
        if( type == "SET_STRESS_XI") { return SET_STRESS_XI ; }
        if( type == "SET_STRESS_ETA") { return SET_STRESS_ETA ; }
        if( type == "SET_STRESS_ZETA") { return SET_STRESS_ZETA ; }
        if( type == "SET_NORMAL_STRESS") { return SET_NORMAL_STRESS ; }
        if( type == "SET_TANGENT_STRESS") { return SET_TANGENT_STRESS ; }
        if( type == "VERTICAL_PLANE_SECTIONS") { return VERTICAL_PLANE_SECTIONS ; }
        if( type == "HORIZONTAL_PLANE_SECTIONS") { return HORIZONTAL_PLANE_SECTIONS ; }
        if( type == "nullptr_CONDITION") { return nullptr_CONDITION ; }
        if( type == "SET_GLOBAL_FORCE_VECTOR") { return SET_GLOBAL_FORCE_VECTOR ; }
        if(ok) { *ok = false ; }
        return GENERAL ;
    }
    static std::string fromLagrangeMultiplierType(LagrangeMultiplierType value)
    {
        switch(value)
        {
            case GENERAL: return "GENERAL" ;
            case FIX_ALONG_ALL: return "FIX_ALONG_ALL" ;
            case FIX_ALONG_XI: return "FIX_ALONG_XI" ;
            case SET_ALONG_XI: return "SET_ALONG_XI" ;
            case INCREMENT_ALONG_XI: return "INCREMENT_ALONG_XI" ;
            case FIX_ALONG_ETA: return "FIX_ALONG_ETA" ;
            case SET_ALONG_ETA: return "SET_ALONG_ETA" ;
            case INCREMENT_ALONG_ETA: return "INCREMENT_ALONG_ETA" ;
            case FIX_ALONG_ZETA: return "FIX_ALONG_ZETA" ;
            case SET_ALONG_ZETA: return "SET_ALONG_ZETA" ;
            case INCREMENT_ALONG_ZETA: return "INCREMENT_ALONG_ZETA" ;
            case FIX_ALONG_XI_ETA: return "FIX_ALONG_XI_ETA" ;
            case SET_ALONG_XI_ETA: return "SET_ALONG_XI_ETA" ;
            case INCREMENT_ALONG_XI_ETA: return "INCREMENT_ALONG_XI_ETA" ;
            case FIX_ALONG_XI_ZETA: return "FIX_ALONG_XI_ZETA" ;
            case SET_ALONG_XI_ZETA: return "SET_ALONG_XI_ZETA" ;
            case INCREMENT_ALONG_XI_ZETA: return "INCREMENT_ALONG_XI_ZETA" ;
            case FIX_ALONG_ETA_ZETA: return "FIX_ALONG_ETA_ZETA" ;
            case SET_ALONG_ETA_ZETA: return "SET_ALONG_ETA_ZETA" ;
            case INCREMENT_ALONG_ETA_ZETA: return "INCREMENT_ALONG_ETA_ZETA" ;
            case FIX_ALONG_INDEXED_AXIS: return "FIX_ALONG_INDEXED_AXIS" ;
            case SET_ALONG_INDEXED_AXIS: return "SET_ALONG_INDEXED_AXIS" ;
            case INCREMENT_ALONG_INDEXED_AXIS: return "INCREMENT_ALONG_INDEXED_AXIS" ;
            case SET_PROPORTIONAL_DISPLACEMENT: return "SET_PROPORTIONAL_DISPLACEMENT" ;
            case SET_PROPORTIONAL_DISPLACEMENT_XI_ETA: return "SET_PROPORTIONAL_DISPLACEMENT_XI_ETA" ;
            case SET_PROPORTIONAL_DISPLACEMENT_XI_ZETA: return "SET_PROPORTIONAL_DISPLACEMENT_XI_ZETA" ;
            case SET_PROPORTIONAL_DISPLACEMENT_ETA_XI: return "SET_PROPORTIONAL_DISPLACEMENT_ETA_XI" ;
            case SET_PROPORTIONAL_DISPLACEMENT_ETA_ZETA: return "SET_PROPORTIONAL_DISPLACEMENT_ETA_ZETA" ;
            case SET_PROPORTIONAL_DISPLACEMENT_ZETA_XI: return "SET_PROPORTIONAL_DISPLACEMENT_ZETA_XI" ;
            case SET_PROPORTIONAL_DISPLACEMENT_ZETA_ETA: return "SET_PROPORTIONAL_DISPLACEMENT_ZETA_ETA" ;
            case FIX_NORMAL_DISPLACEMENT: return "FIX_NORMAL_DISPLACEMENT" ;
            case FIX_TANGENT_DISPLACEMENT: return "FIX_TANGENT_DISPLACEMENT" ;
            case SET_NORMAL_DISPLACEMENT: return "SET_NORMAL_DISPLACEMENT" ;
            case SET_TANGENT_DISPLACEMENT: return "SET_TANGENT_DISPLACEMENT" ;
            case SET_FORCE_XI: return "SET_FORCE_XI" ;
            case SET_FORCE_ETA: return "SET_FORCE_ETA" ;
            case SET_FORCE_ZETA: return "SET_FORCE_ZETA" ;
            case SET_FORCE_INDEXED_AXIS: return "SET_FORCE_INDEXED_AXIS" ;
            case INCREMENT_FORCE_INDEXED_AXIS: return "INCREMENT_FORCE_INDEXED_AXIS" ;
            case SET_FLUX_XI: return "SET_FLUX_XI" ;
            case SET_FLUX_ETA: return "SET_FLUX_ETA" ;
            case SET_FLUX_ZETA: return "SET_FLUX_ZETA" ;
            case SET_VOLUMIC_STRESS_XI: return "SET_VOLUMIC_STRESS_XI" ;
            case SET_VOLUMIC_STRESS_XI_ETA: return "SET_VOLUMIC_STRESS_XI_ETA" ;
            case SET_VOLUMIC_STRESS_ETA: return "SET_VOLUMIC_STRESS_ETA" ;
            case SET_VOLUMIC_STRESS_ZETA: return "SET_VOLUMIC_STRESS_ZETA" ;
            case SET_VOLUMIC_STRESS_XI_ZETA: return "SET_VOLUMIC_STRESS_XI_ZETA" ;
            case SET_VOLUMIC_STRESS_ETA_ZETA: return "SET_VOLUMIC_STRESS_ETA_ZETA" ;
            case SET_STRESS_XI: return "SET_STRESS_XI" ;
            case SET_STRESS_ETA: return "SET_STRESS_ETA" ;
            case SET_STRESS_ZETA: return "SET_STRESS_ZETA" ;
            case SET_NORMAL_STRESS: return "SET_NORMAL_STRESS" ;
            case SET_TANGENT_STRESS: return "SET_TANGENT_STRESS" ;
            case VERTICAL_PLANE_SECTIONS: return "VERTICAL_PLANE_SECTIONS" ;
            case HORIZONTAL_PLANE_SECTIONS: return "HORIZONTAL_PLANE_SECTIONS" ;
            case nullptr_CONDITION: return "nullptr_CONDITION" ;
            case SET_GLOBAL_FORCE_VECTOR: return "SET_GLOBAL_FORCE_VECTOR" ;
        }
        return "GENERAL" ;
    }
   
    // parsed from header file: ../utilities/configuration.h
    static bool getbool(std::string type, bool * ok = 0)
    {
        if(ok) { *ok = true ; }
        if( type == "false") { return false ; }
        if( type == "true") { return true ; }
        if(ok) { *ok = false ; }
        return false ;
    }
    static std::string frombool(bool value) { return ( value ? "true" : "false" ) ; }
   
    // parsed from header file: ../utilities/writer/voxel_writer.h
    static VWFieldType getVWFieldType(std::string type, bool * ok = 0)
    {
        if(ok) { *ok = true ; }
        if( type == "VWFT_PRINCIPAL_ANGLE") { return VWFT_PRINCIPAL_ANGLE ; }
        if( type == "VWFT_STIFFNESS") { return VWFT_STIFFNESS ; }
        if( type == "VWFT_STRAIN") { return VWFT_STRAIN ; }
        if( type == "VWFT_STRESS") { return VWFT_STRESS ; }
        if( type == "VWFT_PRINCIPAL_STRAIN") { return VWFT_PRINCIPAL_STRAIN ; }
        if( type == "VWFT_PRINCIPAL_STRESS") { return VWFT_PRINCIPAL_STRESS ; }
        if( type == "VWFT_STRAIN_AND_STRESS") { return VWFT_STRAIN_AND_STRESS ; }
        if( type == "VWFT_CONCENTRATION") { return VWFT_CONCENTRATION ; }
        if( type == "VWFT_GRADIENT") { return VWFT_GRADIENT ; }
        if( type == "VWFT_FLUX") { return VWFT_FLUX ; }
        if( type == "VWFT_GRADIENT_AND_FLUX") { return VWFT_GRADIENT_AND_FLUX ; }
        if( type == "VWFT_VON_MISES") { return VWFT_VON_MISES ; }
        if( type == "VWFT_ENRICHEMENT") { return VWFT_ENRICHEMENT ; }
        if( type == "VWFT_DAMAGE") { return VWFT_DAMAGE ; }
        if(ok) { *ok = false ; }
        return VWFT_PRINCIPAL_ANGLE ;
    }
    static std::string fromVWFieldType(VWFieldType value)
    {
        switch(value)
        {
            case VWFT_PRINCIPAL_ANGLE: return "VWFT_PRINCIPAL_ANGLE" ;
            case VWFT_STIFFNESS: return "VWFT_STIFFNESS" ;
            case VWFT_STRAIN: return "VWFT_STRAIN" ;
            case VWFT_STRESS: return "VWFT_STRESS" ;
            case VWFT_PRINCIPAL_STRAIN: return "VWFT_PRINCIPAL_STRAIN" ;
            case VWFT_PRINCIPAL_STRESS: return "VWFT_PRINCIPAL_STRESS" ;
            case VWFT_STRAIN_AND_STRESS: return "VWFT_STRAIN_AND_STRESS" ;
            case VWFT_CONCENTRATION: return "VWFT_CONCENTRATION" ;
            case VWFT_GRADIENT: return "VWFT_GRADIENT" ;
            case VWFT_FLUX: return "VWFT_FLUX" ;
            case VWFT_GRADIENT_AND_FLUX: return "VWFT_GRADIENT_AND_FLUX" ;
            case VWFT_VON_MISES: return "VWFT_VON_MISES" ;
            case VWFT_ENRICHEMENT: return "VWFT_ENRICHEMENT" ;
            case VWFT_DAMAGE: return "VWFT_DAMAGE" ;
        }
        return "VWFT_PRINCIPAL_ANGLE" ;
    }
   
    // parsed from header file: ../utilities/writer/triangle_writer.h
    static TWFieldType getTWFieldType(std::string type, bool * ok = 0)
    {
        if(ok) { *ok = true ; }
        if( type == "TWFT_COORDINATE") { return TWFT_COORDINATE ; }
        if( type == "TWFT_DISPLACEMENTS") { return TWFT_DISPLACEMENTS ; }
        if( type == "TWFT_SCALAR") { return TWFT_SCALAR ; }
        if( type == "TWFT_DOH") { return TWFT_DOH ; }
        if( type == "TWFT_PRINCIPAL_ANGLE") { return TWFT_PRINCIPAL_ANGLE ; }
        if( type == "TWFT_TRIANGLE_ANGLE") { return TWFT_TRIANGLE_ANGLE ; }
        if( type == "TWFT_CRACK_ANGLE") { return TWFT_CRACK_ANGLE ; }
        if( type == "TWFT_CRITERION") { return TWFT_CRITERION ; }
        if( type == "TWFT_STIFFNESS") { return TWFT_STIFFNESS ; }
        if( type == "TWFT_SPIN") { return TWFT_LARGE_DEFORMATION_TRANSFORM ; }
        if( type == "TWFT_VISCOSITY") { return TWFT_VISCOSITY ; }
        if( type == "TWFT_STIFFNESS_X") { return TWFT_STIFFNESS_X ; }
        if( type == "TWFT_STIFFNESS_Y") { return TWFT_STIFFNESS_Y ; }
        if( type == "TWFT_STIFFNESS_Z") { return TWFT_STIFFNESS_Z ; }
        if( type == "TWFT_ENRICHMENT") { return TWFT_ENRICHMENT ; }
        if( type == "TWFT_IMPOSED_STRESS_NORM") { return TWFT_IMPOSED_STRESS_NORM ; }
        if( type == "TWFT_PRINCIPAL_STRESS") { return TWFT_PRINCIPAL_STRESS ; }
        if( type == "TWFT_PRINCIPAL_STRAIN") { return TWFT_PRINCIPAL_STRAIN ; }
        if( type == "TWFT_DAMAGE") { return TWFT_DAMAGE ; }
        if( type == "TWFT_GRADIENT") { return TWFT_GRADIENT ; }
        if( type == "TWFT_FLUX") { return TWFT_FLUX ; }
        if( type == "TWFT_GRADIENT_AND_FLUX") { return TWFT_GRADIENT_AND_FLUX ; }
        if( type == "TWFT_VON_MISES") { return TWFT_VON_MISES ; }
        if( type == "TWFT_CRACKS") { return TWFT_CRACKS ; }
        if( type == "TWFT_INTERSECTION") { return TWFT_INTERSECTION ; }
        if( type == "TWFT_FIELD_TYPE") { return TWFT_FIELD_TYPE ; }
        if( type == "TWFT_INTERNAL_VARIABLE") { return TWFT_INTERNAL_VARIABLE ; }
        if(ok) { *ok = false ; }
        return TWFT_COORDINATE ;
    }
    static std::string fromTWFieldType(TWFieldType value)
    {
        switch(value)
        {
            case TWFT_COORDINATE: return "TWFT_COORDINATE" ;
            case TWFT_LARGE_DEFORMATION_TRANSFORM: return "TWFT_LARGE_DEFORMATION_TRANSFORM" ;
            case TWFT_LARGE_DEFORMATION_ANGLE: return "TWFT_LARGE_DEFORMATION_ANGLE" ;
            case TWFT_DISPLACEMENTS: return "TWFT_DISPLACEMENTS" ;
            case TWFT_SCALAR: return "TWFT_SCALAR" ;
            case TWFT_DOH: return "TWFT_DOH" ;
            case TWFT_PRINCIPAL_ANGLE: return "TWFT_PRINCIPAL_ANGLE" ;
            case TWFT_TRIANGLE_ANGLE: return "TWFT_TRIANGLE_ANGLE" ;
            case TWFT_CRACK_ANGLE: return "TWFT_CRACK_ANGLE" ;
            case TWFT_CRITERION: return "TWFT_CRITERION" ;
            case TWFT_STIFFNESS: return "TWFT_STIFFNESS" ;
            case TWFT_VISCOSITY: return "TWFT_VISCOSITY" ;
            case TWFT_STIFFNESS_X: return "TWFT_STIFFNESS_X" ;
            case TWFT_STIFFNESS_Y: return "TWFT_STIFFNESS_Y" ;
            case TWFT_STIFFNESS_Z: return "TWFT_STIFFNESS_Z" ;
            case TWFT_ENRICHMENT: return "TWFT_ENRICHMENT" ;
            case TWFT_IMPOSED_STRESS_NORM: return "TWFT_IMPOSED_STRESS_NORM" ;
            case TWFT_PRINCIPAL_STRESS: return "TWFT_PRINCIPAL_STRESS" ;
            case TWFT_PRINCIPAL_STRAIN: return "TWFT_PRINCIPAL_STRAIN" ;
            case TWFT_DAMAGE: return "TWFT_DAMAGE" ;
            case TWFT_GRADIENT: return "TWFT_GRADIENT" ;
            case TWFT_FLUX: return "TWFT_FLUX" ;
            case TWFT_GRADIENT_AND_FLUX: return "TWFT_GRADIENT_AND_FLUX" ;
            case TWFT_VON_MISES: return "TWFT_VON_MISES" ;
            case TWFT_CRACKS: return "TWFT_CRACKS" ;
            case TWFT_INTERSECTION: return "TWFT_INTERSECTION" ;
            case TWFT_FIELD_TYPE: return "TWFT_FIELD_TYPE" ;
            case TWFT_INTERNAL_VARIABLE: return "TWFT_INTERNAL_VARIABLE" ;
        }
        return "TWFT_COORDINATE" ;
    }
   
    // parsed from header file: ../utilities/parser/config_parser.h
    static XMLTokenType getXMLTokenType(std::string type, bool * ok = 0)
    {
        if(ok) { *ok = true ; }
        if( type == "XML_OPEN_COMMENT") { return XML_OPEN_COMMENT ; }
        if( type == "XML_CLOSE_COMMENT") { return XML_CLOSE_COMMENT ; }
        if( type == "XML_INLINE_COMMENT") { return XML_INLINE_COMMENT ; }
        if( type == "XML_OPEN") { return XML_OPEN ; }
        if( type == "XML_CLOSE") { return XML_CLOSE ; }
        if( type == "XML_INLINE") { return XML_INLINE ; }
        if( type == "XML_VALUE") { return XML_VALUE ; }
        if( type == "XML_INVALID") { return XML_INVALID ; }
        if( type == "XML_HEADER") { return XML_HEADER ; }
        if(ok) { *ok = false ; }
        return XML_OPEN_COMMENT ;
    }
    static std::string fromXMLTokenType(XMLTokenType value)
    {
        switch(value)
        {
            case XML_OPEN_COMMENT: return "XML_OPEN_COMMENT" ;
            case XML_CLOSE_COMMENT: return "XML_CLOSE_COMMENT" ;
            case XML_INLINE_COMMENT: return "XML_INLINE_COMMENT" ;
            case XML_OPEN: return "XML_OPEN" ;
            case XML_CLOSE: return "XML_CLOSE" ;
            case XML_INLINE: return "XML_INLINE" ;
            case XML_VALUE: return "XML_VALUE" ;
            case XML_INVALID: return "XML_INVALID" ;
            case XML_HEADER: return "XML_HEADER" ;
        }
        return "XML_OPEN_COMMENT" ;
    }
   
    // parsed from header file: ../utilities/granulo.h
    static TypeInclusion getTypeInclusion(std::string type, bool * ok = 0)
    {
        if(ok) { *ok = true ; }
        if( type == "CIRCLE_INCLUSION") { return CIRCLE_INCLUSION ; }
        if( type == "SPHERE_INCLUSION") { return SPHERE_INCLUSION ; }
        if( type == "ELLIPSE_INCLUSION") { return ELLIPSE_INCLUSION ; }
        if(ok) { *ok = false ; }
        return CIRCLE_INCLUSION ;
    }
    static std::string fromTypeInclusion(TypeInclusion value)
    {
        switch(value)
        {
            case CIRCLE_INCLUSION: return "CIRCLE_INCLUSION" ;
            case SPHERE_INCLUSION: return "SPHERE_INCLUSION" ;
            case ELLIPSE_INCLUSION: return "ELLIPSE_INCLUSION" ;
        }
        return "CIRCLE_INCLUSION" ;
    }
   
    // parsed from header file: ../utilities/granulo.h
    static PSDSpecificationType getPSDSpecificationType(std::string type, bool * ok = 0)
    {
        if(ok) { *ok = true ; }
        if( type == "CUMULATIVE_PERCENT") { return CUMULATIVE_PERCENT ; }
        if( type == "CUMULATIVE_FRACTION") { return CUMULATIVE_FRACTION ; }
        if( type == "CUMULATIVE_ABSOLUTE") { return CUMULATIVE_ABSOLUTE ; }
        if( type == "CUMULATIVE_PERCENT_REVERSE") { return CUMULATIVE_PERCENT_REVERSE ; }
        if( type == "CUMULATIVE_FRACTION_REVERSE") { return CUMULATIVE_FRACTION_REVERSE ; }
        if( type == "CUMULATIVE_ABSOLUTE_REVERSE") { return CUMULATIVE_ABSOLUTE_REVERSE ; }
        if(ok) { *ok = false ; }
        return CUMULATIVE_PERCENT ;
    }
    static std::string fromPSDSpecificationType(PSDSpecificationType value)
    {
        switch(value)
        {
            case CUMULATIVE_PERCENT: return "CUMULATIVE_PERCENT" ;
            case CUMULATIVE_FRACTION: return "CUMULATIVE_FRACTION" ;
            case CUMULATIVE_ABSOLUTE: return "CUMULATIVE_ABSOLUTE" ;
            case CUMULATIVE_PERCENT_REVERSE: return "CUMULATIVE_PERCENT_REVERSE" ;
            case CUMULATIVE_FRACTION_REVERSE: return "CUMULATIVE_FRACTION_REVERSE" ;
            case CUMULATIVE_ABSOLUTE_REVERSE: return "CUMULATIVE_ABSOLUTE_REVERSE" ;
        }
        return "CUMULATIVE_PERCENT" ;
    }
   
    // parsed from header file: ../utilities/granulo.h
    static VoronoiWeight getVoronoiWeight(std::string type, bool * ok = 0)
    {
        if(ok) { *ok = true ; }
        if( type == "VORONOI_CONSTANT") { return VORONOI_CONSTANT ; }
        if( type == "VORONOI_RADIUS") { return VORONOI_RADIUS ; }
        if( type == "VORONOI_RANDOM") { return VORONOI_RANDOM ; }
        if(ok) { *ok = false ; }
        return VORONOI_CONSTANT ;
    }
    static std::string fromVoronoiWeight(VoronoiWeight value)
    {
        switch(value)
        {
            case VORONOI_CONSTANT: return "VORONOI_CONSTANT" ;
            case VORONOI_RADIUS: return "VORONOI_RADIUS" ;
            case VORONOI_RANDOM: return "VORONOI_RANDOM" ;
        }
        return "VORONOI_CONSTANT" ;
    }
   
    // parsed from header file: ../utilities/tensor.h
    static planeType getplaneType(std::string type, bool * ok = 0)
    {
        if(ok) { *ok = true ; }
        if( type == "PLANE_STRESS") { return PLANE_STRESS ; }
        if( type == "PLANE_STRAIN") { return PLANE_STRAIN ; }
        if( type == "PLANE_STRESS_FREE_G") { return PLANE_STRESS_FREE_G ; }
        if(ok) { *ok = false ; }
        return PLANE_STRESS ;
    }
    static std::string fromplaneType(planeType value)
    {
        switch(value)
        {
            case PLANE_STRESS: return "PLANE_STRESS" ;
            case PLANE_STRAIN: return "PLANE_STRAIN" ;
            case PLANE_STRESS_FREE_G: return "PLANE_STRESS_FREE_G" ;
        }
        return "PLANE_STRESS" ;
    }
   
    // parsed from header file: ../utilities/tensor.h
    static SymmetryType getSymmetryType(std::string type, bool * ok = 0)
    {
        if(ok) { *ok = true ; }
        if( type == "SYMMETRY_CUBIC") { return SYMMETRY_CUBIC ; }
        if( type == "SYMMETRY_HEXAGONAL") { return SYMMETRY_HEXAGONAL ; }
        if( type == "SYMMETRY_MONOCLINIC") { return SYMMETRY_MONOCLINIC ; }
        if( type == "SYMMETRY_ORTHORHOMBIC") { return SYMMETRY_ORTHORHOMBIC ; }
        if( type == "SYMMETRY_TETRAGONAL") { return SYMMETRY_TETRAGONAL ; }
        if( type == "SYMMETRY_TRIGONAL") { return SYMMETRY_TRIGONAL ; }
        if( type == "SYMMETRY_TRICLINIC") { return SYMMETRY_TRICLINIC ; }
        if(ok) { *ok = false ; }
        return SYMMETRY_CUBIC ;
    }
    static std::string fromSymmetryType(SymmetryType value)
    {
        switch(value)
        {
            case SYMMETRY_CUBIC: return "SYMMETRY_CUBIC" ;
            case SYMMETRY_HEXAGONAL: return "SYMMETRY_HEXAGONAL" ;
            case SYMMETRY_MONOCLINIC: return "SYMMETRY_MONOCLINIC" ;
            case SYMMETRY_ORTHORHOMBIC: return "SYMMETRY_ORTHORHOMBIC" ;
            case SYMMETRY_TETRAGONAL: return "SYMMETRY_TETRAGONAL" ;
            case SYMMETRY_TRIGONAL: return "SYMMETRY_TRIGONAL" ;
            case SYMMETRY_TRICLINIC: return "SYMMETRY_TRICLINIC" ;
        }
        return "SYMMETRY_CUBIC" ;
    }
   
    // parsed from header file: ../utilities/tensor.h
    static IsotropicMaterialParameters getIsotropicMaterialParameters(std::string type, bool * ok = 0)
    {
        if(ok) { *ok = true ; }
        if( type == "YOUNG_POISSON") { return YOUNG_POISSON ; }
        if( type == "BULK_SHEAR") { return BULK_SHEAR ; }
        if( type == "YOUNG_SHEAR") { return YOUNG_SHEAR ; }
        if(ok) { *ok = false ; }
        return YOUNG_POISSON ;
    }
    static std::string fromIsotropicMaterialParameters(IsotropicMaterialParameters value)
    {
        switch(value)
        {
            case YOUNG_POISSON: return "YOUNG_POISSON" ;
            case BULK_SHEAR: return "BULK_SHEAR" ;
            case YOUNG_SHEAR: return "YOUNG_SHEAR" ;
        }
        return "YOUNG_POISSON" ;
    }
   
    // parsed from header file: ../utilities/tensor.h
    static CompositionType getCompositionType(std::string type, bool * ok = 0)
    {
        if(ok) { *ok = true ; }
        if( type == "SINGLE_OFF_DIAGONAL_VALUES") { return SINGLE_OFF_DIAGONAL_VALUES ; }
        if( type == "DOUBLE_OFF_DIAGONAL_VALUES") { return DOUBLE_OFF_DIAGONAL_VALUES ; }
        if(ok) { *ok = false ; }
        return SINGLE_OFF_DIAGONAL_VALUES ;
    }
    static std::string fromCompositionType(CompositionType value)
    {
        switch(value)
        {
            case SINGLE_OFF_DIAGONAL_VALUES: return "SINGLE_OFF_DIAGONAL_VALUES" ;
            case DOUBLE_OFF_DIAGONAL_VALUES: return "DOUBLE_OFF_DIAGONAL_VALUES" ;
        }
        return "SINGLE_OFF_DIAGONAL_VALUES" ;
    }
   

} ;

}

#endif // __ENUMERATION_TRANSLATOR_H__
