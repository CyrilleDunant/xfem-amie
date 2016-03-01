//
// C++ Interface: post-processor
//
// Description:
//
//
// Author: Alain Giorla <alain.giorla@gmail.com>, (C) 2008-2016
//
// Copyright: See COPYING file that comes with this distribution
//
//
//

#ifndef __POSTPROCESSOR_H__
#define __POSTPROCESSOR_H__

#include "../features/features.h"
#include "../mesher/delaunay.h"

namespace Amie
{
    struct PostProcessor
    {
        int cacheIndex ;
        double instant ;
        bool all ;

        PostProcessor(int c = -1, double t = 0) : cacheIndex(c), instant(t), all(c<0) { if(c<0) { cacheIndex = -1 ; } } 

        virtual Vector postProcess( FeatureTree * F ) = 0 ;

        static void write( std::string file, FeatureTree * F, std::vector<PostProcessor *> posts, bool console = true, bool erase = false ) ;
    } ;

/*PARSE DoNothing PostProcessor */
    struct DoNothingPostProcessor : public PostProcessor
    {
        DoNothingPostProcessor() : PostProcessor(-1,0) { } ;

        virtual Vector postProcess( FeatureTree * F ) { return Vector(0) ; }
    } ;

    struct FieldPostProcessor : public PostProcessor
    {
        FieldType field ;
        std::string variable ;
        std::string correction = "correction_factor" ;

        FieldPostProcessor( std::string var, int c = -1, double t = 0) ; 
        FieldPostProcessor( FieldType f, int c = -1, double t = 0) : PostProcessor(c,t), field(f), variable(std::string()) { } ;

        virtual Vector postProcess( FeatureTree * F ) = 0 ;

        void setCorrectionFactor( std::string cor ) { correction = cor ; }
    } ;

/*PARSE AverageField PostProcessor
        @string[field] // name of the field to consider
        @value[index] -1 // index of the InclusionFamily on which the post-processing is done
        @value[instant] 1 // instant at which the values are captured
*/
    struct AverageFieldPostProcessor : public FieldPostProcessor
    {
        AverageFieldPostProcessor( std::string var, int c = -1, double t = 0) : FieldPostProcessor( var, c, t ) { } ;
        AverageFieldPostProcessor( FieldType f, int c = -1, double t = 0) : FieldPostProcessor(f,c,t) { } ;

        virtual Vector postProcess( FeatureTree * F ) ;
    } ;

/*PARSE MinimumField PostProcessor
        @string[field] // name of the field to consider
        @value[index] -1 // index of the InclusionFamily on which the post-processing is done
        @value[instant] 1 // instant at which the values are captured
*/
    struct MinimumFieldPostProcessor : public FieldPostProcessor
    {
        MinimumFieldPostProcessor( std::string var, int c = -1, double t = 0) : FieldPostProcessor( var, c, t ) { } ;
        MinimumFieldPostProcessor( FieldType f, int c = -1, double t = 0) : FieldPostProcessor(f,c,t) { } ;

        virtual Vector postProcess( FeatureTree * F ) ;
    } ;
    
/*PARSE MaximumField PostProcessor
        @string[field] // name of the field to consider
        @value[index] -1 // index of the InclusionFamily on which the post-processing is done
        @value[instant] 1 // instant at which the values are captured
*/
    struct MaximumFieldPostProcessor : public FieldPostProcessor
    {
        MaximumFieldPostProcessor( std::string var, int c = -1, double t = 0) : FieldPostProcessor( var, c, t ) { } ;
        MaximumFieldPostProcessor( FieldType f, int c = -1, double t = 0) : FieldPostProcessor(f,c,t) { } ;

        virtual Vector postProcess( FeatureTree * F ) ;
    } ;

/*PARSE LocalField PostProcessor
        @string[field] // name of the field to consider
        @value[x] // X-coordinate of the point where the field is measured
        @value[y] // Y-coordinate of the point where the field is measured
        @value[instant] 1 // instant at which the values are captured
*/
    struct LocalFieldPostProcessor : public FieldPostProcessor
    {
        Point p ; 
        DelaunayTriangle * trg ;

        LocalFieldPostProcessor( std::string var, double x, double y, double t = 0) : FieldPostProcessor( var, -1, t ), p(x,y), trg(nullptr) { } ;
        LocalFieldPostProcessor( FieldType f, double x, double y, double t = 0) : FieldPostProcessor(f,-1,t), p(x,y), trg(nullptr) { } ;

        virtual Vector postProcess( FeatureTree * F ) ;
    } ;

/*PARSE ProfileField PostProcessor
        @string[field] // name of the field to consider
        @value[x1] // X-coordinate of the first point in the line
        @value[y1] // Y-coordinate of the first point in the line
        @value[x2] // X-coordinate of the last point in the line
        @value[y2] // Y-coordinate of the last point in the line
        @value[divisions] 10 // number of points where the value is captured
        @value[instant] 1 // instant at which the values are captured
*/
    struct ProfileFieldPostProcessor : public FieldPostProcessor
    {
        Point start ;
        Point end ;
        size_t div ;
        std::vector<PostProcessor *> gauges ;

        ProfileFieldPostProcessor( std::string var, double x1, double y1, double x2, double y2, size_t d = 10, double t = 0) : FieldPostProcessor( var, -1, t ), start(x1,y1), end(x2,y2), div(d) { } ;
        ProfileFieldPostProcessor( FieldType f, double x1, double y1, double x2, double y2, size_t d = 10, double t = 0) : FieldPostProcessor( f, -1, t ), start(x1,y1), end(x2,y2), div(d) { } ;
        
        virtual Vector postProcess( FeatureTree * F ) ;

        PostProcessor * getPostProcessor(size_t i) ;
    } ;
 
/*PARSE LinearStrainGauge PostProcessor
        @value[x1] // X-coordinate of the first anchor point 
        @value[y1] // Y-coordinate of the first anchor point 
        @value[x2] // X-coordinate of the second anchor point 
        @value[y2] // Y-coordinate of the second anchor point 
        @value[instant] 1 // instant at which the values are captured
*/
    struct LinearStrainGaugePostProcessor : public PostProcessor
    {
        LocalFieldPostProcessor left ;
        LocalFieldPostProcessor right ;

        LinearStrainGaugePostProcessor( double x1, double y1, double x2, double y2, double t = 0 ) : PostProcessor(-1,t), left( DISPLACEMENT_FIELD, x1, y1, t), right( DISPLACEMENT_FIELD, x2, y2, t) { }

        virtual Vector postProcess( FeatureTree * F ) ;
    } ;


/*PARSE MacroscopicStrain PostProcessor
        @value[x] // X-coordinate of the vertical line along which the Y displacement is measured
        @value[y] // Y-coordinate of the horizontal line along which the X displacement is measured 
        @value[instant] 1 // instant at which the values are captured
*/
    struct MacroscopicStrainPostProcessor : public PostProcessor
    {
        Point p ;
        Point dim ;
        int left ;
        int right ;
        int bottom ;
        int top ;
        size_t ndof ;

        MacroscopicStrainPostProcessor(double x, double y, double t = 0 ) : PostProcessor(-1, t), p(x,y), dim(0,0), left(-1), right(-1), bottom(-1), top(-1), ndof(0) { } ;

        virtual Vector postProcess( FeatureTree * F ) ;

        static std::pair<Point*, double> getClosestBoundingPoint( DelaunayTriangle * tri, Point p, double instant, Variable axis = ONE ) ;
        static int getTargetID( FeatureTree * F, Point p, double instant, Variable axis = ONE ) ;
    } ;

/*PARSE MaximumMacroscopicStrain PostProcessor
        @value[divisions] 10 // Number of divisions along which the displacements are measured
        @value[instant] 1 // instant at which the values are captured
*/
    struct MaximumMacroscopicStrainPostProcessor : public PostProcessor
    {
        size_t ndiv ;
        std::vector<MacroscopicStrainPostProcessor *> gauges ;

        MaximumMacroscopicStrainPostProcessor( size_t d = 10, double t = 0) : PostProcessor(-1, t), ndiv(d) {  } ;

        virtual Vector postProcess( FeatureTree * F ) ;
    } ;



}

#endif


