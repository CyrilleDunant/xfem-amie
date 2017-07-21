//
// C++ Interface: voxel writer
//
// Description:
//
//
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2010-2011
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "triangle_writer.h"
#include "voxel_writer.h"
#include "../itoa.h"
#include "../../physics/dual_behaviour.h"
#include "../../physics/stiffness.h"
#include "../../physics/materials/paste_behaviour.h"
#include "../../physics/logarithmic_creep_with_external_parameters.h"
#include "../../elements/generalized_spacetime_viscoelastic_element_state.h"
#include <iostream>
#include <fstream>
#include <iomanip>


namespace Amie
{

TriangleWriter::TriangleWriter( std::string f, FeatureTree *F, int t )
{
    filename = f ;
    source = F ;
    counter = 0 ;
    values.push_back(std::vector< std::vector<std::valarray<double> > >(0));

    if( source )
    {
        layers = source->listLayers() ;

        for( size_t j = 0 ; j < layers.size() ; j++ )
        {

            nTriangles.push_back( source->get2DMesh(layers[j])->begin().size() );
            int count = 0 ;

            for( auto i = source->get2DMesh(layers[j])->begin() ; i != source->get2DMesh(layers[j])->end() ; i++ )
                if( i->getBehaviour() && i->getBehaviour()->type != VOID_BEHAVIOUR )
                    count++ ;

            nTriangles.back() = count ;
            timePlane.push_back( t ) ;

            if( timePlane.back() < 0 )
                timePlane.back() = 0 ;

            if( timePlane.back() >= (int)source->get2DMesh(layers[j])->begin()->timePlanes() )
                timePlane.back() = (int)source->get2DMesh(layers[j])->begin()->timePlanes() - 1 ;

            layerTranslator[layers[j]] = j ;
            values.back().push_back( std::vector<std::valarray<double> >( 0 ) );

        }

        getField( TWFT_COORDINATE, false ) ;
        getField( TWFT_DISPLACEMENTS, false ) ;
    }
}

BinaryTriangleWriter::BinaryTriangleWriter( std::string f, FeatureTree *F, int t ) : TriangleWriter( f, F, t ) { }

MultiTriangleWriter::MultiTriangleWriter(std::string h, std::string b, FeatureTree *F, int t) : TriangleWriter(b, F, t)
{
    counter = 0 ;
    head = h ;
    base = b ;

    std::fstream headstream ;

    headstream.open( head.c_str(), std::ios::out ) ;
    headstream << "MULTI" << std::endl ;
    headstream.close();

}

void TriangleWriter::reset( FeatureTree *F, int t )
{
    fields.clear();
    if(counter > (int)values.size()-1)
    {
        values.push_back(std::vector< std::vector<std::valarray<double> > >(0));
        values.back().clear() ;
    }
    nTriangles.clear();
    layers.clear();
    timePlane.clear() ;
    extraFields.clear() ;
    layerTranslator.clear() ;
    source = F ;
    if( source  )
    {
        layers = source->listLayers() ;
        for( size_t j = 0 ; j < layers.size() ; j++ )
        {

            if(source->get2DMesh(layers[j])->begin().size() == 0)
            {
                continue ;
            }
            nTriangles.push_back( source->get2DMesh(layers[j])->begin().size() );
            int count = 0 ;

            for( auto i = source->get2DMesh(layers[j])->begin() ; i != source->get2DMesh(layers[j])->end() ; i++ )
            {

                if( i->getBehaviour() && i->getBehaviour()->type != VOID_BEHAVIOUR )
                {
                    count++ ;
                }
            }
            nTriangles.back() = count ;
            timePlane.push_back( t ) ;

            if( timePlane.back() < 0 )
                timePlane.back() = 0 ;

            if( timePlane.back() >= (int)source->get2DMesh(layers[j])->begin()->timePlanes() )
                timePlane.back() = (int)source->get2DMesh(layers[j])->begin()->timePlanes() - 1 ;

            layerTranslator[layers[j]] = j ;
            values.back().push_back( std::vector<std::valarray<double> >( 0 ) );

        }
        getField( TWFT_COORDINATE, false ) ;
        getField( TWFT_DISPLACEMENTS, false ) ;
    }
}

void HSVtoRGB( double *r, double *g, double *b, double h, double s, double v )
{
    int i;
    double f, p, q, t;

    if ( s == 0 )
    {
        // achromatic (grey)
        *r = *g = *b = v;
        return;
    }

    h /= 60.;                        // sector 0 to 5

    i = ( int )floor( h );
    f = h - i;                      // factorial part of h
    p = v * ( 1. - s );
    q = v * ( 1. - s * f );
    t = v * ( 1. - s * ( 1. - f ) );

    switch ( i )
    {

    case 0:
        *r = v;
        *g = t;
        *b = p;
        break;

    case 1:
        *r = q;
        *g = v;
        *b = p;
        break;

    case 2:
        *r = p;
        *g = v;
        *b = t;
        break;

    case 3:
        *r = p;
        *g = q;
        *b = v;
        break;

    case 4:
        *r = t;
        *g = p;
        *b = v;
        break;

    default:                // case 5:
        *r = v;
        *g = p;
        *b = q;
        break;
    }
}

void TriangleWriter::writeSvg(double factor, bool incolor)
{
    if(nTriangles.empty())
    {
        std::cout << "not writing empty svg file" << std::endl ;
        return ;
    }
    double maxx = 0 ;
    double maxy = 0 ;
    double minx = 0 ;
    double miny = 0 ;

    bool isScalar = false ;
    for( size_t k = 0 ; k < fields.size() ; k++ )
    {
        if(fields[k] == TWFT_SCALAR)
        {
            isScalar = true ;
            break ;
        }
    }

    for(int j = 0 ; j < counter ; j++)
    {
        for( int i = 0 ; i < nTriangles[0] ; i++ )
        {
            maxx = std::max(std::max(std::max(maxx,values[j][0][0][i]+values[j][0][6][i]*factor) ,values[j][0][2][i]+values[j][0][7][i]*factor),values[j][0][4][i]+values[j][0][8][i]*factor) ;
            maxy = std::max(std::max(std::max(maxy,-values[j][0][1][i]-values[j][0][9][i]*factor) ,-values[j][0][3][i]-values[j][0][10][i]*factor),-values[j][0][5][i]-values[j][0][11][i]*factor) ;
            minx = std::min(std::min(std::min(minx,values[j][0][0][i]+values[j][0][6][i]*factor) ,values[j][0][2][i]+values[j][0][7][i]*factor),values[j][0][4][i]+values[j][0][8][i]*factor) ;
            miny = std::min(std::min(std::min(miny,-values[j][0][1][i]-values[j][0][9][i]*factor) ,-values[j][0][3][i]-values[j][0][10][i]*factor),-values[j][0][5][i]-values[j][0][11][i]*factor) ;
        }
    }
    // assume that we want the total height to be 400 px
    double gfactor = 400./(maxy-miny) ;

    for(int m = 0 ; m < counter ; m++)
    {


        std::string fname = base + std::string("_") + itoa(m+1)+".svg" ;
        std::fstream outfile(fname.c_str(), std::ios::out | std::ios_base::trunc) ;
        if(!outfile.is_open())
        {
            std::cout << "failed opening file for SVG output." << std::endl ;
            return ;
        }
        outfile << "<svg xmlns=\"http://www.w3.org/2000/svg\" width =\""<< (maxx-minx)*gfactor*layers.size()*1.1+330<< "\" height =\"" << (maxy-miny)*gfactor*1.1*(values.back()[0].size()-isScalar*5-5)/3.<<"\" version=\"1.1\">" << std::endl ;
        outfile.flush();

        if(incolor)
        {
            //define a linear hue gradient_flux

            outfile <<"\t<linearGradient id=\"hue\" x1=\"0%\" y1=\"0%\" x2=\"0%\" y2=\"100%\">"<< std::endl ;
            outfile <<"\t\t<stop offset=\"0%\" stop-color=\"#ff0000\"/>"<< std::endl ;
            outfile <<"\t\t<stop offset=\"20.7%\" stop-color=\"#ffff00\"/>"<< std::endl ;
            outfile <<"\t\t<stop offset=\"41.46%\" stop-color=\"#00ff00\"/>"<< std::endl ;
            outfile <<"\t\t<stop offset=\"61%\" stop-color=\"#00ffff\"/>"<< std::endl ;
            outfile <<"\t\t<stop offset=\"80.49%\" stop-color=\"#0000ff\"/>"<< std::endl ;
            outfile <<"\t\t<stop offset=\"100%\" stop-color=\"#ff00ff\"/>"<< std::endl ;
            outfile <<"\t\t<stop offset=\"122%\" stop-color=\"#ff0000\"/>"<< std::endl ;
            outfile <<"\t</linearGradient>"<< std::endl ;
            outfile.flush();
        }
        else
        {

            outfile <<"\t<linearGradient id=\"hue\" x1=\"0%\" y1=\"0%\" x2=\"0%\" y2=\"100%\">"<< std::endl ;
            outfile <<"\t\t<stop offset=\"0%\" stop-color=\"#000000\"/>"<< std::endl ;
            outfile <<"\t\t<stop offset=\"100%\" stop-color=\"#ffffff\"/>"<< std::endl ;
            outfile <<"\t</linearGradient>"<< std::endl ;
            outfile.flush();
        }
        std::vector<double> minval ;
        std::vector<double> maxval ;

        for(size_t j = 6 ; j < values.back()[0].size() ; j+=3)
        {
            minval.push_back(values.back()[0][j][0]) ;
            maxval.push_back(values.back()[0][j][0]) ;

            for( size_t l = 0 ; l < layers.size() ; l++ )
            {
                for( int i = 0 ; i < nTriangles[0] ; i++ )
                {

                    for(int p = 0 ; p < counter ; p++)
                    {
                        maxval.back() = std::max(values[p][l][j][i], maxval.back()) ;
                        minval.back() = std::min(values[p][l][j][i], minval.back()) ;
                    }
                }
            }
        }


        std::vector<size_t> skips ;
        for( size_t k = 0 ; k < layers.size() ; k++ )
        {
            outfile << "\t<g>" << std::endl ;
            int jit = 0 ;
            int fieldCounter = 1;
            int localFieldCounter = 3 ;

            for(size_t j = 6 ; j < values.back()[0].size() ; j+=3)
            {
                std::string currentField = nameOfField(fields[fieldCounter], fieldCounter, fieldsOther, fieldsInternal) ;
                if(isScalar && (fields[fieldCounter] == TWFT_COORDINATE || fields[fieldCounter] == TWFT_DISPLACEMENTS))
                {
                    if(localFieldCounter < numberOfFields(fields[fieldCounter],fieldCounter , fieldsOther, fieldsInternal))
                    {
                        localFieldCounter +=3 ;
                    }
                    else
                    {
                        localFieldCounter = 3 ;
                        fieldCounter++ ;
                    }
                    skips.push_back(j);
                    jit++ ;
                    continue ;
                }

                if(localFieldCounter < numberOfFields(fields[fieldCounter],fieldCounter, fieldsOther, fieldsInternal))
                {
                    localFieldCounter +=3 ;
                }
                else
                {
                    localFieldCounter = 3 ;
                    fieldCounter++ ;
                }


                double dval = (maxval[jit]-minval[jit]) ;
                if (dval < 1e-14)
                    dval = 1 ;
                outfile << "\t\t<!-- maxval = " <<  maxval[jit] << ", minval = " << minval[jit] << "-->" << std::endl ;
                outfile << "\t\t<g>" << std::endl ;
                outfile << "\t\t<text font-size=\"28\" x=\"" << -gfactor*minx+((maxx-minx)*gfactor*1.1)*(k-0.3+0.05) << "\" y=\"" << (maxy-miny)*gfactor*1.1*(j-isScalar*5-6+0.27)/3. << "\">" << currentField << "</text>" << std::endl ;
                for( int i = 0 ; i < nTriangles[0] ; i++ )
                {
                    double val = (values[m][k][j][i]+values[m][k][j+1][i]+values[m][k][j+2][i])/3. ;
                    double r,g,b ;
                    if(incolor)
                        HSVtoRGB(&r,&g,&b,295.2*(maxval[jit]-val)/dval+2, 1., 1.);
                    else
                        HSVtoRGB(&r,&g,&b,180., 0., (maxval[jit]-val)/dval);

                    outfile << "\t\t\t<polygon points =\""
                            <<  gfactor*(-minx+values[m][k][0][i]+values[m][k][6][i]*factor)+((maxx-minx)*gfactor*1.1)*(k+0.05) << "," << (maxy-miny)*gfactor*1.1*(j-isScalar*5-6+0.3)/3.+gfactor*(maxy-values[m][k][1][i]-values[m][k][9][i]*factor) << " "
                            <<  gfactor*(-minx+values[m][k][2][i]+values[m][k][7][i]*factor)+((maxx-minx)*gfactor*1.1)*(k+0.05) << "," << (maxy-miny)*gfactor*1.1*(j-isScalar*5-6+0.3)/3.+gfactor*(maxy-values[m][k][3][i]-values[m][k][10][i]*factor) << " "
                            <<  gfactor*(-minx+values[m][k][4][i]+values[m][k][8][i]*factor)+((maxx-minx)*gfactor*1.1)*(k+0.05) << ","<< (maxy-miny)*gfactor*1.1*(j-isScalar*5-6+0.3)/3.+gfactor*(maxy-values[m][k][5][i]-values[m][k][11][i]*factor)<< " "
                            << "\" style=\"stroke-width:0.0;stroke:darkgray;fill-opacity:1\" fill=\"rgb("<< r*100 <<"%,"<<g*100<<"%," <<b*100 << "%)"<<"\"" << std::flush ;

                    outfile << " />" << std::endl ;
                }
                outfile << "\t\t</g>" << std::endl ;
                outfile.flush();
                jit++ ;
            }
            outfile << "\t</g>" << std::endl ;
        }
        int jit = 0 ;

        for(size_t j = 6 ; j < values.back()[0].size() ; j+=3)
        {
            for(size_t k = 0 ; k < skips.size() ; k++)
                if(j == skips[k])
                    goto postScale ;

            outfile << "\t\t<g>" << std::endl ;
            outfile << "\t\t\t<text font-size=\"28\" x=\"" << -gfactor*minx+((maxx-minx)*gfactor*1.1)*(layers.size()-0.5+0.05)+60
                    << "\" y=\"" << (maxy-miny)*gfactor*1.1*((j-isScalar*5-6+0.5)/3.)+12 << "\">" << std::setprecision(4) << maxval[jit] << "</text>" << std::endl ;
            outfile << "\t\t\t<text font-size=\"28\" x=\"" << -gfactor*minx+((maxx-minx)*gfactor*1.1)*(layers.size()-0.5+0.05)+60
                    << "\" y=\"" <<  (maxy-miny)*gfactor*1.1*((j-isScalar*5-6+0.5)/3.)+0.5*(maxy-miny)*gfactor*.35+12 << "\">" << std::setprecision(4) <<  minval[jit]*0.25+maxval[jit]*0.75 << "</text>" << std::endl ;
            outfile << "\t\t\t<text font-size=\"28\" x=\"" << -gfactor*minx+((maxx-minx)*gfactor*1.1)*(layers.size()-0.5+0.05)+60
                    << "\" y=\"" << (maxy-miny)*gfactor*1.1*((j-isScalar*5-6+0.5)/3.)+(maxy-miny)*gfactor*.35+12  << "\">" << std::setprecision(4) << (minval[jit]+maxval[jit])*.5 << "</text>" << std::endl ;
            outfile << "\t\t\t<text font-size=\"28\" x=\"" << -gfactor*minx+((maxx-minx)*gfactor*1.1)*(layers.size()-0.5+0.05)+60
                    << "\" y=\"" << (maxy-miny)*gfactor*1.1*((j-isScalar*5-6+0.5)/3.)+1.5*(maxy-miny)*gfactor*.35+12  << "\">" << std::setprecision(4) << minval[jit]*0.75+maxval[jit]*0.25 << "</text>" << std::endl ;
            outfile << "\t\t\t<text font-size=\"28\" x=\"" << -gfactor*minx+((maxx-minx)*gfactor*1.1)*(layers.size()-0.5+0.05)+60
                    << "\" y=\"" << (maxy-miny)*gfactor*1.1*((j-isScalar*5-6+0.5)/3.)+2.*(maxy-miny)*gfactor*.35+12   << "\">"<< std::setprecision(4) << minval[jit] << "</text>" << std::endl ;


            outfile << "\t\t\t<rect width=\"30\" height=\"" <<  (maxy-miny)*gfactor*.7
                    << "\" x=\"" << -gfactor*minx+((maxx-minx)*gfactor*1.1)*(layers.size()-0.5+0.05) + 20
                    << "\" y=\"" << (maxy-miny)*gfactor*1.1*((j-isScalar*5-6+0.5)/3.)
                    << "\"  style=\"fill-opacity:1\" fill=\"url(#hue)\"/>" << std::endl ;
            outfile << "\t\t</g>" << std::endl ;

postScale:
            jit++ ;
            outfile.flush();
        }

        outfile << "</svg>"<< std::endl ;
        outfile.flush();
        outfile.close();

    }
// 	wait(1) ;
// 	for(size_t m = 0 ; m < counter ; m++)
// 	  outfile[m].close();

}

void TriangleWriter::write()
{

    writeHeader( layers[0], false ) ;
    std::fstream outfile  ;
    outfile.open( filename.c_str(), std::ios::out | std::ios::app ) ;

    for( int i = 0 ; i < nTriangles[0] ; i++ )
    {
        for( size_t j = 0 ; j < values.back()[0].size() ; j++ )
        {
            outfile << values.back()[0][j][i] << " " ;
        }

        outfile << std::endl ;
    }

    outfile.close();

    std::string filename_orig = filename;

    for( size_t k = 1 ; k < layers.size() ; k++ )
    {
        writeHeader( layers[k], true ) ;
        outfile.open( filename_orig.c_str(), std::ios::out | std::ios::app ) ;

        for( int i = 0 ; i < nTriangles[k] ; i++ )
        {
            for( size_t j = 0 ; j < values.back()[k].size() ; j++ )
            {
                outfile << values.back()[k][j][i] << " " ;
            }

            outfile << std::endl ;
        }

        outfile.close();
    }
}

void BinaryTriangleWriter::write()
{
    writeHeader( false ) ;

    std::vector<std::valarray<unsigned char> > norm ;
    std::valarray<bool> voids( false, values.back()[0][0].size() ) ;

    for( size_t i = 6 ; i < values.back()[0].size() ; i++ )
    {
        std::valarray<unsigned char> n = normalizeArray( values.back()[0][i], voids ) ;
        norm.push_back( n ) ;
    }

    std::fstream outbin ;
    outbin.open( filename.c_str(), std::ios::out | std::ios::app ) ;

    for( int i = 0 ; i < nTriangles[0] ; i++ )
    {
        for( size_t j = 0 ; j < 6 ; j++ )
            outbin << values.back()[0][j][i] << " " ;

        for( size_t j = 6 ; j < values.back()[0].size() ; j++ )
        {
            unsigned char val = norm[j - 6][i] ;
            outbin << ( unsigned short int ) val << " " ;
        }

        outbin << std::endl ;
    }

    outbin.close();
}


void TriangleWriter::append()
{

    for( size_t k = 0 ; k < nTriangles.size() ; k++ )
    {
        writeHeader( true ) ;
        std::fstream outfile  ;
        outfile.open( filename.c_str(), std::ios::out | std::ios::app ) ;

        for( int i = 0 ; i < nTriangles[k] ; i++ )
        {
            for( size_t j = 0 ; j < values.back()[k].size() ; j++ )
            {
                outfile << values.back()[k][j][i] << " " ;
            }

            outfile << std::endl ;
        }

        outfile.close();
    }

}

void BinaryTriangleWriter::append()
{
    writeHeader( true ) ;

    std::vector<std::valarray<unsigned char> > norm ;
    std::valarray<bool> voids( false, values[0].size() ) ;

    for( size_t i = 6 ; i < values.back().size() ; i++ )
    {
        std::valarray<unsigned char> n = normalizeArray( values.back()[0][i], voids ) ;
        norm.push_back( n ) ;
    }

    std::fstream outbin ;
    outbin.open( filename.c_str(), std::ios::out | std::ios::app ) ;

    for( int i = 0 ; i < nTriangles[0] ; i++ )
    {
        for( size_t j = 0 ; j < 6 ; j++ )
            outbin << values.back()[0][j][i] << " " ;

        for( size_t j = 6 ; j < values.back()[0].size() ; j++ )
        {
            unsigned char val = norm[j - 6][i] ;
            outbin << ( unsigned short int ) val << " " ;
        }

        outbin << std::endl ;
    }

    outbin.close();
}

void MultiTriangleWriter::append()
{
    writeHeader( layers[0], true ,false ) ;
    std::fstream outfile  ;
    outfile.open( filename.c_str(), std::ios::out | std::ios::app ) ;

    for( int i = 0 ; i < nTriangles[0] ; i++ )
    {
        for( size_t j = 0 ; j < values.back()[0].size() ; j++ )
        {
            outfile << values.back()[0][j][i] << " " ;
        }

        outfile << std::endl ;
    }
    outfile.close();

    std::string filename_orig = filename;

    for( size_t k = 1 ; k < layers.size() ; k++ )
    {
        writeHeader( layers[k], false, true ) ;
        outfile.open( filename_orig.c_str(), std::ios::out | std::ios::app ) ;

        for( int i = 0 ; i < nTriangles[k] ; i++ )
        {
            for( size_t j = 0 ; j < values.back()[k].size() ; j++ )
            {
                outfile << values.back()[k][j][i] << " " ;
            }

            outfile << std::endl ;
        }

        outfile.close();
    }
    writeIndexFile();
}


void TriangleWriter::getField( TWFieldType field, bool extra, std::string fieldName, double offset )
{

    if(field == TWFT_INTERNAL_VARIABLE)
    {
        fieldsInternal[fields.size()] = fieldName ;
    }
    fields.push_back(field);
    for( size_t j = 0 ; j < layers.size() ; j++ )
    {
        std::vector<std::valarray<double> > val = getDoubleValues( field, fields.size()-1, layers[j], fieldName, offset ) ;
        std::reverse( val.begin(), val.end() );
        values.back()[layerTranslator[layers[j]]].insert( values.back()[layerTranslator[layers[j]]].end(), val.begin(), val.end() ) ;
    }

}


void TriangleWriter::getField( FieldType field, bool extra )
{
    fieldsOther[fields.size()] = field;
    fields.push_back(TWFT_FIELD_TYPE);
    for( size_t j = 0 ; j < layers.size() ; j++ )
    {
        std::vector<std::valarray<double> > val = getDoubleValues( field, layers[j] ) ;
        std::reverse( val.begin(), val.end() );
        values.back()[layerTranslator[layers[j]]].insert( values.back()[layerTranslator[layers[j]]].end(), val.begin(), val.end() ) ;
    }

}

std::vector<std::valarray<double> > TriangleWriter::getDoubleValues( TWFieldType field, size_t index, int layer, std::string fieldName, double offset )
{
    std::vector<std::valarray<double> > ret ;
    int iterator = 0 ;

    for( int i = 0 ; i < numberOfFields( field, index , fieldsOther, fieldsInternal) ; i++ )
    {
        std::valarray<double> reti( nTriangles[layerTranslator[layer]] ) ;
        ret.push_back( reti ) ;
    }

    if( field == TWFT_GRADIENT || field == TWFT_GRADIENT_AND_FLUX || field == TWFT_FLUX )
    {
        std::pair<Vector, Vector> gradient_flux = source->getGradientAndFluxInLayer( layer, false ) ;
// 		std::vector<DelaunayTriangle *> triangles = source->getElements2DInLayer( layer ) ;

        switch( field )
        {


        case TWFT_FLUX:
            gradient_flux.first.resize( 0 ) ;

            for( auto i = source->get2DMesh(layer)->begin() ; i != source->get2DMesh(layer)->end() ; i++ )
            {
                if(  i->getBehaviour() && i->getBehaviour()->type != VOID_BEHAVIOUR && !i->getBehaviour()->fractured() )
                {
                    // j11
                    ret[5][iterator] = gradient_flux.second[i.getPosition() * 3 * 2 + 0] ;
                    ret[4][iterator] = gradient_flux.second[i.getPosition() * 3 * 2 + 2] ;
                    ret[3][iterator] = gradient_flux.second[i.getPosition() * 3 * 2 + 4] ;

                    // j22
                    ret[2][iterator] = gradient_flux.second[i.getPosition() * 3 * 2 + 1] ;
                    ret[1][iterator] = gradient_flux.second[i.getPosition() * 3 * 2 + 3] ;
                    ret[0][iterator++] = gradient_flux.second[i.getPosition() * 3 * 2 + 5] ;
                }
                else if( i->getBehaviour() && i->getBehaviour()->type != VOID_BEHAVIOUR)
                {
                    // j11
                    ret[5][iterator] = 0 ;
                    ret[4][iterator] = 0 ;
                    ret[3][iterator] = 0 ;

                    // j22
                    ret[2][iterator] = 0 ;
                    ret[1][iterator] = 0 ;
                    ret[0][iterator++] = 0 ;
                }
            }

            break ;


        case TWFT_GRADIENT_AND_FLUX:

            for( auto i = source->get2DMesh(layer)->begin() ; i != source->get2DMesh(layer)->end() ; i++ )
            {
                if(  i->getBehaviour() && i->getBehaviour()->type != VOID_BEHAVIOUR && !i->getBehaviour()->fractured() )
                {
                    // d11
                    ret[11][iterator] = gradient_flux.first[i.getPosition() * 3 * 2 + 0] ;
                    ret[10][iterator] = gradient_flux.first[i.getPosition() * 3 * 2 + 2] ;
                    ret[9][iterator] = gradient_flux.first[i.getPosition() * 3 * 2 + 4] ;

                    // d22
                    ret[8][iterator] = gradient_flux.first[i.getPosition() * 3 * 2 + 1] ;
                    ret[7][iterator] = gradient_flux.first[i.getPosition() * 3 * 2 + 3] ;
                    ret[6][iterator] = gradient_flux.first[i.getPosition() * 3 * 2 + 5] ;

                    // j11
                    ret[5][iterator] = gradient_flux.second[i.getPosition() * 3 * 2 + 0] ;
                    ret[4][iterator] = gradient_flux.second[i.getPosition() * 3 * 2 + 2] ;
                    ret[3][iterator] = gradient_flux.second[i.getPosition() * 3 * 2 + 4] ;

                    // j22
                    ret[2][iterator] = gradient_flux.second[i.getPosition() * 3 * 2 + 1] ;
                    ret[1][iterator] = gradient_flux.second[i.getPosition() * 3 * 2 + 3] ;
                    ret[0][iterator++] = gradient_flux.second[i.getPosition() * 3 * 2 + 5] ;
                }
                else if( i->getBehaviour() && i->getBehaviour()->type != VOID_BEHAVIOUR )
                {
                    ret[11][iterator] = 0 ;
                    ret[10][iterator] = 0 ;
                    ret[9][iterator] = 0 ;

                    // d22
                    ret[8][iterator] = 0 ;
                    ret[7][iterator] = 0 ;
                    ret[6][iterator] = 0 ;

                    // j11
                    ret[5][iterator] = 0 ;
                    ret[4][iterator] = 0 ;
                    ret[3][iterator] = 0 ;

                    // j22
                    ret[2][iterator] = 0 ;
                    ret[1][iterator] = 0 ;
                    ret[0][iterator++] =0 ;
                }
            }

            break ;

        case TWFT_GRADIENT:
            gradient_flux.second.resize( 0 ) ;

            for( auto i = source->get2DMesh(layer)->begin() ; i != source->get2DMesh(layer)->end() ; i++ )
            {
                if(  i->getBehaviour() &&i->getBehaviour()->type != VOID_BEHAVIOUR && !i->getBehaviour()->fractured() )
                {
                    // d11
                    ret[5][iterator] = gradient_flux.first[i.getPosition() * 3 * 2 + 0] ;
                    ret[4][iterator] = gradient_flux.first[i.getPosition() * 3 * 2 + 2] ;
                    ret[3][iterator] = gradient_flux.first[i.getPosition() * 3 * 2 + 4] ;

                    // d22
                    ret[2][iterator] = gradient_flux.first[i.getPosition() * 3 * 2 + 1] ;
                    ret[1][iterator] = gradient_flux.first[i.getPosition() * 3 * 2 + 3] ;
                    ret[0][iterator++] = gradient_flux.first[i.getPosition() * 3 * 2 + 5] ;
                }
                else if(i->getBehaviour() &&i->getBehaviour()->type != VOID_BEHAVIOUR )
                {
                    ret[5][iterator] = 0 ;
                    ret[4][iterator] = 0 ;
                    ret[3][iterator] = 0 ;

                    ret[2][iterator] = 0 ;
                    ret[1][iterator] = 0 ;
                    ret[0][iterator++] = 0;
                }
            }

            break ;

        default:
            break ;

        }

    }
    else
    {
        if( field == TWFT_DISPLACEMENTS )
        {
            Vector x = source->getDisplacements(-1, false) ;
// 			std::vector<DelaunayTriangle *> triangles = source->getElements2DInLayer( layer ) ;
            int pointsPerTri = source->get2DMesh(layer)->begin()->getBoundingPoints().size() ;

            int pointsPerTimePlanes = pointsPerTri / source->get2DMesh(layer)->begin()->timePlanes() ;
            int factor = pointsPerTimePlanes / 3 ;
            if( timePlane[layerTranslator[layer]] >= (signed long) source->get2DMesh(layer)->begin()->timePlanes() )
                timePlane[layerTranslator[layer]] = source->get2DMesh(layer)->begin()->timePlanes() - 1 ;

            int time_offset = timePlane[layerTranslator[layer]] * pointsPerTri / source->get2DMesh(layer)->begin()->timePlanes() ;

           
            for( auto i = source->get2DMesh(layer)->begin() ; i != source->get2DMesh(layer)->end() ; i++ )
            {
                if(  i->getBehaviour() && i->getBehaviour()->type != VOID_BEHAVIOUR )
                { 
                    VirtualMachine vm ;
                    size_t dof = i->getBehaviour()->getNumberOfDegreesOfFreedom() ;

                    size_t id1 = i->getBoundingPoint( factor * 0 + time_offset ).getId() ;
                    size_t id2 = i->getBoundingPoint( factor * 1 + time_offset ).getId() ;
                    size_t id3 = i->getBoundingPoint( factor * 2 + time_offset ).getId() ;

                    
//                     ret[5][iterator] = x[id1 * dof] ;
//                     ret[4][iterator] = x[id2 * dof] ;
//                     ret[3][iterator] = x[id3 * dof] ;
//                     ret[2][iterator] = x[id1 * dof + 1] ;
//                     ret[1][iterator] = x[id2 * dof + 1] ;
//                     ret[0][iterator++] = x[id3 * dof + 1] ;
                    
//                     Matrix tran =  i->getState().transform(&vm);//identity(2); ;// inverse2x2Matrix(i->getState().transform(&vm))
                    Vector xy = { x[id1 * dof], x[id1 * dof + 1]} ;
//                     xy = xy*tran ;
                    ret[5][iterator] = xy[0] ;
                    ret[2][iterator] = xy[1] ;
                    
                    xy = { x[id2 * dof], x[id2 * dof + 1]} ;
//                     xy = xy*tran ;
                    ret[4][iterator] = xy[0] ;
                    ret[1][iterator] = xy[1] ;
                    
                    xy = { x[id3 * dof], x[id3 * dof + 1]} ;
//                     xy = xy*tran ;
                    ret[3][iterator] = xy[0] ;
                    ret[0][iterator++] = xy[1] ;
                }
            }
        }
        else if( field == TWFT_LARGE_DEFORMATION_TRANSFORM)
        {
            VirtualMachine vm ;
            for( auto i = source->get2DMesh(layer)->begin() ; i != source->get2DMesh(layer)->end() ; i++ )
            {
                if( i->getBehaviour() && i->getBehaviour()->type != VOID_BEHAVIOUR )
                {
                    if(i->getState().getGlobalTransform().numCols() == 0)
                    {
                            ret[0][iterator] = 0;
                            ret[1][iterator] = 0 ;
                            ret[2][iterator++] = 0 ;
                            continue ;
                    }
                    
                    Matrix m = i->getState().getGlobalTransform() ;
                    double d = sqrt(std::inner_product(&m.array()[0], &m.array()[m.array().size()], &m.array()[0], double(0))) ;
                    
                    ret[0][iterator] = d ;
                    ret[1][iterator] = d ;
                    ret[2][iterator++] = d ;
                    
                }
                else if ( i->getBehaviour() && i->getBehaviour()->type != VOID_BEHAVIOUR )
                {
                    ret[0][iterator] = 0;
                    ret[1][iterator] = 0 ;
                    ret[2][iterator++] = 0 ;
                }
            }
            
        }
        else if( field == TWFT_LARGE_DEFORMATION_ANGLE)
        {
            VirtualMachine vm ;
            for( auto i = source->get2DMesh(layer)->begin() ; i != source->get2DMesh(layer)->end() ; i++ )
            {
                if( i->getBehaviour() && i->getBehaviour()->type != VOID_BEHAVIOUR )
                {
                    double d = 180.*i->getState().getGlobalTransformAngle()/M_PI ;
                    
                    ret[0][iterator] = d ;
                    ret[1][iterator] = d ;
                    ret[2][iterator++] = d ;
                    
                }
                else if ( i->getBehaviour() && i->getBehaviour()->type != VOID_BEHAVIOUR )
                {
                    ret[0][iterator] = 0;
                    ret[1][iterator] = 0 ;
                    ret[2][iterator++] = 0 ;
                }
            }
            
        }
        else if( field == TWFT_SCALAR )
        {
            Vector x = source->getDisplacements(-1, false) ;

            int pointsPerTri = source->get2DMesh(layer)->begin()->getBoundingPoints().size() ;
            int pointsPerTimePlanes = pointsPerTri / source->get2DMesh(layer)->begin()->timePlanes() ;
            int factor = pointsPerTimePlanes / 3 ;
            if( timePlane[layerTranslator[layer]] >= (int)source->get2DMesh(layer)->begin()->timePlanes() )
                timePlane[layerTranslator[layer]] = (int)source->get2DMesh(layer)->begin()->timePlanes() - 1 ;

            int time_offset = timePlane[layerTranslator[layer]] * pointsPerTri / source->get2DMesh(layer)->begin()->timePlanes() ;

            for( auto i = source->get2DMesh(layer)->begin() ; i != source->get2DMesh(layer)->end() ; i++ )
            {
                if(  i->getBehaviour() && i->getBehaviour()->type != VOID_BEHAVIOUR )
                {
                    size_t id1 = i->getBoundingPoint( factor * 0 + time_offset ).getId() ;
                    size_t id2 = i->getBoundingPoint( factor * 1 + time_offset ).getId() ;
                    size_t id3 = i->getBoundingPoint( factor * 2 + time_offset ).getId() ;

                    ret[2][iterator] = x[id1  ] ;
                    ret[1][iterator] = x[id2  ] ;
                    ret[0][iterator++] = x[id3 ] ;
                }
            }
        }
        else if( field == TWFT_DOH )
        {

            for( auto i = source->get2DMesh(layer)->begin() ; i != source->get2DMesh(layer)->end() ; i++ )
            {
                if( i->getBehaviour() && i->getBehaviour()->type != VOID_BEHAVIOUR && dynamic_cast<HydratingDiffusionCementPaste *>(i->getBehaviour()) )
                {
                    double d = dynamic_cast<HydratingDiffusionCementPaste *>(i->getBehaviour())->doh ;

                    ret[0][iterator] = d;
                    ret[1][iterator] = d ;
                    ret[2][iterator++] = d ;

                }
                else if ( i->getBehaviour() && i->getBehaviour()->type != VOID_BEHAVIOUR )
                {
                    ret[0][iterator] = 0;
                    ret[1][iterator] = 0 ;
                    ret[2][iterator++] = 0 ;
                }
            }
        }
// 		else if( field == TWFT_CRACKS )
// 		{
// 			std::vector<DelaunayTriangle *> triangles = source->getElements2DInLayer( layer ) ;
//
// 			for( int i = 0 ; i < triangles.size() ; i++ )
// 			{
// 				if( triangles[i]->getBehaviour() && triangles[i]->getBehaviour()->type != VOID_BEHAVIOUR &&  triangles[i]->getBehaviour()->getDamageModel() &&  triangles[i]->getBehaviour()->getDamageModel()->getState().max() > POINT_TOLERANCE)
// 				{
//
// 					std::pair<double, double> np = triangles[i]->getBehaviour()->getFractureCriterion()->getCrackOpeningAndSlip(triangles[i]->getState()) ;
// 					ret[0][iterator] = np.first;
// 					ret[1][iterator] = np.first ;
// 					ret[2][iterator] = np.first ;
// 					ret[3][iterator] = np.second;
// 					ret[4][iterator] = np.second ;
// 					ret[5][iterator++] = np.second ;
//
// 				}
// 				else if ( triangles[i]->getBehaviour() && triangles[i]->getBehaviour()->type != VOID_BEHAVIOUR )
// 				{
// 					ret[0][iterator] = 0;
// 					ret[1][iterator] = 0 ;
// 					ret[2][iterator] = 0 ;
// 					ret[3][iterator] = 0;
// 					ret[4][iterator] = 0 ;
// 					ret[5][iterator++] = 0 ;
// 				}
// 			}
// 		}
        else if( field == TWFT_DAMAGE)
        {

            for( auto i = source->get2DMesh(layer)->begin() ; i != source->get2DMesh(layer)->end() ; i++ )
            {
                if( i->getBehaviour() && i->getBehaviour()->type != VOID_BEHAVIOUR &&  i->getBehaviour()->getDamageModel() )
                {
                    double d = 0 ;
                    for(size_t m = 0 ;  m <  i->getBehaviour()->getDamageModel()->getState().size() ; m++)
                        d =  std::max(i->getBehaviour()->getDamageModel()->getState()[m], d) ;
                    if( i->getBehaviour()->getDamageModel()->fractured() )
                        d = 1 ;

                    ret[0][iterator] = d;
                    ret[1][iterator] = d ;
                    ret[2][iterator++] = d ;

                }
                else if ( i->getBehaviour() && i->getBehaviour()->type != VOID_BEHAVIOUR )
                {
                    ret[0][iterator] = 0;
                    ret[1][iterator] = 0 ;
                    ret[2][iterator++] = 0 ;
                }
            }
        }
        else if( field == TWFT_IMPOSED_STRESS_NORM )
        {

            Vector x( source->get2DMesh(layer)->begin().size() * 3 ) ;

            for( auto i = source->get2DMesh(layer)->begin() ; i != source->get2DMesh(layer)->end() ; i++ )
            {
                if( i->getBehaviour() && i->getBehaviour()->type != VOID_BEHAVIOUR &&  i->getBehaviour()->hasInducedForces() && i->getBehaviour()->getDamageModel())
                {
                    Vector dv = i->getBehaviour()->getDamageModel()->getImposedStress(i->getCenter());

                    double d = sqrt(std::inner_product(&dv[0], &dv[dv.size()], &dv[0], 0.));
                    if( i->getBehaviour()->getDamageModel()->fractured() )
                        d = 1 ;

                    ret[0][iterator] = d;
                    ret[1][iterator] = d ;
                    ret[2][iterator++] = d ;

                }
                else if ( i->getBehaviour() && i->getBehaviour()->type != VOID_BEHAVIOUR )
                {
                    ret[0][iterator] = 0;
                    ret[1][iterator] = 0 ;
                    ret[2][iterator++] = 0 ;
                }
            }
        }
        else if( field == TWFT_INTERNAL_VARIABLE)
        {
            std::map<std::string, double> empty ;
            for( auto i = source->get2DMesh(layer)->begin() ; i != source->get2DMesh(layer)->end() ; i++ )
            {
                if( fieldName.find("C") == 0 && fieldName.size() == 3 && i->getBehaviour() && i->getBehaviour()->type != VOID_BEHAVIOUR )
                {
                    int r = ctoi( fieldName.at(1) ) -1 ;
                    int c = ctoi( fieldName.at(2) ) -1 ;
                    double v = std::abs((i->getBehaviour())->getTensor( Point(0,0,0) )[r][c]) ;
                    ret[0][iterator] = v ;
                    ret[1][iterator] = v ;
                    ret[2][iterator++] = v ;
                }
                else if( fieldName.find("E") == 0 && fieldName.size() == 3 && i->getBehaviour() && i->getBehaviour()->type != VOID_BEHAVIOUR )
                {
                    int r = ctoi( fieldName.at(1) ) -1 ;
                    int c = ctoi( fieldName.at(2) ) -1 ;
                    double v = (i->getBehaviour())->getViscousTensor( Point(0,0,0) )[r][c] ;
                    ret[0][iterator] = v ;
                    ret[1][iterator] = v ;
                    ret[2][iterator++] = v ;
                }
                else if( i->getBehaviour() && i->getBehaviour()->type != VOID_BEHAVIOUR && dynamic_cast<LogarithmicCreepWithExternalParameters *>(i->getBehaviour()) )
                {
                    double v = dynamic_cast<GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables&>(i->getState()).get(fieldName, empty ) - offset;
                    ret[0][iterator] = v ;
                    ret[1][iterator] = v ;
                    ret[2][iterator++] = v ;
                }
                else if ( i->getBehaviour() && i->getBehaviour()->type != VOID_BEHAVIOUR )
                {
                    ret[0][iterator] = 0;
                    ret[1][iterator] = 0 ;
                    ret[2][iterator++] = 0 ;
                }
            }
        }
        else
        {

            for( auto i = source->get2DMesh(layer)->begin() ; i != source->get2DMesh(layer)->end() ; i++ )
            {
                std::pair<bool, std::vector<double> > val = getDoubleValue( i, field , index) ;

                if(  i->getBehaviour() && i->getBehaviour()->type != VOID_BEHAVIOUR )
                {
                    if( val.first )
                    {
                        for( int j = 0 ; j < numberOfFields( field , index, fieldsOther, fieldsInternal) ; j++ )
                            ret[j][iterator] = val.second[j] ;

                        iterator++ ;
                    }
                }
            }
        }
    }

    return ret ;
}

std::vector<std::valarray<double> > TriangleWriter::getDoubleValues( FieldType field, int layer )
{
    std::vector<std::valarray<double> > ret ;
    int iterator = 0 ;

    int blocks = source->get2DMesh(layer)->begin()->getBehaviour()->getNumberOfDegreesOfFreedom() / 2 ;
    int size = fieldTypeElementarySize( field, SPACE_TWO_DIMENSIONAL, blocks) ;
    int pointsPerTri = source->get2DMesh(layer)->begin()->getBoundingPoints().size() ;
    int pointsPerPlane = pointsPerTri ;
    if(source->get2DMesh(layer)->begin()->timePlanes() > 1)
        pointsPerPlane /= source->get2DMesh(layer)->begin()->timePlanes() ;
// 	int factor = pointsPerPlane/3 ;

    int time_offset = timePlane[layerTranslator[layer]] * pointsPerTri / source->get2DMesh(layer)->begin()->timePlanes() ;

    for( int i = 0 ; i < size*3 ; i++ )
    {
        std::valarray<double> reti( nTriangles[layerTranslator[layer]] ) ;
        ret.push_back( reti ) ;
    }

    VirtualMachine vm ;

    for( auto i = source->get2DMesh(layer)->begin() ; i != source->get2DMesh(layer)->end() ; i++ )
    {
        if(  i->getBehaviour() && i->getBehaviour()->type != VOID_BEHAVIOUR )
        {
            Vector first(size) ;
            Vector second(size) ;
            Vector third(size) ;

            size_t n = i->getBoundingPoints().size()/3 ;
            if(i->timePlanes() > 1)
                n /= i->timePlanes() ;

            i->getState().getField(field,  i->getBoundingPoint(time_offset+0), first, false, &vm);
            i->getState().getField(field,  i->getBoundingPoint(time_offset+n), second, false, &vm);
            i->getState().getField(field,  i->getBoundingPoint(time_offset+2*n), third, false, &vm);
//             i->getState().getField(field,  i->getCenter(), first, false, 0);
//             i->getState().getField(field,  i->getCenter(), second, false, 0);
//             i->getState().getField(field,  i->getCenter(), third, false, 0);

            for(int j = 0 ; j < size ; j++)
            {
                ret[(size-1-j)*3+2][iterator] = first[j] ;
                ret[(size-1-j)*3+1][iterator] = second[j] ;
                ret[(size-1-j)*3+0][iterator] = third[j] ;
            }

            iterator++ ;
        }
    }

    return ret ;
}

std::pair<bool, std::vector<double> > TriangleWriter::getDoubleValue( DelaunayTriangle *tri, TWFieldType field , size_t index)
{
    bool found = false ;
    std::vector<double> ret( numberOfFields( field,  index, fieldsOther, fieldsInternal) ) ;

    if( tri->getBehaviour() && tri->getBehaviour()->type != VOID_BEHAVIOUR )
    {
        switch( field )
        {
        case TWFT_COORDINATE:
        {
            size_t n = tri->getBoundingPoints().size()/3 ;
            if(tri->timePlanes() > 1)
                n /= tri->timePlanes() ;
            ret[5] = tri->getBoundingPoint(0).getX() ;
            ret[4] = tri->getBoundingPoint(0).getY() ;
            ret[3] = tri->getBoundingPoint(n).getX() ;
            ret[2] = tri->getBoundingPoint(n).getY() ;
            ret[1] = tri->getBoundingPoint(2*n).getX() ;
            ret[0] = tri->getBoundingPoint(2*n).getY() ;
            found = true ;
            break ;
        }
        case TWFT_PRINCIPAL_ANGLE:
        {
            Vector v(0., 1) ;
            tri->getState().getField( PRINCIPAL_STRESS_ANGLE_FIELD, *tri->first, v, false , 0) ;
            ret[2] = 180.*v[0]/M_PI ;
            tri->getState().getField( PRINCIPAL_STRESS_ANGLE_FIELD, *tri->second, v, false, 0 ) ;
            ret[1] = 180.*v[0]/M_PI ;
            tri->getState().getField( PRINCIPAL_STRESS_ANGLE_FIELD, *tri->third, v, false , 0) ;
            ret[0] = 180.*v[0]/M_PI ;

            found = true ;
            break ;
        }
        case TWFT_TRIANGLE_ANGLE:
        {
            Vector v(0., 1) ;
            double a = tri->area() ;
            double r = tri->getRadius() ;
            ret[2] = a/(M_PI*r*r) ;
            ret[1] = a/(M_PI*r*r) ;
            ret[0] = a/(M_PI*r*r) ;
// 				tri->getState().getField( PRINCIPAL_STRESS_ANGLE_FIELD, *tri->first, v, false ) ;
// 				ret[2] = 180.*v[0]/M_PI ;
// 				tri->getState().getField( PRINCIPAL_STRESS_ANGLE_FIELD, *tri->second, v, false ) ;
// 				ret[1] = 180.*v[0]/M_PI ;
// 				tri->getState().getField( PRINCIPAL_STRESS_ANGLE_FIELD, *tri->third, v, false ) ;
// 				ret[0] = 180.*v[0]/M_PI ;

            found = true ;
            break ;
        }
        case TWFT_ENRICHMENT:
        {
            ret[2] = tri->getEnrichmentFunctions().size() ;
            ret[1] = tri->getEnrichmentFunctions().size() ;
            ret[0] = tri->getEnrichmentFunctions().size() ;
            found = true ;
            break ;
        }
        case TWFT_INTERSECTION:
        {
            if(intersection)
            {
                bool inter = tri->intersects(intersection) ;
                ret[2] = inter ;
                ret[1] = inter ;
                ret[0] = inter ;
                found = true ;
            }
            else
            {
                ret[2] = 0 ;
                ret[1] = 0 ;
                ret[0] = 0 ;
                found = true ;
            }
            break ;
        }
        case TWFT_CRACK_ANGLE:
        {
            if( tri->getBehaviour()->getFractureCriterion() && tri->getBehaviour()->getDamageModel()->getState().max() > POINT_TOLERANCE)
            {
                double angle = 180.*tri->getBehaviour()->getFractureCriterion()->getSmoothedField(PRINCIPAL_STRESS_ANGLE_FIELD,  tri->getState())[0]/M_PI ;
                while (angle < 0)
                    angle += 90 ;
                ret[2] =  angle ;
                ret[1] =  angle ;
                ret[0] =  angle ;
            }
            else
            {
                ret[2] = 0 ;
                ret[1] = 0 ;
                ret[0] = 0 ;
            }

            found = true ;
            break ;
        }
        case TWFT_CRITERION:
        {
            if( tri->getBehaviour() && tri->getBehaviour()->getFractureCriterion())
            {
                double d = tri->getBehaviour()->getFractureCriterion()->getScoreAtTimeStepEnd() ;
                ret[2] = d ;
                ret[1] = d ;
                ret[0] = d ;
            }
            else if( tri->getBehaviour() )
            {
                double d = -1 ;
                ret[2] = d ;
                ret[1] = d ;
                ret[0] = d ;
            }

            found = true ;
            break ;
        }
        case TWFT_PRINCIPAL_STRESS:
        {
            Vector v(0., 2) ;
            tri->getState().getField( PRINCIPAL_REAL_STRESS_FIELD, *tri->first, v, false, 0) ;
            ret[5] = v[1] ;
            ret[2] = v[0] ;
            tri->getState().getField( PRINCIPAL_REAL_STRESS_FIELD, *tri->second, v, false, 0) ;
            ret[4] = v[1] ;
            ret[1] = v[0] ;
            tri->getState().getField( PRINCIPAL_REAL_STRESS_FIELD, *tri->third, v, false, 0) ;
            ret[3] = v[1] ;
            ret[0] = v[0] ;
            found = true ;
            break ;
        }
        case TWFT_PRINCIPAL_STRAIN:
        {
            Vector v(0., 2) ;
            tri->getState().getField( PRINCIPAL_TOTAL_STRAIN_FIELD, *tri->first, v, false, 0) ;
            ret[5] = v[1] ;
            ret[2] = v[0] ;
            tri->getState().getField( PRINCIPAL_TOTAL_STRAIN_FIELD, *tri->second, v, false, 0) ;
            ret[4] = v[1] ;
            ret[1] = v[0] ;
            tri->getState().getField( PRINCIPAL_TOTAL_STRAIN_FIELD, *tri->third, v, false, 0) ;
            ret[3] = v[1] ;
            ret[0] = v[0] ;
            found = true ;
            break ;
        }
        case TWFT_STIFFNESS:
        {
            if(tri->getBehaviour()->fractured())
            {
                ret[2] = 0. ;
                ret[1] = 0. ;
                ret[0] = 0. ;
                found = true ;
                break ;
            }


            double t = 0 ;
            size_t n = tri->getBoundingPoints().size()/3 ;
            if(tri->timePlanes() > 1)
                n /= tri->timePlanes() ;

            if( tri->timePlanes() > 1 )
                t = -1 + timePlane[0] * 2 / ( (int)tri->timePlanes() - 1 ) ;

            Point A = tri->getBoundingPoint(0) ;
            Point B = tri->getBoundingPoint(n) ;
            Point C = tri->getBoundingPoint(2*n) ;
            Point A_ = tri->inLocalCoordinates(A) ;
            A_.getT() = t ;
            Point B_ = tri->inLocalCoordinates(B) ;
            B_.getT() = t ;
            Point C_ = tri->inLocalCoordinates(C) ;
            C_.getT() = t ;

            ret[2] = (tri->getBehaviour()->getTensor( A_ )[0][0]+tri->getBehaviour()->getTensor( A_ )[1][1])*.5 ;
            ret[1] = (tri->getBehaviour()->getTensor( B_ )[0][0]+tri->getBehaviour()->getTensor( B_ )[1][1])*.5 ;
            ret[0] = (tri->getBehaviour()->getTensor( C_ )[0][0]+tri->getBehaviour()->getTensor( C_ )[1][1])*.5 ;
            found = true ;

            break ;
        }
        case TWFT_VISCOSITY:
        {
            double t = 0 ;
            size_t n = tri->getBoundingPoints().size()/3 ;
            if(tri->timePlanes() > 1)
                n /= tri->timePlanes() ;

            if( tri->timePlanes() > 1 )
                t = -1 + timePlane[0] * 2 / ( (int)tri->timePlanes() - 1 ) ;

            Point A = tri->getBoundingPoint(0) ;
            Point B = tri->getBoundingPoint(n) ;
            Point C = tri->getBoundingPoint(2*n) ;
            Point A_ = tri->inLocalCoordinates(A) ;
            A_.getT() = t ;
            Point B_ = tri->inLocalCoordinates(B) ;
            B_.getT() = t ;
            Point C_ = tri->inLocalCoordinates(C) ;
            C_.getT() = t ;

            ret[2] = (tri->getBehaviour()->getViscousTensor( A_ )[3][3]+tri->getBehaviour()->getViscousTensor( A_ )[4][4])*.5 ;
            ret[1] = (tri->getBehaviour()->getViscousTensor( B_ )[3][3]+tri->getBehaviour()->getViscousTensor( B_ )[4][4])*.5 ;
            ret[0] = (tri->getBehaviour()->getViscousTensor( C_ )[3][3]+tri->getBehaviour()->getViscousTensor( C_ )[4][4])*.5 ;
            found = true ;

            break ;
        }
        case TWFT_STIFFNESS_X:
        {
            double t = 0 ;

            if( tri->timePlanes() > 1 )
                t = -1 + timePlane[0] * 2 / ( (int)tri->timePlanes() - 1 ) ;


            ret[2] = tri->getBehaviour()->getTensor( Point( 0.3333, 0.3333, 0.3333, t ) )[0][0] ;
            ret[1] = ret[2] ;
            ret[0] = ret[2] ;
            found = true ;


            break ;
        }
        case TWFT_STIFFNESS_Y:
        {
            double t = 0 ;

            if( tri->timePlanes() > 1 )
                t = -1 + timePlane[0] * 2 / ( (int)tri->timePlanes() - 1 ) ;


            ret[2] = tri->getBehaviour()->getTensor( Point( 0.3333, 0.3333, 0.3333, t ) )[1][1] ;
            ret[1] = ret[2] ;
            ret[0] = ret[2] ;
            found = true ;


            break ;
        }
        case TWFT_STIFFNESS_Z:
        {
            double t = 0 ;

            if( tri->timePlanes() > 1 )
                t = -1 + timePlane[0] * 2 / ( (int)tri->timePlanes() - 1 ) ;


            ret[2] = tri->getBehaviour()->getTensor( Point( 0.3333,0.3333, 0.3333, t ) )[2][2] ;
            ret[1] = ret[2] ;
            ret[0] = ret[2] ;
            found = true ;

            break ;
        }
        case TWFT_VON_MISES:
        {
            Vector v(0., 1) ;
            tri->getState().getField( VON_MISES_REAL_STRESS_FIELD, *tri->first, v, false , 0) ;
            ret[2] = ret[1] = ret[0] = v[0] ;
            found = true ;
            break ;
        }
        default:
            break ;
        }
    }

    return std::make_pair( found, ret ) ;

}

void TriangleWriter::writeHeader( int layer, bool append )
{
    std::fstream outstream ;

    if( append )
        outstream.open( filename.c_str(), std::ios::out | std::ios::app ) ;
    else
        outstream.open( filename.c_str(), std::ios::out ) ;

    outstream << "TRIANGLES" << std::endl ;
    outstream << ( int ) nTriangles[layerTranslator[layer]] << std::endl ;
    outstream << 3 << std::endl ;
    outstream << ( ( int ) values.back()[layerTranslator[layer]].size() - 6 ) / 3 << std::endl ;
    outstream.close() ;
}

void BinaryTriangleWriter::writeHeader(int layer, bool append )
{
    std::fstream outstream ;

    if( append )
        outstream.open( filename.c_str(), std::ios::out | std::ios::app ) ;
    else
        outstream.open( filename.c_str(), std::ios::out ) ;

    outstream << "BIN_TRIANGLES" << std::endl ;
    outstream << ( int ) values.back()[0].size() << std::endl ;
    outstream << 3 << std::endl ;
    outstream << ( ( int ) values.back().size() - 6 ) / 3 << std::endl ;
    outstream.close() ;
}

void MultiTriangleWriter::writeIndexFile()
{
    std::fstream headstream ;

    headstream.open( head.c_str(), std::ios::out | std::ios::app ) ;
    headstream << filename << std::endl ;
    headstream.close();
}

void MultiTriangleWriter::writeHeader( int layer, bool firstlayer, bool append )
{
    if(firstlayer)
    {
        nextCounter() ;
    }

    std::fstream outstream ;
    if( append )
        outstream.open( filename.c_str(), std::ios::out | std::ios::app ) ;
    else
        outstream.open( filename.c_str(), std::ios::out ) ;
    outstream << "TRIANGLES" << std::endl ;
    outstream << ( int ) nTriangles[layerTranslator[layer]] << std::endl ;
    outstream << 3 << std::endl ;
    outstream << ( ( int ) values.back()[layerTranslator[layer]].size() - 6 ) / 3 << std::endl ;
    outstream.close() ;
}



int numberOfFields( TWFieldType field, size_t index , std::map<size_t, FieldType> & fieldsOther, std::map<size_t, std::string> & external)
{
    switch( field )
    {
    case TWFT_COORDINATE:
        return 6 ;
    case TWFT_DISPLACEMENTS:
        return 6 ;
    case TWFT_SCALAR:
        return 3 ;
    case TWFT_DOH:
        return 3 ;
    case TWFT_PRINCIPAL_ANGLE:
        return 3 ;
    case TWFT_TRIANGLE_ANGLE:
        return 3 ;
    case TWFT_CRACK_ANGLE:
        return 3 ;
    case TWFT_STIFFNESS:
        return 3 ;
    case TWFT_VISCOSITY:
        return 3 ;
    case TWFT_STIFFNESS_X:
        return 3 ;
    case TWFT_STIFFNESS_Y:
        return 3 ;
    case TWFT_STIFFNESS_Z:
        return 3 ;
    case TWFT_INTERSECTION:
        return 3 ;
        /*		case TWFT_STRAIN:
        			return 9 ;*/
    case TWFT_PRINCIPAL_STRAIN:
        return 6 ;
    case TWFT_PRINCIPAL_STRESS:
        return 6 ;
        /*		case TWFT_STRAIN_AND_STRESS:
        			return 18 ;
        		case TWFT_STRESS:
        			return 9 ;*/
    case TWFT_GRADIENT:
        return 6 ;
    case TWFT_GRADIENT_AND_FLUX:
        return 12 ;
    case TWFT_FLUX:
        return 6 ;
    case TWFT_VON_MISES:
        return 3 ;
    case TWFT_ENRICHMENT:
        return 3 ;
    case TWFT_CRITERION:
        return 3 ;
    case TWFT_DAMAGE:
        return 3 ;
    case TWFT_IMPOSED_STRESS_NORM:
        return 3 ;
    case TWFT_LARGE_DEFORMATION_TRANSFORM:
        return 3 ;
    case TWFT_LARGE_DEFORMATION_ANGLE:
        return 3 ;
    case TWFT_CRACKS:
        return 6 ;
    case TWFT_FIELD_TYPE:
        return numberOfFields(fieldsOther[index]) ;
    case TWFT_INTERNAL_VARIABLE:
    {
        return 3 ;
    }
    }

    return 3 ;
}

std::string nameOfField(TWFieldType field, size_t index , std::map<size_t, FieldType> & fieldsOther,  std::map<size_t, std::string> & external)
{
    switch( field )
    {
    case TWFT_COORDINATE:
        return std::string("Coordinates") ;
    case TWFT_DISPLACEMENTS:
        return std::string("Displacements") ;
    case TWFT_SCALAR:
        return std::string("Scalar Field") ;
    case TWFT_DOH:
        return std::string("Degree Of Hydration") ;
    case TWFT_PRINCIPAL_ANGLE:
        return std::string("Angle of Principal Stresses") ;
    case TWFT_TRIANGLE_ANGLE:
        return std::string("Triangle angle") ;
    case TWFT_CRACK_ANGLE:
        return std::string("Cracking Angle") ;
    case TWFT_STIFFNESS:
        return std::string("Stiffness") ;
    case TWFT_VISCOSITY:
        return std::string("Viscosity") ;
    case TWFT_STIFFNESS_X:
        return std::string("Stiffness X") ;
    case TWFT_STIFFNESS_Y:
        return std::string("Stiffness Y") ;
    case TWFT_STIFFNESS_Z:
        return std::string("Stiffness Z") ;
    case TWFT_LARGE_DEFORMATION_TRANSFORM:
        return std::string("spin angle") ;
    case TWFT_ENRICHMENT:
        return std::string("Enrichment") ;
        /*		case TWFT_STRAIN:
        			return std::string("Strain") ;*/
    case TWFT_PRINCIPAL_STRAIN:
        return std::string("Principal Strain") ;
    case TWFT_PRINCIPAL_STRESS:
        return std::string("Principal Stress") ;
    case TWFT_INTERSECTION:
        return std::string("Intersection") ;
        /*		case TWFT_STRAIN_AND_STRESS:
        			return std::string("Strains and Stress") ;
        		case TWFT_STRESS:
        			return std::string("Stress") ;*/
    case TWFT_GRADIENT:
        return std::string("Gradient") ;
    case TWFT_GRADIENT_AND_FLUX:
        return std::string("Gradient and Flux") ;
    case TWFT_FLUX:
        return std::string("Flux") ;
    case TWFT_VON_MISES:
        return std::string("Von Mises Stress") ;
    case TWFT_CRITERION:
        return std::string("Criterion") ;
    case TWFT_DAMAGE:
        return std::string("Damage") ;
    case TWFT_IMPOSED_STRESS_NORM:
        return std::string("Imposed Stress") ;
    case TWFT_CRACKS:
        return std::string("Crack Opening") ;
    case TWFT_LARGE_DEFORMATION_ANGLE:
         return std::string("Large deformation angle") ;
    case TWFT_FIELD_TYPE:
        return nameOfField(fieldsOther[index]) ;
    case TWFT_INTERNAL_VARIABLE:
    {
        if(external.find(index) != external.end())
        {
            return external[index] ;
        }
        return std::string("Internal Variable") ;
    }
    }

    return std::string("Not a Field") ;
}

std::string nameOfField(FieldType field)
{
    switch( field )
    {
    case DISPLACEMENT_FIELD :
        return std::string("Displacements") ;
    case ENRICHED_DISPLACEMENT_FIELD :
        return std::string("Enriched Displacements") ;
    case SPEED_FIELD :
        return std::string("") ;
    case FLUX_FIELD :
        return std::string("") ;
    case GRADIENT_FIELD :
        return std::string("") ;
    case TOTAL_STRAIN_FIELD :
        return std::string("Strain") ;
    case MECHANICAL_STRAIN_FIELD :
        return std::string("Mechanical Strain") ;
    case STRAIN_RATE_FIELD :
        return std::string("Strain Rate") ;
    case EFFECTIVE_STRESS_FIELD :
        return std::string("Stress") ;
    case REAL_STRESS_FIELD :
        return std::string("Stress") ;
    case IMPOSED_STRAIN_FIELD :
        return std::string("Imposed Strain") ;
    case IMPOSED_STRESS_FIELD :
        return std::string("Imposed Stress") ;
    case PRINCIPAL_IMPOSED_STRAIN_FIELD :
        return std::string("Principal Imposed Strain") ;
    case PRINCIPAL_IMPOSED_STRESS_FIELD :
        return std::string("Principal Imposed Stress") ;
    case PRINCIPAL_TOTAL_STRAIN_FIELD :
        return std::string("Principal Strain") ;
    case PRINCIPAL_MECHANICAL_STRAIN_FIELD :
        return std::string("Principal Mechanical Strain") ;
    case PRINCIPAL_EFFECTIVE_STRESS_FIELD :
        return std::string("Principal Effective Stress") ;
    case PRINCIPAL_REAL_STRESS_FIELD :
        return std::string("Principal Stress") ;
    case NON_ENRICHED_STRAIN_FIELD :
        return std::string("") ;
    case NON_ENRICHED_STRAIN_RATE_FIELD :
        return std::string("") ;
    case NON_ENRICHED_EFFECTIVE_STRESS_FIELD :
        return std::string("") ;
    case NON_ENRICHED_REAL_STRESS_FIELD :
        return std::string("") ;
    case VON_MISES_STRAIN_FIELD :
        return std::string("Von Mises Strain") ;
    case VON_MISES_REAL_STRESS_FIELD :
        return std::string("Von Mises Stress") ;
    case VON_MISES_EFFECTIVE_STRESS_FIELD :
        return std::string("") ;
    case PRINCIPAL_STRAIN_ANGLE_FIELD :
        return std::string("Principal Strain Angle") ;
    case PRINCIPAL_STRESS_ANGLE_FIELD :
        return std::string("Principal Stress Angle") ;
    case SCALAR_DAMAGE_FIELD :
        return std::string("Scalar Damage") ;
    case TENSOR_DAMAGE_FIELD :
        return std::string("Tensor Damage") ;
    case INTERNAL_VARIABLE_FIELD :
        return std::string("") ;
    case GENERALIZED_VISCOELASTIC_DISPLACEMENT_FIELD :
        return std::string("") ;
    case GENERALIZED_VISCOELASTIC_ENRICHED_DISPLACEMENT_FIELD :
        return std::string("") ;
    case GENERALIZED_VISCOELASTIC_SPEED_FIELD :
        return std::string("") ;
    case GENERALIZED_VISCOELASTIC_STRAIN_FIELD :
        return std::string("") ;
    case GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD :
        return std::string("") ;
    case GENERALIZED_VISCOELASTIC_EFFECTIVE_STRESS_FIELD :
        return std::string("") ;
    case GENERALIZED_VISCOELASTIC_REAL_STRESS_FIELD :
        return std::string("") ;
    case GENERALIZED_VISCOELASTIC_PRINCIPAL_STRAIN_FIELD :
        return std::string("") ;
    case GENERALIZED_VISCOELASTIC_PRINCIPAL_EFFECTIVE_STRESS_FIELD :
        return std::string("") ;
    case GENERALIZED_VISCOELASTIC_PRINCIPAL_REAL_STRESS_FIELD :
        return std::string("") ;
    case GENERALIZED_VISCOELASTIC_NON_ENRICHED_STRAIN_FIELD :
        return std::string("") ;
    case GENERALIZED_VISCOELASTIC_NON_ENRICHED_STRAIN_RATE_FIELD :
        return std::string("") ;
    case GENERALIZED_VISCOELASTIC_NON_ENRICHED_EFFECTIVE_STRESS_FIELD :
        return std::string("") ;
    case GENERALIZED_VISCOELASTIC_NON_ENRICHED_REAL_STRESS_FIELD :
        return std::string("") ;
    }
    return std::string("Not a Field") ;
}

int numberOfFields( FieldType field )
{
    switch( field )
    {
    case DISPLACEMENT_FIELD :
        return 6 ;
    case ENRICHED_DISPLACEMENT_FIELD :
        return 6 ;
    case SPEED_FIELD :
        return 6 ;
    case FLUX_FIELD :
        return 6 ;
    case GRADIENT_FIELD :
        return 6 ;
    case TOTAL_STRAIN_FIELD :
        return 9 ;
    case MECHANICAL_STRAIN_FIELD :
        return 9 ;
    case STRAIN_RATE_FIELD :
        return 9 ;
    case EFFECTIVE_STRESS_FIELD :
        return 9 ;
    case IMPOSED_STRAIN_FIELD :
        return 9 ;
    case IMPOSED_STRESS_FIELD :
        return 9 ;
    case PRINCIPAL_IMPOSED_STRAIN_FIELD :
        return 6 ;
    case PRINCIPAL_IMPOSED_STRESS_FIELD :
        return 6 ;
    case REAL_STRESS_FIELD :
        return 9 ;
    case PRINCIPAL_TOTAL_STRAIN_FIELD :
        return 6 ;
    case PRINCIPAL_MECHANICAL_STRAIN_FIELD :
        return 6 ;
    case PRINCIPAL_EFFECTIVE_STRESS_FIELD :
        return 6 ;
    case PRINCIPAL_REAL_STRESS_FIELD :
        return 6 ;
    case NON_ENRICHED_STRAIN_FIELD :
        return 9 ;
    case NON_ENRICHED_STRAIN_RATE_FIELD :
        return 9 ;
    case NON_ENRICHED_EFFECTIVE_STRESS_FIELD :
        return 9;
    case NON_ENRICHED_REAL_STRESS_FIELD :
        return 9 ;
    case VON_MISES_STRAIN_FIELD :
        return 3 ;
    case VON_MISES_REAL_STRESS_FIELD :
        return 3 ;
    case VON_MISES_EFFECTIVE_STRESS_FIELD :
        return 3 ;
    case PRINCIPAL_STRAIN_ANGLE_FIELD :
        return 3 ;
    case PRINCIPAL_STRESS_ANGLE_FIELD :
        return 3 ;
    case SCALAR_DAMAGE_FIELD :
        return 3 ;
    case TENSOR_DAMAGE_FIELD :
        return 6 ;
    case INTERNAL_VARIABLE_FIELD :
        return 3 ;
    case GENERALIZED_VISCOELASTIC_DISPLACEMENT_FIELD :
        return 6 ;
    case GENERALIZED_VISCOELASTIC_ENRICHED_DISPLACEMENT_FIELD :
        return 3 ;
    case GENERALIZED_VISCOELASTIC_SPEED_FIELD :
        return 6 ;
    case GENERALIZED_VISCOELASTIC_STRAIN_FIELD :
        return 9 ;
    case GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD :
        return 9 ;
    case GENERALIZED_VISCOELASTIC_EFFECTIVE_STRESS_FIELD :
        return 9;
    case GENERALIZED_VISCOELASTIC_REAL_STRESS_FIELD :
        return 9 ;
    case GENERALIZED_VISCOELASTIC_PRINCIPAL_STRAIN_FIELD :
        return 6 ;
    case GENERALIZED_VISCOELASTIC_PRINCIPAL_EFFECTIVE_STRESS_FIELD :
        return 6 ;
    case GENERALIZED_VISCOELASTIC_PRINCIPAL_REAL_STRESS_FIELD :
        return 9 ;
    case GENERALIZED_VISCOELASTIC_NON_ENRICHED_STRAIN_FIELD :
        return 9 ;
    case GENERALIZED_VISCOELASTIC_NON_ENRICHED_STRAIN_RATE_FIELD :
        return 9 ;
    case GENERALIZED_VISCOELASTIC_NON_ENRICHED_EFFECTIVE_STRESS_FIELD :
        return 9 ;
    case GENERALIZED_VISCOELASTIC_NON_ENRICHED_REAL_STRESS_FIELD :
        return 9 ;
    }
    return 3. ;
}

}



























