// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2009-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "../features/features.h"
#include "../physics/stiffness_with_imposed_deformation.h"
#include "../features/pore.h"
#include "../features/sample.h"
#include "../features/inclusion.h"
#include "../features/crack.h"
#include "../utilities/itoa.h"
#include "../utilities/writer/triangle_writer.h"

#include <fstream>

#include <cmath>
#include <typeinfo>
#include <limits>
#include <time.h> 
#define DEBUG 

#include "../features/sample3d.h"

using namespace Amie ;
using namespace std;

FeatureTree * featureTree ;

double width = 1. ;
double height = 1. ;
Sample sample(nullptr, width, height, 0.5, 0.5) ;	


std::string outfilename("results") ;
std::fstream outresultfile ;
int totit = 0 ;

Vector step(int nInc)
{
	featureTree->setMaxIterationsPerStep(50) ;
	featureTree->step() ;

	std::string filename("triangles") ;
	filename.append(itoa(totit++, 10)) ;
	std::cout << filename << std::endl ;
	std::fstream outfile  ;
	outfile.open(filename.c_str(), std::ios::out) ;
	
	TriangleWriter writer(filename, featureTree) ;
	writer.getField(REAL_STRESS_FIELD ) ;
	writer.getField(STRAIN_FIELD ) ;
	writer.getField(TWFT_VON_MISES) ;
	writer.getField(TWFT_STIFFNESS) ;
	writer.write() ;

	Vector avge = featureTree->getAverageField(STRAIN_FIELD) ;
	Vector avgs = featureTree->getAverageField(REAL_STRESS_FIELD) ;
	Vector results(6) ;
	results[0] = avge[0] ;
	results[1] = avge[1] ;
	results[2] = avge[2] ;
	results[3] = avgs[0] ;
	results[4] = avgs[1] ;
	results[5] = avgs[2] ;
	return results ;
}

Matrix buildStiffnessMatrix(Vector p)
{
	Matrix prop(3,3) ;
	double E = p[1] * 1e9 ;
	double nu = p[2] ;
	prop[0][0] = E/(1.-nu*nu) ; prop[0][1] =E/(1.-nu*nu)*nu ; prop[0][2] = 0 ;
	prop[1][0] = E/(1.-nu*nu)*nu ; prop[1][1] = E/(1.-nu*nu) ; prop[1][2] = 0 ; 
	prop[2][0] = 0 ; prop[2][1] = 0 ; prop[2][2] = E/(1.-nu*nu)*(1.-nu)/2. ; 
	return prop ;  
}

std::vector<Inclusion *> buildInclusion(int nInc, double pInc)
{
	std::vector<Inclusion *> inc ;
	double radius = std::sqrt(pInc / (nInc * nInc * 3.1415926)) ;
	std::cout << radius << std::endl ;
	double dist = (1 - nInc * radius * 2) / (nInc + 1) ;
	std::cout << dist << std::endl ;
	for(int i = 0 ; i < nInc ; i++)
	{
		for(int j = 0 ; j < nInc ; j++)
		{
			inc.push_back(new Inclusion(radius,(i +1) * dist + radius * (2 * i + 1),(j +1) * dist  + radius * (2 * j + 1))) ;
//			std::cout << (i +1) * dist + radius * (2 * i + 1) << ";" << (j +1) * dist  + radius * (2 * j + 1) << std::endl ;
		}
	}
	return inc ;
}


Vector thermalDeformation(Vector matrixProp, Vector inclusionProp, int nInc, double pInc, double temperature)
{
	Vector result(6) ;
	FeatureTree F(&sample) ;
	featureTree = &F ;
	
	Matrix mK = buildStiffnessMatrix(matrixProp) ;
	Matrix iK = buildStiffnessMatrix(inclusionProp) ;
	double mA = matrixProp[0] * 1e-6 ;
	double iA = inclusionProp[0] * 1e-6 ;
	
	Vector mDef(3) ;
	mDef[0] = mA * temperature ;
	mDef[1] = mA * temperature ;
	mDef[2] = mA * temperature ;
	Vector iDef(3) ;
	iDef[0] = iA * temperature ;
	iDef[1] = iA * temperature ;
	iDef[2] = iA * temperature ;
	
	StiffnessWithImposedDeformation * mTDef = new StiffnessWithImposedDeformation(mK,mDef) ;
	StiffnessWithImposedDeformation * iTDef = new StiffnessWithImposedDeformation(iK,iDef) ;
	
	std::vector<Inclusion *> inc = buildInclusion(nInc,pInc) ;
	
	sample.setBehaviour(mTDef) ;
	for(size_t i = 0 ; i < inc.size() ; i++)
	{
		  inc[i]->setBehaviour(iTDef) ;
		  F.addFeature(&sample, inc[i]) ;
	}
	F.setSamplingNumber(1024) ;
	F.setOrder(LINEAR) ;

	result = step(nInc) ;
	return result ;
		
}


int main(int argc, char *argv[])
{
  
	outresultfile.open(outfilename.c_str(), std::ios::out) ;
	Vector exp(6) ;
	
	Vector paste_prop(3) ;
	Vector similipaste_prop(3) ;
	Vector quartz_prop(3) ;
	Vector limestone_prop(3) ;

	paste_prop[0] = 9.69 ;
	paste_prop[1] = 25 ;
	paste_prop[2] = 0.22 ;

	similipaste_prop[0] = 4.34 ;
	similipaste_prop[1] = 25 ;
	similipaste_prop[2] = 0.22 ;

	quartz_prop[0] = 11.38 ;
	quartz_prop[1] = 70 ;
	quartz_prop[2] = 0.2 ;

	limestone_prop[0] = 4.34 ;
	limestone_prop[1] = 50 ;
	limestone_prop[2] = 0.22 ;
    
	exp = thermalDeformation(paste_prop,similipaste_prop,1,0.75,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;
	exp = thermalDeformation(paste_prop,similipaste_prop,3,0.75,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;
	exp = thermalDeformation(paste_prop,similipaste_prop,10,0.75,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;
	exp = thermalDeformation(paste_prop,similipaste_prop,33,0.75,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;
	exp = thermalDeformation(paste_prop,quartz_prop,1,0.75,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;
	exp = thermalDeformation(paste_prop,quartz_prop,3,0.75,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;
	exp = thermalDeformation(paste_prop,quartz_prop,10,0.75,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;	
	exp = thermalDeformation(paste_prop,quartz_prop,33,0.75,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;
 	exp = thermalDeformation(paste_prop,limestone_prop,1,0.75,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;
	exp = thermalDeformation(paste_prop,limestone_prop,3,0.75,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;
	exp = thermalDeformation(paste_prop,limestone_prop,10,0.75,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;	
	exp = thermalDeformation(paste_prop,limestone_prop,33,0.75,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;
	exp = thermalDeformation(paste_prop,similipaste_prop,1,0.5,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;
	exp = thermalDeformation(paste_prop,similipaste_prop,3,0.5,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;
	exp = thermalDeformation(paste_prop,similipaste_prop,10,0.5,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;	
	exp = thermalDeformation(paste_prop,similipaste_prop,33,0.5,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;
	exp = thermalDeformation(paste_prop,quartz_prop,1,0.5,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;
	exp = thermalDeformation(paste_prop,quartz_prop,3,0.5,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;
	exp = thermalDeformation(paste_prop,quartz_prop,10,0.5,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;	
	exp = thermalDeformation(paste_prop,quartz_prop,33,0.5,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;
	exp = thermalDeformation(paste_prop,limestone_prop,1,0.5,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;
	exp = thermalDeformation(paste_prop,limestone_prop,3,0.5,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;
	exp = thermalDeformation(paste_prop,limestone_prop,10,0.5,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;	
	exp = thermalDeformation(paste_prop,limestone_prop,33,0.5,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;
	exp = thermalDeformation(paste_prop,similipaste_prop,1,0.25,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;
	exp = thermalDeformation(paste_prop,similipaste_prop,3,0.25,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;
	exp = thermalDeformation(paste_prop,similipaste_prop,10,0.25,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;	
	exp = thermalDeformation(paste_prop,similipaste_prop,33,0.25,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;
	exp = thermalDeformation(paste_prop,quartz_prop,1,0.25,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;
	exp = thermalDeformation(paste_prop,quartz_prop,3,0.25,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;
	exp = thermalDeformation(paste_prop,quartz_prop,10,0.25,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;	
	exp = thermalDeformation(paste_prop,quartz_prop,33,0.25,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;
	exp = thermalDeformation(paste_prop,limestone_prop,1,0.25,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;
	exp = thermalDeformation(paste_prop,limestone_prop,3,0.25,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;
	exp = thermalDeformation(paste_prop,limestone_prop,10,0.25,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;	
	exp = thermalDeformation(paste_prop,limestone_prop,33,0.25,50) ;
	outresultfile << exp[0] << " ; " << exp[1] << " ; " << exp[2] << " ; " << exp[3] << " ; " << exp[4] << " ; " << exp[5] << "\n" ;
	return 0 ;
}
