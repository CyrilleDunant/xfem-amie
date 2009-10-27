//
// C++ Implementation: main_jerome
//
// Description: 
//
//
// Author:  Jérôme Krebs <jerome.krebs@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include <stdlib.h>
#include <GL/glut.h>
#include "utilities/granulo.h"
#include <vector>
#include <time.h> //utilisation du module de temps pour initialiser la fonction rand
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include "../geometry/geometry_2D.h"
#include "../geometry/geometry_base.h"
#include "../features/sample.h"
#include "../features/inclusion.h"
#include "../utilities/placement.cpp"

using namespace Mu ;

std::vector <Inclusion *> inclusions;
double longueurX = 200.;
double longueurY = 200.;
std::vector<DelaunayTriangle *> triangles ;
int nombreGranulatsPlaces;




void reshape(int width, int height) {

	// Set the new viewport size
// 	glViewport(-250, -250, (GLint)width, (GLint)height);

	glViewport(-900, -900, 1800, 1800);
	// Choose the projection matrix to be the matrix 
	// manipulated by the following calls
//	glMatrixMode(GL_PROJECTION);

	// Set the projection matrix to be the identity matrix
//	glLoadIdentity();

	// Define the dimensions of the Orthographic Viewing Volume
//	glOrtho(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

	// Choose the modelview matrix to be the matrix
	// manipulated by further calls
//	glMatrixMode(GL_MODELVIEW);

	// Clear the window
	glClear(GL_COLOR_BUFFER_BIT);
}

void draw(void) {

	// Set the drawing color
	
	glColor3f(1.0, 1.0, 1.0);
	for (size_t i=0; i<nombreGranulatsPlaces; i++)
	{
		
	
		for(int j =0; j<360; j++)
		{
		
			double x1=inclusions[i]->getCenter().x/longueurX+(inclusions[i]->getRadius()/longueurX)*cos(j);
			//double x1=inclusions/50.+(inclusions[i]->getRadius()/50.)*cos(j);
			double y1=inclusions[i]->getCenter().y/longueurY+(inclusions[i]->getRadius()/longueurY)*sin(j);
			glBegin(GL_POINTS);
			glVertex2f(x1, y1);
		
			glEnd();
		}
		
	}

// 	glColor3f(1.0, 1.0, 1.0);
// 	for(size_t i= 0; i<triangles.size(); i++)
// 	{
// 		glBegin(GL_LINE_LOOP);
// 		glVertex2f(triangles[i]->first->x/longueurX,triangles[i]->first->y/longueurY);
// 		glVertex2f(triangles[i]->second->x/longueurX,triangles[i]->second->y/longueurY);
// 		glVertex2f(triangles[i]->third->x/longueurX,triangles[i]->third->y/longueurY);
// 		glEnd();
// 	}
// 	
	// Flush the buffer to force drawing of all objects thus far
	glFlush();
}



void init(void) 
{
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glShadeModel(GL_SMOOTH);   // Enables Smooth Shading
	glEnable(GL_LINE_SMOOTH) ;
	glEnable(GL_POLYGON_SMOOTH);
	glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST) ;
	
// 	glPointSize(std::max(0.4*(double)width()/(double)columns, 1.));
	glClearColor(0.0f,0.0f,0.0f,0.0f);                                      // Black Background

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA,GL_ONE); 
} 



double pourcentAgregatsSurfaceOccupe(double longueurX, double longueurY, std::vector<Inclusion *> inclusions)
{
	double pourcentAgregats;
	double aireTotaleAgregats=0;
	for(size_t i=0; i<nombreGranulatsPlaces; i++)
	{
		aireTotaleAgregats += inclusions[i]->area();
	}
	pourcentAgregats = 100*aireTotaleAgregats/(longueurX*longueurY);
	std::cout<<"Le pourcentage d'agregat occupé par rapport au volume disponible vaut :"<<pourcentAgregats<<"%"<<std::endl;

return pourcentAgregats;
}


// void distributionGranulatsX(double nombreGranulatsPlaces, std::vector<Inclusion *> inclusions)
// {
// 	file<<"r , ";
// 	for(size_t i=0; i<nombreGranulatsPlaces; i++)
// 	{ 
// 		if(inclusions[i]->getCenter()->y>15 && inclusions[i]->getCenter()->y<35)
// 		{
// // 			std::cout<<inclusions[i]->getRadius()<<"   "<<inclusions[i]->getCenter()->x<<std::endl;
// // 			file << inclusions[i]->getRadius() << inclusions[i]->getCenter()->x<< std::endl ;
// 			file << inclusions[i]->getRadius() << " , " ;
// 		}
// 	}
// 	file<<std::endl<<"x , ";
// 	for(size_t i=0; i<nombreGranulatsPlaces; i++)
// 	{ 
// 		if(inclusions[i]->getCenter()->y>15 && inclusions[i]->getCenter()->y<35)
// 		{
// // 			std::cout<<inclusions[i]->getRadius()<<"   "<<inclusions[i]->getCenter()->x<<std::endl;
// 			file <<  inclusions[i]->getCenter()->x<< " , " ;
// 		}
// 	}
// 	file<<std::endl;
// }



int main(int argc, char **argv) 
{
	
	size_t nombreEssai =1;
	std::ostringstream basename;
	basename << "BOLOME_A_esssss.csv" ;
	std::ofstream file(basename.str().c_str()) ;
	

	for(size_t a =1;a <=nombreEssai; a++)
	{
// 		std::ostringstream basename;
// 		basename << "BOLOME_A_50_def_" ;
// 		basename << a ;
// 		basename << ".csv" ;
// 		std::ofstream file(basename.str().c_str()) ;
	
	
		double masseInitiale = 3000.;//25x25:28; 50x50:221; 100x100:1775; 200x200:14200
		double densite = 0.0268;
		double rayonGranulatMax = 8.;
		int triesMax = 1000;
	
		double pourcentMasseMin =0.10;//BOLOME_A, 200x200 : 0.1
	
	
	
// 		Granulo g(0.001, 1., masseInitiale,densite) ; //(b,n)		Granulo g(0.015, 1.8, masseInitiale,densite) ; //(b,n)

		GranuloBolome g(masseInitiale, densite, BOLOME_D);
	
// 		std::cout<<"Type de granulo : "<<"BOLOME_B"<<std::endl;
// 		std::cout<<"masse initiale : "<<masseInitiale<<" g"<<std::endl;	
// 		std::cout<<"densité granulat : "<<densite<<" g/mm3"<<std::endl;
// 		std::cout<<"rayon granulat max : "<<rayonGranulatMax<<" mm"<<std::endl;
// 		std::cout<<"longueur X : "<<longueurX<<" mm"<<std::endl;
// 		std::cout<<"longueur Y : "<<longueurY<<" mm"<<std::endl;
// 		std::cout<<"pourcentage de masse min : "<<pourcentMasseMin<<" %"<<std::endl;
// 	
		inclusions = g(rayonGranulatMax, pourcentMasseMin);//masseInitiale, densite
		size_t seed = time(NULL) ;
		int nombreGranulatsGeneres = inclusions.size();
		std::cout<<"plus petit granulat généré " << inclusions.back()->getRadius() <<std::endl;


		srand(time(NULL)); //change la fonction rand en fonction du temps
		


// 		std::vector < Inclusion *> inclusionsPasPlacees;
// 		inclusionsPasPlacees=inclusions;
// 		inclusions.clear();
// 		for(size_t i=0; i<nombreGranulatsGeneres;i++)
// 		{
// 			if(inclusionsPasPlacees[i]->getRadius()<0.46)
// 			{
// 				inclusions.push_back(inclusionsPasPlacees[i]);
// 			}
// 		}
		

	
		inclusions=placement(longueurX, longueurY, inclusions, &nombreGranulatsPlaces, triesMax);
	
// 		std::cout << std::endl ;
// 		std::cout<<"nombre granulats places "<<nombreGranulatsPlaces<<std::endl;
		std::cout<<"rayon plus petit granulat "<<inclusions[nombreGranulatsPlaces]->getRadius()<<std::endl;
	
	// 	distributionGranulatsX( nombreGranulatsPlaces, inclusions);//distribution de la taille des granulats le long de x
	
	 	
	
		
		
		Sample s(longueurX, longueurY, longueurX/2, longueurY/2) ; //initialise le carré de base de dimension longueurX et longueurY et de centre 0,0
	
		FeatureTree ft(&s) ;//initialise le problème
		std::vector<Inclusion *> inclusionsPlacees;
		double cimentE = 1e9;
		double cimentNu = 0.2;
		Matrix m0(3,3) ;
		m0[0][0] = cimentE/(1-cimentNu*cimentNu) ; m0[0][1] =cimentE/(1-cimentNu*cimentNu)*cimentNu ; m0[0][2] = 0 ;
		m0[1][0] = cimentE/(1-cimentNu*cimentNu)*cimentNu ; m0[1][1] = cimentE/(1-cimentNu*cimentNu) ; m0[1][2] = 0 ; 
		m0[2][0] = 0 ; m0[2][1] = 0 ; m0[2][2] = cimentE/(1-cimentNu*cimentNu)*(1.-cimentNu)/2. ; 
		s.setBehaviour(new Stiffness(m0)) ; //initialise la rigidité de la matrice
	
		double granulatE = 60e9;
		double granulatNu = 0.2;
		m0[0][0] = granulatE/(1-granulatNu*granulatNu) ; m0[0][1] =granulatE/(1-granulatNu*granulatNu)*granulatNu ; m0[0][2] = 0 ;
		m0[1][0] = granulatE/(1-granulatNu*granulatNu)*granulatNu ; m0[1][1] = granulatE/(1-granulatNu*granulatNu) ; m0[1][2] = 0 ; 
		m0[2][0] = 0 ; m0[2][1] = 0 ; m0[2][2] = granulatE/(1-granulatNu*granulatNu)*(1.-granulatNu)/2. ; 

		inclusionsPlacees.push_back(inclusions[0]);
	
		ft.addFeature(&s, inclusionsPlacees[0]) ;
	
		inclusionsPlacees[0]->setBehaviour(new Stiffness(m0)) ; //initialise la rigidité des granulats
		for(size_t h=1;h<nombreGranulatsPlaces;h++)
		{
			inclusionsPlacees.push_back(inclusions[h]);
	// 		Inclusion *i = new Inclusion(inclusions[h], x[h],y[h]) ;//initialise une inclusion i de rayon 20 et de centre 0,0
			ft.addFeature(inclusionsPlacees[h-1], inclusionsPlacees[h]) ;//ajoute au problème une inclusion dans le carré de base 
			inclusionsPlacees[h]->setBehaviour(new Stiffness(m0)) ;
		}
// 		
		double pourcentAgregatsAireOccupe = pourcentAgregatsSurfaceOccupe(longueurX,longueurY, inclusionsPlacees);
// 
		/**time*/
// 		if(a==1)
// 		{
// 			file << "\"essai en mode\" , " << "\"Type de granulo\" , " << "\"nombre point sur le périmètre\" , " << "\"e/c\" , " << "\"masse de ciment par m3 de béton\" , " <<  "\"nombre d'essai\" , " <<  "\"masse initiale \" , " << "\"densité granulat\" , " <<  "\"rayon granulat max\" , " << "\"module élastique du ciment\" , " <<  "\"coefficient de poisson du ciment\" , " <<  "\"module élastique des granulats\" , " <<  "\"coefficient de poisson des granulats\" , " << "\"longueur X\" , " << "\"longueur Y\" , " <<  "\"pourcentage de masse min\" , " << "\"nombre de granulats générés\" , " << "\"nombre de granulats placés\" , " <<  "\"rayon plus petit granulat\" , " << "\"pourcentage de l'aire occupé par granulats \" , " <<  "\"nombre de répétition d'essai de placement\" , " << "\"contrainte moyenne x\" , " << "\"contrainte moyenne y\" , " <<"\"contrainte moyenne xy\" , " <<  "\"déformation moyenne x\" , " << "\"déformation moyenne y\" , " <<  "\"déformation moyenne xy\" , " << "\"module de poisson\" , " << "\"module élastique réel\" , " <<  "\"module de cisaillement\" , " <<  "\"module élastique effectif\" , " << std::endl ;
// 		}
// 		file << "deformation" << " , "<< "BOLOME_A" << " , "<<	"128" << " , "<< "0.5" << " , " << "350 kg/m3" << " , " << nombreEssai << " , " << masseInitiale<<" , " << densite<<" , " << rayonGranulatMax<<" , " << cimentE<<" , " << cimentNu<<" , " << granulatE<<" , " << granulatNu<<" , " << longueurX<<" , " << longueurY<<" , " << pourcentMasseMin<<" , " << nombreGranulatsGeneres << " , " << nombreGranulatsPlaces << " , " << inclusionsPlacees[nombreGranulatsPlaces-1]->getRadius()<<" , " << pourcentAgregatsAireOccupe <<" , " << triesMax  ;
// 
// return 0;
// }
// }
/** */









/*

ft.sample(64) ;//nombre de point sur le périmètre du carré de base
		
		ft.generateElements() ;//génère les éléments
		
		triangles = ft.getTriangles() ;//stocke les coordonnées des triangles du maillage
	
	
		/**Applcation d'une déformation*/
	
		for(size_t i = 0 ; i< triangles.size() ; i++)//boucle pour les conditions limites et l'application d'une déformation
		{
	
			for(size_t j = 0 ; j < triangles[i]->getBoundingPoints().size() ;j++)//une boucle sur tous les points du triangle
			{
				if(triangles[i]->getBoundingPoint(j).x == longueurX) // condition limite est
				{
					ft.getAssembly()->setPointAlong(XI, -1, triangles[i]->getBoundingPoint(j).id) ;//compression de -1 en x 
				}
				if(triangles[i]->getBoundingPoint(j).x == 0)//condition limite ouest
				{
					ft.getAssembly()->setPointAlong(XI, 1, triangles[i]->getBoundingPoint(j).id) ; 
				}
				if(triangles[i]->getBoundingPoint(j).y == longueurY)//condition limite nord
				{
					ft.getAssembly()->setPointAlong(ETA, 0, triangles[i]->getBoundingPoint(j).id) ;
				}
				if(triangles[i]->getBoundingPoint(j).y == 0)//condition limite sud
				{
					ft.getAssembly()->setPointAlong(ETA, 0, triangles[i]->getBoundingPoint(j).id) ;
				}
			}
		}
	
	
		/**Applcation d'une contrainte*/
	/*
		for(size_t i = 0 ; i< triangles.size() ; i++)//boucle pour les conditions limites et l'application d'une contrainte
		{
	
			for(size_t j = 0 ; j < triangles[i]->getBoundingPoints()->size() ;j++)//une boucle sur tous les points du triangle
			{
				if(triangles[i]->getBoundingPoint(j)->x == longueurX) // condition limite est
				{
					ft.getAssembly()->setForceOn(XI, -1.0e9, triangles[i]->getBoundingPoint(j)->id) ;//compression de -1 en x et 0 en y-1.65e9
				}
				if(triangles[i]->getBoundingPoint(j)->x == 0)//condition limite ouest
				{
					ft.getAssembly()->setPointAlong(XI, 0, triangles[i]->getBoundingPoint(j)->id) ;
				}
// 				if(triangles[i]->getBoundingPoint(j)->y == longueurY)//condition limite nord
// 				{
// 					ft.getAssembly()->setPointAlong(XI, ETA, triangles[i]->getBoundingPoint(j)->id) ;
// 				}
				if(triangles[i]->getBoundingPoint(j)->y == 0)//condition limite sud
				{
					ft.getAssembly()->setPointAlong(ETA,0 , triangles[i]->getBoundingPoint(j)->id) ;
				}
			}
		}*/
	
	
		ft.step(0) ;//caclul les contraintes, déformations...
	// 	Vector stresses(triangles.size()*3);
		Point strainMoyenne(0,0,0);
		Point stressMoyenne(0,0,0);
	
// 		for(size_t i=0; i< triangles.size();i++)
// 		{
// 			
// 			Vector strain = triangles[i]->getState()->getStrain(triangles[i]->getCenter())	;
// 			//stresses[i*3] = stress[0] ;
// 			//stresses[i*3+1] = stress[1] ;
// 	// 		std::cout<<strain[0]<< "  " << strain[1]<< "  " << strain[2] <<std::endl;
// 			strainMoyenne.x += ((strain[0]*triangles[i]->area())/longueurX/longueurY);
// 			strainMoyenne.y += ((strain[1]*triangles[i]->area())/longueurX/longueurY);
// 			strainMoyenne.z += ((strain[2]*triangles[i]->area())/longueurX/longueurY);
// 	
// 			Vector stress = triangles[i]->getState()->getStress(triangles[i]->getCenter())	;
// 	// 		//stresses[i*3] = stress[0] ;
// 	// 		//stresses[i*3+1] = stress[1] ;
// 	// 		std::cout<<stress[0]<< "  " << stress[1]<< "  " << stress[2] <<std::endl;
// 			stressMoyenne.x += ((stress[0]*triangles[i]->area())/longueurX/longueurY);
// 			stressMoyenne.y += ((stress[1]*triangles[i]->area())/longueurX/longueurY);
// 			stressMoyenne.z += ((stress[2]*triangles[i]->area())/longueurX/longueurY);
// 		}

		std::vector <double> deformation(3*triangles.size());
		std::vector <double> contrainte(3*triangles.size());
		for(size_t i=0; i< triangles.size();i++)
		{
			
			Vector strain = triangles[i]->getState()->getStrain(triangles[i]->getCenter())	;
			deformation[i*3] = strain[0] ;
			deformation[i*3+1] = strain[1] ;
			deformation[i*3+2] = strain[2] ;
	// 		std::cout<<strain[0]<< "  " << strain[1]<< "  " << strain[2] <<std::endl;
			strainMoyenne.x += ((strain[0]*triangles[i]->area())/longueurX/longueurY);
			strainMoyenne.y += ((strain[1]*triangles[i]->area())/longueurX/longueurY);
			strainMoyenne.z += ((strain[2]*triangles[i]->area())/longueurX/longueurY);
	
			Vector stress = triangles[i]->getState()->getStress(triangles[i]->getCenter())	;
			contrainte[i*3] = stress[0] ;
			contrainte[i*3+1] = stress[1] ;
			contrainte[i*3+2] = stress[2] ;
	// 		std::cout<<stress[0]<< "  " << stress[1]<< "  " << stress[2] <<std::endl;
			stressMoyenne.x += ((stress[0]*triangles[i]->area())/longueurX/longueurY);
			stressMoyenne.y += ((stress[1]*triangles[i]->area())/longueurX/longueurY);
			stressMoyenne.z += ((stress[2]*triangles[i]->area())/longueurX/longueurY);
		}	

		
// 		std::cout<<"stress moyenne x  "<<stressMoyenne.x<<std::endl;
// 		std::cout<<"stress moyenne y  "<<stressMoyenne.y<<std::endl;
// 		std::cout<<"stress moyenne xy  "<<strainMoyenne.z<<std::endl;
// 		std::cout<<"strain moyenne x  "<<strainMoyenne.x<<std::endl;
// 		std::cout<<"strain moyenne y  "<<strainMoyenne.y<<std::endl;
// 		std::cout<<"strain moyenne xy  "<<strainMoyenne.z<<std::endl;
// 	
		double modulePoisson =(stressMoyenne.y*strainMoyenne.x-stressMoyenne.x*strainMoyenne.y)/(stressMoyenne.x*strainMoyenne.x-stressMoyenne.y*strainMoyenne.y);
	
		double moduleE = (1./strainMoyenne.x)*( stressMoyenne.x-stressMoyenne.y*modulePoisson);
		
		double moduleG = stressMoyenne.z/strainMoyenne.z;
		double moduleEEffectif = stressMoyenne.x/strainMoyenne.x;
// 		std::cout<<"module elastique "<<moduleE<<std::endl;
// 		std::cout<<"module de cisaillement "<<moduleG<<std::endl;
// 		std::cout<<"module de poisson "<<modulePoisson<<std::endl;

// 		std::cout<<"1 inclusion "<<inclusions[0]->getRadius()<<std::endl;
	

// 		if(a==1)
// 		{
// 		file <<std::endl<<std::endl<< "\"essai en mode\" , " << "contrainte" << std::endl ;
// 		file << "\"Type de granulo\" , " << "BOLOME_A" << std::endl ;
// 		file << "\"nombre point sur le périmètre\" , " << "128" << std::endl ;
// 		file << "\"e/c\" , " << "0.5" << std::endl ;
// 		file << "\"masse de ciment par m3 de béton\" , " << "350 kg/m3" << std::endl ;
// 		file << "\"numéro d'essai\" , " << nombreEssai << std::endl ;
// 		file << "\"masse initiale \" , " << masseInitiale<<" , g" << std::endl ;
// 		file << "\"densité granulat\" , " << densite<<" , g/mm3" << std::endl ;
// 		file << "\"rayon granulat max\" , " << rayonGranulatMax<<" , mm" << std::endl ;
// 		file << "\"module élastique du ciment\" , " << cimentE<<" , Pa" << std::endl ;	
// 		file << "\"coefficient de poisson du ciment\" , " << cimentNu<<" , -" << std::endl ;	
// 		file << "\"module élastique des granulats\" , " << granulatE<<" , Pa" << std::endl ;	
// 		file << "\"coefficient de poisson des granulats\" , " << granulatNu<<" , -" << std::endl ;	
// 		file << "\"longueur X\" , " << longueurX<<" , mm" << std::endl ;
// 		file << "\"longueur Y\" , " << longueurY<<" , mm" << std::endl ;
// 		file << "\"pourcentage de masse min\" , " << pourcentMasseMin<<" , %" << std::endl ;
// 		file << "\"nombre de granulats générés\" , " << nombreGranulatsGeneres << std::endl ;
// 		file << "\"nombre de granulats placés\" , " << nombreGranulatsPlaces << std::endl ;
// 		file << "\"rayon plus petit granulat\" , " << inclusionsPlacees[nombreGranulatsPlaces-1]->getRadius()<<" , mm"  << std::endl ;
// 		file << "\"pourcentage de l'aire occupé par granulats \" , " << pourcentAgregatsAireOccupe <<" , %"<< std::endl ;
// 		file << "\"nombre de répétition d'essai de placement\" , " << triesMax << std::endl ;
// 		file << "\"contrainte moyenne x\" , " << stressMoyenne.x<<" , Pa" << std::endl ;
// 		file << "\"contrainte moyenne y\" , " <<stressMoyenne.y <<" , Pa" << std::endl ;
// 		file << "\"contrainte moyenne xy\" , " << stressMoyenne.z<<" , Pa" << std::endl ;
// 		file << "\"déformation moyenne x\" , " << strainMoyenne.x<<" , -" << std::endl ;
// 		file << "\"déformation moyenne y\" , " << strainMoyenne.y<<" , -" << std::endl ;
// 		file << "\"déformation moyenne xy\" , " <<strainMoyenne.z <<" , -" << std::endl ;
// 		file << "\"module de poisson\" , " << modulePoisson << std::endl ;
// 		file << "\"module élastique réel\" , " << moduleE<<" , Pa" << std::endl ;
// 		file << "\"module de cisaillement\" , " << moduleG<<" , Pa" << std::endl ;
// 		file << "\"module élastique effectif\" , " << moduleEEffectif<<" , Pa" << std::endl ;
// 		}







		/**Affichage résultats généraux*/
		if(a==1)
		{
			file << "\"essai en mode\" , " << "\"Type de granulo\" , " << "\"nombre point sur le périmètre\" , " << "\"e/c\" , " << "\"masse de ciment par m3 de béton\" , " <<  "\"nombre d'essai\" , " <<  "\"masse initiale \" , " << "\"densité granulat\" , " <<  "\"rayon granulat max\" , " << "\"module élastique du ciment\" , " <<  "\"coefficient de poisson du ciment\" , " <<  "\"module élastique des granulats\" , " <<  "\"coefficient de poisson des granulats\" , " << "\"longueur X\" , " << "\"longueur Y\" , " <<  "\"pourcentage de masse min\" , " << "\"nombre de granulats générés\" , " << "\"nombre de granulats placés\" , " <<  "\"rayon plus petit granulat\" , " << "\"pourcentage de l'aire occupé par granulats \" , " <<  "\"nombre de répétition d'essai de placement\" , " << "\"contrainte moyenne x\" , " << "\"contrainte moyenne y\" , " <<"\"contrainte moyenne xy\" , " <<  "\"déformation moyenne x\" , " << "\"déformation moyenne y\" , " <<  "\"déformation moyenne xy\" , " << "\"module de poisson\" , " << "\"module élastique réel\" , " <<  "\"module de cisaillement\" , " <<  "\"module élastique effectif\" , " << std::endl ;
		}
		file << "deformation" << " , "<< "BOLOME_A" << " , "<<	"64" << " , "<< "0.5" << " , " << "350 kg/m3" << " , " << nombreEssai << " , " << masseInitiale<<" , " << densite<<" , " << rayonGranulatMax<<" , " << cimentE<<" , " << cimentNu<<" , " << granulatE<<" , " << granulatNu<<" , " << longueurX<<" , " << longueurY<<" , " << pourcentMasseMin<<" , " << nombreGranulatsGeneres << " , " << nombreGranulatsPlaces << " , " << inclusionsPlacees[nombreGranulatsPlaces-1]->getRadius()<<" , " << pourcentAgregatsAireOccupe <<" , " << triesMax << " , " << stressMoyenne.x<<" , " <<stressMoyenne.y <<" , " << stressMoyenne.z<<" , " << strainMoyenne.x<<" , " << strainMoyenne.y<<" , " <<strainMoyenne.z <<" , " << modulePoisson << " , " << moduleE<<" , " << moduleG<<" , " << moduleEEffectif<< std::endl ;



	 	/**Affichage rayon et coordonnée x (horizontal)*/
// 		file<<"r , ";
// 		for(size_t i=0; i<nombreGranulatsPlaces; i++)
// 		{ 
// 			if(inclusions[i]->getCenter()->y>25 && inclusions[i]->getCenter()->y<75)// 15 35 pour 50x50
// 			{
// 	// 			std::cout<<inclusions[i]->getRadius()<<"   "<<inclusions[i]->getCenter()->x<<std::endl;
// 	// 			file << inclusions[i]->getRadius() << inclusions[i]->getCenter()->x<< std::endl ;
// 				file << inclusions[i]->getRadius() << " , " ;
// 			}
// 		}
// 		file<<std::endl<<"x , ";
// 		for(size_t i=0; i<nombreGranulatsPlaces; i++)
// 		{ 
// 			if(inclusions[i]->getCenter()->y>25 && inclusions[i]->getCenter()->y<75)
// 			{
// 	// 			std::cout<<inclusions[i]->getRadius()<<"   "<<inclusions[i]->getCenter()->x<<std::endl;
// 				file <<  inclusions[i]->getCenter()->x<< " , " ;
// 			}
// 		}
// 		file<<std::endl;

		/**Affichage rayon et coordonnée x (vertical)*/
// 		file<<"r , x"<<std::endl;
// 		for(size_t i=0; i<nombreGranulatsPlaces; i++)
// 		{ 
// 			if(inclusions[i]->getCenter()->y>25 && inclusions[i]->getCenter()->y<75)// 15 35 pour 50x50
// 			{
// // 				std::cout<<inclusions[i]->getRadius()<<"   "<<inclusions[i]->getCenter()->x<<std::endl;
// 				file << inclusions[i]->getRadius()<< " , " << inclusions[i]->getCenter()->x<< std::endl ;
// 				
// 			}
// 		}
		



		/**Affichage contrainte et déformation par triangle*/

// 		file << std::endl<<std::endl<<"\"déformation par triangle\"  "  << std::endl ;
// 		for(size_t i=0; i<triangles.size();i++)
// 		{
// 			file << deformation[i*3]<<" , " << deformation[i*3+1]<<" , "<<deformation[i*3+2]<< std::endl ;
// 		}
// 		file << std::endl<<std::endl<<"\"contrainte par triangle\"  "  << std::endl ;
// 		for(size_t i=0; i<triangles.size();i++)
// 		{
// 			file << contrainte[i*3]<<" , " << contrainte[i*3+1]<<" , "<<contrainte[i*3+2]<< std::endl ;
// 		}



		/**Affichage rayon et %masse réelle*/
// 		file<<"r , %masse réelle"<<std::endl;
// 		double pourcentMasseGranulatsPlaces=1;
// 		double masseTotaleGranulatsPlaces=0;
// 		for(size_t i=0; i<nombreGranulatsPlaces;i++)
// 		{
// 			masseTotaleGranulatsPlaces += 4.*3.14*densite*(inclusions[i]->getRadius()*inclusions[i]->getRadius()*inclusions[i]->getRadius())/3.;
// 		}
// 		double masse=0;
// 		for(size_t i=0; i<nombreGranulatsPlaces;i++)
// 		{
// 			file<<inclusions[i]->getRadius()<<" , "<<pourcentMasseGranulatsPlaces<<std::endl;
// 			masse += 4.*3.14*densite*(inclusions[i]->getRadius()*inclusions[i]->getRadius()*inclusions[i]->getRadius())/3.;
// 			pourcentMasseGranulatsPlaces = (masseTotaleGranulatsPlaces-masse)/masseTotaleGranulatsPlaces; 
// 			
// 		}


		/**Affichage rayon et %masse initiale*/
		file<<"r , %masse initiale"<<std::endl;
		double pourcentMasseGranulatsPlaces=1;
		double masseTotaleGranulatsPlaces=0;

		double masse=0;
		for(size_t i=0; i<nombreGranulatsPlaces;i++)
		{
			file<<inclusions[i]->getRadius()<<" , "<<pourcentMasseGranulatsPlaces<<std::endl;
			masse += 4.*3.14*densite*(inclusions[i]->getRadius()*inclusions[i]->getRadius()*inclusions[i]->getRadius())/3.;
			pourcentMasseGranulatsPlaces = (masseInitiale-masse)/masseInitiale; 
			
		}


		/**Contrainte x et %V*/
// 		file<<"Contrainte x, %V"<<std::endl;
// 		std::vector <double> intervalle;
// 		double min=contrainte[0];
// 		double max=contrainte[0];
// 		for(size_t i=0;i<contrainte.size()/3; i++)
// 		{
// 			for (size_t j=0; j<i; j++)
// 			{
// 				if(max<contrainte[j*3])
// 				{
// 					max=contrainte[j*3];
// 				}
// 				if(min>contrainte[j*3])
// 				{
// 					min = contrainte[j*3];
// 				}
// 			}
// 		}		
// 		double inter = (max-min)/100.;
// 		for(size_t i=0; i<100;i++)
// 		{
// 			intervalle.push_back(i*inter+min);
// 		}	
// 		std::vector <double> volume;
// 		for(size_t i=0; i<intervalle.size()-1;i++)
// 		{
// 			volume.push_back(0.);
// 		}
// 		for(size_t i=0; i<triangles.size();i++)
// 		{
// 			for(size_t j=0; j<triangles.size();j++)
// 			{
// 				if(contrainte[3*i]>=intervalle[j] && contrainte[3*i]<intervalle[j+1])
// 				{
// 					
// 					volume[j]+= triangles[i]->area()/(longueurX*longueurY);
// 					break;
// 				}
// 			}
// 		}
// 		for(size_t i=0; i<100; i++)
// 		
// 		{
// 			file<<intervalle[i]<<" , "<<volume[i]<<std::endl;
// 		}

		/**Déformation x et %V*/
// 		file<<"déformation x, %V"<<std::endl;
// 		std::vector <double> defintervalle;
// 		double defmin=deformation[0];
// 		double defmax=deformation[0];
// 		for(size_t i=0;i<deformation.size()/3; i++)
// 		{
// 			for (size_t j=0; j<i; j++)
// 			{
// 				if(defmax<deformation[j*3])
// 				{
// 					defmax=deformation[j*3];
// 				}
// 				if(defmin>deformation[j*3])
// 				{
// 					defmin = deformation[j*3];
// 				}
// 			}
// 		}		
// 		double definter = (defmax-defmin)/100.;
// 		for(size_t i=0; i<100;i++)
// 		{
// 			defintervalle.push_back(i*definter+defmin);
// 		}	
// 		std::vector <double> defvolume;
// 		for(size_t i=0; i<defintervalle.size()-1;i++)
// 		{
// 			defvolume.push_back(0.);
// 		}
// 		for(size_t i=0; i<triangles.size();i++)
// 		{
// 			for(size_t j=0; j<triangles.size();j++)
// 			{
// 				if(deformation[3*i]>=defintervalle[j] && deformation[3*i]<defintervalle[j+1])
// 				{
// 					
// 					defvolume[j]+= triangles[i]->area()/(longueurX*longueurY);
// 					break;
// 				}
// 			}
// 		}
// 		for(size_t i=0; i<100; i++)
// 		
// 		{
// 			file<<defintervalle[i]<<" , "<<defvolume[i]<<std::endl;
// 		}*/

		
		if(nombreEssai==(a))
		{
			file.close() ;
			glutInitWindowSize(900, 900) ;
			glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA);
			glutInit(&argc, argv) ;	
			glutInitWindowPosition(0,0);
		
			// Open a window, name it "Hello World"
			if (glutCreateWindow("Hello World") == GL_FALSE) {
				return 0;
			}
		
			// Set the clear color to black
			glClearColor(0.0, 0.0, 0.0, 0.0);
		
			// Assign reshape() to be the function called whenever 
			// a reshape event occurs
			glutReshapeFunc(reshape);
		
			// Assign draw() to be the function called whenever a display
			// event occurs, generally after a resize or expose event
			glutDisplayFunc(draw);
		
			// Pass program control to tk's event handling code
			// In other words, loop forever
			glutMainLoop();
		}
	inclusions.clear();
	inclusionsPlacees.clear();
	triangles.clear();
	}

	return 0 ;
}

