
// Author:  Jérôme Krebs <jerome.krebs@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//

#include "placement.h"
#include "features/inclusion.h"


using namespace Mu ;


double chiffreAleatoire(double longueur) // fonction qui retourne une valeur aléatoire
{
	double chiffreAleatoire = (double) rand()/RAND_MAX*longueur;
	return chiffreAleatoire;
}



bool bord(double r, double longueurX, double longueurY, double x, double y)//fonction du problème de bord
{
	if(r>x ||r>(longueurX-x) ||r>y||r>(longueurY-y))
	{			
		return true;	
	}
	return false;
}


// double espaceMinGranulat(int j,std::vector <Inclusion *> inclusions)
// {
// 	double rayonGranulat;
// 	if(inclusions[j]->getRadius()>4.0)//espace entre 2 granulats
// 	{	
// 		rayonGranulat = inclusions[j]->getRadius()+0.1;
// 	}
// 	else rayonGranulat =inclusions[j]->getRadius()+0.1*inclusions[j]->getRadius();
// 	return rayonGranulat;
// }



// Point calculNewCoordGranulat(Point coordNewGranulat, Point coordOldGranulat,std::vector <Inclusion *> inclusions,int i, int j)
// {
// 	
// 	Point vecteur2Point = coordNewGranulat-coordOldGranulat;
// 	Point vecteur2PointUnitaire = vecteur2Point/vecteur2Point.norm();
// 	double distance = inclusions[i]->getRadius()+inclusions[j]->getRadius()+0;// 0.05=1 rayon du plus petit granulat
// 	coordNewGranulat = coordOldGranulat + (vecteur2PointUnitaire*distance);
// 	j=-1;
// 	return coordNewGranulat;
// }


std::vector<Inclusion *> placement(double longueurX, double longueurY, std::vector<Inclusion *> inclusions, int *nombreGranulatsPlaces, int triesMax)
{
	bool suite =true;
	bool sortie;
	int tries = 0 ;
	double x;
	double y;

	for(int i=0; i<inclusions.size()&& tries< triesMax; i++) //boucle qui calcul les coordonn�s al�toires des granulats  
		{
			
			*nombreGranulatsPlaces = i;
			sortie = false; 
			tries++ ;
			x=chiffreAleatoire(longueurX);

			y=chiffreAleatoire(longueurY);
			//std::cout<< "\r"<< i<<"/" << inclusions.size() << "  "<<inclusions[i] << "    "<<x<<"     "  << y<< "         " <<std::flush; 	

			Point coordNewGranulat(x,y);

			if (bord(inclusions[i]->getRadius(), longueurX, longueurY, coordNewGranulat.x,coordNewGranulat.y))
			{
				i--;
				//std::cout<<"granulat trop pres du bord"<<std::endl;
// 				tries++ ;
				suite=false;
			}

			if (suite)//probl�e du superposage des granulats
			{
				
				Circle c1(inclusions[i]->getRadius(),coordNewGranulat.x,coordNewGranulat.y);
				for(int j=0; j<i; ++j)
				{
// 					tries++ ;
					Point coordOldGranulat(inclusions[j]->getCenter().x,inclusions[j]->getCenter().y);	
					double rayonGranulat;

					//rayonGranulat=espaceMinGranulat(j,inclusions);
					if(inclusions[j]->getRadius()>4.0)//espace entre 2 granulats
					{	
						rayonGranulat = inclusions[j]->getRadius()+0.04;
					}
					else rayonGranulat =inclusions[j]->getRadius()+0.01*inclusions[j]->getRadius();

					Circle c2(rayonGranulat,coordOldGranulat.x,coordOldGranulat.y);
					bool inter = c1.intersects(&c2);
					
					/** placement normal */
					if(inter==true) //s'il y a intersection et que les nouvelles coordonnées ont déjà été calculées une fois
					{
						i--;
						tries++ ;
						sortie =false;
						suite=false;
						break;
					}					



					/** placement amélioré */
// 					if(inter==true && sortie ==true) //s'il y a intersection et que les nouvelles coordonnées ont déjà été calculées une fois			
// 					
// 					{
// 						i--;
// // 						tries++ ;
// 						sortie =false;
// 						suite=false;
// 						break;
// 					}
// 					
// 					if(inter==true) //s'il y a intersection, on cherche un nouveau centre qui est r1+r2+0.1 
// 					{
// 						//coordNewGranulat= calculNewCoordGranulat(coordNewGranulat, coordOldGranulat, inclusions,i,j);
// 						tries++;
// 						Point vecteur2Point = coordNewGranulat-coordOldGranulat;
// 						Point vecteur2PointUnitaire = vecteur2Point/vecteur2Point.norm();
// 						double distance = inclusions[i]->getRadius()+inclusions[j]->getRadius()+0;// 0.05=1 rayon du plus petit granulat
// 						coordNewGranulat = coordOldGranulat + (vecteur2PointUnitaire*distance);
// 						j=-1;
// 						
// 						if (bord(inclusions[i]->getRadius(), longueurX, longueurY, coordNewGranulat.x,coordNewGranulat.y))
// 						{
// 							//std::cout<<"granulat trop pres du bord"<<std::endl;
// 							i--;
// // 							tries++ ;
// 							suite=false;
// 							break;
// 						}
// 						sortie = true;	
// 					}
				}
			}
		if(suite)
			{	
				//std::cout<<coordNewGranulat.x<<std::endl;
				tries = 0 ;
				inclusions[i]->getCenter().x = coordNewGranulat.x;
				inclusions[i]->getCenter().y = coordNewGranulat.y;	
			}
		suite =true;

		}
	return inclusions;
		
}

