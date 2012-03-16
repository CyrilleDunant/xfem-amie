
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2009-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "main.h"
#include "../physics/homogenization/composite.h"
#include "../physics/homogenization/phase.h"
#include "../physics/materials/gel_behaviour.h"
#include "../physics/materials/aggregate_behaviour.h"
#include "../physics/materials/paste_behaviour.h"

#include <fstream>

using namespace Mu ;

int main(int argc, char *argv[])
{
	std::fstream writer ;
	writer.open("elastic.csv", std::ios::out) ;
	for(int i = 0 ; i < 1000 ; i++)
	{
		double f_gel = ((double) i)*0.001 ;

		Phase gel(new GelBehaviour(), f_gel) ;
		Phase aggregate(new ElasticOnlyAggregateBehaviour(), (1.-f_gel)) ;
		MoriTanakaMatrixInclusionComposite lowerAggregate(aggregate, gel) ;
		lowerAggregate.apply() ;

		double f_paste = 0.36 ;
		lowerAggregate.volume = 1.-f_paste ;

		Phase cement(new ElasticOnlyPasteBehaviour(), 0.36) ;
		MoriTanakaMatrixInclusionComposite lowerConcrete(cement, lowerAggregate) ;
		lowerConcrete.apply() ;
		Matrix S = lowerConcrete.C ;
		Composite::invertTensor(S) ;
		Vector lowerAlpha = S*lowerConcrete.beta ;

		InverseMoriTanakaMatrixInclusionComposite upperAggregate(aggregate, gel) ;
		upperAggregate.apply() ;

		upperAggregate.volume = 1.-f_paste ;

		InverseMoriTanakaMatrixInclusionComposite upperConcrete(cement, lowerAggregate) ;
		upperConcrete.apply() ;
		S = upperConcrete.C ;
		Composite::invertTensor(S) ;
		Vector upperAlpha = S*upperConcrete.beta ;
	
		writer << f_gel << "," << lowerConcrete.C[0][0] << "," << lowerAlpha[0] << "," << upperConcrete.C[0][0] << "," << upperAlpha[0]<< std::endl ;
	}
	
	return 0 ;
}
