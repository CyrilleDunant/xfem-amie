
#include "optimizer.h"
#include "itoa.h"

namespace Mu
{

double normalRand(double mean, double var)
{
	double m = (double)rand()/(double)(RAND_MAX);
	double c = (double)rand()/(double)(RAND_MAX);
	return (sqrt(-2.*log(m))*cos(2.*M_PI*c))*var+mean ;

}
	
GeneticAlgorithmOptimizer::GeneticAlgorithmOptimizer(std::vector<double> vars, std::vector<std::pair<double, double> >  bounds, std::vector<double *> lowLevelVars, std::vector<std::pair<double, double> >  lowLevelbounds, const Function & objectiveFunction) : vars(vars), bounds(bounds), lowLevelVars(lowLevelVars), lowLevelbounds(lowLevelbounds), objectiveFunction(objectiveFunction), lowLevelFunction(NULL)
{

}
GeneticAlgorithmOptimizer::GeneticAlgorithmOptimizer(const std::vector<double> & vars, const std::vector<std::pair<double, double> >  & bounds, const Function & objectiveFunction): vars(vars), bounds(bounds), objectiveFunction(objectiveFunction), lowLevelFunction(NULL)
{
	
}

GeneticAlgorithmOptimizer::GeneticAlgorithmOptimizer(std::vector<double *> lowLevelVars,  std::vector<std::pair<double, double> >  lowLevelbounds, double (*lowLevelFunction)()) : lowLevelVars(lowLevelVars), lowLevelbounds(lowLevelbounds), lowLevelFunction(lowLevelFunction)
{
}

GeneticAlgorithmOptimizer::GeneticAlgorithmOptimizer(const std::vector<std::pair<std::string, double> > & nvars, const std::vector<std::pair<double, double> >  & nbounds, const Function & objectiveFunction) : namedVars(nvars), namedVarsBounds(nbounds), objectiveFunction(objectiveFunction), lowLevelFunction(NULL)
{
}

double GeneticAlgorithmOptimizer::applyFunction(const std::vector<double> & vals) const
{
	// x, y, z, t, u, v, w
	
	switch (vals.size())
	{
		case 0 :
			return VirtualMachine().eval(objectiveFunction, namedVars) ;
		case 1 :
			return VirtualMachine().eval(objectiveFunction, namedVars, vals[0]) ;
		case 2 :
			return VirtualMachine().eval(objectiveFunction, namedVars, vals[0], vals[1]) ;
		case 3 :
			return VirtualMachine().eval(objectiveFunction, namedVars, vals[0], vals[1], vals[2]) ;
		case 4 :
			return VirtualMachine().eval(objectiveFunction, namedVars, vals[0], vals[1], vals[2], vals[3]) ;
		case 5 :
			return VirtualMachine().eval(objectiveFunction, namedVars, vals[0], vals[1], vals[2], vals[3], vals[4]) ;
		case 6 :
			return VirtualMachine().eval(objectiveFunction, namedVars, vals[0], vals[1], vals[2], vals[3], vals[4], vals[5]) ;
		case 7 :
			return VirtualMachine().eval(objectiveFunction, namedVars, vals[0], vals[1], vals[2], vals[3], vals[4], vals[5], vals[6]) ;
	}
	return 0. ;
}

double GeneticAlgorithmOptimizer::optimize(double eps, int Maxit, int population, double elitism, double factor)
{
	if(lowLevelFunction)
		return lowLevelOptimize(eps, Maxit, population, elitism, factor) ;
	
	int it = 0 ;
	int maxit = Maxit ;
	if(maxit < 0)
		maxit = 100 ;
	
	if(population < 0)
		population = 10 ;
	
	std::vector< std::vector<double> > individuals ;
	std::vector< std::vector<double> > nindividuals ;
	std::vector< std::vector<double> > llindividuals ;
	for(size_t i = 0 ; i < population ; i++)
	{
		std::vector<double> newindividual ;
		for(size_t j = 0 ;  j < vars.size() ; j++)
		{
			newindividual.push_back((bounds[j].second-bounds[j].first)*random()/RAND_MAX + bounds[j].first);
		}
		individuals.push_back(newindividual);
	}
	
	for(size_t i = 0 ; i < population ; i++)
	{
		std::vector<double> newindividual ;
		for(size_t j = 0 ;  j < namedVars.size() ; j++)
		{
			newindividual.push_back((namedVarsBounds[j].second-namedVarsBounds[j].first)*random()/RAND_MAX + namedVarsBounds[j].first);
		}
		nindividuals.push_back(newindividual);
	}
	
	for(size_t i = 0 ; i < population ; i++)
	{
		std::vector<double> newllindividual ;
		for(size_t j = 0 ;  j < lowLevelVars.size() ; j++)
		{
			newllindividual.push_back((lowLevelbounds[j].second-lowLevelbounds[j].first)*random()/RAND_MAX + lowLevelbounds[j].first);
		}
		llindividuals.push_back(newllindividual);
	}
	
	double err = eps*100 ;
	std::map<double, std::vector<std::vector<double> > > sorted ;
	while (err > eps && it < maxit )
	{
		//evaluate
		for(size_t i = 0 ; i < population ; i++)
		{
			for(size_t j = 0 ;  j < lowLevelVars.size() ; j++)
				*lowLevelVars[j] = llindividuals[i][j] ;
			for(size_t j = 0 ;  j < namedVars.size() ; j++)
				namedVars[j].second = nindividuals[i][j] ;
			
			std::vector<std::vector<double> >  curr;
			curr.push_back(individuals[i]);
			curr.push_back(nindividuals[i]);
			curr.push_back(llindividuals[i]);
			sorted[applyFunction(individuals[i])] = curr ;
		}
		
		std::vector< std::vector<double> > newindividuals ;
		std::vector< std::vector<double> > newnindividuals ;
		std::vector< std::vector<double> > newllindividuals ;
				
		//reproduce - elite individuals are kept
		for(size_t i = 0 ; i < std::min((int)(population*elitism), (int)sorted.size()) ; i++)
		{
			std::map<double, std::vector<std::vector<double> > >::iterator iter = sorted.begin() ;
			for(size_t j = 0 ; j < i ; j++)
				iter++ ;
			newindividuals.push_back(iter->second[0]);
			newnindividuals.push_back(iter->second[1]);
			newllindividuals.push_back(iter->second[2]);
		}
		
		//reproduce - the rest fills the available slots and mutates
		int i = 0 ;
		while( newindividuals.size() < individuals.size())
		{
			for(size_t n = 0 ; n < 4 && newindividuals.size() <= individuals.size(); n++)
			{
				std::map<double, std::vector<std::vector<double> > >::iterator iter = sorted.begin() ;
				for(size_t j = 0 ; j < i ; j++)
					iter++ ;
				newindividuals.push_back(iter->second[0]);
				newnindividuals.push_back(iter->second[1]);
				newllindividuals.push_back(iter->second[2]);
				double factor = 0.85;
				for(size_t j = 0 ; j < newindividuals.back().size() ; j++)
				{
					double sigma = err*(namedVarsBounds[j].second-namedVarsBounds[j].first)/population*factor ;
					sigma = std::min(sigma, (bounds[j].second-bounds[j].first)/2) ;
					newindividuals.back()[j] = normalRand(newindividuals.back()[j], sigma) ; //newindividuals.back()[j]*(1.-err*population*.5) + err*population*.5*(bounds[j].second-bounds[j].first)*random()/RAND_MAX ;
					newindividuals.back()[j] = std::max(newindividuals.back()[j], bounds[j].first) ;
					newindividuals.back()[j] = std::min(newindividuals.back()[j], bounds[j].second) ;
				}
				for(size_t j = 0 ; j < newnindividuals.back().size() ; j++)
				{
					double sigma = err*(namedVarsBounds[j].second-namedVarsBounds[j].first)/population*factor ;
					sigma = std::min(sigma, (namedVarsBounds[j].second-namedVarsBounds[j].first)/2) ;
					newnindividuals.back()[j] = normalRand(newnindividuals.back()[j], sigma) ;//newnindividuals.back()[j]*(1.-err*population*.5) + err*population*.5*(namedVarsBounds[j].second-namedVarsBounds[j].first)*random()/RAND_MAX ;
					newnindividuals.back()[j] = std::max(newnindividuals.back()[j], namedVarsBounds[j].first) ;
					newnindividuals.back()[j] = std::min(newnindividuals.back()[j], namedVarsBounds[j].second) ;
				}
				for(size_t j = 0 ; j < newllindividuals.back().size() ; j++)
				{
					double sigma = err*(lowLevelbounds[j].second-lowLevelbounds[j].first)/population*factor ;
					sigma = std::min(sigma, (lowLevelbounds[j].second-lowLevelbounds[j].first)/2) ;
					newllindividuals.back()[j] = normalRand(newllindividuals.back()[j], sigma) ;//newllindividuals.back()[j]*(1.-err*population*.5) + err*population*.5*(lowLevelbounds[j].second-lowLevelbounds[j].first)*random()/RAND_MAX ;
					newllindividuals.back()[j] = std::max(newllindividuals.back()[j], lowLevelbounds[j].first) ;
					newllindividuals.back()[j] = std::min(newllindividuals.back()[j], lowLevelbounds[j].second) ;
				}
			}
			
			i++ ;
			i = std::min(i, (int)sorted.size()-1) ;
		}
// 		std::cout << sorted.begin()->first << std::endl ; 
		individuals = newindividuals ;
		nindividuals = newnindividuals ;
		llindividuals = newllindividuals ;
		
		//iterate
		it++ ;
		err = sorted.begin()->first ;
	}
	std::cout << it << std::endl;
	vars = sorted.begin()->second[0] ;
	for(size_t i = 0 ; i < namedVars.size() ; i++)
		namedVars[i].second = sorted.begin()->second[1][i] ;
	for(size_t i = 0 ; i < lowLevelVars.size() ; i++)
		*lowLevelVars[i] = sorted.begin()->second[2][i] ;
	
	return err ;
}

double GeneticAlgorithmOptimizer::lowLevelOptimize(double eps, int Maxit, int population, double elitism, double factor)
{
	int it = 0 ;
	int maxit = Maxit ;
	if(maxit < 0)
		maxit = 100 ;
	
	if(population < 0)
		population = 10 ;
	
	std::vector< std::vector<double> > llindividuals ;

	for(size_t i = 0 ; i < population ; i++)
	{
		std::vector<double> newllindividual ;
		for(size_t j = 0 ;  j < lowLevelVars.size() ; j++)
		{
			newllindividual.push_back((lowLevelbounds[j].second-lowLevelbounds[j].first)*random()/RAND_MAX + lowLevelbounds[j].first);
		}
		llindividuals.push_back(newllindividual);
	}
	
	double err = eps*100 ;
	std::map<double, std::vector<double> > sorted ;
	while (err > eps && it < maxit )
	{
		//evaluate
		for(size_t i = 0 ; i < population ; i++)
		{
			for(size_t j = 0 ;  j < lowLevelVars.size() ; j++)
				*lowLevelVars[j] = llindividuals[i][j] ;
			double llf = lowLevelFunction() ;
			sorted[llf] = llindividuals[i] ;
		}
		
		std::vector< std::vector<double> > newllindividuals ;
				
		//reproduce - elite individuals are kept
		for(size_t i = 0 ; i < std::min((int)(population*elitism), (int)sorted.size()) ; i++)
		{
			std::map<double, std::vector<double> >::iterator iter = sorted.begin() ;
			for(size_t j = 0 ; j < i ; j++)
				iter++ ;
			newllindividuals.push_back(iter->second);
		}
		
		//reproduce - the rest fills the available slots and mutates
		int i = 0 ;
		while( newllindividuals.size() < llindividuals.size())
		{
			for(size_t n = 0 ; n < 4 && newllindividuals.size() <= llindividuals.size(); n++)
			{
				std::map<double, std::vector<double> >::iterator iter = sorted.begin() ;
				for(size_t j = 0 ; j < i ; j++)
					iter++ ;
				newllindividuals.push_back(iter->second);
				double factor = 0.85;

				for(size_t j = 0 ; j < newllindividuals.back().size() ; j++)
				{
					double sigma = err*(lowLevelbounds[j].second-lowLevelbounds[j].first)/population*factor ;
					sigma = std::min(sigma, (lowLevelbounds[j].second-lowLevelbounds[j].first)/2) ;
					newllindividuals.back()[j] = normalRand(newllindividuals.back()[j], sigma) ;//newllindividuals.back()[j]*(1.-err*population*.5) + err*population*.5*(lowLevelbounds[j].second-lowLevelbounds[j].first)*random()/RAND_MAX ;
					newllindividuals.back()[j] = std::max(newllindividuals.back()[j], lowLevelbounds[j].first) ;
					newllindividuals.back()[j] = std::min(newllindividuals.back()[j], lowLevelbounds[j].second) ;
				}
			}
			
			i++ ;
			i = std::min(i, (int)sorted.size()-1) ;
		}
// 		std::cout << sorted.begin()->first << std::endl ; 

		llindividuals = newllindividuals ;
		for(size_t j = 0 ;  j < lowLevelVars.size() ; j++)
			*lowLevelVars[j] = sorted.begin()->second[j] ;
		//iterate
		it++ ;
		err = sorted.begin()->first ;
		std::cout << getValues()[0].first << " " << getValues()[0].second <<std::endl ;
		std::cout << getValues()[1].first << " " << getValues()[1].second <<std::endl ;
		std::cout << getValues()[2].first << " " << getValues()[2].second <<std::endl ;
		std::cout << getValues()[3].first << " " << getValues()[3].second <<std::endl ;
		std::cout << err << std::endl ;
	}
	std::cout << it << std::endl;

	for(size_t i = 0 ; i < lowLevelVars.size() ; i++)
		*lowLevelVars[i] = sorted.begin()->second[i] ;
	
	return err ;
}

std::vector<std::pair<std::string, double> > GeneticAlgorithmOptimizer::getValues() const
{
	std::vector<std::pair<std::string, double> > ret ;
	
	for(size_t i = 0 ; i < vars.size() ; i++)
	{
		if(i == 0)
			ret.push_back(std::make_pair("x", vars[i]));
		if(i == 1)
			ret.push_back(std::make_pair("y", vars[i]));
		if(i == 2)
			ret.push_back(std::make_pair("z", vars[i]));
		if(i == 3)
			ret.push_back(std::make_pair("t", vars[i]));
		if(i == 4)
			ret.push_back(std::make_pair("u", vars[i]));
		if(i == 5)
			ret.push_back(std::make_pair("v", vars[i]));
		if(i == 6)
			ret.push_back(std::make_pair("w", vars[i]));
	}
	for(size_t i = 0 ; i < namedVars.size() ; i++)
	{
		ret.push_back(namedVars[i]);
	}
	
	for(size_t i = 0 ; i < lowLevelVars.size() ; i++)
	{
		std::string name("lowLevel") ;
		name += itoa(i) ;
		ret.push_back(std::make_pair(name, *lowLevelVars[i]));
	}
	
	return ret ;
}

LeastSquaresApproximation::LeastSquaresApproximation(const Vector & measures, const Matrix & linearModel) : measures(measures), linearModel(linearModel), linearModelChanged(true), X0t(linearModel.numCols(), linearModel.numRows()), X0tX0(linearModel.numRows(), linearModel.numRows()), parameters(0., linearModel.numRows())
{
	X0t   = linearModel.transpose() ;
	X0tX0 = linearModel*X0t ;
	
}

void LeastSquaresApproximation::setParameterValue(int i, double v)
{
	for(size_t j = 0 ; j < linearModel.numCols() ; j++)
	{
		if(i == j)
			linearModel[i][j] = 1 ;
		else
		{
			linearModel[i][j] = 0 ;
		}
	}
	
	linearModelChanged = true ;
}

void LeastSquaresApproximation::setLinearModel(int i, int j, double v)
{
	linearModel[i][j] = v ;
	linearModelChanged = true ;
}

double LeastSquaresApproximation::optimize()
{
	if(linearModelChanged)
	{
		X0t = linearModel.transpose() ;
		X0tX0 = linearModel*X0t ;
	}
	
	Vector lb0 = linearModel*measures ;
	parameters = solveSystem(X0tX0, lb0, parameters) ;

	return std::abs((Vector)(X0tX0*parameters)-lb0).max() ;
	linearModelChanged = false ;
}

} ;
