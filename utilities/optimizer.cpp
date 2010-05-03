#include "optimizer.h"

namespace Mu
{

GeneticAlgorithmOptimizer::GeneticAlgorithmOptimizer(std::vector<double> vars, std::vector<std::pair<double, double> >  bounds, std::vector<double *> lowLevelVars,  
							  std::vector<std::pair<double, double> >  lowLevelbounds, const Function & objectiveFunction) : vars(vars), bounds(bounds), lowLevelVars(lowLevelVars), lowLevelbounds(lowLevelbounds), objectiveFunction(objectiveFunction)
{

}
GeneticAlgorithmOptimizer::GeneticAlgorithmOptimizer(const std::vector<double> & vars, const std::vector<std::pair<double, double> >  & bounds, const Function & objectiveFunction): vars(vars), bounds(bounds), objectiveFunction(objectiveFunction)
{
	
}

double GeneticAlgorithmOptimizer::applyFunction(const std::vector<double> & vals) const
{
	// x, y, z, t, u, v, w
	

	switch (vals.size())
	{
		case 1 :
			return VirtualMachine().eval(objectiveFunction, vals[0]) ;
		case 2 :
			return VirtualMachine().eval(objectiveFunction, vals[0], vals[1]) ;
		case 3 :
			return VirtualMachine().eval(objectiveFunction, vals[0], vals[1], vals[2]) ;
		case 4 :
			return VirtualMachine().eval(objectiveFunction, vals[0], vals[1], vals[2], vals[3]) ;
		case 5 :
			return VirtualMachine().eval(objectiveFunction, vals[0], vals[1], vals[2], vals[3], vals[4]) ;
		case 6 :
			return VirtualMachine().eval(objectiveFunction, vals[0], vals[1], vals[2], vals[3], vals[4], vals[5]) ;
		case 7 :
			return VirtualMachine().eval(objectiveFunction, vals[0], vals[1], vals[2], vals[3], vals[4], vals[5], vals[6]) ;
	}

	return 0. ;
}

double GeneticAlgorithmOptimizer::optimize(double eps, int Maxit, int population, double elitism)
{
	int it = 0 ;
	int maxit = Maxit ;
	if(maxit < 0)
		maxit = 100 ;
	
	if(population < 0)
		population = 10 ;
	
	std::vector< std::vector<double> > individuals ;
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
		std::vector<double> newllindividual ;
		for(size_t j = 0 ;  j < lowLevelVars.size() ; j++)
		{
			newllindividual.push_back((lowLevelbounds[j].second-lowLevelbounds[j].first)*random()/RAND_MAX + lowLevelbounds[j].first);
		}
		llindividuals.push_back(newllindividual);
	}
	
	double err = eps*100 ;
	std::map<double, std::pair<std::vector<double>, std::vector<double> > > sorted ;
	while (err > eps && it < maxit )
	{
		//evaluate
		
		for(size_t i = 0 ; i < population ; i++)
		{
			for(size_t j = 0 ;  j < lowLevelVars.size() ; j++)
				*lowLevelVars[j] = llindividuals[i][j] ;
			
			sorted[applyFunction(individuals[i])] = std::make_pair(individuals[i], llindividuals[i]) ;
		}
		
		std::vector< std::vector<double> > newindividuals ;
		std::vector< std::vector<double> > newllindividuals ;
				
		//reproduce - elite individuals are kept
		for(size_t i = 0 ; i < std::min((int)(population*elitism), (int)sorted.size()) ; i++)
		{
			std::map<double, std::pair<std::vector<double>, std::vector<double> > >::iterator iter = sorted.begin() ;
			for(size_t j = 0 ; j < i ; j++)
				iter++ ;
			newindividuals.push_back(iter->second.first);
			newllindividuals.push_back(iter->second.second);
		}
		
		//reproduce - the rest fills the available slots and mutates
		int i = 0 ;
		while( newindividuals.size() < individuals.size())
		{
			for(size_t n = 0 ; n < 4 && newindividuals.size() <= individuals.size(); n++)
			{
				std::map<double, std::pair<std::vector<double>, std::vector<double> > >::iterator iter = sorted.begin() ;
				for(size_t j = 0 ; j < i ; j++)
					iter++ ;
				newindividuals.push_back(iter->second.first);
				newllindividuals.push_back(iter->second.second);
				for(size_t j = 0 ; j < newindividuals.back().size() ; j++)
				{
					newindividuals.back()[j] = newindividuals.back()[j]*(1.-err*population) + err*population*(bounds[j].second-bounds[j].first)*random()/RAND_MAX ;
					newindividuals.back()[j] = std::max(newindividuals.back()[j], bounds[j].first) ;
					newindividuals.back()[j] = std::min(newindividuals.back()[j], bounds[j].second) ;
				}
				for(size_t j = 0 ; j < newllindividuals.back().size() ; j++)
				{
					newllindividuals.back()[j] = newllindividuals.back()[j]*(1.-err*population) + err*population*(lowLevelbounds[j].second-lowLevelbounds[j].first)*random()/RAND_MAX ;
					newllindividuals.back()[j] = std::max(newllindividuals.back()[j], lowLevelbounds[j].first) ;
					newllindividuals.back()[j] = std::min(newllindividuals.back()[j], lowLevelbounds[j].second) ;
				}
			}
			
			i++ ;
			i = std::min(i, (int)sorted.size()-1) ;
		}
		std::cout << sorted.begin()->first << std::endl ; 
		individuals = newindividuals ;
		llindividuals = newllindividuals ;
		
		//iterate
		it++ ;
		err = sorted.begin()->first ;
	}
	vars = sorted.begin()->second.first ;
	return err ;
}
	
std::vector<double> GeneticAlgorithmOptimizer::getValues() const
{
	std::vector<double> ret ;
	
	for(size_t i = 0 ; i < vars.size() ; i++)
		ret.push_back(vars[i]);
	for(size_t i = 0 ; i < lowLevelVars.size() ; i++)
		ret.push_back(*lowLevelVars[i]);
	
	return ret ;
}



} ;
