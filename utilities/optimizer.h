#ifndef GA_OPTIMIZER_H
#define GA_OPTIMIZER_H

#include <vector>
#include <complex>
#include "../geometry/geometry_base.h"
#include "../polynomial/vm_base.h"

namespace Mu
{
class GeneticAlgorithmOptimizer 
{
protected:
	friend class Segment ;
	friend class Geometry ;
	std::vector<double> vars ;
	std::vector<std::pair<double, double> > bounds ; 
	std::vector<double *> lowLevelVars ;
	std::vector<std::pair<double, double> >  lowLevelbounds ;
	
	Function objectiveFunction ;
	
	double applyFunction(const std::vector< double >& vals) const;
public:
	GeneticAlgorithmOptimizer(std::vector<double> vars, std::vector<std::pair<double, double> >  bounds, std::vector<double *> lowLevelVars,  
							  std::vector<std::pair<double, double> >  lowLevelbounds, const Function & objectiveFunction) ;
	GeneticAlgorithmOptimizer(const std::vector<double> & vars, const std::vector<std::pair<double, double> >  & bounds, const Function & objectiveFunction) ;
							  
	double optimize(double eps = 0.0001, int Maxit = -1, int population = -1, double elitism = .1) ;
	
	std::vector<double> getValues() const ;
} ;



} ;

#endif
