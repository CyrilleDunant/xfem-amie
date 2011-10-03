// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef GA_OPTIMIZER_H
#define GA_OPTIMIZER_H

#include <vector>
#include <complex>
#include "../geometry/geometry_base.h"
#include "../polynomial/vm_base.h"
#include "../features/microstructuregenerator.h"

namespace Mu
{

class GeneticAlgorithmOptimizer 
{
protected:

	std::vector<double> vars ;
	std::vector<std::pair<double, double> > bounds ; 
	std::vector<std::pair<std::string, double> > namedVars ;
	std::vector<std::pair<double, double> > namedVarsBounds ; 
	std::vector<double *> lowLevelVars ;
	std::vector<std::pair<double, double> >  lowLevelbounds ;
	
	Function objectiveFunction ;
	double (*lowLevelFunction)() ;
	MicrostructureGenerator * generator ;
	
	double applyFunction(const std::vector< double >& vals) const;
	double lowLevelOptimize(double eps = 0.0001, int Maxit = -1, int population = -1, double elitism = .1, double factor = .65) ;
	double generatorOptimize(double eps = 0.0001, int Maxit = -1, int population = -1, double elitism = .1, double factor = .65) ;
	
public:
	GeneticAlgorithmOptimizer(std::vector<double> vars, std::vector<std::pair<double, double> >  bounds, std::vector<double *> lowLevelVars,  std::vector<std::pair<double, double> >  lowLevelbounds, const Function & objectiveFunction) ;
		
	GeneticAlgorithmOptimizer(std::vector<double *> lowLevelVars,  std::vector<std::pair<double, double> >  lowLevelbounds, double (*lowLevelFunction)()) ;
	
	GeneticAlgorithmOptimizer(std::vector<double *> lowLevelVars,  std::vector<std::pair<double, double> >  lowLevelbounds, MicrostructureGenerator * generator) ;
	
	GeneticAlgorithmOptimizer(const std::vector<std::pair<std::string, double> > & nvars, const std::vector<std::pair<double, double> >  & nbounds, const Function & objectiveFunction) ;
							  
	double optimize(double eps = 0.0001, int Maxit = -1, int population = -1, double elitism = .1, double factor = .65) ;
	
	std::vector<std::pair<std::string, double> > getValues() const ;
} ;

class LeastSquaresApproximation
{
protected:
	Vector measures ;
	Matrix linearModel ;
	Vector parameters ;
	
	std::vector<std::pair<size_t, double> > fixedValues ;
	
	bool linearModelChanged ;
	Matrix X0t ;
	Matrix X0tX0 ;
	Matrix * Q ;
public:
	
	LeastSquaresApproximation(const Vector & measures, const Matrix & linearModel) ;
	~LeastSquaresApproximation() ;
	
	
	void setLinearModel(int i, int j, double v) ;
	const Matrix & getLinearModel() const { return linearModel ;} ;
	const Vector & getMeasures() const { return measures ; } ;
	void setMeasures(const Vector &m)  
	{ 
		measures.resize(m.size()) ;
		measures = m; 
	} ;
	
	double optimize();
	Vector getApproximation() const ;
	void setParameterValue(size_t i, double v) ;
	void clearParameterValues() ;
	
	const Vector & getParameters() const { return parameters ; } ;
    void printParameters() const;
	
	
} ;

} ;

#endif
