// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "optimizer.h"
#include "itoa.h"

namespace Amie
{

GeneticAlgorithmOptimizer::GeneticAlgorithmOptimizer(std::vector<double> vars, std::vector<std::pair<double, double> >  bounds, std::vector<double *> lowLevelVars, std::vector<std::pair<double, double> >  lowLevelbounds, const Function & objectiveFunction) : vars(vars), bounds(bounds), lowLevelVars(lowLevelVars), lowLevelbounds(lowLevelbounds), objectiveFunction(objectiveFunction), lowLevelFunction(nullptr), generator(nullptr)
{

}

GeneticAlgorithmOptimizer::GeneticAlgorithmOptimizer(std::vector<double *> lowLevelVars,  std::vector<std::pair<double, double> >  lowLevelbounds, double (*lowLevelFunction)()) : lowLevelVars(lowLevelVars), lowLevelbounds(lowLevelbounds), lowLevelFunction(lowLevelFunction), generator(nullptr)
{
}

GeneticAlgorithmOptimizer::GeneticAlgorithmOptimizer(std::vector<double *> lowLevelVars,  std::vector<std::pair<double, double> >  lowLevelbounds, MicrostructureGenerator * generator) : lowLevelVars(lowLevelVars), lowLevelbounds(lowLevelbounds), lowLevelFunction(nullptr), generator(generator)
{
}

GeneticAlgorithmOptimizer::GeneticAlgorithmOptimizer(const std::vector<std::pair<std::string, double> > & nvars, const std::vector<std::pair<double, double> >  & nbounds, const Function & objectiveFunction) : namedVars(nvars), namedVarsBounds(nbounds), objectiveFunction(objectiveFunction), lowLevelFunction(nullptr), generator(nullptr)
{
}

std::vector<double> GeneticAlgorithmOptimizer::getRoots( double a, double b, double c) 
{
    std::vector<double> ret ;
    double delta = b*b-4.*a*c ;
    if(std::abs(delta) < POINT_TOLERANCE)
       ret.push_back( -b/(2.*a) ) ;
    else if(delta > POINT_TOLERANCE)
    {
        ret.push_back( (-b+std::sqrt(delta))/(2.*a) ) ;
        ret.push_back( (-b-std::sqrt(delta))/(2.*a) ) ;
    }
    return ret ;
}


double GeneticAlgorithmOptimizer::applyFunction(const std::vector<double> & vals) const
{
    // x, y, z, t, u, v, w

    switch (vals.size())
    {
    case 0 :
        return VirtualMachine().eval(objectiveFunction) ;
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

double GeneticAlgorithmOptimizer::optimize(double eps, int Maxit, int population, double elitism, double factor)
{
    if(lowLevelFunction)
        return lowLevelOptimize(eps, Maxit, population, elitism, factor) ;

    if(generator)
        return generatorOptimize(eps, Maxit, population, elitism, factor) ;

    int it = 0 ;
    int maxit = Maxit ;
    if(maxit < 0)
        maxit = 100 ;

    if(population < 0)
        population = 10 ;

    std::vector< std::vector<double> > individuals ;
    std::vector< std::vector<double> > nindividuals ;
    std::vector< std::vector<double> > llindividuals ;
    for(int i = 0 ; i < population ; i++)
    {
        std::vector<double> newindividual ;
        for(size_t j = 0 ;  j < vars.size() ; j++)
        {
            newindividual.push_back((bounds[j].second-bounds[j].first)*rand()/RAND_MAX + bounds[j].first);
        }
        individuals.push_back(newindividual);
    }

    for(int i = 0 ; i < population ; i++)
    {
        std::vector<double> newindividual ;
        for(size_t j = 0 ;  j < namedVars.size() ; j++)
        {
            newindividual.push_back((namedVarsBounds[j].second-namedVarsBounds[j].first)*rand()/RAND_MAX + namedVarsBounds[j].first);
        }
        nindividuals.push_back(newindividual);
    }

    for(int i = 0 ; i < population ; i++)
    {
        std::vector<double> newllindividual ;
        for(size_t j = 0 ;  j < lowLevelVars.size() ; j++)
        {
            newllindividual.push_back((lowLevelbounds[j].second-lowLevelbounds[j].first)*rand()/RAND_MAX + lowLevelbounds[j].first);
        }
        llindividuals.push_back(newllindividual);
    }

    double err = eps*100 ;
    std::map<double, std::vector<std::vector<double> > > sorted ;
    while (err > eps && it < maxit )
    {
        sorted.clear();
        //evaluate
        for(int i = 0 ; i < population ; i++)
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
        for(int i = 0 ; i < std::min((int)(population*elitism), (int)sorted.size()) ; i++)
        {
            auto iter = sorted.begin() ;
            for(int j = 0 ; j < i ; j++)
                iter++ ;
            newindividuals.push_back(iter->second[0]);
            newnindividuals.push_back(iter->second[1]);
            newllindividuals.push_back(iter->second[2]);
        }
        std::default_random_engine generator;

        //reproduce - the rest fills the available slots and mutates
        int i = 0 ;
        while( newindividuals.size() < individuals.size())
        {
            for(size_t n = 0 ; n < 4 && newindividuals.size() <= individuals.size(); n++)
            {
                auto iter = sorted.begin() ;
                for(int j = 0 ; j < i ; j++)
                    iter++ ;
                newindividuals.push_back(iter->second[0]);
                newnindividuals.push_back(iter->second[1]);
                newllindividuals.push_back(iter->second[2]);
                double factor = 0.85;
                for(size_t j = 0 ; j < newindividuals.back().size() ; j++)
                {
                    double sigma = err*(namedVarsBounds[j].second-namedVarsBounds[j].first)/population*factor ;
                    sigma = std::min(sigma, (bounds[j].second-bounds[j].first)/2) ;
                    std::normal_distribution< double > distribution(newindividuals.back()[j], sigma) ;
                    newindividuals.back()[j] = distribution(generator) ; //newindividuals.back()[j]*(1.-err*population*.5) + err*population*.5*(bounds[j].second-bounds[j].first)*random()/RAND_MAX ;
                    newindividuals.back()[j] = std::max(newindividuals.back()[j], bounds[j].first) ;
                    newindividuals.back()[j] = std::min(newindividuals.back()[j], bounds[j].second) ;
                }
                for(size_t j = 0 ; j < newnindividuals.back().size() ; j++)
                {
                    double sigma = err*(namedVarsBounds[j].second-namedVarsBounds[j].first)/population*factor ;
                    sigma = std::min(sigma, (namedVarsBounds[j].second-namedVarsBounds[j].first)/2) ;
                    std::normal_distribution< double > distribution(newnindividuals.back()[j], sigma) ;
                    newnindividuals.back()[j] = distribution(generator) ; //newindividuals.back()[j]*(1.-err*population*.5) + err*population*.5*(bounds[j].second-bounds[j].first)*random()/RAND_MAX ;
                    newnindividuals.back()[j] = std::max(newnindividuals.back()[j], namedVarsBounds[j].first) ;
                    newnindividuals.back()[j] = std::min(newnindividuals.back()[j], namedVarsBounds[j].second) ;
                }
                for(size_t j = 0 ; j < newllindividuals.back().size() ; j++)
                {
                    double sigma = err*(lowLevelbounds[j].second-lowLevelbounds[j].first)/population*factor ;
                    sigma = std::min(sigma, (lowLevelbounds[j].second-lowLevelbounds[j].first)/2) ;
                    std::normal_distribution< double > distribution(newllindividuals.back()[j], sigma) ;
                    newllindividuals.back()[j] = distribution(generator) ; //newindividuals.back()[j]*(1.-err*population*.5) + err*population*.5*(bounds[j].second-bounds[j].first)*random()/RAND_MAX ;
                    newllindividuals.back()[j] = std::max(newllindividuals.back()[j], lowLevelbounds[j].first) ;
                    newllindividuals.back()[j] = std::min(newllindividuals.back()[j], lowLevelbounds[j].second) ;
                }
            }

            i++ ;
            i = std::min(i, (int)sorted.size()-1) ;
        }
        std::cout << sorted.begin()->first << std::endl ;
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
    std::default_random_engine generator ;
    std::uniform_real_distribution< double > uniform(0,1) ;

    for(int i = 0 ; i < population ; i++)
    {
        std::vector<double> newllindividual ;
        for(size_t j = 0 ;  j < lowLevelVars.size() ; j++)
        {
            newllindividual.push_back((lowLevelbounds[j].second-lowLevelbounds[j].first)*uniform(generator) + lowLevelbounds[j].first);
        }
        llindividuals.push_back(newllindividual);
    }

    double err = eps*100 ;
    double err0 = eps*100 ;
    std::map<double, std::vector<double> > sorted ;
    
 
    while (err > eps && it < maxit )
    {
        sorted.clear();
        //evaluate
        for(int i = 0 ; i < population ; i++)
        {
            while(true)
            {
                for(size_t j = 0 ;  j < lowLevelVars.size() ; j++)
                {
                    *lowLevelVars[j] = llindividuals[i][j] ;
                }

                double llf = lowLevelFunction() ;
                if(!std::isnan(llf))
                {
                    sorted[llf] = llindividuals[i] ;
                    break ;
                }

                for(size_t j = 0 ;  j < lowLevelVars.size() ; j++)
                    llindividuals[i][j] = (lowLevelbounds[j].second-lowLevelbounds[j].first)*uniform(generator) + lowLevelbounds[j].first;

            }
        }

        std::vector< std::vector<double> > newllindividuals ;

        //reproduce - elite individuals are kept
        for(int i = 0 ; i < std::max(1, (int)round(elitism*population)) ; i++)
        {
            auto iter = sorted.begin() ;
            auto iterend = sorted.end() ;
            if(iter == iterend)
                break ;
            for(size_t j = 0 ; (int)j < i-1 && j < sorted.size()-1 ; j++)
                iter++ ;
            newllindividuals.push_back(iter->second);
        }

        //reproduce - the rest fills the available slots and mutates
        int i = 0 ;
        while( newllindividuals.size() != llindividuals.size())
        {
            double test = uniform(generator) ;
            if(sorted.size() && test >= (double)i/((double)sorted.size()))
            {
                auto iter = sorted.begin() ;
                for(size_t j = 0 ; (int)j < i-1 && j < sorted.size()-1 ; j++)
                    iter++ ;
                newllindividuals.push_back(iter->second);

                for(size_t j = 0 ; j < newllindividuals.back().size() ; j++)
                {
                    double sigma = (lowLevelbounds[j].second-lowLevelbounds[j].first)*factor/err0 ;
                    sigma = std::min(sigma, (lowLevelbounds[j].second-lowLevelbounds[j].first)*.5) ;
                    std::normal_distribution< double > distribution(newllindividuals.back()[j], sigma) ;
                    newllindividuals.back()[j] = distribution(generator) ; //newindividuals.back()[j]*(1.-err*population*.5) + err*population*.5*(bounds[j].second-bounds[j].first)*random()/RAND_MAX ;
                    newllindividuals.back()[j] = std::max(newllindividuals.back()[j], lowLevelbounds[j].first) ;
                    newllindividuals.back()[j] = std::min(newllindividuals.back()[j], lowLevelbounds[j].second) ;
                }

            }
            else
            {
                auto iter = sorted.begin() ;
                auto iterend = sorted.end() ;

                for(int j = 0 ; j < i && iter != iterend ; j++)
                    iter++ ;
                if(iter != iterend)
                    newllindividuals.push_back(iter->second);
            }
            i++ ;
            if(i >= (int)sorted.size())
                i = 0 ;
        }

        llindividuals = newllindividuals ;
        for(size_t j = 0 ;  j < lowLevelVars.size() ; j++)
            *lowLevelVars[j] = sorted.begin()->second[j] ;

        // we call the function, because it might have some border effects
        lowLevelFunction() ;


        for(size_t i = 0 ; i < lowLevelVars.size() ; i++)
            std::cout << *lowLevelVars[i] << " ("<<exp(*lowLevelVars[i]) << "), "<< std::flush ;
        //iterate
        err = sorted.begin()->first ;
        if(it == 0)
            err0 = err ;
        std::cout << err <<  std::endl ;
        it++ ;
    }
    std::cout << it << std::endl;

// 	for(size_t i = 0 ; i < lowLevelVars.size() ; i++)
// 		*lowLevelVars[i] = sorted.begin()->second[i] ;

    return err ;
}

double GeneticAlgorithmOptimizer::generatorOptimize(double eps, int Maxit, int population, double elitism, double factor)
{
    int it = 0 ;
    int maxit = Maxit ;
    if(maxit < 0)
        maxit = 100 ;

    if(population < 0)
        population = 10 ;

    std::vector< std::vector<double> > llindividuals ;

    std::default_random_engine random_generator ;
    std::uniform_real_distribution< double > uniform(0,1) ;
    for(int i = 0 ; i < population ; i++)
    {
        std::vector<double> newllindividual ;
        for(size_t j = 0 ;  j < lowLevelVars.size() ; j++)
        {
            newllindividual.push_back((lowLevelbounds[j].second-lowLevelbounds[j].first)*uniform(random_generator) + lowLevelbounds[j].first);
        }
        llindividuals.push_back(newllindividual);
    }

    double err = eps*100 ;
    std::map<double, std::vector<double> > sorted ;
    while (err > eps && it < maxit )
    {
        //evaluate
        for(int i = (it != 0)*(std::max(1, (int)round(elitism*population))-1) ; i < population ; i++)
        {
            while(true)
            {
                for(size_t j = 0 ;  j < lowLevelVars.size() ; j++)
                {
                    *lowLevelVars[j] = llindividuals[i][j] ;
                }

                double llf = generator->score() ;
                if(!std::isnan(llf))
                {
                    sorted[llf] = llindividuals[i] ;
                    break ;
                }

                for(size_t j = 0 ;  j < lowLevelVars.size() ; j++)
                    llindividuals[i][j] = (lowLevelbounds[j].second-lowLevelbounds[j].first)*uniform(random_generator) + lowLevelbounds[j].first;

            }
        }
        auto iter = sorted.begin() ;
        for(int i = 0 ;  i < population && iter !=sorted.end() ; i++)
            iter++ ;
        if(it != 0)
            sorted.erase(iter, sorted.end());
        std::vector< std::vector<double> > newllindividuals ;

        //reproduce - elite individuals are kept
        for(int i = 0 ; i < std::max(1, (int)round(elitism*population)) ; i++)
        {
            auto iter = sorted.begin() ;
            auto iterend = sorted.end() ;
            if(iter == iterend)
                break ;
            for(size_t j = 0 ; (int)j < i-1 && j < sorted.size()-1 ; j++)
                iter++ ;
            newllindividuals.push_back(iter->second);
        }

        //reproduce - the rest fills the available slots and mutates
        size_t i = 0 ;
        while( newllindividuals.size() != llindividuals.size())
        {
            double test = uniform(random_generator) ;
            if(sorted.size() && test >= (double)i/((double)sorted.size()))
            {
                auto iter = sorted.begin() ;
                for(size_t j = 0 ; j < i-1 && j < sorted.size()-1 ; j++)
                    iter++ ;
                newllindividuals.push_back(iter->second);

                for(size_t j = 0 ; j < newllindividuals.back().size() ; j++)
                {
                    double sigma = (lowLevelbounds[j].second-lowLevelbounds[j].first)*.1 ;
// 					sigma = std::min(sigma, (lowLevelbounds[j].second-lowLevelbounds[j].first)*.5) ;
                    std::normal_distribution< double > distribution(newllindividuals.back()[j], sigma) ;
                    newllindividuals.back()[j] = distribution(random_generator) ; //newindividuals.back()[j]*(1.-err*population*.5) + err*population*.5*(bounds[j].second-bounds[j].first)*random()/RAND_MAX ;
                    newllindividuals.back()[j] = std::max(newllindividuals.back()[j], lowLevelbounds[j].first) ;
                    newllindividuals.back()[j] = std::min(newllindividuals.back()[j], lowLevelbounds[j].second) ;
                }

            }
            else
            {
                auto iter = sorted.begin() ;
                auto iterend = sorted.end() ;

                for(size_t j = 0 ; j < i && iter != iterend ; j++)
                    iter++ ;
                if(iter != iterend)
                    newllindividuals.push_back(iter->second);
            }
            i++ ;
            if(i >= sorted.size())
                i = 0 ;
        }

        llindividuals = newllindividuals ;
        for(size_t j = 0 ;  j < lowLevelVars.size() ; j++)
            *lowLevelVars[j] = sorted.begin()->second[j] ;

        // we call the function, because it might have some border effects
        generator->score() ;
        generator->print() ;

        for(size_t i = 0 ; i < lowLevelVars.size() ; i++)
            std::cout << *lowLevelVars[i] << " ("<<exp(*lowLevelVars[i]) << "), "<< std::flush ;
        //iterate
        err = sorted.begin()->first ;
        std::cout << err <<  std::endl ;
        it++ ;
    }
    std::cout << it << std::endl;

// 	for(size_t i = 0 ; i < lowLevelVars.size() ; i++)
// 		*lowLevelVars[i] = sorted.begin()->second[i] ;

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

LeastSquaresApproximation::~LeastSquaresApproximation()
{
    delete Q ;
}

LeastSquaresApproximation::LeastSquaresApproximation(const Vector & measures, const Matrix & linearModel) : measures(measures), linearModel(linearModel), parameters(0., linearModel.numRows()), linearModelChanged(true), X0t(linearModel.numCols(), linearModel.numRows()), X0tX0(linearModel.numRows(), linearModel.numRows())
{
    Q = nullptr ;
    X0t   = linearModel.transpose() ;
    X0tX0 = linearModel*X0t ;

}

Vector LeastSquaresApproximation::getApproximation() const
{
    Vector ret(linearModel.numCols()) ;
    for(size_t i = 0 ; i < linearModel.numCols() ; ++i)
    {
        ret[i] = 0 ;
        for(size_t j = 0 ; j < linearModel.numRows() ; ++j)
        {
            ret[i] += linearModel[j][i]*parameters[j] ;
        }
    }

    return ret ;
}



void LeastSquaresApproximation::printParameters() const
{

    for(size_t i = 0 ; i < parameters.size() ; i++)
        std::cout << parameters[i] << std::endl ;
}

void LeastSquaresApproximation::setParameterValue(size_t i, double v)
{
    fixedValues.push_back(std::make_pair(i, v));
}

void LeastSquaresApproximation::clearParameterValues()
{
    fixedValues.clear();
    delete Q ;
    Q = nullptr ;
}

void LeastSquaresApproximation::setLinearModel(int i, int j, double v)
{
    linearModel[i][j] = v ;
    linearModelChanged = true ;
}

double LeastSquaresApproximation::optimize()
{
    Matrix * A = nullptr;
    Vector consts(fixedValues.size()) ;
    for(size_t i = 0 ; i < consts.size() ; i++)
        consts[i] = fixedValues[i].second ;

    if(!fixedValues.empty())
    {
        //setting the constrained values.
        delete Q ;
        Q = new Matrix(parameters.size(), parameters.size()) ;

        //The Q matrix is the line-swap matrix which puts all
        //the constrained values together
        //first, it is initialised to the identity.
        for( size_t i = 0 ; i < parameters.size() ; i++ )
        {
            (*Q)[i][i] = 1 ;
        }

        //the lines are swapped as we find new constraints
        int c = 0 ;
        Matrix S(parameters.size(), parameters.size()) ;
        for( size_t j = 0 ; j < parameters.size() ; j++ )
        {
            S[j][j] = 1 ;
        }
        for( size_t i = 0 ; i < fixedValues.size() ; i++ )
        {
            for( size_t j = 0 ; j < parameters.size() ; j++ )
            {
                std::swap(S[fixedValues[i].first][j], S[c][j]);
            }
            c++ ;
        }
        S *= *Q ;
        (*Q) = S ;

        // The problem changes and becomes:
        A = new Matrix((*Q)*linearModel) ;
        Matrix A1(fixedValues.size(), linearModel.numCols() ) ;
        for(size_t i = 0 ; i < fixedValues.size() ; i++)
        {
            for(size_t j = 0 ; j < linearModel.numCols() ; j++)
            {
                A1[i][j] = (*A)[i][j] ;
            }
        }

        Matrix A2(linearModel.numRows()-fixedValues.size(), linearModel.numCols() ) ;
        for(size_t i = 0 ; i < linearModel.numRows()-fixedValues.size() ; i++)
        {
            for(size_t j = 0 ; j < linearModel.numCols() ; j++)
            {
                A2[i][j] = (*A)[i+fixedValues.size()][j] ;
            }
        }

        measures = measures-consts*A1 ;
        linearModel.resize(A2.numRows(), A2.numCols());
        linearModel = A2 ;
        parameters.resize(A2.numRows());
        linearModelChanged = true ;
    }


    if(linearModelChanged)
    {
        X0t.resize(linearModel.numCols(), linearModel.numRows());
        X0t = linearModel.transpose() ;
        X0tX0.resize(linearModel.numRows(), linearModel.numRows());
        X0tX0 = linearModel*X0t ;
    }
    Vector lb0 = linearModel*measures ;
    parameters = 0 ;
    parameters = solveSystem(X0tX0, lb0, parameters) ;

    if(Q)
    {
        Vector newp(parameters.size()+consts.size()) ;
        for(size_t i = 0 ; i < consts.size() ; i++)
            newp[i] = consts[i] ;
        for(size_t i = 0 ; i < parameters.size() ; i++)
            newp[i+consts.size()] = parameters[i] ;

        newp = newp*(*Q) ;

        parameters.resize(newp.size());
        parameters = newp;
        linearModel.resize(A->numRows(), A->numCols()) ;
        linearModel = *A ;
        X0t = linearModel.transpose() ;
        X0tX0.resize(A->numRows(), A->numRows());
        X0tX0 = linearModel*X0t ;
        delete A ;
    }

    return std::abs((Vector)(X0tX0*parameters)-lb0).max() ;
    linearModelChanged = false ;
}

}
