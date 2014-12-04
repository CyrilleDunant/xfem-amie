#include "crackmanager.h"
    
using namespace Amie ;

CrackManager::CrackManager(BranchedCrack * first, double criticalEnergy, double minCRadius, double maxExpansion) : EnrichmentManager(first), criticalEnergy(criticalEnergy), minCRadius(1e-3), findRadius(false), findExtension(false), iteration(0), maxExpansion(maxExpansion) {} ;
    
bool CrackManager::step(double dt, Vector * v, Mesh< DelaunayTriangle, DelaunayTreeItem >* dtree) 
{ 
    if(stable) //checkpoint: start iteration
    {
        stable = false ;
        findRadius = true ;
        
        double maxScore = -1;
        movingSet.clear();
        iteration = 0 ;
        for(auto & i : featureSet)
        {
            for(auto & j : static_cast<BranchedCrack *>(i)->getTips() )
            {
                std::pair<double, double> score = static_cast<BranchedCrack *>(i)->computeJIntegralAtTip(j, dtree) ;
                double scalarScore = (score.first*score.first+score.second*score.second)/criticalEnergy -1 ;
                if(scalarScore > maxScore)
                {
                    findRadius = true ;
                    Vector dr(2); dr[0] = 1 ; dr[1] = 0. ;
                    Matrix rotate0(2,2) ;
                    rotate0[0][0] = cos(j.second) ; rotate0[0][1] = sin(j.second) ;
                    rotate0[1][0] = -sin(j.second) ; rotate0[1][1] = cos(j.second) ;
                    
                    dr = rotate0*dr ;
                    
                    maxScore = scalarScore ;
                    movingSet.clear();
                    movingSet.push_back(std::make_pair(j.first,static_cast<BranchedCrack *>(i)));
                    
                    startAngles.clear();
                    startAngles.push_back(static_cast<BranchedCrack *>(i)->propagationAngleFromTip(j, dtree));
                    
                    
                    //there is a minimum radius of curvature.
                    //max is infty.
                    
                    downRadius = minCRadius ; ;
                    upRadius = minCRadius*1e4 ;
                    currentRadius = (downRadius+upRadius)*.5 ; 
                    Point center(-dr[1], dr[0]) ;
                    centers.clear();
                    centers.push_back(center);
                    center /=center.norm() ;
                    center *= currentRadius ;
                    center += *j.first ;
                    
                    upExtension = maxExpansion ;
                    downExtension = 0 ;
                    currentExtension = maxExpansion*2. ;
                    
                    double alpha = 0.5*upExtension/currentRadius ;
                    rotate0[0][0] = cos(alpha) ; rotate0[0][1] = sin(alpha) ;
                    rotate0[1][0] = -sin(alpha) ; rotate0[1][1] = cos(alpha) ;
                    
                    Point e0 = *j.first - center ;
                    e0 *=rotate0 ;
                    e0 += center ;
                    Point e1 = e0-center ;
                    e1 *=rotate0 ;
                    e1+= center ;
                    for(auto k : movingTips)
                    {
                        delete k.first.first ;
                        delete k.second.first ;
                    }
                    movingTips.clear();
                    movingTips.push_back(std::make_pair(std::make_pair(new Point(e0), 0.), std::make_pair(new Point(e1), startAngles.back()+2.*alpha)));
                }
            }
        }
        if(movingSet.empty())
            stable = true ;
        else
        {
            for(size_t i = 0 ; i < movingSet.size() ; i++)
            {
                movingSet[i].second->grow(movingSet[i].first, movingTips[i].first.first) ;
                movingSet[i].second->grow(movingTips[i].first.first, movingTips[i].second.first) ;
            }
        }
    }
    else if(findRadius)
    {
        for(size_t i = 0 ; i < movingSet.size() ; i++)
        {
            double apparentAngle = movingSet[i].second->propagationAngleFromTip(movingTips[i].second, dtree) ;
            double effectiveAngle = 0.5*upExtension/currentRadius ;
            
            if(apparentAngle > effectiveAngle)
                upRadius = currentRadius ;
            else
                downRadius = currentRadius ;
            
            currentRadius = (downRadius+upRadius)*.5 ; 
            double alpha = 0.5*upExtension/currentRadius ;
            Matrix rotate0(2,2) ;
            rotate0[0][0] = cos(alpha) ; rotate0[0][1] = sin(alpha) ;
            rotate0[1][0] = -sin(alpha) ; rotate0[1][1] = cos(alpha) ;
            
            
            Point center = centers[i];
            center /=center.norm() ;
            center *= currentRadius ;
            center += *movingSet[i].first ;
                    
            Point e0 = *movingSet[i].first - center ;
            e0 *=rotate0 ;
            e0 += center ;
            Point e1 = e0-center ;
            e1 *=rotate0 ;
            e1+= center ;
            movingSet[i].second->move(movingTips[i].first.first, e0) ;
            movingSet[i].second->move(movingTips[i].second.first, e1) ;
        }
        if(iteration++ > 8)
        {
            findRadius = false ;
            iteration = 0 ;
            findExtension = true ;
        }
    }
    else if (findExtension)
    {
        double maxScore = -1 ;
        for(auto & i : featureSet)
        {
            for(auto & j : static_cast<BranchedCrack *>(i)->getTips() )
            {
                std::pair<double, double> score = static_cast<BranchedCrack *>(i)->computeJIntegralAtTip(j, dtree) ;
                double scalarScore = (score.first*score.first+score.second*score.second)/criticalEnergy -1 ;
                if(scalarScore > maxScore)
                    maxScore = scalarScore ;
            }
        }
        
        std::pair<double, double>  cscore = movingSet[0].second->computeJIntegralAtTip(movingTips[0].second, dtree) ;
        double currentScore = (cscore.first*cscore.first+cscore.second*cscore.second)/criticalEnergy -1 ;
        
        if(currentScore < maxScore - POINT_TOLERANCE_3D)
        {
            upExtension = currentExtension ;
            currentExtension = (upExtension+downExtension)*.5 ;
        }
        else
        {
            downExtension = currentExtension ;
            currentExtension = (upExtension+downExtension)*.5 ; 
        }
        

        double alpha = 0.5*currentExtension/currentRadius ;
        Matrix rotate0(2,2) ;
        rotate0[0][0] = cos(alpha) ; rotate0[0][1] = sin(alpha) ;
        rotate0[1][0] = -sin(alpha) ; rotate0[1][1] = cos(alpha) ;
        
        
        Point center = centers[0];
        center /=center.norm() ;
        center *= currentRadius ;
        center += *movingSet[0].first ;
                
        Point e0 = *movingSet[0].first - center ;
        e0 *=rotate0 ;
        e0 += center ;
        Point e1 = e0-center ;
        e1 *=rotate0 ;
        e1+= center ;
        movingSet[0].second->move(movingTips[0].first.first, e0) ;
        movingSet[0].second->move(movingTips[0].second.first, e1) ;
        
        if(iteration++ > 8)
        {
            findExtension = false ;
            iteration = 0 ;
            stable = true ;
        }
    }

    
};

void CrackManager::setCriticalEnergy(double e)
{
    criticalEnergy = e ;
}

void CrackManager::setMaximumRadius(double r)
{
    minCRadius = r*1e4 ;
}

void CrackManager::setMinimumRadius(double r)
{
    minCRadius = r ;
}

void CrackManager::setMaximumExpansion(double m)
{
    maxExpansion = m ;
}

bool CrackManager::step(double dt, Vector * v, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * dtree) 
{ 

};

bool CrackManager::converged() { return stable ; } ;