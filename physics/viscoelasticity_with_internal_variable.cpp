
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
//
// Copyright: See COPYING file that comes with this distribution
//

#include "viscoelasticity_with_internal_variable.h"
#include "../mesher/delaunay.h" 
#include "../polynomial/vm_base.h" 
#include "../features/boundarycondition.h"

using namespace Mu ;




ViscoElasticity::ViscoElasticity( double _tau_k, double _tau_g, Vector g, Vector k) : LinearForm(Matrix(6,6)), tau_g(g.size()), tau_k(k.size()), g(g), k(k), a_g(6,g.size() ), a_k(k.size()), average_delta_sigma(6)
{
	num_dof = 3 ;
	
	tau_g[0] = 1 ;
	for(size_t i = 1 ;  i < g.size() ; i++)
	{
		tau_g[i] = pow(10, i-1)*_tau_g ;
	}
	
	tau_k[0] = 1 ;
	for(size_t i = 1 ;  i < k.size() ; i++)
	{
		tau_k[i] = pow(10, i-1)*_tau_k ;
	}
	
	
	double E = (9.*k[0]*g[0])/(3.*k[0]+g[0]);
	double nu = (3.*k[0]-2.*g[0])/(6.*k[0]+2.*g[0]);
	
	param[0][0]  = 1.-nu ; param[0][1] = nu   ; param[0][2] = nu ;
	param[1][0]  = nu   ; param[1][1] = 1.-nu ; param[1][2] = nu ;
	param[2][0]  = nu   ; param[2][1] = nu   ; param[2][2] = 1.-nu  ;
	param[3][3]  = 0.5 - nu ;
	param[4][4]  = 0.5 - nu ;
	param[5][5]  = 0.5 - nu ;
	param *= E/((1.+nu)*(1.-2.*nu)) ;
	
// 		std::cout << "E = " << E << ", nu = " << nu << std::endl; 
		
	time_d = true;
	
}
	
ViscoElasticity::~ViscoElasticity() { } ;


void ViscoElasticity::apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix &ret, VirtualMachine * vm) const
{
	std::vector<Variable> v ;
	v.push_back(XI);
	v.push_back(ETA);
	v.push_back(ZETA);
	vm->ieval(Gradient(p_i) * param * Gradient(p_j, true), gp, Jinv,v,ret) ;
}

void ViscoElasticity::step(double timestep, ElementState & currentState)
{
	
		// computing G
	Vector lambda_g(1., g.size()) ;
	Vector gama_g(0., g.size()) ;
	for(size_t i = 1 ; i < lambda_g.size() ; i++)
	{
		
		if ((timestep/tau_g[i]) < 1e-5)
		{
			gama_g[i]=0;
			lambda_g[i]=1;
		}
		
		gama_g[i]= 1.-exp(-timestep/tau_g[i]);
		lambda_g[i] = gama_g[i] * tau_g[i]/timestep ;
	}
	
	double G = 0;
	for(size_t i = 0 ; i < g.size() ; i++)
	{
		G += g[i]*lambda_g[i] ;
	}
	
	Vector lambda_k(1., k.size()) ;
	Vector gama_k(0., k.size()) ;
	for(size_t i = 1 ; i < lambda_k.size() ; i++)
	{
		
		if ((timestep/tau_k[i]) < 1e-5)
		{	
			gama_k[i]=0;
			lambda_k[i]=1;
		}
		
		gama_k[i] = 1.-exp(-timestep/tau_k[i]);
		lambda_k[i] = gama_k[i] * tau_k[i]/timestep ;
	}
	
	double K = 0;
	for(size_t i = 0 ; i < k.size() ; i++)
	{
		K += k[i]*lambda_k[i] ;
	}
	
	Vector Deltastrain_g(3);
	Vector previousStrain_g(3);
	Vector Deltastrain = currentState.getDeltaStrain(currentState.getParent()->getCenter()) ;
	Vector previousStrain = currentState.getStrain(currentState.getParent()->getCenter())-Deltastrain ;
	double compoundDeltaStrain = std::accumulate(&Deltastrain[0], &Deltastrain[2], (double)(0.)) ;
	double previousCompoundStrain = std::accumulate(&previousStrain[0], &previousStrain[2], (double)(0.)) ;
	Deltastrain_g-= std::accumulate(&Deltastrain[0], &Deltastrain[2], (double)(0.))/3. ;
	previousStrain_g -=std::accumulate(&previousStrain[0], &previousStrain[2], (double)(0.))/3. ;
	
	for(size_t i = 1 ; i < tau_k.size() ; i++)
	{
#warning this is only valid for linear elements
			
		a_k[i] = (1.-gama_k[i])*a_k[i]+compoundDeltaStrain*(1.-lambda_k[i])+previousCompoundStrain*gama_k[i] ;
	}
	
			
	for(size_t j = 0 ; j < 3 ; j++)
	{
#warning this is only valid for linear elements
			
		for(size_t i = 1 ; i < tau_g.size() ; i++)
			
			a_g[j][i] = a_g[j][i]*(1.-gama_g[i])+Deltastrain_g[j]*(1.-lambda_g[i])+previousStrain_g[j]*gama_g[i] ;
	}
	
	
	
	for(size_t j = 3 ; j < 6 ; j++)
	{
#warning this is only valid for linear elements
			
		for(size_t i = 1 ; i < tau_g.size() ; i++)
			a_g[j][i] = a_g[j][i]*(1.-gama_g[i])+Deltastrain[j]*(1.-lambda_g[i])+previousStrain[j]*gama_g[i] ;
	}
		
		
	param[0][0] = K+4./3.*G ; param[0][1] = K-2./3.*G ; param[0][2] = K-2./3.*G ;
	param[1][0] = K-2./3.*G ; param[1][1] = K+4./3.*G; param[1][2] =  K-2./3.*G ;
	param[2][0] = K-2./3.*G ; param[2][1] = K-2./3.*G; param[2][2] =  K+4./3.*G ;
	param[3][3] = G ;
	param[4][4] = G ;
	param[5][5] = G ;
		
		//update average_delta_sigma
	Vector average_delta_sigma_g(6) ;
	double average_delta_sigma_k(6) ;
// 		average_delta_segma_g(g.size());
		
		
	for(size_t i = 1 ; i < tau_k.size() ; i++)
	{
		
		average_delta_sigma_k += k[i]*gama_k[i]*(previousCompoundStrain-a_k[i]) ;
	}
	
	
	
	
	for(size_t j = 0 ; j < 3 ; j++)
	{
		
// 			for(size_t j = 0 ; j < a_g.numCols() ; j++)
		for(size_t i = 1 ; i < tau_g.size() ; i++)
		{
			average_delta_sigma_g[j] += g[i]*gama_g[i]*(previousStrain_g[j]-a_g[j][i]) ;
		}
	}
		
		
		
	for(size_t j = 3 ; j < 6 ; j++)
	{
		
// 			for(size_t k = 1 ; k < previousStrain.size() ; k++)
// 				std::cout << previousStrain[k] << std::endl ;
		for(size_t i = 1 ; i < tau_g.size() ; i++)
		{
// 				std::cout << previousStrain[j] << std::endl ;
			average_delta_sigma_g[j] += g[i]*gama_g[i]*(previousStrain[j]-a_g[j][i]) ;
		}
	}
	
		
	for(size_t j = 0 ; j <6 ; j++)
	{
		if ( j <3 )
		{
			// 			for(size_t k = 1 ; k < previousStrain.size() ; k++)
			// 				std::cout << previousStrain[k] << std::endl ;	
			average_delta_sigma[j] = -(average_delta_sigma_k+average_delta_sigma_g[j]);
			std::cout << j  << average_delta_sigma[j]<< std::endl ;
		}

		else
		{
			average_delta_sigma[j] = -(average_delta_sigma_g[j]);
			std::cout << j  << average_delta_sigma[j]<< std::endl ;
		}
	}
}
	
std::vector<BoundaryCondition * > ViscoElasticity::getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const
{
	std::vector<Variable> v ;
	v.push_back(XI);
	v.push_back(ETA);
	v.push_back(ZETA);
	Vector 	f =  VirtualMachine().ieval(Gradient(p_i, true)*average_delta_sigma, gp, Jinv,v) ;

	std::vector<BoundaryCondition * > ret ;
	if(f.size() == 2)
	{
		ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_XI, dynamic_cast<ElementarySurface *>(s.getParent()), id, f[0]));
		ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_ETA, dynamic_cast<ElementarySurface *>(s.getParent()), id, f[1]));
	}
	if(f.size() == 3)
	{
		ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_XI, dynamic_cast<ElementaryVolume *>(s.getParent()), id, f[0]));
		ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_ETA, dynamic_cast<ElementaryVolume *>(s.getParent()), id, f[1]));
		ret.push_back(new DofDefinedBoundaryCondition(SET_FORCE_ZETA, dynamic_cast<ElementaryVolume *>(s.getParent()), id, f[2]));
	}
	return ret ;
}

bool ViscoElasticity::hasInducedForces()
{
	return true ;
}

Form * ViscoElasticity::getCopy() const 
{
	return new ViscoElasticity(*this) ;
}




