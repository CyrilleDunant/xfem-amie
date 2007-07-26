
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
//
// Copyright: See COPYING file that comes with this distribution
//

#include "physics.h"
#include "../delaunay.h" 

using namespace Mu ;

BimaterialInterface::BimaterialInterface(Geometry * in, const Matrix & inT,  const Matrix & outTensor) : LinearForm(outTensor, false, true, 2), inTensor(inT), inGeometry(in)  { }

BimaterialInterface::~BimaterialInterface() { } ;

void BimaterialInterface::transform(const Function & x, const Function & y)
{
	xtransform = x ;
	ytransform = y ;
}

Matrix BimaterialInterface::getTensor(const Point * p) const
{
	VirtualMachine vm ;
	
	FunctionMatrix C(3,3) ;
	
	Function domain(inGeometry, xtransform, ytransform) ;
	
	for(size_t i = 0 ; i < 3 ; i++)
	{
		for(size_t j = 0 ; j < 3 ; j++)
		{
			C[i][j] = f_positivity(domain)*param[i][j] + f_negativity(domain)*inTensor[i][j];
		}
	}
	
	return vm.eval(C, p->x, p->y) ;
}

Matrix BimaterialInterface::apply(const Function & p_i, const Function & p_j, const IntegrableEntity *e) const
{
	VirtualMachine vm ;
	Function domain(inGeometry, xtransform, ytransform) ;
	FunctionMatrix C(3,3) ;
	std::vector<Variable> v ;
	v.push_back(XI);
	v.push_back(ETA);
	
	for(size_t i = 0 ; i < 3 ; i++)
	{
		for(size_t j = 0 ; j < 3 ; j++)
		{
			C[i][j] = f_positivity(domain)*param[i][j] + f_negativity(domain)*inTensor[i][j];
		}
	}
	
	return vm.ieval(Gradient(p_i) * C * Gradient(p_j, true), e,v) ;
}

Matrix BimaterialInterface::apply(const Function & p_i, const Function & p_j, const std::valarray< std::pair<Point,double> > &gp, const std::valarray<Matrix> &Jinv) const
{
	
	VirtualMachine vm ;
	
	Function domain(inGeometry) ;
	Function position = xtransform ;
	position = xtransform(ytransform) ;
	position = domain(position) ;
	
	FunctionMatrix C(3,3) ;
	std::vector<Variable> v ;
	v.push_back(XI);
	v.push_back(ETA);
	
	for(size_t i = 0 ; i < 3 ; i++)
	{
		for(size_t j = 0 ; j < 3 ; j++)
		{
			C[i][j] = f_positivity(position)*param[i][j] + f_negativity(position)*inTensor[i][j];
		}
	}
	
	return vm.ieval(Gradient(p_i) * C * Gradient(p_j, true), gp, Jinv,v) ;
}

bool BimaterialInterface::fractured() const
{
	return false;
}

Form * BimaterialInterface::getCopy() const 
{
	return new BimaterialInterface(*this) ;
}


Vector BimaterialInterface::getForces(const ElementState * s, const Function & p_i, const Function & p_j, const std::valarray< std::pair<Point, double> > &gp, const std::valarray<Matrix> &Jinv) const 
{
	return Vector(0) ;
}


RadialStiffnessGradient::RadialStiffnessGradient(double E_int, double nu_int, double rint, double E_ext, double nu_ext, double rext, Point c) : LinearForm(Matrix(3,3), false, true, 2), paramAlt(3,3), r_ext(rext), r_int(rint), dr(r_ext-r_int), centre(c)
{

	param[0][0] = E_int/(1-nu_int*nu_int) ; param[0][1] =E_int/(1-nu_int*nu_int)*nu_int ; param[0][2] = 0 ;
	param[1][0] = E_int/(1-nu_int*nu_int)*nu_int ; param[1][1] = E_int/(1-nu_int*nu_int) ; param[1][2] = 0 ; 
	param[2][0] = 0 ; param[2][1] = 0 ; param[2][2] = E_int/(1-nu_int*nu_int)*(1.-nu_int)/2. ; 
	
	paramAlt[0][0] = E_ext/(1-nu_ext*nu_ext) ;        paramAlt[0][1] = E_ext/(1-nu_ext*nu_ext)*nu_ext ; paramAlt[0][2] = 0 ;
	paramAlt[1][0] = E_ext/(1-nu_ext*nu_ext)*nu_ext ; paramAlt[1][1] = E_ext/(1-nu_ext*nu_ext) ;        paramAlt[1][2] = 0 ; 
	paramAlt[2][0] = 0 ; paramAlt[2][1] = 0 ;         paramAlt[2][2] = E_ext/(1-nu_ext*nu_ext)*(1.-nu_ext)/2. ; 
	this->space_d = true ;
}

RadialStiffnessGradient::~RadialStiffnessGradient() { } ;

void RadialStiffnessGradient::transform(const Function & x, const Function & y)
{
	r = f_sqrt(((x-centre.x)^2)+((y-centre.y)^2)) ;
}

Matrix RadialStiffnessGradient::apply(const Function & p_i, const Function & p_j, const IntegrableEntity *e) const
{
	VirtualMachine vm ;
	
	FunctionMatrix C(3,3) ;
	std::vector<Variable> v ;
	v.push_back(XI);
	v.push_back(ETA);

	for(size_t i = 0 ; i < 3 ; i++)
	{
		for(size_t j = 0 ; j < 3 ; j++)
		{
// 				double val =  (paramAlt[i][j] - param[i][j])/dr ;
			C[i][j] = ((r-r_int)/dr)*paramAlt[i][j] - ((r-r_ext)/dr)*param[i][j];
		}
	}
	
	return vm.ieval(Gradient(p_i) * C * Gradient(p_j, true), e,v) ;
}

bool RadialStiffnessGradient::fractured() const
{
	return false ;
}

Matrix RadialStiffnessGradient::getTensor(const Point * p) const
{
	VirtualMachine vm ;
	
	FunctionMatrix C(3,3) ;
	
	for(size_t i = 0 ; i < 3 ; i++)
	{
		for(size_t j = 0 ; j < 3 ; j++)
		{
// 				double val =  (paramAlt[i][j] - param[i][j])/dr ;
			C[i][j] = ((r-r_int)/dr)*paramAlt[i][j] - ((r-r_ext)/dr)*param[i][j];
		}
	}
	
	return vm.eval(C, p->x, p->y) ;
}

Matrix RadialStiffnessGradient::apply(const Function & p_i, const Function & p_j, const std::valarray< std::pair<Point,double> > &gp, const std::valarray<Matrix> &Jinv) const
{
	
	FunctionMatrix C(3,3) ;
	std::vector<Variable> v ;
	v.push_back(XI);
	v.push_back(ETA);

	for(size_t i = 0 ; i < 3 ; i++)
	{
		for(size_t j = 0 ; j < 3 ; j++)
		{
			C[i][j] = ((r-r_int)/dr)*paramAlt[i][j] - ((r-r_ext)/dr)*param[i][j];
		}
	}
	
	return VirtualMachine().ieval(Gradient(p_i) * C * Gradient(p_j, true), gp, Jinv,v) ;
}

Form * RadialStiffnessGradient::getCopy() const 
{
	return new RadialStiffnessGradient(*this) ;
}

Vector RadialStiffnessGradient::getForces(const ElementState * s, const Function & p_i, const Function & p_j, const std::valarray< std::pair<Point, double> > &gp, const std::valarray<Matrix> &Jinv) const 
{
	return Vector(0) ;
}

// Diffusion2D::Diffusion2D(double capacity, double conductivityX, double conductivityY) : LinearForm(Matrix(2,2), false, true, 1)
// {
// 	
// 	c = capacity ;
// 	
// 	alpha =  0.5;
// 	tau = .01;
// 	
// 	param[0][0] = conductivityX ; 
// 	param[1][1] = conductivityY ; 
// 
// }
// 	
// Diffusion2D::~Diffusion2D() { } ;
// 
// 
// Matrix Diffusion2D::apply(const Function & p_i, const Function & p_j, const IntegrableEntity *e) const
// {
// 	Matrix ret(1,1) ;
// 	
// 	ret[0][0] = VirtualMachine().ieval(VectorGradient(p_i) * param * VectorGradient(p_j, true), e)/(alpha*tau) + c*VirtualMachine().ieval(p_i * p_j, e);
// 	
// 	return ret ;
// }
// 	
// Matrix Diffusion2D::apply(const Function & p_i, const Function & p_j, const std::valarray< std::pair<Point,double> > &gp, const std::valarray<Matrix> &Jinv) const
// {
// 	Matrix ret(1,1) ;
// 	
// 	ret[0][0] =  VirtualMachine().ieval(VectorGradient(p_i) * param * VectorGradient(p_j, true) , gp, Jinv)/(alpha*tau) + c*VirtualMachine().ieval(p_i * p_j, gp) ;
// 	
// 	return ret ;
// }
// 	
// void Diffusion2D::step(double timestep, ElementState * currentState) 
// {
// 	Vector du = currentState->getBuffer() ;
// 	currentState->getDisplacements() +=  tau*currentState->getPreviousDisplacements()+ du;
// 	currentState->getPreviousDisplacements() += du/(alpha*tau) ;
// }
// 	
// bool Diffusion2D::fractured() const
// {
// 	return false;
// }
// 	
// 
// Form * Diffusion2D::getCopy() const 
// {
// 	return new Diffusion2D(*this) ;
// }
// 
// bool Diffusion2D::hasInducedForces() const 
// {
// 	return true ;
// }
// 
// Vector Diffusion2D::getForces(const ElementState * s, const Function & p_i, const Function & p_j, const std::valarray< std::pair<Point, double> > &gp, const std::valarray<Matrix> &Jinv) const 
// {
// 	
// 	Vector ret(1) ;
// 
// 	size_t index = 0 ;
// 	for(size_t i = 0 ; i < s->getParent()->getBoundingPoints().size() ; i++)
// 	{
// 		if(VirtualMachine().eval(p_j, s->getParent()->inLocalCoordinates( s->getParent()->getBoundingPoint(i))) > .5)
// 		{
// 			index = i ;
// 			break ;
// 		}
// 	}
// 	
// 	Vector v = s->getPreviousDisplacements() ;
// 	Vector u = s->getDisplacements() ;
// 	
// 	ret[0] = -c*VirtualMachine().ieval(p_i*p_j, gp)*(v[index]) -
// 			VirtualMachine().ieval(VectorGradient(p_i) * param * VectorGradient(p_j, true), gp, Jinv)*
// 						(u[index]+tau*v[index]);
// 	
// 	return ret ;
// }
	
NonLinearStiffness::NonLinearStiffness(Function f, double n) 
{
	E = f ;
	nu = n ;
	this->time_d = false ;
	this->type = NON_LINEAR ;
	this->state = NULL ;
	FunctionMatrix m0(3,3) ;
	VirtualMachine vm ;
	
	m0[0][0] = E/(1-nu*nu) ; m0[0][1] =E/(1-nu*nu)*nu ; m0[0][2] = Function() ;
	m0[1][0] = E/(1-nu*nu)*nu ; m0[1][1] = E/(1-nu*nu) ; m0[1][2] = Function() ; 
	m0[2][0] = Function() ; m0[2][1] = Function() ; m0[2][2] = E/(1-nu*nu)*(1.-nu)/2. ; 
	
	param = vm.eval(m0, Point(0.333333333, 0.3333333333)) ;
}

NonLinearStiffness::~NonLinearStiffness() { } ;

Matrix NonLinearStiffness::apply(const Function & p_i, const Function & p_j, const IntegrableEntity *e) const
{
	VirtualMachine vm ;
	std::valarray< std::pair<Point,double> > gp = e->getGaussPoints() ;
	std::valarray<Point> pts(gp.size()) ;
	for(size_t i = 0; i < gp.size() ; i++)
	{
		pts[i] = gp[i].first ;
	}
	
	double E_ = 0;
	Matrix m0(3,3) ;
	Vector displacements = state->getDisplacements( pts) ;
	for(size_t i = 0 ; i < gp.size() ; i++)
	{
		E_ += vm.eval(E, displacements[i*2],displacements[i*2+1])*gp[i].second ;
	}
	
	m0[0][0] = E_/(1-nu*nu) ; m0[0][1] =E_/(1-nu*nu)*nu ; m0[0][2] = 0 ;
	m0[1][0] = E_/(1-nu*nu)*nu ; m0[1][1] = E_/(1-nu*nu) ; m0[1][2] = 0 ; 
	m0[2][0] = 0 ; m0[2][1] = 0 ; m0[2][2] = E_/(1-nu*nu)*(1.-nu)/2. ; 
	
	std::vector<Variable> v ;
	v.push_back(XI);
	v.push_back(ETA);
	
	return vm.ieval(Gradient(p_i) * m0 * Gradient(p_j, true), e,v) ;
}

bool NonLinearStiffness::hasInducedForces() const
{
	return true ;
} ;

bool NonLinearStiffness::hasInducedMatrix() const 
{
	return true ;
}

Matrix NonLinearStiffness::apply(const Function & p_i, const Function & p_j, const std::valarray< std::pair<Point, double> > &gp, const std::valarray<Matrix> &Jinv) const
{
	
	VirtualMachine vm ;
	Matrix m0(3,3) ;
	
	std::valarray<Point> pts(gp.size()) ;
	for(size_t i = 0; i < gp.size() ; i++)
	{
		pts[i] = gp[i].first ;
	}
	
	Vector displacements = state->getDisplacements( pts) ;
	double E_ = 0;
	
	for(size_t i = 0 ; i < gp.size() ; i++)
	{
		E_ += vm.eval(E, displacements[i*2],displacements[i*2+1])*gp[i].second ;
	}
	
	m0[0][0] = E_/(1-nu*nu) ; m0[0][1] =E_/(1-nu*nu)*nu ; m0[0][2] = 0 ;
	m0[1][0] = E_/(1-nu*nu)*nu ; m0[1][1] = E_/(1-nu*nu) ; m0[1][2] = 0 ; 
	m0[2][0] = 0 ; m0[2][1] = 0 ; m0[2][2] = E_/(1-nu*nu)*(1.-nu)/2. ; 
	
	std::vector<Variable> v ;
	v.push_back(XI);
	v.push_back(ETA);
	
	return vm.ieval(Gradient(p_i) * m0 * Gradient(p_j, true), gp, Jinv,v) ;
}

Vector NonLinearStiffness::getForces(const ElementState * s, const Function & p_i, const Function & p_j, const std::valarray< std::pair<Point, double> > &gp, const std::valarray<Matrix> &Jinv) const 
{
	Vector stress = state->getStress(gp) ; 

	std::vector<Variable> v ;
	v.push_back(XI);
	v.push_back(ETA);
	
	return VirtualMachine().ieval(Gradient(p_j, true)*stress, gp, Jinv,v) ;
}

bool NonLinearStiffness::isActive() const 
{
	if(!state)
		return false ;
	
	std::valarray< std::pair<Point, double> > gp = state->getParent()->getGaussPoints() ;
	std::valarray<Point> pts(gp.size()) ;
	for(size_t i = 0; i < gp.size() ; i++)
	{
		pts[i] = gp[i].first ;
	}
	
	Vector displacements = state->getDisplacements( pts) ;
	double E_ = 0;
	VirtualMachine vm ;
	
	for(size_t i = 0 ; i < gp.size() ; i++)
	{
		E_ += vm.eval(E, displacements[i*2],displacements[i*2+1])*gp[i].second ;
	}
	
	return E_ > 1e-6 ;
}

Form * NonLinearStiffness::getCopy() const 
{
	return new NonLinearStiffness(*this) ;
}
	



TwoDCohesiveForces::TwoDCohesiveForces(const IntegrableEntity *s, const IntegrableEntity *t, const SegmentedLine * sl) 
{
	this->time_d = false ;
	this->type = NON_LINEAR ;
	this->target = t ;
	this->source = s ;
	active = false ;
	startArea = s->area() ;
	for(size_t i = 0 ; i < sl->getBoundingPoints().size()-1 ; i++)
	{
		Segment test(sl->getBoundingPoint(i), sl->getBoundingPoint(i+1)) ;
		normals.push_back(test.normal()) ;
	}
	
}
	
Form * TwoDCohesiveForces::getCopy() const
{
	return new TwoDCohesiveForces(*this) ;
}

TwoDCohesiveForces::~TwoDCohesiveForces() { } ;

Matrix TwoDCohesiveForces::apply(const Function & p_i, const Function & p_j, const IntegrableEntity *e) const
{
	return Matrix() ;
}
	
Matrix TwoDCohesiveForces::apply(const Function & p_i, const Function & p_j, const std::valarray< std::pair<Point, double> > &gp, const std::valarray<Matrix> &Jinv) const
{
	return Matrix() ;
}
	
bool TwoDCohesiveForces::hasInducedForces() const
{
	return true ;
}

bool TwoDCohesiveForces::hasInducedMatrix() const 
{
	return false ;
}
	
	
Vector TwoDCohesiveForces::getForces(const ElementState * s, const Function & p_i, const Function & p_j, const std::valarray< std::pair<Point, double> > &gp, const std::valarray<Matrix> &Jinv) const 
{

	Vector ret(0., 2) ;
	
	bool enrichedDof = !(p_j.getDofID() == -1) ;
	
	if(!enrichedDof)
		return ret ;
	

	Vector apparentStress = source->getState()->getNonEnrichedStress(gp,Jinv) ; 
	for(size_t i = 0 ; i < gp.size() ; i++)
	{
		if(enrichedDof)
		{
			Vector stress(3) ;
			for(size_t j = 0 ; j < 3 ; j++)
				stress[j] = apparentStress[i*3+j] ;
			std::vector<Variable> v ;
			v.push_back(XI);
			v.push_back(ETA);
			
			Matrix grad = VirtualMachine().geval(p_j, Jinv[i],v, gp[i].first, true) ;
			Vector force = (Vector)(grad*stress) ;
				
			double normalAmplitude = force[0]*normals[0].x + force[1]*normals[0].y;
			double tangeantAmplitude = -force[0]*normals[0].y + force[1]*normals[0].x;
				
			Vector normalForce(2) ;
			normalForce[0] = normals[0].x*normalAmplitude ;
			normalForce[1] = normals[0].y*normalAmplitude ;
				
			Vector tangeantForce(2) ;
			tangeantForce[0] = -normals[0].y*tangeantAmplitude ;
			tangeantForce[1] = normals[0].x*tangeantAmplitude ;

			ret += normalForce*gp[i].second ;
			ret += tangeantForce*gp[i].second ;
			
		}

	}
	
	return ret ;

}
	
void TwoDCohesiveForces::step(double timestep, ElementState * s) 
{
// 		this->state = currentState ;
	
	
	s->getPreviousPreviousDisplacements() = s->getPreviousDisplacements() ;
	s->getPreviousDisplacements() = s->getDisplacements() ;
	s->getPreviousPreviousEnrichedDisplacements() = s->getPreviousEnrichedDisplacements() ;
	s->getPreviousEnrichedDisplacements() = s->getEnrichedDisplacements() ;
	
	size_t ndofs = s->getParent()->getBehaviour()->getNumberOfDegreesOfFreedom() ;
	int offset = ndofs-1 ;
	
	if(s->getEnrichedDisplacements().size() != s->getParent()->getEnrichmentFunctions().size()*ndofs)
		s->getEnrichedDisplacements().resize(s->getParent()->getEnrichmentFunctions().size()*ndofs) ;
	
	for(size_t i = 0 ; i < s->getParent()->getBoundingPoints().size() ; i++)
	{
		s->getDisplacements()[i*ndofs] = s->getBuffer()[i*ndofs] ;
		s->getDisplacements()[i*ndofs+offset] = s->getBuffer()[i*ndofs+offset] ;
	}
	
	for(size_t i = 0 ; i < s->getParent()->getEnrichmentFunctions().size() ; i++)
	{
		std::cout << i*ndofs << " vs " << s->getEnrichedDisplacements().size() << std::endl ;
		s->getEnrichedDisplacements()[i*ndofs] = s->getBuffer()[(i+s->getParent()->getBoundingPoints().size())*ndofs] ;
		s->getEnrichedDisplacements()[i*ndofs+offset] = s->getBuffer()[(i+s->getParent()->getBoundingPoints().size())*ndofs+offset] ;
	}
}
	

bool TwoDCohesiveForces::isActive() const 
{
	Point a = source->getBoundingPoint(0) + source->getState()->getDisplacements(source->getBoundingPoint(0), false);
	Point b = source->getBoundingPoint(source->getBoundingPoints().size()/3) + source->getState()->getDisplacements(source->getBoundingPoint(source->getBoundingPoints().size()/3), false);
	Point c = source->getBoundingPoint(2*source->getBoundingPoints().size()/3) + source->getState()->getDisplacements(source->getBoundingPoint(2*source->getBoundingPoints().size()/3), false);
	
	double newArea = Triangle(a,b,c).area() ;
	
	if(newArea< startArea)
		return true ;
	
	return false ;
	
}


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
		
	this->time_d = true;
	
}
	
ViscoElasticity::~ViscoElasticity() { } ;

Matrix ViscoElasticity::apply(const Function & p_i, const Function & p_j, const IntegrableEntity *e) const
{
	std::vector<Variable> v ;
	v.push_back(XI);
	v.push_back(ETA);
	v.push_back(ZETA);
	return VirtualMachine().ieval(Gradient(p_i) * param * Gradient(p_j, true), e,v) ;
}

Matrix ViscoElasticity::apply(const Function & p_i, const Function & p_j, const std::valarray< std::pair<Point,double> > &gp, const std::valarray<Matrix> &Jinv) const
{
	std::vector<Variable> v ;
	v.push_back(XI);
	v.push_back(ETA);
	v.push_back(ZETA);
	return VirtualMachine().ieval(Gradient(p_i) * param * Gradient(p_j, true), gp, Jinv,v) ;
}
void ViscoElasticity::step(double timestep, ElementState * currentState)
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
	Vector Deltastrain = currentState->getDeltaStrain(currentState->getParent()->getCenter()) ;
	Vector previousStrain = currentState->getStrain(currentState->getParent()->getCenter())-Deltastrain ;
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
	
Vector ViscoElasticity::getForces(const ElementState * s, const Function & p_i, const Function & p_j, const std::valarray< std::pair<Point, double> > &gp, const std::valarray<Matrix> &Jinv) const
{
	std::vector<Variable> v ;
	v.push_back(XI);
	v.push_back(ETA);
	v.push_back(ZETA);
	return VirtualMachine().ieval(Gradient(p_j, true)*average_delta_sigma, gp, Jinv,v) ;
}

bool ViscoElasticity::hasInducedForces()
{
	return true ;
}

Form * ViscoElasticity::getCopy() const 
{
	return new ViscoElasticity(*this) ;
}



