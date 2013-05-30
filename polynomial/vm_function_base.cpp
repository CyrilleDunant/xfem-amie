//
// C++ Implementation: vm_function_base
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2013
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "vm_function_base.h"
#include "../elements/elements.h"
#include <string.h>

using namespace Mu ;


void concatenateFunctions(const Function & src0, const Function & src1, Function & dst)
{
	Function tmpdst;

	tmpdst.hasGeoOp = src0.hasGeoOp || src1.hasGeoOp || dst.hasGeoOp;
	if(tmpdst.hasGeoOp )
		tmpdst.geo_op.resize(HEAP_SIZE, (GeometryOperation*)nullptr);

	tmpdst.adress_a = src0.adress_a ;
	tmpdst.byteCode = src0.byteCode ;
	tmpdst.values = src0.values ;

	for(size_t i = 0 ; i < src1.values.size() ; i++)
		tmpdst.values.push_back(src1.values[i]) ;
	
	if(src0.hasGeoOp)
	{
		for(size_t i = 0 ; i < src0.byteCode.size() ; i++)
		{
			if(src0.geo_op[i])
				tmpdst.geo_op[i] = src0.geo_op[i]->getCopy() ;
		}
	}
	
 	tmpdst.adress_t.resize( src0.adress_t.size() + src1.adress_t.size(), 0) ;
 	for(size_t i = 0 ; i < tmpdst.adress_t.size() ; i++)
 		tmpdst.adress_t[i] = HEAP_VARIABLE_TRANSFORM_OFFSET + i ;

	if(src0.adress_t.size() != src0.transformed.size() )
	{
		std::cout << "la???" << std::endl ;
//		VirtualMachine().print(src0) ;
	}
	
	for(size_t i = 0 ; i < src0.adress_t.size() ; i++)
	{
		Function f = src0.transform(i) ;
		tmpdst.setVariableTransform( src0.transformed[i], f, false);
	}

	for(size_t i = 0 ; i < src1.adress_t.size() ; i++)
	{
		Function f = src1.transform(i) ;
		tmpdst.setVariableTransform( src1.transformed[i], f , false);
	}
	
	for(size_t i = 0 ; i < src1.byteCode.size() ; i++)
	{
		tmpdst.byteCode.push_back(src1.byteCode[i])  ;
		tmpdst.adress_a.push_back(0);
		tmpdst.adress_a.push_back(0);
		tmpdst.adress_a.push_back(0);
		tmpdst.adress_a.push_back(0);
		
		if(src1.hasGeoOp)
			if(src1.geo_op[i])
				tmpdst.geo_op[i+src0.byteCode.size()] = src1.geo_op[i]->getCopy() ;
			
		if(src1.adress_a[i*4] >= HEAP_SIZE-src1.values.size()-1)
		{
			tmpdst.adress_a[(i+src0.byteCode.size())*4] = src1.adress_a[i*4]-src0.values.size() ;
		}
		
		if(src1.adress_a[i*4+1] >= HEAP_SIZE-src1.values.size()-1)
		{
			tmpdst.adress_a[(i+src0.byteCode.size())*4+1] = src1.adress_a[i*4+1]-src0.values.size() ;
		}
		
		if(src1.adress_a[i*4+2] >= HEAP_SIZE-src1.values.size()-1)
		{
			tmpdst.adress_a[(i+src0.byteCode.size())*4+2] = src1.adress_a[i*4+2]-src0.values.size() ;
		}
		

		if(src1.adress_a[i*4] >= 8 && src1.adress_a[i*4] < HEAP_SIZE-src1.values.size())
			tmpdst.adress_a[(i+src0.byteCode.size())*4] = src1.adress_a[i*4]+1 ;
		if(src1.adress_a[i*4+1] >= 8 && src1.adress_a[i*4+1] < HEAP_SIZE-src1.values.size())
			tmpdst.adress_a[(i+src0.byteCode.size())*4+1] = src1.adress_a[i*4+1]+1 ;
		if(src1.adress_a[i*4+2] >= 8 && src1.adress_a[i*4+2] < HEAP_SIZE-src1.values.size())
			tmpdst.adress_a[(i+src0.byteCode.size())*4+2] = src1.adress_a[i*4+2]+1 ;
		
		if(src1.adress_a[i*4] >= HEAP_VARIABLE_TRANSFORM_OFFSET && src1.adress_a[i*4] < HEAP_SIZE - src1.values.size()-1)
		{
			tmpdst.adress_a[ (i+src0.byteCode.size())*4 ] = src1.adress_a[i*4]+src0.adress_t.size() ;
		}
		if(src1.adress_a[i*4+1] >= HEAP_VARIABLE_TRANSFORM_OFFSET && src1.adress_a[i*4+1] < HEAP_SIZE - src1.values.size()-1)
		{
			tmpdst.adress_a[ (i+src0.byteCode.size())*4+1 ] = src1.adress_a[i*4+1]+src0.adress_t.size() ;
		}
		if(src1.adress_a[i*4+2] >= HEAP_VARIABLE_TRANSFORM_OFFSET && src1.adress_a[i*4+2] < HEAP_SIZE - src1.values.size()-1)
		{
			tmpdst.adress_a[ (i+src0.byteCode.size())*4+2 ] = src1.adress_a[i*4+2]+src0.adress_t.size() ;
		}
		
		if(src1.adress_a[i*4] < 8 )
			tmpdst.adress_a[(i+src0.byteCode.size())*4] = src1.adress_a[i*4] ;
		if(src1.adress_a[i*4+1] < 8 )
			tmpdst.adress_a[(i+src0.byteCode.size())*4+1] = src1.adress_a[i*4+1] ;
		if(src1.adress_a[i*4+2] < 8 )
			tmpdst.adress_a[(i+src0.byteCode.size())*4+2] = src1.adress_a[i*4+2] ;
	}
	
	dst = tmpdst ;
}

void concatenateFunctions(const Function & src0_, const Function & src1_, const Function & src2_, Function & dst)
{
	Function src0(src0_) ;
	Function src1(src1_) ;
	Function src2(src2_) ;
	dst.adress_a = src0.adress_a ;
	dst.byteCode = src0.byteCode ;
	dst.hasGeoOp = src0.hasGeoOp || src1.hasGeoOp || src2.hasGeoOp  || dst.hasGeoOp;
	if(dst.hasGeoOp )
		dst.geo_op.resize(HEAP_SIZE, (GeometryOperation*)nullptr);
	
	for(size_t i = 0 ; i < src0.values.size() ; i++)
		dst.values.push_back(src0.values[i]) ;
	for(size_t i = 0 ; i < src1.values.size() ; i++)
		dst.values.push_back( src1.values[i]) ;
	for(size_t i = 0 ; i < src2.values.size() ; i++)
		dst.values.push_back( src2.values[i]) ;
	
	if(src0.hasGeoOp)
		for(size_t i = 0 ; i < src0.byteCode.size() ; i++)
		{
			if(src0.geo_op[i])
				dst.geo_op[i] = src0.geo_op[i]->getCopy() ;
		}

 	dst.adress_t.resize( src0.adress_t.size() + src1.adress_t.size() + src2.adress_t.size(), 0) ;
 	for(size_t i = 0 ; i < dst.adress_t.size() ; i++)
 		dst.adress_t[i] = HEAP_VARIABLE_TRANSFORM_OFFSET + i ;

	for(size_t i = 0 ; i < src0.adress_t.size() ; i++)
	{
		Function f = src0.transform(i) ;
		dst.setVariableTransform( src0.transformed[0], f, false);
	}

	for(size_t i = 0 ; i < src1.adress_t.size() ; i++)
	{
		Function f = src1.transform(i) ;
		dst.setVariableTransform( src1.transformed[i], f , false);
	}
		
	for(size_t i = 0 ; i < src2.adress_t.size() ; i++)
	{
		Function f = src2.transform(i) ;
		dst.setVariableTransform( src2.transformed[i], f , false);
	}
		
	for(size_t i = 0 ; i < src1.byteCode.size() ; i++)
	{
		dst.byteCode.push_back(src1.byteCode[i]) ;
		dst.adress_a.push_back(0);
		dst.adress_a.push_back(0);
		dst.adress_a.push_back(0);
		dst.adress_a.push_back(0);
		
		if(src1.hasGeoOp)
			if(src1.geo_op[i])
				dst.geo_op[i+src0.byteCode.size()] = src1.geo_op[i]->getCopy() ;
			
		if(src1.adress_a[i*4] >= HEAP_SIZE-src1.values.size())
		{
			dst.adress_a[(i+src0.byteCode.size())*4] = src1.adress_a[i*4]-src0.values.size() ;
		}
		
		if(src1.adress_a[i*4+1] >= HEAP_SIZE-src1.values.size())
		{
			dst.adress_a[(i+src0.byteCode.size())*4+1] = src1.adress_a[i*4+1]-src0.values.size() ;
		}
		
		if(src1.adress_a[i*4+2] >= HEAP_SIZE-src1.values.size())
		{
			dst.adress_a[(i+src0.byteCode.size())*4+2] = src1.adress_a[i*4+2]-src0.values.size() ;
		}
		
		if(src1.adress_a[i*4] >= 8 && src1.adress_a[i*4] < HEAP_SIZE-src1.values.size())
			dst.adress_a[(i+src0.byteCode.size())*4] = src1.adress_a[i*4]+1 ;
		
		if(src1.adress_a[i*4+1] >= 8 && src1.adress_a[i*4+1] < HEAP_SIZE-src1.values.size())
			dst.adress_a[(i+src0.byteCode.size())*4+1] = src1.adress_a[i*4+1]+1 ;
		
		if(src1.adress_a[i*4+2] >= 8 && src1.adress_a[i*4+2] < HEAP_SIZE-src1.values.size())
			dst.adress_a[(i+src0.byteCode.size())*4+2] = src1.adress_a[i*4+2]+1 ;

		if(src1.adress_a[i*4] >= HEAP_VARIABLE_TRANSFORM_OFFSET && src1.adress_a[i*4] < HEAP_SIZE - src1.values.size()-1)
		{
			dst.adress_a[ (i+src0.byteCode.size())*4 ] = src1.adress_a[i*4]+src0.adress_t.size() ;
		}
		if(src1.adress_a[i*4+1] >= HEAP_VARIABLE_TRANSFORM_OFFSET && src1.adress_a[i*4+1] < HEAP_SIZE - src1.values.size()-1)
		{
			dst.adress_a[ (i+src0.byteCode.size())*4+1 ] = src1.adress_a[i*4+1]+src0.adress_t.size() ;
		}
		if(src1.adress_a[i*4+2] >= HEAP_VARIABLE_TRANSFORM_OFFSET && src1.adress_a[i*4+2] < HEAP_SIZE - src1.values.size()-1)
		{
			dst.adress_a[ (i+src0.byteCode.size())*4+2 ] = src1.adress_a[i*4+2]+src0.adress_t.size() ;
		}
		
		if(src1.adress_a[i*4] < 8 )
			dst.adress_a[(i+src0.byteCode.size())*4] = src1.adress_a[i*4] ;
		
		if(src1.adress_a[i*4+1] < 8 )
			dst.adress_a[(i+src0.byteCode.size())*4+1] = src1.adress_a[i*4+1] ;
		
		if(src1.adress_a[i*4+2] < 8 )
			dst.adress_a[(i+src0.byteCode.size())*4+2] = src1.adress_a[i*4+2] ;
	}
	
	for(size_t i = 0 ; i < src2.byteCode.size() ; i++)
	{
		dst.byteCode.push_back(src2.byteCode[i]) ;
		dst.adress_a.push_back(0);
		dst.adress_a.push_back(0);
		dst.adress_a.push_back(0);
		dst.adress_a.push_back(0);
		
		if(src2.hasGeoOp)
			if(src2.geo_op[i])
				dst.geo_op[i+src0.byteCode.size()+src1.byteCode.size()] = src2.geo_op[i]->getCopy() ;
			
		if(src2.adress_a[i*4] >= HEAP_SIZE-src2.values.size())
		{
			dst.adress_a[(i+src0.byteCode.size()+src1.byteCode.size())*4] = src2.adress_a[i*4]-src0.values.size()-src1.values.size() ;
		}
		
		if(src2.adress_a[i*4+1] >= HEAP_SIZE-src2.values.size())
		{
			dst.adress_a[(i+src0.byteCode.size()+src1.byteCode.size())*4+1] = src2.adress_a[i*4+1]-src0.values.size()-src1.values.size() ;
		}
		
		if(src2.adress_a[i*4+2] >= HEAP_SIZE-src2.values.size())
		{
			dst.adress_a[(i+src0.byteCode.size()+src1.byteCode.size())*4+2] = src2.adress_a[i*4+2]-src0.values.size()-src1.values.size() ;
		}
		
		if(src2.adress_a[i*4] >= 8 && src2.adress_a[i*4] < HEAP_SIZE-src2.values.size())
			dst.adress_a[(i+src0.byteCode.size()+src1.byteCode.size())*4] = src2.adress_a[i*4]+2 ;
		
		if(src2.adress_a[i*4+1] >= 8 && src2.adress_a[i*4+1] < HEAP_SIZE-src2.values.size())
			dst.adress_a[(i+src0.byteCode.size()+src1.byteCode.size())*4+1] = src2.adress_a[i*4+1]+2 ;
		
		if(src2.adress_a[i*4+2] >= 8 && src2.adress_a[i*4+2] < HEAP_SIZE-src2.values.size())
			dst.adress_a[(i+src0.byteCode.size()+src1.byteCode.size())*4+2] = src2.adress_a[i*4+2]+2 ;

		if(src2.adress_a[i*4] >= HEAP_VARIABLE_TRANSFORM_OFFSET && src2.adress_a[i*4] < HEAP_SIZE - src2.values.size()-1)
		{
			dst.adress_a[ (i+src0.byteCode.size()+src1.byteCode.size())*4 ] = src2.adress_a[i*4]+src0.adress_t.size()+src1.adress_t.size() ;
		}
		if(src2.adress_a[i*4+1] >= HEAP_VARIABLE_TRANSFORM_OFFSET && src2.adress_a[i*4+1] < HEAP_SIZE - src2.values.size()-1)
		{
			dst.adress_a[ (i+src0.byteCode.size()+src1.byteCode.size())*4+1 ] = src2.adress_a[i*4+1]+src0.adress_t.size()+src1.adress_t.size() ;
		}
		if(src2.adress_a[i*4+2] >= HEAP_VARIABLE_TRANSFORM_OFFSET && src2.adress_a[i*4+2] < HEAP_SIZE - src2.values.size()-1)
		{
			dst.adress_a[ (i+src0.byteCode.size()+src1.byteCode.size())*4+2 ] = src1.adress_a[i*4+2]+src0.adress_t.size()+src1.adress_t.size() ;
		}
		
		if(src2.adress_a[i*4] < 8 )
			dst.adress_a[(i+src0.byteCode.size()+src1.byteCode.size())*4] = src2.adress_a[i*4] ;
		
		if(src2.adress_a[i*4+1] < 8 )
			dst.adress_a[(i+src0.byteCode.size()+src1.byteCode.size())*4+1] = src2.adress_a[i*4+1] ;
		
		if(src2.adress_a[i*4+2] < 8 )
			dst.adress_a[(i+src0.byteCode.size()+src1.byteCode.size())*4+2] = src2.adress_a[i*4+2] ;
	}

  
}


void  Function::preCalculate(const GaussPointArray & gp , std::vector<Variable> & var, const double eps)
{
	//this needs to be in to stages, otherwise memory gets accessed too early.
	Vector * newVal = new Vector(VirtualMachine().eval(*this, gp)) ;
	if(precalc.find(gp.id) != precalc.end())
		delete precalc[gp.id] ;
	precalc[gp.id] = newVal ;
	std::map<Variable, Vector *> val ;
	for(size_t i = 0 ; i < var.size() ; i++)
	{
		if(val.find(var[i]) == val.end())
			val[var[i]] = new Vector(VirtualMachine().deval(*this,var[i] ,gp, eps)) ;
	}

	if(dprecalc.find(gp.id) != dprecalc.end())
	{
		for(auto i = dprecalc[gp.id].begin() ;  i != dprecalc[gp.id].end() ; ++i)
			delete i->second ;
	}

	dprecalc[gp.id] = val ;
}

void  Function::preCalculate(const GaussPointArray & gp )
{
	if(precalc.find(gp.id) != precalc.end())
		delete precalc[gp.id] ;
	precalc[gp.id] = new Vector(VirtualMachine().eval(*this, gp)) ;
}

const Vector & Function::getPrecalculatedValue(const GaussPointArray &gp ) const
{
	return *(precalc.find(gp.id)->second) ;
}

const Vector & Function::getPrecalculatedValue(const GaussPointArray &gp, Variable v) const
{
	return *dprecalc.find(gp.id)->second.find(v)->second ;
}

bool Function::precalculated(const GaussPointArray & gp) const
{
	return (precalc.find(gp.id) != precalc.end()) ;
}

bool Function::precalculated(const GaussPointArray & gp, Variable v) const
{
	auto precalculatedDerivative = dprecalc.find(gp.id) ;
	return precalculatedDerivative != dprecalc.end() && precalculatedDerivative->second.find(v) != precalculatedDerivative->second.end() ;
}


GtM Gradient::operator*( const Matrix & f) const
{
	return GtM(*this, f) ;
}

GDDtM GradientDotDot::operator*( const Matrix & f) const
{
	return GDDtM(*this, f) ;
}

GDDtML GradientDotDot::operator*( const std::vector<Matrix> & f) const
{
	return GDDtML(*this, f) ;
}

GDDtMtG GDDtM::operator*( const Gradient & f) const
{
	return GDDtMtG(*this, f) ;
}

GDDtMLtG GDDtML::operator*( const Gradient & f) const
{
	return GDDtMLtG(*this, f) ;
}

GtML Gradient::operator*( const std::vector<Matrix> & f) const
{
	return GtML(*this, f) ;
}

GDtML GradientDot::operator*( const std::vector<Matrix> & f) const
{
	return GDtML(*this, f) ;
}

GDtMLtG GDtML::operator*( const Gradient & f) const
{
	return GDtMLtG(*this, f) ;
}

GDtMLtGD GDtML::operator*( const GradientDot & f) const
{
	return GDtMLtGD(*this, f) ;
}


GtV Gradient::operator*( const Vector & f) const
{
	return GtV(*this, f) ;
}

GtVL Gradient::operator*( const std::vector<Vector> & f) const
{
	return GtVL(*this, f) ;
}

GDtVL GradientDot::operator*( const std::vector<Vector> & f) const
{
	return GDtVL(*this, f) ;
}

GtMtG GtM::operator*(const Gradient & f) const
{
	return GtMtG(this->first, this->second, f) ;
}

GtMLtG GtML::operator*(const Gradient & f) const
{
	return GtMLtG(this->first, this->second, f) ;
}

GtMLtGD GtML::operator*(const GradientDot & f) const
{
	return GtMLtGD(*this, f) ;
}

GtMtGD GtM::operator*(const GradientDot & f) const
{
	return GtMtGD(*this, f) ;
}

DtF Differential::operator *(const Function & f) const
{
	return DtF(*this, f) ;
}

DtD Differential::operator *(const Differential & f) const
{
	return DtD(*this, f) ;
}

DtV Differential::operator *(const Vector & f) const
{
	return DtV(*this, f) ;
}

DtVL Differential::operator *(const std::vector<Vector> & f) const
{
	return DtVL(*this, f) ;
}

DtGtMtG Differential::operator *(const GtMtG & g) const
{
	return DtGtMtG(*this, g) ;
}

GDtMtGD GDtM::operator*(const Mu::GradientDot & f) const
{
	return GDtMtGD(*this, f) ;
}

GDtMtG GDtM::operator*(const Mu::Gradient & f) const
{
	return GDtMtG(*this, f) ;
}

GDtM GradientDot::operator *(const Matrix & f) const
{
	return GDtM(*this, f) ;
}

GDtV GradientDot::operator *(const Vector & v) const
{
	return GDtV(*this, v) ;
}

VGtM VectorGradient::operator *(const Matrix & f) const 
{ 
	return VGtM(*this, f) ; 
}

VGtV VectorGradient::operator *(const Vector & f) const 
{ 
	return VGtV(*this, f) ; 
}


VGtMtVG VGtM::operator*(const Mu::VectorGradient & f) const 
{
	return VGtMtVG(first,second, f) ;
}


Function::Function() : derivative(nullptr),
		transforms(nullptr),
		e_diff(false),
		geo_op((GeometryOperation *)nullptr,0),
		dofID(-1),
		ptID (nullptr),
		hasGeoOp(false)
{
	byteCode.push_back(TOKEN_OPERATION_CONSTANT),
	values.push_back(0) ;
	adress_a.push_back(0);
	adress_a.push_back(0);
	adress_a.push_back(0);
	adress_a.push_back(0);
	initialiseAdresses();
}

Function & Function::operator=(const Function &f)
{
	if(derivative)
		for(size_t i = 0 ; i < derivative->size() ; i++)
			delete (*derivative)[i] ;
		
	delete derivative ;

	if(transforms)
		for(size_t i = 0 ; i < transforms->size() ; i++)
			delete (*transforms)[i] ;
		
	delete transforms ;
	adress_t.clear() ;
	transformed.clear() ;
	
	for(auto i = precalc.begin() ; i != precalc.end() ; ++i)
		delete i->second ;
	for(auto i = dprecalc.begin() ; i != dprecalc.end() ; ++i)
	{
		for(auto j = i->second.begin() ; j != i->second.end() ; ++j)
		{
			delete j->second ;
		}
	}
		
	derivative = nullptr ;
	e_diff = false ;
	if(f.derivative)
	{
		derivative = new std::valarray<Function *>((Function *)nullptr, f.derivative->size()) ;
		e_diff = true ;
		for(size_t i = 0 ; i < f.derivative->size() ; i++)
		{
			if((*f.derivative)[i])
				(*derivative)[i] = new Function(*(*f.derivative)[i]) ;
			else
				(*derivative)[i] = nullptr ;
		}
	}

	transforms = nullptr ;
	if(f.transforms)
	{
		transforms = new std::vector<Function *>() ;
		for(size_t i = 0 ; i < f.transforms->size() ; i++)
		{
			transforms->push_back(new Function(*(*f.transforms)[i])) ;
			transformed.push_back(f.transformed[i]);
			adress_t.push_back(f.adress_t[i]);
		}
	}

	if(hasGeoOp)
	{
		for(size_t i = 0 ; i < byteCode.size() ; i++)
		{
			delete geo_op[i] ;
		}
	}
	
	if(f.hasGeoOp)
	{
		if(!hasGeoOp)
			geo_op.resize(HEAP_SIZE,(GeometryOperation *)nullptr);
		for(size_t i = 0 ; i < f.byteCode.size() ; i++)
		{
			if(f.geo_op[i])
				geo_op[i] = f.geo_op[i]->getCopy() ;
		}
	}
	
	ptID = f.ptID ;
	dofID = f.dofID ;
	iPoint = f.iPoint ;
	byteCode = f.byteCode ;
	values = f.values ;
	adress_a = f.adress_a ;
	hasGeoOp = f.hasGeoOp ;

// 	initialiseAdresses();

	for(auto i = f.precalc.begin() ; i != f.precalc.end() ; ++i)
	{
		precalc[i->first] = new Vector(*i->second) ;
	}
	for(auto i = f.dprecalc.begin() ; i != f.dprecalc.end() ; ++i)
	{
		
		dprecalc[i->first] =  std::map<Variable, Vector *>() ;
		for(auto j = i->second.begin() ; j != i->second.end() ; ++j)
		{
			dprecalc[i->first][j->first] = new Vector(*j->second) ;
		}
	}
	return *this ;
}

int Function::getDofID() const
{
	return this->dofID ;
}
void Function::setDofID(size_t id)
{
	this->dofID = id ;
}

Point * Function::getPoint() const
{
	return this->ptID ;
}
void Function::setPoint(Point * id)
{
	this->ptID = id ;
}

boost::tuple<TokenOperationType, double, std::string> Function::toToken(const std::string & str, int iter, std::vector<double> & val) const
{
	if(isNumeral(str[0]) || (str[0] == '-' && isNumeral(str[1])))
	{
		return boost::make_tuple(TOKEN_OPERATION_CONSTANT, atof(str.c_str()), "") ;
	}
	else if(str == std::string("cst"))
	{
		if(iter < val.size())
			return boost::make_tuple(TOKEN_OPERATION_CONSTANT, val[iter], "cst") ; 
		return boost::make_tuple(TOKEN_OPERATION_CONSTANT, val[val.size()-1], "cst") ;		
	}
	else if(str == std::string("sin"))
	{
		return boost::make_tuple(TOKEN_OPERATION_SIN, 0., "") ;
	}
	else if(str == std::string("cos"))
	{
		return boost::make_tuple(TOKEN_OPERATION_COS, 0., "") ;
	}
	else if(str == std::string("tan"))
	{
		return boost::make_tuple(TOKEN_OPERATION_TAN, 0., "") ;
	}
	else if(str == std::string("sinh"))
	{
		return boost::make_tuple(TOKEN_OPERATION_SINH, 0., "") ;
	}
	else if(str == std::string("cosh"))
	{
		return boost::make_tuple(TOKEN_OPERATION_COSH, 0., "") ;
	}
	else if(str == std::string("tanh"))
	{
		return boost::make_tuple(TOKEN_OPERATION_TANH, 0., "") ;
	}
	else if(str == std::string("exp"))
	{
		return boost::make_tuple(TOKEN_OPERATION_EXP, 0., "") ;
	}
	else if(str == std::string("abs"))
	{
		return boost::make_tuple(TOKEN_OPERATION_ABS, 0., "") ;
	}
	else if(str == std::string("log"))
	{
		return boost::make_tuple(TOKEN_OPERATION_LOG, 0., "") ;
	}
	else if(str == std::string("sqrt"))
	{
		return boost::make_tuple(TOKEN_OPERATION_SQRT, 0., "") ;
	}
	else if(str == std::string("sign"))
	{
		return boost::make_tuple(TOKEN_OPERATION_SIGN, 0., "") ;
	}
	else if(str == std::string("atan2"))
	{
		return boost::make_tuple(TOKEN_OPERATION_ATAN2, 0., "") ;
	}
	else if(str == std::string("sign"))
	{
		return boost::make_tuple(TOKEN_OPERATION_SIGN, 0., "") ;
	}
	else if(str == std::string("positive"))
	{
		return boost::make_tuple(TOKEN_OPERATION_POSITIVITY, 0., "") ;
	}
	else if(str == std::string("negative"))
	{
		return boost::make_tuple(TOKEN_OPERATION_NEGATIVITY, 0., "") ;
	}
	else if(str == std::string("+"))
	{
		return boost::make_tuple(TOKEN_OPERATION_PLUS, 0., "") ;
	}
	else if(str == std::string("-"))
	{
		return boost::make_tuple(TOKEN_OPERATION_MINUS, 0., "") ;
	}
	else if(str == std::string("*"))
	{
		return boost::make_tuple(TOKEN_OPERATION_TIMES, 0., "") ;
	}
	else if(str == std::string("/"))
	{
		return boost::make_tuple(TOKEN_OPERATION_DIVIDES, 0., "") ;
	}
	else if(str == std::string("^"))
	{
		return boost::make_tuple(TOKEN_OPERATION_POWER, 0., "") ;
	}
	else if(str == std::string("x"))
	{
		return boost::make_tuple(TOKEN_OPERATION_X, 0., "") ;
	}
	else if(str == std::string("y"))
	{
		return boost::make_tuple(TOKEN_OPERATION_Y, 0., "") ;
	}
	else if(str == std::string("z"))
	{
		return boost::make_tuple(TOKEN_OPERATION_Z, 0., "") ;
	}
	else if(str == std::string("t"))
	{
		return boost::make_tuple(TOKEN_OPERATION_T, 0., "") ;
	}
	else if(str == std::string("u"))
	{
		return boost::make_tuple(TOKEN_OPERATION_U, 0., "") ;
	}
	else if(str == std::string("v"))
	{
		return boost::make_tuple(TOKEN_OPERATION_V, 0., "") ;
	}
	else 
	{
		return boost::make_tuple(TOKEN_OPERATION_W, 0., "") ;
	}
}
	
bool Function::isOperator(const char c) const
{
	if(
		c == '+' ||
		c == '-' ||
		c == '/' ||
		c == '*' ||
		c == '^'
		)
		return true ;
	return false ;
} ;
	
bool Function::isNumeral(const char c) const
{
	if(
		c == '0' ||
		c == '1' ||
		c == '2' ||
		c == '3' ||
		c == '4' ||
		c == '5' ||
		c == '6' ||
		c == '7' ||
		c == '8' ||
		c == '9' 
		)
		return true ;
	return false ;
} ;
	
bool Function::isSeparator(const char c) const
{
	if(
		c == ' ' ||
		isOperator(c)
		)
		return true ;
	return false ;
} ;

std::pair<size_t, boost::tuple< TokenOperationType, double, std::string>> Function::getNext(size_t init, const std::string & form,  int cstIterator, std::vector<double> & val)
{
	return getNext(init, form.c_str(), cstIterator, val) ;
}

std::pair<size_t, boost::tuple<TokenOperationType, double, std::string>> Function::getNext(size_t init, const char * form, int cstIterator, std::vector<double> & val)
{
	std::string cur_token("") ;
	char cur_char  = form[init] ;
	while(isSeparator(cur_char) && !isOperator(cur_char))
	{
		init++ ;
		cur_char  = form[init] ;
	}
	
	if(isOperator(cur_char))
	{
		if(cur_char != '-')
		{
			cur_token+=cur_char ;
			init++ ;
			return std::make_pair(init, toToken(cur_token, cstIterator, val)) ;
		}
		else if (isNumeral(form[init+1]))
		{
			cur_token+=cur_char ;
			init++ ;
			cur_char  = form[init] ;
			while(!isSeparator(cur_char))
			{
				if(init < strlen(form))
				{
					cur_token+=cur_char ;
					init++ ;
					cur_char  = form[init] ;
				}
				else
					break ;
			}
			return std::make_pair(init, toToken(cur_token, cstIterator, val)) ;
		}
		else
		{
			cur_token+=cur_char ;
			init++ ;
			return std::make_pair(init, toToken(cur_token, cstIterator, val)) ;
		}
	}
	while(!isSeparator(cur_char))
	{
		if(init < strlen(form))
		{
			cur_token+=cur_char ;
			init++ ;
			cur_char  = form[init] ;
		}
		else
			break ;
	}
	return std::make_pair(init, toToken(cur_token, cstIterator, val)) ;
}

void Function::initialiseAdresses(size_t offset)
{
	unsigned short int da = 8 + (offset>0);
	int counter = offset ;
	for(size_t i = offset ; i < byteCode.size() ; i++)
	{
		switch(byteCode[i])
		{
			case TOKEN_OPERATION_CONSTANT:
			{
				adress_a[4*counter+2] = da ;
				adress_a[4*counter++] = da++ ;
				break ;
			}
			case TOKEN_OPERATION_X:
			{
				adress_a[4*counter+2] = da ;
				adress_a[4*counter++] = da++ ;
				break ;
			}
			case TOKEN_OPERATION_Y:
			{
				adress_a[4*counter+2] = da ;
				adress_a[4*counter++] = da++ ;
				break ;
			}
			case TOKEN_OPERATION_Z:
			{
				adress_a[4*counter+2] = da ;
				adress_a[4*counter++] = da++ ;
				break ;
			}
			case TOKEN_OPERATION_T:
			{
				adress_a[4*counter+2] = da ;
				adress_a[4*counter++] = da++ ;
				break ;
			}
			case TOKEN_OPERATION_U:
			{
				adress_a[4*counter+2] = da ;
				adress_a[4*counter++] = da++ ;
				break ;
			}
			case TOKEN_OPERATION_V:
			{
				adress_a[4*counter+2] = da ;
				adress_a[4*counter++] = da++ ;
				break ;
			}
			case TOKEN_OPERATION_W:
			{
				adress_a[4*counter+2] = da ;
				adress_a[4*counter++] = da++ ;
				break ;
			}
			case TOKEN_OPERATION_COS:
			{
				adress_a[4*counter+2] = da-1 ;
				adress_a[4*counter++] = da-1 ;
				break ;
			}
			case TOKEN_OPERATION_ABS:
			{
				adress_a[4*counter+2] = da-1 ;
				adress_a[4*counter++] = da-1 ;
				break ;
			}
			case TOKEN_OPERATION_TAN:
			{
				adress_a[4*counter+2] = da-1 ;
				adress_a[4*counter++] = da-1 ;
				break ;
			}
			case TOKEN_OPERATION_SIN:
			{
				adress_a[4*counter+2] = da-1 ;
				adress_a[4*counter++] = da-1 ;
				break ;
			}
			case TOKEN_OPERATION_EXP:
			{
				adress_a[4*counter+2] = da-1 ;
				adress_a[4*counter++] = da-1 ;
				break ;
			}
			case TOKEN_OPERATION_SIGN:
			{
				adress_a[4*counter+2] = da-1 ;
				adress_a[4*counter++] = da-1 ;
				break ;
			}
			case TOKEN_OPERATION_POSITIVITY:
			{
				adress_a[4*counter+2] = da-1 ;
				adress_a[4*counter++] = da-1 ;
				break ;
			}
			case TOKEN_OPERATION_NEGATIVITY:
			{
				adress_a[4*counter+2] = da-1 ;
				adress_a[4*counter++] = da-1 ;
				break ;
			}
			case TOKEN_OPERATION_LOG:
			{
				adress_a[4*counter+2] = da-1 ;
				adress_a[4*counter++] = da-1 ;
				break ;
			}
			case TOKEN_OPERATION_COSH:
			{
				adress_a[4*counter+2] = da-1 ;
				adress_a[4*counter++] = da-1 ;
				break ;
			}
			case TOKEN_OPERATION_SINH:
			{
				adress_a[4*counter+2] = da-1 ;
				adress_a[4*counter++] = da -1;
				break ;
			}
			case TOKEN_OPERATION_TANH:
			{
				adress_a[4*counter+2] = da-1 ;
				adress_a[4*counter++] = da-1 ;
				break ;
			}
			case TOKEN_OPERATION_SQRT:
			{
				adress_a[4*counter+2] = da-1 ;
				adress_a[4*counter++] = da-1 ;
				break ;
			}
			case TOKEN_OPERATION_BESSEL:
			{
				adress_a[counter+2] = da-1 ;
				adress_a[counter++] = da-1 ;
				break ;
			}
			case TOKEN_OPERATION_ATAN2:
			{
				adress_a[4*counter] = da-2 ;
				adress_a[4*counter+1] = da-1 ;
				adress_a[4*counter++] = da-2 ;
				--da ;
				break ;
			}
			case TOKEN_OPERATION_INTERPOLATE:
			{
				adress_a[4*counter+2] = da-2 ;
				adress_a[4*counter+1] = da-1 ;
				adress_a[4*counter++] = da-2 ;
				--da ;
				break ;
			}
			case TOKEN_OPERATION_PLUS:
			{
				adress_a[4*counter+2] = da-2 ;
				adress_a[4*counter+1] = da-1 ;
				adress_a[4*counter++] = da-2 ;
				--da ;
				break ;
			}
			case TOKEN_OPERATION_MINUS:
			{
				adress_a[4*counter+2] = da-2 ;
				adress_a[4*counter+1] = da-1 ;
				adress_a[4*counter++] = da-2 ;
				--da ;
				break ;
			}
			case TOKEN_OPERATION_TIMES:
			{
				adress_a[4*counter+2] = da-2 ;
				adress_a[4*counter+1] = da-1 ;
				adress_a[4*counter++] = da-2 ;
				--da ;
				break ;
			}
			case TOKEN_OPERATION_DIVIDES:
			{
				adress_a[4*counter+2] = da-2 ;
				adress_a[4*counter+1] = da-1 ;
				adress_a[4*counter++] = da-2 ;
				--da ;
				break ;
			}
			case TOKEN_OPERATION_POWER:
			{
				adress_a[4*counter+2] = da-2 ;
				adress_a[4*counter+1] = da-1 ;
				adress_a[4*counter++] = da-2 ;
				--da ;
				break ;
			}
			case TOKEN_OPERATION_GEO_OPERATION:
			{
				int offset = geo_op[i]->adressOffset() ;
				if(offset == -1) 
				{
					adress_a[4*counter+2] = da-2 ;
					adress_a[4*counter+1] = da-1 ;
					adress_a[4*counter++] = da-2 ;
					--da ;
				}
				
				if(offset == -2)
				{
					adress_a[4*counter+2] = da-3 ;
					adress_a[4*counter+1] = da-2 ;
					adress_a[4*counter++] = da-1 ;
					da -= 2 ;
				}
				if(offset == 0) 
				{
					adress_a[4*counter+2] = da-3 ;
					adress_a[4*counter+1] = da-2 ;
					adress_a[4*counter++] = da-1 ;
				}
			}
			default:
			break ;
		}
	}
	
	std::vector<TokenOperationType> newbyteCode ;
	std::valarray<GeometryOperation *> newgeo_op(geo_op) ;
	std::vector<double> newvalues ;
	std::vector<short unsigned int> newadress_a(adress_a) ;
	
	size_t constIterator = 0 ;
	for(size_t i = offset ; i < byteCode.size() ; i++)
	{
		switch(byteCode[i])
		{
			case TOKEN_OPERATION_X:
			{
				bool found = false ;
				for(size_t j = i+1 ; j < byteCode.size() ; j++)
				{
					if(adress_a[i*4] == adress_a[j*4])
					{
						adress_a[j*4] = 1 ;
						found = true ;
						break ;
					}
					if(adress_a[i*4] == adress_a[j*4+1])
					{
						adress_a[j*4+1] = 1 ;
						found = true ;
						break ;
					}
					if(adress_a[i*4] == adress_a[j*4+2])
					{
						adress_a[j*4+2] = 1 ;
						found = true ;
						break ;
					}
				}
				if(!found)
				{
					adress_a[i*4] = 1 ;
				}
				break ;
			}
			case TOKEN_OPERATION_Y:
			{
				bool found = false ;
				for(size_t j = i+1 ; j < byteCode.size() ; j++)
				{
					if(adress_a[i*4] == adress_a[j*4])
					{
						adress_a[j*4] = 2 ;
						found = true ;
						break ;
					}
					if(adress_a[i*4] == adress_a[j*4+1])
					{
						adress_a[j*4+1] = 2 ;
						found = true ;
						break ;
					}
					if(adress_a[i*4] == adress_a[j*4+2])
					{
						adress_a[j*4+2] = 2 ;
						found = true ;
						break ;
					}
				}
				if(!found)
				{
					adress_a[i*4] = 2 ;
				}
				break ;
			}
			case TOKEN_OPERATION_Z:
			{
				bool found = false ;
				for(size_t j = i+1 ; j < byteCode.size() ; j++)
				{
					if(adress_a[i*4] == adress_a[j*4])
					{
						adress_a[j*4] = 3 ;
						found = true ;
						break ;
					}
					if(adress_a[i*4] == adress_a[j*4+1])
					{
						adress_a[j*4+1] = 3 ;
						found = true ;
						break ;
					}
					if(adress_a[i*4] == adress_a[j*4+2])
					{
						adress_a[j*4+2] = 3 ;
						found = true ;
						break ;
					}
				}
				if(!found)
				{
					adress_a[i*4] = 3 ;
				}
				break ;
			}
			case TOKEN_OPERATION_T:
			{
				bool found = false ;
				for(size_t j = i+1 ; j < byteCode.size() ; j++)
				{
					if(adress_a[i*4] == adress_a[j*4])
					{
						adress_a[j*4] = 4 ;
						found = true ;
						break ;
					}
					if(adress_a[i*4] == adress_a[j*4+1])
					{
						adress_a[j*4+1] = 4 ;
						found = true ;
						break ;
					}
					if(adress_a[i*4] == adress_a[j*4+2])
					{
						adress_a[j*4+2] = 4 ;
						found = true ;
						break ;
					}
				}
				if(!found)
				{
					adress_a[i*4] = 4 ;
				}
				break ;
			}
			case TOKEN_OPERATION_U:
			{
				bool found = false ;
				for(size_t j = i+1 ; j < byteCode.size() ; j++)
				{
					if(adress_a[i*4] == adress_a[j*4])
					{
						adress_a[j*4] = 5 ;
						found = true ;
						break ;
					}
					if(adress_a[i*4] == adress_a[j*4+1])
					{
						adress_a[j*4+1] = 5 ;
						found = true ;
						break ;
					}
					if(adress_a[i*4] == adress_a[j*4+2])
					{
						adress_a[j*4+2] = 5 ;
						found = true ;
						break ;
					}
				}
				if(!found)
				{
					adress_a[i*4] = 5 ;
				}
				break ;
			}
			case TOKEN_OPERATION_V:
			{
				bool found = false ;
				for(size_t j = i+1 ; j < byteCode.size() ; j++)
				{
				if(adress_a[i*4] == adress_a[j*4])
					{
						adress_a[j*4] = 6 ;
						found = true ;
						break ;
					}
					if(adress_a[i*4] == adress_a[j*4+1])
					{
						adress_a[j*4+1] = 6 ;
						found = true ;
						break ;
					}
					if(adress_a[i*4] == adress_a[j*4+2])
					{
						adress_a[j*4+2] = 6 ;
						found = true ;
						break ;
					}
				}
				if(!found)
				{
					adress_a[i*4] = 6 ;
				}
				break ;
			}
			case TOKEN_OPERATION_W:
			{
				bool found = false ;
				for(size_t j = i+1 ; j < byteCode.size() ; j++)
				{
					if(adress_a[i*4] == adress_a[j*4])
					{
						adress_a[j*4] = 7 ;
						found = true ;
						break ;
					}
					if(adress_a[i*4] == adress_a[j*4+1])
					{
						adress_a[j*4+1] = 7 ;
						found = true ;
						break ;
					}
					if(adress_a[i*4] == adress_a[j*4+2])
					{
						adress_a[j*4+2] = 7 ;
						found = true ;
						break ;
					}
				}
				if(!found)
				{
					adress_a[i*4] = 7 ;
				}
				
				break ;
			}
			case TOKEN_OPERATION_CONSTANT:
			{
				bool found = false ;
				for(size_t j = i+1 ; j < byteCode.size() ; j++)
				{
					if(adress_a[i*4] == adress_a[j*4])
					{
						newvalues.push_back(values[constIterator++]) ;
						adress_a[j*4] = HEAP_SIZE-newvalues.size() ;
						found = true ;
						break ;
					}
					else if(adress_a[i*4] == adress_a[j*4+1])
					{
						newvalues.push_back(values[constIterator++]) ;
						adress_a[j*4+1] = HEAP_SIZE-newvalues.size() ;
						found = true ;
						break ;
					}
					else if(adress_a[i*4] == adress_a[j*4+2])
					{
						newvalues.push_back(values[constIterator++]) ;
						adress_a[j*4+2] = HEAP_SIZE-newvalues.size() ;
						found = true ;
						break ;
					}
				}
				if(!found)
				{
					newvalues.push_back(values[constIterator++]) ;
					adress_a[i*4] = HEAP_SIZE-newvalues.size() ;
				}
				
				break ;
			}
		}
	}
	
	counter = offset ;
	if(byteCode.size() > 1)
	{
		for(size_t i = offset ; i < byteCode.size() ; i++)
		{
			if(byteCode[i] == TOKEN_OPERATION_X ||
				byteCode[i] == TOKEN_OPERATION_Y ||
				byteCode[i] == TOKEN_OPERATION_Z ||
				byteCode[i] == TOKEN_OPERATION_T ||
				byteCode[i] == TOKEN_OPERATION_W ||
				byteCode[i] == TOKEN_OPERATION_U ||
				byteCode[i] == TOKEN_OPERATION_V ||
				byteCode[i] == TOKEN_OPERATION_CONSTANT 
			)
			{
				continue ;
			}
			if(adress_a[i*4] == adress_a[i*4+2] && i < byteCode.size()-1 && i < newbyteCode.size())
			{
				if(byteCode[i] == TOKEN_OPERATION_PLUS)
					byteCode[i] = TOKEN_OPERATION_INPLACE_PLUS ;
				if(byteCode[i] == TOKEN_OPERATION_MINUS)
					byteCode[i] = TOKEN_OPERATION_INPLACE_MINUS ;
				if(byteCode[i] == TOKEN_OPERATION_TIMES)
					byteCode[i] = TOKEN_OPERATION_INPLACE_TIMES ;
				if(byteCode[i] == TOKEN_OPERATION_DIVIDES)
					byteCode[i] = TOKEN_OPERATION_INPLACE_DIVIDES ;
				if(byteCode[i] == TOKEN_OPERATION_POWER)
					byteCode[i] = TOKEN_OPERATION_INPLACE_POWER ;
			}

			newbyteCode.push_back(byteCode[i]); 
			newadress_a.push_back(0);newadress_a.push_back(0);newadress_a.push_back(0);newadress_a.push_back(0);
			if(hasGeoOp)
				newgeo_op[counter]=geo_op[i]; 
			newadress_a[counter*4]= adress_a[i*4];
			newadress_a[counter*4+1]= adress_a[i*4+1];
			newadress_a[counter*4+2]= adress_a[i*4+2];
			counter++ ;

		}
		byteCode = newbyteCode; 
		geo_op= newgeo_op; 
		values = newvalues ;
		adress_a= newadress_a;
	}
}

Function::Function(const char *f): derivative(nullptr),
		transforms(nullptr),
		e_diff(false),
		geo_op((GeometryOperation *)nullptr,0),
		dofID(-1),
		ptID (nullptr),
		hasGeoOp(false)
{

	size_t init = 0 ;
	std::vector<double> val ;
	int iter = 0 ;
	while(init < strlen(f))
	{
		std::pair<size_t, boost::tuple<TokenOperationType, double, std::string> > temp = getNext(init, f, iter, val) ;
		byteCode.push_back(boost::get<0>(temp.second)) ;
		adress_a.push_back(0);
		adress_a.push_back(0);
		adress_a.push_back(0);
		adress_a.push_back(0);
		if(byteCode.back() == TOKEN_OPERATION_CONSTANT)
			values.push_back(boost::get<1>(temp.second));
		
		

		init = temp.first ;
	}
	initialiseAdresses();
}

Function::Function(const char *f, double v): derivative(nullptr),
		transforms(nullptr),
		e_diff(false),
		geo_op((GeometryOperation *)nullptr,0),
		dofID(-1),
		ptID (nullptr),
		hasGeoOp(false)
{

	size_t init = 0 ;
	std::vector<double> val ;
	val.push_back(v);
	while(init < strlen(f))
	{
		std::pair<size_t, boost::tuple<TokenOperationType, double, std::string> > temp = getNext(init, f, 0, val) ;
		byteCode.push_back(boost::get<0>(temp.second)) ;
		adress_a.push_back(0);
		adress_a.push_back(0);
		adress_a.push_back(0);
		adress_a.push_back(0);
		if(byteCode.back() == TOKEN_OPERATION_CONSTANT)
			values.push_back(boost::get<1>(temp.second));

		init = temp.first ;
	}
	initialiseAdresses();
}

Function::Function(const char *f, std::vector<double> & val): derivative(nullptr),
		transforms(nullptr),
		e_diff(false),
		geo_op((GeometryOperation *)nullptr,0),
		dofID(-1),
		ptID (nullptr),
		hasGeoOp(false)
{

	size_t init = 0 ;
	while(init < strlen(f))
	{
		std::pair<size_t, boost::tuple<TokenOperationType, double, std::string> > temp = getNext(init, f, 0, val) ;
		byteCode.push_back(boost::get<0>(temp.second)) ;
		adress_a.push_back(0);
		adress_a.push_back(0);
		adress_a.push_back(0);
		adress_a.push_back(0);
		if(byteCode.back() == TOKEN_OPERATION_CONSTANT)
			values.push_back(boost::get<1>(temp.second));

		init = temp.first ;
	}
	initialiseAdresses();
}

Function::Function(const std::string &f) : derivative(nullptr),
		transforms(nullptr),
		e_diff(false),
		geo_op((GeometryOperation *)nullptr,0),
		dofID(-1),
		ptID (nullptr),
		hasGeoOp(false)
{
	size_t init = 0 ;
	std::vector<double> val ;
	while(init < f.length())
	{
		std::pair<size_t, boost::tuple<TokenOperationType, double, std::string> > temp = getNext(init, f,0, val) ;
		byteCode.push_back(boost::get<0>(temp.second)) ;
		adress_a.push_back(0);
		adress_a.push_back(0);
		adress_a.push_back(0);
		adress_a.push_back(0);
		if(byteCode.back() == TOKEN_OPERATION_CONSTANT)
			values.push_back(boost::get<1>(temp.second));
	}
	initialiseAdresses();
}

Function::Function(const std::string &f, double v) : derivative(nullptr),
		transforms(nullptr),
		e_diff(false),
		geo_op((GeometryOperation *)nullptr,0),
		dofID(-1),
		ptID (nullptr),
		hasGeoOp(false)
{
	size_t init = 0 ;
	std::vector<double> val ;
	val.push_back(v);
	while(init < f.length())
	{
		std::pair<size_t, boost::tuple<TokenOperationType, double, std::string> > temp = getNext(init, f,0, val) ;
		byteCode.push_back(boost::get<0>(temp.second)) ;
		adress_a.push_back(0);
		adress_a.push_back(0);
		adress_a.push_back(0);
		adress_a.push_back(0);
		if(byteCode.back() == TOKEN_OPERATION_CONSTANT)
			values.push_back(boost::get<1>(temp.second));
	}
	initialiseAdresses();
}

Function::Function(const std::string &f, std::vector<double> & val) : derivative(nullptr),
		transforms(nullptr),
		e_diff(false),
		geo_op((GeometryOperation *)nullptr,0),
		dofID(-1),
		ptID (nullptr),
		hasGeoOp(false)
{
	size_t init = 0 ;
	while(init < f.length())
	{
		std::pair<size_t, boost::tuple<TokenOperationType, double, std::string> > temp = getNext(init, f,0, val) ;
		byteCode.push_back(boost::get<0>(temp.second)) ;
		adress_a.push_back(0);
		adress_a.push_back(0);
		adress_a.push_back(0);
		adress_a.push_back(0);
		if(byteCode.back() == TOKEN_OPERATION_CONSTANT)
			values.push_back(boost::get<1>(temp.second));
	}
	initialiseAdresses();
}

Function::Function(const Line & l, ElementarySurface * s) : derivative(nullptr),
		transforms(nullptr),
		e_diff(false),
		geo_op((GeometryOperation *)nullptr,HEAP_SIZE),
		dofID(-1),
		ptID (nullptr),
		hasGeoOp(true)
{
	Function g = s->getXTransform() ;
	Function f = s->getYTransform() ;
	concatenateFunctions(g, f, *this);
	hasGeoOp = true ;
	byteCode.push_back(TOKEN_OPERATION_GEO_OPERATION);
	adress_a.push_back(0);
	adress_a.push_back(0);
	adress_a.push_back(0);
	adress_a.push_back(0);
	adress_a[(byteCode.size()-1)*4+2] = 8 ;
	adress_a[(byteCode.size()-1)*4+1] = 9 ;
	adress_a[(byteCode.size()-1)*4] = 8 ;
	geo_op[byteCode.size()-1] = new LineDistanceOperation(l) ;
	initialiseAdresses();
}

Function::Function(const Point & l,  ElementarySurface * s) : derivative(nullptr),
		transforms(nullptr),
		e_diff(false),
		geo_op((GeometryOperation *)nullptr,HEAP_SIZE),
		dofID(-1),
		ptID (nullptr)
{
	Function g = s->getXTransform() ;
	Function f = s->getYTransform() ;
	concatenateFunctions(g, f, *this);
	byteCode.push_back(TOKEN_OPERATION_GEO_OPERATION);
	adress_a.push_back(0);
	adress_a.push_back(0);
	adress_a.push_back(0);
	adress_a.push_back(0);
	adress_a[(byteCode.size()-1)*4+2] = 8 ;
	adress_a[(byteCode.size()-1)*4+1] = 9 ;
	adress_a[(byteCode.size()-1)*4] = 8 ;
	hasGeoOp = true ;
	geo_op[byteCode.size()-1] = new PointDistanceBinaryOperation(l) ;
}

Function::Function(const Point & l,  ElementaryVolume * s) : derivative(nullptr),
		transforms(nullptr),
		e_diff(false),
		geo_op((GeometryOperation *)nullptr,HEAP_SIZE),
		dofID(-1),
		ptID (nullptr)
{
	Function g = s->getXTransform() ;
	Function f = s->getYTransform() ;
	Function h = s->getZTransform() ;
	concatenateFunctions(g, f, h,  *this);
	hasGeoOp = true ;
	byteCode.push_back(TOKEN_OPERATION_GEO_OPERATION);
	adress_a.push_back(0);
	adress_a.push_back(0);
	adress_a.push_back(0);
	adress_a.push_back(0);
	adress_a[(byteCode.size()-1)*4+2] = 8 ;
	adress_a[(byteCode.size()-1)*4+1] = 9 ;
	adress_a[(byteCode.size()-1)*4] = 10 ;
	
	geo_op[byteCode.size()-1] = new PointDistanceTrinaryOperation(l) ;
}

Function::Function(double a,  ElementarySurface * s) : derivative(nullptr),
		transforms(nullptr),
		e_diff(false),
		geo_op((GeometryOperation *)nullptr,HEAP_SIZE),
		adress_a(HEAP_SIZE*4),
		dofID(-1),
		ptID (nullptr),
		hasGeoOp(true)
{
	byteCode.push_back(TOKEN_OPERATION_X) ;
	byteCode.push_back(TOKEN_OPERATION_Y) ;
	byteCode.push_back(TOKEN_OPERATION_GEO_OPERATION) ;
	adress_a.push_back(0); adress_a.push_back(0); adress_a.push_back(0); adress_a.push_back(0);
	adress_a.push_back(0); adress_a.push_back(0); adress_a.push_back(0); adress_a.push_back(0);
	adress_a.push_back(0); adress_a.push_back(0); adress_a.push_back(0); adress_a.push_back(0);
	geo_op[byteCode.size()-1] = new RotationBinaryOperation(a) ;
}


Function::Function(double a,const Point & p,   ElementarySurface * s): derivative(nullptr),
		transforms(nullptr),
		e_diff(false),
		geo_op((GeometryOperation *)nullptr,HEAP_SIZE),
		dofID(-1),
		ptID (nullptr),
		hasGeoOp(true)
{
	Function g = s->getXTransform() ;
	Function f = s->getYTransform() ;
	concatenateFunctions(g, f, *this);
	byteCode.push_back(TOKEN_OPERATION_GEO_OPERATION); 
	adress_a.push_back(0); adress_a.push_back(0); adress_a.push_back(0); adress_a.push_back(0);
	adress_a[(byteCode.size()-1)*4+2] = 8 ;
	adress_a[(byteCode.size()-1)*4+1] = 9 ;
	adress_a[(byteCode.size()-1)*4] = 8 ;
	
	geo_op[byteCode.size()-1] = new AngleBinaryOperation(a,p) ;
}

Function::Function( const Geometry * geo, const ElementarySurface * s) : derivative(nullptr),
		transforms(nullptr),
		e_diff(false),
		geo_op((GeometryOperation *)nullptr,HEAP_SIZE),
		dofID(-1),
		ptID (nullptr),
		hasGeoOp(true)
{
	Function g = s->getXTransform() ;
	Function f = s->getYTransform() ;
	concatenateFunctions(g, f, *this);
	byteCode.push_back(TOKEN_OPERATION_GEO_OPERATION);
	adress_a.push_back(0); adress_a.push_back(0); adress_a.push_back(0); adress_a.push_back(0);
	adress_a[(byteCode.size()-1)*4+2] = 8 ;
	adress_a[(byteCode.size()-1)*4+1] = 9 ;
	adress_a[(byteCode.size()-1)*4] = 8 ;
	
	geo_op[byteCode.size()-1] = new DomainBinaryOperation(geo) ;
}

Function f_project(const Geometry *g, const Function &x, const Function &y)
{
	Function ret ;
	if(x.hasGeoOp || y.hasGeoOp)
		ret.geo_op.resize(HEAP_SIZE,(GeometryOperation *)nullptr);		
	ret.byteCode = y.byteCode;
	ret.values = y.values;
	if(y.hasGeoOp)
		for(size_t i = 0 ; i < y.byteCode.size() ; i++)
		{
			if(y.geo_op[i])
				ret.geo_op[i] = y.geo_op[i]->getCopy() ;
		}
	for(size_t i = 0 ; i < x.values.size() ; i++)
		ret.values.push_back( x.values[i]);
	
	for(size_t i = 0 ; i < x.byteCode.size() ; i++)
	{
		ret.byteCode.push_back(x.byteCode[i]);

		if(x.hasGeoOp)
		{
			if(x.geo_op[i])
				ret.geo_op[i+y.byteCode.size()] = x.geo_op[i]->getCopy() ;
		}
	}
	ret.byteCode.push_back(TOKEN_OPERATION_GEO_OPERATION); 
	ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0);
	ret.geo_op[ret.byteCode.size()-1] = new ProjectionBinaryOperation(g) ;
	ret.hasGeoOp = true ;
	ret.initialiseAdresses();
	return ret ;
}

Function::Function(const std::vector<Segment> s , ElementarySurface * u, PositionTokenType t) : derivative(nullptr),
		transforms(nullptr),
		e_diff(false),
		geo_op((GeometryOperation *)nullptr,HEAP_SIZE),
		dofID(-1),
		ptID (nullptr),
		hasGeoOp(true)
{

	switch(t)
	{
	case POSITION_TOKEN :
		{
			byteCode.push_back(TOKEN_OPERATION_X); 
			byteCode.push_back(TOKEN_OPERATION_Y) ;
			byteCode.push_back(TOKEN_OPERATION_GEO_OPERATION) ;
			adress_a.push_back(0); adress_a.push_back(0); adress_a.push_back(0); adress_a.push_back(0);
			adress_a.push_back(0); adress_a.push_back(0); adress_a.push_back(0); adress_a.push_back(0);
			adress_a.push_back(0); adress_a.push_back(0); adress_a.push_back(0); adress_a.push_back(0);
			geo_op[byteCode.size()-1] = new PositionOperation(s) ;

			
			break ;
		}
	case PROJECTION_TOKEN :
		{

			for(size_t i = 0 ; i < s.size() ; i++)
			{
				byteCode.push_back(TOKEN_OPERATION_GEO_OPERATION) ;
				adress_a.push_back(0); adress_a.push_back(0); adress_a.push_back(0); adress_a.push_back(0);
				geo_op[byteCode.size()-1] = new ProjectionOperation2D(s[i]) ;
				byteCode.push_back(TOKEN_OPERATION_PLUS)  ;
				adress_a.push_back(0); adress_a.push_back(0); adress_a.push_back(0); adress_a.push_back(0);
			}
				
			
			break ;
		}
	}
	initialiseAdresses();
}

Function::Function(const Function &f) : e_diff(f.e_diff), byteCode(f.byteCode),transforms(nullptr),
 geo_op((GeometryOperation*)nullptr,f.geo_op.size()),values(f.values),adress_a(f.adress_a),ptID(f.ptID),hasGeoOp(f.hasGeoOp),dofID(f.dofID)
{
	if(f.derivative)
	{
		derivative = new std::valarray<Function *>((Function *)nullptr,f.derivative->size()) ;
		for(size_t i = 0 ; i < f.derivative->size() ; i++)
		{
			if((*f.derivative)[i])
				(*derivative)[i] = new Function(*(*f.derivative)[i]) ;
			else
				(*derivative)[i] = nullptr ;
		}
		e_diff = true ;
	}
	else
	{
		derivative = nullptr ;
	}
	
	if(f.transforms)
	{
		transforms = new std::vector<Function *>() ;
		for(size_t i = 0 ; i < f.transforms->size() ; i++)
		{
			transforms->push_back(new Function(*(*f.transforms)[i]));
			adress_t.push_back(f.adress_t[i]);
			transformed.push_back(f.transformed[i]);
		}
	}
		

	if(hasGeoOp)
	{
		for(size_t i = 0 ; i < f.byteCode.size() ; i++)
		{
			if(f.geo_op[i])
				geo_op[i] = f.geo_op[i]->getCopy() ;
		}
	}
	
	for(auto i = f.precalc.begin() ; i != f.precalc.end() ; ++i)
		precalc[i->first] = new Vector(*i->second) ;
	for(auto i = f.dprecalc.begin() ; i != f.dprecalc.end() ; ++i)
	{
		dprecalc[i->first] =  std::map<Variable, Vector *>() ;
		for(auto j = i->second.begin() ; j != i->second.end() ; ++j)
		{
			dprecalc[i->first][j->first] = new Vector(*j->second) ;
		}
	}
	
}


Function::~Function()
{
	
	if(hasGeoOp)
	{
		for(size_t i = 0 ; i < byteCode.size() ; i++)
		{
			delete geo_op[i] ;
		}
	}
	
	if(derivative)
	{
		for(size_t i = 0 ; i < derivative->size() ; i++)
		{
			delete (*derivative)[i] ;
		}
		delete derivative ;
	}

	for(auto i = precalc.begin() ; i != precalc.end() ; ++i)
		delete i->second ;
	for(auto i = dprecalc.begin() ; i != dprecalc.end() ; ++i)
	{
		for(auto j = i->second.begin() ; j != i->second.end() ; ++j)
		{
			delete j->second ;
		}
	}
}


bool Function::isDifferentiable() const 
{
	return e_diff && derivative;
}

bool Function::isDifferentiable(const Variable v) const 
{
	if(derivative && derivative->size() > v) 
	{
		return (*derivative)[v] != nullptr ;
	}
	return false ;
}

bool Function::isDifferentiable(size_t v) const 
{
	if(derivative && derivative->size() > v) 
	{
		return (*derivative)[v] != nullptr ;
	}
	return false ;
}

int Function::getNumberOfDerivatives() const 
{
	if(derivative)
	{
		return derivative->size() ;
	}
	return 0 ;
}

void Function::setVariableTransform( const Variable v, Function & f, bool replaceToken) 
{
#warning you have to manually call makeVariableTransformDerivative() once all variable transforms have been set
  
	if( derivativeTransformed )
	{
		  std::cout << "try to add variable transform to derived function - do nothing instead" << std::endl ;
		  return ;
	}
    
	if(!transforms)
		transforms = new std::vector<Function *>() ;
	
	transforms->push_back(new Function(f));
	transformed.push_back(v) ;
	
	if(replaceToken)
	{
		adress_t.push_back( HEAP_VARIABLE_TRANSFORM_OFFSET+adress_t.size() );
		size_t last = adress_t.size()-1 ;
		for(size_t i = 0 ; i < byteCode.size() ; i++)
		{
//			std::cout << adress_a[i*4] << "\t" << adress_a[i*4+1] << "\t" << adress_a[i*4+2] << "\t" << adress_a[i*4+3] << "\n" ;
			if(adress_a[ i*4 ] == (size_t) v+1)
				  adress_a[ i*4 ] = adress_t[last] ;
			if(adress_a[ i*4 +1 ] == (size_t) v+1)
				  adress_a[ i*4 +1 ] = adress_t[last] ;
			if(adress_a[ i*4 +2 ] == (size_t) v+1)
				  adress_a[ i*4 +2 ] = adress_t[last] ;
			if(adress_a[ i*4 +3 ] == (size_t) v+1)
				  adress_a[ i*4 +3 ] = adress_t[last] ;
		}
	}
}

void Function::makeVariableTransformDerivative() 
{
	derivativeTransformed = true ;
  
	if(!derivative)
	{
		return ;
	}
	
	for(size_t j = 0 ; j <derivative->size() ; j++)
	{
		for(size_t i = 0 ; i < transforms->size() ; i++)
		{
			if(! (*transforms)[i]->isDifferentiable(j) )
			{
				for(size_t k = 0 ; k < derivative->size() ; k++)
					delete (*derivative)[k] ;
				delete derivative ;
				derivative = nullptr ;
				e_diff = false ;
//				std::cout << "no replacement derivative - out" << std::endl ;
				return ;
			}
		}
	}
	
	std::vector<Function *> original ;
	for(size_t i = 0 ; i < derivative->size() ; i++)
		original.push_back((*derivative)[i]);
  
	std::vector<Function *> replacements ;
	for(size_t i = 0 ; i < original.size() ; i++)
		replacements.push_back( nullptr );
	
	for(size_t i = 0 ; i < transforms->size() ; i++)
	{
		if( !replacements[ transformed[i] ] )
			replacements[ transformed[i] ] = (*transforms)[i] ;
	}
	
	for(size_t i = 0 ; i < derivative->size() ; i++)
	{
// 		std::cout << "  << " << i << " >>" << std::endl ;
		Function newderivative ;
		bool initialized = false ;
		for(size_t j = 0 ; j < replacements.size() ; j++)
		{
// 			std::cout << "  [[ " << j << " ]]" << std::endl ;
			if(replacements[j] ) 
			{
				Function oldderivative = *(original[j]) ;
				for(size_t k = 0 ; k < replacements.size() ; k++)
				{
					if(replacements[k])
						oldderivative.setVariableTransform( (const Variable) k, *(replacements[k]) );
				}
				oldderivative.makeVariableTransformDerivative() ;
				if(oldderivative.isDifferentiable(TIME_VARIABLE))
					  std::cout << "b" ;
				oldderivative *= replacements[j]->d( (const Variable) i ) ;
				if(! initialized )
				{
					newderivative = oldderivative ;
					initialized = true ;
				}
				else
				{
// 					std::cout << newderivative.adress_t.size() << " --- " << newderivative.transformed.size() << std::endl ;
					newderivative += oldderivative ;
				}
			}
			else if( i == j )
			{
				if(! initialized )
				{
					newderivative = *(original[i]) ;
					initialized = true ;
				}
				else
					newderivative += *(original[i]) ;
			}
		}
		
		(*derivative)[i] = new Function(newderivative) ;

	  
	}
	
	for(size_t i = 0 ; i < original.size() ; i++)
		delete original[i] ;
  
}

size_t Function::derivationDepth() const
{
	size_t depth = 0 ;
	if(isDifferentiable())
	{
		bool init = false ;
		for(size_t i = 0 ; i < derivative->size() ; i++)
		{
			if(init)
				depth = std::min( depth, (*derivative)[i]->derivationDepth() ) ;			
			else
			{
				depth = (*derivative)[i]->derivationDepth() ;
				init = true ;
			}
		}
		depth++ ;
	}
	return depth ;
	
}

Function Function::operator*(const Function &f) const
{

	Function ret ;
	concatenateFunctions(*this, f, ret);
	ret.byteCode.push_back(TOKEN_OPERATION_TIMES)  ;
	ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0);
	ret.adress_a[(ret.byteCode.size()-1)*4+2] = 8 ;
	ret.adress_a[(ret.byteCode.size()-1)*4+1] = 9 ;
	ret.adress_a[(ret.byteCode.size()-1)*4] = 8 ;
	
	
	int n = std::max( getNumberOfDerivatives(), f.getNumberOfDerivatives() ) ;
	if(n > 0)
	{
		ret.setNumberOfDerivatives(n) ;
		for(size_t i = 0 ; i < n ; i++)
		{
			if( isDifferentiable(i) && f.isDifferentiable(i) )
			{
				Function diff = (*this) * f.d((const Variable) i) + d((const Variable) i)*f ;
				ret.setDerivative( (const Variable) i, diff ) ;
			}
		}
	}
	ret.hasGeoOp = f.hasGeoOp ;
	
	
	return ret ;
}
	
Function Function::operator/(const Function &f) const 
{

	Function ret ;
	concatenateFunctions(*this, f, ret);
	ret.byteCode.push_back(TOKEN_OPERATION_DIVIDES)  ;
	ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0);
	ret.adress_a[(ret.byteCode.size()-1)*4+2] = 8 ;
	ret.adress_a[(ret.byteCode.size()-1)*4+1] = 9 ;
	ret.adress_a[(ret.byteCode.size()-1)*4] = 8 ;
	
	int n = std::max( getNumberOfDerivatives(), f.getNumberOfDerivatives() ) ;
	if(n > 0)
	{
		ret.setNumberOfDerivatives(n) ;
		for(size_t i = 0 ; i < n ; i++)
		{
			if( isDifferentiable(i) && f.isDifferentiable(i) )
			{
				Function diff =  (d((const Variable) i)*f - (*this) * f.d((const Variable) i))/(f*f) ;
				ret.setDerivative( (const Variable) i, diff ) ;
			}
		}
	}
	ret.hasGeoOp = f.hasGeoOp ;
	
	return ret ;
}
	
Function Function::operator+(const Function &f) const
{
	int n = std::min( getNumberOfDerivatives(), f.getNumberOfDerivatives() ) ;
	std::vector<Function *> newderivatives ;
	for(size_t i = 0 ; i < n ; i++)
	{
		if( isDifferentiable(i) && f.isDifferentiable(i) )
		{
			newderivatives.push_back( new Function( f.d((const Variable) i) + d((const Variable) i))) ;
		}
	}
	
	Function ret ;
	concatenateFunctions(*this, f, ret);
	ret.byteCode.push_back(TOKEN_OPERATION_PLUS);
	ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0);
	ret.adress_a[(ret.byteCode.size()-1)*4+2] = 8 ;
	ret.adress_a[(ret.byteCode.size()-1)*4+1] = 9 ;
	ret.adress_a[(ret.byteCode.size()-1)*4] = 8 ;
	
	ret.setNumberOfDerivatives(newderivatives.size());
	
	for(size_t i = 0 ; i < newderivatives.size() ; i++)
		(*(ret.derivative))[i] = newderivatives[i] ;
	
	ret.hasGeoOp = f.hasGeoOp ;
	return ret ;
}
	
Function operator-(const double & a, const Function &f)
{
	Function ret(f) ;
	ret.values.clear();
	ret.values.push_back(a) ;
	ret.values.insert(ret.values.end(), f.values.begin(), f.values.end()) ;
	
	for(size_t i = 0 ; i < f.byteCode.size() ; i++)
	{
		if(f.hasGeoOp)
			if(f.geo_op[i])
			{
				ret.geo_op[i] = f.geo_op[i]->getCopy() ;
			}
			
		if(f.adress_a[i*4] >= HEAP_SIZE-f.values.size()-2)
			ret.adress_a[(i)*4] = f.adress_a[i*4]-1 ;
		
		if(f.adress_a[i*4+1] >= HEAP_SIZE-f.values.size()-2)
			ret.adress_a[(i)*4+1] = f.adress_a[i*4+1]-1 ;
		
		if(f.adress_a[i*4+2] >= HEAP_SIZE-f.values.size()-2)
			ret.adress_a[(i)*4+2] = f.adress_a[i*4+2]-1 ;
		
		if(f.adress_a[i*4] >= 8 && f.adress_a[i*4] < HEAP_SIZE-1-f.values.size())
			ret.adress_a[(i)*4] = f.adress_a[i*4] ;
		
		if(f.adress_a[i*4+1] >= 8 && f.adress_a[i*4+1] < HEAP_SIZE-1-f.values.size())
			ret.adress_a[(i)*4+1] = f.adress_a[i*4+1] ;
		
		if(f.adress_a[i*4+2] >= 8 && f.adress_a[i*4+2] < HEAP_SIZE-1-f.values.size())
			ret.adress_a[(i)*4+2] = f.adress_a[i*4+2] ;
		
		if(f.adress_a[i*4] < 8 )
			ret.adress_a[(i)*4] = f.adress_a[i*4] ;
		
		if(f.adress_a[i*4+1] < 8 )
			ret.adress_a[(i)*4+1] = f.adress_a[i*4+1] ;
		
		if(f.adress_a[i*4+2] < 8 )
			ret.adress_a[(i)*4+2] = f.adress_a[i*4+2] ;

	}
	ret.byteCode.push_back(TOKEN_OPERATION_MINUS);
	ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0);
	ret.adress_a[(ret.byteCode.size()-1)*4+2] = 8 ;
	ret.adress_a[(ret.byteCode.size()-1)*4+1] = 8 ;
	ret.adress_a[(ret.byteCode.size()-1)*4] = HEAP_SIZE-1 ;
	
	
	int n = f.getNumberOfDerivatives() ;
	if(n > 0)
	{
		ret.setNumberOfDerivatives(n) ;
		for(size_t i = 0 ; i < n ; i++)
		{
			if( f.isDifferentiable(i) )
			{
				Function diff =  f.d((const Variable) i) * -1;
				ret.setDerivative( (const Variable) i, diff ) ;
			}
		}
	}
	ret.hasGeoOp = f.hasGeoOp ;
	return ret ;
	
}

Function operator*(const double & a, const Function &f)
{
	Function ret(f) ;
	ret.values.clear();
	ret.values.push_back(a) ;
	ret.values.insert(ret.values.end(), f.values.begin(), f.values.end()) ;
	
	for(size_t i = 0 ; i < f.byteCode.size() ; i++)
	{
		if(f.hasGeoOp)
			if(f.geo_op[i])
			{
				ret.geo_op[i] = f.geo_op[i]->getCopy() ;
			}
			
		if(f.adress_a[i*4] >= HEAP_SIZE-2-f.values.size())
			ret.adress_a[(i)*4] = f.adress_a[i*4]-1 ;
		
		if(f.adress_a[i*4+1] >= HEAP_SIZE-2-f.values.size())
			ret.adress_a[(i)*4+1] = f.adress_a[i*4+1]-1 ;
		
		if(f.adress_a[i*4+2] >= HEAP_SIZE-2-f.values.size())
			ret.adress_a[(i)*4+2] = f.adress_a[i*4+2]-1 ;
		
		if(f.adress_a[i*4] >= 8 && f.adress_a[i*4] < HEAP_SIZE-f.values.size()-1)
			ret.adress_a[(i)*4] = f.adress_a[i*4]+1 ;
		if(f.adress_a[i*4+1] >= 8 && f.adress_a[i*4+1] < HEAP_SIZE-f.values.size()-1)
			ret.adress_a[(i)*4+1] = f.adress_a[i*4+1]+1 ;
		if(f.adress_a[i*4+2] >= 8 && f.adress_a[i*4+2] < HEAP_SIZE-f.values.size()-1)
			ret.adress_a[(i)*4+2] = f.adress_a[i*4+2]+1 ;
		if(f.adress_a[i*4] < 8 )
			ret.adress_a[(i)*4] = f.adress_a[i*4] ;
		if(f.adress_a[i*4+1] < 8 )
			ret.adress_a[(i)*4+1] = f.adress_a[i*4+1] ;
		if(f.adress_a[i*4+2] < 8 )
			ret.adress_a[(i)*4+2] = f.adress_a[i*4+2] ;

	}
	ret.byteCode.push_back(TOKEN_OPERATION_TIMES);
	ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0);
	ret.adress_a[(ret.byteCode.size()-1)*4+2] = 8 ;
	ret.adress_a[(ret.byteCode.size()-1)*4+1] = 8 ;
	ret.adress_a[(ret.byteCode.size()-1)*4] = HEAP_SIZE-1 ;
	
	
	int n = f.getNumberOfDerivatives() ;
	if(n > 0)
	{
		ret.setNumberOfDerivatives(n) ;
		for(size_t i = 0 ; i < n ; i++)
		{
			if(  f.isDifferentiable(i) )
			{
				Function diff =  f.d((const Variable) i) * a;
				ret.setDerivative( (const Variable) i, diff ) ;
			}
		}
	}
	
	return ret ;
}

Function operator+(const double & a, const Function &f)
{
	Function ret(f) ;

	ret.values.clear();
	ret.values.push_back(a) ;
	ret.values.insert(ret.values.end(), f.values.begin(), f.values.end()) ;
	
	for(size_t i = 0 ; i < f.byteCode.size() ; i++)
	{

		if(f.hasGeoOp)
			if(f.geo_op[i])
				ret.geo_op[i] = f.geo_op[i]->getCopy() ;
			
		if(f.adress_a[i*4] >= HEAP_SIZE-2-f.values.size())
			ret.adress_a[(i)*4] = f.adress_a[i*4]-1 ;
		
		if(f.adress_a[i*4+1] >= HEAP_SIZE-2-f.values.size())
			ret.adress_a[(i)*4+1] = f.adress_a[i*4+1]-1 ;
		
		if(f.adress_a[i*4+2] >= HEAP_SIZE-2-f.values.size())
			ret.adress_a[(i)*4+2] = f.adress_a[i*4+2]-1 ;
		
		if(f.adress_a[i*4] >= 8 && f.adress_a[i*4] < HEAP_SIZE-f.values.size()-1)
			ret.adress_a[(i)*4] = f.adress_a[i*4]+1 ;
		if(f.adress_a[i*4+1] >= 8 && f.adress_a[i*4+1] < HEAP_SIZE-f.values.size()-1)
			ret.adress_a[(i)*4+1] = f.adress_a[i*4+1]+1 ;
		if(f.adress_a[i*4+2] >= 8 && f.adress_a[i*4+2] < HEAP_SIZE-f.values.size()-1)
			ret.adress_a[(i)*4+2] = f.adress_a[i*4+2]+1 ;
		if(f.adress_a[i*4] < 8 )
			ret.adress_a[(i)*4] = f.adress_a[i*4] ;
		if(f.adress_a[i*4+1] < 8 )
			ret.adress_a[(i)*4+1] = f.adress_a[i*4+1] ;
		if(f.adress_a[i*4+2] < 8 )
			ret.adress_a[(i)*4+2] = f.adress_a[i*4+2] ;

	}
	ret.byteCode.push_back(TOKEN_OPERATION_PLUS);
	ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0);
	ret.adress_a[(ret.byteCode.size()-1)*4+2] = 8 ;
	ret.adress_a[(ret.byteCode.size()-1)*4+1] = 8 ;
	ret.adress_a[(ret.byteCode.size()-1)*4] = HEAP_SIZE-1 ;
	
	
	int n = f.getNumberOfDerivatives() ;
	if(n > 0)
	{
		ret.setNumberOfDerivatives(n) ;
		for(size_t i = 0 ; i < n ; i++)
		{
			if(  f.isDifferentiable(i) )
			{
				Function diff =  f.d((const Variable) i);
				ret.setDerivative( (const Variable) i, diff ) ;
			}
		}
	}
	
	return ret ;
}

Function operator/(const double & a, const Function &f)
{
	Function ret(f) ;
	ret.values.clear();
	ret.values.push_back(a) ;
	ret.values.insert(ret.values.end(), f.values.begin(), f.values.end()) ;
	
	for(size_t i = 0 ; i < f.byteCode.size() ; i++)
	{

		if(f.hasGeoOp)
			if(f.geo_op[i])
				ret.geo_op[i] = f.geo_op[i]->getCopy() ;
			
		if(f.adress_a[i*4] >= HEAP_SIZE-2-f.values.size())
			ret.adress_a[(i)*4] = f.adress_a[i*4]-1 ;
		
		if(f.adress_a[i*4+1] >= HEAP_SIZE-2-f.values.size())
			ret.adress_a[(i)*4+1] = f.adress_a[i*4+1]-1 ;
		
		if(f.adress_a[i*4+2] >= HEAP_SIZE-2-f.values.size())
			ret.adress_a[(i)*4+2] = f.adress_a[i*4+2]-1 ;
		
		if(f.adress_a[i*4] >= 8 && f.adress_a[i*4] < HEAP_SIZE-f.values.size()-1)
			ret.adress_a[(i)*4] = f.adress_a[i*4]+1 ;
		if(f.adress_a[i*4+1] >= 8 && f.adress_a[i*4+1] < HEAP_SIZE-f.values.size()-1)
			ret.adress_a[(i)*4+1] = f.adress_a[i*4+1]+1 ;
		if(f.adress_a[i*4+2] >= 8 && f.adress_a[i*4+2] < HEAP_SIZE-f.values.size()-1)
			ret.adress_a[(i)*4+2] = f.adress_a[i*4+2]+1 ;
		if(f.adress_a[i*4] < 8 )
			ret.adress_a[(i)*4] = f.adress_a[i*4] ;
		if(f.adress_a[i*4+1] < 8 )
			ret.adress_a[(i)*4+1] = f.adress_a[i*4+1] ;
		if(f.adress_a[i*4+2] < 8 )
			ret.adress_a[(i)*4+2] = f.adress_a[i*4+2] ;

	}
	ret.byteCode.push_back(TOKEN_OPERATION_DIVIDES); 
	ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0);
	ret.adress_a[(ret.byteCode.size()-1)*4+2] = 8 ;
	ret.adress_a[(ret.byteCode.size()-1)*4+1] = 8 ;
	ret.adress_a[(ret.byteCode.size()-1)*4] = HEAP_SIZE-1 ;
	
	
	int n = f.getNumberOfDerivatives() ;
	if(n > 0)
	{
		ret.setNumberOfDerivatives(n) ;
		for(size_t i = 0 ; i < n ; i++)
		{
			if(  f.isDifferentiable(i) )
			{
				Function diff =  (-a) * f.d((const Variable) i) / (f*f) ;
				ret.setDerivative( (const Variable) i, diff ) ;
			}
		}
	}
	
	return ret ;
	
}

Function Function::operator-(const Function &f) const 
{
	int n = std::min( getNumberOfDerivatives(), f.getNumberOfDerivatives() ) ;
	std::vector<Function *> newderivatives ;
	for(size_t i = 0 ; i < n ; i++)
	{
		if( isDifferentiable(i) && f.isDifferentiable(i) )
		{
			newderivatives.push_back( new Function( d((const Variable) i) - f.d((const Variable) i))) ;
		}
	}
	
	Function ret ;
	concatenateFunctions(*this, f, ret);
	ret.byteCode.push_back(TOKEN_OPERATION_MINUS);
	ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0);
	ret.adress_a[(ret.byteCode.size()-1)*4+2] = 8 ;
	ret.adress_a[(ret.byteCode.size()-1)*4+1] = 9 ;
	ret.adress_a[(ret.byteCode.size()-1)*4] = 8 ;
	
	ret.setNumberOfDerivatives(newderivatives.size());
	
	for(size_t i = 0 ; i < newderivatives.size() ; i++)
		(*(ret.derivative))[i] = newderivatives[i] ;
	
	ret.hasGeoOp = f.hasGeoOp ;
	return ret ;

}

Function Function::operator*(const double a) const
{

	Function ret(*this) ;
	
	

	if(hasGeoOp)
	{
		ret.geo_op.resize(HEAP_SIZE, (GeometryOperation*)nullptr);
		for(size_t i = 0 ; i < byteCode.size() ; i++)
		{
			if(geo_op[i])
				ret.geo_op[i] = geo_op[i]->getCopy() ;
		}
	}
	
	ret.values.push_back(a) ;
	ret.byteCode.push_back(TOKEN_OPERATION_TIMES);
	ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0);
	ret.adress_a[(ret.byteCode.size()-1)*4+2] = 8 ;
	ret.adress_a[(ret.byteCode.size()-1)*4+1] = HEAP_SIZE-values.size()-1 ;
	ret.adress_a[(ret.byteCode.size()-1)*4] = 8 ;
	
	int n =  getNumberOfDerivatives() ;

	for(size_t i = 0 ; i < n ; i++)
	{
		if( isDifferentiable(i) )
		{
			ret.d( (const Variable)i) *= a;
		}
	}
	ret.hasGeoOp = hasGeoOp ;
	return ret ;
}

Function Function::operator/(const double a) const 
{
	
	Function ret(*this) ;
	ret.values.push_back(a) ;

	if(hasGeoOp)
	{
		if(hasGeoOp)
			ret.geo_op.resize(HEAP_SIZE, (GeometryOperation*)nullptr);
		for(size_t i = 0 ; i < byteCode.size() ; i++)
		{
			if(geo_op[i])
				ret.geo_op[i] = geo_op[i]->getCopy() ;
		}
	}
	ret.byteCode.push_back(TOKEN_OPERATION_DIVIDES);
	ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0);
	ret.adress_a[(ret.byteCode.size()-1)*4+2] = 8 ;
	ret.adress_a[(ret.byteCode.size()-1)*4+1] = HEAP_SIZE-values.size()-1 ;
	ret.adress_a[(ret.byteCode.size()-1)*4] = 8 ;
	
	int n =  getNumberOfDerivatives() ;

	for(size_t i = 0 ; i < n ; i++)
	{
		if( isDifferentiable(i)  )
		{
			ret.d((const Variable) i) /= a ;
		}
	}
	ret.hasGeoOp = hasGeoOp ;
	return ret ;
}

Function Function::operator+(const double a) const
{
	Function ret(*this) ;
	ret.values.push_back(a) ;
	if(hasGeoOp)
		ret.geo_op.resize(HEAP_SIZE, (GeometryOperation*)nullptr);
	if(hasGeoOp)
	{
		for(size_t i = 0 ; i < byteCode.size() ; i++)
		{
			if(geo_op[i])
				ret.geo_op[i] = geo_op[i]->getCopy() ;
		}
	}
	ret.byteCode.push_back(TOKEN_OPERATION_PLUS);
	ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0);
	ret.adress_a[(ret.byteCode.size()-1)*4+2] = 8 ;
	ret.adress_a[(ret.byteCode.size()-1)*4+1] = HEAP_SIZE-values.size()-1 ;
	ret.adress_a[(ret.byteCode.size()-1)*4] = 8 ;
	

	return ret ;
}

Function Function::operator-(const double a) const 
{
	Function ret(*this) ;

	ret.values.push_back(a) ;

	if(hasGeoOp)
	{
		ret.geo_op.resize(HEAP_SIZE, (GeometryOperation*)nullptr);
		for(size_t i = 0 ; i < byteCode.size() ; i++)
		{
			if(geo_op[i])
				ret.geo_op[i] = geo_op[i]->getCopy() ;
		}
	}
	ret.byteCode.push_back(TOKEN_OPERATION_MINUS); 
	ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0);
	ret.adress_a[(ret.byteCode.size()-1)*4+2] = 8 ;
	ret.adress_a[(ret.byteCode.size()-1)*4+1] = HEAP_SIZE-values.size()-1 ;
	ret.adress_a[(ret.byteCode.size()-1)*4] = 8 ;
	

	return ret ;
}

Function  Function::operator^(const int a) const
{
	Function ret(*this) ;

	ret.values.push_back(a) ;

	if(hasGeoOp)
	{
		ret.geo_op.resize(HEAP_SIZE, (GeometryOperation*)nullptr);
		for(size_t i = 0 ; i < byteCode.size() ; i++)
		{
			if(geo_op[i])
				ret.geo_op[i] = geo_op[i]->getCopy() ;
		}
	}
	ret.byteCode.push_back(TOKEN_OPERATION_POWER);
	ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0);
	ret.adress_a[(ret.byteCode.size()-1)*4+2] = 8 ;
	ret.adress_a[(ret.byteCode.size()-1)*4+1] = HEAP_SIZE-values.size()-1 ;
	ret.adress_a[(ret.byteCode.size()-1)*4] = 8 ;
	
	int n =  getNumberOfDerivatives() ;
	if(n > 0)
	{
		ret.setNumberOfDerivatives(n) ;
		for(size_t i = 0 ; i < n ; i++)
		{
			if( isDifferentiable(i)  && a > 0)
			{
				Function diff =  a* d((const Variable) i) * ((*this) ^ (a-1) ) ;
				ret.setDerivative( (const Variable) i, diff ) ;
			}
		}
	}
	return ret ;
}

void Function::operator*=(const Function &f) 
{
	if(f.hasGeoOp && !hasGeoOp)
		geo_op.resize(HEAP_SIZE, (GeometryOperation*)nullptr);
	int n = std::min( getNumberOfDerivatives(), f.getNumberOfDerivatives() ) ;
	std::vector<Function *> newderivatives ;
	if(n > 0)
	{
		for(size_t i = 0 ; i < n ; i++)
		{
			if( isDifferentiable(i) && f.isDifferentiable(i) )
			{
				newderivatives.push_back( new Function( (*this) * f.d((const Variable) i) + d((const Variable) i)*f)) ;
			}
		}
	}

	concatenateFunctions(*this, f, *this);
	byteCode.push_back(TOKEN_OPERATION_TIMES);
	adress_a.push_back(0); adress_a.push_back(0); adress_a.push_back(0); adress_a.push_back(0);
	adress_a[(byteCode.size()-1)*4+2] = 8 ;
	adress_a[(byteCode.size()-1)*4+1] = 9 ;
	adress_a[(byteCode.size()-1)*4] = 8 ;
	 
	hasGeoOp = hasGeoOp  || f.hasGeoOp ;	
	
	setNumberOfDerivatives(newderivatives.size());
	for(size_t i = 0 ; i < newderivatives.size() ; i++)
		(*derivative)[i] = newderivatives[i] ;
}


void Function::operator/=(const Function &f)  
{
	if(f.hasGeoOp && !hasGeoOp)
		geo_op.resize(HEAP_SIZE, (GeometryOperation*)nullptr);
	std::vector<Function *> newderivatives ;
	int n = std::min( getNumberOfDerivatives(), f.getNumberOfDerivatives() ) ;
	if(n > 0)
	{
		for(size_t i = 0 ; i < n ; i++)
		{
			if( isDifferentiable(i) && f.isDifferentiable(i) )
			{
				newderivatives.push_back (new Function((d((const Variable)i)*f-*this*f.d((const Variable)i))/(f*f))) ;
			}
		}
	}
	
	concatenateFunctions(*this, f, *this);
	byteCode.push_back(TOKEN_OPERATION_DIVIDES)  ;
	adress_a.push_back(0); adress_a.push_back(0); adress_a.push_back(0); adress_a.push_back(0);
	adress_a[(byteCode.size()-1)*4+2] = 8 ;
	adress_a[(byteCode.size()-1)*4+1] = 9 ;
	adress_a[(byteCode.size()-1)*4] = 8 ;
	
	setNumberOfDerivatives(newderivatives.size());
	for(size_t i = 0 ; i < newderivatives.size() ; i++)
		(*derivative)[i] = newderivatives[i] ;
	hasGeoOp = hasGeoOp  || f.hasGeoOp ;
}

void Function::operator+=(const Function &f) 
{
	if(f.hasGeoOp && !hasGeoOp)
		geo_op.resize(HEAP_SIZE, (GeometryOperation*)nullptr);
	std::vector<Function *> newderivatives ;
	int n = std::min( getNumberOfDerivatives(), f.getNumberOfDerivatives() ) ;
	if(n > 0)
	{
		for(size_t i = 0 ; i < n ; i++)
		{
			if( isDifferentiable(i) && f.isDifferentiable(i) )
			{
				newderivatives.push_back (new Function(d((const Variable)i) +f.d((const Variable)i))) ;
			}
		}
	}
	
	concatenateFunctions(*this, f, *this);
	byteCode.push_back(TOKEN_OPERATION_PLUS) ;
	adress_a.push_back(0); adress_a.push_back(0); adress_a.push_back(0); adress_a.push_back(0);
	adress_a[(byteCode.size()-1)*4+2] = 8 ;
	adress_a[(byteCode.size()-1)*4+1] = 9 ;
	adress_a[(byteCode.size()-1)*4] = 8 ;
	

	setNumberOfDerivatives(newderivatives.size());
	for(size_t i = 0 ; i < newderivatives.size() ; i++)
		(*derivative)[i] = newderivatives[i] ;
	
  hasGeoOp = hasGeoOp  || f.hasGeoOp ;
}

void Function::operator-=(const Function &f)  
{
	if(f.hasGeoOp && !hasGeoOp)
		geo_op.resize(HEAP_SIZE, (GeometryOperation*)nullptr);
	std::vector<Function *> newderivatives ;
	int n = std::min( getNumberOfDerivatives(), f.getNumberOfDerivatives() ) ;
	if(n > 0)
	{
		for(size_t i = 0 ; i < n ; i++)
		{
			if( isDifferentiable(i) && f.isDifferentiable(i) )
			{
				newderivatives.push_back (new Function(d((const Variable)i) - f.d((const Variable)i))) ;
			}
		}
	}
	
	concatenateFunctions(*this, f, *this);
	byteCode.push_back(TOKEN_OPERATION_MINUS) ;
	adress_a.push_back(0); adress_a.push_back(0); adress_a.push_back(0); adress_a.push_back(0);
	adress_a[(byteCode.size()-1)*4+2] = 8 ;
	adress_a[(byteCode.size()-1)*4+1] = 9 ;
	adress_a[(byteCode.size()-1)*4] = 8 ;
	
	setNumberOfDerivatives(newderivatives.size());
	for(size_t i = 0 ; i < newderivatives.size() ; i++)
		(*derivative)[i] = newderivatives[i] ;
	
	hasGeoOp = hasGeoOp  || f.hasGeoOp ;
}

void Function::operator*=(const double a) 
{

	byteCode.push_back(TOKEN_OPERATION_TIMES) ;
	adress_a.push_back(0); adress_a.push_back(0); adress_a.push_back(0); adress_a.push_back(0);
	values.push_back(a) ;
	adress_a[(byteCode.size()-1)*4+2] = 8 ;
	adress_a[(byteCode.size()-1)*4+1] = HEAP_SIZE-values.size() ;
	adress_a[(byteCode.size()-1)*4] = 8 ;
	
	if(derivative)
	{
		for(size_t i = 0 ; i < derivative->size() ; i++)
		{
			*(*derivative) [i] *= a ;
		}
	}
}

void Function::operator/=(const double a)  
{
	if(std::abs(a-1) < POINT_TOLERANCE_2D)
		return ;
	if(std::abs(a) < POINT_TOLERANCE_2D)
	{
		std::cout << "Divide By Zero" << std::endl ;
		exit(0) ;
	}
	byteCode.push_back(TOKEN_OPERATION_DIVIDES);
	adress_a.push_back(0); adress_a.push_back(0); adress_a.push_back(0); adress_a.push_back(0);
	values.push_back(a) ;
	adress_a[(byteCode.size()-1)*4+2] = 8 ;
	adress_a[(byteCode.size()-1)*4+1] = HEAP_SIZE-values.size() ;
	adress_a[(byteCode.size()-1)*4] = 8 ;
	
	
	if(derivative)
	{
		for(size_t i = 0 ; i < derivative->size() ; i++)
		{
			if(isDifferentiable(i))
			{
				*(*derivative) [(const Variable) i] /= a ;
			}
		}
	}
	
}

void Function::operator+=(const double a) 
{
	if(std::abs(a) < POINT_TOLERANCE_2D)
		return ;
	byteCode.push_back(TOKEN_OPERATION_PLUS);
	adress_a.push_back(0); adress_a.push_back(0); adress_a.push_back(0); adress_a.push_back(0);
	values.push_back(a) ;
	adress_a[(byteCode.size()-1)*4+2] = 8 ;
	adress_a[(byteCode.size()-1)*4+1] = HEAP_SIZE-values.size() ;
	adress_a[(byteCode.size()-1)*4] = 8 ;
	
}

void Function::operator-=(const double a)  
{
	if(std::abs(a) < POINT_TOLERANCE_2D)
		return ;
	byteCode.push_back(TOKEN_OPERATION_MINUS);
	adress_a.push_back(0); adress_a.push_back(0); adress_a.push_back(0); adress_a.push_back(0);
	values.push_back(a) ;
	adress_a[(byteCode.size()-1)*4+2] = 8 ;
	adress_a[(byteCode.size()-1)*4+1] = HEAP_SIZE-values.size() ;
	adress_a[(byteCode.size()-1)*4] = 8 ;
	
}

void Function::setNumberOfDerivatives(int n) 
{
	if(derivative)
	{
		for(size_t i = 0 ; i < derivative->size() ; i++)
			delete (*derivative)[i] ;
		derivative->resize(n,(Function *)nullptr) ;
	}
	else 
		derivative = new std::valarray<Function *>((Function *)nullptr, n) ;

	e_diff = n > 0 ;
	
	if (n == 0)
	{
		if(derivative)
			for(size_t i = 0 ; i < derivative->size() ; i++)
				delete (*derivative)[i] ;
		delete derivative ;
		derivative = nullptr ;
	}
}

void Function::setDerivative( const Variable v, Function & f) 
{
	if(derivative)
	{
		if(derivative->size() < v+1)
		{
			std::vector<Function *> oldfuncts ;
			for(size_t i = 0 ; i < getNumberOfDerivatives() ; i++)
				oldfuncts.push_back((*derivative)[i]);
			
			derivative->resize(v+1,(Function *)nullptr);
			for(size_t i = 0 ; i < oldfuncts.size() ; i++)
				(*derivative)[i] = oldfuncts[i] ;
			for(size_t i = oldfuncts.size() ; i < derivative->size() ; i++)
				(*derivative)[i] = new Function("0") ;
		}
		delete (*derivative)[v] ;
		(*derivative)[v] = new Function(f) ;
	}
}

const Function & Function::d(const Variable v) const
{
	return *(*derivative)[v] ;
}

Function & Function::d(const Variable v) 
{
	return *(*derivative)[v] ;
}

const Function & Function::transform(size_t i) const
{
	return *(*transforms)[i] ;
}

Function & Function::transform(size_t i) 
{
	return *(*transforms)[i] ;
}

std::valarray<Function *> & Function::getDerivatives() const
{
	return *derivative ;
}

std::valarray<Function *> & Function::getDerivatives()
{
	return *derivative ;
}


Function f_exp(const Function &f)
{
	Function ret = f ;
	ret.byteCode.push_back(TOKEN_OPERATION_EXP);
	ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0);
	ret.adress_a[(ret.byteCode.size()-1)*4+2] = 8 ;
	ret.adress_a[(ret.byteCode.size()-1)*4+1] = 9 ;
	ret.adress_a[(ret.byteCode.size()-1)*4] = 8 ;
	
	return ret ;
}

Function f_abs(const Function &f, bool differentiate)
{
	Function ret = f ;
	ret.byteCode.push_back(TOKEN_OPERATION_ABS); 
	ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0);
	ret.adress_a[(ret.byteCode.size()-1)*4+2] = 8 ;
	ret.adress_a[(ret.byteCode.size()-1)*4+1] = 9 ;
	ret.adress_a[(ret.byteCode.size()-1)*4] = 8 ;
	
	
	ret.setNumberOfDerivatives(f.getNumberOfDerivatives()*differentiate);
	
	if(differentiate && f.isDifferentiable())
	{
		for(size_t i = 0 ; i < f.getNumberOfDerivatives() ; i++)
		{
			Function d = f.d((const Variable)i)*f_positivity(f,false)-f.d((const Variable)i)*f_negativity(f,false) ;
			if(f.d((const Variable) i).isDifferentiable())
			{
				d.setNumberOfDerivatives( f.d((const Variable) i).getNumberOfDerivatives() );
				for(size_t k = 0 ; k < d.getNumberOfDerivatives() ; k++)
				{
					Function g = f.d((const Variable)i).d((const Variable) k)*f_positivity(f,false)-f.d((const Variable)i).d((const Variable) k)*f_negativity(f,false) ;
					d.setDerivative((const Variable) k, g);
				}
			}
			ret.getDerivatives()[i] =  new Function(d) ;
		}
	}
	return ret ;
}

Function f_log(const Function &f)
{
	Function ret(f) ;
	ret.byteCode.push_back(TOKEN_OPERATION_LOG);
	ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0);
	ret.adress_a[(ret.byteCode.size()-1)*4+2] = 8 ;
	ret.adress_a[(ret.byteCode.size()-1)*4+1] = 9 ;
	ret.adress_a[(ret.byteCode.size()-1)*4] = 8 ;
	
	return ret ;
}

Function f_sqrt(const Function &f, bool differentiate)
{
	Function ret(f) ;
	ret.byteCode.push_back(TOKEN_OPERATION_SQRT);
	ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0);
	ret.adress_a[(ret.byteCode.size()-1)*4+2] = 8 ;
	ret.adress_a[(ret.byteCode.size()-1)*4+1] = 9 ;
	ret.adress_a[(ret.byteCode.size()-1)*4] = 8 ;
	

 	ret.setNumberOfDerivatives(0);//f.getNumberOfDerivatives()*differentiate) ;
// 	
// 	if(differentiate && f.isDifferentiable())
// 	{
// 		for(size_t i = 0 ; i < f.getNumberOfDerivatives() ; i++)
// 		{
// 			Function d = f.d((const Variable)i)*0.5/f_sqrt(f,f.isDifferentiable() && f.d((const Variable)i).isDifferentiable()) ;
// 			ret.getDerivatives()[i] = new Function(d) ;
// 		}
// 	}
	return ret ;
}

Function f_atan2(const Function &f0, const Function &f1)
{
	Function ret ;
	concatenateFunctions(f0, f1, ret);
	ret.byteCode.push_back(TOKEN_OPERATION_ATAN2);
	ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0);
	ret.adress_a[(ret.byteCode.size()-1)*4+2] = 8 ;
	ret.adress_a[(ret.byteCode.size()-1)*4+1] = 9 ;
	ret.adress_a[(ret.byteCode.size()-1)*4] = 8 ;
	return ret ;
}

Function f_sin(const Function &f)
{
	Function ret = f ;
	ret.byteCode.push_back(TOKEN_OPERATION_SIN);
	ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0);
	ret.adress_a[(ret.byteCode.size()-1)*4+2] = 8 ;
	ret.adress_a[(ret.byteCode.size()-1)*4+1] = 9 ;
	ret.adress_a[(ret.byteCode.size()-1)*4] = 8 ;
	
	return ret ;
}


Mu::Function f_sign(const Mu::Function &f)
{
	Function ret = f ;
	ret.byteCode.push_back(TOKEN_OPERATION_SIGN);
	ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0);
	ret.adress_a[(ret.byteCode.size()-1)*4+2] = 8 ;
	ret.adress_a[(ret.byteCode.size()-1)*4+1] = 9 ;
	ret.adress_a[(ret.byteCode.size()-1)*4] = 8 ;
	
	return ret ;
}



Mu::Function f_positivity(const Mu::Function &f, bool differentiate)
{
	Function ret(f) ;
	ret.byteCode.push_back(TOKEN_OPERATION_POSITIVITY);
	ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0);
	ret.adress_a[(ret.byteCode.size()-1)*4+2] = 8 ;
	ret.adress_a[(ret.byteCode.size()-1)*4+1] = 9 ;
	ret.adress_a[(ret.byteCode.size()-1)*4] = 8 ;
	
	
	if(ret.isDifferentiable() && !differentiate)
	{
		ret.setNumberOfDerivatives(0);
	}
	else if(ret.isDifferentiable())
	{
		Function zero("0") ;
		
		size_t  j = 0 ;
		while(j < ret.getDerivatives().size())
		{
			if( f.isDifferentiable((const Variable) j) )
			{
				Function z("0") ;
				zero.setNumberOfDerivatives( ret.getNumberOfDerivatives() );
				for(size_t k = 0 ; k < ret.getDerivatives().size() ; k++)
				{
					zero.setDerivative((const Variable) k, z);
				}
				j = ret.getDerivatives().size() ;
			}
		}
		
		
		for(size_t i = 0 ; i < ret.getDerivatives().size() ; i++)
			ret.setDerivative((const Variable)i,  zero) ;
	}
	return ret ;
}


Mu::Function f_negativity(const Mu::Function &f, bool differentiate)
{
	Function ret = f ;
	ret.byteCode.push_back(TOKEN_OPERATION_NEGATIVITY);
	ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0);
	ret.adress_a[(ret.byteCode.size()-1)*4+2] = 8 ;
	ret.adress_a[(ret.byteCode.size()-1)*4+1] = 9 ;
	ret.adress_a[(ret.byteCode.size()-1)*4] = 8 ;
	
	if(ret.isDifferentiable() && !differentiate)
	{
		ret.setNumberOfDerivatives(0);
	}
	else
	{
		Function zero("0") ;

		size_t  j = 0 ;
		while(j < ret.getDerivatives().size())
		{
			if( f.isDifferentiable((const Variable) j) )
			{
				Function z("0") ;
				zero.setNumberOfDerivatives( ret.getNumberOfDerivatives() );
				for(size_t k = 0 ; k < ret.getDerivatives().size() ; k++)
				{
					zero.setDerivative((const Variable) k, z);
				}
				j = ret.getDerivatives().size() ;
			}
		}
		
		for(size_t i = 0 ; i < ret.getDerivatives().size() ; i++)
			ret.setDerivative((const Variable)i,  zero) ;
	}
	return ret ;
}


Function f_cos(const Function &f)
{
	Function ret = f ;
	ret.byteCode.push_back(TOKEN_OPERATION_COS) ;
	ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0); ret.adress_a.push_back(0);
	ret.adress_a[(ret.byteCode.size()-1)*4+2] = 8 ;
	ret.adress_a[(ret.byteCode.size()-1)*4+1] = 9 ;
	ret.adress_a[(ret.byteCode.size()-1)*4] = 8 ;
	
	return ret ;
} 

const std::vector< Point > & Function::getIntegrationHint() const
{
	return iPoint ;
}

const Point & Function::getIntegrationHint(size_t i) const
{
	return iPoint[i] ;
}

void Function::setIntegrationHint(const std::vector< Point > v)
{
	iPoint.clear() ;
	iPoint = v ;
}

void Function::addIntegrationHint(const Point  p)
{
	iPoint.push_back(p) ;
}

bool Function::hasIntegrationHint() const
{
	return iPoint.size() > 0 ;
}





bool Function::isBinaryOperator(TokenOperationType t) const
{
	return ( t== TOKEN_OPERATION_PLUS || 
	         t== TOKEN_OPERATION_MINUS || 
	         t== TOKEN_OPERATION_TIMES || 
	         t== TOKEN_OPERATION_DIVIDES || 
	         t== TOKEN_OPERATION_POWER ||
	         t== TOKEN_OPERATION_INTERPOLATE ||
	         t== TOKEN_OPERATION_ATAN2) ;
}

bool Function::isBinaryVectorisableOperator(TokenOperationType t) const
{
	return ( t == TOKEN_OPERATION_PLUS || 
	         t == TOKEN_OPERATION_MINUS || 
	         t == TOKEN_OPERATION_DIVIDES || 
	         t == TOKEN_OPERATION_POWER ) ;
}

bool Function::isUnaryOperator(TokenOperationType t) const
{
	return ( t == TOKEN_OPERATION_SIN || 
	         t == TOKEN_OPERATION_COS || 
	         t == TOKEN_OPERATION_SQRT || 
	         t == TOKEN_OPERATION_COSH || 
	         t == TOKEN_OPERATION_SINH ||
	         t == TOKEN_OPERATION_TANH || 
	         t == TOKEN_OPERATION_EXP || 
	         t == TOKEN_OPERATION_LOG || 
	         t == TOKEN_OPERATION_TAN  
	       ) ;
}
