//
// C++ Implementation: mechanical analytic homogenization
//
// Description:
//
//
// Author:  Alain Giorla, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "elastic_homogenization.h"
#include "../../utilities/matrixops.h"

namespace Mu
{


Matrix cauchyGreen(std::pair<double,double> prop, bool hooke, SpaceDimensionality dim) 
{
	double E = prop.first ;
	double nu = prop.second ;
	if(!hooke)
	{
		double k = prop.first ;
		double mu = prop.second ;
		E = 9*k*mu / (3*k+mu) ;
		nu = (3*k-2*mu) / (6*k+2*mu) ;
	}
	switch(dim)
	{
		case SPACE_TWO_DIMENSIONAL:
		{
			Matrix cg(3,3) ;
			cg[0][0] = E/(1.-nu*nu) ; cg[0][1] =E/(1.-nu*nu)*nu ; cg[0][2] = 0 ;
			cg[1][0] = E/(1.-nu*nu)*nu ; cg[1][1] = E/(1.-nu*nu) ; cg[1][2] = 0 ; 
			cg[2][0] = 0 ; cg[2][1] = 0 ; cg[2][2] = .99*E/(1.-nu*nu)*(1.-nu)/2. ; 
			return cg ;
		}
		case SPACE_THREE_DIMENSIONAL:
		{
			Matrix cgg(6,6) ;
			cgg[0][0] = 1. - nu ; cgg[0][1] = nu ; cgg[0][2] = nu ;
			cgg[1][0] = nu ; cgg[1][1] = 1. - nu ; cgg[1][2] = nu ;
			cgg[2][0] = nu ; cgg[2][1] = nu ; cgg[2][2] = 1. - nu ;
			cgg[3][3] = 0.5 - nu ;
			cgg[4][4] = 0.5 - nu ;
			cgg[5][5] = 0.5 - nu ;
			cgg *= E/((1.+nu)*(1.-2.*nu)) ;
			return cgg ;
		}
	}
	return Matrix(0,0) ;
}






ElasticHomogenizationScheme::ElasticHomogenizationScheme(int i) : Scheme(i)
{
	input.push_back(TAG_VOLUME_FRACTION) ;
	input.push_back(TAG_BULK_MODULUS) ;
	input.push_back(TAG_SHEAR_MODULUS) ;

	output.push_back(TAG_BULK_MODULUS) ;
	output.push_back(TAG_SHEAR_MODULUS) ;
}





DilutedScheme::DilutedScheme() : ElasticHomogenizationScheme(2)
{
}

Vector DilutedScheme::process(const Matrix & data)
{
	double k_mat = data[0][1] ;
	double mu_mat = data[0][2] ;

	double f_inc = data[1][0] ;
	double k_inc = data[1][1] ;
	double mu_inc = data[1][2] ;

	double k_hom = k_mat + f_inc * (k_inc - k_mat) * (3*k_mat + 4*mu_mat) / (3*k_inc + 4*mu_mat) ;
	double mu_hom = mu_mat + f_inc * (mu_inc - mu_mat) * (5*mu_mat*(3*k_mat + 4*mu_mat)) / (mu_mat*(9*k_mat + 8*mu_mat) + 6*mu_inc*(k_mat + 2*mu_mat)) ;

	Vector processed(2) ;
	processed[0]=k_hom ;
	processed[1]=mu_hom ;

	return processed ;
}





GeneralizedDilutedScheme::GeneralizedDilutedScheme() : ElasticHomogenizationScheme(-1)
{
}

Vector GeneralizedDilutedScheme::process(const Matrix & data)
{
	double k_mat = data[0][1] ;
	double mu_mat = data[0][2] ;

	double k_hom = k_mat ;
	double mu_hom = mu_mat ;

	for(size_t i = 1 ; i < data.size() ; i++)
	{
		double f_inc = data[i][0] ;
		double k_inc = data[i][1] ;
		double mu_inc = data [i][2] ;

		k_hom += f_inc * (k_inc - k_mat) * (3*k_mat + 4*mu_mat) / (3*k_inc + 4*mu_mat) ;
		mu_hom += f_inc * (mu_inc - mu_mat) * (5*mu_mat*(3*k_mat + 4*mu_mat)) / (mu_mat*(9*k_mat + 8*mu_mat) + 6*mu_inc*(k_mat + 2*mu_mat)) ;
	}

	Vector processed(2) ;
	processed[0]=k_hom ;
	processed[1]=mu_hom ;

	return processed ;
}


IncrementalScheme::IncrementalScheme(double d) : ElasticHomogenizationScheme(2)
{
	dalpha = d;
}

Vector IncrementalScheme::process(const Matrix & data)
{
	double k_mat = data[0][1] ;
	double mu_mat = data[0][2] ;

	double k_hom = k_mat ;
	double mu_hom = mu_mat ;

	double k_inc = data[1][1] ;
	double mu_inc = data[1][2] ;

	double f = 0 ;
	while(f < data[1][0])
	{
		k_hom += (k_inc-k_hom) * ((3*k_hom + 4*mu_hom) / (3*k_inc + 4*mu_hom)) * dalpha/(1-f) ;
		mu_hom += (mu_inc - mu_hom) * (5*mu_hom*(3*k_hom + 4*mu_hom) / (mu_hom*(9*k_hom+8*mu_hom) + 6*mu_inc*(k_hom+2*mu_hom))) * dalpha/(1-f) ;
		f += dalpha ;
	}

	double dlast = data[1][0] - (f-dalpha) ;
	k_hom += (k_inc-k_hom) * ((3*k_hom + 4*mu_hom) / (3*k_inc + 4*mu_hom)) * dlast/(1-f) ;
	mu_hom += (mu_inc - mu_hom) * (5*mu_hom*(3*k_hom + 4*mu_hom) / (mu_hom*(9*k_hom+8*mu_hom) + 6*mu_inc*(k_hom+2*mu_hom))) * dlast/(1-f) ;

	Vector processed(2) ;
	processed[0]=k_hom ;
	processed[1]=mu_hom ;

	return processed ;

}






MoriTanaka::MoriTanaka() : ElasticHomogenizationScheme(2)
{
}

Vector MoriTanaka::process(const Matrix & data)
{
	double k_mat = data[0][1] ;
	double mu_mat = data[0][2] ;

	double k_hom = k_mat ;
	double mu_hom = mu_mat  ;

	double f_inc = data[1][0] ;
	double k_inc = data[1][1] ;
	double mu_inc = data[1][2] ;

	double Sk = (4*mu_mat) / (3*k_mat) ;
	double Smu = (3 + 2*Sk) / (2 + 3*Sk) ;
	double K = k_inc / k_mat ;
	double Mu = mu_inc / mu_mat ;

	k_hom *= (K/Sk + 1 + f_inc*(K-1)) ;
	k_hom /= (1 + (K + f_inc*(1-K))/Sk) ;
	mu_hom *= (Mu/Smu + 1 + f_inc*(Mu-1)) ;
	mu_hom /= (1 + (Mu + f_inc*(1-Mu))/Smu) ;

	Vector processed(2) ;
	processed[0]=k_hom ;
	processed[1]=mu_hom ;

	return processed ;
}



GeneralizedMoriTanaka::GeneralizedMoriTanaka() : ElasticHomogenizationScheme(-1)
{
}

Vector GeneralizedMoriTanaka::process(const Matrix & data)
{
	double k_mat = data[0][1] ;
	double mu_mat = data[0][2] ;

	double k_hom = k_mat * data[0][0] ;
	double mu_hom = mu_mat * data[0][0] ;

	for(size_t i = 1 ; i < data.size() ; i++)
	{
		Matrix idata(2,3) ;
		idata[0][0] = data[0][0] ;
		idata[0][1] = data[0][1] ;
		idata[0][2] = data[0][2] ;
		idata[1][0] = data[i][0] ;
		idata[1][1] = data[i][1] ;
		idata[1][2] = data[i][2] ;
		idata[0][0] = 1 - idata[1][0] ;
		Vector mt_res(2) ;
		mt_res = MoriTanaka().process(idata) ;

		k_hom += data[i][0] * (mt_res[0]) ;
		mu_hom += data[i][0] * (mt_res[1]) ;
	}

	Vector processed(2) ;
	processed[0]=k_hom ;
	processed[1]=mu_hom ;

	return processed ;
}



SelfConsistent::SelfConsistent() : ElasticHomogenizationScheme(2)
{
}

Vector SelfConsistent::process(const Matrix & data)
{
	double & f_mat = data[0][0] ;
	double & k_mat = data[0][1] ;
	double & mu_mat = data[0][2] ;
	double & f_inc = data[1][0] ;
	double & k_inc = data[1][1] ;
	double & mu_inc = data[1][2] ;

	Vector data_hom(3) ;
	data_hom[1] = (k_mat*f_mat+k_inc*f_inc) ;
	data_hom[2] = (mu_mat*mu_mat+mu_inc*mu_inc) ;
	double & k_hom = data_hom[1] ;
	double & mu_hom = data_hom[2] ;


	double error = 1. ;
	Vector gmt_mat(2) ;
	Vector gmt_inc(2) ;

	Matrix data_mat(2,3) ;
	Matrix data_inc(2,3) ;

	double & f_mat_hom = data_mat[0][0] ;
	double & k_mat_hom = data_mat[0][1] ;
	double & mu_mat_hom = data_mat[0][2] ;
	double & f_inc_hom = data_inc[0][0] ;
	double & k_inc_hom = data_inc[0][1] ;
	double & mu_inc_hom = data_inc[0][2] ;

//	for(size_t i = 0 ; i < data_hom.size() ; i++)
//	{
		data_mat[0][1] = k_hom ;
		data_inc[0][1] = k_hom ;
		data_mat[0][2] = mu_hom ;
		data_inc[0][2] = mu_hom ;

		data_mat[1][0] = f_mat ;
		data_inc[1][0] = f_inc ;
		data_mat[1][1] = k_mat ;
		data_inc[1][1] = k_inc ;
		data_mat[1][2] = mu_mat ;
		data_inc[1][2] = mu_inc ;
//	}


	f_inc_hom = 1.- f_inc ;
	f_mat_hom = 1.- f_mat ;

	double & k_mat_mt = gmt_mat[0] ;
	double & mu_mat_mt = gmt_mat[1] ;
	double & k_inc_mt = gmt_inc[0] ;
	double & mu_inc_mt = gmt_inc[1] ;

	double Ak_mat_mt = 0. ;
	double Ak_inc_mt = 0. ;
	double Amu_mat_mt = 0. ;
	double Amu_inc_mt = 0. ;

	while(error > 1e-6)
	{
		gmt_mat = MoriTanaka().process(data_mat) ;
		gmt_inc = MoriTanaka().process(data_inc) ;

		Ak_mat_mt = (k_mat_mt-k_mat_hom)/(k_mat-k_mat_hom) /f_mat ;
		Ak_inc_mt = (k_inc_mt-k_inc_hom)/(k_inc-k_inc_hom) /f_inc ;

		Amu_mat_mt = (mu_mat_mt-mu_mat_hom)/(mu_mat-mu_mat_hom) /f_mat ;
		Amu_inc_mt = (mu_inc_mt-mu_inc_hom)/(mu_inc-mu_inc_hom) /f_inc ;

		double ksc = k_mat + f_inc*(k_inc-k_mat)*Ak_inc_mt/(f_inc*Ak_inc_mt+f_mat*Ak_mat_mt) ;
		double musc = mu_mat + f_inc*(mu_inc-mu_mat)*Amu_inc_mt/(f_inc*Amu_inc_mt+f_mat*Amu_mat_mt) ;

		error = std::max(std::abs(k_mat_hom - ksc)/k_mat_hom,std::abs(mu_mat_hom - musc)/mu_mat_hom) ;

		k_mat_hom = ksc ;
		k_inc_hom = ksc ;

		mu_mat_hom = musc ;
		mu_inc_hom = musc ;

	}


	Vector processed(2) ;
	processed[0]=k_mat_hom ;
	processed[1]=mu_mat_hom ;

	return processed ;
}


GeneralizedSelfConsistent::GeneralizedSelfConsistent() : ElasticHomogenizationScheme(-1)
{
}

Vector GeneralizedSelfConsistent::process(const Matrix & data)
{
	Vector kmu_mean(2) ;
	std::vector<Tag> avg ;
	avg.push_back(TAG_BULK_MODULUS) ;
	avg.push_back(TAG_SHEAR_MODULUS) ;
	MeanScheme mean(true,true,avg) ;
	kmu_mean = mean.process(data) ;
	if(!mean.isOK())
	{
		mean.print() ;
		s = STATUS_BAD_HOMOGENIZATION ;
	}
	Vector hom(3) ;
	hom[1] = kmu_mean[0] ;
	hom[2] = kmu_mean[1] ;
	Vector next_hom(3) ;
	next_hom = hom ;
	double error = 1. ;
	int nit = 0 ;
	Matrix iscm(2,3) ;		
	Matrix isc(2,3) ;
	Vector res_mt(2) ;
	Vector resm_mt(2) ;
	while(error > 1e-6)
	{
		error = 1. ;
		Vector dk(data.numRows()) ;
		Vector dmu(data.numRows()) ;

		iscm.array() = 0 ;
		iscm[0][0] = 1. - data[0][0] ;
		iscm[0][1] = hom[1] ;
		iscm[0][2] = hom[2] ;
		iscm[1][0] = data[0][0] ;
		iscm[1][1] = data[0][1] ;
		iscm[1][2] = data[0][2] ;
		resm_mt = 0 ;
		
		resm_mt = MoriTanaka().process(iscm) ;

		double fm = iscm[1][0] ;

		double khomm = iscm[0][1] ;
		double km = iscm[1][1] ;
		double kmtm = resm_mt[0] ;

		double muhomm = iscm[0][2] ;
		double mum = iscm[1][2] ;
		double mumtm = resm_mt[1] ;


		double Akm = (kmtm-khomm)/(km-khomm)/fm ;
		double Amum = (mumtm-muhomm)/(mum-muhomm)/fm ;

		dk[0] = data[0][1] ;
		dmu[0] = data[0][2] ;

		for(size_t i = 1 ; i < data.numRows() ; i++)
		{
			isc.array() = 0 ;
			res_mt = 0 ;
			isc[0][0] = 1. - data[i][0] ;
			isc[0][1] = hom[1] ;
			isc[0][2] = hom[2] ;
			isc[1][0] = data[i][0] ;
			isc[1][1] = data[i][1] ;
			isc[1][2] = data[i][2] ;

			
			res_mt = MoriTanaka().process(isc) ;

			double fi = isc[1][0] ;

			double khom = isc[0][1] ;
			double ki = isc[1][1] ;
			double kmt = res_mt[0] ;

			double muhom = isc[0][2] ;
			double mui = isc[1][2] ;
			double mumt = res_mt[1] ;


			double Ak = (kmt-khom)/(ki-khom)/fi ;
			double Amu = (mumt-muhom)/(mui-muhom)/fi ;
			
			dk[i] = fi*(ki-km)*Ak/(fi*Ak+(1.-fi)*Akm) ;
			dmu[i] = fi*(mui-mum)*Amu/(fi*Amu+(1.-fi)*Amum) ;

		}

		next_hom[1] = 0. ;
		next_hom[2] = 0. ;

		for(size_t i = 0 ; i < dk.size() ; i++)
		{
			next_hom[1] += dk[i] ;
			next_hom[2] += dmu[i] ;
		}

		error = std::max(std::abs(hom[1]-next_hom[1])/hom[1],std::abs(hom[2]-next_hom[2])/hom[2]) ;

		hom = next_hom ;

	}


	Vector processed(2) ;
	processed[0]=hom[1] ;
	processed[1]=hom[2] ;

	return processed ;
}



}


