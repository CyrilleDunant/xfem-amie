#include "scheme_template.h"
#include "homogenization_base.h"

using namespace Mu ;

PhaseTemplate::PhaseTemplate(int n)
{
	max = n ;
	done = 0 ;
	prop.clear() ;
}


PhaseTemplate::PhaseTemplate(int n, std::vector<PType> t)
{
	max = n ;
	done = 0 ;
	prop.clear() ;
	for(int j = 0 ; j < (int) t.size() ; j++)
	{
		if(!hasProperties(t[j]))
		{
			prop.push_back(Properties(t[j],0.)) ;
		}
	}
}

bool PhaseTemplate::hasProperties(PType t) const
{
	return indexOfProperties(t) != -1 ;
}

int PhaseTemplate::indexOfProperties(PType t) const
{
	for(int i = 0 ; i < (int) prop.size() ; i++)
	{
		if(prop[i].is(t))
		{
			return i ;
		}
	}
	return -1 ;
}

	
bool PhaseTemplate::cast(Material m, Material father)
{
	if(filled() || (m.sizeProperties() < (int) prop.size()))
	{
		return false ;
	}
	else
	{
		bool valid = true ;
		int i = 0 ;
		while(valid && i < (int) prop.size())
		{
			valid = cast(i, m.searchProperties(prop[i].type(), SEARCH_TRY_ANYTHING, true, father)) ;
			i++ ;
		}
		done++ ;
		return valid ;
	}
}

bool PhaseTemplate::cast(int i, Properties p)
{
	if(p.is(prop[i].type()))
	{
		prop[i].set(p.val()) ;
		return true ;
	}
	else
	{
		return false ;
	}
}

std::vector<double> PhaseTemplate::val() const
{
	std::vector<double> v ;
	for(int i = 0 ; i < (int) prop.size() ; i++)
	{
		v.push_back(prop[i].val()) ;
	}
	return v ;
}

int PhaseTemplate::size() const 
{
	return prop.size() ;
}

int PhaseTemplate::phases(int count) const 
{
	if(max == -1)
	{
		return count ;
	}
	else
	{
		return max ;
	}
}
	
bool PhaseTemplate::filled() const
{
	if(max == -1)
	{
		return false ;
	}
	else
	{
		return done >= max ;
	}
}

void PhaseTemplate::force(std::vector<double> data)
{
	for(int i = 0 ; i < std::min((int) data.size(), (int) prop.size()) ; i++)
	{
		prop[i].set(data[i]) ;
	}
}

Properties PhaseTemplate::get(int i) const
{
	if(i < 0 || i >= (int) prop.size())
	{
		return Properties(P_BAD_INDEX,0.) ;
	}
	else
	{
		return prop[i] ;
	}
}

PhaseTemplate PhaseTemplate::makeBulkShearTemplate(int n) 
{
	std::vector<PType> t ;
	t.push_back(P_BULK_MODULUS) ;
	t.push_back(P_SHEAR_MODULUS) ;
	return PhaseTemplate(n, t) ;
}

PhaseTemplate PhaseTemplate::makeVolumeBulkShearTemplate(int n) 
{
	std::vector<PType> t ;
	t.push_back(P_VOLUME_FRACTION) ;
	t.push_back(P_BULK_MODULUS) ;
	t.push_back(P_SHEAR_MODULUS) ;
	return PhaseTemplate(n, t) ;
}

PhaseTemplate PhaseTemplate::makeExpansionTemplate(int n)
{
	std::vector<PType> t ;
	t.push_back(P_EXPANSION_COEFFICIENT) ;
	return PhaseTemplate(n, t) ;
}

PhaseTemplate PhaseTemplate::makeBulkExpansionTemplate(int n)
{
	std::vector<PType> t ;
	t.push_back(P_BULK_MODULUS) ;
	t.push_back(P_EXPANSION_COEFFICIENT) ;
	return PhaseTemplate(n, t) ;
}

PhaseTemplate PhaseTemplate::makeVolumeBulkExpansionTemplate(int n)
{
	std::vector<PType> t ;
	t.push_back(P_VOLUME_FRACTION) ;
	t.push_back(P_BULK_MODULUS) ;
	t.push_back(P_EXPANSION_COEFFICIENT) ;
	return PhaseTemplate(n, t) ;
}

PhaseTemplate PhaseTemplate::makeVolumeBulkShearExpansionTemplate(int n)
{
	std::vector<PType> t ;
	t.push_back(P_VOLUME_FRACTION) ;
	t.push_back(P_BULK_MODULUS) ;
	t.push_back(P_SHEAR_MODULUS) ;
	t.push_back(P_EXPANSION_COEFFICIENT) ;
	return PhaseTemplate(n, t) ;
}





SchemeTemplate::SchemeTemplate(PhaseTemplate r, std::vector<PhaseTemplate> p)
{
	phases.clear() ;
	data.clear() ;
	for(int i = 0 ; i < (int) p .size() ; i++)
	{
		phases.push_back(p[i]) ;
	}
	scheme = SCHEME_DO_NOTHING ;
}

SchemeTemplate::SchemeTemplate(HomogenizationScheme s)
{
	phases.clear() ;
	data.clear() ;
	scheme = s ;
	if(s == ELASTICITY_DILUTED
		|| s == ELASTICITY_INCREMENTAL
		|| s == ELASTICITY_MORI_TANAKA
		|| s == ELASTICITY_SELF_CONSISTENT)
	{
		PhaseTemplate matrix = PhaseTemplate::makeBulkShearTemplate(1) ; 
		PhaseTemplate inclusion = PhaseTemplate::makeVolumeBulkShearTemplate(1) ; 
		phases.push_back(matrix) ;
		phases.push_back(inclusion) ;
		result = PhaseTemplate::makeBulkShearTemplate(1) ;
	}
		
	if(s == ELASTICITY_GENERALIZED_DILUTED
		|| s == ELASTICITY_GENERALIZED_MORI_TANAKA
		|| s == ELASTICITY_GENERALIZED_SELF_CONSISTENT)
	{
		PhaseTemplate matrix = PhaseTemplate::makeBulkShearTemplate(1) ; 
		PhaseTemplate inclusion = PhaseTemplate::makeVolumeBulkShearTemplate(-1) ; 
		phases.push_back(matrix) ;
		phases.push_back(inclusion) ;
		result = PhaseTemplate::makeBulkShearTemplate(1) ;
	}

	if(s == EXPANSION_HOBBS || s == EXPANSION_TURNER)
	{
		PhaseTemplate matrix = PhaseTemplate::makeBulkExpansionTemplate(1) ; 
		PhaseTemplate inclusion = PhaseTemplate::makeVolumeBulkExpansionTemplate(1) ; 
		phases.push_back(matrix) ;
		phases.push_back(inclusion) ;
		result = PhaseTemplate::makeExpansionTemplate(1) ;
	}

	if(s == EXPANSION_KERNER)
	{
		PhaseTemplate matrix = PhaseTemplate::makeVolumeBulkShearExpansionTemplate(1) ; 
		PhaseTemplate inclusion = PhaseTemplate::makeBulkExpansionTemplate(1) ; 
		phases.push_back(matrix) ;
		phases.push_back(inclusion) ;
		result = PhaseTemplate::makeExpansionTemplate(1) ;
	}	
	

}

int SchemeTemplate::sizeResult() const
{
	return result.size() ;
}

Properties SchemeTemplate::getResult(int i) const 
{
	return result.get(i) ;
}
	
bool SchemeTemplate::cast(Material m)
{
	if(scheme == SCHEME_DO_NOTHING)
	{
		return true ;
	}

	int count = m.sizePhase() ;
	int cast = 0 ;
	for(int i = 0 ; i < (int) phases.size() ; i++)
	{
		int here = phases[i].phases(count) ;
		cast += here ;
		count -= here ;
		if(count < 0)
		{
			return false ;
		}
	}
	int j = 0 ;
	for(int i = 0 ; i < (int) phases.size() ; i++)
	{
		while(!phases[i].filled() && j < m.sizePhase())
		{
			if(phases[i].cast(m.getPhase(j), m))
			{
				data.push_back(phases[i].val()) ;
				j++ ;
			}
			else
			{
				return false ;
			}
		}
	}
	result.force(makeScheme(scheme, data)) ;
	return true ;
}

std::vector<double> SchemeTemplate::makeScheme(HomogenizationScheme s, std::vector<std::vector<double> > d) 
{
	std::vector<double> empty ;
	switch(s)
	{
		case ELASTICITY_DILUTED:
			return SchemeTemplate::elasticityDilutedScheme(d) ;
		case ELASTICITY_GENERALIZED_DILUTED:
			return SchemeTemplate::elasticityGeneralizedDilutedScheme(d) ;
		case ELASTICITY_INCREMENTAL:
			return SchemeTemplate::elasticityIncrementalScheme(d) ;
		case ELASTICITY_MORI_TANAKA:
			return SchemeTemplate::elasticityMoriTanakaScheme(d) ;
		case ELASTICITY_GENERALIZED_MORI_TANAKA:
			return SchemeTemplate::elasticityGeneralizedMoriTanakaScheme(d) ;
		case ELASTICITY_SELF_CONSISTENT:
			return SchemeTemplate::elasticitySelfConsistentScheme(d) ;
		case ELASTICITY_GENERALIZED_SELF_CONSISTENT:
			return SchemeTemplate::elasticityGeneralizedSelfConsistentScheme(d) ;
		case EXPANSION_HOBBS:
			return SchemeTemplate::expansionHobbsScheme(d) ;
		case EXPANSION_KERNER:
			return SchemeTemplate::expansionKernerScheme(d) ;
		case EXPANSION_TURNER:
			return SchemeTemplate::expansionTurnerScheme(d) ;

	}
	return empty ;
}
	
std::vector<double> SchemeTemplate::elasticityDilutedScheme(std::vector<std::vector<double> > d)
{
	std::vector<double> processed ;

	double k_mat = d[0][0] ;
	double mu_mat = d[0][1] ;

	double f_inc = d[1][0] ;
	double k_inc = d[1][1] ;
	double mu_inc = d[1][2] ;

	double k_hom = k_mat + 
		f_inc * (k_inc - k_mat) * (3*k_mat + 4*mu_mat) / 
		(3*k_inc + 4*mu_mat) ;
	double mu_hom = mu_mat + 
		f_inc * (mu_inc - mu_mat) * (5*mu_mat*(3*k_mat + 4*mu_mat)) / 
		(mu_mat*(9*k_mat + 8*mu_mat) + 6*mu_inc*(k_mat + 2*mu_mat)) ;

	processed.push_back(k_hom) ;
	processed.push_back(mu_hom) ;

	return processed ;
}

std::vector<double> SchemeTemplate::elasticityGeneralizedDilutedScheme(std::vector<std::vector<double> > d) 
{
	std::vector<double> processed ;
	
	double k_mat = d[0][0] ;
	double mu_mat = d[0][1] ;

	double k_hom = k_mat ;
	double mu_hom = mu_mat ;

	for(size_t i = 1 ; i < d.size() ; i++)
	{
		double f_inc = d[i][0] ;
		double k_inc = d[i][1] ;
		double mu_inc = d [i][2] ;

		k_hom += f_inc * (k_inc - k_mat) * (3*k_mat + 4*mu_mat) / 
			(3*k_inc + 4*mu_mat) ;
		mu_hom += f_inc * (mu_inc - mu_mat) * (5*mu_mat*(3*k_mat + 4*mu_mat)) / 
			(mu_mat*(9*k_mat + 8*mu_mat) + 6*mu_inc*(k_mat + 2*mu_mat)) ;
	}

	processed.push_back(k_hom) ;
	processed.push_back(mu_hom) ;

	return processed ;
}

std::vector<double> SchemeTemplate::elasticityIncrementalScheme(std::vector<std::vector<double> > d, double dalpha) 
{
	std::vector<double> processed ;

	double k_mat = d[0][0] ;
	double mu_mat = d[0][1] ;

	double f_inc = d[1][0] ;
	double k_inc = d[1][1] ;
	double mu_inc = d[1][2] ;

	double k_hom = k_mat ;
	double mu_hom = mu_mat ;

	double f = 0 ;
	while(f < f_inc)
	{
		k_hom += (k_inc-k_hom) * ((3*k_hom + 4*mu_hom) / 
			(3*k_inc + 4*mu_hom)) * dalpha/(1-f) ;
		mu_hom += (mu_inc - mu_hom) * (5*mu_hom*(3*k_hom + 4*mu_hom) / 
			(mu_hom*(9*k_hom+8*mu_hom) + 6*mu_inc*(k_hom+2*mu_hom))) * dalpha/(1-f) ;
		f += dalpha ;
	}

	double dlast = d[1][0] - (f-dalpha) ;
	k_hom += (k_inc-k_hom) * ((3*k_hom + 4*mu_hom) / 
		(3*k_inc + 4*mu_hom)) * dlast/(1-f) ;
	mu_hom += (mu_inc - mu_hom) * (5*mu_hom*(3*k_hom + 4*mu_hom) / 
		(mu_hom*(9*k_hom+8*mu_hom) + 6*mu_inc*(k_hom+2*mu_hom))) * dlast/(1-f) ;

	processed.push_back(k_hom) ;
	processed.push_back(mu_hom) ;

	return processed ;
}

std::vector<double> SchemeTemplate::elasticityMoriTanakaScheme(std::vector<std::vector<double> > d) 
{
	std::vector<double> processed ;

	double k_mat = d[0][0] ;
	double mu_mat = d[0][1] ;

	double k_hom = k_mat ;
	double mu_hom = mu_mat  ;

	double f_inc = d[1][0] ;
	double k_inc = d[1][1] ;
	double mu_inc = d[1][2] ;

	double Sk = (4*mu_mat) / (3*k_mat) ;
	double Smu = (3 + 2*Sk) / (2 + 3*Sk) ;
	double K = k_inc / k_mat ;
	double Mu = mu_inc / mu_mat ;

	k_hom *= (K/Sk + 1 + f_inc*(K-1)) ;
	k_hom /= (1 + (K + f_inc*(1-K))/Sk) ;
	mu_hom *= (Mu/Smu + 1 + f_inc*(Mu-1)) ;
	mu_hom /= (1 + (Mu + f_inc*(1-Mu))/Smu) ;

	processed.push_back(k_hom) ;
	processed.push_back(mu_hom) ;

	return processed ;
}

std::vector<double> SchemeTemplate::elasticityGeneralizedMoriTanakaScheme(std::vector<std::vector<double> > d) 
{
	std::vector<double> processed ;

	double k_mat = d[0][0] ;
	double mu_mat = d[0][1] ;

	double f_mat = 1. ;
	for(size_t i = 1 ; i < d.size() ; i++)
	{
		f_mat -= d[i][0] ;
	}

	double k_hom = k_mat * f_mat ;
	double mu_hom = mu_mat * f_mat ;

	std::vector<std::vector<double> > iter ;
	std::vector<double> mat ;
	std::vector<double> inc ;
	std::vector<double> mt_results ;

	for(size_t i = 1 ; i < d.size() ; i++)
	{
		iter.clear() ;
		mat.clear() ;
		inc.clear() ;
		
		mat.push_back(d[0][0]) ;
		mat.push_back(d[0][1]) ;

		inc.push_back(d[i][0]) ;
		inc.push_back(d[i][1]) ;
		inc.push_back(d[i][0]) ;
		
		iter.push_back(mat) ;
		iter.push_back(inc) ;

		mt_results = SchemeTemplate::elasticityMoriTanakaScheme(iter) ;

		k_hom += d[i][0] * (mt_results[0]) ;
		mu_hom += d[i][0] * (mt_results[1]) ;
	}

	processed.push_back(k_hom) ;
	processed.push_back(mu_hom) ;

	return processed ;
}

std::vector<double> SchemeTemplate::elasticitySelfConsistentScheme(std::vector<std::vector<double> > d) 
{
	std::vector<double> processed ;

	double f_mat = 1. - d[1][0] ;
	double & k_mat = d[0][0] ;
	double & mu_mat = d[0][1] ;
	double & f_inc = d[1][0] ;
	double & k_inc = d[1][1] ;
	double & mu_inc = d[1][2] ;

	std::vector<double> data_hom ;
	data_hom.push_back(k_mat*f_mat+k_inc*f_inc) ;
	data_hom.push_back(mu_mat*mu_mat+mu_inc*mu_inc) ;
	double & k_hom = data_hom[0] ;
	double & mu_hom = data_hom[1] ;


	double error = 1. ;
	std::vector<double> gmt_mat ;
	gmt_mat.push_back(0) ;
	gmt_mat.push_back(0) ;
	std::vector<double> gmt_inc ;
	gmt_inc.push_back(0) ;
	gmt_inc.push_back(0) ;

	std::vector<std::vector<double> > data_mat ;
	std::vector<std::vector<double> > data_inc ;
	{
		std::vector<double> empt2 ;
		empt2.push_back(0) ;
		empt2.push_back(0) ;
		std::vector<double> empt3 ;
		empt3.push_back(0) ;
		empt3.push_back(0) ;
		empt3.push_back(0) ;
		data_mat.push_back(empt2) ;
		data_mat.push_back(empt3) ;
		data_inc.push_back(empt2) ;
		data_inc.push_back(empt3) ;
	}

	double f_mat_hom = 1. - data_mat[1][0] ;
	double & k_mat_hom = data_mat[0][0] ;
	double & mu_mat_hom = data_mat[0][1] ;
	double f_inc_hom = 1. - data_inc[1][0] ;
	double & k_inc_hom = data_inc[0][0] ;
	double & mu_inc_hom = data_inc[0][1] ;

	data_mat[0][0] = k_hom ;
	data_inc[0][0] = k_hom ;
	data_mat[0][1] = mu_hom ;
	data_inc[0][1] = mu_hom ;

	data_mat[1][0] = f_mat ;
	data_inc[1][0] = f_inc ;
	data_mat[1][1] = k_mat ;
	data_inc[1][1] = k_inc ;
	data_mat[1][2] = mu_mat ;
	data_inc[1][2] = mu_inc ;

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
		gmt_mat = SchemeTemplate::elasticityMoriTanakaScheme(data_mat) ;
		gmt_inc = SchemeTemplate::elasticityMoriTanakaScheme(data_inc) ;

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


	processed.push_back(k_mat_hom) ;
	processed.push_back(mu_mat_hom) ;

	return processed ;
}

std::vector<double> SchemeTemplate::elasticityGeneralizedSelfConsistentScheme(std::vector<std::vector<double> > d) 
{
	std::vector<double> processed ;

	std::vector<double> kmu_mean ;
	kmu_mean.push_back(0.) ;
	kmu_mean.push_back(0.) ;
	double f_mat = 1. ;
	for(size_t i = 1 ; i < d.size() ; i++)
	{
		kmu_mean[0] += d[i][1] * d[i][0] ;
		kmu_mean[1] += d[i][2] * d[i][0] ;
		f_mat -= d[i][0] ;
	}
	kmu_mean[0] += d[0][0] * f_mat ;
	kmu_mean[1] += d[0][1] * f_mat ;
	
	
	std::vector<double> hom ;
	hom.push_back(kmu_mean[0]) ;
	hom.push_back(kmu_mean[1]) ;
	std::vector<double> next_hom ;
	next_hom.push_back(kmu_mean[0]) ;
	next_hom.push_back(kmu_mean[1]) ;
	
	double error = 1. ;
	int nit = 0 ;
	
	std::vector<std::vector<double> > iscm ;
	std::vector<std::vector<double> > isc ;
	{
		std::vector<double> empt2 ;
		empt2.push_back(0) ;
		empt2.push_back(0) ;
		std::vector<double> empt3 ;
		empt3.push_back(0) ;
		empt3.push_back(0) ;
		empt3.push_back(0) ;
		iscm.push_back(empt2) ;
		iscm.push_back(empt3) ;
		isc.push_back(empt2) ;
		isc.push_back(empt3) ;
	}
	std::vector<double> mt_mat_results ;
	mt_mat_results.push_back(0) ;
	mt_mat_results.push_back(0) ;
	std::vector<double> mt_results ;
	mt_results.push_back(0) ;
	mt_results.push_back(0) ;
	while(error > 1e-6)
	{
		error = 1. ;
		std::vector<double> dk ;
		std::vector<double> dmu ;

		iscm[0][0] = hom[0] ;
		iscm[0][1] = hom[1] ;
		iscm[1][0] = f_mat ;
		iscm[1][1] = d[0][0] ;
		iscm[1][2] = d[0][1] ;

		mt_mat_results = SchemeTemplate::elasticityMoriTanakaScheme(iscm) ;

		double fm = iscm[1][0] ;

		double khomm = iscm[0][0] ;
		double km = iscm[1][1] ;
		double kmtm = mt_mat_results[0] ;

		double muhomm = iscm[0][1] ;
		double mum = iscm[1][2] ;
		double mumtm = mt_mat_results[1] ;


		double Akm = (kmtm-khomm)/(km-khomm)/fm ;
		double Amum = (mumtm-muhomm)/(mum-muhomm)/fm ;

		dk.push_back(d[0][0]) ;
		dmu.push_back(d[0][1]) ;

		for(size_t i = 1 ; i < d.size() ; i++)
		{
			isc[0][0] = hom[0] ;
			isc[0][1] = hom[1] ;
			isc[1][0] = d[i][0] ;
			isc[1][1] = d[i][1] ;
			isc[1][2] = d[i][2] ;

			
			mt_results = SchemeTemplate::elasticityMoriTanakaScheme(isc) ;

			double fi = isc[1][0] ;

			double khom = isc[0][0] ;
			double ki = isc[1][1] ;
			double kmt = mt_results[0] ;

			double muhom = isc[0][1] ;
			double mui = isc[1][2] ;
			double mumt = mt_results[1] ;


			double Ak = (kmt-khom)/(ki-khom)/fi ;
			double Amu = (mumt-muhom)/(mui-muhom)/fi ;
			
			dk.push_back(fi*(ki-km)*Ak/(fi*Ak+(1.-fi)*Akm)) ;
			dmu.push_back(fi*(mui-mum)*Amu/(fi*Amu+(1.-fi)*Amum)) ;

		}

		next_hom[0] = 0. ;
		next_hom[1] = 0. ;

		for(size_t i = 0 ; i < dk.size() ; i++)
		{
			next_hom[0] += dk[i] ;
			next_hom[1] += dmu[i] ;
		}

		error = std::max(std::abs(hom[0]-next_hom[0])/hom[0],std::abs(hom[1]-next_hom[1])/hom[1]) ;

		hom[0] = next_hom[0] ;
		hom[1] = next_hom[1] ;
		
		std::cout << hom[0] << ";" << hom[1] << std::endl ;

	}

	processed.push_back(hom[0]) ;
	processed.push_back(hom[1]) ;

	return processed ;
}

std::vector<double> SchemeTemplate::expansionHobbsScheme(std::vector<std::vector<double> > d)
{
	double kmat = d[0][0] ;
	double amat = d[0][1] ;
	
	double finc = d[1][0] ;
	double kinc = d[1][1] ;
	double ainc = d[1][2] ;
	
	std::vector<double> processed ;
	processed.push_back(amat + 2*finc*kinc*(ainc-amat) / (kmat+kinc+finc*(kinc-kmat))) ;

	return processed ;

}

std::vector<double> SchemeTemplate::expansionKernerScheme(std::vector<std::vector<double> > d)
{
	double fmat = 1 - d[1][0] ;
	double kmat = d[0][0] ;
	double amat = d[0][1] ;
	
	double finc = d[1][0] ;
	double kinc = d[1][1] ;
	double ainc = d[1][2] ;
	
	std::vector<double> processed ;
	processed.push_back((amat*kmat*fmat + ainc*kinc*finc) / (kmat*fmat + kinc*finc)) ;

	return processed ;

}

std::vector<double> SchemeTemplate::expansionTurnerScheme(std::vector<std::vector<double> > d)
{
	double fmat = d[0][0] ;
	double kmat = d[0][1] ;
	double mumat= d[0][2] ;
	double amat = d[0][3] ;
	
	double finc = 1 - d[0][0] ;
	double kinc = d[1][0] ;
	double ainc = d[1][1] ;
	
	std::vector<double> processed ;
	processed.push_back(amat*fmat + ainc*finc) ;
	processed[0] += (fmat*finc)*(ainc - amat)*(kinc-kmat)/(kmat*fmat+kinc*finc+3*kmat*kinc/(4*mumat)) ;

	return processed ;

}



