#ifndef SCHEME_TEMPLATE_H
#define SCHEME_TEMPLATE_H

#include "homogenization_base.h"

namespace Mu
{

class PhaseTemplate
{
	std::vector<Properties> prop ;
	int max ;
	int done ;
	
public:
	PhaseTemplate(int n = 1) ;
	PhaseTemplate(int n, std::vector<PType> t) ;
	
	bool cast(Material m) ;
	bool cast(int i, Properties p) ;

	std::vector<double> val() const ;

	bool hasProperties(PType t) const ;
	int indexOfProperties(PType t) const ;	
	int size() const ;
	int phases(int count = -1) const ;
	bool filled() const ;
	void force(std::vector<double> data) ;
	Properties get(int i) const ;
	
	static PhaseTemplate makeBulkShearTemplate(int n = -1) ;
	static PhaseTemplate makeVolumeBulkShearTemplate(int n = -1) ;

} ;

class SchemeTemplate
{
	HomogenizationScheme scheme ;
	std::vector<PhaseTemplate> phases ;
	PhaseTemplate result ;
	std::vector<std::vector<double> > data ;
	
public:
	SchemeTemplate(PhaseTemplate r, std::vector<PhaseTemplate> p) ;
	SchemeTemplate(HomogenizationScheme s) ;
	
	bool cast(Material m) ;
	
	int sizeResult() const ;
	Properties getResult(int i) const ;

	static std::vector<double> makeScheme(HomogenizationScheme s, std::vector<std::vector<double> > d) ;	
	
	static std::vector<double> elasticityDilutedScheme(std::vector<std::vector<double> > d) ;
	static std::vector<double> elasticityGeneralizedDilutedScheme(std::vector<std::vector<double> > d) ;
	static std::vector<double> elasticityIncrementalScheme(std::vector<std::vector<double> > d, double dalpha = 1e-5) ;
	static std::vector<double> elasticityMoriTanakaScheme(std::vector<std::vector<double> > d) ;
	static std::vector<double> elasticityGeneralizedMoriTanakaScheme(std::vector<std::vector<double> > d) ;
	static std::vector<double> elasticitySelfConsistentScheme(std::vector<std::vector<double> > d) ;
	static std::vector<double> elasticityGeneralizedSelfConsistentScheme(std::vector<std::vector<double> > d) ;
	

} ;
	
	
} ;

#endif
