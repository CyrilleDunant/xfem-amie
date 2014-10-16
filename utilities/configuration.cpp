
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//

#include "configuration.h"
#include "../physics/stiffness.h"
#include "../physics/stiffness_and_fracture.h"
#include "../physics/stiffness_with_imposed_deformation.h"
#include "../physics/logarithmic_creep_with_external_parameters.h"
#include "../physics/viscoelasticity.h"
#include "../physics/viscoelasticity_and_fracture.h"
#include "../physics/viscoelasticity_and_imposed_deformation.h"
#include "../physics/material_laws/temperature_material_laws.h"
#include "../physics/materials/aggregate_behaviour.h"
#include "../physics/materials/paste_behaviour.h"
#include "../physics/materials/gel_behaviour.h"
#include "../physics/materials/rebar_behaviour.h"
#include "../physics/materials/steel_behaviour.h"
#include "../physics/materials/concrete_behaviour.h"
#include "../physics/damagemodels/spacetimefiberbasedisotropiclineardamage.h"
#include "../physics/damagemodels/isotropiclineardamage.h"
#include "../physics/damagemodels/plasticstrain.h"
#include "../physics/fracturecriteria/maxstrain.h"
#include "../physics/fracturecriteria/mohrcoulomb.h"
#include "../physics/fracturecriteria/nonlocalvonmises.h"
#include "writer/triangle_writer.h"
#include "iostream"
#include "fstream"


using namespace Amie ;

ConfigTreeItem::ConfigTreeItem() : father(nullptr), label("trunk"), data(0.), str("")
{

}

ConfigTreeItem::ConfigTreeItem(ConfigTreeItem * f, std::string l) : father(f), label(l), data(0.), str("")
{
	if(f)
		f->addChild(this) ;
}

ConfigTreeItem::ConfigTreeItem(ConfigTreeItem * f, std::string l, double d) : father(f), label(l), data(d), str("")
{
	if(f)
		f->addChild(this) ;
}

ConfigTreeItem::ConfigTreeItem(ConfigTreeItem * f, std::string l, std::string s) : father(f), label(l), data(0.), str(s)
{
	if(f)
		f->addChild(this) ;
}

ConfigTreeItem * ConfigTreeItem::getFather() const 
{
	return father ;
}

ConfigTreeItem * ConfigTreeItem::getChild(std::string childLabel)  const
{
	for(size_t i = 0 ; i < children.size() ; i++)
	{
		if(children[i]->is(childLabel))
			return children[i] ;
	}
	return nullptr ;
}

ConfigTreeItem * ConfigTreeItem::getLastChild() const
{
	if(children.size() > 0)
		return children[children.size()-1] ;
	return nullptr ;
}

std::vector<ConfigTreeItem *> ConfigTreeItem::getAllChildren(std::string childLabel) const
{
	std::vector<ConfigTreeItem *> ret ;
	for(size_t i = 0 ; i < children.size() ; i++)
	{
		if(children[i]->is(childLabel))
			ret.push_back(children[i]) ; 
	}
	return ret ;
}

std::vector<ConfigTreeItem *> ConfigTreeItem::getAllChildren() const
{
	return children ;
}

ConfigTreeItem * ConfigTreeItem::getChild(std::vector<std::string> childLabelDecomposed) const
{
	if(childLabelDecomposed.size() == 0)
		return nullptr ;
	if(childLabelDecomposed.size() == 1)
		return getChild(childLabelDecomposed[0]) ;

	ConfigTreeItem * child = getChild(childLabelDecomposed[0]) ;
	if(child)
	{
		std::vector<std::string> newDecomposition ;
		for(size_t i = 1 ; i < childLabelDecomposed.size() ; i++)
		{
			newDecomposition.push_back(childLabelDecomposed[i]) ;
		}
		return child->getChild( newDecomposition ) ;
	}
	return nullptr ;
}

ConfigTreeItem * ConfigTreeItem::getChildFromFullLabel(std::string childLabel)  const
{
	return getChild( ConfigTreeItem::decompose( childLabel ) ) ;
}

void ConfigTreeItem::addChild(ConfigTreeItem * c) 
{
	if(c)
	{
		children.push_back(c) ;
		c->setFather(this) ;
	}
}

void ConfigTreeItem::setFather(ConfigTreeItem * f) 
{
	if(father && f)
		father = f ;
}

ConfigTreeItem * ConfigTreeItem::getRoot() const 
{
	if(isTrunk())
		return nullptr ;
	ConfigTreeItem * root = this->getFather() ;
	while(!root->isTrunk())
		root = root->getFather() ;
	return root ;
}

bool ConfigTreeItem::is(std::string l) const { return label == l ; }

bool ConfigTreeItem::isTrunk() const { return father == nullptr ; }

bool ConfigTreeItem::isLeaf() const { return children.size() == 0 ; }

std::string ConfigTreeItem::getFullLabel() const 
{
	if(isTrunk())
		return label ;

	std::string path = father->getFullLabel() ;
	path.append(".") ;
	path.append(label) ;

	return path ;
}

void ConfigTreeItem::setLabel(std::string l) 
{
	if(isTrunk())
		return ;

	label = l ;
}
std::vector<std::string> ConfigTreeItem::decompose(std::string full)
{
	std::vector<std::string> decomposition ;
	std::string current = full ;
	size_t found = current.find(".") ;
	while(found != std::string::npos)
	{
		decomposition.push_back(current.substr(0, found)) ;
		current = current.substr(found+1) ;
		found = current.find(".") ;
	}
	decomposition.push_back(current) ;
	return decomposition ;
}

void ConfigTreeItem::print() const
{
	std::cout << getFullLabel() ;
	if(isLeaf())
	{
		std::cout << " = " ;
		if(str.size() > 0)
			std::cout << str ;
		else
			std::cout << data ;
	}
	std::cout << std::endl ;
} ;

void ConfigTreeItem::printTree() const
{
	print() ;
	for(size_t i = 0 ; i < children.size() ; i++)
		children[i]->printTree() ;
} ;

bool ConfigTreeItem::hasChild( std::string l )  const
{
	for(size_t i = 0 ; i < children.size() ; i++)
	{
		if(children[i]->is(l))
			return true ;
	}
	return false ;
}

bool ConfigTreeItem::hasChildFromFullLabel( std::string childLabel )  const
{
	return hasChild( ConfigTreeItem::decompose(childLabel) ) ;
}

bool ConfigTreeItem::hasChild( std::vector<std::string> childLabelDecomposed )  const
{
	if(childLabelDecomposed.size() == 0)
		return false ;
	if(childLabelDecomposed.size() == 1)
		return hasChild(childLabelDecomposed[0]) ;

	if(hasChild(childLabelDecomposed[0]))
	{
		std::vector<std::string> newDecomposition ;
		for(size_t i = 1 ; i < childLabelDecomposed.size() ; i++)
		{
			newDecomposition.push_back(childLabelDecomposed[i]) ;
		}
		return getChild(childLabelDecomposed[0])->hasChild( newDecomposition ) ;

	}
	else
		return false ;
	
}

double ConfigTreeItem::getData(std::string path, double defaultValue) const 
{
	if(hasChildFromFullLabel(path))
		return getChildFromFullLabel(path)->getData() ;
	return defaultValue ;
}

std::string ConfigTreeItem::getStringData(std::string path, std::string defaultValue) const 
{
	if(hasChildFromFullLabel(path))
		return getChildFromFullLabel(path)->getStringData() ;
	return defaultValue ;
}

Sample * ConfigTreeItem::getSample() const
{
	bool spaceTime = (ConfigTreeItem::translateOrder(getRoot()->getStringData("discretization.order","LINEAR")) >= CONSTANT_TIME_LINEAR) ;
	double sampleWidth = getData("width",1) ;
	double sampleHeight = getData("height",1) ;
	double sampleCenterX = getData("center.x",0) ;
	double sampleCenterY = getData("center.y",0) ;
	Sample * ret = new Sample(nullptr, sampleWidth, sampleHeight, sampleCenterX, sampleCenterY) ;
	if(hasChild("behaviour"))
	{
		ret->setBehaviour( getChild("behaviour")->getBehaviour( SPACE_TWO_DIMENSIONAL, spaceTime ) ) ;
	}

	return ret ;
}

Matrix ConfigTreeItem::getStiffnessMatrix(SpaceDimensionality dim, planeType pt) const 
{
	if(hasChild("bulk_modulus"))
	{
		double k = getData("bulk_modulus",1e9) ;
		double mu = getData("shear_modulus", 1e9) ;
	        return Material::cauchyGreen( k, mu, false, dim, pt) ;
	} else {
		double E = getData("young_modulus",1e9) ;
		double nu = getData("poisson_ratio", 0.2) ;
	        return Material::cauchyGreen( E, nu, true, dim, pt) ;
	}
}

Vector ConfigTreeItem::getImposedDeformation(SpaceDimensionality dim) const 
{
	double alpha = getData("imposed_deformation", 0.) ;
	Vector imp(3+3*(dim == SPACE_THREE_DIMENSIONAL)) ;
	imp = 0 ;
	for(size_t i = 0 ; i < 2+(dim == SPACE_THREE_DIMENSIONAL) ; i++)
		imp[i] = alpha ;
	return imp ;
}

Function ConfigTreeItem::getFunction() const 
{
	std::string function = getStringData() ;
	Function f(function.c_str()) ;
	return f ;
}

ExternalMaterialLaw * ConfigTreeItem::getExternalMaterialLaw() const 
{
	ExternalMaterialLaw * ret = nullptr ;
	std::string type = getStringData("type","ABSTRACT") ;
	if(type == "ABSTRACT")
		return ret ;

	if(type == "THERMAL_EXPANSION")
	{
		ret = new ThermalExpansionMaterialLaw("temperature = 293") ;
		ret->setDefaultValue("temperature", getData("reference_temperature", 293)) ;
		return ret ;
	}

	if(type == "RADIATION_INDUCED_VOLUMETRIC_EXPANSION")
	{
		return new RadiationInducedExpansionMaterialLaw() ;
	}

	if(type == "DRYING_SHRINKAGE")
	{
		ret = new RadiationInducedExpansionMaterialLaw() ;
		ret->setDefaultValue("relative_humidity", getData("reference_relative_humidity", 1.)) ;
		return ret ;
	}

	if(type == "ARRHENIUS")
	{
		std::string affected = getStringData("parameter_affected","none") ;
		ret = new ArrheniusMaterialLaw(affected, "temperature = 293") ;
		ret->setDefaultValue("temperature", getData("reference_temperature", 293)) ;
		return ret ;
	}

	if(type == "CREEP_ARRHENIUS")
	{
		ret = new CreepArrheniusMaterialLaw("temperature = 293") ;
		ret->setDefaultValue("temperature", getData("reference_temperature", 293)) ;
		return ret ;
	}

	if(type == "CREEP_HUMIDITY")
	{
		ret = new CreepRelativeHumidityMaterialLaw() ;
		ret->setDefaultValue("creep_humidity_coefficient",5.) ;
		return ret ;
	}

	if(type == "CONSTANT")
	{
		ret = new ConstantExternalMaterialLaw(std::string()) ;
		std::vector<ConfigTreeItem *> parameters = getAllChildren() ;
		for(size_t i = 1 ; i < parameters.size() ; i++)
		{
			ret->setDefaultValue( parameters[i]->getLabel(), parameters[i]->getData() ) ;
		}
		return ret ;
	}

	if(type == "SPACE_TIME_DEPENDENT")
	{
		if(!hasChild("function"))
		{
			std::cout << "no function found while creating SpaceTimeDependentMaterialLaw" << std::endl ;
			return nullptr ;
		}
		Function f = getChild("function")->getFunction() ;
		ret = new SpaceTimeDependentExternalMaterialLaw( getStringData("output_parameter","none"), f, getStringData("additive", "FALSE") == "TRUE") ;
		return ret ;
	}

	if(type == "SIMPLE_DEPENDENT")
	{
		if(!hasChild("function"))
		{
			std::cout << "no function found while creating SimpleDependentExternalMaterialLaw" << std::endl ;
			return nullptr ;
		}
		Function f = getChild("function")->getFunction() ;
		ret = new SimpleDependentExternalMaterialLaw( getStringData("output_parameter","none"), getStringData("input_parameter","x"), f, getStringData("additive", "FALSE") == "TRUE") ;
		return ret ;
	}

	if(type == "VARIABLE_DEPENDENT")
	{
		if(!hasChild("function"))
		{
			std::cout << "no function found while creating SimpleDependentExternalMaterialLaw" << std::endl ;
			return nullptr ;
		}
		Function f = getChild("function")->getFunction() ;
		ret = new VariableDependentExternalMaterialLaw( getStringData("output_parameter","none"), f, getStringData("additive", "FALSE") == "TRUE") ;
		if(hasChild("x"))
			dynamic_cast<VariableDependentExternalMaterialLaw*>(ret)->setAsX( getStringData("x","x")) ;
		if(hasChild("y"))
			dynamic_cast<VariableDependentExternalMaterialLaw*>(ret)->setAsY( getStringData("y","y")) ;
		if(hasChild("z"))
			dynamic_cast<VariableDependentExternalMaterialLaw*>(ret)->setAsZ( getStringData("z","z")) ;
		if(hasChild("t"))
			dynamic_cast<VariableDependentExternalMaterialLaw*>(ret)->setAsT( getStringData("t","t")) ;
		if(hasChild("u"))
			dynamic_cast<VariableDependentExternalMaterialLaw*>(ret)->setAsU( getStringData("u","u")) ;
		if(hasChild("v"))
			dynamic_cast<VariableDependentExternalMaterialLaw*>(ret)->setAsV( getStringData("v","v")) ;
		if(hasChild("w"))
			dynamic_cast<VariableDependentExternalMaterialLaw*>(ret)->setAsW( getStringData("w","w")) ;
		return ret ;
	}

	if(type == "LINEAR_INTERPOLATED")
	{
		ret = new LinearInterpolatedExternalMaterialLaw( std::make_pair(getStringData("input_parameter","t"),getStringData("output_parameter","none")), getStringData("file_name","file_name")) ;
		return ret ;
	}

	if(type == "TIME_DERIVATIVE")
	{
		ret = new TimeDerivativeMaterialLaw( getStringData("input_parameter","none"), getStringData("output_parameter", getStringData("input_parameter","none")+"_rate")) ;
		return ret ;
	}

	if(type == "TIME_INTEGRAL")
	{
		ret = new TimeIntegralMaterialLaw( getStringData("input_parameter","none"), getStringData("output_parameter", getStringData("input_parameter","none")+"_integral")) ;
		return ret ;
	}

	if(type == "MINIMUM" || type == "MAXIMUM")
	{
		std::vector<ConfigTreeItem *> all = getAllChildren("input_parameter") ;
		std::vector<std::string> coord ;
		for(size_t i = 0 ; i < all.size() ; i++)
			coord.push_back( all[i]->getStringData() ) ;
		bool add = getStringData("additive", "FALSE") == "TRUE" ;
		if(type == "MINIMUM")
			ret = new MinimumMaterialLaw(  getStringData("output_parameter","none"), coord, add ) ;
		else
			ret = new MaximumMaterialLaw(  getStringData("output_parameter","none"), coord, add ) ;
	}

	if(type == "GET_FIELD")
	{
		bool isFieldType = false ;
		std::string field =  getStringData("field", "NOT_A_FIELD") ;
		FieldType f = ConfigTreeItem::translateFieldType( getStringData("field", "NOT_A_FIELD"), isFieldType ) ;
		if(isFieldType)
		{
			std::transform(field.begin(), field.end(), field.begin(), ::tolower);
			ret = new GetFieldMaterialLaw( f, field ) ;
		}
	}

	return ret ;
}

FractureCriterion * ConfigTreeItem::getFractureCriterion(bool spaceTime) const 
{
	FractureCriterion * ret = nullptr ;
	std::string type = getStringData("type","ABSTRACT") ;
	if(spaceTime)
	{
		if(type == "MAXIMUM_TENSILE_STRAIN")
		{
			ret = new SpaceTimeNonLocalMaximumStrain( getData("limit_tensile_strain",0.001) ) ;
		}
		if(type == "MAXIMUM_TENSILE_STRESS")
		{
			ret = new SpaceTimeNonLocalMaximumStress( getData( "limit_tensile_stress",1e6) ) ;
		}
		if(type == "LINEAR_SOFTENING_MAXIMUM_TENSILE_STRAIN")
		{
			ret = new SpaceTimeNonLocalLinearSofteningMaximumStrain( getData("limit_tensile_strain",0.001), getData("limit_tensile_stress",1e6), getData("maximum_tensile_strain",0.002) ) ;
		}
		if(type == "ELLIPSOIDAL_SOFTENING_MAXIMUM_TENSILE_STRESS")
		{
			ret = new SpaceTimeNonLocalEllipsoidalMixedCriterion( getData("limit_tensile_strain",0.001), getData("limit_tensile_stress",1e6), getData("instantaneous_modulus",10e9), getData( "relaxed_modulus",1e9) ) ;
		}
	}
	else
	{
		if(type == "MOHR_COULOMB")
		{
			ret = new NonLocalMohrCoulomb( getData("limit_tensile_strain",0.001), getData("limit_compressive_strain",-0.001), getFather()->getData("young_modulus",1e9) ) ;
		}
		if(type == "LINEAR_SOFTENING_MOHR_COULOMB")
		{
			ret = new NonLocalLinearlyDecreasingMohrCoulomb( getData("limit_tensile_strain",0.001), getData("limit_compressive_strain",-0.001), getData("maximum_tensile_strain",0.002), getData("maximum_compressive_strain",-0.002), getFather()->getData("young_modulus",1e9) ) ;
		}
		if(type == "EXPONENTIAL_SOFTENING_MOHR_COULOMB")
		{
			ret = new NonLocalExponentiallyDecreasingMohrCoulomb( getData("limit_tensile_strain",0.001), getData("limit_compressive_strain",-0.001), getData("maximum_tensile_strain",0.002), getData("maximum_compressive_strain",-0.002), getFather()->getData("young_modulus",1e9) ) ;
		}
		if(type == "MCFT")
		{
			ret = new NonLocalMCFT( getData("limit_compressive_strain",-0.001), getFather()->getData("young_modulus",1e9), getData("material_characteristic_radius",0.001)) ;
			std::vector<ConfigTreeItem *> rebarItems = getAllChildren("rebar") ;
			std::vector<std::pair<double, double> > rebars ;
			for(size_t i = 0 ; i < rebarItems.size() ; i++)
			{
				rebars.push_back( std::make_pair ( rebarItems[i]->getData("location",0.001), rebarItems[i]->getData("diameter",0.0001) ) ) ;
			}
		}
		if(type == "VON_MISES")
		{
			ret = new NonLocalVonMises( getData("limit_tensile_stress",1e6),  getFather()->getData("young_modulus",1e9), getData("material_characteristic_radius",0.001) ) ;
		}
	}
	if(ret && hasChild("material_characteristic_radius"))
		ret->setMaterialCharacteristicRadius(getData("material_characteristic_radius",0.001)) ;

	return ret ;
}

DamageModel * ConfigTreeItem::getDamageModel(bool spaceTime) const 
{
	DamageModel * ret = nullptr ;
	std::string type = getStringData("type","ABSTRACT") ;
	if(type == "ISOTROPIC_INCREMENTAL_LINEAR_DAMAGE")
	{
		if(spaceTime)
			ret = new SpaceTimeFiberBasedIsotropicLinearDamage( getData("damage_increment",0.1), getData("time_tolerance",1e-9) ) ;
		else
			ret = new FiberBasedIsotropicLinearDamage( getData("damage_increment",0.1) ) ;
	}
	if(type == "ISOTROPIC_LINEAR_DAMAGE")
	{
		if(spaceTime)
			std::cout << "warning: damage model incompatible with finite element discretization" << std::endl ;
		ret = new IsotropicLinearDamage() ;
	}
	if(type == "PLASTIC_STRAIN")
	{
		if(spaceTime)
			std::cout << "warning: damage model incompatible with finite element discretization" << std::endl ;
		ret = new PlasticStrain() ;
	}

	if(ret)
	{
		if(hasChild("maximum_damage"))
			ret->setThresholdDamageDensity(getData("maximum_damage",1.)) ;
	}

	return ret ;
}

Form * ConfigTreeItem::getBehaviour(SpaceDimensionality dim, bool spaceTime) const
{
	std::string type = getStringData("type","VOID_BEHAVIOUR") ;

	if(type == std::string("VOID_BEHAVIOUR"))
		return new VoidForm() ;

	planeType pt = PLANE_STRESS ;
	if(hasChild("plane_type"))
		pt = ConfigTreeItem::translatePlaneType( getStringData("plane_type", "PLANE_STRESS") ) ;

	if(type == std::string("ELASTICITY"))
	{
		if(spaceTime)
		{
			int blocks = getData("additional_viscoelastic_variables",0) ;
			return new Viscoelasticity( PURE_ELASTICITY, getStiffnessMatrix(dim, pt), blocks) ;
		}
		return new Stiffness( getStiffnessMatrix(dim, pt) ) ;
	}

	if(type == std::string("ELASTICITY_AND_FRACTURE"))
	{
		FractureCriterion * frac = getChild("fracture_criterion")->getFractureCriterion(spaceTime) ;
		DamageModel * dam = getChild("damage_model")->getDamageModel(spaceTime) ;
		if(frac && dam)
		{
			if(spaceTime)
			{
				int blocks = getData("additional_viscoelastic_variables",0) ;
				return new ViscoelasticityAndFracture( PURE_ELASTICITY, getStiffnessMatrix(dim, pt), frac, dam, blocks) ;
			}
			return new StiffnessAndFracture( getStiffnessMatrix(dim, pt), frac, dam ) ;
		}
		else
		{
			if(!frac)
				std::cout << "fracture criterion undefined!" << std::endl ;
			if(!dam)
				std::cout << "damage model undefined!" << std::endl ;

			if(spaceTime)
			{
				int blocks = getData("additional_viscoelastic_variables",0) ;
				return new Viscoelasticity( PURE_ELASTICITY, getStiffnessMatrix(dim, pt), blocks) ;
			}
			return new Stiffness( getStiffnessMatrix(dim, pt) ) ;
		}
	}

	if(type == std::string("ELASTICITY_AND_IMPOSED_DEFORMATION"))
	{
		if(spaceTime)
		{
			int blocks = getData("additional_viscoelastic_variables",0) ;
			Vector imp = getImposedDeformation(dim) ;
			return new ViscoelasticityAndImposedDeformation( PURE_ELASTICITY, getStiffnessMatrix(dim, pt), imp, blocks) ;
		}
		return new StiffnessWithImposedDeformation( getStiffnessMatrix(dim, pt), getImposedDeformation(dim) ) ;
	}

	if(type == std::string("LOGARITHMIC_CREEP"))
			
	{
		if(!spaceTime)
			return nullptr ;


		LogarithmicCreepWithExternalParameters * log = nullptr;

		if(hasChild("fracture_criterion") && hasChild("damage_model"))
		{
			FractureCriterion * frac = getChild("fracture_criterion")->getFractureCriterion(true) ;
			DamageModel * dam = getChild("damage_model")->getDamageModel(true) ;
			if(frac && dam)
				log = new LogarithmicCreepWithExternalParameters(std::string(), frac, dam) ;
		}

		if(!log)
			log = new LogarithmicCreepWithExternalParameters(std::string()) ;

		if(!hasChild("parameters"))
		{
			std::cout << "log-creep law without parameters!" << std::endl ;
			return nullptr ;
		}

		std::vector<ConfigTreeItem *> parameters = getChild("parameters")->getAllChildren() ;
		for(size_t i = 0 ; i < parameters.size() ; i++)
		{
//			std::cout << parameters[i]->getLabel() << "\t" << parameters[i]->getData() << std::endl ;
			log->addMaterialParameter( parameters[i]->getLabel(), parameters[i]->getData() ) ;
		}

		std::vector<ConfigTreeItem *> laws = getAllChildren("material_law") ;
		for(size_t i = 0 ; i < laws.size() ; i++)
		{
			log->addMaterialLaw( laws[i]->getExternalMaterialLaw() ) ;
		}

//		std::cout << getFullLabel() << std::endl ;

		return log ;
	}


	if(type == std::string("PASTE_BEHAVIOUR"))
	{
		if(getStringData("damage","TRUE") == "FALSE")
		{
			if(spaceTime)
				return new ViscoElasticOnlyPasteBehaviour( getData("young_modulus", 12e9) , getData( "poisson_ratio", 0.3), getData("short_term_creep_modulus", 0.3), getData("long_term_creep_modulus", 0.37), dim) ;
			return new ElasticOnlyPasteBehaviour( getData("young_modulus", 12e9) , getData( "poisson_ratio", 0.3), dim) ;
		}
		if(spaceTime)
			return new ViscoDamagePasteBehaviour( getData("young_modulus", 12e9) , getData( "poisson_ratio", 0.3), getData( "tensile_strain_limit", 0.0003), 0.00025, dim) ;
		return new PasteBehaviour( getData("young_modulus", 12e9) , getData( "poisson_ratio", 0.3), getData( "tensile_strain_limit", 0.0016666), 0.001666, 12000, dim) ;
	}

	if(type == std::string("AGGREGATE_BEHAVIOUR"))
	{
		if(getStringData("damage","TRUE") == "FALSE")
		{
			if(spaceTime)
				return new ViscoElasticOnlyAggregateBehaviour( getData("young_modulus", 59e9) , getData( "poisson_ratio", 0.3), dim) ;
			return new ElasticOnlyAggregateBehaviour( getData("young_modulus", 59e9) , getData( "poisson_ratio", 0.3), dim) ;
		}
		if(spaceTime)
			return new ViscoDamageAggregateBehaviour( getData("young_modulus", 59e9) , getData( "poisson_ratio", 0.3), getData( "tensile_strain_limit", 0.0012), 0.00025, dim) ;
		return new AggregateBehaviour( getData("young_modulus", 59e9) , getData( "poisson_ratio", 0.3), getData( "tensile_strain_limit", 0.0012), 0.00044, 12000, dim) ;
	}

	if(type == std::string("ASR_GEL_BEHAVIOUR"))
	{
		if(spaceTime)
			return new ViscoElasticOnlyGelBehaviour( getData("young_modulus", 22e9) , getData( "poisson_ratio", 0.18), getData( "imposed_deformation", 0.22), dim) ;
		return new GelBehaviour( getData("young_modulus", 22e9) , getData( "poisson_ratio", 0.18), getData( "imposed_deformation", 0.22), dim) ;
	}

	if(type == std::string("CONCRETE_BEHAVIOUR"))
	{
		if(spaceTime)
		{
			std::cout << "cannot make a space-time concrete behaviour" << std::endl ;
			return nullptr ;
		}
		return new ConcreteBehaviour( getData("young_modulus", 37e9) , getData( "poisson_ratio", 0.3), getData( "compressive_strength", -37e6), PLANE_STRESS, UPPER_BOUND, dim) ;
	}

	if(type == std::string("REBAR_BEHAVIOUR"))
	{
		if(spaceTime)
		{
			std::cout << "cannot make a space-time rebar behaviour" << std::endl ;
			return nullptr ;
		}
		return new RebarBehaviour( getData("young_modulus", 59e9) , getData( "poisson_ratio", 0.3), getData( "tensile_strength", 57000), dim) ;
	}

	if(type == std::string("STEEL_BEHAVIOUR"))
	{
		if(spaceTime)
		{
			std::cout << "cannot make a space-time steel behaviour" << std::endl ;
			return nullptr ;
		}
		return new SteelBehaviour( getData("young_modulus", 210e9) , getData( "poisson_ratio", 0.3), getData( "tensile_strength", 276e6), dim) ;
	}

	if(type == std::string("VISCOSITY"))
	{
		int blocks = getData("additional_viscoelastic_variables",0) ;
		return new Viscoelasticity( PURE_VISCOSITY, getStiffnessMatrix(dim, pt)*getData("characteristic_time",1.), blocks) ;
	}

	if(type == std::string("KELVIN_VOIGT"))
	{
		int blocks = getData("additional_viscoelastic_variables",0) ;
		return new Viscoelasticity( KELVIN_VOIGT, getStiffnessMatrix(dim, pt), getStiffnessMatrix(dim, pt)*getData("characteristic_time",1.), blocks) ;
	}

	if(type == std::string("MAXWELL"))
	{
		int blocks = getData("additional_viscoelastic_variables",0) ;
		return new Viscoelasticity( MAXWELL, getStiffnessMatrix(dim, pt), getStiffnessMatrix(dim, pt)*getData("characteristic_time",1.), blocks) ;
	}

	if(type == std::string("BURGER"))
	{
		int blocks = getData("additional_viscoelastic_variables",0) ;
		return new Viscoelasticity( BURGER, getChildFromFullLabel("kelvin_voigt")->getStiffnessMatrix(dim, pt), getChildFromFullLabel("kelvin_voigt")->getStiffnessMatrix(dim, pt)*getData("kelvin_voigt.characteristic_time",1.), getChildFromFullLabel("maxwell")->getStiffnessMatrix(dim, pt), getChildFromFullLabel("maxwell")->getStiffnessMatrix(dim, pt)*getData("maxwell.characteristic_time",1.),blocks) ;
	}

	if(type == std::string("GENERALIZED_KELVIN_VOIGT"))
	{
		int blocks = getData("additional_viscoelastic_variables",0) ;
		std::vector<std::pair<Matrix, Matrix> > branches ;
		std::vector<ConfigTreeItem *> config = getAllChildren("branch") ;
		for(size_t i = 0 ; i < config.size() ; i++)
			branches.push_back(std::make_pair( config[i]->getStiffnessMatrix(dim, pt), config[i]->getStiffnessMatrix(dim, pt)*config[i]->getData("characteristic_time",1.) ) ) ;
		return new Viscoelasticity( GENERALIZED_KELVIN_VOIGT, getChild("first_branch")->getStiffnessMatrix(dim, pt), branches, blocks) ;
	}

	if(type == std::string("GENERALIZED_MAXWELL"))
	{
		int blocks = getData("additional_viscoelastic_variables",0) ;
		std::vector<std::pair<Matrix, Matrix> > branches ;
		std::vector<ConfigTreeItem *> config = getAllChildren("branch") ;
		for(size_t i = 0 ; i < config.size() ; i++)
			branches.push_back(std::make_pair( config[i]->getStiffnessMatrix(dim, pt), config[i]->getStiffnessMatrix(dim, pt)*config[i]->getData("characteristic_time",1.) ) ) ;
		return new Viscoelasticity( GENERALIZED_MAXWELL, getChild("first_branch")->getStiffnessMatrix(dim, pt), branches, blocks) ;
	}

	return new VoidForm() ;
}

Vector ConfigTreeItem::readVectorFromFile() const 
{
	std::fstream in ;
	in.open( getStringData().c_str(), std::ios::in ) ;
	if(in.fail())
	{
		std::cout << "unable to read vector from file " << getStringData() << std::endl ;
		return Vector() ;
	}
	std::vector<double> tmp ;
	while(!in.eof())
	{
		double buff ;
		in >> buff ;
		tmp.push_back(buff) ;
	}
	Vector ret(tmp.size()-1) ;
	for(size_t i = 0 ; i < ret.size() ; i++)
		ret[i] = tmp[i] ;
	return ret ;
}

Order ConfigTreeItem::translateOrder(std::string order)
{
	if(order == "CONSTANT")
		return CONSTANT ;
	if(order == "LINEAR")
		return LINEAR ;
	if(order == "QUADRATIC")
		return QUADRATIC ;
	if(order == "CUBIC")
		return CUBIC ;
	if(order == "QUADRIC")
		return QUADRIC ;
	if(order == "QUINTIC")
		return QUINTIC ;
	if(order == "CONSTANT_TIME_LINEAR")
		return CONSTANT_TIME_LINEAR ;
	if(order == "CONSTANT_TIME_QUADRATIC")
		return CONSTANT_TIME_QUADRATIC ;
	if(order == "LINEAR_TIME_LINEAR")
		return LINEAR_TIME_LINEAR ;
	if(order == "LINEAR_TIME_QUADRATIC")
		return LINEAR_TIME_QUADRATIC ;
	if(order == "QUADRATIC_TIME_LINEAR")
		return QUADRATIC_TIME_LINEAR ;
	if(order == "QUADRATIC_TIME_QUADRATIC")
		return QUADRATIC_TIME_QUADRATIC ;
	if(order == "CUBIC_TIME_LINEAR")
		return CUBIC_TIME_LINEAR ;
	if(order == "CUBIC_TIME_QUADRATIC")
		return CUBIC_TIME_QUADRATIC ;
	if(order == "QUADRIC_TIME_LINEAR")
		return QUADRIC_TIME_LINEAR ;
	if(order == "QUADRIC_TIME_QUADRATIC")
		return QUADRIC_TIME_QUADRATIC ;
	if(order == "QUINTIC_TIME_LINEAR")
		return QUINTIC_TIME_LINEAR ;
	if(order == "QUINTIC_TIME_QUADRATIC")
		return QUINTIC_TIME_QUADRATIC ;
	if(order == "QUADTREE_REFINED")
		return QUADTREE_REFINED ;
	if(order == "REGULAR_GRID")
		return REGULAR_GRID ;
	return LINEAR ;
}

SamplingRestrictionType ConfigTreeItem::translateSamplingRestrictionType( std::string restriction ) 
{
	if(restriction == "SAMPLE_RESTRICT_4")
		return SAMPLE_RESTRICT_4 ;
	if(restriction == "SAMPLE_RESTRICT_8")
		return SAMPLE_RESTRICT_8 ;
	if(restriction == "SAMPLE_RESTRICT_16")
		return SAMPLE_RESTRICT_16 ;
	return SAMPLE_NO_RESTRICTION ;
}

planeType ConfigTreeItem::translatePlaneType( std::string type ) 
{
	if(type == "PLANE_STRAIN")
		return PLANE_STRAIN ;
	if(type == "PLANE_STRESS")
		return PLANE_STRESS ;
	if(type == "PLANE_STRESS_FREE_G")
		return PLANE_STRESS_FREE_G ;
	return PLANE_STRESS ;
}

PSDSpecificationType ConfigTreeItem::translatePSDSpecificationType( std::string specification ) 
{
	if(specification == "CUMULATIVE_PERCENT")
		return CUMULATIVE_PERCENT ;
	if(specification == "CUMULATIVE_FRACTION")
		return CUMULATIVE_FRACTION ;
	if(specification == "CUMULATIVE_ABSOLUTE")
		return CUMULATIVE_ABSOLUTE ;
	if(specification == "CUMULATIVE_PERCENT_REVERSE")
		return CUMULATIVE_PERCENT_REVERSE ;
	if(specification == "CUMULATIVE_FRACTION_REVERSE")
		return CUMULATIVE_FRACTION_REVERSE ;
	if(specification == "CUMULATIVE_ABSOLUTE_REVERSE")
		return CUMULATIVE_ABSOLUTE ;
	return CUMULATIVE_PERCENT ;
}

ParticleSizeDistribution * ConfigTreeItem::getParticleSizeDistribution() const 
{
	std::string type = getStringData("type","CONSTANT") ;

	if(type == std::string("CONSTANT"))
	{
		return new ConstantSizeDistribution() ;
	}

	if(type == std::string("BOLOME_A"))
	{
		return new PSDBolomeA() ;
	}

	if(type == std::string("BOLOME_B"))
	{
		return new PSDBolomeB() ;
	}

	if(type == std::string("BOLOME_C"))
	{
		return new PSDBolomeC() ;
	}

	if(type == std::string("BOLOME_D"))
	{
		return new PSDBolomeD() ;
	}
	
	if(type == std::string("FROM_CUMULATIVE_FILE"))
	{
		std::string filename = getStringData("file_name","file_not_found") ;
		if(filename == std::string("file_not_found"))
		{
			std::cout << "no psd file specified, falling back to constant particle size distribution" << std::endl ;
			return new ConstantSizeDistribution() ;
		}
		std::string specification = getStringData("psd_specification_type","CUMULATIVE_PERCENT") ;
		double factor = getData("factor",1) ;
		double cutOffUp = getData("cutoff.up", -1.) ;
		double cutOffDown = getData("cutoff.down", -1.) ;
		return new GranuloFromCumulativePSD(filename, ConfigTreeItem::translatePSDSpecificationType( specification ), factor, cutOffUp, cutOffDown) ;
	}

	return new ConstantSizeDistribution() ;
}

TypeInclusion ConfigTreeItem::translateInclusionType(std::string type) 
{
	if(type == std::string("CIRCLE"))
		return CIRCLE_INCLUSION ;
	if(type == std::string("SPHERE"))
		return SPHERE_INCLUSION ;
	if(type == std::string("ELLIPSE"))
		return ELLIPSE_INCLUSION ;
	return CIRCLE_INCLUSION ;
}

GeometryType ConfigTreeItem::translateGeometryType(std::string type) 
{
	if(type == std::string("CIRCLE"))
		return CIRCLE ;
	if(type == std::string("LAYERED_CIRCLE"))
		return LAYERED_CIRCLE ;
	if(type == std::string("TRIANGLE"))
		return TRIANGLE ;
	if(type == std::string("RECTANGLE"))
		return RECTANGLE ;
	if(type == std::string("PARALLELOGRAMME"))
		return PARALLELOGRAMME ;
	if(type == std::string("CONVEX_POLYGON"))
		return CONVEX_POLYGON ;
	if(type == std::string("SEGMENTED_LINE"))
		return SEGMENTED_LINE ;
	if(type == std::string("ORIENTABLE_CIRCLE"))
		return ORIENTABLE_CIRCLE ;
	if(type == std::string("CLOSED_NURB"))
		return CLOSED_NURB ;
	if(type == std::string("TETRAHEDRON"))
		return TETRAHEDRON ;
	if(type == std::string("HEXAHEDRON"))
		return HEXAHEDRON ;
	if(type == std::string("SPHERE"))
		return SPHERE ;
	if(type == std::string("LAYERED_SPHERE"))
		return LAYERED_SPHERE ;
	if(type == std::string("REGULAR_OCTAHEDRON"))
		return REGULAR_OCTAHEDRON ;
	if(type == std::string("ELLIPSE"))
		return ELLIPSE ;
	if(type == std::string("LEVEL_SET"))
		return LEVEL_SET ;
	if(type == std::string("TIME_DEPENDENT_CIRCLE"))
		return TIME_DEPENDENT_CIRCLE ;
	return CIRCLE ;
}


std::vector<Feature *> ConfigTreeItem::getInclusions(FeatureTree * F, std::vector<Feature *> base, std::vector<Geometry *> brothers) const 
{
	std::vector<Feature *> ret ;
	ConfigTreeItem * psdConfig = getChild("particle_size_distribution") ;
	if(!psdConfig)
		return ret ;
	std::string type = psdConfig->getStringData("type","CONSTANT") ;

	Form * behaviour = nullptr ;
	if(hasChild("behaviour"))
		behaviour = getChild("behaviour")->getBehaviour( F->is2D() ? SPACE_TWO_DIMENSIONAL : SPACE_THREE_DIMENSIONAL , F->getOrder() >= CONSTANT_TIME_LINEAR ) ;

	if(type == "FROM_INCLUSION_FILE")
	{
		std::vector<std::string> columns ;
		std::vector<ConfigTreeItem *> tree = psdConfig->getAllChildren("column") ;
		for(size_t i = 0 ; i < tree.size() ; i++)
			columns.push_back(tree[i]->getStringData()) ;
		std::string filename = psdConfig->getStringData("file_name","file_not_found") ;
		if(filename == std::string("file_not_found"))
		{
			std::cout << "no inclusion file specified" << std::endl ;
			return ret ;
		}
		GranuloFromFile granulo( filename, columns ) ;
		std::string inclusion = getStringData("geometry.type","CIRCLE") ;
		int n = psdConfig->getData("number",1) ;
		ret = granulo.getFeatures( ConfigTreeItem::translateInclusionType( inclusion ), n ) ;
		if(behaviour)
		{
			for(size_t i = 0 ; i < ret.size() ; i++)
				ret[i]->setBehaviour( behaviour ) ;
		}
		for(size_t i = 0 ; i < ret.size() ; i++)
		{
			bool placed = false ;
			for(size_t j = 0 ; j < base.size() ; j++)
			{
				if(base[j]->in(ret[i]->getCenter()))
				{
					F->addFeature( base[j], ret[i]) ;
					placed = true ;
					break ; 
				}
			}
			if(!placed)
				F->addFeature( F->getFeature(0), ret[i]) ;
		}
	}
	else if(type == "FROM_PARENT_DISTRIBUTION")
	{
		Function f("x") ;
		if(psdConfig->hasChild("layer_thickness_function"))
			f = psdConfig->getChild("layer_thickness_function")->getFunction() ;
		else
			f -= psdConfig->getData("layer_thickness",0.0001) ;
		
		VirtualMachine vm ;
		for(size_t i = 0 ; i < base.size() ; i++)
		{
			double newr = vm.eval(f, base[i]->getRadius()) ;
			if(newr > POINT_TOLERANCE_2D && newr < base[i]->getRadius() - POINT_TOLERANCE_2D)
			{
				Inclusion * inc = new Inclusion( newr, base[i]->getCenter().getX(), base[i]->getCenter().getY()) ;
				inc->setBehaviour(behaviour) ;
				F->addFeature( base[i], inc) ;
				ret.push_back(inc) ;
			}
		}

	}
	else
	{
		ParticleSizeDistribution * psd = psdConfig->getParticleSizeDistribution() ;
		int n = psdConfig->getData("number",1) ;
		double rmax = psdConfig->getData("rmax",0.1) ;
		double fraction = psdConfig->getData("fraction", 0.8) ;
		std::string geometry = getStringData("geometry.type","CIRCLE") ;
		double aspectRatio = getData("geometry.aspect_ratio",1.) ;
		double orientation = getData("geometry.orientation",M_PI) ;
		Sample * placement = nullptr ;
		int tries = getData("placement.tries", 1e6) ;
		double spacing = getData("placement.spacing",1e-5) ;
		if(hasChildFromFullLabel("placement.box"))
		{
			placement = getChildFromFullLabel("placement.box")->getSample() ;
		}
		if(base.size() == 0)
		{
			ret = PSDGenerator::get2DConcrete( F, behaviour, n, rmax, spacing, psd, ConfigTreeItem::translateGeometryType(geometry), aspectRatio, orientation, tries, fraction, dynamic_cast<Rectangle *>(placement), brothers) ;
			std::cout << ret.size() << std::endl ;
		}
		else
		{
			ret = PSDGenerator::get2DEmbeddedInclusions( F, behaviour, base, n, rmax, spacing, psd, ConfigTreeItem::translateGeometryType(geometry), aspectRatio, orientation, tries, fraction, dynamic_cast<Rectangle *>(placement), brothers) ;
		}
	}

	if(hasChild("sampling_factor"))
	{
		for(size_t i = 0 ; i < ret.size() ; i++)
			F->setSamplingFactor( ret[i], getData("sampling_factor",1.) ) ;
	}

	if(hasChild("intersection_sampling_factor"))
	{
		for(size_t i = 0 ; i < ret.size() ; i++)
		{
			if(F->getFeature(0)->intersects(ret[i]))
				F->setSamplingFactor( ret[i], getData("intersection_sampling_factor",1.) ) ;
		}
	}

	if(hasChild("inclusions"))
	{
		std::vector<Feature *> newbase = ret ;
		std::vector<Geometry *> newbrothers ;
		std::vector<ConfigTreeItem *> newInclusions = getAllChildren("inclusions") ;
		for(size_t i = 0 ; i < newInclusions.size() ; i++)
		{
			std::vector<Feature *> tmp = newInclusions[i]->getInclusions( F, newbase, newbrothers ) ;
			for(size_t j = 0 ; j < tmp.size() ; j++)
				newbrothers.push_back( dynamic_cast<Geometry *>(tmp[j]) ) ;
		}
	}

	return ret ;
}

LagrangeMultiplierType ConfigTreeItem::translateLagrangeMultiplierType(std::string type) 
{
	if(type == "GENERAL")
		return GENERAL ;
	if(type == "FIX_ALONG_ALL")
		return FIX_ALONG_ALL ;
	if(type == "FIX_ALONG_XI")
		return FIX_ALONG_XI ;
	if(type == "FIX_ALONG_ETA")
		return FIX_ALONG_ETA ;
	if(type == "FIX_ALONG_ZETA")
		return FIX_ALONG_ZETA ;
	if(type == "FIX_ALONG_XI_ETA")
		return FIX_ALONG_XI_ETA ;
	if(type == "FIX_ALONG_XI_ZETA")
		return FIX_ALONG_XI_ZETA ;
	if(type == "FIX_ALONG_ETA_ZETA")
		return FIX_ALONG_ETA_ZETA ;
	if(type == "SET_ALONG_XI")
		return SET_ALONG_XI ;
	if(type == "SET_ALONG_ETA")
		return SET_ALONG_ETA ;
	if(type == "SET_ALONG_ZETA")
		return SET_ALONG_ZETA ;
	if(type == "SET_ALONG_XI_ETA")
		return SET_ALONG_XI_ETA ;
	if(type == "SET_ALONG_XI_ZETA")
		return SET_ALONG_XI_ZETA ;
	if(type == "SET_ALONG_ETA_ZETA")
		return SET_ALONG_ETA_ZETA ;
	if(type == "FIX_ALONG_INDEXED_AXIS")
		return FIX_ALONG_INDEXED_AXIS ;
	if(type == "SET_ALONG_INDEXED_AXIS")
		return SET_ALONG_INDEXED_AXIS ;
	if(type == "SET_FORCE_XI")
		return SET_FORCE_XI ;
	if(type == "SET_FORCE_ETA")
		return SET_FORCE_ETA ;
	if(type == "SET_FORCE_ZETA")
		return SET_FORCE_ZETA ;
	if(type == "SET_FORCE_INDEXED_AXIS")
		return SET_FORCE_INDEXED_AXIS ;
	if(type == "SET_FLUX_XI")
		return SET_FLUX_XI ;
	if(type == "SET_FLUX_ETA")
		return SET_FLUX_ETA ;
	if(type == "SET_FLUX_ZETA")
		return SET_FLUX_ZETA ;
	if(type == "SET_STRESS_XI")
		return SET_STRESS_XI ;
	if(type == "SET_STRESS_ETA")
		return SET_STRESS_ETA ;
	if(type == "SET_STRESS_ZETA")
		return SET_STRESS_ZETA ;
	if(type == "SET_VOLUMIC_STRESS_XI")
		return SET_VOLUMIC_STRESS_XI ;
	if(type == "SET_VOLUMIC_STRESS_ETA")
		return SET_VOLUMIC_STRESS_ETA ;
	if(type == "SET_VOLUMIC_STRESS_ZETA")
		return SET_VOLUMIC_STRESS_ZETA ;
	if(type == "SET_NORMAL_STRESS")
		return SET_NORMAL_STRESS ;
	if(type == "SET_TANGENT_STRESS")
		return SET_TANGENT_STRESS ;
	if(type == "VERTICAL_PLANE_SECTIONS")
		return VERTICAL_PLANE_SECTIONS ;
	if(type == "HORIZONTAL_PLANE_SECTIONS")
		return HORIZONTAL_PLANE_SECTIONS ;
	if(type == "SET_GLOBAL_FORCE_VECTOR")
		return SET_GLOBAL_FORCE_VECTOR ;
	if(type == "nullptr_CONDITION")
		return nullptr_CONDITION ;
	return GENERAL ;
}

BoundingBoxPosition ConfigTreeItem::translateBoundingBoxPosition( std::string position ) 
{
	if(position == "TOP")
		return TOP ;
	if(position == "TOP_LEFT")
		return TOP_LEFT ;
	if(position == "TOP_RIGHT")
		return TOP_RIGHT ;
	if(position == "TOP_BACK")
		return TOP_BACK ;
	if(position == "TOP_LEFT_FRONT")
		return TOP_LEFT_FRONT ;
	if(position == "TOP_LEFT_BACK")
		return TOP_LEFT_BACK ;
	if(position == "TOP_RIGHT_FRONT")
		return TOP_RIGHT_FRONT ;
	if(position == "TOP_RIGHT_BACK")
		return TOP_RIGHT_BACK ;
	if(position == "BOTTOM")
		return BOTTOM ;
	if(position == "BOTTOM_LEFT")
		return BOTTOM_LEFT ;
	if(position == "BOTTOM_RIGHT")
		return BOTTOM_RIGHT ;
	if(position == "BOTTOM_BACK")
		return BOTTOM_BACK ;
	if(position == "BOTTOM_LEFT_FRONT")
		return BOTTOM_LEFT_FRONT ;
	if(position == "BOTTOM_LEFT_BACK")
		return BOTTOM_LEFT_BACK ;
	if(position == "BOTTOM_RIGHT_FRONT")
		return BOTTOM_RIGHT_FRONT ;
	if(position == "BOTTOM_RIGHT_BACK")
		return BOTTOM_RIGHT_BACK ;
	if(position == "LEFT")
		return LEFT ;
	if(position == "RIGHT")
		return RIGHT ;
	if(position == "FRONT")
		return FRONT ;
	if(position == "FRONT_LEFT")
		return FRONT_LEFT ;
	if(position == "FRONT_RIGHT")
		return FRONT_RIGHT ;
	if(position == "FRONT_TOP")
		return FRONT_TOP ;
	if(position == "FRONT_BOTTOM")
		return FRONT_BOTTOM ;
	if(position == "BACK")
		return BACK ;
	if(position == "BACK_LEFT")
		return BACK_LEFT ;
	if(position == "BACK_RIGHT")
		return BACK_RIGHT ;
	if(position == "BEFORE")
		return BEFORE ;
	if(position == "NOW")
		return NOW ;
	if(position == "AFTER")
		return AFTER ;
	if(position == "TOP_BEFORE")
		return TOP_BEFORE ;
	if(position == "TOP_LEFT_BEFORE")
		return TOP_LEFT_BEFORE ;
	if(position == "TOP_RIGHT_BEFORE")
		return TOP_RIGHT_BEFORE ;
	if(position == "TOP_BACK_BEFORE")
		return TOP_BACK_BEFORE ;
	if(position == "TOP_LEFT_FRONT_BEFORE")
		return TOP_LEFT_FRONT_BEFORE ;
	if(position == "TOP_LEFT_BACK_BEFORE")
		return TOP_LEFT_BACK_BEFORE ;
	if(position == "TOP_RIGHT_FRONT_BEFORE")
		return TOP_RIGHT_FRONT_BEFORE ;
	if(position == "TOP_RIGHT_BACK_BEFORE")
		return TOP_RIGHT_BACK_BEFORE ;
	if(position == "BOTTOM_BEFORE")
		return BOTTOM_BEFORE ;
	if(position == "BOTTOM_LEFT_BEFORE")
		return BOTTOM_LEFT_BEFORE ;
	if(position == "BOTTOM_RIGHT_BEFORE")
		return BOTTOM_RIGHT_BEFORE ;
	if(position == "BOTTOM_BACK_BEFORE")
		return BOTTOM_BACK_BEFORE ;
	if(position == "BOTTOM_LEFT_FRONT_BEFORE")
		return BOTTOM_LEFT_FRONT_BEFORE ;
	if(position == "BOTTOM_LEFT_BACK_BEFORE")
		return BOTTOM_LEFT_BACK_BEFORE ;
	if(position == "BOTTOM_RIGHT_FRONT_BEFORE")
		return BOTTOM_RIGHT_FRONT_BEFORE ;
	if(position == "BOTTOM_RIGHT_BACK_BEFORE")
		return BOTTOM_RIGHT_BACK_BEFORE ;
	if(position == "LEFT_BEFORE")
		return LEFT_BEFORE ;
	if(position == "RIGHT_BEFORE")
		return RIGHT_BEFORE ;
	if(position == "FRONT_BEFORE")
		return FRONT_BEFORE ;
	if(position == "FRONT_LEFT_BEFORE")
		return FRONT_LEFT_BEFORE ;
	if(position == "FRONT_RIGHT_BEFORE")
		return FRONT_RIGHT_BEFORE ;
	if(position == "FRONT_TOP_BEFORE")
		return FRONT_TOP_BEFORE ;
	if(position == "FRONT_BOTTOM_BEFORE")
		return FRONT_BOTTOM_BEFORE ;
	if(position == "BACK_BEFORE")
		return BACK_BEFORE ;
	if(position == "BACK_LEFT_BEFORE")
		return BACK_LEFT_BEFORE ;
	if(position == "BACK_RIGHT_BEFORE")
		return BACK_RIGHT_BEFORE ;
	if(position == "TOP_NOW")
		return TOP_NOW ;
	if(position == "TOP_LEFT_NOW")
		return TOP_LEFT_NOW ;
	if(position == "TOP_RIGHT_NOW")
		return TOP_RIGHT_NOW ;
	if(position == "TOP_LEFT_FRONT_NOW")
		return TOP_LEFT_FRONT_NOW ;
	if(position == "TOP_LEFT_BACK_NOW")
		return TOP_LEFT_BACK_NOW ;
	if(position == "TOP_RIGHT_FRONT_NOW")
		return TOP_RIGHT_FRONT_NOW ;
	if(position == "TOP_RIGHT_BACK_NOW")
		return TOP_RIGHT_BACK_NOW ;
	if(position == "BOTTOM_NOW")
		return BOTTOM_NOW ;
	if(position == "BOTTOM_LEFT_NOW")
		return BOTTOM_LEFT_NOW ;
	if(position == "BOTTOM_RIGHT_NOW")
		return BOTTOM_RIGHT_NOW ;
	if(position == "BOTTOM_LEFT_FRONT_NOW")
		return BOTTOM_LEFT_FRONT_NOW ;
	if(position == "BOTTOM_LEFT_BACK_NOW")
		return BOTTOM_LEFT_BACK_NOW ;
	if(position == "BOTTOM_RIGHT_FRONT_NOW")
		return BOTTOM_RIGHT_FRONT_NOW ;
	if(position == "BOTTOM_RIGHT_BACK_NOW")
		return BOTTOM_RIGHT_BACK_NOW ;
	if(position == "LEFT_NOW")
		return LEFT_NOW ;
	if(position == "RIGHT_NOW")
		return RIGHT_NOW ;
	if(position == "FRONT_NOW")
		return FRONT_NOW ;
	if(position == "FRONT_LEFT_NOW")
		return FRONT_LEFT_NOW ;
	if(position == "FRONT_RIGHT_NOW")
		return FRONT_RIGHT_NOW ;
	if(position == "FRONT_TOP_NOW")
		return FRONT_TOP_NOW ;
	if(position == "FRONT_BOTTOM_NOW")
		return FRONT_BOTTOM_NOW ;
	if(position == "BACK_NOW")
		return BACK_NOW ;
	if(position == "BACK_LEFT_NOW")
		return BACK_LEFT_NOW ;
	if(position == "BACK_RIGHT_NOW")
		return BACK_RIGHT_NOW ;
	if(position == "TOP_AFTER")
		return TOP_AFTER ;
	if(position == "TOP_LEFT_AFTER")
		return TOP_LEFT_AFTER ;
	if(position == "TOP_RIGHT_AFTER")
		return TOP_RIGHT_AFTER ;
	if(position == "TOP_BACK_AFTER")
		return TOP_BACK_AFTER ;
	if(position == "TOP_LEFT_FRONT_AFTER")
		return TOP_LEFT_FRONT_AFTER ;
	if(position == "TOP_LEFT_BACK_AFTER")
		return TOP_LEFT_BACK_AFTER ;
	if(position == "TOP_RIGHT_FRONT_AFTER")
		return TOP_RIGHT_FRONT_AFTER ;
	if(position == "TOP_RIGHT_BACK_AFTER")
		return TOP_RIGHT_BACK_AFTER ;
	if(position == "BOTTOM_AFTER")
		return BOTTOM_AFTER ;
	if(position == "BOTTOM_LEFT_AFTER")
		return BOTTOM_LEFT_AFTER ;
	if(position == "BOTTOM_RIGHT_AFTER")
		return BOTTOM_RIGHT_AFTER ;
	if(position == "BOTTOM_BACK_AFTER")
		return BOTTOM_BACK_AFTER ;
	if(position == "BOTTOM_LEFT_FRONT_AFTER")
		return BOTTOM_LEFT_FRONT_AFTER ;
	if(position == "BOTTOM_LEFT_BACK_AFTER")
		return BOTTOM_LEFT_BACK_AFTER ;
	if(position == "BOTTOM_RIGHT_FRONT_AFTER")
		return BOTTOM_RIGHT_FRONT_AFTER ;
	if(position == "BOTTOM_RIGHT_BACK_AFTER")
		return BOTTOM_RIGHT_BACK_AFTER ;
	if(position == "LEFT_AFTER")
		return LEFT_AFTER ;
	if(position == "RIGHT_AFTER")
		return RIGHT_AFTER ;
	if(position == "FRONT_AFTER")
		return FRONT_AFTER ;
	if(position == "FRONT_LEFT_AFTER")
		return FRONT_LEFT_AFTER ;
	if(position == "FRONT_RIGHT_AFTER")
		return FRONT_RIGHT_AFTER ;
	if(position == "FRONT_TOP_AFTER")
		return FRONT_TOP_AFTER ;
	if(position == "FRONT_BOTTOM_AFTER")
		return FRONT_BOTTOM_AFTER ;
	if(position == "BACK_AFTER")
		return BACK_AFTER ;
	if(position == "BACK_LEFT_AFTER")
		return BACK_LEFT_AFTER ;
	if(position == "BACK_RIGHT_AFTER")
		return BACK_RIGHT_AFTER ;	
	return NOW ;
}


BoundaryCondition * ConfigTreeItem::getBoundaryCondition() const 
{
	if( !hasChild("rate") )
		return new BoundingBoxDefinedBoundaryCondition( ConfigTreeItem::translateLagrangeMultiplierType( getStringData( "condition", "GENERAL" ) ),
			ConfigTreeItem::translateBoundingBoxPosition( getStringData( "position", "NOW" ) ),
			getData( "value", 0 ), getData( "axis", 0 ) ) ;
	Function t("t") ;
	t *= getData( "rate", 0 ) ;
	return new BoundingBoxDefinedBoundaryCondition( ConfigTreeItem::translateLagrangeMultiplierType( getStringData( "condition", "GENERAL" ) ),
		ConfigTreeItem::translateBoundingBoxPosition( getStringData( "position", "NOW" ) ),
		t, getData( "axis", 0 ) ) ;
}

bool ConfigTreeItem::isAtTimeStep(int i, int nsteps) const
{
	if(getStringData("at","NONE") == "ALL")
		return true ;
	if(getStringData("at","NONE") == "LAST")
		return i == nsteps-1 ;
	if(getStringData("at","NONE") == "REGULAR")
		return i%((int) getData("every", 2)) == 0 ;
	return false ;
}

TWFieldType ConfigTreeItem::translateTriangleWriterFieldType( std::string field, bool & ok)
{
	ok = true ;
	if(field == "CRITERION")
		return TWFT_CRITERION ;
	if(field == "STIFFNESS")
		return TWFT_STIFFNESS ;
	if(field == "VISCOSITY")
		return TWFT_VISCOSITY ;
	ok = false ;
	return TWFT_COORDINATE ;
}


FieldType ConfigTreeItem::translateFieldType( std::string field, bool & ok)
{
	ok = true ;
	if(field == "DISPLACEMENT_FIELD")
		return DISPLACEMENT_FIELD ;
	if(field == "ENRICHED_DISPLACEMENT_FIELD")
		return ENRICHED_DISPLACEMENT_FIELD ;
	if(field == "SPEED_FIELD")
		return SPEED_FIELD ;
	if(field == "FLUX_FIELD")
		return FLUX_FIELD ;
	if(field == "GRADIENT_FIELD")
		return GRADIENT_FIELD ;
	if(field == "STRAIN_RATE_FIELD")
		return STRAIN_RATE_FIELD ;
	if(field == "STRAIN_FIELD")
		return STRAIN_FIELD ;
	if(field == "EFFECTIVE_STRESS_FIELD")
		return EFFECTIVE_STRESS_FIELD ;
	if(field == "REAL_STRESS_FIELD")
		return REAL_STRESS_FIELD ;
	if(field == "PRINCIPAL_STRAIN_FIELD")
		return PRINCIPAL_STRAIN_FIELD ;
	if(field == "PRINCIPAL_EFFECTIVE_STRESS_FIELD")
		return PRINCIPAL_EFFECTIVE_STRESS_FIELD ;
	if(field == "PRINCIPAL_REAL_STRESS_FIELD")
		return PRINCIPAL_REAL_STRESS_FIELD ;
	if(field == "NON_ENRICHED_STRAIN_RATE_FIELD")
		return NON_ENRICHED_STRAIN_RATE_FIELD ;
	if(field == "NON_ENRICHED_STRAIN_FIELD")
		return NON_ENRICHED_STRAIN_FIELD ;
	if(field == "NON_ENRICHED_EFFECTIVE_STRESS_FIELD")
		return NON_ENRICHED_EFFECTIVE_STRESS_FIELD ;
	if(field == "NON_ENRICHED_REAL_STRESS_FIELD")
		return NON_ENRICHED_REAL_STRESS_FIELD ;
	if(field == "VON_MISES_STRAIN_FIELD")
		return VON_MISES_STRAIN_FIELD ;
	if(field == "VON_MISES_EFFECTIVE_STRESS_FIELD")
		return VON_MISES_EFFECTIVE_STRESS_FIELD ;
	if(field == "VON_MISES_REAL_STRESS_FIELD")
		return VON_MISES_REAL_STRESS_FIELD ;
	if(field == "PRINCIPAL_ANGLE_FIELD")
		return PRINCIPAL_ANGLE_FIELD ;
	if(field == "INTERNAL_VARIABLE_FIELD")
		return INTERNAL_VARIABLE_FIELD ;
	if(field == "GENERALIZED_VISCOELASTIC_DISPLACEMENT_FIELD")
		return GENERALIZED_VISCOELASTIC_DISPLACEMENT_FIELD ;
	if(field == "GENERALIZED_VISCOELASTIC_ENRICHED_DISPLACEMENT_FIELD")
		return GENERALIZED_VISCOELASTIC_ENRICHED_DISPLACEMENT_FIELD ;
	if(field == "GENERALIZED_VISCOELASTIC_SPEED_FIELD")
		return GENERALIZED_VISCOELASTIC_SPEED_FIELD ;
	if(field == "GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD")
		return GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD ;
	if(field == "GENERALIZED_VISCOELASTIC_STRAIN_FIELD")
		return GENERALIZED_VISCOELASTIC_STRAIN_FIELD ;
	if(field == "GENERALIZED_VISCOELASTIC_EFFECTIVE_STRESS_FIELD")
		return GENERALIZED_VISCOELASTIC_EFFECTIVE_STRESS_FIELD ;
	if(field == "GENERALIZED_VISCOELASTIC_REAL_STRESS_FIELD")
		return GENERALIZED_VISCOELASTIC_REAL_STRESS_FIELD ;
	if(field == "GENERALIZED_VISCOELASTIC_PRINCIPAL_STRAIN_FIELD")
		return GENERALIZED_VISCOELASTIC_PRINCIPAL_STRAIN_FIELD ;
	if(field == "GENERALIZED_VISCOELASTIC_PRINCIPAL_EFFECTIVE_STRESS_FIELD")
		return GENERALIZED_VISCOELASTIC_PRINCIPAL_EFFECTIVE_STRESS_FIELD ;
	if(field == "GENERALIZED_VISCOELASTIC_PRINCIPAL_REAL_STRESS_FIELD")
		return GENERALIZED_VISCOELASTIC_PRINCIPAL_REAL_STRESS_FIELD ;
	if(field == "GENERALIZED_VISCOELASTIC_NON_ENRICHED_STRAIN_RATE_FIELD")
		return GENERALIZED_VISCOELASTIC_NON_ENRICHED_STRAIN_RATE_FIELD ;
	if(field == "GENERALIZED_VISCOELASTIC_NON_ENRICHED_STRAIN_FIELD")
		return GENERALIZED_VISCOELASTIC_NON_ENRICHED_STRAIN_FIELD ;
	if(field == "GENERALIZED_VISCOELASTIC_NON_ENRICHED_EFFECTIVE_STRESS_FIELD")
		return GENERALIZED_VISCOELASTIC_NON_ENRICHED_EFFECTIVE_STRESS_FIELD ;
	if(field == "GENERALIZED_VISCOELASTIC_NON_ENRICHED_REAL_STRESS_FIELD")
		return NON_ENRICHED_REAL_STRESS_FIELD ;
	ok = false ;
	return DISPLACEMENT_FIELD ;
}


void ConfigTreeItem::writeOutput(FeatureTree * F, int i, int nsteps) 
{
	if(getStringData("file_name","file_not_found") == "file_not_found")
		return ;

	if(i == 1)
	{
		std::fstream out ;
		out.open(getStringData("file_name","output").c_str(), std::ios::out) ;
		out.close() ;
	}
	if(getChild("time_step")->isAtTimeStep(i, nsteps))
	{
		std::fstream out ;
		out.open(getStringData("file_name","output").c_str(), std::ios::out | std::ios::app) ;
		out << F->getCurrentTime() << "\t" ;
		std::vector<ConfigTreeItem *> fields = getAllChildren("field") ;
		std::string instant = getStringData("instant","NOW") ;
		for(size_t i = 0 ; i < fields.size() ; i++)
		{
			bool isFieldType = true ;
			Vector f = F->getAverageField( ConfigTreeItem::translateFieldType( fields[i]->getStringData(), isFieldType ), -1, (instant == "AFTER") - (instant == "BEFORE") ) ;
			for(size_t j = 0 ; j < f.size() ; j++)
				out << f[j] << "\t" ;
		}
		out << std::endl ;
		out.close() ;
	}
}

void ConfigTreeItem::exportTriangles(FeatureTree * F, int i, int nsteps) 
{
	if(getStringData("file_name","file_not_found") == "file_not_found")
		return ;

	if(getChild("time_step")->isAtTimeStep(i, nsteps))
	{
		std::string instant = getStringData("instant","NOW") ;
		TriangleWriter trg(getStringData("file_name","output_").append(itoa(i)), F, (instant == "AFTER") - (instant == "BEFORE") ) ;
		std::vector<ConfigTreeItem *> fields = getAllChildren("field") ;
		for(size_t i = 0 ; i < fields.size() ; i++)
		{
			bool isFieldType = true ;
			FieldType f = ConfigTreeItem::translateFieldType( fields[i]->getStringData(), isFieldType ) ;
			bool isTWFieldType = false ;
			TWFieldType g = ConfigTreeItem::translateTriangleWriterFieldType( fields[i]->getStringData(), isTWFieldType ) ;
			if(isFieldType)
				trg.getField( f ) ;
			else if(isTWFieldType)
				trg.getField( g ) ;
			else
				trg.getField(  fields[i]->getStringData() ) ;
		}
		trg.write() ;
	}
}

void ConfigTreeItem::exportSvgTriangles(MultiTriangleWriter * trg, FeatureTree * F, int i, int nsteps) 
{
	if(!trg)
		return ;
	
	if(getChild("time_step")->isAtTimeStep(i, nsteps))
	{
		std::string instant = getStringData("instant","NOW") ;
		trg->reset(F, (instant == "AFTER") - (instant == "BEFORE")) ;
		std::vector<ConfigTreeItem *> fields = getAllChildren("field") ;
		for(size_t i = 0 ; i < fields.size() ; i++)
		{
			bool isFieldType = true ;
			FieldType f = ConfigTreeItem::translateFieldType( fields[i]->getStringData(), isFieldType ) ;
			bool isTWFieldType = false ;
			TWFieldType g = ConfigTreeItem::translateTriangleWriterFieldType( fields[i]->getStringData(), isTWFieldType ) ;
			if(isFieldType)
			{
				trg->getField( f ) ;
			}
			else if(isTWFieldType)
			{
				trg->getField( g ) ;
			}
			else
				trg->getField(  fields[i]->getStringData() ) ;
		}
		trg->append() ;
		trg->writeSvg() ;
	}
}

#ifdef __WIN32
void ConfigTreeItem::makeWindowsPath() 
{
	for(size_t i = 0 ; i < children.size() ; i++)
		children[i]->makeWindowsPath() ;

	if(label != "function")
	{
		std::string path = str ;
		size_t backslash = path.find("/") ;
		while(backslash < std::string::npos)
		{
			std::string left = path.substr(0, backslash) ;
			std::string right = path.substr(backslash+1) ;
			path = left+"\\"+right ;
			backslash = path.find("/") ;
		}
		str = path ;
		std::cout << path << std::endl ;
	}

}
#endif

