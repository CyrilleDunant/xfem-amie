#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include "../utilities/parser.h"
#include "../physics/physics.h"
#include "../solvers/assembly.h"
#include "../elements/elements.h"
#include "../polynomial/variable.h"
#include "../geometry/geometry_3D.h"
#include "../sparse/sparse_matrix.h"

using namespace Mu ;

typedef enum
{
	NONE,
	DISPLACEMENT,
	LOAD,
	HOMOGENEOUS_X_DISPLACEMENT,
	BLOCKED_ALONG_X,
	BLOCKED_ALONG_Y,
	BLOCKED_ALONG_Z,
	PERIODIC_BOUNDARY_CONDITION,
} BoundaryCondition ;

typedef enum
{
	 CH = 1,
	CSH = 2,
	C3S = 3,
	ACIER = 5
} MaterialType ;

std::valarray<double> strainFromDisplacements(const std::valarray<double> & disp, const std::vector< Hexahedron*> & hex)
{
	
	HexahedralElement * father = new  HexahedralElement(LINEAR);
	
	std::valarray<double> strain(0., 6*hex.size()) ;
	
	for(size_t i  = 0 ; i < hex.size() ; i++)
	{
		HexahedralElement* e = new HexahedralElement(father, hex[i]) ;	
		
// 		std::valarray<Matrix> xi(2) ;
// 		xi[1][0][0] = 1 ;
// 		std::valarray<Matrix> eta(2) ;
// 		eta[0][1][0] = 1 ;
// 		std::valarray<Matrix> zeta(2) ;
// 		zeta[0][0][1] = 1 ;
// 		std::valarray<Matrix> one(2) ;
// 		one[0][0][0] = 1 ;
// 		
// 		Function x(xi) ;
// 		Function y(eta) ;
// 		Function z(zeta) ;
		
		VirtualMachine vm ;
		double s0_x = vm.eval(e->getShapeFunction(0).d(XI),Point()) ;
		double s0_y = vm.eval(e->getShapeFunction(0).d(ETA),Point()) ;
		double s0_z = vm.eval(e->getShapeFunction(0).d(ZETA),Point()) ;
		double s1_x = vm.eval(e->getShapeFunction(1).d(XI),Point()) ;
		double s1_y = vm.eval(e->getShapeFunction(1).d(ETA),Point()) ;
		double s1_z = vm.eval(e->getShapeFunction(1).d(ZETA),Point()) ;
		double s2_x = vm.eval(e->getShapeFunction(2).d(XI),Point()) ;
		double s2_y = vm.eval(e->getShapeFunction(2).d(ETA),Point()) ;
		double s2_z = vm.eval(e->getShapeFunction(2).d(ZETA),Point()) ;
		double s3_x = vm.eval(e->getShapeFunction(3).d(XI),Point()) ;
		double s3_y = vm.eval(e->getShapeFunction(3).d(ETA),Point()) ;
		double s3_z = vm.eval(e->getShapeFunction(3).d(ZETA),Point()) ;
		double s4_x = vm.eval(e->getShapeFunction(4).d(XI),Point()) ;
		double s4_y = vm.eval(e->getShapeFunction(4).d(ETA),Point()) ;
		double s4_z = vm.eval(e->getShapeFunction(4).d(ZETA),Point()) ;
		double s5_x = vm.eval(e->getShapeFunction(5).d(XI),Point()) ;
		double s5_y = vm.eval(e->getShapeFunction(5).d(ETA),Point()) ;
		double s5_z = vm.eval(e->getShapeFunction(5).d(ZETA),Point()) ;
		double s6_x = vm.eval(e->getShapeFunction(6).d(XI),Point()) ;
		double s6_y = vm.eval(e->getShapeFunction(6).d(ETA),Point()) ;
		double s6_z = vm.eval(e->getShapeFunction(6).d(ZETA),Point()) ;
		double s7_x = vm.eval(e->getShapeFunction(7).d(XI),Point()) ;
		double s7_y = vm.eval(e->getShapeFunction(7).d(ETA),Point()) ;
		double s7_z = vm.eval(e->getShapeFunction(7).d(ZETA),Point()) ;
		
		
		double xdxi = s0_x*disp[e->Hexahedron::getBoundingPoint(0).id*3] +  s1_x*disp[e->Hexahedron::getBoundingPoint(1).id*3] + s2_x*disp[e->Hexahedron::getBoundingPoint(2).id*3] + s3_x*disp[e->Hexahedron::getBoundingPoint(3).id*3]+ s4_x*disp[e->Hexahedron::getBoundingPoint(4).id*3]+ s5_x*disp[e->Hexahedron::getBoundingPoint(5).id*3]+ s6_x*disp[e->Hexahedron::getBoundingPoint(6).id*3]+ s7_x*disp[e->Hexahedron::getBoundingPoint(7).id*3];
		
		double xdeta =  s0_y*disp[e->Hexahedron::getBoundingPoint(0).id*3] +  s1_y*disp[e->Hexahedron::getBoundingPoint(1).id*3] + s2_y*disp[e->Hexahedron::getBoundingPoint(2).id*3] + s3_y*disp[e->Hexahedron::getBoundingPoint(3).id*3]+ s4_y*disp[e->Hexahedron::getBoundingPoint(4).id*3]+ s5_y*disp[e->Hexahedron::getBoundingPoint(5).id*3]+ s6_y*disp[e->Hexahedron::getBoundingPoint(6).id*3]+ s7_y*disp[e->Hexahedron::getBoundingPoint(7).id*3];
		
		double xdzeta = s0_z*disp[e->Hexahedron::getBoundingPoint(0).id*3] +  s1_z*disp[e->Hexahedron::getBoundingPoint(1).id*3] + s2_z*disp[e->Hexahedron::getBoundingPoint(2).id*3] + s3_z*disp[e->Hexahedron::getBoundingPoint(3).id*3]+ s4_z*disp[e->Hexahedron::getBoundingPoint(4).id*3]+ s5_z*disp[e->Hexahedron::getBoundingPoint(5).id*3]+ s6_z*disp[e->Hexahedron::getBoundingPoint(6).id*3]+ s7_z*disp[e->Hexahedron::getBoundingPoint(7).id*3];
		
		double ydxi = s0_x*disp[e->Hexahedron::getBoundingPoint(0).id*3+1] +  s1_x*disp[e->Hexahedron::getBoundingPoint(1).id*3+1] + s2_x*disp[e->Hexahedron::getBoundingPoint(2).id*3+1] + s3_x*disp[e->Hexahedron::getBoundingPoint(3).id*3+1]+ s4_x*disp[e->Hexahedron::getBoundingPoint(4).id*3+1]+ s5_x*disp[e->Hexahedron::getBoundingPoint(5).id*3+1]+ s6_x*disp[e->Hexahedron::getBoundingPoint(6).id*3+1]+ s7_x*disp[e->Hexahedron::getBoundingPoint(7).id*3+1];
		
		
		double ydeta =  s0_y*disp[e->Hexahedron::getBoundingPoint(0).id*3+1] +  s1_y*disp[e->Hexahedron::getBoundingPoint(1).id*3+1] + s2_y*disp[e->Hexahedron::getBoundingPoint(2).id*3+1] + s3_y*disp[e->Hexahedron::getBoundingPoint(3).id*3+1]+ s4_y*disp[e->Hexahedron::getBoundingPoint(4).id*3+1]+ s5_y*disp[e->Hexahedron::getBoundingPoint(5).id*3+1]+ s6_y*disp[e->Hexahedron::getBoundingPoint(6).id*3+1]+ s7_y*disp[e->Hexahedron::getBoundingPoint(7).id*3+1];
		
		
		double ydzeta =s0_z*disp[e->Hexahedron::getBoundingPoint(0).id*3+1] +  s1_z*disp[e->Hexahedron::getBoundingPoint(1).id*3+1] + s2_z*disp[e->Hexahedron::getBoundingPoint(2).id*3+1] + s3_z*disp[e->Hexahedron::getBoundingPoint(3).id*3+1]+ s4_z*disp[e->Hexahedron::getBoundingPoint(4).id*3+1]+ s5_z*disp[e->Hexahedron::getBoundingPoint(5).id*3+1]+ s6_z*disp[e->Hexahedron::getBoundingPoint(6).id*3+1]+ s7_z*disp[e->Hexahedron::getBoundingPoint(7).id*3+1]; 
		
		
		double zdxi =  s0_x*disp[e->Hexahedron::getBoundingPoint(0).id*3+2] +  s1_x*disp[e->Hexahedron::getBoundingPoint(1).id*3+2] + s2_x*disp[e->Hexahedron::getBoundingPoint(2).id*3+2] + s3_x*disp[e->Hexahedron::getBoundingPoint(3).id*3+2]+ s4_x*disp[e->Hexahedron::getBoundingPoint(4).id*3+2]+ s5_x*disp[e->Hexahedron::getBoundingPoint(5).id*3+2]+ s6_x*disp[e->Hexahedron::getBoundingPoint(6).id*3+2]+ s7_x*disp[e->Hexahedron::getBoundingPoint(7).id*3+2];
		
		
		double zdeta = s0_y*disp[e->Hexahedron::getBoundingPoint(0).id*3+2] +  s1_y*disp[e->Hexahedron::getBoundingPoint(1).id*3+2] + s2_y*disp[e->Hexahedron::getBoundingPoint(2).id*3+2] + s3_y*disp[e->Hexahedron::getBoundingPoint(3).id*3+2]+ s4_y*disp[e->Hexahedron::getBoundingPoint(4).id*3+2]+ s5_y*disp[e->Hexahedron::getBoundingPoint(5).id*3+2]+ s6_y*disp[e->Hexahedron::getBoundingPoint(6).id*3+2]+ s7_y*disp[e->Hexahedron::getBoundingPoint(7).id*3+2];
		
		 
		double zdzeta = s0_z*disp[e->Hexahedron::getBoundingPoint(0).id*3+2] +  s1_z*disp[e->Hexahedron::getBoundingPoint(1).id*3+2] + s2_z*disp[e->Hexahedron::getBoundingPoint(2).id*3+2] + s3_z*disp[e->Hexahedron::getBoundingPoint(3).id*3+2]+ s4_z*disp[e->Hexahedron::getBoundingPoint(4).id*3+2]+ s5_z*disp[e->Hexahedron::getBoundingPoint(5).id*3+2]+ s6_z*disp[e->Hexahedron::getBoundingPoint(6).id*3+2]+ s7_z*disp[e->Hexahedron::getBoundingPoint(7).id*3+2]; 
		
		
		
		
		Matrix Jinv ( e->getInverseJacobianMatrix(Point())) ; //any point would do : is constant on the element
		strain[i*6  ] = xdxi*Jinv[0][0] + xdeta*Jinv[0][1] + xdzeta*Jinv[0][2]  ;   // epsilon_xx
		strain[i*6+1] = ydxi*Jinv[1][0] + ydeta*Jinv[1][1] + ydzeta*Jinv[1][2]  ;   // epsilon_yy
		strain[i*6+2] = zdxi*Jinv[2][0] + zdeta*Jinv[2][1] + zdzeta*Jinv[2][2]  ;   // epsilon_zz
		strain[i*6+3] = (xdxi*Jinv[1][0] + xdeta*Jinv[1][1] + xdzeta*Jinv[1][2])*0.5 + 
			(ydxi*Jinv[0][0] + ydeta*Jinv[0][1] + ydzeta*Jinv[0][2])*0.5 ;      // epsilon_xy
		strain[i*6+4] = (xdxi*Jinv[2][0] + xdeta*Jinv[2][1] + xdzeta*Jinv[2][2])*0.5 + 
			(zdxi*Jinv[0][0] + zdeta*Jinv[0][1] + zdzeta*Jinv[0][2])*0.5 ;      // epsilon_xz
		strain[i*6+5] =  (ydxi*Jinv[2][0] + ydeta*Jinv[2][1] + ydzeta*Jinv[2][2])*0.5 + 
			(zdxi*Jinv[1][0] + zdeta*Jinv[1][1] + zdzeta*Jinv[1][2])*0.5 ;      // epsilon_yz
		
	}
	
	return strain ;
}

std::valarray<double> stressFromStrain(const std::valarray<double> & strain, const  Matrix & cg)
{
	return strain*cg ;
	
}

std::valarray<double> stressFromStrain(const std::valarray<double>*  strain, std::map<MaterialType, Matrix> cg,const  std::vector<std::pair<std::valarray<int>, size_t> > * t_data)
{
	std::valarray<double> ret(0., strain->size()) ;
	
	for(size_t i = 0 ; i < t_data->size() ; i++)
	{
		std::valarray<double> temp(0.,6) ;
		std::copy(&(*strain)[i*6], &(*strain)[i*6+6], &temp[0]) ;
		temp[3] *= 2 ;
		temp[4] *= 2 ;
		temp[5] *= 2 ;
		temp = stressFromStrain(temp, cg[(MaterialType)(*t_data)[i].second]);
		std::copy(&temp[0], &temp[6],&ret[i*6]) ;
	}
	
	return ret ;
}

std::valarray<double> averageStrain(const std::valarray<double> *displacements, const std::vector<Point *> * points, const std::vector < std::vector<std::valarray<int> > * > & boundPoints)
{
	//(*boundPoints[face_id][Rectangle_id])[point_id]
	
	std::valarray<double> ret(0., 6) ;
	
	for(size_t i = 0 ;  i < 6 ; i++)
	{
		Point normal ;
		
		if(i == 0)
			normal.set(-1, 0, 0) ;
		if(i==1)
			normal.set(1, 0, 0) ;
		if(i==2)
			normal.set(0, -1, 0) ;
		if(i==3)
			normal.set(0, 1, 0) ;
		if(i==4)
			normal.set(0, 0, -1) ;
		if(i==5)
			normal.set(0, 0, 1) ;
		
		for(size_t j = 0 ; j < boundPoints[i]->size() ; j++)
		{
			double area = Parallelogramme((*points)[(*boundPoints[i])[j][0]], (*points)[(*boundPoints[i])[j][1]], (*points)[(*boundPoints[i])[j][2]], (*points)[(*boundPoints[i])[j][3]]).area() ;
//			std::cout << area << std::endl;	
			for(size_t k = 0 ; k < 4 ; k++)
			{
				size_t this_point_id = (*boundPoints[i])[j][k] ;
				ret[0] += 2./4.*(*displacements)[this_point_id*3]*normal.x*area ; //00 -> xx
				ret[1] += 2./4.*(*displacements)[this_point_id*3+1]*normal.y*area ; //11 -> yy
				ret[2] += 2./4.*(*displacements)[this_point_id*3+2]*normal.z*area ; //22 -> zz
				ret[3] += 1./4.*((*displacements)[this_point_id*3]*normal.y + (*displacements)[this_point_id*3+1]*normal.x)*area;//01 ->xy
				ret[4] += 1./4.*((*displacements)[this_point_id*3]*normal.z  + (*displacements)[this_point_id*3+2]*normal.x)*area;//02 ->xz
				ret[5] += 1./4.*((*displacements)[this_point_id*3+1]*normal.z  + (*displacements)[this_point_id*3+2]*normal.y)*area;//12 ->yz
			}
		}
	}
	return ret ;
}

std::valarray<double> averageStress(const std::valarray<double> & stress, const std::vector<Hexahedron *> & hex)
{
	std::valarray<double> ret(0., 6) ;
	
	for(size_t i = 0 ; i < hex.size() ; i++)
	{
		for(size_t j = 0 ; j < 6 ; j++)
			ret[j] += stress[i*6+j]*hex[i]->volume() ;
		
	}
	
	return ret ;
}

double solidVolume(const std::vector< Hexahedron*> hex)
{
	double ret = 0 ;
	for(size_t i = 0 ; i < hex.size() ; i++)
	{
		ret += hex[i]->volume() ;
	}
	
	return ret ;
}

double totalVolume(const std::vector<Point *> * points)
{
	double max_x = (*points)[0]->x ;
	double min_x = (*points)[0]->x ;
	double max_y = (*points)[0]->y ;
	double min_y = (*points)[0]->y ;
	double max_z = (*points)[0]->z ;
	double min_z = (*points)[0]->z ;
	for(size_t i = 1 ; i < points->size() ; i++)
	{
		if((*points)[i]->x > max_x )
			max_x = (*points)[i]->x ;
		if((*points)[i]->x < min_x )
			min_x = (*points)[i]->x ;
		if((*points)[i]->y > max_y )
			max_y = (*points)[i]->y ;
		if((*points)[i]->y < min_y )
			min_y = (*points)[i]->y ;
		if((*points)[i]->z > max_z )
			max_z = (*points)[i]->z ;
		if((*points)[i]->z < min_z )
			min_z = (*points)[i]->z ;
	}
	
	return (max_x-min_x)*(max_y-min_y)*(max_z-min_z) ;
}

int main(size_t argc, char * argv[])
{
		
	double E_csh = 31000 ;      // in MPa
	double nu_csh = 0.28 ;
	
	double E_ch = 40000 ;      // in MPa
	double nu_ch = 0.3 ;
	
	double E_c3s = 135000 ;      // in MPa
	double nu_c3s = 0.31 ;
	
	double E_acier = 350000 ;      // in MPa
	double nu_acier = 0.18 ;
	double E_eff= 350000;
	double nu_eff=0.18 ;
	double K_eff_mat,K_acier;
	double G_acier,G_eff_mat;

	std::string pointFile("") ;
	bool got_p_file = false ;
	std::string hexFile("") ;
	bool got_t_file = false ;
	std::string boundFile("") ;
	bool got_b_file = false ;

//	Added here.
	bool check_periodic = false ;
//	Added here.

//	std::string periodFile("") ;
//	bool got_i_file = false ;
	for(size_t i = 1 ; i< argc ; i++)
	{
		if(std::string(argv[i]) == std::string("--point-file"))
		{
			pointFile = argv[++i] ;
			got_p_file = true ;
		}
		else if(std::string(argv[i]) == std::string("--hex-file"))
		{
			hexFile = argv[++i] ;
			got_t_file = true ;
		}
		else if(std::string(argv[i])== ("--cshrix-modulus"))
		{
			E_csh = atof(argv[++i]) ;
		}
		else if(std::string(argv[i])== ("--c3slusion-modulus"))
		{
			E_c3s = atof(argv[++i]) ;
		}
		else if(std::string(argv[i])== ("--cshrix-nu"))
		{
			nu_csh = atof(argv[++i]) ;
		}
		else if(std::string(argv[i])== ("--c3slusion-nu"))
		{
			nu_c3s = atof(argv[++i]) ;
		}
		else if(std::string(argv[i])== ("--boundary-file"))
		{
			boundFile = argv[++i] ;
			got_b_file = true ;
		}
//		else if(std::string(argv[i])== ("--periodicity-file"))
//		{
//			periodFile = argv[++i] ;
//			got_i_file = true ;
//		}
		else if(std::string(argv[i]) == std::string("--help"))
		{
			std::cout << "USAGE:" << std::endl ;
			std::cout << "  * --point-file FILE        set the file containing point coordinates and ids" << std::endl ;
			std::cout << "  * --hex-file FILE          set the file containing the neighbourhood inforcshion" << std::endl ;
			std::cout << "  * --boundary-file FILE     set the file containing the list of boundary points" << std::endl ;
//			std::cout << "  * --Periodic-file FILE     set the file containing the list of identity changement" << std::endl ;
			std::cout << "  * --cshrix-modulus VAL     set the cshrix modulus (default 20 GPa)" << std::endl ;
			std::cout << "  * --c3slusion-modulus VAL  set the c3slusions modulus (default 60 GPa)" << std::endl ;
			std::cout << "  * --cshrix-nu VAL          set cshrix poisson ratio (default 0.3)" << std::endl ;
			std::cout << "  * --c3slusion-nu VAL       set the c3slusions poisson ratio (default 0.18)" << std::endl ;
			std::cout <<  std::endl ;
			std::cout << "  * --help                   print this help and exit" << std::endl ;
			std::cout <<  std::endl ;
			return 0 ;
		}
		else
		{
			std::cerr << "invalid argument, use \"--help\" for more inforcshion" << std::endl ;
			return 1 ;
		}
	}
	
	if(!got_p_file)
	{
		std::cout << "please enter points file :" << std::flush ;
		std::cin >> pointFile ;
	}
	if(!got_t_file)
	{
		std::cout << "please enter hex file :" << std::flush ;
		std::cin >> hexFile ;
	}
	if(!got_b_file)
	{
		std::cout << "please enter boundary file :" << std::flush ;
		std::cin >> boundFile ;
	}
//	if(!got_i_file)
//	{
//		std::cout << "please enter periodicity file :" << std::flush ;
//		std::cin >> periodFile ;
//	}
	std::cout << "=== summary ===" << std::endl ;
	std::cout << "points     : " << pointFile << std::endl;
	std::cout << "hex       : " << hexFile << std::endl;
	std::cout << "boundary   : " << boundFile << std::endl;
//	std::cout << "periodicity   : " << periodFile << std::endl;

	//end cooking. Do stuff, now.
	
	PointParser3D pointsData(pointFile.c_str()) ;
	pointsData.readData() ;
	HexahedronParser hexData(hexFile.c_str()) ;
	hexData.readData() ;
	HexahedronBoundaryParser boundData(boundFile.c_str()) ;
	boundData.readData() ;
//	PeriodParser periodSetsData(periodFile.c_str()) ;
//	periodSetsData.readData() ;
	
	Matrix cgStressC3S(6,6) ;
	cgStressC3S[0][0] = 1. - nu_c3s ; cgStressC3S[0][1] = nu_c3s ; cgStressC3S[0][2] = nu_c3s ;
	cgStressC3S[1][0] = nu_c3s ; cgStressC3S[1][1] = 1. - nu_c3s ; cgStressC3S[1][2] = nu_c3s ;
	cgStressC3S[2][0] = nu_c3s ; cgStressC3S[2][1] = nu_c3s ; cgStressC3S[2][2] = 1. - nu_c3s ;
	cgStressC3S[3][3] = 0.5 - nu_c3s ;
	cgStressC3S[4][4] = 0.5 - nu_c3s ;
	cgStressC3S[5][5] = 0.5 - nu_c3s ;
	cgStressC3S *= E_c3s/((1.+nu_c3s)*(1.-2.*nu_c3s)) ;
	
	Matrix cgStressCSH(6,6) ;
	cgStressCSH[0][0] = 1. - nu_csh ; cgStressCSH[0][1] = nu_csh ; cgStressCSH[0][2] = nu_csh ;
	cgStressCSH[1][0] = nu_csh ; cgStressCSH[1][1] = 1. - nu_csh ; cgStressCSH[1][2] = nu_csh ;
	cgStressCSH[2][0] = nu_csh ; cgStressCSH[2][1] = nu_csh ; cgStressCSH[2][2] = 1. - nu_csh ;
	cgStressCSH[3][3] = 0.5 - nu_csh ;
	cgStressCSH[4][4] = 0.5 - nu_csh ;
	cgStressCSH[5][5] = 0.5 - nu_csh ;
	cgStressCSH *= E_csh/((1.+nu_csh)*(1.-2.*nu_csh)) ;
	
	Matrix cgStressCH(6,6) ;
	cgStressCH[0][0] = 1. - nu_ch ; cgStressCH[0][1] = nu_ch ; cgStressCH[0][2] = nu_ch ;
	cgStressCH[1][0] = nu_ch ; cgStressCH[1][1] = 1. - nu_ch ; cgStressCH[1][2] = nu_ch ;
	cgStressCH[2][0] = nu_ch ; cgStressCH[2][1] = nu_ch ; cgStressCH[2][2] = 1. - nu_ch ;
	cgStressCH[3][3] = 0.5 - nu_ch ;
	cgStressCH[4][4] = 0.5 - nu_ch ;
	cgStressCH[5][5] = 0.5 - nu_ch ;
	cgStressCH *= E_ch/((1.+nu_ch)*(1.-2.*nu_ch)) ;
	
	//Matrix cgStressACIER(6,6) ;
	//cgStressACIER[0][0] = 1 - nu_acier ; cgStressACIER[0][1] = nu_acier ; cgStressACIER[0][2] = nu_acier ;
	//cgStressACIER[1][0] = nu_acier ; cgStressACIER[1][1] = 1 - nu_acier ; cgStressACIER[1][2] = nu_acier ;
	//cgStressACIER[2][0] = nu_acier ; cgStressACIER[2][1] = nu_acier ; cgStressACIER[2][2] = 1 - nu_acier ;
	//cgStressACIER[3][3] = 0.5 - nu_acier ;
	//cgStressACIER[4][4] = 0.5 - nu_acier ;
	//cgStressACIER[5][5] = 0.5 - nu_acier ;
	//cgStressACIER *= E_acier/((1+nu_acier)*(1-2*nu_acier)) ;
	
// 	Stiffness csh_stiffness(cgStressCSH) ;
// 	Stiffness c3s_stiffness(cgStressC3S) ;
// 	Stiffness ch_stiffness(cgStressCH) ;
	//Stiffness acier_stiffness(cgStressACIER) ;
	
	//let us get the data
	std::vector<Point *> * points = pointsData.getData() ;
	std::vector<std::pair<std::valarray<int>, size_t> > * t_data = hexData.getData() ;
	std::vector<Hexahedron *> hex ;
	
	for(size_t i = 0 ;  i < t_data->size() ; i++)
	{
		hex.push_back(new Hexahedron((*points)[(*t_data)[i].first[0]], (*points)[(*t_data)[i].first[1]], (*points)[(*t_data)[i].first[2]], (*points)[(*t_data)[i].first[3]], (*points)[(*t_data)[i].first[4]], (*points)[(*t_data)[i].first[5]], (*points)[(*t_data)[i].first[6]], (*points)[(*t_data)[i].first[7]])) ;
	}
	
	
	std::cout << "=== summary ===" << std::endl ;
	std::cout << "points     : " << points->size() << std::endl;
	std::cout << "hex       : " << hex.size() << std::endl;
	
	HexahedralElement * father = new  HexahedralElement(LINEAR);

	
	//ask the user what to do
	
	std::valarray<BoundaryCondition> bc(NONE, 7) ;
	std::map<size_t, std::valarray<double> *> face_load ;
	std::map<size_t, std::valarray<double> * > face_displacement ;
	std::valarray<double> b ;
	std::valarray<double> x ;
	double displacement_factor = 1 ;
	
	//bool quit = false ; 
	
	//while(!quit)
	//{

 //     We need some physical data, so we set the default
	
// 	Hexahedron toto ;
// 	
// 	std::cout << toto.volume() << std::endl ;
// 	std::cout << toto.area() << std::endl ;
// 	return 0 ;

		
		//Assembly3D assembly(&csh_stiffness) ;
		
		for(int i = 0 ; i < 6 ; i++)
		{
			std::cout << "instruction for face " << i+1 << ": " << std::endl ;
			std::cout << "      [0]  No cond " << std::endl ;
			std::cout << "      [1]  Apply displacement " << std::endl ;
			std::cout << "      [2]  Apply load " << std::endl ;
			std::cout << "      [4]  Set X " << std::endl ;
			std::cout << "      [5]  Set Y " << std::endl ;
			std::cout << "      [6]  Set Z " << std::endl ;
			std::cout << "      [7]  Periodic Boundary Condition " << std::endl ;
			std::cout << "                      " << std::endl ;
			std::cout << "      [3]  homogeneous x displacement " << std::endl ;
			std::cout << " choice : " << std::flush ;
			
			size_t choice = DISPLACEMENT;
			std::cin >> choice ;
			if(choice == HOMOGENEOUS_X_DISPLACEMENT)
			{
				bc[i] = (BoundaryCondition)choice ;
				std::cout << " factor ? " << std::flush ;
				std::cin >> displacement_factor ;
				break ;
			}
			
			bc[i] = (BoundaryCondition)choice ;

//	Added here.
			if(choice == PERIODIC_BOUNDARY_CONDITION)
			{
				check_periodic = true ;
				std::valarray<double> val(0., 3) ;
				
			std::cout << " in x : value : " << std::flush ;
				std::cin >> val[0] ;
			std::cout << " in y : value : " << std::flush ;
				std::cin >> val[1] ;
			std::cout << " in z : value : " << std::flush ;
				std::cin >> val[2] ;
				
				face_load[i] =new std::valarray<double>(val) ;
				face_load[i+1] =new std::valarray<double>(val) ;
				i = i + 1 ;
			}
//	Added here.
			
			if(choice == LOAD)
			{
				std::valarray<double> val(0., 3) ;
				
			std::cout << " in x : value : " << std::flush ;
				std::cin >> val[0] ;
			std::cout << " in y : value : " << std::flush ;
				std::cin >> val[1] ;
			std::cout << " in z : value : " << std::flush ;
				std::cin >> val[2] ;
				
				face_load[i] =new std::valarray<double>(val) ;
			}
			if(choice == DISPLACEMENT)
			{
				
				std::valarray<double> val(0., 3) ;
				
			std::cout << " in x : value : " << std::flush ;
				std::cin >> val[0] ;
			std::cout << " in y : value : " << std::flush ;
				std::cin >> val[1] ;
			std::cout << " in z : value : " << std::flush ;
				std::cin >> val[2] ;
				
				face_displacement[i] =new std::valarray<double>(val) ;
			}
			
			if(choice == BLOCKED_ALONG_X || choice == BLOCKED_ALONG_Y || choice == BLOCKED_ALONG_Z)
			{
				
				std::valarray<double> val(0., 3) ;
				
			std::cout << " value : " << std::flush ;
				std::cin >> val[0] ;
				
				face_displacement[i] =new std::valarray<double>(val) ;
			}
			
			if(i == 5)
			{
				std::cout << "=== summary ===" << std::endl ;
				
				if(bc[0] == LOAD)
				{
					std::cout << "* face 1 *" << std::endl;
					std::cout << "load = " << (*face_load[0])[0] << ", "<<  (*face_load[0])[1] << ", "<<   (*face_load[0])[2] << std::endl ;
				}

//	Added here.
				else if(bc[0] == PERIODIC_BOUNDARY_CONDITION)
				{
					bc[1] = PERIODIC_BOUNDARY_CONDITION ;
					std::cout << "* face 1 *" << std::endl;
					std::cout << "load = " << (*face_load[0])[0] << ", "<<  (*face_load[0])[1] << ", "<<   (*face_load[0])[2] << std::endl ;

					std::cout << "* face 2 *" << std::endl;
					std::cout << "load = " << (*face_load[1])[0] << ", "<<  (*face_load[1])[1] << ", "<<   (*face_load[1])[2] << std::endl ;

				}
//	Added here.

				else if(bc[0] == DISPLACEMENT)
				{
					std::cout << "* face 1 *" << std::endl;
					std::cout << "displacement = " << (*face_displacement[0])[0] << ", "<<  (*face_displacement[0])[1] << ", "<<   (*face_displacement[0])[2] << std::endl ;
				}
				else if(bc[0] == BLOCKED_ALONG_X)
				{
					std::cout << "* face 1 *" << std::endl;
					std::cout << "displacement x = " << (*face_displacement[0])[0] << std::endl ;
				}
				else if(bc[0] == BLOCKED_ALONG_Y)
				{
					std::cout << "* face 1 *" << std::endl;
					std::cout << "displacement y = " << (*face_displacement[0])[0] << std::endl ;
				}
				else if(bc[0] == BLOCKED_ALONG_Z)
				{
					std::cout << "* face 1 *" << std::endl;
					std::cout << "displacement z = " << (*face_displacement[0])[0] << std::endl ;
				}
				
				if(bc[1] == LOAD)
				{
					std::cout << "* face 2 *" << std::endl;
					std::cout << "load = " << (*face_load[1])[0] << ", "<<  (*face_load[1])[1] << ", "<<   (*face_load[1])[2] << std::endl ;
				}
				else if(bc[1] == DISPLACEMENT)
				{
					std::cout << "* face 2 *" << std::endl;
					std::cout << "displacement = " <<  (*face_displacement[1])[0] << ", "<<  (*face_displacement[1])[1] << ", "<<   (*face_displacement[1])[2] << std::endl ;
				}
				else if(bc[1] == BLOCKED_ALONG_X)
				{
					std::cout << "* face 2 *" << std::endl;
					std::cout << "displacement x= " << (*face_displacement[1])[0] << std::endl ;
				}
				else if(bc[1] == BLOCKED_ALONG_Y)
				{
					std::cout << "* face 2 *" << std::endl;
					std::cout << "displacement y = " << (*face_displacement[1])[0] << std::endl ;
				}
				else if(bc[1] == BLOCKED_ALONG_Z)
				{
					std::cout << "* face 2 *" << std::endl;
					std::cout << "displacement z = " << (*face_displacement[1])[0] << std::endl ;
				}
				
				if(bc[2] == LOAD)
				{
					std::cout << "* face 3 *" << std::endl;
					std::cout << "load = " << (*face_load[2])[0] << ", "<<  (*face_load[2])[1] << ", "<<   (*face_load[2])[2] << std::endl ;
				}

//	Added here.
				else if(bc[2] == PERIODIC_BOUNDARY_CONDITION)
				{
					bc[3] = PERIODIC_BOUNDARY_CONDITION ;
					std::cout << "* face 3 *" << std::endl;
					std::cout << "load = " << (*face_load[2])[0] << ", "<<  (*face_load[2])[1] << ", "<<   (*face_load[2])[2] << std::endl ;

					std::cout << "* face 4 *" << std::endl;
					std::cout << "load = " << (*face_load[3])[0] << ", "<<  (*face_load[3])[1] << ", "<<   (*face_load[3])[2] << std::endl ;

				}
//	Added here.

				else if(bc[2] == DISPLACEMENT)
				{
					std::cout << "* face 3 *" << std::endl;
					std::cout << "displacement = " <<  (*face_displacement[2])[0] << ", "<<  (*face_displacement[2])[1] << ", "<<   (*face_displacement[2])[2] << std::endl ;
				}
				else if(bc[2] == BLOCKED_ALONG_X)
				{
					std::cout << "* face 3 *" << std::endl;
					std::cout << "displacement x = " << (*face_displacement[2])[0] << std::endl ;
				}
				else if(bc[2] == BLOCKED_ALONG_Y)
				{
					std::cout << "* face 3 *" << std::endl;
					std::cout << "displacement y = " << (*face_displacement[2])[0] << std::endl ;
				}
				else if(bc[2] == BLOCKED_ALONG_Z)
				{
					std::cout << "* face 3 *" << std::endl;
					std::cout << "displacement z = " << (*face_displacement[2])[0] << std::endl ;
				}
				
				if(bc[3] == LOAD)
				{
					std::cout << "* face 4 * " << std::endl;
					std::cout << "load = " << (*face_load[3])[0] << ", "<<  (*face_load[3])[1] << ", "<<   (*face_load[3])[2] << std::endl ;
				}
				else if(bc[3] == DISPLACEMENT)
				{
					std::cout << "* face 4 * " << std::endl;
					std::cout << "displacement = " <<  (*face_displacement[3])[0] << ", "<<  (*face_displacement[3])[1] << ", "<<   (*face_displacement[3])[2] << std::endl ;
				}
				else if(bc[3] == BLOCKED_ALONG_X)
				{
					std::cout << "* face 4 *" << std::endl;
					std::cout << "displacement x = " << (*face_displacement[3])[0] << std::endl ;
				}
				else if(bc[3] == BLOCKED_ALONG_Y)
				{
					std::cout << "* face 4 *" << std::endl;
					std::cout << "displacement y = " << (*face_displacement[3])[0] << std::endl ;
				}
				else if(bc[3] == BLOCKED_ALONG_Z)
				{
					std::cout << "* face 4 *" << std::endl;
					std::cout << "displacement z = " << (*face_displacement[3])[0] << std::endl ;
				}
				
				if(bc[4] == LOAD)
				{
					std::cout << "* face 5 * " << std::endl;
					std::cout << "load = " << (*face_load[4])[0] << ", "<<  (*face_load[4])[1] << ", "<<   (*face_load[4])[2] << std::endl ;
				}

//	Added here.
				else if(bc[4] == PERIODIC_BOUNDARY_CONDITION)
				{
					bc[5] = PERIODIC_BOUNDARY_CONDITION ;
					std::cout << "* face 5 *" << std::endl;
					std::cout << "load = " << (*face_load[4])[0] << ", "<<  (*face_load[4])[1] << ", "<<   (*face_load[4])[2] << std::endl ;

					std::cout << "* face 6 *" << std::endl;
					std::cout << "load = " << (*face_load[5])[0] << ", "<<  (*face_load[5])[1] << ", "<<   (*face_load[5])[2] << std::endl ;

				}
//	Added here.

				else if(bc[4] == DISPLACEMENT)
				{
					std::cout << "* face 5 * " << std::endl;
					std::cout << "displacement = " <<  (*face_displacement[4])[0] << ", "<<  (*face_displacement[4])[1] << ", "<<   (*face_displacement[4])[2] << std::endl ;
				}
				else if(bc[4] == BLOCKED_ALONG_X)
				{
					std::cout << "* face 5 *" << std::endl;
					std::cout << "displacement x = " << (*face_displacement[4])[0] << std::endl ;
				}
				else if(bc[4] == BLOCKED_ALONG_Y)
				{
					std::cout << "* face 5 *" << std::endl;
					std::cout << "displacement y = " << (*face_displacement[4])[0] << std::endl ;
				}
				else if(bc[4] == BLOCKED_ALONG_Z)
				{
					std::cout << "* face 5 *" << std::endl;
					std::cout << "displacement z = " << (*face_displacement[4])[0] << std::endl ;
				}
				
				if(bc[5] == LOAD)
				{
					std::cout << "* face 6 * " << std::endl;
					std::cout << "load = " << (*face_load[5])[0] << ", "<<  (*face_load[5])[1] << ", "<<   (*face_load[5])[2] << std::endl ;
				}
				else if(bc[5] == DISPLACEMENT)
				{
					std::cout << "* face 6 * " << std::endl;
					std::cout << "displacement = " <<  (*face_displacement[5])[0] << ", "<<  (*face_displacement[5])[1] << ", "<<   (*face_displacement[5])[2] << std::endl ;
				}
				else if(bc[5] == BLOCKED_ALONG_X)
				{
					std::cout << "* face 6 *" << std::endl;
					std::cout << "displacement x = " << (*face_displacement[5])[0] << std::endl ;
				}
				else if(bc[5] == BLOCKED_ALONG_Y)
				{
					std::cout << "* face 6 *" << std::endl;
					std::cout << "displacement y = " << (*face_displacement[5])[0] << std::endl ;
				}
				else if(bc[5] == BLOCKED_ALONG_Z)
				{
					std::cout << "* face 6 *" << std::endl;
					std::cout << "displacement z = " << (*face_displacement[5])[0] << std::endl ;
				}
				
				bool yes = true ;
				
				std::cout << "is this correct ?" << std::flush ;
				std::cin >> yes ;
				if(!yes)
					i = -1 ;
			}
		}
	
	
	int l=0;
	double eps = 0.01;
	double Av_Stress_mat_xx;
	double Av_Stress_mat_yy;
	double Av_Stress_mat_zz;
	double Av_Strain_mat_xx;
	double Av_Strain_mat_yy;
	double Av_Strain_mat_zz;
	std::cout<< " eps "<< eps <<std::endl;
	do
	{
	
	 E_acier= E_eff;
	nu_acier = nu_eff;
	Matrix cgStressACIER(6,6) ;
	cgStressACIER[0][0] = 1. - nu_acier ; cgStressACIER[0][1] = nu_acier ; cgStressACIER[0][2] = nu_acier ;
	cgStressACIER[1][0] = nu_acier ; cgStressACIER[1][1] = 1. - nu_acier ; cgStressACIER[1][2] = nu_acier ;
	cgStressACIER[2][0] = nu_acier ; cgStressACIER[2][1] = nu_acier ; cgStressACIER[2][2] = 1. - nu_acier ;
	cgStressACIER[3][3] = 0.5 - nu_acier ;
	cgStressACIER[4][4] = 0.5 - nu_acier ;
	cgStressACIER[5][5] = 0.5 - nu_acier ;
	cgStressACIER *= E_acier/((1.+nu_acier)*(1.-2.*nu_acier)) ;
// 	Stiffness acier_stiffness(cgStressACIER) ;
	
	Assembly assembly ;
	
//	Added the if condition.

	if(check_periodic == true)
	{
		for(size_t i = 0 ; i < pointsData.getPeriodicIds().size() ; i++)
		{
			(*points)[i]->id = pointsData.getPeriodicIds()[i] ;
		}
	}

//	Added the if condition.
		
		for(size_t i = 0 ; i < hex.size() ; i++)
		{
	
			HexahedralElement*e = new HexahedralElement(father, hex[i]) ;
			
			if((*t_data)[i].second == CSH)
			{
				e->setBehaviour(new Stiffness(cgStressCSH)) ;
				
			}
			else if ((*t_data)[i].second == C3S)
			{
				e->setBehaviour(new Stiffness(cgStressC3S)) ;
				
			}
			else if ((*t_data)[i].second == CH)
			{
				e->setBehaviour(new Stiffness(cgStressCH)) ;
				
			}
			else if ((*t_data)[i].second == ACIER)
			{
				e->setBehaviour(new Stiffness(cgStressACIER)) ;
				
			}
			else
			{
				std::cout << "Arrrgh we do not know our material, we are doomed! we will die!" << std::endl ;
				return 0 ;
			}
			assembly.add(e) ;
			
// 			std::cout << e->volume() << std::endl ;
			std::cout << "\r assembling " << i+1 << "/" << hex.size() << std::flush ;
		}
		std::cout << std::endl;
		
		// //now we apply the boundary conditions
		
		//first, we get the list of point id's in a convenient forcsh.
		std::vector < std::vector<std::valarray<int> > * > points_on_face_id ;
		
		for(size_t i = 0 ; i < 6 ; i++)
		{
			points_on_face_id.push_back(boundData.getData(i)) ;
		}
		 
//  		b.resize(points->size()*3, 0.) ;
		x.resize(points->size()*3, 0.) ;
		
		
		if(bc[0] == HOMOGENEOUS_X_DISPLACEMENT )
		{
			std::vector<size_t > point_list(6) ;
			for(size_t i = 0 ; i < points_on_face_id.size() ; i++) //face
			{
			   for(size_t j = 0 ; j < points_on_face_id[i]->size() ; j++) //carré
				{
				  for(size_t k = 0 ; k < (*points_on_face_id[i])[j].size() ; k++) //point
					{
					 point_list.push_back((*points_on_face_id[i])[j][k]) ;
				
					}
				}
			}
			
// 			std::stable_sort(point_list.begin(), point_list.end()) ;
// 			std::vector<size_t>::iterator e = std::unique(point_list.begin(), point_list.end()) ;
// 			point_list.erase(e, point_list.end()) ;
			
// 			for(size_t i = 0 ; i < point_list.size() ; i++)
// 			{	
// 				point_list[i] = pointsData.getPeriodicIds()[point_list[i]];
// 			}
			
			std::stable_sort(point_list.begin(), point_list.end()) ;
			std::vector<size_t>::iterator e = std::unique(point_list.begin(), point_list.end()) ;
			point_list.erase(e, point_list.end()) ;

 			for(size_t i = 0 ; i < point_list.size() ; i++)
 			{	
 				point_list[i] = pointsData.getPeriodicIds()[point_list[i]];
// 				std::cout << i << std::endl ;
 //				std::cout << "point_list[i]"<<  point_list[i]<<std::endl ;
			}
			size_t point_counter = 0 ;
			
			// CECI ESt IMPORTANT
			for(size_t i = 0 ; i < point_list.size() ; i++) //face
			{
				assembly.setPoint(
				                   (*points)[point_list[i]]->x*(displacement_factor-1) ,                   // x
				                   0,                                                                      // y
				                   0,                                                                      // z
// 				                   &b,                                                                     // forces
				                   point_list[i]                                                          // point ID
				                 ) ;
				std::cout << "\r setting displacement BC. point " << ++point_counter << std::flush ;
// 				std::cout << "point_list[i]" << point_list[i] << std::flush ;
// 				std::cout << "&b" << &b << std::flush ;
			}
			std::cout << std::endl;
			
		}
		else
		{
			
			//reminder : we have (*points_on_face_id[FACE])[Parallelogramme][POINT]
			
			
			//we do two passes, one to build the point list, the other to apply the conditions.
			
			//the variables ;
			std::vector<std::vector<size_t> > point_list(6) ;
			size_t point_counter = 0 ;

			//first pass
			
			for(size_t i = 0 ; i < points_on_face_id.size() ; i++) //face
			{
				for(size_t j = 0 ; j < points_on_face_id[i]->size() ; j++) //carré
				{
					for(size_t k = 0 ; k < (*points_on_face_id[i])[j].size() ; k++) //point
					{
						if(bc[i] == DISPLACEMENT || 
						   bc[i] == BLOCKED_ALONG_X || 
						   bc[i] == BLOCKED_ALONG_Y || 
						   bc[i] == BLOCKED_ALONG_Z )
						{
							point_list[i].push_back((*points_on_face_id[i])[j][k]) ;
						}
					}
				}
			}
			//making unique
				
 			for(size_t i = 0 ; i < point_list.size() ; i++)
			{
				std::stable_sort(point_list[i].begin(), point_list[i].end()) ;
				std::vector<size_t>::iterator e = std::unique(point_list[i].begin(), point_list[i].end()) ;
 				point_list[i].erase(e, point_list[i].end()) ;
			}

	

// 				for(size_t i = 0 ; i < point_list.size() ; i++)//face
// 				{
// 						std::cout << "point_list[i].size = " << point_list[i].size() << std::endl ;
// 						for(size_t j = 0 ; j < point_list[i].size() ; j++)//point
// 						{
// 							point_list[i][j] = pointsData.getPeriodicIds()[point_list[i][j]];
// 						}
// 				}




			//second pass -- displacements
			for(size_t i = 0 ; i < point_list.size() ; i++) //face
			{
				for(size_t j = 0 ; j < point_list[i].size() ; j++) //point
				{
					if(bc[i] == DISPLACEMENT)
					{
						assembly.setPoint(
											(*face_displacement[i])[0] ,
											(*face_displacement[i])[1],
											(*face_displacement[i])[2], 
											
// 											&b,
										point_list[i][j]
											) ;
						std::cout << "\r setting displacement BC. point " << ++point_counter << std::flush ;
					}
					if(bc[i] == BLOCKED_ALONG_X)
					{
						assembly.setPointAlong(XI, (*face_displacement[i])[0],  point_list[i][j]) ;
						std::cout << "\r setting displacement BC. point " << ++point_counter << std::flush ;
					}
					if(bc[i] == BLOCKED_ALONG_Y)
					{
						assembly.setPointAlong(ETA, (*face_displacement[i])[0], point_list[i][j]) ;
						std::cout << "\r setting displacement BC. point " << ++point_counter << std::flush ;
					}
					if(bc[i] == BLOCKED_ALONG_Z)
					{
						assembly.setPointAlong(ZETA, (*face_displacement[i])[0],  point_list[i][j]) ;
						std::cout << "\r setting displacement BC. point " << ++point_counter << std::flush ;
					}
				}
			}
			
			std::cout << std::endl ;
			
			//second pass -- loads
//	Added here.

			if(check_periodic == true)
			{
				assembly.setPoint(
								0.,                                                                     // x
								0.,                                                                      // y
								0.,                                                                      // z
								435                                                         // point ID
								) ;
				std::cout << "Zero Displ. point is = 435 "  << std::endl ;
			}

//	Added here.
			for(size_t i = 0 ; i < points_on_face_id.size() ; i++) //face
			{
				if(bc[i] == LOAD)
				{
					for(size_t j = 0 ; j < points_on_face_id[i]->size() ; j++) //carré
					{
						double area = Parallelogramme((*points)[(*points_on_face_id[i])[j][0]], (*points)[(*points_on_face_id[i])[j][1]], (*points)[(*points_on_face_id[i])[j][2]], (*points)[(*points_on_face_id[i])[j][3]]).area() ;
						
						for(size_t k = 0 ; k < (*points_on_face_id[i])[j].size() ; k++) //point
						{
							assembly.setForceOn(XI,(*face_load[i])[0]*area/4,(*points_on_face_id[i])[j][k]) ;
							assembly.setForceOn(ETA,(*face_load[i])[1]*area/4,(*points_on_face_id[i])[j][k]) ;
							assembly.setForceOn(ZETA,(*face_load[i])[2]*area/4,(*points_on_face_id[i])[j][k]) ;
// 							b[(*points_on_face_id[i])[j][k]*3] += (*face_load[i])[0]*area/4 ;
// 							b[(*points_on_face_id[i])[j][k]*3+1] += (*face_load[i])[1]*area/4 ;
// 							b[(*points_on_face_id[i])[j][k]*3+2] += (*face_load[i])[2]*area/4 ;
							
								std::cout << "\r setting load BC. point " << ++point_counter << std::flush ;
						}
					}
				}
//	Added here.
				else if(bc[i] == PERIODIC_BOUNDARY_CONDITION)
				{
					for(size_t j = 0 ; j < points_on_face_id[i]->size() ; j++) //carré
					{
						double area = Parallelogramme((*points)[(*points_on_face_id[i])[j][0]], (*points)[(*points_on_face_id[i])[j][1]], (*points)[(*points_on_face_id[i])[j][2]], (*points)[(*points_on_face_id[i])[j][3]]).area() ;
						
						for(size_t k = 0 ; k < (*points_on_face_id[i])[j].size() ; k++) //point
						{
							assembly.setForceOn(XI,(*face_load[i])[0]*area/4,pointsData.getPeriodicIds()[(*points_on_face_id[i])[j][k]]) ;
							assembly.setForceOn(ETA,(*face_load[i])[1]*area/4,pointsData.getPeriodicIds()[(*points_on_face_id[i])[j][k]]) ;
							assembly.setForceOn(ZETA,(*face_load[i])[2]*area/4,pointsData.getPeriodicIds()[(*points_on_face_id[i])[j][k]]) ;
// 							b[(*points_on_face_id[i])[j][k]*3] += (*face_load[i])[0]*area/4 ;
// 							b[(*points_on_face_id[i])[j][k]*3+1] += (*face_load[i])[1]*area/4 ;
// 							b[(*points_on_face_id[i])[j][k]*3+2] += (*face_load[i])[2]*area/4 ;
							
								std::cout << "\r setting periodic boundary condition load BC. point " << ++point_counter << std::flush ;
						}
					}
				}
//	Added here.
			}
			
			std::cout << std::endl ;
		}
		
//  	    std::cout << "coucou" << std::endl;
		x = assembly.cgsolve(x, x.size());

// (b, x, x.size()) ;
		
		std::valarray<double> strains = strainFromDisplacements(x, hex) ;
		
		std::map<MaterialType, Matrix> cg ;
		cg[C3S] = cgStressC3S ;
		cg[CSH] = cgStressCSH ;
		cg[ACIER] = cgStressACIER ;
		cg[CH] = cgStressCH ;
		std::valarray<double> stress = stressFromStrain(&strains, cg, t_data) ;
		
		std::valarray<double> avgStrain = averageStrain(&x, points, points_on_face_id) ;
		std::valarray<double> avgStress = averageStress(stress, hex) ;
		std::valarray<double> avgStrainAcier(0., 6) ;
		std::valarray<double> avgStressAcier(0., 6) ;
		
		double volume_ch =0;
		double volume_csh =0;
		double volume_c3s =0;
		double volume_acier =0;
		
		for(size_t i = 0 ; i < hex.size() ; i++)
		{
			if((*t_data)[i].second == CH)
			{
				volume_ch += hex[i]->volume() ;
			}
			else if ((*t_data)[i].second == C3S)
			{
				volume_c3s += hex[i]->volume() ;
			}
			else if ((*t_data)[i].second == CSH)
			{
				volume_csh += hex[i]->volume() ;
			}
			else if ((*t_data)[i].second == ACIER)
			{
				volume_acier += hex[i]->volume() ;
				avgStressAcier[0] += stress[i*6]*hex[i]->volume() ;
				avgStrainAcier[0] += strains[i*6]*hex[i]->volume() ;
				avgStressAcier[1] += stress[i*6+1]*hex[i]->volume() ;
				avgStrainAcier[1] += strains[i*6+1]*hex[i]->volume() ;
				avgStressAcier[2] += stress[i*6+2]*hex[i]->volume() ;
				avgStrainAcier[2] += strains[i*6+2]*hex[i]->volume() ;
				avgStressAcier[3] += stress[i*6+3]*hex[i]->volume() ;
				avgStrainAcier[3] += strains[i*6+3]*hex[i]->volume() ;
				avgStressAcier[4] += stress[i*6+4]*hex[i]->volume() ;
				avgStrainAcier[4] += strains[i*6+4]*hex[i]->volume() ;
				avgStressAcier[5] += stress[i*6+5]*hex[i]->volume() ;
				avgStrainAcier[5] += strains[i*6+5]*hex[i]->volume() ;
			}
		}
		
		double V_solid = solidVolume(hex) ;
		double V_tot = totalVolume(points) ;
		avgStrain /= 2*V_tot ;
		avgStress /= V_tot ;
		avgStrainAcier /= volume_acier ;
		avgStressAcier /= volume_acier ;
		
		std::cout << "volume CH = " << volume_ch 
			<< ", volume C3S = " <<  volume_c3s 
			<< ", volume CSH = " << volume_csh 
			<< ", volume Acier = " << volume_acier 
			<< ", volume totale = " << V_tot 
			<< std::endl ;
		std::cout << "concentration CH = " << volume_ch/(V_tot-volume_acier) 
			<< ", concentration C3S = " << volume_c3s/(V_tot-volume_acier) 
			<< ", concentration CSH = " << volume_csh/(V_tot-volume_acier) 
			<< ", concentration Acier = " << volume_acier/V_tot 
			<< std::endl ;
		std::cout << "Concentration Solide = " << (V_solid-volume_acier)/(V_tot-volume_acier) << std::endl ;
		std::cout << std::endl ;
		std::cout << "average strain = " 
			<< " xx = " << avgStrain[0]
			<< ", yy = " << avgStrain[1] 
			<< ", zz = " << avgStrain[2] 
			<< std::endl ;
		std::cout << "average strain = " 
			<< " xy = " << avgStrain[3] 
			<< ", xz = " << avgStrain[4] 
			<< ", yz = " << avgStrain[5] 
			<< std::endl ;
		std::cout << std::endl ;
		std::cout << "average stress = " 
			<< " xx = " << avgStress[0] 
			<< ", yy = " << avgStress[1]
			<< ", zz = " << avgStress[2] 
			<< std::endl ;
		std::cout << "average stress = " 
			<< " xy = " << avgStress[3] 
			<< ", xz = " << avgStress[4] 
			<< ", yz = " << avgStress[5] 
			<< std::endl ;
		std::cout << std::endl ;


		std::cout << std::endl ;
		std::cout << std::endl ;
		
		// Modulus calculate

		std::cout << "average strain (ohne stahl) = " 
			<< " xx = " << (avgStrain[0] - (volume_acier/V_tot)*avgStrainAcier[0])/(1-(volume_acier/V_tot))
			<< ", yy = " <<( avgStrain[1] - (volume_acier/V_tot)*avgStrainAcier[1])/(1-(volume_acier/V_tot))
			<< ", zz = " << (avgStrain[2] -  (volume_acier/V_tot)*avgStrainAcier[2])/(1-(volume_acier/V_tot))
			<< std::endl ;
		std::cout << "average strain (ohne stahl) = " 
			<< " xy = " << (avgStrain[3] - (volume_acier/V_tot)*avgStrainAcier[3])/(1-(volume_acier/V_tot))
			<< ", xz = " << (avgStrain[4] - (volume_acier/V_tot)*avgStrainAcier[4])/(1-(volume_acier/V_tot))
			<< ", yz = " << (avgStrain[5] - (volume_acier/V_tot)*avgStrainAcier[5])/(1-(volume_acier/V_tot))
			<< std::endl ;
		std::cout << std::endl ;
		std::cout << "average stress (ohne stahl) = " 
			<< " xx = " << (avgStress[0] - (volume_acier/V_tot)*avgStressAcier[0])/(1-(volume_acier/V_tot)) 
			<< ", yy = " <<( avgStress[1] - (volume_acier/V_tot)*avgStressAcier[1])/(1-(volume_acier/V_tot))
			<< ", zz = " << (avgStress[2] -  (volume_acier/V_tot)*avgStressAcier[2])/(1-(volume_acier/V_tot))
			<< std::endl ;
		std::cout << "average stress (ohne stahl) = " 
			<< " xy = " << (avgStress[3] - (volume_acier/V_tot)*avgStressAcier[3])/(1-(volume_acier/V_tot))
			<< ", xz = " << (avgStress[4] - (volume_acier/V_tot)*avgStressAcier[4])/(1-(volume_acier/V_tot))
			<< ", yz = " << (avgStress[5] - (volume_acier/V_tot)*avgStressAcier[5])/(1-(volume_acier/V_tot))
			<< std::endl ;
		std::cout << std::endl ;
		std::cout << std::endl ;
		std::cout << "average strain acier= " << " xx = " << avgStrainAcier[0] << ", yy = " << avgStrainAcier[1] << ", zz = " << avgStrainAcier[2] << std::endl ;
		std::cout << "average strain acier= " << " xy = " << avgStrainAcier[3] << ", xz = " << avgStrainAcier[4] << ", yz = " << avgStrainAcier[5] << std::endl ;
		std::cout << std::endl ;
		std::cout << "average stress acier= " << " xx = " << avgStressAcier[0] << ", yy = " << avgStressAcier[1] << ", zz = " << avgStressAcier[2] << std::endl ;
		std::cout << "average stress acier= " << " xy = " << avgStressAcier[3] << ", xz = " << avgStressAcier[4] << ", yz = " << avgStressAcier[5] << std::endl ;

		// Properties calculate
		Av_Stress_mat_xx = (avgStress[0] - (volume_acier/V_tot)*avgStressAcier[0])/(1-(volume_acier/V_tot)) ;
		Av_Stress_mat_yy= (avgStress[1] - (volume_acier/V_tot)*avgStressAcier[1])/(1-(volume_acier/V_tot));
		Av_Stress_mat_zz= (avgStress[2] -  (volume_acier/V_tot)*avgStressAcier[2])/(1-(volume_acier/V_tot));
		Av_Strain_mat_xx =(avgStrain[0] - (volume_acier/V_tot)*avgStrainAcier[0])/(1-(volume_acier/V_tot));
		Av_Strain_mat_yy=( avgStrain[1] - (volume_acier/V_tot)*avgStrainAcier[1])/(1-(volume_acier/V_tot));
		Av_Strain_mat_zz= (avgStrain[2] -  (volume_acier/V_tot)*avgStrainAcier[2])/(1-(volume_acier/V_tot));
		
		nu_acier = (avgStressAcier[1]*avgStrainAcier[0]-avgStressAcier[0]*avgStrainAcier[1])/                                                    	(((avgStressAcier[0]+avgStressAcier[1])*(avgStrainAcier[0]-avgStrainAcier[1]))+((avgStressAcier[0]-avgStressAcier[1])* avgStrainAcier[2] ));
		
		nu_eff =((Av_Stress_mat_yy*Av_Strain_mat_xx)-(Av_Stress_mat_xx*Av_Strain_mat_yy))/                                           (((Av_Stress_mat_xx+Av_Stress_mat_yy)*(Av_Strain_mat_xx-Av_Strain_mat_yy))+((Av_Stress_mat_xx-Av_Stress_mat_yy)*Av_Strain_mat_zz));
		
		E_eff =(Av_Stress_mat_yy*(1+nu_eff)*(1-2*nu_eff))/(nu_eff*Av_Strain_mat_xx+(1-nu_eff)*Av_Strain_mat_yy+nu_eff*Av_Strain_mat_zz);
		
		K_eff_mat = (Av_Stress_mat_xx+Av_Stress_mat_yy+Av_Stress_mat_zz)/(3*(Av_Strain_mat_xx+Av_Strain_mat_yy+Av_Strain_mat_zz));
		
		G_eff_mat =(Av_Stress_mat_xx-((1/3)*(Av_Stress_mat_xx+Av_Stress_mat_yy+Av_Stress_mat_zz)))/                                                (2*((Av_Strain_mat_xx)-((1/3)*(Av_Strain_mat_xx+Av_Strain_mat_yy+Av_Strain_mat_zz))));   
		
		
		
		E_acier = fabs((avgStressAcier[1]*(1+nu_acier)*(1-2*nu_acier))/(nu_acier*avgStrainAcier[0]+                                               (1-nu_acier)*avgStrainAcier[1]+ nu_acier*avgStrainAcier[2]));
		
		K_acier =(avgStressAcier[0]+avgStressAcier[1]+avgStressAcier[2])/(3*(avgStrainAcier[0]+avgStrainAcier[1]+avgStrainAcier[2]));
		
		G_acier = (avgStressAcier[0]-(1/3)*(avgStressAcier[0]+avgStressAcier[1]+avgStressAcier[2]))/(2*((avgStrainAcier[0])-(1/3)*  (avgStrainAcier[0]+avgStrainAcier[1]+avgStrainAcier[2])));
		
		std::cout << std::endl ;
		std::cout << std::endl ;
		
		std::cout << "Modulus (ohne stahl) = "
		          << E_eff
			<< std::endl ;
	   	std::cout << "Poisson's ratio (ohne stahl) = "
		         << nu_eff
			<< std::endl ;
		std::cout << "Bulk Modulus (ohne stahl)="
		           <<  K_eff_mat
			  << std::endl ;
		
		std::cout << "Shear Modulus (ohne stahl)="
		           << G_eff_mat
			  << std::endl ;
		
		std::cout << "Bulk Modulus (ohne stahl) with formulation="
		           << E_eff/(3*(1-2*nu_eff))
			  << std::endl ;
			  
		std::cout << "Shear Modulus (ohne stahl) with formulation="
		           << E_eff/(2*(1+nu_eff))
			  << std::endl ;
		
		  
		std::cout << std::endl ;
		std::cout << std::endl;
			  
		std::cout << "Modulus (Acier)="
		         << E_acier 
			  << std::endl ;
		std::cout << "Poisson's ratio (ohne stahl) = "
		        << nu_acier
			<< std::endl ;
		std::cout << " Bulk Modulus(Acier)="
		          << K_acier
			  << std::endl ;
		std::cout << "Shear Modulus (Acier)="
		           << G_acier
			   << std::endl ;
			
		std::cout << " Bulk Modulus(Acier) with formulation="
		           << E_acier/(3*(1-2*nu_acier))
			  << std::endl ;
		std::cout << "Shear Modulus (Acier) with formulation="
		           <<E_acier/(2*(1+nu_acier))
			   << std::endl ;
			
			
	     //std::cout<< "quit?" << std::flush;
		//std::cin>>quit;
		
		
		std::fstream file ;
		
		file.open("strains.txt", std::ios::out) ;
		for(size_t i = 0 ; i < strains.size() ; i=i+6)
		{
			file << i/6+1 << "\t" << strains[i] << "\t" << strains[i+1] << "\t" << strains[i+2] << "\t" << strains[i+3] << "\t" << strains[i+4] << "\t" << strains[i+5] << std::endl ;
		}
		file.close() ;
		
		file.open("stress.txt", std::ios::out) ;
		for(size_t i = 0 ; i < strains.size() ; i=i+6)
		{
			file << i/6+1 << "\t" << stress[i] << "\t" << stress[i+1] << "\t" << stress[i+2] << "\t" << stress[i+3] << "\t" << stress[i+4] << "\t" << stress[i+5] << std::endl ;
		}
		file.close() ;
		
		file.open("displacements.txt", std::ios::out) ;
		for(size_t i = 0 ; i < x.size() ; i=i+3)
		{
			file << i/3+1 << "\t" << x[i] << "\t" << x[i+1] << "\t" << x[i+2] << std::endl ;
		}
		file.close() ;
		
		   l=l+1;
		
	}
	while( (((sqrt(((E_eff/E_acier)-1)*((E_eff/E_acier)-1))+((nu_eff/nu_acier)-1)*((nu_eff/nu_acier)-1)))>= eps));

	std::cout<< "converged after "<< l <<" étèrations" << std::endl;
	return 0 ;
	
}
