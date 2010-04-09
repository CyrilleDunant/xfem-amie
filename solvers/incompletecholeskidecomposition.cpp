
//
// C++ Implementation: incompletecholeskidecomposition
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "incompletecholeskidecomposition.h"
#include "choleskidecomposed.h"

using namespace Mu ;

// InCompleteCholesky::InCompleteCholesky(const CoordinateIndexedSparseMatrix &A_) : A(A_)
// {
// 	stable = true ;
// 	
// 	Vector diag = A_.diagonal() ;
// 	double fac =  1 ; 
// 	
// 	std::cout << "froebenius Norm = " << A_.froebeniusNorm() << ", infinity Norm = "<< A_.infinityNorm()<< std::endl ;
// 	
// 	size_t tries = 0 ;
// 	for(int i = 0 ; i < (int)A.row_size.size() && (int)tries < 16; i++)
// 	{
// 		if(i%1000 == 0)
// 			std::cout << "\r try " << tries <<  ", computing sparse Choleski line " << i << "/" << A.row_size.size() << std::flush ;
// 		for(size_t j = 0 ; A.column_index[A.accumulated_row_size[i]+j] < (size_t)i ; j++)
// 		{
// 			size_t a_i = A.accumulated_row_size[i]+j ;
// 			size_t col = A.column_index[a_i] ;
// 			
// 			if(std::abs((int)col-i) < 24)
// // 			if(std::abs(A[col][col]) > 1e-6)
// 			{
// 				A.array[a_i] = (A_.array[a_i] - innerProduct( A[i], A[col], col))/(A[col][col]*fac) ;
// // 			else
// // 				A.array[a_i] = 0 ;
// 				A[col][i] = A.array[a_i] ;
// 			}
// 		}	
// 		
// 		double sum = innerProduct( A[i], A[i], i) ; 
// 		double val = diag[i] - sum;
// 		if (val <= 0)
// 		{
// 			std::cout << val << " = " << diag[i] <<" - "<< sum << std::endl ;
// // 			if(diag[i] < 0)
// // 			{
// // 				val = 2 ;
// // 				A[i][i] = sqrt(val);
// // // 				std::cout << "\n oops, diag = "<< diag[i] << ", sum = "<< sum  << std::endl ;
// // 			}
// // 			else
// // 			{
// 			fac+=0.1 ;
// 			tries++ ;
// 			i = -1 ;
// // 			}
// 		}
// 		else
// 			A[i][i] = sqrt(val);
// 		
// 	}
// 	
// 	if (tries == 4)
// 	{
// 		stable = false ;
// 		std::valarray<unsigned int> rs(A_.row_size.size()) ;
// 		std::valarray<unsigned int> ci(A_.row_size.size()) ;
// 		
// 		rs = 1 ;
// 		
// 		for(size_t i = 0 ; i < rs.size() ; i++)
// 			ci[i] = i ;
// 		
// 		A = CoordinateIndexedSparseMatrix(rs, ci, 1) ;
// 		
// 		for(size_t i = 0 ; i < rs.size() ; i++)
// 		{
// 			A.array[i] = diag[i] ;
// 		}
// 	}
// 	
// 	std::cout << " ...done" << std::endl ;
// 	d.resize(diag.size()) ;
// 	for(size_t i = 0 ; i < diag.size() ; i++)
// 	{
// 		d = 1./diag[i] ;
// 	}
// }
// 
// void InCompleteCholesky::precondition(const Vector &v, Vector &t) const
// {
// // 	if(stable)
// // 	{
// //		Vector y = UpperTriangular(A, v).solve(Vector(0), NULL) ;
// //		Vector ret= LowerTriangular(A, y).solve(Vector(0), NULL) ;
// 	
// 	t = CholeskiDecomposed(A,v,d).solve(Vector(0), NULL) ;
// // 	}
// // 	else
// // 		return UpperTriangular(A, v).solve(Vector(0), NULL) ;
// }
