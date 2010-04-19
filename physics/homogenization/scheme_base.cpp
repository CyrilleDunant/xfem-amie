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

#include "scheme_base.h"

using namespace Mu ;




/*Vector Mu::getRawData(const Material & mat, PropertiesType p)
{
	Vector out(0) ;
	size_t i = getFirstIndex(mat,p) ;
	if((i+1) > 0)
		return mat[i].getValues() ;
	return out ;
}

Matrix Mu::getRawVectorData(const std::vector<Material> & mat, PropertiesType p)
{
	Matrix out() ;
	for(size_t i = 0 ; i < mat.size() ; i++)
	{
		Vector iout = getRawData(mat[i],p) ;
		out.push_back(iout) ;
	}
	size_t z = out[0].size() ;
	for(size_t i = 0 ; i < out.size() ; i++)
		z = std::min(z,out[i].size()) ;
	for(size_t i = 0 ; i < out.size() ; i++)
	{
		while(out[i].size() > z)
			out[i].pop_back() ;
	}
	return out ;
}

std::vector<std::vector<double> > Mu::getRawVectorData(const std::vector<Material> & mat, const std::vector<PropertiesType> & p)
{
	std::vector<std::vector<double> > out ;
	for(size_t i = 0 ; i < p.size() ; i++)
	{
		std::vector<std::vector<double> > iout = getRawVectorData(mat,p[i]) ;
		if(i == 0)
			out = iout ;
		else
		{
			for(size_t j = 0 ; j < out.size() ; j++)
			{
				for(size_t k = 0 ; k < iout[j].size() ; k++)
				{
					out[j].push_back(iout[j][k]) ;
				}
			}
		}
	}
	size_t z = out[0].size() ;
	for(size_t i = 0 ; i < out.size() ; i++)
		z = std::min(z,out[i].size()) ;
	for(size_t i = 0 ; i < out.size() ; i++)
	{
		while(out[i].size() > z)
			out[i].pop_back() ;
	}
	return out ;	
}

std::vector<std::vector<double> > Mu::separateByPropertiesType(const std::vector<double> & val, const std::vector<PropertiesType> & p)
{
	std::vector<std::vector<double> > separate ;
	size_t count = 0 ;
	for(size_t i = 0 ; i < p.size() ; i++)
	{
		size_t s = standardNVal(p[i]) ;
		size_t next_count = count ;
		std::vector<double> iseparate ;
		if(s == -1)
			next_count = val.size() ;
		if(s > 0)
			next_count = std::min(val.size(), count+s) ;
		for(size_t j = count ; j < next_count ; j++)
			iseparate.push_back(val[j]) ;
		separate.push_back(iseparate) ;
		count = next_count ;
	}
	return separate ;
}

*/

HomogenizationScheme::HomogenizationScheme(size_t n, PropertiesType in)
{
	nPhases = n ;
	input.push_back(in) ;
	output.push_back(in) ;
}

HomogenizationScheme::HomogenizationScheme(size_t n, PropertiesType in, PropertiesType out)
{
	nPhases = n ;
	input.push_back(in) ;
	output.push_back(out) ;
}

HomogenizationScheme::HomogenizationScheme(size_t n, std::vector<PropertiesType> & in)
{
	nPhases = n ;
	for(size_t i = 0 ; i < in.size() ; i++)
	{
		input.push_back(in[i]) ;
		output.push_back(in[i]) ;
	}
}

HomogenizationScheme::HomogenizationScheme(size_t n, std::vector<PropertiesType> & in, std::vector<PropertiesType> & out)
{
	nPhases = n ;
	for(size_t i = 0 ; i < in.size() ; i++)
		input.push_back(in[i]) ;
	for(size_t j = 0 ; j < out.size() ; j++)
		output.push_back(out[j]) ;
}

bool HomogenizationScheme::verify(const Material & mat)
{
	std::vector<bool> found ;
	for(size_t i = 0 ; i < input.size() ; i++)
		found.push_back(false) ;
	for(size_t i = 0 ; i < mat.size() ; i++)
	{
		for(size_t j = 0 ; j < input.size() ; j++)
		{
			if(mat[i].getPropertiesType() == input[j])
				found[j] = true ;
		}
	}
	for(size_t i = 0 ; i < found.size() ; i++)
	{
		if(!found[i])
			return false ;
	}
	return true ;
}

std::vector<std::vector<size_t> > HomogenizationScheme::getAllPositions(const std::vector<Material> & mat)
{
	std::vector<std::vector<size_t> > index_by_input_by_mat ;
	for(size_t i = 0 ; i < mat.size() ; i++)
	{
		std::vector<size_t> index_by_mat ;
		for(size_t j = 0 ; j < input.size() ; j++)
			index_by_mat.push_back(mat[i].getFirstIndex(input[j])) ;
		index_by_input_by_mat.push_back(index_by_mat) ;
	}
	return index_by_input_by_mat ;
}

std::vector<size_t> HomogenizationScheme::getSchemeSize(const std::vector<Material> & mat)
{
	std::vector<size_t> sz ;
	sz.push_back(nPhases) ;
	if(sz[0]+1 == 0)
		sz[0] = mat.size() ;
	for(size_t i = 0 ; i < input.size() ; i++)
	{
		size_t index = mat[0].getFirstIndex(input[i]) ;
		size_t to_add = mat[0][index].getNVal() ;
		for(size_t j = 0 ; j < mat.size() ; j++)
		{
			index = mat[j].getFirstIndex(input[i]) ;
			to_add = std::min(to_add,mat[j][index].getNVal()) ;
		}
		sz.push_back(to_add) ;
	}
	return sz ;
}

Matrix HomogenizationScheme::getRawData(const std::vector<Material> & mat)
{
	std::vector<size_t> sz = getSchemeSize(mat) ;
	size_t param_sz = 0 ;
	for(size_t i = 1 ; i < sz.size() ; i++)
		param_sz += sz[i] ;

	Matrix data(sz[0],param_sz) ;

	size_t param_i = 0 ;
	for(size_t i = 0 ; i < input.size() ; i++)
	{
		for(size_t k = 0 ; k < mat.size() ; k++)
		{
			size_t index = mat[k].getFirstIndex(input[i]) ;
			for(size_t j = 0 ; j < sz[i+1] ; j++)
			{
				data[k][param_i+j] = mat[k][index].getValue(j) ;
			}
		}
		param_i += sz[i+1] ;
	}

	return data ;
}

std::pair<bool, Material> HomogenizationScheme::apply(const std::vector<Material> & mat)
{
	bool v = true ;
	size_t i = 0 ;
	while(v && i < mat.size())
	{
		v = verify(mat[i]) ;
		i++ ;
	}
	if(!v)
		return std::make_pair(v, mat[0]) ;

	std::vector<Material> mat_reduced ;
	for(size_t i = 0 ; i < mat.size() ; i++)
	{
		size_t combine = -1 ;
		for(size_t j = 0 ; j < mat_reduced.size() ; j++)
		{
			for(size_t k = 0 ; k < output.size() ; k++)
			{
				if(mat[i].equals(mat_reduced[j],output[k]))
					combine = j ;
			}
		}
		if(combine+1 > 0)
			mat_reduced[combine].combine(mat[i],input[0]) ;
		else
			mat_reduced.push_back(mat[i]) ;
	}


	Matrix raw = getRawData(mat_reduced) ;

//	raw.print() ;

/*	for(size_t i = 0 ; i < raw.numRows() ; i++)
		std::cout << raw[i][0] << ";" ;
	std::cout << "   =>    " ; */

	Vector processed = processData(raw) ;

/*	for(size_t i = 0 ; i < processed.size() ; i++)
		std::cout << processed[i] << ";" ;
	std::cout << std::endl ;*/

	Material m_out ;
	size_t count = 0 ;
	size_t next_count = 0 ;
	for(size_t i = 0 ; i < output.size() ; i++)
	{
		if(next_count < processed.size())
		{
			std::vector<double> prop ;
			next_count = count + standardNVal(output[i]) ;
			if(next_count < count)
				next_count = processed.size() ;
	
			for(size_t j = count ; j < next_count ; j++)
				prop.push_back(processed[j]) ;
	
			m_out.push_back(Properties(output[i],prop)) ;
		} else {
			m_out.push_back(Properties()) ;
		}
		count = next_count ;
	}

//	m_out[0].print() ;

	return std::make_pair(v, m_out) ;
}

Vector HomogenizationScheme::processData(const Matrix & data) 
{
	Vector processed(data.numCols()) ;

	for(size_t i = 0 ; i < data.numRows() ; i++)
		processed[i] = data[0][i] ;

	return processed ;
}






MeanSeries::MeanSeries() : HomogenizationScheme(-1,FRACTION,ABSTRACT)
{
	input.push_back(ABSTRACT) ;
}

MeanSeries::MeanSeries(PropertiesType p) : HomogenizationScheme(-1,FRACTION,p)
{
	input.push_back(p) ;
}

Vector MeanSeries::processData(const Matrix & data) 
{
	Vector mean(data.numCols()-1) ;
	if(data.numCols() > 1)
	{
		for(size_t i = 0 ; i < data.numRows() ; i++)
		{
			for(size_t j = 1 ; j < data.numCols() ; j++)
				mean[j-1] += data[i][0] * data[i][j] ;
		}
	}
	return mean ;
}



MeanParallel::MeanParallel() : HomogenizationScheme(-1,FRACTION,ABSTRACT)
{
	input.push_back(ABSTRACT) ;
}

MeanParallel::MeanParallel(PropertiesType p) : HomogenizationScheme(-1,FRACTION,p)
{
	input.push_back(p) ;
}

Vector MeanParallel::processData(const Matrix & data) 
{
	Vector mean(data.numCols()-1) ;
	if(data.numCols() > 1)
	{
		for(size_t i = 0 ; i < data.numRows() ; i++)
		{
			for(size_t j = 1 ; j < data.numCols() ; j++)
				mean[j-1] += data[i][0] / data[i][j] ;
		}
		for(size_t i = 0 ; i < mean.size() ; i++)
			mean[i] = 1 / mean[i] ;
	}
	return mean ;
}




