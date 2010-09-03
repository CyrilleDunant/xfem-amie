#ifndef __TRIANGLE_DATA_READER_H__
#define __TRIANGLE_DATA_READER_H__

#include <QString>
#include <QFile>
#include <QTextStream>

class TriangleDataReader
{
protected:
	QFile file ;
	QTextStream stream;
	
	quint64 m_numberOfTriangles ;
	quint64 m_numberOfPointsPerTriangle ;
	quint64 m_numberOfExtraFields ;
	
public:
	TriangleDataReader(const QString f) : file(f) { 
	
		if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
		{
			m_numberOfExtraFields = 0 ;
			m_numberOfTriangles = 0 ;
			m_numberOfPointsPerTriangle = 0;
			return;
		}
		stream.setDevice(&file);
		QString type ;
		stream >> type ;
		if(type != QString("TRIANGLES"))
		{
			m_numberOfExtraFields = 0 ;
			m_numberOfTriangles = 0 ;
			m_numberOfPointsPerTriangle = 0;
			return;
		}
		stream >> m_numberOfTriangles ;
		stream >> m_numberOfPointsPerTriangle ;
		stream >> m_numberOfExtraFields ;

	}
	
	size_t numberOfTriangles() const { return m_numberOfTriangles ;}
	size_t numberOfPointsPerTriangle() const { return m_numberOfPointsPerTriangle ;}
	size_t numberOfExtraFields() const { return m_numberOfExtraFields ;}
	
	virtual std::vector< std::valarray<float> > * data()
	{
		std::vector< std::valarray<float> > * d = new  std::vector< std::valarray<float> >();
		
		for(size_t j = 0 ; j < m_numberOfPointsPerTriangle ; j++)
		{
			d->push_back(std::valarray<float>(0., m_numberOfTriangles)) ;
			d->push_back(std::valarray<float>(0., m_numberOfTriangles)) ;
		}
		
		for (size_t i = 0 ; i < m_numberOfExtraFields*m_numberOfPointsPerTriangle ; i++)
			d->push_back(std::valarray<float>(0., m_numberOfTriangles)) ;
		
		for(size_t i = 0 ; i < m_numberOfTriangles ; i++)
		{
			for(size_t j = 0 ; j < m_numberOfPointsPerTriangle ; j++)
			{
				stream >> (*d)[j*2][i] >> (*d)[j*2+1][i];
			}
			
			for (size_t j = 0 ; j < m_numberOfExtraFields*m_numberOfPointsPerTriangle ; j++)
			{
				stream >> (*d)[m_numberOfPointsPerTriangle*2+j][i] ;
			}
		}
		
		return d ;
	}
	
	virtual ~TriangleDataReader() { file.close() ; }
} ;

#endif
