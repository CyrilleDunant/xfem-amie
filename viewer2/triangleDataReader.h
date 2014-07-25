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
	
	std::vector<quint64> m_positions ;
	quint64 lastKnownPosition ;
	
	bool binary ;
	
public:
	TriangleDataReader(const QString f) : file(f) { 
	
		lastKnownPosition = 0;
		if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
		{
			m_numberOfExtraFields = 0 ;
			m_numberOfTriangles = 0 ;
			m_numberOfPointsPerTriangle = 0;
			return;
		}
		stream.setDevice(&file);
		
		m_positions.clear() ;
		m_positions.push_back((quint64) 0) ;
		prepareFirstHeader() ;

	}
	
	size_t numberOfTriangles() const { return m_numberOfTriangles ;}
	size_t numberOfPointsPerTriangle() const { return m_numberOfPointsPerTriangle ;}
	size_t numberOfExtraFields() const { return m_numberOfExtraFields ;}
	size_t numberOfTimePlanes() const { return m_positions.size() ;}
	
	bool atEnd() const { return stream.atEnd() ; }
	QString getFileName() { return file.fileName() ; }
	void reload()
	{
	    stream.setDevice(&file);
	}
	
	virtual std::vector< std::valarray<float> > * data()
	{
	    std::cerr << file.fileName().toStdString() << std::endl ;
		std::vector< std::valarray<float> > * d = new  std::vector< std::valarray<float> >();
		
		for(size_t j = 0 ; j < m_numberOfPointsPerTriangle ; j++)
		{
			d->push_back(std::valarray<float>(0., m_numberOfTriangles)) ;
			d->push_back(std::valarray<float>(0., m_numberOfTriangles)) ;
		}
		
		for (size_t i = 0 ; i < m_numberOfExtraFields*m_numberOfPointsPerTriangle ; i++)
			d->push_back(std::valarray<float>(0., m_numberOfTriangles)) ;
		
		if(!binary)
		{
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
		}
		else
		{
			double scale = 1./255. ;
			QChar  c ;
			for(size_t i = 0 ; i < m_numberOfTriangles ; i++)
			{
				for(size_t j = 0 ; j < m_numberOfPointsPerTriangle ; j++)
				{
					stream >> (*d)[j*2][i] >> (*d)[j*2+1][i];
/*					stream >> c ;
					(*d)[j*2][i] = scale * c.unicode() ;
					stream >> c ;
					(*d)[j*2+1][i] = scale * c.unicode() ;*/
//					std::cerr << (*d)[j*2][i] << "\t" << (*d)[j*2+1][i] << std::endl ;
				}
			
				for (size_t j = 0 ; j < m_numberOfExtraFields*m_numberOfPointsPerTriangle ; j++)
				{
//					stream >> c ;
//					(*d)[m_numberOfPointsPerTriangle*2+j][i] = scale * c.unicode() ;
					stream >> (*d)[m_numberOfPointsPerTriangle*2+j][i] ;
				}
//				std::cerr << c.unicode() << std::endl ;
			}		
		}
		
		lastKnownPosition += m_numberOfTriangles * (m_numberOfPointsPerTriangle * ( 2 + m_numberOfExtraFields)) ;

		return d ;
	}
	
	virtual std::vector< std::valarray<float> > * dataNext()
	{
		goToPosition(lastKnownPosition) ;
		
		std::vector< std::valarray<float> > * d = new  std::vector< std::valarray<float> >();
		if(!prepareNextHeader())
			return d ;
		
		for(size_t j = 0 ; j < m_numberOfPointsPerTriangle ; j++)
		{
			d->push_back(std::valarray<float>(0., m_numberOfTriangles)) ;
			d->push_back(std::valarray<float>(0., m_numberOfTriangles)) ;
		}
		
		for (size_t i = 0 ; i < m_numberOfExtraFields*m_numberOfPointsPerTriangle ; i++)
			d->push_back(std::valarray<float>(0., m_numberOfTriangles)) ;
		
		if(!binary)
		{
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
		}
		else
		{
			double scale = 1./255. ;
			QChar  c ;
			for(size_t i = 0 ; i < m_numberOfTriangles ; i++)
			{
				for(size_t j = 0 ; j < m_numberOfPointsPerTriangle ; j++)
				{
					stream >> (*d)[j*2][i] >> (*d)[j*2+1][i];
/*					stream >> c ;
					(*d)[j*2][i] = scale * c.unicode() ;
					stream >> c ;
					(*d)[j*2+1][i] = scale * c.unicode() ;*/
//					std::cerr << (*d)[j*2][i] << "\t" << (*d)[j*2+1][i] << std::endl ;
				}
			
				for (size_t j = 0 ; j < m_numberOfExtraFields*m_numberOfPointsPerTriangle ; j++)
				{
//					stream >> c ;
//					(*d)[m_numberOfPointsPerTriangle*2+j][i] = scale * c.unicode() ;
					stream >> (*d)[m_numberOfPointsPerTriangle*2+j][i] ;
				}
//				std::cerr << c.unicode() << std::endl ;
			}		
		}
		
		lastKnownPosition += m_numberOfTriangles * (m_numberOfPointsPerTriangle * ( 2 + m_numberOfExtraFields)) ;

		return d ;
	}
	
	virtual std::vector< std::valarray<float> > * dataAtPos(size_t timePos)
	{
		std::vector< std::valarray<float> > * d = new  std::vector< std::valarray<float> >();
		if(timePos < m_positions.size() && prepareHeaderAt(m_positions[timePos]))
		{
			for(size_t j = 0 ; j < m_numberOfPointsPerTriangle ; j++)
			{
				d->push_back(std::valarray<float>(0., m_numberOfTriangles)) ;
				d->push_back(std::valarray<float>(0., m_numberOfTriangles)) ;
			}
		
			for (size_t i = 0 ; i < m_numberOfExtraFields*m_numberOfPointsPerTriangle ; i++)
				d->push_back(std::valarray<float>(0., m_numberOfTriangles)) ;
		
			if(!binary)
			{
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
			}
			else
			{
				double scale = 1./255. ;
				QChar  c ;
				for(size_t i = 0 ; i < m_numberOfTriangles ; i++)
				{
					for(size_t j = 0 ; j < m_numberOfPointsPerTriangle ; j++)
					{
						stream >> (*d)[j*2][i] >> (*d)[j*2+1][i];
	/*					stream >> c ;
						(*d)[j*2][i] = scale * c.unicode() ;
						stream >> c ;
						(*d)[j*2+1][i] = scale * c.unicode() ;*/
	//					std::cerr << (*d)[j*2][i] << "\t" << (*d)[j*2+1][i] << std::endl ;
					}
			
					for (size_t j = 0 ; j < m_numberOfExtraFields*m_numberOfPointsPerTriangle ; j++)
					{
	//					stream >> c ;
	//					(*d)[m_numberOfPointsPerTriangle*2+j][i] = scale * c.unicode() ;
						stream >> (*d)[m_numberOfPointsPerTriangle*2+j][i] ;
					}
	//				std::cerr << c.unicode() << std::endl ;
				}
			}		
		}
		
		if(m_positions[timePos] + 4 + m_numberOfTriangles * (m_numberOfPointsPerTriangle * ( 2 + m_numberOfExtraFields)) > lastKnownPosition)
			lastKnownPosition = m_positions[timePos] + 4 + m_numberOfTriangles * (m_numberOfPointsPerTriangle * ( 2 + m_numberOfExtraFields)) ;
		
		return d ;
	}
	

	virtual ~TriangleDataReader() { file.close() ; }
	
protected:
	bool addPosition(quint64 pos)
	{
		bool found = false ;
		for(size_t i = 0 ; i < m_positions.size() ; i++)
		{
			found |= (pos == m_positions[i]) ;
		}
		if(!found)
			m_positions.push_back(pos) ;		
		return !found ;
	}
	
	void goToPosition(quint64 position)
	{
		stream.seek(0) ;
		quint64 i = 0 ;
		QString buffer ;
		while(i < position)
		{
			stream >> buffer ;
			i++ ;
		}
	}
	

	virtual bool prepareFirstHeader()
	{
		QString type ;
		stream >> type ;
		if(type != QString("TRIANGLES") && type != QString("BIN_TRIANGLES"))
		{
			m_numberOfExtraFields = 0 ;
			m_numberOfTriangles = 0 ;
			m_numberOfPointsPerTriangle = 0;
			return false;
		}
		binary = (type == QString("BIN_TRIANGLES")) ;
		stream >> m_numberOfTriangles ;
		stream >> m_numberOfPointsPerTriangle ;
		stream >> m_numberOfExtraFields ;
		
		lastKnownPosition = 4 ;
		
		return true ;
	}
	
	virtual bool prepareNextHeader()
	{
		if(stream.atEnd())
		{
			return false ;
		}
		
		QString type ;
		stream >> type ;
		if(type != QString("TRIANGLES") && type != QString("BIN_TRIANGLES"))
		{
			m_numberOfExtraFields = 0 ;
			m_numberOfTriangles = 0 ;
			m_numberOfPointsPerTriangle = 0;
			return false;
		}
		
		binary = (type == QString("BIN_TRIANGLES")) ;
		stream >> m_numberOfTriangles ;
		stream >> m_numberOfPointsPerTriangle ;
		stream >> m_numberOfExtraFields ;

		if(addPosition(lastKnownPosition))
			lastKnownPosition += 4 ;

				
		return true ;
	}

	virtual bool prepareHeaderAt(quint64 pos)
	{
		goToPosition(pos) ;
		
		QString type ;
		stream >> type ;
		if(type != QString("TRIANGLES") && type != QString("BIN_TRIANGLES"))
		{
			m_numberOfExtraFields = 0 ;
			m_numberOfTriangles = 0 ;
			m_numberOfPointsPerTriangle = 0;
			return false;
		}
		
		binary = (type == QString("BIN_TRIANGLES")) ;
		stream >> m_numberOfTriangles ;
		stream >> m_numberOfPointsPerTriangle ;
		stream >> m_numberOfExtraFields ;
		
		if(addPosition(pos))
			lastKnownPosition += 4 ;

		
		return true ;
	}


} ;

#endif
