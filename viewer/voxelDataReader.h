#ifndef __VOXEL_DATA_READER_H__
#define __VOXEL_DATA_READER_H__

#include <QString>
#include <QFile>
#include <QDataStream>
#include <QTextStream>

#include <vector>
#include <valarray>
#include <iostream>

class VoxelDataReader :public QObject
{
Q_OBJECT      // must include this if you use Qt signals/slots
		
protected:
	QFile file ;
	QDataStream streamBin;
	QTextStream streamAsc;
	quint64 numberOfFields ;
	quint64 r ;
	quint64 c ;
	quint64 s ;
	
	qint64 position ;
	
public:
	
	VoxelDataReader(const QString f) : file(f) { 
	
		if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
		{
			r = 0 ;
			c = 0 ;
			s = 0 ;
			numberOfFields = 0 ;
			return;
		}

		streamAsc.setDevice(&file) ;
		
		QString type ;
		streamAsc >> type ;
		if(type != "VOXELS")
		{
			r = 0 ;
			c = 0 ;
			s = 0 ;
			numberOfFields = 0 ;
			return;
		}
		streamAsc >> numberOfFields ;
		streamAsc >> r ;
		streamAsc >> c ;
		streamAsc >> s ;
		
		std::cout << "i = " << r << ", j = " << c << ", k =" << s << std::endl ;
// 		file.setOpenMode(QIODevice::ReadOnly) ;
		position = file.pos() ; //save position
		file.close() ;
	}


	size_t rows() const { return r ;}
	size_t columns() const { return c ;}
	size_t strips() const { return s ;}
	
	virtual std::vector< std::valarray<quint8> > * data()  =  0 ;
	
	virtual ~VoxelDataReader() { file.close() ; }
	
signals:
	void progressed(int) ;
} ;

#endif
