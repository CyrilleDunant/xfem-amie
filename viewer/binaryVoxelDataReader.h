#ifndef __BINARY_VOXEL_DATA_READER_H__
#define __BINARY_VOXEL_DATA_READER_H__

#include "voxelDataReader.h"

#include <QEventLoop>

class BinaryVoxelDataReader : public VoxelDataReader
{
	
public:
	BinaryVoxelDataReader(const QString file) : VoxelDataReader(file) {
	} ;
	
	virtual std::vector< std::valarray<quint8> > * data() 
	{
		file.open(QIODevice::ReadOnly) ;
// 		if(file.isSequential())
// 			std::cout << "ouch" << std::endl ;
// 		else
// 			std::cout << position << std::endl ;
// 		file.seek(position) ;
		streamBin.setDevice(&file);
		streamBin.skipRawData(21) ;
		
		std::vector< std::valarray<quint8> > * d = new  std::vector< std::valarray<quint8> >(numberOfFields,
			std::valarray<quint8>((quint8)0, r*c*s+3)) ;

		
		std::cout << "fields = " << numberOfFields << std::endl ;
		std::cout << "datapoints = " << r*c*s << std::endl ;
		
		QEventLoop eloop ;
		for(size_t i = 0 ; i < r*c*s ; i++)
		{
			if(i % (r*c*s/100) == 0)
			{
				int percent = (int)((float)i/(float)(r*c*s)*100.f) ;
				emit progressed(percent) ;
				eloop.wakeUp() ;
				eloop.processEvents() ;
			}
			
			for(size_t j = 0 ; j < numberOfFields ; j++)
			{
				streamBin >> (quint8&)((*d)[j])[i] ;
// 				std::cout << i << " -> "<< (unsigned int)((*d)[j][i]) << std::endl ;
			}
		}
		emit progressed(100) ;
		eloop.processEvents() ;
		eloop.exit() ;
		return d ;
	}
	
	virtual ~BinaryVoxelDataReader() { };
} ;

#endif
