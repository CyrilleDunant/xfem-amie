#ifndef __ASCII_VOXEL_DATA_READER_H__
#define __ASCII_VOXEL_DATA_READER_H__

#include "voxelDataReader.h"

class AsciiVoxelDataReader : public VoxelDataReader
{
public:
    AsciiVoxelDataReader(const QString file) : VoxelDataReader(file) {
    } ;

    virtual std::vector< std::valarray<quint8> > * data()
    {
        file.open(QIODevice::ReadOnly | QIODevice::Text) ;
        file.seek(position) ;
        streamAsc.setDevice(&file);

        std::vector< std::valarray<quint8> > * d = new  std::vector< std::valarray<quint8> >(numberOfFields,
                std::valarray<quint8>((quint8)0, r*c*s+3)) ;


        std::cout << "fields = " << numberOfFields << std::endl ;
        std::cout << "datapoints = " << r*c*s << std::endl ;
        for(size_t i = 0 ; i < r*c*s ; i++)
        {
            if(i%10000 == 0)
                std::cout << "\r reading ... (%) " << (int)((float)i/(float)(r*c*s)*100.f) << std::flush ;

            for(size_t j = 0 ; j < numberOfFields ; j++)
            {
                streamAsc >> (int&)((*d)[j])[i] ;
            }
        }
        std::cout << "\r reading ... (%) 100        " << std::endl ;
        return d ;
    }

    virtual ~AsciiVoxelDataReader() { };
} ;

#endif
