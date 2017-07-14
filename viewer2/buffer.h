// #define nullptr NULL ;

#ifndef BUFFER_H
#define BUFFER_H

#include "triangleDataReader.h"

struct LayerBuffer
{
	public:
		int numberOfPointsPerTriangle ;
		std::vector<std::valarray<float > > * values ;

		LayerBuffer( int n, std::vector<std::valarray<float> > * v )
		{
			numberOfPointsPerTriangle = n ;
			values = v ;
		}

		LayerBuffer()
		{
			numberOfPointsPerTriangle = 0 ;
			values = nullptr;
		}

};

class FileBuffer : public QList<LayerBuffer>
{
	public:
		QString file ;
		QString short_name ;

		FileBuffer() : QList<LayerBuffer>() { }

		FileBuffer( QString n ) : QList<LayerBuffer>()
		{
			file = n ;
			short_name = n.section( '/', -1 ) ;
			TriangleDataReader reader( n ) ;

			int np = reader.numberOfPointsPerTriangle() ;
			push_back( LayerBuffer( np, reader.data() ) );
			std::vector<std::valarray<float> > * d = reader.dataNext() ;

			while( d->size() > 0 )
			{
				push_back( LayerBuffer( np, d ) );
				d = reader.dataNext() ;
			}
		}

		bool is( QString n )
		{
			return file.compare( n ) == 0 || short_name.compare( n ) == 0 ;
		}
};

class Buffer : public QList<FileBuffer>
{
	QString master ;
	std::vector<QString> files ;
	public:
		Buffer() : QList<FileBuffer>() {  }

		Buffer( const QString & n ) : QList<FileBuffer>()
		{
			master = n ;
			QFile file( n ) ;

			if( !file.open( QIODevice::ReadOnly | QIODevice::Text ) )
				return ;

			QTextStream stream(&file);
			
			QString type ;
			stream >> type ;
			

			if( type != "MULTI" )
				return ;

			bool go_on = true ;
			QString path = n.section( '/', 0, -2 ) + "/" ;

			while( go_on )
			{
				stream >> type ;

				if( type.size() > 0 )
                {
					files.push_back( path + type ) ;
                    push_back(FileBuffer(files.back())) ;
                }
				else
					go_on = false ;
			}
			file.close();
			
		}
		
		bool update()
		{
			QFile file( master ) ;

			if( !file.open( QIODevice::ReadOnly | QIODevice::Text ) )
				return false;

			QTextStream stream;
			stream.setDevice( &file );
			QString type ;
			stream >> type ;
			file.close();

			if( type != "MULTI" )
				return false;
			int originalSize = files.size() ;
			files.clear() ;
			bool go_on = true ;
			QString path = master.section( '/', 0, -2 ) + "/" ;
			int newSize = 0 ;
			while( go_on )
			{
				newSize++ ;
				stream >> type ;

				if( type.size() > 0 )
					files.push_back( path + type ) ;
				else
					go_on = false ;
			}
			
			for(int i = originalSize ; i < newSize-1 ; i++)
				push_back(FileBuffer(files[i])) ;
			
			return originalSize != newSize ;
		}

		int fileInBuffer( QString file )
		{
			for( int i = 0 ; i < size() ; i++ )
			{
				if( value( i ).is( file ) )
					return i ;
			}

			return -1 ;
		}

};

#endif // BUFFER_H
