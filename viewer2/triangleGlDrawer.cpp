#include "triangleGlDrawer.h"

void TriangleGLDrawer::computeDisplay( ) const
{

	glMatrixMode( GL_MODELVIEW );

	glColor4f( 1, 1, 1, 0.5 ) ;
	glBegin( GL_LINE_LOOP ) ;
	glVertex2f( 0.50 * .8 - .15, 0.50 * .8 + .05) ;
	glVertex2f( -0.50 * .8 - .15, 0.50 * .8 + .05) ;
	glVertex2f( -0.50 * .8 - .15, -0.50 * .8 + .05) ;
	glVertex2f( 0.50 * .8 - .15, -0.50 * .8 + .05) ;
	glEnd() ;
	glCallList( currentDisplayList ) ;
	glCallList( wireFrameList ) ;
}


void TriangleGLDrawer::mouseReleaseEvent( QMouseEvent *event )
{
	moving = false ;
	leftDown = false ;

	if( event->button() == Qt::LeftButton )
	{
		QImage fb = this->grabFrameBuffer() ;
		QRgb color = fb.pixel(event->x(), event->y());
		
		
// 		double prop = (((float)QColor(color).lightness()-15.56339)*1.065)/255. ;
// 		valUnderCursor = prop * ( limits[currentSet].first - ( limits[currentSet].first - limits[currentSet].second ) * ( 1. - fracup / 10000. ) ) + ( 1. - prop ) * ( limits[currentSet].second + ( limits[currentSet].first - limits[currentSet].second ) * fracdown / 10000. ) ;
		
		double downcut =  ( limits[currentSet].first - limits[currentSet].second ) * ( fracdown / 10000. ) ;
		double upcut =  ( limits[currentSet].second - limits[currentSet].first ) * ( fracup / 10000.-1. ) ;
		valUnderCursor = 1.049*limits[currentSet].first + downcut + (float)QColor(color).lightness()/255.*( limits[currentSet].second- limits[currentSet].first ) - ((float)QColor(color).lightness()/255.)*0.049 - downcut;
		mousePosOnLeftClick = QPoint( 0, 0 ) ;
	}

	paintGL() ;
}

void TriangleGLDrawer::mouseMoveEvent( QMouseEvent *event )
{
	if( leftDown )
	{
		xtransleft -= ( int )round( ( float )( mousePosOnLeftClick.x() - event->x() ) * ( 100000. / ( float )zoom ) ) ;
		ytransleft += ( int )round( ( float )( mousePosOnLeftClick.y() - event->y() ) * ( 100000. / ( float )zoom ) ) ;
		emit xTranslationChanged( xtransleft ) ;
		emit yTranslationChanged( ytransleft ) ;
		mousePosOnLeftClick = event->pos() ;
	}

	paintGL() ;
}

void TriangleGLDrawer::mousePressEvent( QMouseEvent *event )
{
// 	grabKeyboard() ;

	moving = true ;

	if( event->button() == Qt::LeftButton )
	{
		leftDown = true ;
		mousePosOnLeftClick = event->pos() ;
	}
}

void TriangleGLDrawer::keyPressEvent( QKeyEvent *event )
{
	if( event->key() == Qt::Key_Plus )
	{
		zoom = ( int )round( zoom * 1.1 ) ;
		setZoom( zoom ) ;
	}
	else if( event->key() == Qt::Key_Minus )
	{
		zoom = ( int )round( zoom / 1.1 ) ;
		setZoom( zoom ) ;
	}
	else if( event->key() == Qt::Key_D )
	{
		currentDisplayList++ ;
		currentDisplayList -= displayList ;
		currentDisplayList %= numberOfExtraFields ;
		currentDisplayList += displayList ;
	}
	else
	{
		event->ignore() ;
	}

	paintGL() ;
}

void TriangleGLDrawer::wheelEvent( QWheelEvent *event )
{

	int newZoom = zoom ;

	if( event->delta() > 0 )
		newZoom = ( float )zoom * 1.1 ;
	else
		newZoom = ( float )zoom * 0.9 ;

	zoom = newZoom ;
	setZoom( newZoom );
}

void TriangleGLDrawer::setSet( int set )
{
// 	set -= displayList ;
	if( set < 0 )
		set = numberOfExtraFields - 1 ;

	if( set == numberOfExtraFields )
		set = 0 ;

	currentDisplayList = set + displayList ;
	currentSet = set ;
	paintGL() ;

	emit setChanged( set ) ;
}

void TriangleGLDrawer::setXTranslation( int trans )
{
	xtransleft = trans ;
	emit xTranslationChanged( xtransleft ) ;
}

void TriangleGLDrawer::setYTranslation( int trans )
{
	ytransleft = trans ;
	emit yTranslationChanged( ytransleft ) ;
}

void TriangleGLDrawer::paintGL()
{

	QTime startTime = QTime::currentTime();
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT ); // Clear Screen And Depth Buffer

	glMatrixMode( GL_MODELVIEW );                           // Select The Modelview Matrix

	glLoadIdentity();
	glLookAt( getXtrans(), getYtrans(), zpos, getXtrans(), getYtrans(), 0, 0, 1, 0 ) ;
	computeDisplay() ;

	glColor3f( 0., 0., 0. ) ;
	int elapsedTime = startTime.msecsTo( QTime::currentTime() );
	renderText( 10, 20, QString( "%0 fps" ).arg( 1000.0f / ( float )elapsedTime ) );
	size_t r, g, b ;

	renderText( 10, 40, QString( "%0 []" ).arg( ( float )valUnderCursor) );
	
	
	if( !limits.empty() )
	{
		double rel =  (1.-(valUnderCursor-limits[currentSet].first)/(limits[currentSet].second - limits[currentSet].first)) ;
		glBegin( GL_QUAD_STRIP ) ;

// 		for(double i = 1. ; i > fracup/10000. ; i -= 0.001)
// 		{
// 			HSVtoRGB(&r, &g, &b, 180., 0., 0.1) ;
// 			glColor4ub(r, g, b, 255) ;
// 			glVertex2f((.805-0.5) , (i-0.5)*.7 ) ;
// 			glVertex2f((.835-0.5) , (i-0.5)*.7 ) ;
// 		}
		for( double i = ( double )fracup / 10000. ; i > ( double )fracdown / 10000. ; i -= 0.001 )
		{
			double where = ( 1. - ( ( double )fracup / 10000. - i ) / ( ( double )fracup / 10000. - ( double )fracdown / 10000. ) ) ;
			HSVtoRGB( &r, &g, &b, 180., 0., ( 1. - where ) ) ;
			if(std::abs(rel-where) < .5e-2 )
			{
				r = 255 ;
				g = 0 ;
				b = 0 ;
			}
			
			glColor4ub( r, g, b, 255 ) ;
			glVertex2f( ( .805 - 0.5 ) , ( where - 0.5 )*.7 ) ;
			glVertex2f( ( .835 - 0.5 ) , ( where - 0.5 )*.7 ) ;
		}

// 		for(double i = fracdown/10000. ; i > 0 ; i -= 0.001)
// 		{
// 			HSVtoRGB(&r, &g, &b, 180., 0., 0.9) ;
// 			glColor4ub(r, g, b, 255) ;
// 			glVertex2f((.805-0.5) , (i-0.5)*.7 ) ;
// 			glVertex2f((.835-0.5) , (i-0.5)*.7 ) ;
// 		}
		glEnd() ;
		glColor4ub( 0, 0, 0, 255 ) ;

		for( double i = 1. ; i > 0 ; i -= 0.1 )
		{

			double v = i * ( limits[currentSet].first - ( limits[currentSet].first - limits[currentSet].second ) * ( 1. - fracup / 10000. ) ) + ( 1. - i ) * ( limits[currentSet].second + ( limits[currentSet].first - limits[currentSet].second ) * fracdown / 10000. ) ;
			renderText( ( .805 - 0.5 ) + 0.05,
			            ( i - 0.5 )*.7,
			            0.,
			            QString::number( v, 'f', 2) );
		}
	}

	glFinish();
	swapBuffers() ;
}

void TriangleGLDrawer::resizeGL( int w, int h )
{

	if( h == 0 )                                     // Prevent A Divide By Zero By
		h = 1;                                         // Making Height Equal One

	if( h > w )
		glViewport( 0, 0, h, h );
	else
		glViewport( 0, 0, w, w );



// 	glPointSize(std::max(0.5*(float)width()/(float)columns, 1.));

	setZoom( zoom ) ;
}

void TriangleGLDrawer::setZoom( int percent )
{
	size_t p ;

	if( percent == 0 )
		p = 1 ;
	else
		p = percent ;

	float fact = 45.*( ( float )100 / ( float )p ) ;
	glMatrixMode( GL_PROJECTION );
	glLoadIdentity();
	glPerspective( fact, 1., 1.f, 30.0f );
	glMatrixMode( GL_MODELVIEW );
	glLoadIdentity();

	paintGL() ;

	emit zoomChanged( p ) ;
}

void TriangleGLDrawer::openFile( const QString f )
{

	delete valuesAtPoint ;
	delete reader ;

	reader = new TriangleDataReader( f ) ;

	numberOfTriangles = reader->numberOfTriangles() ;
	numberOfPointsPerTriangle = reader->numberOfPointsPerTriangle() ;
	numberOfExtraFields = reader->numberOfExtraFields() ;
	numberOfExtraTimePlanes = reader->numberOfTimePlanes() - 1 ;
	currentTimePlane = 0 ;

	valuesAtPoint = reader->data() ;

}

void TriangleGLDrawer::openFile( TriangleDataReader *f )
{

	delete valuesAtPoint ;
	delete reader ;

	reader = f ;

	numberOfTriangles = reader->numberOfTriangles() ;
	numberOfPointsPerTriangle = reader->numberOfPointsPerTriangle() ;
	numberOfExtraFields = reader->numberOfExtraFields() ;
	numberOfExtraTimePlanes = reader->numberOfTimePlanes() - 1 ;
	currentTimePlane = 0 ;

	valuesAtPoint = reader->data() ;

}


void TriangleGLDrawer::initializeGL()
{


	displayList = glGenLists( numberOfExtraFields+1 ) ;
	currentDisplayList = currentSet + 1 ;
	glViewport( 0, 0, 600, 600 ) ;
//	glEnable( GL_BLEND );
	glShadeModel( GL_SMOOTH ); // Enables Smooth Shading
	glEnable( GL_LINE_SMOOTH ) ;
	glEnable( GL_POLYGON_SMOOTH );
	glDisable( GL_DEPTH_TEST );
	glHint( GL_POLYGON_SMOOTH_HINT, GL_NICEST ) ;

// 	glPointSize(std::max(0.4*(float)width()/(float)columns, 1.));
	glClearColor( 1.0f, 1.0f, 1.0f, 0.0f );                                 // Black Background
	glClearDepth( 1.0f );                                                   // Depth Buffer Setup
	glDisable( GL_DEPTH_TEST );
//	glBlendFunc( GL_SRC_ALPHA, GL_ONE );

	computeDisplayList() ;
}

QSize TriangleGLDrawer::minimumSizeHint() const
{
	return QSize( 50, 50 );
}

QSize TriangleGLDrawer::sizeHint() const
{
	return QSize( 600, 600 );
}

void TriangleGLDrawer::grab()
{

	paintGL() ;

	QImage grab = grabFrameBuffer() ;
	QString s( fileName ) ;

	s += QString( ".png" ) ;

	grab.save( s, "png" ) ;

}

void TriangleGLDrawer::computeDisplayList()
{

	if( !minmaxinit )
	{
		max_x = ( *valuesAtPoint )[0][0] ;
		min_x = ( *valuesAtPoint )[0][0] ;
		max_y = ( *valuesAtPoint )[1][0] ;
		min_y = ( *valuesAtPoint )[1][0] ;
		minmaxinit = true ;
	}

	for( size_t i = 0 ; i < numberOfTriangles ; i++ )
	{
		for( size_t j = 0 ; j < numberOfPointsPerTriangle ; j++ )
		{
			if( ( *valuesAtPoint )[2 * j][i] > max_x )
				max_x = ( *valuesAtPoint )[2 * j][i] ;

			if( ( *valuesAtPoint )[2 * j + 1][i] > max_y )
				max_y = ( *valuesAtPoint )[2 * j + 1][i] ;

			if( ( *valuesAtPoint )[2 * j][i] < min_x )
				min_x = ( *valuesAtPoint )[2 * j][i];

			if( ( *valuesAtPoint )[2 * j + 1][i] < min_y )
				min_y = ( *valuesAtPoint )[2 * j + 1][i] ;
		}
	}

	for( size_t N = 0 ; N < numberOfExtraFields ; N++ )
	{
		glNewList( displayList + N, GL_COMPILE ) ;

		float max_val = ( *valuesAtPoint )[( 2 + N ) * numberOfPointsPerTriangle][0] ;
		float min_val = ( *valuesAtPoint )[( 2 + N ) * numberOfPointsPerTriangle][0] ;
		std::vector<float> vals ;

		for( size_t i = 0 ; i < numberOfTriangles ; i++ )
		{
			for( size_t j = 0 ; j < numberOfPointsPerTriangle ; j++ )
			{
				vals.push_back( ( *valuesAtPoint )[( 2 + N )*numberOfPointsPerTriangle + j][i] ) ;

				if( vals.back() > max_val )
					max_val = vals.back() ;

				if( vals.back() < min_val )
					min_val = vals.back();

			}
		}

		if( limits.size() <= N )
		{
			limits.push_back( std::make_pair( max_val, min_val ) ) ;
		}

		if( limits[N].first < max_val )
			limits[N].first = max_val ;

		if( limits[N].second > min_val )
			limits[N].second = min_val ;

		if( std::abs( limits[N].first + limits[N].second ) / ( limits[N].first - limits[N].second ) < 1e-6 )
		{
			limits[N].first = 1e-6 ;
			limits[N].second = -1e-6 ;
		}

		max_val = limits[N].first ;
		min_val = limits[N].second ;
		if(max_val > 0 && min_val > 0)
		{
			min_val = 0 ;
			limits[N].second = 0;
		}
		if(max_val < 0 && min_val < 0)
		{
			max_val = 0 ;
			limits[N].first = 0;
		}
		
		double maxdelta = std::max( max_x - min_x, max_y - min_y ) / 0.8 ;
		double mindelta = std::min( max_x - min_x, max_y - min_y ) / 0.8 ;
		double cx = 0.7 ;
		double cy = 0.7 ;
		double mag = scale ;

		if( max_x - min_x > max_y - min_y )
			cy = 0.7 * mindelta / maxdelta ;
		else
			cx = 0.7 * mindelta / maxdelta ;

		size_t r, g, b ;

		for( size_t i = 0 ; i < numberOfTriangles ; i++ )
		{
			glBegin( GL_TRIANGLES ) ;

			for( size_t j = 0 ; j < numberOfPointsPerTriangle ; j++ )
			{
				float v = ( ( *valuesAtPoint )[( 2 + N ) * numberOfPointsPerTriangle + j][i] - min_val ) / ( max_val - min_val );
				if( std::abs(min_val) > 1e-4  &&  abs(max_val / min_val - 1) < 1e-6)
					v = 0.5 ;
				
				if( v < ( double )fracdown / 10000. )
					v = 0 ;
				else if( v > ( double )fracup / 10000. )
					v = 1 ;
				else
					v = ( v - ( double )fracdown / 10000. ) / ( ( ( double )fracup - ( double )fracdown ) / 10000. ) ;

				HSVtoRGB( &r, &g, &b, 180., 0., 1.-v ) ;
				glColor4ub( std::min(r,(size_t)240), std::min(g,(size_t)240), std::min(b,(size_t)240), 255 ) ;

				double dx = ( *valuesAtPoint )[( 2 ) * numberOfPointsPerTriangle + j][i] * mag/std::max((max_x-min_x), (max_y-min_y)) ;
				double dy = ( *valuesAtPoint )[( 3 ) * numberOfPointsPerTriangle + j][i] * mag/std::max((max_x-min_x), (max_y-min_y)) ;
				glVertex2f( ( ( *valuesAtPoint )[j * 2][i] - min_x ) / maxdelta - 0.5 * cx + dx - 0.2, ( ( *valuesAtPoint )[j * 2 + 1][i] - min_y ) / maxdelta - 0.5 * cy + dy ) ;
			}

			glEnd() ;
		}

		
// 		glLineWidth(1) ;
// 		for( size_t i = 0 ; i < numberOfTriangles ; i++ )
// 		{
// 			glBegin( GL_LINE_LOOP ) ;
// 
// 			for( size_t j = 0 ; j < numberOfPointsPerTriangle ; j++ )
// 			{
// 
// 	//			HSVtoRGB( &r, &g, &b, 180., 0., 0. ) ;
// 
// 				glColor4ub( 0, 0, 0, 0 ) ;
// 
// 				double dx = ( *valuesAtPoint )[( 2 ) * numberOfPointsPerTriangle + j][i] * mag/std::max((max_x-min_x), (max_y-min_y)) ;
// 				double dy = ( *valuesAtPoint )[( 3 ) * numberOfPointsPerTriangle + j][i] * mag/std::max((max_x-min_x), (max_y-min_y)) ;
// 				glVertex2f( ( ( *valuesAtPoint )[j * 2][i] - min_x ) / maxdelta - 0.5 * cx + dx - 0.2, ( ( *valuesAtPoint )[j * 2 + 1][i] - min_y ) / maxdelta - 0.5 * cy + dy ) ;
// 			}
// 
// 			glEnd() ;
// 		}
// 			
		glEndList() ;
	}
	
	glNewList( displayList+numberOfExtraFields, GL_COMPILE ) ;
	wireFrameList = displayList+numberOfExtraFields ;
	double maxdelta = std::max( max_x - min_x, max_y - min_y ) / 0.8 ;
	double mindelta = std::min( max_x - min_x, max_y - min_y ) / 0.8 ;
	double cx = 0.7 ;
	double cy = 0.7 ;
	double mag = scale ;

	if( max_x - min_x > max_y - min_y )
		cy = 0.7 * mindelta / maxdelta ;
	else
		cx = 0.7 * mindelta / maxdelta ;

	size_t r, g, b ;

	glLineWidth(1) ;
/*	for( size_t i = 0 ; i < numberOfTriangles ; i++ )
	{
		glBegin( GL_LINE_LOOP ) ;

		for( size_t j = 0 ; j < numberOfPointsPerTriangle ; j++ )
		{

//			HSVtoRGB( &r, &g, &b, 180., 0., 0. ) ;

			glColor4ub( 0, 0, 0, 0 ) ;

			double dx = ( *valuesAtPoint )[( 2 ) * numberOfPointsPerTriangle + j][i] * mag/std::max((max_x-min_x), (max_y-min_y)) ;
			double dy = ( *valuesAtPoint )[( 3 ) * numberOfPointsPerTriangle + j][i] * mag/std::max((max_x-min_x), (max_y-min_y)) ;
			glVertex2f( ( ( *valuesAtPoint )[j * 2][i] - min_x ) / maxdelta - 0.5 * cx + dx - 0.2, ( ( *valuesAtPoint )[j * 2 + 1][i] - min_y ) / maxdelta - 0.5 * cy + dy ) ;
		}

		glEnd() ;
	}*/

	glEndList() ;
}


void TriangleGLDrawer::setSegmentDown( int v )
{
	if( v < 0 )
		v = 0 ;

	if( v > fracup )
		v = fracup - 1 ;

	fracdown = v ;
	computeDisplayList() ;
	paintGL() ;
	emit segmentDownChanged( v ) ;

}

void TriangleGLDrawer::setScale( int v )
{
	if( v < 0 )
		v = 0 ;

	scale = v ;
	computeDisplayList() ;
	paintGL() ;
	emit scaleChanged( v ) ;
}

void TriangleGLDrawer::setTimePlane( int tp )
{
	if(reader)
	{
		delete valuesAtPoint ;

		if( tp <= numberOfExtraTimePlanes )
		{
			valuesAtPoint = reader->dataAtPos( tp ) ;
		}
		else
		{
			valuesAtPoint = reader->dataNext() ;

			if( valuesAtPoint->size() == 0 )
			{
				delete valuesAtPoint ;
				tp = 0 ;
				valuesAtPoint = reader->dataAtPos( tp ) ;
			}
			else
			{
				numberOfExtraTimePlanes = tp ;
			}
		}

		currentTimePlane = tp ;
		numberOfTriangles = reader->numberOfTriangles() ;
		numberOfPointsPerTriangle = reader->numberOfPointsPerTriangle() ;
		numberOfExtraFields = reader->numberOfExtraFields() ;

		computeDisplayList() ;
		paintGL() ;
	}

//	emit timePlaneChanged(tp) ;
}


void TriangleGLDrawer::setSegmentUp( int v )
{
	if( v < fracdown )
		v = fracdown + 1 ;

	if( v > 10000 )
		v = 10000 ;

	if( fracup != v )
	{
		fracup = v ;

		computeDisplayList() ;
		paintGL() ;
		emit segmentUpChanged( v ) ;
	}

//
}

void TriangleGLDrawer::HSVtoRGB( size_t *r, size_t *g, size_t *b, float h, float s, float v ) const
{
	int i;
	float f, p, q, t;

	if( s == 0 )
	{
		// achromatic (grey)
		*r = *g = *b = static_cast<size_t>( v * 255. );
		return;
	}

	h /= 60.;                        // sector 0 to 5
	i = ( int )floor( h );
	f = h - i;                      // factorial part of h
	p = v * ( 1. - s );
	q = v * ( 1. - s * f );
	t = v * ( 1. - s * ( 1. - f ) );

	switch( i )
	{
		case 0:
			*r = static_cast<size_t>( round( v * 255. ) );
			*g = static_cast<size_t>( round( t * 255. ) );
			*b = static_cast<size_t>( round( p * 255. ) );
			break;
		case 1:
			*r = static_cast<size_t>( round( q * 255. ) );
			*g = static_cast<size_t>( round( v * 255. ) );
			*b = static_cast<size_t>( round( p * 255. ) );
			break;
		case 2:
			*r = static_cast<size_t>( round( p * 255. ) );
			*g = static_cast<size_t>( round( v * 255. ) );
			*b = static_cast<size_t>( round( t * 255. ) );
			break;
		case 3:
			*r = static_cast<size_t>( round( p * 255. ) );
			*g = static_cast<size_t>( round( q * 255. ) );
			*b = static_cast<size_t>( round( v * 255. ) );
			break;
		case 4:
			*r = static_cast<size_t>( round( t * 255. ) );
			*g = static_cast<size_t>( round( p * 255. ) );
			*b = static_cast<size_t>( round( v * 255. ) );
			break;
		default:                // case 5:
			*r = static_cast<size_t>( round( v * 255. ) );
			*g = static_cast<size_t>( round( p * 255. ) );
			*b = static_cast<size_t>( round( q * 255. ) );
			break;
	}
}

void TriangleGLDrawer::reset( QString f, const std::vector<std::pair<float, float> > & l )
{
	glDeleteLists( displayList, numberOfExtraFields + 2 ) ;
	limits = l ;
	valuesAtPoint = nullptr ;
	reader = nullptr ;

	mousePosOnLeftClick = QPoint( 0, 0 );

	leftDown = false;
	moving = false ;

	xtransleft = 0;
	ytransleft  = 0;

	zoom = 100 ;

	zpos = 1.5 ;
	currentSet = 0 ;

	openFile( f ) ;
	fracup = 10000;
	fracdown = 0;

	scale = 1 ;

	computeDisplayList();
}
TriangleGLDrawer::TriangleGLDrawer( QString f, const std::vector<std::pair<float, float> > & limits, QWidget *parent ) : QGLWidget( parent ), limits( limits )
{

	valuesAtPoint = nullptr ;
	reader = nullptr ;

	mousePosOnLeftClick = QPoint( 0, 0 );

	leftDown = false;
	moving = false ;

	xtransleft = 0;
	ytransleft  = 0;

	zoom = 100 ;

	zpos = 1.5 ;
	currentSet = 0 ;

	openFile( f ) ;
	fracup = 10000;
	fracdown = 0;

	scale = 1 ;

	minmaxinit = false ;


}

void TriangleGLDrawer::reset( std::vector<std::valarray<float> > * v, int np, int set, const std::vector<std::pair<float, float> > & l )
{
	glDeleteLists( displayList, numberOfExtraFields + 2 ) ;
	limits = l ;
	valuesAtPoint = v ;
	reader = nullptr ;

	mousePosOnLeftClick = QPoint( 0, 0 );

	leftDown = false;
	moving = false ;

	xtransleft = 0;
	ytransleft  = 0;

	zoom = 100 ;

	zpos = 1.5 ;
	currentSet = set ;

//	openFile(f) ;
	fracup = 10000;
	fracdown = 0;

// 	scale = 1 ;

	numberOfTriangles = ( *v )[0].size() ;
	numberOfPointsPerTriangle = np ;
	numberOfExtraFields = ( ( v->size() ) - 2 * np ) / np ;
	numberOfExtraTimePlanes = 0 ;
	currentTimePlane = 0 ;
	computeDisplayList();
	valUnderCursor = 0 ;
}

TriangleGLDrawer::TriangleGLDrawer( std::vector<std::valarray<float> > * v, int np, int set, const std::vector<std::pair<float, float> > & limits, QWidget *parent ) : QGLWidget( parent ), limits( limits )
{
	valUnderCursor = 0 ;
	valuesAtPoint = v ;
	reader = nullptr ;

	mousePosOnLeftClick = QPoint( 0, 0 );

	leftDown = false;
	moving = false ;

	xtransleft = 0;
	ytransleft  = 0;

	zoom = 100 ;

	zpos = 1.5 ;
	currentSet = set ;

//	openFile(f) ;
	fracup = 10000;
	fracdown = 0;

	scale = 1 ;

	minmaxinit = false ;

	numberOfTriangles = ( *v )[0].size() ;
	numberOfPointsPerTriangle = np ;
	numberOfExtraFields = ( ( v->size() ) - 2 * np ) / np ;
	numberOfExtraTimePlanes = 0 ;
	currentTimePlane = 0 ;
}


void TriangleGLDrawer::reset( TriangleDataReader *f, int set, const std::vector<std::pair<float, float> > & l )
{
	glDeleteLists( displayList, numberOfExtraFields + 2 ) ;
	limits = l ;
	valuesAtPoint = nullptr ;
	reader = nullptr ;

	mousePosOnLeftClick = QPoint( 0, 0 );

	leftDown = false;
	moving = false ;

	xtransleft = 0;
	ytransleft  = 0;

	zoom = 100 ;

	zpos = 1.5 ;
	currentSet = set ;
	currentDisplayList = set + 1 ;

	openFile( f ) ;
	fracup = 10000;
	fracdown = 0;

// 	scale = 1 ;

	computeDisplayList();
}

TriangleGLDrawer::TriangleGLDrawer( TriangleDataReader *f, int set, const std::vector<std::pair<float, float> > & limits, QWidget *parent ): QGLWidget( parent ), limits( limits )
{
	valuesAtPoint = nullptr ;
	reader = nullptr ;
	valUnderCursor = 0 ;
	mousePosOnLeftClick = QPoint( 0, 0 );

	leftDown = false;
	moving = false ;

	xtransleft = 0;
	ytransleft  = 0;

	zoom = 100 ;

	zpos = 1.5 ;
	currentSet = set ;
	currentDisplayList = set + 1 ;

	openFile( f ) ;
	fracup = 10000;
	fracdown = 0;

	scale = 1 ;

	minmaxinit = false ;
}


void TriangleGLDrawer::reset( QString f, int set, const std::vector<std::pair<float, float> > & l )
{
	glDeleteLists( displayList, numberOfExtraFields + 2 ) ;
	limits = l ;
	valuesAtPoint = nullptr ;
	reader = nullptr ;

	mousePosOnLeftClick = QPoint( 0, 0 );

	leftDown = false;
	moving = false ;

	xtransleft = 0;
	ytransleft  = 0;

	zoom = 100 ;

	zpos = 1.5 ;
	currentSet = set ;

	openFile( f ) ;
	fracup = 10000;
	fracdown = 0;

// 	scale = 1 ;

	computeDisplayList();
}

TriangleGLDrawer::TriangleGLDrawer( QString f, int set, const std::vector<std::pair<float, float> > & limits, QWidget *parent ) : QGLWidget( parent ), limits( limits )
{
	valUnderCursor = 0 ;
	valuesAtPoint = nullptr ;
	reader = nullptr ;

	mousePosOnLeftClick = QPoint( 0, 0 );

	leftDown = false;
	moving = false ;

	xtransleft = 0;
	ytransleft  = 0;

	zoom = 100 ;

	zpos = 1.5 ;
	currentSet = set ;

	openFile( f ) ;
	fracup = 10000;
	fracdown = 0;

	scale = 1 ;

	minmaxinit = false ;

}

void TriangleGLDrawer::reset()
{
	glDeleteLists( displayList, numberOfExtraFields + 2 ) ;

	valuesAtPoint = new std::vector< std::valarray<float> >( 0 ) ;
	reader = nullptr ;

	mousePosOnLeftClick = QPoint( 0, 0 );

	leftDown = false;
	moving = false ;

	numberOfExtraFields = 0;
	numberOfTriangles = 0 ;
	numberOfPointsPerTriangle = 0;

	xtransleft = 0;
	ytransleft  = 0;

	zoom = 100 ;

	zpos = 1.5 ;
	currentSet = 0 ;

	fracup = 10000;
	fracdown = 0;

// 	scale = 1 ;
	computeDisplayList();
}

TriangleGLDrawer::TriangleGLDrawer( QWidget *parent ) : QGLWidget( parent )
{
	valUnderCursor = 0 ;
	minmaxinit = false ;

	valuesAtPoint = new std::vector< std::valarray<float> >( 0 ) ;
	reader = nullptr ;

	mousePosOnLeftClick = QPoint( 0, 0 );

	leftDown = false;
	moving = false ;

	numberOfExtraFields = 0;
	numberOfTriangles = 0 ;
	numberOfPointsPerTriangle = 0;

	xtransleft = 0;
	ytransleft  = 0;

	zoom = 100 ;

	zpos = 1.5 ;
	currentSet = 0 ;

	fracup = 10000;
	fracdown = 0;

	scale = 1 ;
}

TriangleGLDrawer::~TriangleGLDrawer()
{
	glDeleteLists( displayList, numberOfExtraFields + 2 ) ;
	makeCurrent();
}


float TriangleGLDrawer::getXtrans() const
{
	return  -( float )xtransleft / ( ( float )width() * 1000. );
}

float TriangleGLDrawer::getYtrans() const
{
	return -( float )ytransleft / ( ( float )height() * 1000. );
}
