#include "triangleGlDrawer.h"

void TriangleGLDrawer::computeDisplay( ) const {
	
	glMatrixMode(GL_MODELVIEW);
	
	glColor4f(1, 1, 1, 0.5) ;
	glBegin(GL_LINE_LOOP) ;
		glVertex2f( 0.50, 0.50) ;
		glVertex2f(-0.50, 0.50) ;
		glVertex2f(-0.50, -0.50) ;
		glVertex2f( 0.50, -0.50) ;
	glEnd() ;
	glCallList( currentDisplayList) ;
}

void TriangleGLDrawer::mouseReleaseEvent(QMouseEvent *event) {
	moving = false ;
	leftDown = false ;
	if(event->button() == Qt::LeftButton)
		mousePosOnLeftClick = QPoint(0,0) ;

	paintGL() ;
}

void TriangleGLDrawer::mouseMoveEvent(QMouseEvent *event) {
	if(leftDown)
	{
		xtransleft -= (int)round((float)(mousePosOnLeftClick.x()-event->x())*(100000./(float)zoom)) ;
		ytransleft += (int)round((float)(mousePosOnLeftClick.y()-event->y())*(100000./(float)zoom)) ;
		emit xTranslationChanged(xtransleft) ;
		emit yTranslationChanged(ytransleft) ;
		mousePosOnLeftClick = event->pos() ;
	}
	
	paintGL() ;
}

void TriangleGLDrawer::mousePressEvent(QMouseEvent *event) {
// 	grabKeyboard() ;
	
	moving = true ;
	if(event->button() == Qt::LeftButton)
	{
		leftDown = true ;
		mousePosOnLeftClick = event->pos() ;
	}
}

void TriangleGLDrawer::keyPressEvent ( QKeyEvent * event ) {
	if(event->key() == Qt::Key_Plus)
	{
		zoom = (int)round(zoom * 1.1) ;
		setZoom(zoom) ;
	}
	else if (event->key() == Qt::Key_Minus)
	{
		zoom = (int)round(zoom / 1.1) ;
		setZoom(zoom) ;
	}
	else if(event->key() == Qt::Key_D)
	{
		currentDisplayList++ ;
		currentDisplayList-=displayList ;
		currentDisplayList%=numberOfExtraFields ;
		currentDisplayList+=displayList ;
	}
	else
	{
		event->ignore() ;
	}
	
	paintGL() ;
}

void TriangleGLDrawer::wheelEvent(QWheelEvent * event ) {
	
	int newZoom = zoom ;
	if(event->delta() > 0)
		newZoom = (float)zoom * 1.1 ;
	else
		newZoom = (float)zoom * 0.9 ;
	zoom = newZoom ;
	setZoom(newZoom);
}

void TriangleGLDrawer::setSet(int set) {
// 	set -= displayList ;
	set %= numberOfExtraFields ;
	set += displayList ;
	currentDisplayList = set ;
	paintGL() ;
	
	emit setChanged(set-displayList) ;
}

void TriangleGLDrawer::setXTranslation(int trans){
	xtransleft = trans ;
	emit xTranslationChanged(xtransleft) ;
}

void TriangleGLDrawer::setYTranslation(int trans){
	ytransleft = trans ;
	emit yTranslationChanged(ytransleft) ;
}

void TriangleGLDrawer::paintGL() {
	
	QTime startTime = QTime::currentTime();
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);  // Clear Screen And Depth Buffer
	
	glMatrixMode(GL_MODELVIEW);                             // Select The Modelview Matrix
	
	glLoadIdentity();
	gluLookAt(getXtrans(), getYtrans(), zpos, getXtrans(), getYtrans(), 0, 0, 1, 0) ;
	computeDisplay() ;
	
	glColor3f(1., 1., 1.) ;
	int elapsedTime = startTime.msecsTo(QTime::currentTime());
	renderText(10, 20, QString("%0 fps").arg(1000.0f / (float)elapsedTime));
	
	glFlush();
	swapBuffers() ;
}

void TriangleGLDrawer::resizeGL(int w, int h) {
	
	if (h==0)                                        // Prevent A Divide By Zero By
		h=1;                                           // Making Height Equal One
	
	if(h > w)
		glViewport(0,0, h, h);
	else
		glViewport(0,0, w, w);


	
// 	glPointSize(std::max(0.5*(float)width()/(float)columns, 1.));
	
	setZoom(zoom) ;	
}

void TriangleGLDrawer::setZoom(int percent) {
	size_t p ;
	if(percent == 0)
		p = 1 ;
	else
		p = percent ;
	
	float fact = 45.*((float)100/(float)p) ;
	glMatrixMode (GL_PROJECTION);
	glLoadIdentity ();
	gluPerspective(fact,1.,1.f,30.0f);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity(); 
		
	paintGL() ;
	
	emit zoomChanged(p) ;
}

void TriangleGLDrawer::openFile(const QString f) {

	delete valuesAtPoint ;
	TriangleDataReader reader(f) ;
	
	numberOfTriangles = reader.numberOfTriangles() ;	
	numberOfPointsPerTriangle = reader.numberOfPointsPerTriangle() ;	
	numberOfExtraFields = reader.numberOfExtraFields() ;
	
	valuesAtPoint = reader.data() ;
}

void TriangleGLDrawer::initializeGL() {
	
	//get the extended GL functions
	//glBlendEquationSeparate = (void (*)(GLenum, GLfloat))glXGetProcAddressARB((const GLubyte*)"glBlendEquationSeparate") ;
	//glPointParameterfARB = (void (*)(GLenum, GLfloat))glXGetProcAddressARB((const GLubyte*)"glPointParameterfARB") ;
	//glPointParameterfvARB = (void (*)(GLenum, GLfloat *))glXGetProcAddressARB((const GLubyte*)"glPointParameterfvARB") ;
	//glActiveTextureARB = (void (*)(GLuint))glXGetProcAddressARB((const GLubyte*)"glActiveTextureARB") ;
	//glTexEnvf = (void (*)(GLenum, GLenum, GLboolean))glXGetProcAddressARB((const GLubyte*)"glTexEnvf") ;
	
	displayList = glGenLists(numberOfExtraFields) ;
	std::cout << "numberOfExtraFields = " << numberOfExtraFields << std::endl ;
	std::cout << "displayList = " << displayList << std::endl ;
	currentDisplayList = displayList ;
	glViewport(0, 0, 600, 600) ;
	glShadeModel(GL_SMOOTH);   // Enables Smooth Shading
	glEnable(GL_LINE_SMOOTH) ;
	glEnable(GL_POLYGON_SMOOTH);
	glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST) ;
	
// 	glPointSize(std::max(0.4*(float)width()/(float)columns, 1.));
	glClearColor(0.0f,0.0f,0.0f,0.0f);                                      // Black Background
	glClearDepth(1.0f);                                                     // Depth Buffer Setup
	glDisable(GL_DEPTH_TEST);
	
	
	
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA,GL_ONE); 
	
	computeDisplayList() ;
}

QSize TriangleGLDrawer::minimumSizeHint() const{
	return QSize(50, 50);
}

QSize TriangleGLDrawer::sizeHint() const{
	return QSize(600, 600);
}
 
void TriangleGLDrawer::grab()
{

	paintGL() ;
	
	QImage grab = grabFrameBuffer() ;
	QString s(fileName) ;

	s+= QString(".png") ;
	
	grab.save(s, "png") ;
		
}
 
 void TriangleGLDrawer::computeDisplayList() {
	
	 float max_x = (*valuesAtPoint)[0][0] ;
	 float min_x = (*valuesAtPoint)[0][0] ; 
	 float max_y = (*valuesAtPoint)[1][0] ;
	 float min_y = (*valuesAtPoint)[1][0] ;

	 for(size_t i =0 ; i< numberOfTriangles ; i++)
	 {
	   for(size_t j = 0 ; j < numberOfPointsPerTriangle ; j++)
	   {
	     if( (*valuesAtPoint)[2*j][i] > max_x)
	       max_x = (*valuesAtPoint)[2*j][i] ;
	     if( (*valuesAtPoint)[2*j+1][i] > max_y)
	       max_y = (*valuesAtPoint)[2*j+1][i] ;
	     if( (*valuesAtPoint)[2*j][i] < min_x)
	       min_x = (*valuesAtPoint)[2*j][i];
	     if( (*valuesAtPoint)[2*j+1][i] < min_y)
	       min_y = (*valuesAtPoint)[2*j+1][i] ;
	   }
	 }
	 
	 for(size_t N = 0 ; N < numberOfExtraFields ;N++)
	 {
		glNewList( displayList+N, GL_COMPILE) ;
		 
		float max_val = (*valuesAtPoint)[(2+N)*numberOfPointsPerTriangle][0] ; 
		float min_val = (*valuesAtPoint)[(2+N)*numberOfPointsPerTriangle][0] ;
		 
		std::vector<float> vals ;
		for(size_t i =0 ; i< numberOfTriangles ; i++)
		{
		  for(size_t j = 0 ; j < numberOfPointsPerTriangle ; j++)
		  {
			vals.push_back((*valuesAtPoint)[(2+N)*numberOfPointsPerTriangle+j][i]) ;
		    if( vals.back() > max_val)
		      max_val = vals.back() ;
		    if( vals.back() < min_val)
		      min_val = vals.back();

		  }
		}
		 
		 
		 std::sort(vals.begin(), vals.end()) ;
		 vals.erase(std::unique(vals.begin(), vals.end()), vals.end()) ;
		 if(limits.size() <= N)
		 {
		   limits.push_back(std::make_pair(vals[.001*(vals.size()-1)], vals[.999*(vals.size()-1)])) ;
		   if(vals.front() > 0 &&  vals.back() > 0)
		      limits[N].first = 0 ;
		 }
		 
		 if(limits[N].first > vals[.001*(vals.size()-1)])
		   limits[N].first = vals[.001*(vals.size()-1)] ;
		 if(limits[N].second < vals[.999*(vals.size()-1)])
		   limits[N].second = vals[.999*(vals.size()-1)] ;
		 if(vals.front() > 0 &&  vals.back() > 0)
		      limits[N].first = 0 ;
		 min_val = limits[N].first ;
		 max_val = limits[N].second ;
		 bool logplot = false ;
		 if(min_val >= 0 &&  max_val >= 0)
		 {
				min_val = 0 ;
				logplot = true ;
				std::cout <<"logplot" << std::endl ;
		 }
		 
		 double halfval = (min_val + max_val)*.5 ;
		 double maxdelta = std::max(max_x-min_x, max_y-min_y) ;
		 double mindelta = std::min(max_x-min_x, max_y-min_y) ;
		 double cx = 1 ;
		 double cy = 1 ;
		 double mag = 10 ;
		 if(max_x-min_x > max_y-min_y)
			 cy = mindelta/maxdelta ;
		 else
			 cx = mindelta/maxdelta ;
		for(size_t i = 0 ; i< numberOfTriangles ; i++)
		{
			
			size_t r, g, b ;

			if((*valuesAtPoint)[(2+N)*numberOfPointsPerTriangle][i] <= halfval)
			{
				if((*valuesAtPoint)[(2+N)*numberOfPointsPerTriangle][i] >= min_val && (*valuesAtPoint)[(2+N)*numberOfPointsPerTriangle][i] <= max_val)
				{
					glBegin(GL_TRIANGLES) ;
					for(size_t j = 0 ; j < numberOfPointsPerTriangle ; j++)
					{
						float v = (std::max(std::min((*valuesAtPoint)[(2+N)*numberOfPointsPerTriangle+j][i], max_val), min_val) - min_val)/(max_val-min_val);
						HSVtoRGB(&r, &g, &b, 330.-v*360., 1./*-0.5*exp(v+1)/exp(2)*/, 1.) ;
						glColor4ub(r, g, b, 255) ;
						double dx = (*valuesAtPoint)[(2)*numberOfPointsPerTriangle+j][i]*mag ;
						double dy = (*valuesAtPoint)[(3)*numberOfPointsPerTriangle+j][i]*mag ;
						glVertex2f(((*valuesAtPoint)[j*2][i]-min_x)/maxdelta-0.5*cx + dx, ((*valuesAtPoint)[j*2+1][i]-min_y)/maxdelta-0.5*cy +dy) ;
		// 				std::cout << ((*valuesAtPoint)[j*2][i]-min_x)/(max_x-min_x)-0.5 << ", "<< ((*valuesAtPoint)[j*2+1][i]-min_y)/(max_y-min_y)-0.5 << std::endl ;
					}
					glEnd() ;
				}
				else if((*valuesAtPoint)[(2+N)*numberOfPointsPerTriangle][i] < min_val)
				{
					glBegin(GL_TRIANGLES) ;
					for(size_t j = 0 ; j < numberOfPointsPerTriangle ; j++)
					{
						HSVtoRGB(&r, &g, &b, 0. , 1., .5) ;
						glColor4ub(r, g, b, 255) ;
						
						double dx = (*valuesAtPoint)[(2)*numberOfPointsPerTriangle+j][i]*mag ;
						double dy = (*valuesAtPoint)[(3)*numberOfPointsPerTriangle+j][i]*mag ;
						
						glVertex2f(((*valuesAtPoint)[j*2][i]-min_x)/maxdelta-0.5*cx +dx, ((*valuesAtPoint)[j*2+1][i]-min_y)/maxdelta-0.5*cy +dy) ;
		// 				std::cout << ((*valuesAtPoint)[j*2][i]-min_x)/(max_x-min_x)-0.5 << ", "<< ((*valuesAtPoint)[j*2+1][i]-min_y)/(max_y-min_y)-0.5 << std::endl ;
					}
					glEnd() ;
				}
				else if((*valuesAtPoint)[(2+N)*numberOfPointsPerTriangle][i] > max_val)
				{
					glBegin(GL_TRIANGLES) ;
					for(size_t j = 0 ; j < numberOfPointsPerTriangle ; j++)
					{
						HSVtoRGB(&r, &g, &b, 1. , 1., .5) ;
						glColor4ub(r, g, b, 255) ;
						
						double dx = (*valuesAtPoint)[(2)*numberOfPointsPerTriangle+j][i]*mag ;
						double dy = (*valuesAtPoint)[(3)*numberOfPointsPerTriangle+j][i]*mag ;
						
						glVertex2f(((*valuesAtPoint)[j*2][i]-min_x)/maxdelta-0.5*cx+dx, ((*valuesAtPoint)[j*2+1][i]-min_y)/maxdelta-0.5*cy+dy) ;
		// 				std::cout << ((*valuesAtPoint)[j*2][i]-min_x)/(max_x-min_x)-0.5 << ", "<< ((*valuesAtPoint)[j*2+1][i]-min_y)/(max_y-min_y)-0.5 << std::endl ;
					}
					glEnd() ;
				}
			}
			else
			{
				if((*valuesAtPoint)[(2+N)*numberOfPointsPerTriangle][i] >= min_val && (*valuesAtPoint)[(2+N)*numberOfPointsPerTriangle][i] <= max_val)
				{
					glBegin(GL_TRIANGLES) ;
					for(size_t j = 0 ; j < numberOfPointsPerTriangle ; j++)
					{
						float v = (std::max(std::min((*valuesAtPoint)[(2+N)*numberOfPointsPerTriangle+j][i], max_val), min_val) - min_val)/(max_val-min_val);
						HSVtoRGB(&r, &g, &b, 330.-v*360., 1./*-0.5*exp(v+1)/exp(2)*/, 1.) ;
						glColor4ub(r, g, b, 255) ;
						
						double dx = (*valuesAtPoint)[(2)*numberOfPointsPerTriangle+j][i]*mag ;
						double dy = (*valuesAtPoint)[(3)*numberOfPointsPerTriangle+j][i]*mag ;
						
						glVertex2f(((*valuesAtPoint)[j*2][i]-min_x)/maxdelta-0.5*cx +dx, ((*valuesAtPoint)[j*2+1][i]-min_y)/maxdelta-0.5*cy+dy) ;
		// 				std::cout << ((*valuesAtPoint)[j*2][i]-min_x)/(max_x-min_x)-0.5 << ", "<< ((*valuesAtPoint)[j*2+1][i]-min_y)/(max_y-min_y)-0.5 << std::endl ;
					}
					glEnd() ;
				}
				else if((*valuesAtPoint)[(2+N)*numberOfPointsPerTriangle][i] < min_val)
				{
					glBegin(GL_TRIANGLES) ;
					for(size_t j = 0 ; j < numberOfPointsPerTriangle ; j++)
					{
						HSVtoRGB(&r, &g, &b, 0. , 1., .5) ;
						glColor4ub(r, g, b, 255) ;
						
						double dx = (*valuesAtPoint)[(2)*numberOfPointsPerTriangle+j][i]*mag ;
						double dy = (*valuesAtPoint)[(3)*numberOfPointsPerTriangle+j][i]*mag ;
						
						glVertex2f(((*valuesAtPoint)[j*2][i]-min_x)/maxdelta-0.5*cx +dx, ((*valuesAtPoint)[j*2+1][i]-min_y)/maxdelta-0.5*cy +dy) ;
		// 				std::cout << ((*valuesAtPoint)[j*2][i]-min_x)/(max_x-min_x)-0.5 << ", "<< ((*valuesAtPoint)[j*2+1][i]-min_y)/(max_y-min_y)-0.5 << std::endl ;
					}
					glEnd() ;
				}
				else if((*valuesAtPoint)[(2+N)*numberOfPointsPerTriangle][i] > max_val)
				{
					glBegin(GL_TRIANGLES) ;
					for(size_t j = 0 ; j < numberOfPointsPerTriangle ; j++)
					{
						HSVtoRGB(&r, &g, &b, 1. , 1., .5) ;
						glColor4ub(r, g, b, 255) ;
						
						double dx = (*valuesAtPoint)[(2)*numberOfPointsPerTriangle+j][i]*mag ;
						double dy = (*valuesAtPoint)[(3)*numberOfPointsPerTriangle+j][i]*mag ;
						
						glVertex2f(((*valuesAtPoint)[j*2][i]-min_x)/maxdelta-0.5*cx +dx, ((*valuesAtPoint)[j*2+1][i]-min_y)/maxdelta-0.5*cy +dy) ;
		// 				std::cout << ((*valuesAtPoint)[j*2][i]-min_x)/(max_x-min_x)-0.5 << ", "<< ((*valuesAtPoint)[j*2+1][i]-min_y)/(max_y-min_y)-0.5 << std::endl ;
					}
					glEnd() ;
				}
				
			}
		}
		
		glEndList() ;
	 }
}

void TriangleGLDrawer::HSVtoRGB( size_t *r, size_t *g, size_t *b, float h, float s, float v ) const {
	int i;
	float f, p, q, t;
	if( s == 0 ) {
                // achromatic (grey)
		*r = *g = *b = static_cast<size_t>(v*255.);
		return;
	}
	h /= 60.;                        // sector 0 to 5
	i = (int)floor( h );
	f = h - i;                      // factorial part of h
	p = v * ( 1. - s );
	q = v * ( 1. - s * f );
	t = v * ( 1. - s * ( 1. - f ) );
	switch( i ) {
	case 0:
		*r = static_cast<size_t>(round(v*255.));
		*g = static_cast<size_t>(round(t*255.));
		*b = static_cast<size_t>(round(p*255.));
		break;
	case 1:
		*r = static_cast<size_t>(round(q*255.));
		*g = static_cast<size_t>(round(v*255.));
		*b = static_cast<size_t>(round(p*255.));
		break;
	case 2:
		*r = static_cast<size_t>(round(p*255.));
		*g = static_cast<size_t>(round(v*255.));
		*b = static_cast<size_t>(round(t*255.));
		break;
	case 3:
		*r = static_cast<size_t>(round(p*255.));
		*g = static_cast<size_t>(round(q*255.));
		*b = static_cast<size_t>(round(v*255.));
		break;
	case 4:
		*r = static_cast<size_t>(round(t*255.));
		*g = static_cast<size_t>(round(p*255.));
		*b = static_cast<size_t>(round(v*255.));
		break;
	default:                // case 5:
		*r = static_cast<size_t>(round(v*255.));
		*g = static_cast<size_t>(round(p*255.));
		*b = static_cast<size_t>(round(q*255.));
		break;
	}
}

TriangleGLDrawer::TriangleGLDrawer(QString f, QWidget *parent) : QGLWidget(parent) {
	
	valuesAtPoint = NULL ;
	
	mousePosOnLeftClick = QPoint(0,0);
	
	leftDown = false;
	moving = false ;
	
	xtransleft = 0;
	ytransleft  = 0;
	
	zoom = 100 ;
	
	zpos = 1.5 ;
	
	openFile(f) ;

}

TriangleGLDrawer::TriangleGLDrawer(QWidget *parent) : QGLWidget(parent) {
	
	valuesAtPoint = new std::vector< std::valarray<float> >(0) ;
	
	mousePosOnLeftClick = QPoint(0,0);
	
	leftDown = false;
	moving = false ;
	
	numberOfExtraFields = 0;
	numberOfTriangles = 0 ;
	numberOfPointsPerTriangle = 0;
	
	xtransleft = 0;
	ytransleft  = 0;
	
	zoom = 100 ;
	
	zpos = 1.5 ;	
	
}

TriangleGLDrawer::~TriangleGLDrawer()  {
	glDeleteLists(displayList, numberOfExtraFields+2) ;
	makeCurrent();
}


float TriangleGLDrawer::getXtrans() const {
	return  -(float)xtransleft/((float)width()*1000.);
}

float TriangleGLDrawer::getYtrans() const {
	return -(float)ytransleft/((float)height()*1000.);
}
