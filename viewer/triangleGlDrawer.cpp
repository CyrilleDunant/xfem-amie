#include "triangleGlDrawer.h"

void TriangleGLDrawer::computeDisplay( ) const {
	
	glMatrixMode(GL_MODELVIEW);
	
	glColor4f(1, 1, 1, 0.5) ;
	glBegin(GL_LINE_LOOP) ;
		glVertex2f( 0.50*.8-.15, 0.50*.8) ;
		glVertex2f(-0.50*.8-.15, 0.50*.8) ;
		glVertex2f(-0.50*.8-.15, -0.50*.8) ;
		glVertex2f( 0.50*.8-.15, -0.50*.8) ;
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
	if(set < 0)
		set = numberOfExtraFields-1 ;
	if(set == numberOfExtraFields)
		set = 0 ;
	
	currentDisplayList = set+displayList ;
	currentSet = set ;
	paintGL() ;
	
	emit setChanged(set) ;
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
	size_t r, g, b ;
	
	if(!limits.empty())
	{
		glBegin(GL_QUAD_STRIP) ;
		for(double i = 1. ; i > 0 ; i -= 0.01)
		{
			HSVtoRGB(&r, &g, &b, 180., 0., 0.1+0.8*(1.-i)) ;
			glColor4ub(r, g, b, 255) ;
			glVertex2f((.805-0.5) , (i-0.5)*.7 ) ;
			glVertex2f((.835-0.5) , (i-0.5)*.7 ) ;
		}
		glEnd() ;
		glColor4ub(255, 0, 0, 255) ;
		for(double i = 1. ; i > 0 ; i -= 0.1)
		{
			renderText((.805-0.5)+0.05, 
									(i-0.5)*.7,  
									0.,
									QString("%0").arg(i*(limits[currentSet].first-limits[currentSet].second)+limits[currentSet].second));
		}
	}
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
		

		 
		 if(limits.size() <= N)
		 {
		   limits.push_back(std::make_pair(max_val, min_val)) ;
		 }
		 if(limits[N].first < max_val)
		   limits[N].first = max_val ;
		 if(limits[N].second > min_val)
		   limits[N].second = min_val ;
		 
		 if(std::abs(limits[N].first+limits[N].second)/(limits[N].first- limits[N].second)< 1e-6 || isnan(std::abs(limits[N].first+limits[N].second)/(limits[N].first- limits[N].second)))
		 {
			 limits[N].first = 1e-6 ;
			 limits[N].second = -1e-6 ;
		 }

		 max_val = limits[N].first ;
		 min_val = limits[N].second ;
		 
		 double maxdelta = std::max(max_x-min_x, max_y-min_y)/0.8 ;
		 double mindelta = std::min(max_x-min_x, max_y-min_y)/0.8 ;
		 double cx = 0.7 ;
		 double cy = 0.7 ;
		 double mag = 10 ;
		 if(max_x-min_x > max_y-min_y)
			 cy = 0.7*mindelta/maxdelta ;
		 else
			 cx = 0.7*mindelta/maxdelta ;
		 size_t r, g, b ;
		for(size_t i = 0 ; i< numberOfTriangles ; i++)
		{
					glBegin(GL_TRIANGLES) ;
					for(size_t j = 0 ; j < numberOfPointsPerTriangle ; j++)
					{
						float v = ((*valuesAtPoint)[(2+N)*numberOfPointsPerTriangle+j][i] - min_val)/(max_val-min_val);
						HSVtoRGB(&r, &g, &b, 180., 0./*-0.5*exp(v+1)/exp(2)*/, 0.9-v*0.8) ;
						glColor4ub(r, g, b, 255) ;
						double dx = (*valuesAtPoint)[(2)*numberOfPointsPerTriangle+j][i]*mag ;
						double dy = (*valuesAtPoint)[(3)*numberOfPointsPerTriangle+j][i]*mag ;
						glVertex2f(((*valuesAtPoint)[j*2][i]-min_x)/maxdelta-0.5*cx + dx-0.2, ((*valuesAtPoint)[j*2+1][i]-min_y)/maxdelta-0.5*cy +dy) ;
					}
					glEnd() ;
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

TriangleGLDrawer::TriangleGLDrawer(QString f, const std::vector<std::pair<float, float> > & limits, QWidget *parent) : QGLWidget(parent), limits(limits) {
	
	valuesAtPoint = NULL ;
	
	mousePosOnLeftClick = QPoint(0,0);
	
	leftDown = false;
	moving = false ;
	
	xtransleft = 0;
	ytransleft  = 0;
	
	zoom = 100 ;
	
	zpos = 1.5 ;
	currentSet = 0 ;
	
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
	currentSet = 0 ;
	
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
