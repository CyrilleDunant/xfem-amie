#include "voxelGlDrawer.h"

void VoxelGLDrawer::computeDisplay( ) {
	glMatrixMode(GL_MODELVIEW);

	GLint targetList = currentDisplayList ;
	
	if(alpha == 255)
		targetList = displayList ;
	else
	{
		if(xangle >= 15 && xangle <= 180-15 )
			targetList = displayList ;
		else if(xangle > 180+15 && xangle < 360-15)
			targetList =displayList+1 ;
		else if(xangle > 180-15 && xangle <= 180+15 )
		{
			if(yangle < 90-15 || yangle > 270+15)
				targetList =displayList+3 ;
			else if(yangle > 90+15 && yangle < 270-15)
				targetList =displayList+4 ;
			else if(yangle >= 90-15 && yangle <= 90+15)
				targetList =displayList+8 ;
			else
				targetList =displayList+6 ;
		}
		else
		{
			if(yangle < 90-15 || yangle > 270+15)
				targetList =displayList+2 ;
			else if(yangle > 90+15 && yangle < 270-15)
				targetList =displayList+5 ;
			else if(yangle >= 90-15 && yangle <= 90+15)
				targetList = displayList+7 ;
			else 
				targetList = displayList+9 ;
		}
	}
	if(targetList != currentDisplayList)
	{
		currentDisplayList = targetList ;
	}
	glCallList(targetList) ;
	std::cout << "\r"<< xangle << "  " << yangle << std::flush ;
// glEnable(GL_DEPTH_TEST);
// 	glDisable(GL_TEXTURE_2D) ;
// 	glColor4f(1, 1, 1, 0.5f) ;
// 	glBegin(GL_LINE_LOOP) ;
// 	glVertex3f(0.5f, size_y*0.5f, size_z*0.5f) ;
// 	glVertex3f(0.5f,-size_y*0.5f, size_z*0.5f) ;
// 	glVertex3f(0.5f,-size_y*0.5f, -size_z*0.5f) ;
// 	glVertex3f(0.5f, size_y*0.5f, -size_z*0.5f) ;
// 	glEnd() ;
// 	glBegin(GL_LINE_LOOP) ;
// 	glVertex3f(-0.5f, size_y*0.5f, size_z*0.5f) ;
// 	glVertex3f(-0.5f,-size_y*0.5f, size_z*0.5f) ;
// 	glVertex3f(-0.5f,-size_y*0.5f, -size_z*0.5f) ;
// 	glVertex3f(-0.5f, size_y*0.5f, -size_z*0.5f) ;
// 	glEnd() ;
// 	glBegin(GL_LINES) ;
// 	glVertex3f(0.5f, size_y*0.5f, size_z*0.5f) ;
// 	glVertex3f(-0.5f,size_y*0.5f, size_z*0.5f) ;
// 	glVertex3f(0.5f, -size_y*0.5f, size_z*0.5f) ;
// 	glVertex3f(-0.5f,-size_y*0.5f, size_z*0.5f) ;
// 	glVertex3f(0.5f, -size_y*0.5f, -size_z*0.5f) ;
// 	glVertex3f(-0.5f,-size_y*0.5f, -size_z*0.5f) ;
// 	glVertex3f(0.5f, size_y*0.5f, -size_z*0.5f) ;
// 	glVertex3f(-0.5f,size_y*0.5f, -size_z*0.5f) ;
// 	glEnd() ;
// 	glDisable(GL_DEPTH_TEST);
// 	for(size_t i = 0 ; i < 64 ; i++)
// 		glCallList(displayList+i) ;

}

void VoxelGLDrawer::grab()
{
	for (size_t i = 0 ; i < 360 ; i++)
	{
		xangle++ ;
		xangle%=360 ;
		emit xAngleChanged(xangle) ;
		paintGL() ;
		
		QImage grab = grabFrameBuffer() ;
		QString s("shots/grab") ;
		if(i < 10)
			s+= QString("00") ;
		else if(i < 100)
			s+= QString("0") ;
		s+= QString("%1").arg(i) ;
		s+= QString(".png") ;
		
		grab.save(s, "png") ;
		
	}
}

bool VoxelGLDrawer::isVisible(const size_t i, const size_t j, const size_t k) const {
	
	return true ;
	
	float x = (float)i/7.f-0.5f ; float y = (float)j/7.f-0.5f ; float z = (float)k/7.f-0.5f ;
	if(std::abs(x) > .4 || 
	   std::abs(y) > .4 ||
	   std::abs(z) > .4
	  )
		return true ;
	return false ;
	
	GLfloat mat[4][4] ;
	
	GLfloat position[] = { getXtrans(), getYtrans(), zpos, 1.0f };
	GLfloat pos[] = { 0, 0, 0, 1.0f };
	
	glGetFloatv(GL_MODELVIEW_MATRIX, (GLfloat *)mat);
	for(size_t l = 0 ; l< 4 ; l++)
	{
		for(size_t m = 0 ; m< 4 ; m++)
		{
			pos[l] += position[m] * mat[l][m] ;
		}
	}
	
	std::valarray<float> dists(0.f, 8) ;
	dists[0] = (pos[0]-0.5f)*(pos[0]-0.5f) + (pos[1]-0.5f)*(pos[1]-0.5f) + (pos[2]-0.5f)*(pos[2]-0.5f) ;
	dists[1] = (pos[0]-0.5f)*(pos[0]-0.5f) + (pos[1]-0.5f)*(pos[1]-0.5f) + (pos[2]+0.5f)*(pos[2]+0.5f) ;
	dists[2] = (pos[0]-0.5f)*(pos[0]-0.5f) + (pos[1]+0.5f)*(pos[1]+0.5f) + (pos[2]-0.5f)*(pos[2]-0.5f) ;
	dists[3] = (pos[0]-0.5f)*(pos[0]-0.5f) + (pos[1]+0.5f)*(pos[1]+0.5f) + (pos[2]+0.5f)*(pos[2]+0.5f) ;
	dists[4] = (pos[0]+0.5f)*(pos[0]+0.5f) + (pos[1]-0.5f)*(pos[1]-0.5f) + (pos[2]-0.5f)*(pos[2]-0.5f) ;
	dists[5] = (pos[0]+0.5f)*(pos[0]+0.5f) + (pos[1]-0.5f)*(pos[1]-0.5f) + (pos[2]+0.5f)*(pos[2]+0.5f) ;
	dists[6] = (pos[0]+0.5f)*(pos[0]+0.5f) + (pos[1]+0.5f)*(pos[1]+0.5f) + (pos[2]-0.5f)*(pos[2]-0.5f) ;
	dists[7] = (pos[0]+0.5f)*(pos[0]+0.5f) + (pos[1]+0.5f)*(pos[1]+0.5f) + (pos[2]+0.5f)*(pos[2]+0.5f) ;
		
	return (pos[0]-x)*(pos[0]-x) + (pos[1]-y)*(pos[1]-y) + (pos[2]-z)*(pos[2]-z) > dists[2];
}

void VoxelGLDrawer::mouseReleaseEvent(QMouseEvent *event) {
	moving = false ;
	leftDown = false ;
	rightDown = false ;
	if(event->button() == Qt::LeftButton)
		mousePosOnLeftClick = QPoint(0,0) ;
	else if (event->button() == Qt::RightButton)
		mousePosOnRightClick = QPoint(0,0) ;
	paintGL() ;
}

void VoxelGLDrawer::mouseMoveEvent(QMouseEvent *event) {
	if(leftDown)
	{
		xtransleft -= mousePosOnLeftClick.x()-event->x() ;
		ytransleft += mousePosOnLeftClick.y()-event->y() ;
		emit xTranslationChanged(xtransleft) ;
		emit yTranslationChanged(ytransleft) ;
		mousePosOnLeftClick = event->pos() ;
	}
	else if (rightDown)
	{
		
		xangle += (int)round(180.*(float)(mousePosOnRightClick.x()-event->x())/(float)width()) + 360;
		yangle += (int)round(180.*(float)(mousePosOnRightClick.y()-event->y())/(float)height()) + 360;

		xangle%=360 ;
		yangle%=360 ;
		
		emit xAngleChanged(xangle) ;
		emit yAngleChanged(yangle) ;
		
		mousePosOnRightClick = event->pos() ;
	}
	
	paintGL() ;
}

void VoxelGLDrawer::mousePressEvent(QMouseEvent *event) {
// 	grabKeyboard() ;
	
	moving = true ;
	if(event->button() == Qt::LeftButton)
	{
		leftDown = true ;
		mousePosOnLeftClick = event->pos() ;
	}
	else if (event->button() == Qt::RightButton)
	{
		rightDown = true ;
		mousePosOnRightClick = event->pos() ;
	}
}

void VoxelGLDrawer::wheelEvent(QWheelEvent * event ) {
	setZoom((int)((float)m_zoom * (1.0f + (float)event->delta()/(120.0f*8.0f))));
}

void VoxelGLDrawer::keyPressEvent ( QKeyEvent * event ) {
	if(event->key() == Qt::Key_Plus)
	{
		int alpha = 2;
		for(size_t i = 0 ; i < sz ; i++)
		{
			if(restriction[i])
			{
				size_t __alpha__[4] = {i2RGBA(colour[i])} ;
				alpha = __alpha__[3] ;
				break ;
			}
		}
		alpha++ ;
		setAlpha(alpha) ;
	}
	else if (event->key() == Qt::Key_Minus)
	{
		int alpha = 2;
		for(size_t i = 0 ; i < sz ; i++)
		{
			if(restriction[i])
			{
				size_t __alpha__[4] = {i2RGBA(colour[i])} ;
				alpha = __alpha__[3] ;
				break ;
			}
		}
		alpha-- ;
		setAlpha(alpha) ;
	}
	else if(event->key() == Qt::Key_S)
	{
		if(!slice)
		{
			slice = true ;
			slice_dir = X_SLICE_DIRECTION ;
		}
		else if(slice_dir == X_SLICE_DIRECTION)
		{
			slice_pos = 0 ;
			slice_dir = Y_SLICE_DIRECTION ;
		}
		else if (slice_dir == Y_SLICE_DIRECTION)
		{
			slice_pos = 0 ;
			slice_dir = Z_SLICE_DIRECTION ;
		}
		else if(slice_dir == Z_SLICE_DIRECTION)
		{
			slice_dir = X_SLICE_DIRECTION ;
			slice_pos = 0 ;
			slice = false ;
		}
		computeDisplayList() ;
		paintGL() ;
	}
	else if(event->key() == Qt::Key_O)
	{
		computeDisplayList() ;
		paintGL() ;
	}
	else if(event->key() == Qt::Key_P)
	{
		slice_pos+=10 ;
		if(slice_dir == X_SLICE_DIRECTION)
		{
			slice_pos%=rows ;	
		}
		if(slice_dir == Y_SLICE_DIRECTION)
		{
			slice_pos%=columns ;
			
		}
		if(slice_dir == Z_SLICE_DIRECTION)
		{
			slice_pos%=strips ;
		}
		computeDisplayList() ;
		paintGL() ;
	}
	else
	{
		event->ignore() ;
	}
}

void VoxelGLDrawer::updateRestrictions() {
	restriction = true ;
	size_t del = 0 ;
	size_t c = 0;
	QEventLoop eloop ;
	if(pbar)
	{
		pbar->reset() ;
		connect(this, SIGNAL(progressed(int)), pbar, SLOT(setValue(int)));
	}
	
	
	for(size_t i = 0 ; i < rows ; i++)
	{
		for(size_t j = 0 ; j < columns ; j++)
		{
			for(size_t k = 0 ; k < strips ; k++)
			{
				size_t index = toIndex(i,j,k) ;
				c++ ;
				
				if(c % (sz/100) == 0)
				{
					int percent = (int)((float)c/(float)(sz)*100.f) ;
					emit progressed(percent) ;
					eloop.wakeUp() ;
					eloop.processEvents() ;
				}
				
				if(restriction[index])
				{
					if(isBoundary(i,j,k))
					{
						std::valarray<qint8> n = normalFromSurroundings(i,j,k) ;
						normal[index*3] = n[0] ;
						normal[index*3+1] = n[1] ;
						normal[index*3+2] = n[2] ;
					}
					else
					{
						restriction[index] = false ;
						del++ ;
					}
				}
			}
		}
	}
	emit progressed(100) ;
	eloop.processEvents() ;
	eloop.exit() ;
	std::cout << " Added " << del << " restrictions." << std::endl ;
}

bool VoxelGLDrawer::isBoundary(const size_t& x, const size_t& y, const size_t& z) const {
	
	
	if( x == 0     || 
	x == rows-1    ||
	y == 0         ||
	y == columns-1 ||
	z == 0         ||
	z == strips -1
	)
		return true ;
		
	
// 	return false ;
			
// 	qint8 v = (*valuesAtPoint)[m_currentField][toIndex(x,y,z)] ;

	size_t index = toIndex(x,y,z) ;
	if(!isInRange(index))
		return false ;
	
	for(int i = std::max((int)(x)-1, 0) ; i < (int)std::min(x+2, rows) ; i++)
	{
		if(i == (int)x)
			i++ ;
		for(int j = std::max((int)(y)-1,0) ; j < (int)std::min(y+2, columns) ; j++)
		{
			if(j == (int)y)
				j++ ;
			for(int k = std::max((int)(z)-1, 0) ; k < (int)std::min(z+2, strips) ; k++)
			{
				if(k == (int)z)
					k++ ;

				if(!isInRange(toIndex(i,j,k)))
					return true ;
				
			}
		}
	}
	return false ;
	
}

void VoxelGLDrawer::setAlpha(int a){
	alpha = std::max(1,std::min(a, 255)) ;
	if(alpha == 255)
		glEnable(GL_DEPTH_TEST);
	else
		glDisable(GL_DEPTH_TEST);
	
	emit alphaChanged(alpha) ;
	computeDisplayList() ;
	paintGL() ;
}

void VoxelGLDrawer::setXAngle(int angle) {
	xangle = angle + 360;
	xangle%=360 ;
	emit xAngleChanged(xangle) ;
	paintGL() ;
}

void VoxelGLDrawer::setYAngle( int angle) {
	yangle = angle +360;
	yangle%=360 ;
	emit yAngleChanged(yangle) ;
	paintGL() ;
}

void VoxelGLDrawer::setZAngle( int angle) {
	zangle = angle +360;
	zangle%=360 ;
	emit zAngleChanged(zangle) ;
	paintGL() ;
}

void VoxelGLDrawer::setXTranslation(int trans){
	xtransleft = trans ;
	emit xTranslationChanged(xtransleft) ;
}

void VoxelGLDrawer::setYTranslation(int trans){
	ytransleft = trans ;
	emit yTranslationChanged(ytransleft) ;
}

void VoxelGLDrawer::paintGL() {
	
// 	QTime startTime = QTime::currentTime();
	GLfloat light_position[] = { 0.0, 0.0, zpos, 0 };
	
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);  // Clear Screen And Depth Buffer
	
	glMatrixMode(GL_MODELVIEW);                             // Select The Modelview Matrix
	glLoadIdentity();
	
	gluLookAt(getXtrans(), getYtrans(), zpos, getXtrans(), getYtrans(), 0, 0, 1, 0) ;
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	glRotatef(getXAngle() , 0, -1, 0) ;
	glRotatef(getYAngle() , -1, 0,0) ;

	computeDisplay()  ;
	
	glColor3f(1., 1., 1.) ;
// 	int elapsedTime = startTime.msecsTo(QTime::currentTime());
// 	renderText(10, 20, QString("%0 fps").arg(1000.0f / (float)elapsedTime));
	
	glFlush();
	swapBuffers() ;
}

void VoxelGLDrawer::resizeGL(int w, int h) {
	
	if (h==0)                                          // Prevent A Divide By Zero By
		h=1;   // Making Height Equal One
	
	glViewport(0, 0, w, h);   // Reset The Current Viewport

	
// 	glPointSize(2.0f*(float)width()/(float)columns);
// 	glPointParameterfARB( GL_POINT_SIZE_MIN_ARB, .85*(float)width()/(float)columns );
	
	float fact = 45.f*(100.f/(float)m_zoom) ;
	glPointSize(45.f/fact*(float)width()/(float)columns);
	glMatrixMode (GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(fact,(float)width()/(float)height(),1.,30.0f);
	// 	glOrtho(-1., 1., -1., 1., 0.1, 20.) ;
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity(); 
	
	paintGL() ;
	
// 	setZoom(100) ;	
}

void VoxelGLDrawer::setZoom(int percent) {
	size_t p = percent;
	if(percent < 10)
		p = 10 ;
	else if(percent > 1000)
		p = 1000 ;
	
	m_zoom = p ;
	float fact = 45.*(100.f/(float)p) ;
	glMatrixMode (GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(fact,(float)width()/(float)height(),1.,30.0f);
	// 	glOrtho(-1., 1., -1., 1., 0.1, 20.) ;
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity(); 
	
// 	glPointParameterfARB( GL_POINT_SIZE_MIN_ARB, .85*(45./fact)*(float)width()/(float)columns );
	glPointSize(40.f/fact*(float)width()/(float)columns);
	
// 	computeDisplayList() ;		
	paintGL() ;
	
	emit zoomChanged(p) ;
}

void VoxelGLDrawer::openFile(const QString f) {

	delete valuesAtPoint ;
	slice = false ;
	slice_pos = 0 ;
	slice_dir = X_SLICE_DIRECTION;
	
	BinaryVoxelDataReader reader(f) ;
	if(pbar)
		connect(&reader, SIGNAL(progressed(int)), pbar, SLOT(setValue(int)));
	
	rows = reader.rows() ;
	columns = reader.columns() ;
	strips = reader.strips() ;
	
	sz=rows*columns*strips ;
	colour.resize(sz) ;
	
// 	rpos.resize(sz*3) ;
// 	for(size_t i = 0 ; i < rpos.size() ; i++)
// 	{
// 		rpos[i] = (0.25/(float)columns)*(float)random()/(float)(RAND_MAX + 1) ;
// 	}
	
	restriction.resize(sz) ; 
	restriction = true ;
	
	valuesAtPoint = new std::vector< std::valarray<quint8> >(*reader.data()) ;
	
	size_x = 1 ;
	size_y = (double)columns/(double)rows ;
	size_z = (double)strips/(double)rows ;

	recalculate() ;
}

void VoxelGLDrawer::initializeGL() {
	
	//get the extended GL functions
	//glBlendEquationSeparate = (void (*)(GLenum, GLfloat))glXGetProcAddressARB((const GLubyte*)"glBlendEquationSeparate") ;
	//glPointParameterfARB = (void (*)(GLenum, GLfloat))glXGetProcAddressARB((const GLubyte*)"glPointParameterfARB") ;
	//glPointParameterfvARB = (void (*)(GLenum, GLfloat *))glXGetProcAddressARB((const GLubyte*)"glPointParameterfvARB") ;
	//glActiveTextureARB = (void (*)(GLuint))glXGetProcAddressARB((const GLubyte*)"glActiveTextureARB") ;
	//glTexEnvf = (void (*)(GLenum, GLenum, GLboolean))glXGetProcAddressARB((const GLubyte*)"glTexEnvf") ;
	
	glViewport(0, 0, 600, 600) ;
	LoadGLTextures(&texture[0]) ;
	glShadeModel(GL_SMOOTH);   // Enables Smooth Shading
	glEnable(GL_LINE_SMOOTH) ;
// 	glEnable(GL_POINT_SMOOTH) ;
	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST) ;
// 	glEnable(GL_POLYGON_SMOOTH) ;
// 	glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST) ;
	
// 	glDepthMask(GL_TRUE);
	
	glEnable(GL_COLOR_MATERIAL) ;
	float mat_specular[] = {0.01f, 0.01f, 0.01f, 1.0f};
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
	glColorMaterial(GL_AMBIENT_AND_DIFFUSE, GL_FRONT_AND_BACK) ;
	glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);
	
	glEnable(GL_NORMALIZE) ;
	glClearColor(0.0f,0.0f,0.0f,1.0f);                                      // Black Background
	glClearDepth(1.0f);                                                     // Depth Buffer Setup
	glEnable(GL_DEPTH_TEST);     
// 	glEnable(GL_ALPHA_TEST); 
	
	
	glEnable(GL_LIGHTING) ;
	glEnable(GL_LIGHT0) ;
	
	GLfloat light_ambient[] = { 0.2, 0.2, 0.2, 1.0 };
	GLfloat light_diffuse[] = { 0.6, 0.6, 0.6, 1.0 };
	GLfloat light_specular[] = { .8, .8, .8, 1.0 };
	GLfloat light_position[] = { 0, 0, zpos, 0 };
	
	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	
// 	glEnable(GL_TEXTURE_2D);                                                // Enable Texture Mapping
// 	glBindTexture(GL_TEXTURE_2D,texture[0]) ;
	
	
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); 
// 	glBlendFunc(GL_ONE_MINUS_DST_ALPHA,GL_DST_ALPHA); 
	glEnable(GL_BLEND);
	glPointSize(1.0f*(float)width()/(float)columns);
	
// 	glEnable( GL_POINT_SPRITE_ARB );
// 	glPointParameterfARB( GL_POINT_SIZE_MIN_ARB, .85*(float)width()/(float)columns );
// 	glPointParameterfARB( GL_POINT_SIZE_MAX_ARB, 32.f);
	
// 	glTexEnvf(GL_POINT_SPRITE_ARB, GL_COORD_REPLACE_ARB, GL_TRUE);
	
	computeDisplayList() ;
	setZoom(100) ;	
}

QSize VoxelGLDrawer::minimumSizeHint() const{
	return QSize(50, 50);
}

QSize VoxelGLDrawer::sizeHint() const{
	return QSize(600, 600);
}

void VoxelGLDrawer::setSegmentDown(int v){
	if (v < 0)
		v = 0 ;
	if (v > m_segmentUp)
		v = m_segmentUp-1 ;
	
	m_segmentDown = v ;
	updateRestrictions() ;
	computeDisplayList() ;
	paintGL() ;
	emit segmentDownChanged(v) ;
	
}

void VoxelGLDrawer::setSegmentUp(int v){
	if (v < m_segmentDown)
		v = m_segmentDown+1 ;
	if (v > 255)
		v = 255 ;
	m_segmentUp = v ;
	
	updateRestrictions() ;
	computeDisplayList() ;
	paintGL() ;
	emit segmentUpChanged(v) ;
// 	
}

void VoxelGLDrawer::setField(int f){
	
	if(f < 0)
	{
		f = 0 ;
		emit fieldChanged(f) ;
		return ;
		
	}
	if (f > (int)valuesAtPoint->size()-1)
	{
		f = valuesAtPoint->size()-1 ;
		emit fieldChanged(f) ;
		return ;
	}
	m_currentField = f ;
	m_min = (*valuesAtPoint)[m_currentField].min() ;
	m_max = (*valuesAtPoint)[m_currentField].max() ;
	computeDisplayList() ;
	paintGL() ;
	emit fieldChanged(f) ;
}

bool VoxelGLDrawer::isInRange(int i) const {
// 	std::cout << "m_segmentDown = " << m_segmentDown << "  m_segmentUp = " << m_segmentUp << " val = " << (int)(*valuesAtPoint)[m_currentField][i] << std::endl ;
	if(m_segmentDown == 0 && m_segmentUp == 255)
		return true ;
	
	return (double)((*valuesAtPoint)[m_currentField][i]-m_min)/(double)(m_max-m_min) >= (double)m_segmentDown/255. && 
		(double)((*valuesAtPoint)[m_currentField][i]-m_min)/(double)(m_max-m_min) <= (double)m_segmentUp/255. ;
}

void VoxelGLDrawer::initializePalette(){
	
// 	colors[][] =
// 	{
// 		{1.0f, 0.5f, 0.5f}, {1.0f, 0.75f, 0.5f}, {1.0f, 1.0f, 0.5f}, {0.75f, 1.0f, 0.5f},
// 		{0.5f, 1.0f, 0.5f}, {0.5f, 1.0f, 0.75f}, {0.5f, 1.0f, 1.0f}, {0.5f, 0.75f, 1.0f},
// 		{0.5f, 0.5f, 1.0f}, {0.75f, 0.5f, 1.0f}, {1.0f, 0.5f, 1.0f}, {1.0f, 0.5f, 0.75f}
// 	};
	
// 	for(size_t i = 0 ; i < 256 ; i++)
// 		palette.push_back(RGBA2i(i,i,i,255)) ;
// 	
// 	for(size_t i = 0 ; i < 127 ; i++)
// 		palette[i] = (RGBA2i(i,i,255,255)) ;
	
    int min_val = (*valuesAtPoint)[m_currentField].min() ;
    int max_val = (*valuesAtPoint)[m_currentField].min() ;
	for(int i = 0 ; i < 256 ; i++)
	{
		size_t r; 
		size_t g; 
		size_t b; 
		HSVtoRGB( &r, &g, &b, 180.*(double)i/256., 1, 1 )  ;
		palette.push_back(RGBA2i(r,g,b,255)) ;
	}
	
// 	for(size_t i = 0 ; i < 127 ; i++)
// 		palette[i] = (RGBA2i(i,i,255,255)) ;
	
//  	palette.push_back(RGBA2i(255,127,127,255)) ;
//  	palette.push_back(RGBA2i(255,191,127,255)) ;
//  	palette.push_back(RGBA2i(255,255,127,255)) ;
//  	palette.push_back(RGBA2i(191,255,127,255)) ;
//  	palette.push_back(RGBA2i(127,255,127,255)) ;
//  	palette.push_back(RGBA2i(127,255,191,255)) ;
//  	palette.push_back(RGBA2i(127,255,255,255)) ;
//  	palette.push_back(RGBA2i(127,191,127,255)) ;
//  	palette.push_back(RGBA2i(127,127,255,255)) ;
//  	palette.push_back(RGBA2i(191,127,255,255)) ;
//  	palette.push_back(RGBA2i(255,127,255,255)) ;
//  	palette.push_back(RGBA2i(255,127,191,255)) ;
}

void VoxelGLDrawer::displayPoints(const std::valarray<size_t> & index, int offset, int mult)
{
// // 	delete (vertex0+offset) ;
// // 	delete (color0+offset) ;
// // 	delete (normal0+offset) ;
	std::valarray<float>  tmpvert(3*index.size());
	
	std::valarray<quint8>  tmpcol(4*index.size()) ;
	
	std::valarray<qint8>  tmpnorm(3*index.size()) ;
	
	for(size_t i = 0; i < index.size() ; i++)
	{
		size_t c[4] ={ i2RGBA(palette[colour[index[i]]]) };
		tmpcol[i*4  ] = c[0];
		tmpcol[i*4+1] = c[1];
		tmpcol[i*4+2] = c[2];
		tmpcol[i*4+3] = alpha;
		
		tmpnorm[i*3  ] = normal[index[i]*3  ]*mult;
		tmpnorm[i*3+1] = normal[index[i]*3+1]*mult;
		tmpnorm[i*3+2] = normal[index[i]*3+2]*mult;
		
		float randx = ((float)rand()/RAND_MAX*2.-1.)*0.2/rows ;
		float randy = ((float)rand()/RAND_MAX*2.-1.)*0.2/columns ;
		float randz = ((float)rand()/RAND_MAX*2.-1.)*0.2/strips ;
		std::valarray<size_t> xyz = toArrayPos(index[i]) ;
		tmpvert[i*3  ] = size_x*(float)(xyz[0]+1)/(float)rows-0.5f+randx;
		tmpvert[i*3+1] = size_y*(float)(xyz[1]+1)/(float)columns-size_y*0.5f+randy;
		tmpvert[i*3+2] = size_z*(float)(xyz[2]+1)/(float)strips-size_z*0.5f+randz;
	}
	glNewList(displayList+ offset, GL_COMPILE) ;
	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_NORMAL_ARRAY);
	glEnableClientState(GL_COLOR_ARRAY);
	glColorPointer(4, GL_UNSIGNED_BYTE, 0, &tmpcol[0]);
	glVertexPointer(3, GL_FLOAT, 0, &tmpvert[0]);
	glNormalPointer(GL_BYTE,0, &tmpnorm[0]);
	glDrawArrays(GL_POINTS, 0, index.size() );
	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);
	glDisableClientState(GL_COLOR_ARRAY);
	glEndList();

}

void VoxelGLDrawer::computeDisplayList() {
// 	if(alpha == 255)
// 		glEnable(GL_POINT_SMOOTH) ;
// 	else
// 		glDisable(GL_POINT_SMOOTH) ;
	
	if(slice)
	{
		computeSlice() ;
		return;
	}
	
	 if(sz == 0)
		 return ;
	 
	
	glDeleteLists (displayList , 10 );
	displayList = glGenLists(10) ;

	std::vector< size_t > zordered ;

	for(int k = 0; k < rows ; k++)
	{
		for(int l = 0; l < columns ; l++)
		{
			for(int m = 0 ; m < strips ; m++)
			{
				int curr = toIndex(k,l,m) ;
				if(restriction[curr] && isInRange(curr))
				{
					zordered.push_back(curr) ;
				}
			}
		}
	}
	int l = 0 ;
	std::valarray<size_t> index((size_t)0, zordered.size()) ;
	for(std::vector< size_t >::const_iterator k = zordered.begin() ; k != zordered.end() ; ++k)
		index[l++] = *k ;
	
	displayPoints(index, 0, 1) ;
// 	std::reverse(&index[0], &index[index.size()]);
	
	
	zordered.clear() ;

	for(int k = rows-1; k >= 0 ; k--)
	{
		for(int l = 0; l < columns ; l++)
		{
			for(int m = 0 ; m < strips ; m++)
			{
				int curr = toIndex(k,l,m) ;
				if(restriction[curr] && isInRange(curr))
				{
					zordered.push_back(curr) ;
				}
			}
		}
	}
	l = 0 ;
	for(std::vector< size_t >::const_iterator k = zordered.begin() ; k != zordered.end() ; ++k)
		index[l++] = *k ;
	
	displayPoints(index, 1, 1) ;
	
	zordered.clear() ;
	for(int m = 0 ; m < strips ; m++)
	{
		for(int l = 0; l < columns ; l++)
		{
			for(int k = 0; k < rows ; k++)
			{
				int curr = toIndex(k,l,m) ;
				if(restriction[curr] && isInRange(curr))
				{
					zordered.push_back(curr) ;
				}
			}
		}
	}
	
	l = 0 ;
	for(std::vector< size_t >::const_iterator k = zordered.begin() ; k != zordered.end() ; ++k)
		index[l++] = *k ;
	displayPoints(index, 2, 1) ;
	
	zordered.clear() ;
	for(int m = strips-1 ; m >=0 ; m--)
	{
		for(int l = 0; l < columns ; l++)
		{
			for(int k = 0; k < rows ; k++)
			{
				int curr = toIndex(k,l,m) ;
				if(restriction[curr] && isInRange(curr))
				{
					zordered.push_back(curr) ;
				}
			}
		}
	}
	
	l = 0 ;
	for(std::vector< size_t >::const_iterator k = zordered.begin() ; k != zordered.end() ; ++k)
		index[l++] = *k ;
	displayPoints(index, 3, 1) ;
	
	
	zordered.clear() ;
	for(int m = strips-1 ; m >=0 ; m--)
	{
		for(int l = columns-1; l >= 0 ; l--)
		{
			for(int k = rows-1; k >=0 ; k--)
			{
				int curr = toIndex(k,l,m) ;
				if(restriction[curr] && isInRange(curr))
				{
					zordered.push_back(curr) ;
				}
			}
		}
	}
	
	l = 0 ;
	for(std::vector< size_t >::const_iterator k = zordered.begin() ; k != zordered.end() ; ++k)
		index[l++] = *k ;
	std::reverse(&index[0], &index[index.size()]);
	displayPoints(index, 4, 1) ;
	
	zordered.clear() ;
	for(int m = 0 ; m <strips ; m++)
	{
		for(int l = 0; l < columns ; l++)
		{
			for(int k = 0; k < rows ; k++)
			{
				int curr = toIndex(k,l,m) ;
				if(restriction[curr] && isInRange(curr))
				{
					zordered.push_back(curr) ;
				}
			}
		}
	}
	
	l = 0 ;
	for(std::vector< size_t >::const_iterator k = zordered.begin() ; k != zordered.end() ; ++k)
		index[l++] = *k ;
	std::reverse(&index[0], &index[index.size()]);
	displayPoints(index, 5, 1) ;
	
	zordered.clear() ;
	for(int l = columns-1; l >= 0 ; l--)
	{
		for(int m = 0 ; m < strips ; m++)
		{
			for(int k = rows-1; k >= 0 ; k--)
			{
				int curr = toIndex(k,l,m) ;
				if(restriction[curr] && isInRange(curr))
				{
					zordered.push_back(curr) ;
				}
			}
		}
	}
	l = 0 ;
	for(std::vector< size_t >::const_iterator k = zordered.begin() ; k != zordered.end() ; ++k)
		index[l++] = *k ;
	displayPoints(index, 6, 1) ;
	
	zordered.clear() ;
	for(int l = columns-1; l >= 0 ; l--)
	{
		for(int m = strips-1 ; m >=0  ; m--)
		{
			for(int k = rows-1 ; k >=0  ; k--)
			{
				int curr = toIndex(k,l,m) ;
				if(restriction[curr] && isInRange(curr))
				{
					zordered.push_back(curr) ;
				}
			}
		}
	}
	l = 0 ;
	for(std::vector< size_t >::const_iterator k = zordered.begin() ; k != zordered.end() ; ++k)
		index[l++] = *k ;
	displayPoints(index, 7, 1) ;
	
	zordered.clear() ;
	for(int l = columns-1; l >= 0 ; l--)
	{
		for(int m = 0 ; m < strips ; m++)
		{
			for(int k = rows-1; k >= 0 ; k--)
			{
				int curr = toIndex(k,l,m) ;
				if(restriction[curr] && isInRange(curr))
				{
					zordered.push_back(curr) ;
				}
			}
		}
	}
	l = 0 ;
	for(std::vector< size_t >::const_iterator k = zordered.begin() ; k != zordered.end() ; ++k)
		index[l++] = *k ;
	std::reverse(&index[0], &index[index.size()]);
	displayPoints(index, 8, 1) ;
	
	zordered.clear() ;
	for(int l = columns-1; l >= 0 ; l--)
	{
		for(int m = strips-1 ; m >=0  ; m--)
		{
			for(int k = rows-1 ; k >=0  ; k--)
			{
				int curr = toIndex(k,l,m) ;
				if(restriction[curr] && isInRange(curr))
				{
					zordered.push_back(curr) ;
				}
			}
		}
	}
	l = 0 ;
	for(std::vector< size_t >::const_iterator k = zordered.begin() ; k != zordered.end() ; ++k)
		index[l++] = *k ;
	std::reverse(&index[0], &index[index.size()]);
	displayPoints(index, 9, 1) ;
	
	
}

quint8 VoxelGLDrawer::valueAverage(int delta, const size_t &x, const size_t &y, const size_t &z) const {

// 	int v = 1 ; //(*valuesAtPoint)[m_currentField][index] ;
	
	std::cout << "vai" << std::endl ;
	
	if(!valuesAtPoint || (*valuesAtPoint)[m_currentField].size() == 0)
		return 0;
	
	
	double value = 0 ;
	size_t num = 0 ;

	for(int i = std::max((int)x-delta,0) ; i < (int)x+delta+1 && i <(int)rows  ; i++)
	{
		std::cout << i << std::endl;
		for(int j = std::max((int)y-delta,0) ; j < (int)y+delta+1 && j <(int)columns; j++)
		{

			std::cout << j << std::endl;
			for(int k = std::max((int)z-delta,0) ; k < (int)z+delta+1 && k <(int)strips  ; k++)
			{
				std::cout << k << std::endl;
				
				std::cout << toIndex(i,j,k) << std::endl ;
				std::cout << (*valuesAtPoint)[m_currentField].size() << std::endl ;
				value += (*valuesAtPoint)[m_currentField][toIndex(i,j,k)] ;
				num++ ;
			}
		}
	}
	std::cout << "vao" << std::endl ;
	return (quint8)round(value/num) ;
}

std::valarray<qint8> VoxelGLDrawer::normalFromSurroundings(const size_t &x, const size_t &y, const size_t &z) const {
	
	std::valarray<qint8> n((qint8)0, 3) ;

	size_t index = toIndex(x,y,z) ;
// 	int v = 1 ; //(*valuesAtPoint)[m_currentField][index] ;
		
	if (isInRange(index))
	{
		for(int i = (int)x-3 ; i < (int)x+4 ; i++)
		{
			if(i == (int)x)
				i++ ;
			
			for(int j = (int)y-3 ; j < (int)y+4 ; j++)
			{
				if(j == (int)y)
					j++ ;
				
				for(int k = (int)z-3 ; k < (int)z+4; k++)
				{
					if(k == (int)z)
						k++ ;
					
					int x_ = 0 ;
					if( (int)x > i )
						x_ = -1 ;
					else if( (int)x < i )
						x_ = 1 ;
					
					int y_ = 0 ;
					if( (int)y>j )
						y_ = -1 ;
					else if( (int)y < j )
						y_ = 1 ;
					
					int z_ = 0 ;
					if( (int)z > k )
						z_ = -1 ;
					else if( (int)z < k )
						z_ = 1 ;
					
	// 				float n_ = 1.f/sqrt(x*x + y*y + k*k) ;
					int val ;
					if(i >= 0 && i < (int)rows    &&
					j >= 0 && j < (int)columns &&
					k >= 0 && k < (int)strips 
					)
					{
						val =  !isInRange(toIndex(i,j,k)); //!restriction[toIndex(i,j,k)];
					}
					else
						val = 1 ;
					
	// 				int v_ = std::abs(val - v) ;
					
	// 				if(v_ != 0)
	// 					v_ = 2 ;
					
					n[0] = (qint8)(x_*val+(int)n[0]);
					n[1] = (qint8)(y_*val+(int)n[1]) ;
					n[2] = (qint8)(z_*val+(int)n[2]) ;
				}
			}
		}
	}
	
	if(x == 0 || y == 0 || z == 0 || x == rows-1 || y == columns-1 || z == strips-1)
		n = 0 ;
	
	if(x == 0 && n[0] == 0)
		n[0] = -25 ;
	if(x == rows-1 && n[0] == 0)
		n[0] = 25 ;
	
	if(y == 0 && n[1] == 0)
		n[1] = -25 ;
	if(y == columns-1 && n[1] == 0)
		n[1] = 25 ;
	
	if(z == 0 && n[2] == 0)
		n[2] = -25 ;
	if(z == strips-1 && n[2] == 0)
		n[2] = 25 ;
	
	return n ;
}

bool VoxelGLDrawer::singlePixel(const size_t& x,  const size_t& y, const size_t& z) const {
	size_t index = toIndex(x, y, z) ;
	
	if(!isInRange(index))
		return false ;
	
	size_t count = 0 ;
	

	for(int i = std::max((int)(x)-3, 0) ; i < (int)std::min(x+4, rows) ; i++)
	{
		if(i == (int)x && i < (int)rows-1)
			i++ ;
		for(int j = std::max((int)(y)-3,0) ; j < (int)std::min(y+4, columns) ; j++)
		{
			if(j == (int)y && j < (int)columns-1)
				j++ ;
			for(int k = std::max((int)(z)-3, 0) ; k < (int)std::min(z+4, strips) ; k++)
			{
				if(k == (int)z && k < (int)strips-1)
					k++ ;
					
				if(isInRange(toIndex(i,j,k)))
					count++ ;
				
			}
		}
	}

	
	return count < 10 ;

}

void VoxelGLDrawer::computeDisplayList(size_t x_0, size_t x_1, size_t y_0, size_t y_1, size_t z_0, size_t z_1) {
	
	glEnable(GL_TEXTURE_2D) ;
	glNewList(displayList, GL_COMPILE) ;
	glBegin(GL_POINTS) ;
	
	for(size_t i =x_0 ; i < std::min(x_1, rows) ; i++)
	{
		for(size_t j =y_0 ; j < std::min(y_1,columns) ; j++)
		{
			for(size_t k =z_0 ; k < std::min(z_1, strips) ; k++)
			{
				if(restriction[toIndex(i, j, k)] && isInRange(toIndex(i, j, k)))
				{
					if(slice)
					{
						unsigned char c[4] ={ i2RGBA(colour[toIndex(i, j, k)]) };
						glColor4ub(c[0], c[1], c[2], (unsigned char)(round((double)c[3]*4))) ;
					}
					else
					{
						glColor4ubv((const GLubyte*)(&colour[toIndex(i, j, k)])) ;
					}
					
					glVertex3f((float)(i+1)/(float)(rows)   -0.5f /*+ rpos[ toIndex(i, j, k) *3]*/, 
					           (float)(j+1)/(float)(columns)-0.5f /*+ rpos[ toIndex(i, j, k) *3+1]*/, 
					           (float)(k+1)/(float)(strips) -0.5f /*+ rpos[ toIndex(i, j, k) *3+2]*/) ;
				}
			}

		}
	}
	glEnd() ;
	glEndList() ;
	glDisable(GL_TEXTURE_2D) ;
}

void VoxelGLDrawer::LoadGLTextures(GLuint  * texture) {
// 	GLfloat tex[64] = 
// 	{ 		.8,  .5,  .6,  .8,  .8,  .6, .5, .8,
// 			.5,  .6, .65,  .7,  .7,  .65,  .6, .5,
// 			.6, .65, .7,  .8, .8, .6,  .65, .6,
// 			.8, .7, .8, 1.,  1.,  .8, .7, .8,
// 			.8, .7, .8, 1.,  1.,  .8, .7, .8,
// 			.6, .65, .7,  .8, .8, .7,  .65, .6,
// 			.5,  .5, .65,  .7,  .7,  .65,  .5, .5,
// 			.8,  .5,  .6,  .8,  .8, .6,  .4,  .8 
// 	} ;
	
// 	GLfloat tex[36] = 
// 	{ 		
// 		1.-0.97, 1.-0.90,  1.-0.83,  1.-0.83,  1.-0.90,  1.-0.97, 
// 		1.-0.90, 1.-0.71,  1.-0.502, 1.-0.502, 1.-0.71,  1.-0.90, 
// 		1.-0.83, 1.-0.502, 1.-0.18,  1.-0.18,  1.-0.502, 1.-0.83, 
// 		1.-0.83, 1.-0.502, 1.-0.18,  1.-0.18,  1.-0.502, 1.-0.83, 
// 		1.-0.90, 1.-0.71,  1.-0.502, 1.-0.502, 1.-0.71,  1.-0.90, 
// 		1.-0.97, 1.-0.90,  1.-0.83,  1.-0.83,  1.-0.90,  1.-0.97, 
// 	} ;
	
	GLfloat tex[16] = 
	{ 		
		.9, 1., 1., .9, 
		1., 1., 1., 1., 
		1., 1., 1., 1., 
		.9, 1., 1., .9, 
	} ;
	
	glGenTextures(1, &texture[0]);                  // Create One Textures
	
	glBindTexture(GL_TEXTURE_2D, texture[0]);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);	
	glTexImage2D(GL_TEXTURE_2D, 0, 1, 4, 4, 0, GL_LUMINANCE, GL_FLOAT, tex) ;
	
// 	GLfloat tex2[2] = {1.0, 1.0};
// 	glTexImage2D(GL_TEXTURE_2D, 0, 2, 1, 1, 0, GL_LUMINANCE_ALPHA, GL_FLOAT, tex2) ;
}

void VoxelGLDrawer::setRestriction(size_t i) {
	restriction[i] = false ;
}

VoxelGLDrawer::VoxelGLDrawer(const size_t r, const size_t c, const size_t s, trsFunc tc , QMainWindow *parent ) : QGLWidget(parent), sz(r*c*s), restriction(true, r*c*s) {
	
	valuesAtPoint = nullptr ;
	
	mousePosOnLeftClick = QPoint(0,0);
	mousePosOnRightClick = QPoint(0,0);
	
	rightDown = false;
	leftDown = false;
	moving = false ;
	
	xangle = 20 ;
	yangle = 20 ;
	zangle = 0 ;
	
	xtransleft = 0;
	ytransleft  = 0;
	
	zpos = 2.3 ;
	
	slice_dir = X_SLICE_DIRECTION;
	slice_pos = 0;
	slice = false;
	initializePalette() ;
	rows = r ;
	columns  = c ;
	strips  = s ;
	toColour = tc ;
	
	size_x = 1 ;
	size_y = (double)columns/(double)rows ;
	size_z = (double)strips/(double)rows ;
	
	if(parent)
	{
		pbar = new QProgressBar(parent->statusBar()) ;
		parent->statusBar()->addPermanentWidget(pbar) ;
	}
	else
		pbar = nullptr ;
	
// 	rpos.resize(sz*3) ;
// 	for(size_t i = 0 ; i < rpos.size() ; i++)
// 	{
// 		rpos[i] = (0.25/(float)columns)*(float)random()/(float)(RAND_MAX + 1) ;
// 	}
	
	
	
	colour.resize(sz) ;
	normal.resize(sz*3) ;
	alpha = 255 ;
	
	displayList = 1 ; //glGenLists(1);
	currentDisplayList = -1 ;
	
	start_offset = 0 ;
	
	m_segmentDown  = 0 ;
	m_segmentUp  = 255 ;
	m_currentField = 0 ;
	m_zoom = 100 ;

	recalculate() ;
	updateRestrictions() ;

	m_min = 0;//(*valuesAtPoint)[m_currentField].min() ;
	m_max = 0;//(*valuesAtPoint)[m_currentField].max() ;
	


}

VoxelGLDrawer::VoxelGLDrawer(QString f, QMainWindow *mainwin) : QGLWidget(mainwin) {
	
	valuesAtPoint = nullptr ;
	
	if(mainwin)
	{
		pbar = new QProgressBar(mainwin->statusBar()) ;
		mainwin->statusBar()->addPermanentWidget(pbar) ;
	}
	else 
		pbar = nullptr ;
	
	mousePosOnLeftClick = QPoint(0,0);
	mousePosOnRightClick = QPoint(0,0);
	
	rightDown = false;
	leftDown = false;
	moving = false ;
	slice = false;
	
	xangle = 20 ;
	yangle = 20 ;
	zangle = 0 ;
	
	xtransleft = 0;
	ytransleft  = 0;
	
	zpos = 2.3 ;
	
	toColour = phaseInfo ;

	openFile(f) ;
	
	for(size_t i = 0 ; i < valuesAtPoint->size() ;i++)
	{
		delta.push_back((*valuesAtPoint)[i].max() - (*valuesAtPoint)[i].min()) ;
		min.push_back((*valuesAtPoint)[i].min()) ;
	}
	
	start_offset = 0 ;
	
	m_segmentDown  = 0 ;
	m_segmentUp  = 255 ;
	m_currentField = 0 ;
	
	displayList = 1 ;
	currentDisplayList = -1 ;
	normal.resize(sz*3) ;
	initializePalette() ;
	m_min = (*valuesAtPoint)[m_currentField].min() ;
	m_max = (*valuesAtPoint)[m_currentField].max() ;
	updateRestrictions() ;
	computeDisplayList() ;
	
	m_zoom = 85 ;
	setZoom(85) ;
	
	alpha = 255 ;


	size_x = 1 ;
	size_y = (double)columns/(double)rows ;
	size_z = (double)strips/(double)rows ;

}


VoxelGLDrawer::VoxelGLDrawer(QMainWindow *parent) : QGLWidget(parent) {
	
	if(parent)
	{
		pbar = new QProgressBar(parent->statusBar()) ;
		parent->statusBar()->addPermanentWidget(pbar) ;
	}
	else
		pbar = nullptr ;
	
	valuesAtPoint = new std::vector< std::valarray<quint8> >(0) ;
	
	mousePosOnLeftClick = QPoint(0,0);
	mousePosOnRightClick = QPoint(0,0);
	
	rightDown = false;
	leftDown = false;
	moving = false ;
	slice = false;
	
	xangle = 20 ;
	yangle = 20 ;
	zangle = 0 ;
	
	xtransleft = 0;
	ytransleft  = 0;
	
	rows = 0 ;
	columns = 0;
	strips = 0 ;
	sz = 0 ;
	
	zpos = 2.3 ;
	
	toColour = phaseInfo ;
	
	start_offset = 0 ;
	
	colour.resize(0) ;
	
	alpha = 255 ;

	
	m_segmentDown  = 0 ;
	m_segmentUp  = 255 ;
	m_currentField = 0 ;
	
	displayList = 1 ; 
	currentDisplayList = -1 ;
	
	delta.push_back(1) ; 
	min.push_back(0) ;
	initializePalette() ;
	
	m_zoom = 100 ;
	
	m_min = 0 ;//(*valuesAtPoint)[m_currentField].min() ;
	m_max = 0 ;//(*valuesAtPoint)[m_currentField].max() ;
	
	size_x = 1 ;
	size_y = (double)columns/(double)rows ;
	size_z = (double)strips/(double)rows ;
	

	
}

VoxelGLDrawer::~VoxelGLDrawer()  {
	if(dynamic_cast<QMainWindow *>(parent()))
		dynamic_cast<QMainWindow *>(parent())->statusBar()->removeWidget(pbar) ;

	
	glDeleteLists(displayList, 8) ;
	makeCurrent();
}

void VoxelGLDrawer::computeSlice() {
	if(!slice)
	{
		computeDisplayList() ;
		return ;
	}
	
	slice_pos++ ; 
	
	if(slice_dir == X_SLICE_DIRECTION)
	{
		slice_pos%=rows ;
		computeDisplayList(slice_pos, slice_pos+rows/16, 0, columns, 0, strips) ;		
	}
	if(slice_dir == Y_SLICE_DIRECTION)
	{
		slice_pos%=columns ;
		computeDisplayList(0, rows, slice_pos, slice_pos+columns/16, 0, strips) ;		
		
	}
	if(slice_dir == Z_SLICE_DIRECTION)
	{
		slice_pos%=strips ;
		computeDisplayList(0, rows, 0, columns, slice_pos, slice_pos+strips/16) ;		
	}
}

void VoxelGLDrawer::recalculate() {
	
	
	if(valuesAtPoint->empty() && valuesAtPoint[0].size() != 0)
		return ;
		
// 	std::vector<std::valarray<quint8> > vals(*valuesAtPoint) ;
// 		
// 	for(size_t i = 0 ; i < rows ; i++)
// 	{
// 		for(size_t j = 0 ; j< columns ; j++)
// 		{
// 			for(size_t k = 0 ; k< strips ; k++)
// 			{
// 				vals[0][toIndex(i,j,k)] = valueAverage(1,i,j,k) ;
// 			}
// 		}
// 	}
	
	toColour(valuesAtPoint, &colour, &restriction, rows, columns, strips) ;
}

size_t VoxelGLDrawer::toIndex(const size_t i, const size_t j , const size_t k) const {
	return i + j*columns + k*rows*columns ;
}

std::valarray<size_t> VoxelGLDrawer::toArrayPos( const size_t where) const {
	std::valarray<size_t> ret(3) ;
	size_t planeSize = rows*columns ;
	
	ret[2] = where/planeSize ;
	ret[1] = (where%planeSize)/columns ;
	ret[0] = ((where%planeSize)%columns)%rows ;
	
	return ret ;
}

float VoxelGLDrawer::getXtrans() const {
	return  -(float)xtransleft/(float)width();
}

float VoxelGLDrawer::getYtrans() const {
	return -(float)ytransleft/(float)height();
}

float VoxelGLDrawer::getXAngle() const {
	return (float)xangle ;
}

float VoxelGLDrawer::getYAngle() const {
	return (float)yangle ;
}

float VoxelGLDrawer::getZAngle() const {
	return (float)zangle ;
}

void inSphereColour(const std::vector< std::valarray<qint8> > *d, std::valarray<quint8> * c, std::valarray<bool> * res , const size_t r, const size_t co, const size_t s) {
	for(size_t i = 0 ; i< c->size() ; i++)
	{
		std::valarray<float> pos(3) ;
		size_t planeSize = s*co ;
		
		pos[2] = (float)(i/planeSize)/(float)s -0.5f;
		pos[1] = (float)((i%planeSize)/co)/(float)co -0.5f;
		pos[0] = (float)(((i%planeSize)%co)%r)/(float)r -0.5f;
		
		if((pos[0]+0.2)*(pos[0]+0.2) + (pos[1]+0.2)*(pos[1]+0.2)+ (pos[2])*(pos[2]) < 0.05 ||
		   (pos[0]-0.2)*(pos[0]-0.2) + (pos[1]-0.2)*(pos[1]-0.2)+ (pos[2]+0.2)*(pos[2]+0.2) < 0.02 
		  )
		{
			(*c)[i] = 8 ; 
		}
		else if((pos[0]+0.2)*(pos[0]+0.2) + (pos[1]-0.2)*(pos[1]-0.2)+ (pos[2]-0.3)*(pos[2]-0.3) < 0.015)
		{
			(*c)[i] = 0 ; 
		}
		else
		{
			(*c)[i] = 5 ; 
		}
	}
	*res = true ;
} 

std::valarray<bool> inSphere(const std::vector< std::valarray<qint8> > *d, const size_t r, const size_t c , const size_t s) {
	std::valarray<bool> ret(true, r*c*s) ;
	
	
	for(size_t i = 0 ; i< r*c*s ; i++)
	{
		std::valarray<float> pos(3) ;
		size_t planeSize = s*c ;
		
		pos[2] = (float)(i/planeSize)/(float)s -0.5f;
		pos[1] = (float)((i%planeSize)/c)/(float)c -0.5f;
		pos[0] = (float)(((i%planeSize)%c)%r)/(float)r -0.5f;
		
		
		if((pos[0]+0.2)*(pos[0]+0.2) + (pos[1]+0.2)*(pos[1]+0.2)+ (pos[2])*(pos[2]) < 0.05 ||
		   (pos[0]-0.2)*(pos[0]-0.2) + (pos[1]-0.2)*(pos[1]-0.2)+ (pos[2]+0.2)*(pos[2]+0.2) < 0.02 
		  )
		{
			ret[i] = true ; 
		}
		else if((pos[0]+0.2)*(pos[0]+0.2) + (pos[1]-0.2)*(pos[1]-0.2)+ (pos[2]-0.3)*(pos[2]-0.3) < 0.015)
		{
			ret[i] = true ; 
		}
		else
		{
			ret[i] = false ; 
		}
	}
	
	return ret ;
	
}

std::valarray<bool> outOfSphere(const std::vector< std::valarray<quint8> > *d, const size_t r, const size_t c, const size_t s) {
	std::valarray<bool> ret(true, r*c*s) ;
	
	for(size_t i = 0 ; i < r*c*s ; i++)
	{
		
		std::valarray<float> pos(3) ;
		size_t planeSize = s*c ;
		
		pos[2] = (float)(i/planeSize)/(float)s -0.5f;
		pos[1] = (float)((i%planeSize)/c)/(float)c -0.5f;
		pos[0] = (float)(((i%planeSize)%c)%r)/(float)r -0.5f;
		
		if((pos[0]+0.2)*(pos[0]+0.2) + (pos[1]+0.2)*(pos[1]+0.2)+ (pos[2])*(pos[2]) < 0.05 ||
		   (pos[0]-0.2)*(pos[0]-0.2) + (pos[1]-0.2)*(pos[1]-0.2)+ (pos[2]+0.2)*(pos[2]+0.2) < 0.02 
		  )
		{
			ret[i] = false ; 
		}
		else if((pos[0]+0.2)*(pos[0]+0.2) + (pos[1]-0.2)*(pos[1]-0.2)+ (pos[2]-0.3)*(pos[2]-0.3) < 0.015)
		{
			ret[i] = false ; 
		}
		else
		{
			ret[i] = true ; 
		}
	}
	
	return ret ;
	
}

void phaseInfo (const std::vector< std::valarray<quint8> > *d, std::valarray<quint8> * c, std::valarray<bool> * res, const size_t r, const size_t co, const size_t s){
	
	quint8 min = (*d)[0].min() ;
	quint8 max = (*d)[0].max() ;
	
	
	for(size_t i = 0 ; i< c->size() ; i++)
	{
		(*c)[i] = (quint8)round(255.*(double)((*d)[0][i]-min)/(double)(max-min)) ; (*d)[0][i] ;
// 		if(abs((double)(*c)[i]- 0  ) <=1 || 
// 			abs((double)(*c)[i] - 20 ) <=1 || 
// 			abs((double)(*c)[i] - 40 ) <=1 ||
// 			abs((double)(*c)[i] - 60 ) <=1 ||
// 			abs((double)(*c)[i] - 80 ) <=1 ||
// 			abs((double)(*c)[i] - 100) <=1 ||
// 			abs((double)(*c)[i] - 120) <=1 ||
// 			abs((double)(*c)[i] - 140) <=1 ||
// 			abs((double)(*c)[i] - 160) <=1 ||
// 			abs((double)(*c)[i] - 180) <=1 ||
// 			abs((double)(*c)[i] - 200) <=1 ||
// 			abs((double)(*c)[i] - 220) <=1 ||
// 			abs((double)(*c)[i] - 240) <=1
// 		)
			(*res)[i] = true ;
// 		else
// 			(*res)[i] = false ;
		
	}
// 	for(size_t i = 0 ; i< c->size() ; i++)
// 	{
// // 		std::valarray<float> pos(3) ;
// // 		size_t planeSize = s*co ;
// // 		
// // 		pos[2] = (float)(i/planeSize)/(float)s -0.5f;
// // 		pos[1] = (float)((i%planeSize)/co)/(float)co -0.5f;
// // 		pos[0] = (float)(((i%planeSize)%co)%r)/(float)r -0.5f;
// // 		
// 		if((*d)[0][i] == 0 )
// 		{
// 			(*c)[i] = 8 ; 
// // 			(*res)[i] = false ;
// 		}
// 		else if((*d)[0][i] == 1 )
// 		{
// 			(*c)[i] = 7 ; 
// // 			(*res)[i] = false ;
// 		}
// 		else if ((*d)[0][i] == 2 )
// 		{
// 			(*c)[i] = 2 ; 
// // 			(*res)[i] = false ;
// 		}
// 		else if ((*d)[0][i] == 3 )
// 		{
// 			(*c)[i] = 5 ; 
// // 			(*res)[i] = false ;
// 		}
// 		else if ((*d)[0][i] == 4 )
// 		{
// 			(*c)[i] = 1 ; 
// // 			(*res)[i] = false ;
// 		}
// 		else if ((*d)[0][i] == 5 )
// 		{
// 			(*c)[i] = 3 ; 
// // 			(*res)[i] = false ;
// 		}
// 		else if ((*d)[0][i] == 6 )
// 		{
// 			(*c)[i] = 11 ; 
// // 			(*res)[i] = false ;
// 		}
// 		else
// 		{
// 			(*c)[i] = 10 ; 
// // 			(*res)[i] = false ;
// 		}
// 	}
}
