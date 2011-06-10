#ifndef __Triangle_DRAWER__H__
#define __Triangle_DRAWER__H__

#include <QPushButton>
#include <QWidget>
#include <QMouseEvent>
#include <QKeyEvent>
#include <QTime>
#include <QString>
#include <QtOpenGL/QGLWidget>

#include <valarray>
#include <map>
#include <vector>
#include <string>
#include <iostream>

#include "triangleDataReader.h"

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

#include <GL/glx.h>

#include <GL/glxext.h>
#include <GL/glext.h>

#define i2RGBA( __val_32__ )   (size_t)__val_32__  & 0x000000FF , ((size_t)__val_32__ & 0x0000FF00)>>8 , ((size_t)__val_32__ & 0x00FF0000)>>16 , ((size_t)__val_32__ & 0xFF000000)>>24 

#define RGBA2i( __r__, __g__, __b__ , __a__)  ((size_t)__r__ +  (size_t)((size_t) __g__ << 8)  + (size_t)((size_t) __b__ << 16)  +  (size_t)((size_t) __a__ << 24))



class TriangleGLDrawer : public QGLWidget
{
	Q_OBJECT        // must include this if you use Qt signals/slots

		
public slots:
	void setXTranslation(int)  ; //
	void setYTranslation(int)  ; //
	void setZoom(int)  ; //
	void openFile(const QString f)  ; //
	void openFile(TriangleDataReader * f)  ; //
	void setSet(int set)  ; //
	void grab() ;
	
signals:
	void xTranslationChanged(int)  ; //
	void yTranslationChanged(int)  ; //
	void zoomChanged(int)  ; //
	void setChanged(int) ;
	
	
protected:
	void ( * APIENTRY glBlendEquationSeparate )(GLenum, GLfloat);
	void ( * APIENTRY glPointParameterfARB )(GLenum, GLfloat);
	void ( * APIENTRY glPointParameterfvARB)(GLenum, GLfloat *) ;
	void ( * APIENTRY glTexEnvf )(GLenum, GLenum, GLboolean) ;
	void ( * APIENTRY glActiveTextureARB)(GLuint) ;
	
	
protected:
	QPoint mousePosOnLeftClick ;
	
	bool moving ;
	bool leftDown ;
		
	float zpos ;
	int zoom ;
	

	float getXtrans() const ; //
	float getYtrans() const ; //
	
	size_t currentDisplayList ;
	int currentSet ;
	int fracup ; 
	int fracdown ;
	
	int scale ;
	
	float max_x ;
	float max_y ;
	float min_x ;
	float min_y ;
	
	bool minmaxinit ;
	
public slots:
	void setSegmentDown(int) ;
	void setSegmentUp(int) ;
	void setScale(int) ;
	void setTimePlane(int) ;
	
signals:
	void segmentUpChanged(int) ;
	void segmentDownChanged(int) ;
	void scaleChanged(int) ;
	void timePlaneChanged(int) ;
	
protected:
	std::vector< std::valarray<float> > * valuesAtPoint ;
	TriangleDataReader * reader ;
	
	quint64 numberOfTriangles ;
	quint64 numberOfPointsPerTriangle ;
	quint64 numberOfExtraFields ;
	quint64 numberOfExtraTimePlanes ;
	quint64 currentTimePlane ;
	
	std::vector< std::valarray<size_t> > colour ;

	void computeDisplayList() ; //
	GLint displayList ;
	
public:
	
	int xtransleft ;
	int ytransleft ;
	
	TriangleGLDrawer(std::vector<std::valarray<float> > * v, int np, int set, QWidget * parent = 0 ) ;
	TriangleGLDrawer(TriangleDataReader * f, int set, const std::vector<std::pair<float, float> > & limits, QWidget *parent = 0) ; //
	TriangleGLDrawer(QString f, const std::vector<std::pair<float, float> > & limits, QWidget *parent = 0) ; //
	TriangleGLDrawer(QString f, int set, const std::vector<std::pair<float, float> > & limits, QWidget *parent = 0) ; //
	TriangleGLDrawer(QWidget *parent = 0) ; //
	~TriangleGLDrawer() ; //
	
	QString fileName ;
	
	size_t numberOfElements() const { return numberOfTriangles ; } //
	std::vector<std::pair<float, float> > limits ;

	std::vector< std::valarray<float> > * getValuesAtPoint() { return valuesAtPoint ; }
	
protected:
	GLuint texture[1] ;
	
protected:
	void initializeGL() ; //
	void resizeGL(int w, int h) ; //
	void paintGL() ; //
	QSize minimumSizeHint() const ;
	QSize sizeHint() const ;
	void mousePressEvent(QMouseEvent *event) ; //
	void mouseMoveEvent(QMouseEvent *event) ; //
	void wheelEvent(QWheelEvent * event ) ;
	void mouseReleaseEvent(QMouseEvent *event) ; //
	void keyPressEvent ( QKeyEvent * event ) ; //
	void computeDisplay( ) const ; //
	
	void HSVtoRGB( size_t *r, size_t *g, size_t *b, float h, float s, float v ) const  ;
	
};




#endif 
