// #define nullptr NULL ;



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

#include "glutils.h"

/*
   Return a RGB colour value given a scalar v in the range [vmin,vmax]
   In this case each colour component ranges from 0 (no contribution) to
   1 (fully saturated), modifications for other ranges is trivial.
   The colour is clipped at the end of the scales if v is outside
   the range [vmin,vmax]
*/




typedef struct {
    double r,g,b;
} COLOUR;

COLOUR GetColour(double v) ;

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
	void setMeshDisplay(int state) ;
	
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
	bool meshOnly = false ;
		
	float zpos ;
	int zoom ;
	

	float getXtrans() const ; //
	float getYtrans() const ; //
	
	size_t currentDisplayList ;
	size_t wireFrameList ;
	int currentSet ;
	int fracup ; 
	int fracdown ;
	
	
	
	float max_x ;
	float max_y ;
	float min_x ;
	float min_y ;
	
	bool minmaxinit ;
	
	float valUnderCursor ;
	
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
	int scale ;
	
	TriangleGLDrawer(std::vector<std::valarray<float> > * v, int np, int set, const std::vector<std::pair<float, float> > & limits, QWidget * parent = 0 ) ;
	TriangleGLDrawer(TriangleDataReader * f, int set, const std::vector<std::pair<float, float> > & limits, QWidget *parent = 0) ; //
	TriangleGLDrawer(QString f, const std::vector<std::pair<float, float> > & limits, QWidget *parent = 0) ; //
	TriangleGLDrawer(QString f, int set, const std::vector<std::pair<float, float> > & limits, QWidget *parent = 0) ; //
	TriangleGLDrawer(QWidget *parent = 0) ;
	
	void reset(std::vector<std::valarray<float> > * v, int np, int set, const std::vector<std::pair<float, float> > & limits) ;
	void reset(TriangleDataReader * f, int set, const std::vector<std::pair<float, float> > & limits) ;
	void reset(QString f, const std::vector<std::pair<float, float> > & limits) ;
	void reset(QString f, int set, const std::vector<std::pair<float, float> > & limits) ;
	void reset() ;
	
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
	void computeDisplay( ) ; //
	
	void HSVtoRGB( size_t *r, size_t *g, size_t *b, float h, float s, float v ) const  ;
	
};




#endif 
