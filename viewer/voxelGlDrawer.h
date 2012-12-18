#define nullptr NULL ;



#ifndef __VOXEL_DRAWER__H__
#define __VOXEL_DRAWER__H__

#include <QPushButton>
#include <QWidget>
#include <QMouseEvent>
#include <QKeyEvent>
#include <QTime>
#include <QString>
#include <QWheelEvent>
#include <QtOpenGL/QGLWidget>

#include <QMainWindow>
#include <QProgressBar>
#include <QStatusBar>

#include <valarray>
#include <map>
#include <vector>
#include <string>
#include <iostream>

#include "asciiVoxelDataReader.h"
#include "binaryVoxelDataReader.h"

#include <GL/gl.h>
#include <GL/glu.h>

// #define glXGetProcAddress glXGetProcAddressARB
#include <GL/glx.h>

#include <GL/glxext.h>
#include <GL/glext.h>

#define i2RGBA( __val_32__ )  (unsigned char) ((size_t)__val_32__  & 0x000000FF ),(unsigned char) (((size_t)__val_32__ & 0x0000FF00)>>8 ),(unsigned char)( ((size_t)__val_32__ & 0x00FF0000)>>16 ), (unsigned char)(((size_t)__val_32__ & 0xFF000000)>>24)

#define RGBA2i( __r__, __g__, __b__ , __a__)  ((size_t)__r__ +  (size_t)((size_t) __g__ << 8)  + (size_t)((size_t) __b__ << 16)  +  (size_t)((size_t) __a__ << 24))


typedef void (*trsFunc)(const std::vector< std::valarray<quint8> > *, std::valarray<quint8> *, std::valarray<bool> * res, const size_t, const size_t, const size_t) ;

typedef std::valarray<bool> (*restrictionFunction)(const std::vector< std::valarray<float> > *) ;

void inSphereColour(const std::vector< std::valarray<quint8> > *d, std::valarray<quint8> * c, std::valarray<bool> * res , const size_t, const size_t, const size_t) ;

void phaseInfo (const std::vector< std::valarray<quint8> > *d, std::valarray<quint8> * c , std::valarray<bool> * res, const size_t r, const size_t co, const size_t s);

std::valarray<bool> inSphere(const std::vector< std::valarray<qint8> > *d, const size_t, const size_t, const size_t) ;

std::valarray<bool> outOfSphere(const std::vector< std::valarray<qint8> > *d, const size_t, const size_t, const size_t) ;


typedef enum
{
	X_SLICE_DIRECTION,
	Y_SLICE_DIRECTION,
	Z_SLICE_DIRECTION
} SliceDirection ;

typedef enum
{
	X_SHIFT_DIRECTION,
	Y_SHIFT_DIRECTION,
	Z_SHIFT_DIRECTION
} ShiftDirection ;

class VoxelGLDrawer : public QGLWidget
{
	Q_OBJECT        // must include this if you use Qt signals/slots

		
	QProgressBar * pbar ;
	
public slots:
	void setXAngle(int)  ; //
	void setYAngle(int)  ; //
	void setZAngle(int)  ; //
	void setXTranslation(int)  ; //
	void setYTranslation(int)  ; //
	void setZoom(int)  ; //
	void setAlpha(int)  ; //
	void openFile(const QString f)  ; //
	void setSegmentDown(int) ;
	void setSegmentUp(int) ;
	void setField(int) ;
	void grab() ;
	
signals:
	void xAngleChanged(int)  ; //
	void yAngleChanged(int)  ; //
	void zAngleChanged(int)  ; //
	void xTranslationChanged(int)  ; //
	void yTranslationChanged(int)  ; //
	void zoomChanged(int)  ; //
	void alphaChanged(int)  ; //
	void segmentUpChanged(int) ;
	void segmentDownChanged(int) ;
	void fieldChanged(int) ;
	void progressed(int) ;
	
	
protected:
	void ( * APIENTRY glBlendEquationSeparate )(GLenum, GLfloat);
	void ( * APIENTRY glPointParameterfARB )(GLenum, GLfloat);
	void ( * APIENTRY glPointParameterfvARB)(GLenum, GLfloat *) ;
	void ( * APIENTRY glTexEnvf )(GLenum, GLenum, GLboolean) ;
	void ( * APIENTRY glActiveTextureARB)(GLuint) ;
	void displayPoints(const std::valarray<size_t> & index, int offset, int mult) ;

	
protected:
	QPoint mousePosOnLeftClick ;
	QPoint mousePosOnRightClick ;
	
	bool moving ;
	bool rightDown ;
	bool leftDown ;
		
	float zpos ;
	
	int xtransleft ;
	int ytransleft ;
	float getXtrans() const ; //
	float getYtrans() const ; //
		
	int xangle ;
	int yangle ;
	int zangle ;
	int alpha ;
	
	float getXAngle() const ; //
	float getYAngle() const ; //
	float getZAngle() const ; //
	
	SliceDirection slice_dir ;
	size_t slice_pos ;
	bool slice ;
	size_t start_offset ;
	
	int m_segmentDown ;
	int m_segmentUp ;
	int m_currentField ;
	int m_zoom ;
	quint8 m_min ;
	quint8 m_max ;
	
	bool isInRange(int i) const ;
	std::vector<float> delta ;
	std::vector<float> min ;
	
	double size_x ;
	double size_y ;
	double size_z ;
	
protected:
	size_t sz ;
	std::valarray<bool> restriction ;
	std::vector< std::valarray<quint8> > * valuesAtPoint ;
	std::valarray<quint8> colour ;
	std::valarray<qint8> normal ;
	std::vector<size_t> palette ;
	size_t rows ;
	size_t columns ;
	size_t strips ;
	
	void HSVtoRGB( size_t *r, size_t *g, size_t *b, float h, float s, float v ) const {
	int i;
	float f, p, q, t;
	if( s == 0 ) {
                // achromatic (grey)
		*r = *g = *b = static_cast<size_t>(v*255);
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
		*r = static_cast<size_t>(v*255.);
		*g = static_cast<size_t>(t*255.);
		*b = static_cast<size_t>(p*255.);
		break;
	case 1:
		*r = static_cast<size_t>(q*255.);
		*g = static_cast<size_t>(v*255.);
		*b = static_cast<size_t>(p*255.);
		break;
	case 2:
		*r = static_cast<size_t>(p*255.);
		*g = static_cast<size_t>(v*255.);
		*b = static_cast<size_t>(t*255.);
		break;
	case 3:
		*r = static_cast<size_t>(p*255.);
		*g = static_cast<size_t>(q*255.);
		*b = static_cast<size_t>(v*255.);
		break;
	case 4:
		*r = static_cast<size_t>(t*255.);
		*g = static_cast<size_t>(p*255.);
		*b = static_cast<size_t>(v*255.);
		break;
	default:                // case 5:
		*r = static_cast<size_t>(v*255.);
		*g = static_cast<size_t>(p*255.);
		*b = static_cast<size_t>(q*255.);
		break;
	}
}

	
// 	std::valarray<float> rpos ;
	
	float mat[4][4];
	
	trsFunc toColour ;
	
	void recalculate() ; //
	void initializePalette() ;
	
	void updateRestrictions() ;
	void computeSlice() ; //
	void computeDisplayList() ; //
	void computeDisplayList(size_t x_0, size_t x_1, size_t y_0, size_t y_1, size_t z_0, size_t z_1) ; //
	GLint displayList ;
	GLint currentDisplayList ;
	
	bool isVisible(const size_t i, const size_t j, const size_t k) const ;
	
	quint8 valueAverage(int delta, const size_t &x, const size_t &y, const size_t &z) const ;
	
public:
	VoxelGLDrawer(const size_t r, const size_t c, const size_t s, trsFunc tc = inSphereColour, QMainWindow *parent = 0) ; //
	VoxelGLDrawer(QString f, QMainWindow *parent = 0) ; //
	VoxelGLDrawer(QMainWindow *parent = 0) ; //
	~VoxelGLDrawer() ; //
	
	void setRestriction(size_t i) ;
	size_t numberOfElements() const { return sz ; } //
	size_t toIndex(const size_t i, const size_t j , const size_t k) const ; //
	std::valarray<size_t> toArrayPos( const size_t where) const ; //
	
	std::valarray<qint8> normalFromSurroundings(const size_t&, const size_t&, const size_t&) const ;
	bool isBoundary(const size_t& x, const size_t& y, const size_t& z) const ;
	bool singlePixel(const size_t& x,  const size_t& y, const size_t& z) const ;
	
protected:
	GLuint texture[1] ;
	
protected:
	void LoadGLTextures(GLuint  * texture) ; //
	void initializeGL() ; //
	void resizeGL(int w, int h) ; //
	void paintGL() ; //
	QSize minimumSizeHint() const ;
	QSize sizeHint() const ;
	void mousePressEvent(QMouseEvent *event) ; //
	void mouseMoveEvent(QMouseEvent *event) ; //
	void mouseReleaseEvent(QMouseEvent *event) ; //
	void wheelEvent(QWheelEvent * event ) ;
	void keyPressEvent ( QKeyEvent * event ) ; //
	void computeDisplay( ) ; //

};




#endif 
