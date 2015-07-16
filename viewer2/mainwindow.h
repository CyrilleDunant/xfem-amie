#ifndef __MAINWINDOW_H__
#define __MAINWINDOW_H__

#include <QMainWindow>
#include <QToolButton>
#include <QSlider>
#include <QAction>
#include <QString>
#include <QMenuBar>
#include <QStatusBar>
#include <QHBoxLayout>
#include <QFileDialog>
#include <QToolBar>
#include <QSpinBox>
#include <QCheckBox>
#include <string>

#include "voxelGlDrawer.h"
#include "triangleGlDrawer.h"
#include "triangleDataReader.h"
#include "buffer.h"



class MainWindow : public QMainWindow
{
	Q_OBJECT
		
public:
	MainWindow();
	
private:
	QSlider *createSlider();
	
	VoxelGLDrawer *voxeldisplay;
	TriangleGLDrawer *triangledisplay;
	QSlider *zoomSlider;
	
	QAction *openAct ;
	QAction *exitAct ;
	
	QMenu *fileMenu;
	QToolBar *fileToolBar;
	QSpinBox * zoom ;
	QSpinBox * alpha ;
	QSpinBox * field ;
	QSpinBox * layer ;

	QSlider * time ;
	QSlider * downSlider ;
	QSlider * upSlider ;

	QCheckBox * meshDisplay ;

	QStatusBar *statusbar ;
	QToolButton * printButton ;
	
	Buffer buffer ;

	void createActions() ;
	void createMenus() ;
	void createToolBars() ;
	void createStatusBar() ;
	
	bool triangles(const QString s) const ;
	bool multi(const QString s) const ;
	bool voxels(const QString s) const ;

	void disconnect(TriangleGLDrawer * display) ;
	
public slots:
	void open() ;
	void open(const QString &) ;
	void getFile(int) ;
	void getLayer(int) ;

signals:
	void prepareFile(const QString &) ;
};

#endif
