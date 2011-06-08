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
#include <string>

#include "voxelGlDrawer.h"
#include "triangleGlDrawer.h"

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
	QSpinBox * time ;
	QSlider * downSlider ;
	QSlider * upSlider ;
	QStatusBar *statusbar ;
	QToolButton * printButton ;
	
	QStringList files ;

	void createActions() ;
	void createMenus() ;
	void createToolBars() ;
	void createStatusBar() ;
	
	bool triangles(const QString s) const ;
	bool multi(const QString s) const ;
	bool voxels(const QString s) const ;
	
public slots:
	void open() ;
	void open(const QString &) ;
	void getFile(int) ;

signals:
	void prepareFile(const QString &) ;
};

#endif
