#include <QApplication>
// #define nullptr NULL ;


/*#include "voxelGlDrawer.h"*/
#include "mainwindow.h"
#include "triangleGlDrawer.h"

int main(int argc, char **argv) {
	QApplication app(argc, argv);

	
// 	QWidget window ;
// 	
// 	QPushButton quit("Quit !", &window);
// 	
// 	window.resize(200, 120);
// 	quit.resize(90,40) ;
// 	quit.setGeometry((window.width()-quit.width())>>1,
// 	                 (window.height()-quit.height())>>1, 
// 	                 quit.width(), 
// 	                 quit.height());
// 	quit.setFont(QFont("Bitstream Vera Sans", 24, QFont::Bold)) ;
// 	
// 	QObject::connect(&quit, SIGNAL(clicked()), &app, SLOT(quit())) ;
// 
// 	window.show();
	MainWindow vld ;
// 	
	if(app.arguments().size() > 1)
		vld.open(app.arguments().at(1));

// 	VoxelGLDrawer vld(128, 128, 128) ;
	
	vld.show() ;
	
// 	TriangleGLDrawer t("results") ;
// 	t.show() ;
	return app.exec();
}
