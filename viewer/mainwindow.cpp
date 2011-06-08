#include "mainwindow.h"

MainWindow::MainWindow()
{
	voxeldisplay = new VoxelGLDrawer() ;
	triangledisplay = NULL ;
	setCentralWidget(voxeldisplay);
	
	createActions() ;
	createMenus() ;
	createToolBars() ;
	createStatusBar() ;
	
// 	zoomSlider = createSlider();
// 	
// 	
// 	QHBoxLayout *mainLayout = new QHBoxLayout;
// 	mainLayout->addWidget(voxeldisplay);
// 	mainLayout->addWidget(zoomSlider);
// 	centralWidget()->setLayout(mainLayout);
	
	setWindowTitle("empty");
}

 void MainWindow::createStatusBar()
 {
	 statusbar = statusBar() ;
     statusBar()->showMessage(tr("Ready"));
 }


QSlider *MainWindow::createSlider()
{
	QSlider *slider = new QSlider(Qt::Vertical);
	slider->setRange(1, 400);
	slider->setSingleStep(2);
	slider->setPageStep(10);
	slider->setTickInterval(40);
	slider->setTickPosition(QSlider::TicksRight);
	return slider;
}

void MainWindow::createActions()
{
	
	openAct = new QAction(QIcon::fromTheme("document-open"), tr("&Open..."), this);
	openAct->setShortcut(tr("Ctrl+O"));
	openAct->setStatusTip(tr("Open an existing file"));
	connect(openAct, SIGNAL(triggered()), this, SLOT(open()));
	
	exitAct = new QAction(tr("&Quit"), this);
	exitAct->setShortcut(tr("Ctrl+Q"));
	exitAct->setStatusTip(tr("Exit the application"));
	connect(exitAct, SIGNAL(triggered()), this, SLOT(close()));
}

void MainWindow::createMenus()
{
	fileMenu = menuBar()->addMenu(tr("&File"));
	fileMenu->addAction(openAct);
	fileMenu->addSeparator();
	fileMenu->addAction(exitAct);
	
}

void MainWindow::createToolBars()
{
	fileToolBar = addToolBar(tr("File"));
	fileToolBar->addAction(openAct);
	fileToolBar->addSeparator() ;
	
	printButton = new QToolButton(fileToolBar) ;
	connect(printButton, SIGNAL(released()), voxeldisplay, SLOT(grab()));
	fileToolBar->addWidget(printButton) ;
	printButton->setIcon(QIcon::fromTheme("camera-photo")) ;
	
	zoom  = new QSpinBox(fileToolBar) ;
	zoom->setRange ( 1, 9600 ) ;
	zoom->setValue ( 100 ) ;
	zoom->setSuffix(" %") ;
	connect(zoom, SIGNAL(valueChanged(int)), voxeldisplay, SLOT(setZoom(int)));
	connect(voxeldisplay, SIGNAL(zoomChanged(int)), zoom, SLOT(setValue(int)));
	fileToolBar->addWidget(zoom);
	fileToolBar->addSeparator() ;
	
	alpha  = new QSpinBox(fileToolBar) ;
	alpha->setRange ( 0, 255 ) ;
	alpha->setValue ( 0 ) ;
	connect(alpha, SIGNAL(valueChanged(int)), voxeldisplay, SLOT(setAlpha(int)));
	connect(voxeldisplay, SIGNAL(alphaChanged(int)), alpha, SLOT(setValue(int)));
	fileToolBar->addWidget(alpha);
	
	field  = new QSpinBox(fileToolBar) ;
	field->setRange ( 0, 255 ) ;
	field->setValue ( 0 ) ;
	connect(field, SIGNAL(valueChanged(int)), voxeldisplay, SLOT(setField(int)));
	connect(voxeldisplay, SIGNAL(fieldChanged(int)), field, SLOT(setValue(int)));
	fileToolBar->addWidget(field);
	
	time  = new QSpinBox(fileToolBar) ;
	time->setRange ( 0, 255 ) ;
	time->setValue ( 0 ) ;
//	connect(field, SIGNAL(valueChanged(int)), voxeldisplay, SLOT(setField(int)));
//	connect(voxeldisplay, SIGNAL(fieldChanged(int)), field, SLOT(setValue(int)));
	fileToolBar->addWidget(time);
	
	downSlider = new QSlider(Qt::Horizontal,fileToolBar);
	downSlider->setRange(0, 254);
	downSlider->setSingleStep(5);
	downSlider->setPageStep(10);
	downSlider->setTickInterval(1);
	downSlider->setValue(0) ;
	downSlider->setTracking ( false ) ;
// 	downSlider->setTickPosition(QSlider::TicksDown);
	connect(downSlider, SIGNAL(valueChanged(int)), voxeldisplay, SLOT(setSegmentDown(int)));
	connect(voxeldisplay, SIGNAL(segmentDownChanged(int)), downSlider, SLOT(setValue(int)));
	fileToolBar->addWidget(downSlider);
	
	upSlider = new QSlider(Qt::Horizontal,fileToolBar);
	upSlider->setRange(1, 255);
	upSlider->setSingleStep(1);
	upSlider->setPageStep(10);
	upSlider->setTickInterval(1);
	upSlider->setValue(255) ;
	upSlider->setTracking ( false ) ;
// 	upSlider->setTickPosition(QSlider::TicksDown);
	connect(upSlider, SIGNAL(valueChanged(int)), voxeldisplay, SLOT(setSegmentUp(int)));
	connect(voxeldisplay, SIGNAL(segmentUpChanged(int)), upSlider, SLOT(setValue(int)));
	fileToolBar->addWidget(upSlider);
	
	
// 	fileToolBar->addAction(exitAct);
}

void MainWindow::open()
{
	QString fileName = QFileDialog::getOpenFileName(this);
	if (!fileName.isEmpty())
	{
		if(triangles(fileName))
		{
			delete voxeldisplay ;
			std::vector<std::pair<float, float> > limits ;
			int xpos = 0 ;
			int ypos = 0 ;
			int zoomval = zoom->value(); 
			if(triangledisplay)
			{
			  limits = triangledisplay->limits ;
				xpos = triangledisplay->xtransleft ;
				ypos = triangledisplay->ytransleft ;
				
			}
			delete triangledisplay ;
			voxeldisplay = NULL ;
			triangledisplay = new TriangleGLDrawer(fileName, limits) ;
			triangledisplay->xtransleft = xpos ;
			triangledisplay->ytransleft = ypos ;
			
			setCentralWidget(triangledisplay);
			QFileInfo pathInfo( fileName );
			setWindowTitle(pathInfo.fileName());
			triangledisplay->fileName = fileName ;
			connect(zoom, SIGNAL(valueChanged(int)), triangledisplay, SLOT(setZoom(int)));
			connect(triangledisplay, SIGNAL(zoomChanged(int)), zoom, SLOT(setValue(int)));
			
			connect(alpha, SIGNAL(valueChanged(int)), triangledisplay, SLOT(setSet(int)));
			connect(triangledisplay, SIGNAL(setChanged(int)), alpha, SLOT(setValue(int)));
			
			connect(field, SIGNAL(valueChanged(int)), triangledisplay, SLOT(setScale(int)));
			connect(triangledisplay, SIGNAL(scaleChanged(int)), field, SLOT(setValue(int)));
			
//			connect(time, SIGNAL(valueChanged(int)), triangledisplay, SLOT(setTimePlane(int)));
//			connect(triangledisplay, SIGNAL(timePlaneChanged(int)), time, SLOT(setValue(int)));
			
			connect(downSlider, SIGNAL(valueChanged(int)), triangledisplay, SLOT(setSegmentDown(int)));
			connect(triangledisplay, SIGNAL(segmentDownChanged(int)), downSlider, SLOT(setValue(int)));
			
			connect(upSlider, SIGNAL(valueChanged(int)), triangledisplay, SLOT(setSegmentUp(int)));
			connect(triangledisplay, SIGNAL(segmentUpChanged(int)), upSlider, SLOT(setValue(int)));
			
			connect(printButton, SIGNAL(released()), triangledisplay, SLOT(grab()));
			triangledisplay->setZoom(zoomval) ;
			
			alpha->setValue(0) ;
			alpha->setRange(-1, 1000);
			
			downSlider->setRange(0, 9999);
			downSlider->setValue(0) ;
			
			field->setValue(1) ;
			field->setRange(-1, 1000);
			
			time->setValue(0) ;
			time->setRange(-1, 1000) ;

			upSlider->setRange(1, 10000);
			upSlider->setValue(10000) ;
			
		}
		else if(voxels(fileName))
		{
			delete voxeldisplay ;
			delete triangledisplay ;
			statusBar()->clearMessage() ;
			voxeldisplay = new VoxelGLDrawer(fileName, this) ;
			triangledisplay = NULL ;
			setCentralWidget(voxeldisplay);
			QFileInfo pathInfo( fileName );
			setWindowTitle(pathInfo.fileName());
			
			connect(printButton, SIGNAL(released()), voxeldisplay, SLOT(grab()));
			
			connect(zoom, SIGNAL(valueChanged(int)), voxeldisplay, SLOT(setZoom(int)));
			connect(voxeldisplay, SIGNAL(zoomChanged(int)), zoom, SLOT(setValue(int)));
			
			connect(alpha, SIGNAL(valueChanged(int)), voxeldisplay, SLOT(setAlpha(int)));
			connect(voxeldisplay, SIGNAL(alphaChanged(int)), alpha, SLOT(setValue(int)));
			
			connect(field, SIGNAL(valueChanged(int)), voxeldisplay, SLOT(setField(int)));
			connect(voxeldisplay, SIGNAL(fieldChanged(int)), field, SLOT(setValue(int)));
			
			connect(downSlider, SIGNAL(valueChanged(int)), voxeldisplay, SLOT(setSegmentDown(int)));
			connect(voxeldisplay, SIGNAL(segmentDownChanged(int)), downSlider, SLOT(setValue(int)));
			
			connect(upSlider, SIGNAL(valueChanged(int)), voxeldisplay, SLOT(setSegmentUp(int)));
			connect(voxeldisplay, SIGNAL(segmentUpChanged(int)), upSlider, SLOT(setValue(int)));
			
			alpha->setValue(255) ;
		}
		else if(multi(fileName))
		{
		    QFile list(fileName) ;
		    if (!list.open(QIODevice::ReadOnly | QIODevice::Text))
		    {
			    return ;
		    }
		    QTextStream stream ;
		    stream.setDevice(&list);
		    QString buff ;
		    stream >> buff ; // should be "MULTI"
		    bool go_on = true ;
		    while(go_on)
		    {
			stream >> buff ;
			if(buff.length() > 0)
			{
			    buff = fileName.section('/',0,-2) + "/" + buff ;
			    files << buff ;
			}
			else
			    go_on = false ;
		    }
		    list.close();

		    connect(time, SIGNAL(valueChanged(int)), this, SLOT(getFile(int))) ;
		    connect(this, SIGNAL(prepareFile(QString)), this, SLOT(open(QString))) ;

		    open(files.at(0)) ;
		}
	}
}

void MainWindow::open(const QString &fileName)
{
	if (!fileName.isEmpty())
	{
	    if(triangles(fileName))
		{
		    delete voxeldisplay ;
		    std::vector<std::pair<float, float> > limits ;
		    int xpos = 0 ;
		    int ypos = 0 ;
		    int zoomval = zoom->value();
		    if(triangledisplay)
		    {
		      limits = triangledisplay->limits ;
			    xpos = triangledisplay->xtransleft ;
			    ypos = triangledisplay->ytransleft ;

		    }
		    delete triangledisplay ;
		    voxeldisplay = NULL ;
		    triangledisplay = new TriangleGLDrawer(fileName, alpha->value(), limits) ;
		    triangledisplay->xtransleft = xpos ;
		    triangledisplay->ytransleft = ypos ;

		    setCentralWidget(triangledisplay);
		    QFileInfo pathInfo( fileName );
		    setWindowTitle(pathInfo.fileName());
		    triangledisplay->fileName = fileName ;
		    connect(zoom, SIGNAL(valueChanged(int)), triangledisplay, SLOT(setZoom(int)));
		    connect(triangledisplay, SIGNAL(zoomChanged(int)), zoom, SLOT(setValue(int)));

		    connect(alpha, SIGNAL(valueChanged(int)), triangledisplay, SLOT(setSet(int)));
		    connect(triangledisplay, SIGNAL(setChanged(int)), alpha, SLOT(setValue(int)));

		    connect(field, SIGNAL(valueChanged(int)), triangledisplay, SLOT(setScale(int)));
		    connect(triangledisplay, SIGNAL(scaleChanged(int)), field, SLOT(setValue(int)));

//			connect(time, SIGNAL(valueChanged(int)), triangledisplay, SLOT(setTimePlane(int)));
//			connect(triangledisplay, SIGNAL(timePlaneChanged(int)), time, SLOT(setValue(int)));

		    connect(downSlider, SIGNAL(valueChanged(int)), triangledisplay, SLOT(setSegmentDown(int)));
		    connect(triangledisplay, SIGNAL(segmentDownChanged(int)), downSlider, SLOT(setValue(int)));

		    connect(upSlider, SIGNAL(valueChanged(int)), triangledisplay, SLOT(setSegmentUp(int)));
		    connect(triangledisplay, SIGNAL(segmentUpChanged(int)), upSlider, SLOT(setValue(int)));

		    connect(printButton, SIGNAL(released()), triangledisplay, SLOT(grab()));
		    triangledisplay->setZoom(zoomval) ;

//		    alpha->setValue(0) ;
		    alpha->setRange(-1, 1000);

		    downSlider->setRange(0, 9999);
		    downSlider->setValue(0) ;

//		    field->setValue(1) ;
//		    field->setRange(-1, 1000);

//		    time->setValue(0) ;
//		    time->setRange(-1, 1000) ;

		    upSlider->setRange(1, 10000);
		    upSlider->setValue(10000) ;


		}
		else if(voxels(fileName))
		{
			delete voxeldisplay ;
			delete triangledisplay ;
			voxeldisplay = new VoxelGLDrawer(fileName) ;
			triangledisplay = NULL ;
			setCentralWidget(voxeldisplay);
			QFileInfo pathInfo( fileName );
			setWindowTitle(pathInfo.fileName());
			connect(zoom, SIGNAL(valueChanged(int)), voxeldisplay, SLOT(setZoom(int)));
			connect(voxeldisplay, SIGNAL(zoomChanged(int)), zoom, SLOT(setValue(int)));
			
			connect(alpha, SIGNAL(valueChanged(int)), voxeldisplay, SLOT(setAlpha(int)));
			connect(voxeldisplay, SIGNAL(alphaChanged(int)), alpha, SLOT(setValue(int)));

			connect(field, SIGNAL(valueChanged(int)), voxeldisplay, SLOT(setField(int)));
			connect(voxeldisplay, SIGNAL(fieldChanged(int)), field, SLOT(setValue(int)));
			
			connect(downSlider, SIGNAL(valueChanged(int)), voxeldisplay, SLOT(setSegmentDown(int)));
			connect(voxeldisplay, SIGNAL(segmentDownChanged(int)), downSlider, SLOT(setValue(int)));
			
			connect(upSlider, SIGNAL(valueChanged(int)), voxeldisplay, SLOT(setSegmentUp(int)));
			connect(voxeldisplay, SIGNAL(segmentUpChanged(int)), upSlider, SLOT(setValue(int)));
			
		}
		else if(multi(fileName))
		{
		    QFile list(fileName) ;
		    std::cerr << fileName.toStdString() << std::endl ;
		    if (!list.open(QIODevice::ReadOnly | QIODevice::Text))
		    {
			    return ;
		    }
		    QTextStream stream ;
		    stream.setDevice(&list);
		    QString buff ;
		    stream >> buff ; // should be "MULTI"
		    bool go_on = true ;
		    while(go_on)
		    {
			stream >> buff ;
			if(buff.length() > 0)
			    files << buff ;
			else
			    go_on = false ;
			std::cerr << buff.toStdString() << std::endl ;
		    }
		    list.close();

		    connect(time, SIGNAL(valueChanged(int)), this, SLOT(getFile(int))) ;
		    connect(this, SIGNAL(prepareFile(QString)), this, SLOT(open(QString))) ;

//		    open(files.at(0)) ;
		}
	}
}

void MainWindow::getFile(int i)
{
    int index = i ;
    if(index >= files.size() || index < 0)
    {
	index = 0 ;
	time->setValue(0);
    }
    QString name = files.at(index) ;
    emit prepareFile(name);
}

bool MainWindow::triangles(const QString s) const {
	
	QFile file(s) ;
	
	if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
	{
		return false;
	}
	
	QTextStream stream;
	stream.setDevice(&file);
	QString type ;
	stream >> type ;
	file.close();

	std::cerr << s.toStdString() << std::endl ;
	
	return (type == "TRIANGLES" || type == "BIN_TRIANGLES") ;
}

bool MainWindow::voxels(const QString s) const {
	
	QFile file(s) ;
	
	if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
	{
		return false;
	}
	
	QTextStream stream;
	stream.setDevice(&file);
	QString type ;
	stream >> type ;
	file.close();
	
	return (type == "VOXELS") ;
}

bool MainWindow::multi(const QString s) const {

	QFile file(s) ;

	if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
	{
		return false;
	}

	QTextStream stream;
	stream.setDevice(&file);
	QString type ;
	stream >> type ;
	file.close();

	return (type == "MULTI") ;
}
