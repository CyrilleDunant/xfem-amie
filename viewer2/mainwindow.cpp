#include "mainwindow.h"
#include "triangleDataReader.h"

MainWindow::MainWindow()
{
	voxeldisplay = new VoxelGLDrawer() ;
	triangledisplay = nullptr ;
	setCentralWidget( voxeldisplay );

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

	setWindowTitle( "empty" );

	connect( time, SIGNAL( valueChanged( int ) ), this, SLOT( getFile( int ) ) ) ;
	connect( this, SIGNAL( prepareFile( QString ) ), this, SLOT( open( QString ) ) ) ;

	connect( layer, SIGNAL( valueChanged( int ) ), this, SLOT( getLayer( int ) ) ) ;

}

void MainWindow::createStatusBar()
{
	statusbar = statusBar() ;
	statusBar()->showMessage( tr( "Ready" ) );
}


QSlider *MainWindow::createSlider()
{
	QSlider *slider = new QSlider( Qt::Vertical );
	slider->setRange( 1, 400 );
	slider->setSingleStep( 2 );
	slider->setPageStep( 10 );
	slider->setTickInterval( 40 );
	slider->setTickPosition( QSlider::TicksRight );
	return slider;
}

void MainWindow::createActions()
{

	openAct = new QAction( QIcon::fromTheme( "document-open" ), tr( "&Open..." ), this );
	openAct->setShortcut( tr( "Ctrl+O" ) );
	openAct->setStatusTip( tr( "Open an existing file" ) );
	connect( openAct, SIGNAL( triggered() ), this, SLOT( open() ) );

	exitAct = new QAction( tr( "&Quit" ), this );
	exitAct->setShortcut( tr( "Ctrl+Q" ) );
	exitAct->setStatusTip( tr( "Exit the application" ) );
	connect( exitAct, SIGNAL( triggered() ), this, SLOT( close() ) );
}

void MainWindow::createMenus()
{
	fileMenu = menuBar()->addMenu( tr( "&File" ) );
	fileMenu->addAction( openAct );
	fileMenu->addSeparator();
	fileMenu->addAction( exitAct );

}

void MainWindow::createToolBars()
{
	fileToolBar = addToolBar( tr( "File" ) );
	fileToolBar->addAction( openAct );
	fileToolBar->addSeparator() ;

	printButton = new QToolButton( fileToolBar ) ;
	connect( printButton, SIGNAL( released() ), voxeldisplay, SLOT( grab() ) );
	fileToolBar->addWidget( printButton ) ;
	printButton->setIcon( QIcon::fromTheme( "camera-photo" ) ) ;

	zoom  = new QSpinBox( fileToolBar ) ;
	zoom->setRange( 1, 9600 ) ;
	zoom->setValue( 100 ) ;
	zoom->setSuffix( " %" ) ;
	connect( zoom, SIGNAL( valueChanged( int ) ), voxeldisplay, SLOT( setZoom( int ) ) );
	connect( voxeldisplay, SIGNAL( zoomChanged( int ) ), zoom, SLOT( setValue( int ) ) );
	fileToolBar->addWidget( zoom );
	fileToolBar->addSeparator() ;

	alpha  = new QSpinBox( fileToolBar ) ;
	alpha->setRange( 0, 255 ) ;
	alpha->setValue( 0 ) ;
	connect( alpha, SIGNAL( valueChanged( int ) ), voxeldisplay, SLOT( setAlpha( int ) ) );
	connect( voxeldisplay, SIGNAL( alphaChanged( int ) ), alpha, SLOT( setValue( int ) ) );
	fileToolBar->addWidget( alpha );

	field  = new QSpinBox( fileToolBar ) ;
	field->setRange( 0, 255 ) ;
	field->setValue( 0 ) ;
	connect( field, SIGNAL( valueChanged( int ) ), voxeldisplay, SLOT( setField( int ) ) );
	connect( voxeldisplay, SIGNAL( fieldChanged( int ) ), field, SLOT( setValue( int ) ) );
	fileToolBar->addWidget( field );

	layer  = new QSpinBox( fileToolBar ) ;
	layer->setRange( 0, 255 ) ;
	layer->setValue( 0 ) ;
//	connect(field, SIGNAL(valueChanged(int)), voxeldisplay, SLOT(setField(int)));
//	connect(voxeldisplay, SIGNAL(fieldChanged(int)), field, SLOT(setValue(int)));
	fileToolBar->addWidget( layer );

	time = new QSlider( Qt::Horizontal, fileToolBar ) ;
	time->setRange( 0, 1 );
	time->setSingleStep( 1 );
	time->setPageStep( 1 );
	time->setTickInterval( 1 );
	time->setValue( 0 ) ;
	time->setTracking( false ) ;
// 	downSlider->setTickPosition(QSlider::TicksDown);
	fileToolBar->addWidget( time );

	downSlider = new QSlider( Qt::Horizontal, fileToolBar );
	downSlider->setRange( 0, 254 );
	downSlider->setSingleStep( 5 );
	downSlider->setPageStep( 10 );
	downSlider->setTickInterval( 1 );
	downSlider->setValue( 0 ) ;
	downSlider->setTracking( false ) ;
// 	downSlider->setTickPosition(QSlider::TicksDown);
	connect( downSlider, SIGNAL( valueChanged( int ) ), voxeldisplay, SLOT( setSegmentDown( int ) ) );
	connect( voxeldisplay, SIGNAL( segmentDownChanged( int ) ), downSlider, SLOT( setValue( int ) ) );
	fileToolBar->addWidget( downSlider );

	upSlider = new QSlider( Qt::Horizontal, fileToolBar );
	upSlider->setRange( 1, 255 );
	upSlider->setSingleStep( 1 );
	upSlider->setPageStep( 10 );
	upSlider->setTickInterval( 1 );
	upSlider->setValue( 255 ) ;
	upSlider->setTracking( false ) ;
// 	upSlider->setTickPosition(QSlider::TicksDown);
	connect( upSlider, SIGNAL( valueChanged( int ) ), voxeldisplay, SLOT( setSegmentUp( int ) ) );
	connect( voxeldisplay, SIGNAL( segmentUpChanged( int ) ), upSlider, SLOT( setValue( int ) ) );
	fileToolBar->addWidget( upSlider );


// 	fileToolBar->addAction(exitAct);
}

void MainWindow::open()
{
	QString fileName = QFileDialog::getOpenFileName( this );

	if( !fileName.isEmpty() )
	{
		if( triangles( fileName ) )
		{
			delete voxeldisplay ;
			std::vector<std::pair<float, float> > limits ;
			int xpos = 0 ;
			int ypos = 0 ;
			int zoomval = zoom->value();
			int scale =  1 ;

			if( triangledisplay )
			{
				scale = triangledisplay->scale ;
				limits = triangledisplay->limits ;
				xpos = triangledisplay->xtransleft ;
				ypos = triangledisplay->ytransleft ;
			}

			voxeldisplay = nullptr ;

			int i = buffer.fileInBuffer( fileName ) ;

			if( i < 0 )
			{
				FileBuffer fbuff( fileName ) ;
				buffer.push_back( fbuff );
				i = buffer.fileInBuffer( fileName ) ;
			}

			if(!triangledisplay)
			{
				if( i >= 0 )
				{
					LayerBuffer lbuff = buffer.value( i ).value( layer->value() ) ;
					triangledisplay = new TriangleGLDrawer( lbuff.values, lbuff.numberOfPointsPerTriangle, alpha->value(), limits ) ;
				}
				else
				{
					triangledisplay = new TriangleGLDrawer( fileName, alpha->value(), limits ) ;
				}
			}
			else
			{
				if( i >= 0 )
				{
					LayerBuffer lbuff = buffer.value( i ).value( layer->value() ) ;
					triangledisplay->reset( lbuff.values, lbuff.numberOfPointsPerTriangle, alpha->value(), limits ) ;
				}
				else
				{
					triangledisplay->reset( fileName, alpha->value(), limits ) ;
				}
			}

			triangledisplay->xtransleft = xpos ;
			triangledisplay->ytransleft = ypos ;
			triangledisplay->scale = scale ;
			
			setCentralWidget( triangledisplay );
			QFileInfo pathInfo( fileName );
			setWindowTitle( pathInfo.fileName() );
			triangledisplay->fileName = fileName ;
                        zoom->setRange( 1, 9600 ) ;
			connect( zoom, SIGNAL( valueChanged( int ) ), triangledisplay, SLOT( setZoom( int ) ) );
			connect( triangledisplay, SIGNAL( zoomChanged( int ) ), zoom, SLOT( setValue( int ) ) );

                     
			connect( alpha, SIGNAL( valueChanged( int ) ), triangledisplay, SLOT( setSet( int ) ) );
			connect( triangledisplay, SIGNAL( setChanged( int ) ), alpha, SLOT( setValue( int ) ) );

			connect( field, SIGNAL( valueChanged( int ) ), triangledisplay, SLOT( setScale( int ) ) );
			connect( triangledisplay, SIGNAL( scaleChanged( int ) ), field, SLOT( setValue( int ) ) );

			connect( downSlider, SIGNAL( valueChanged( int ) ), triangledisplay, SLOT( setSegmentDown( int ) ) );
			connect( triangledisplay, SIGNAL( segmentDownChanged( int ) ), downSlider, SLOT( setValue( int ) ) );

			connect( upSlider, SIGNAL( valueChanged( int ) ), triangledisplay, SLOT( setSegmentUp( int ) ) );
			connect( triangledisplay, SIGNAL( segmentUpChanged( int ) ), upSlider, SLOT( setValue( int ) ) );

			connect( printButton, SIGNAL( released() ), triangledisplay, SLOT( grab() ) );

			field->setValue( 1 ) ;
			field->setRange( 0, 4096 );

			downSlider->setRange( 0, 9999 );
			downSlider->setValue( 0 ) ;

// 			field->setValue( 1 ) ;
			field->setRange( -1, 1000 );

			layer->setValue( 0 ) ;
			layer->setRange( -1, 1000 ) ;

			upSlider->setRange( 1, 10000 );
			upSlider->setValue( 10000 ) ;
			
			triangledisplay->setZoom( zoomval ) ;
			triangledisplay->setScale( scale ); 

		}
		else if( voxels( fileName ) )
		{
			delete voxeldisplay ;
			delete triangledisplay ;
			statusBar()->clearMessage() ;
			voxeldisplay = new VoxelGLDrawer( fileName, this ) ;
			triangledisplay = nullptr ;
			setCentralWidget( voxeldisplay );
			QFileInfo pathInfo( fileName );
			setWindowTitle( pathInfo.fileName() );

			connect( printButton, SIGNAL( released() ), voxeldisplay, SLOT( grab() ) );

			connect( zoom, SIGNAL( valueChanged( int ) ), voxeldisplay, SLOT( setZoom( int ) ) );
			connect( voxeldisplay, SIGNAL( zoomChanged( int ) ), zoom, SLOT( setValue( int ) ) );

			connect( alpha, SIGNAL( valueChanged( int ) ), voxeldisplay, SLOT( setAlpha( int ) ) );
			connect( voxeldisplay, SIGNAL( alphaChanged( int ) ), alpha, SLOT( setValue( int ) ) );

			connect( field, SIGNAL( valueChanged( int ) ), voxeldisplay, SLOT( setField( int ) ) );
			connect( voxeldisplay, SIGNAL( fieldChanged( int ) ), field, SLOT( setValue( int ) ) );

			connect( downSlider, SIGNAL( valueChanged( int ) ), voxeldisplay, SLOT( setSegmentDown( int ) ) );
			connect( voxeldisplay, SIGNAL( segmentDownChanged( int ) ), downSlider, SLOT( setValue( int ) ) );

			connect( upSlider, SIGNAL( valueChanged( int ) ), voxeldisplay, SLOT( setSegmentUp( int ) ) );
			connect( voxeldisplay, SIGNAL( segmentUpChanged( int ) ), upSlider, SLOT( setValue( int ) ) );

			alpha->setValue( 255 ) ;
		}
		else if( multi( fileName ) )
		{
                        zoom->setRange( 1, 9600 ) ;
			connect( zoom, SIGNAL( valueChanged( int ) ), triangledisplay, SLOT( setZoom( int ) ) );
			connect( triangledisplay, SIGNAL( zoomChanged( int ) ), zoom, SLOT( setValue( int ) ) );

			connect( alpha, SIGNAL( valueChanged( int ) ), triangledisplay, SLOT( setSet( int ) ) );
			connect( triangledisplay, SIGNAL( setChanged( int ) ), alpha, SLOT( setValue( int ) ) );

			connect( field, SIGNAL( valueChanged( int ) ), triangledisplay, SLOT( setScale( int ) ) );
			connect( triangledisplay, SIGNAL( scaleChanged( int ) ), field, SLOT( setValue( int ) ) );

			connect( downSlider, SIGNAL( valueChanged( int ) ), triangledisplay, SLOT( setSegmentDown( int ) ) );
			connect( triangledisplay, SIGNAL( segmentDownChanged( int ) ), downSlider, SLOT( setValue( int ) ) );

			connect( upSlider, SIGNAL( valueChanged( int ) ), triangledisplay, SLOT( setSegmentUp( int ) ) );
			connect( triangledisplay, SIGNAL( segmentUpChanged( int ) ), upSlider, SLOT( setValue( int ) ) );

			connect( printButton, SIGNAL( released() ), triangledisplay, SLOT( grab() ) );

			field->setValue( 1 ) ;
			field->setRange( 0, 4096 );

			downSlider->setRange( 0, 9999 );
			downSlider->setValue( 0 ) ;

			layer->setValue( 0 ) ;
			layer->setRange( -1, 1000 ) ;

			upSlider->setRange( 1, 10000 );
			upSlider->setValue( 10000 ) ;

			buffer.clear() ;
			buffer = Buffer( fileName ) ;

			time->setRange( 0, buffer.size() - 1 );

			open( buffer.at( 0 ).file ) ;
		}
	}
}

void MainWindow::open( const QString &fileName )
{
	if( !fileName.isEmpty() )
	{
		if( triangles( fileName ) )
		{
			delete voxeldisplay ;
			std::vector<std::pair<float, float> > limits ;
			int xpos = 0 ;
			int ypos = 0 ;
			int zoomval = zoom->value();
			int scale = 1 ;
			voxeldisplay = nullptr ;

			if( triangledisplay )
			{
				scale = triangledisplay->scale ;
				limits = triangledisplay->limits ;
				xpos = triangledisplay->xtransleft ;
				ypos = triangledisplay->ytransleft ;
			}

			int i = buffer.fileInBuffer( fileName ) ;

			if( i < 0 )
			{
				FileBuffer fbuff( fileName ) ;
				buffer.push_back( fbuff );
				i = buffer.fileInBuffer( fileName ) ;
			}
			if(!triangledisplay)
			{
				if( i >= 0 )
				{
					LayerBuffer lbuff = buffer.value( i ).value( layer->value() ) ;
					triangledisplay = new TriangleGLDrawer( lbuff.values, lbuff.numberOfPointsPerTriangle, alpha->value(), limits ) ;
				}
				else
				{
					triangledisplay = new TriangleGLDrawer( fileName, alpha->value(), limits ) ;
				}
			}
			else
			{
				if( i >= 0 )
				{
					LayerBuffer lbuff = buffer.value( i ).value( layer->value() ) ;
					triangledisplay->reset( lbuff.values, lbuff.numberOfPointsPerTriangle, alpha->value(), limits ) ;
				}
				else
				{
					triangledisplay->reset( fileName, alpha->value(), limits ) ;
				}
			}

			triangledisplay->xtransleft = xpos ;
			triangledisplay->ytransleft = ypos ;
			triangledisplay->scale = scale ;
			
			if( triangledisplay == NULL )
				std::cerr << "nullptr" << std::endl ;

			setCentralWidget( triangledisplay );
			QFileInfo pathInfo( fileName );
			setWindowTitle( pathInfo.fileName() );
			triangledisplay->fileName = fileName ;
                        zoom->setRange( 1, 9600 ) ;
                        field->setRange(0, 4096) ;
			connect( zoom, SIGNAL( valueChanged( int ) ), triangledisplay, SLOT( setZoom( int ) ) );
			connect( triangledisplay, SIGNAL( zoomChanged( int ) ), zoom, SLOT( setValue( int ) ) );

                        field->setRange(0, 4096) ;
			connect( alpha, SIGNAL( valueChanged( int ) ), triangledisplay, SLOT( setSet( int ) ) );
			connect( triangledisplay, SIGNAL( setChanged( int ) ), alpha, SLOT( setValue( int ) ) );

			connect( field, SIGNAL( valueChanged( int ) ), triangledisplay, SLOT( setScale( int ) ) );
			connect( triangledisplay, SIGNAL( scaleChanged( int ) ), field, SLOT( setValue( int ) ) );

			connect( downSlider, SIGNAL( valueChanged( int ) ), triangledisplay, SLOT( setSegmentDown( int ) ) );
			connect( triangledisplay, SIGNAL( segmentDownChanged( int ) ), downSlider, SLOT( setValue( int ) ) );

			connect( upSlider, SIGNAL( valueChanged( int ) ), triangledisplay, SLOT( setSegmentUp( int ) ) );
			connect( triangledisplay, SIGNAL( segmentUpChanged( int ) ), upSlider, SLOT( setValue( int ) ) );

			connect( printButton, SIGNAL( released() ), triangledisplay, SLOT( grab() ) );

                        
			downSlider->setRange( 0, 9999 );
			downSlider->setValue( 0 ) ;

			upSlider->setRange( 1, 10000 );
			upSlider->setValue( 10000 ) ;
			
			triangledisplay->setZoom( zoomval ) ;
			triangledisplay->setScale( scale ) ;

		}
		else if( voxels( fileName ) )
		{
			delete voxeldisplay ;
			delete triangledisplay ;
			voxeldisplay = new VoxelGLDrawer( fileName ) ;
			triangledisplay = NULL ;
			setCentralWidget( voxeldisplay );
			QFileInfo pathInfo( fileName );
			setWindowTitle( pathInfo.fileName() );
			connect( zoom, SIGNAL( valueChanged( int ) ), voxeldisplay, SLOT( setZoom( int ) ) );
			connect( voxeldisplay, SIGNAL( zoomChanged( int ) ), zoom, SLOT( setValue( int ) ) );

			connect( alpha, SIGNAL( valueChanged( int ) ), voxeldisplay, SLOT( setAlpha( int ) ) );
			connect( voxeldisplay, SIGNAL( alphaChanged( int ) ), alpha, SLOT( setValue( int ) ) );

			connect( field, SIGNAL( valueChanged( int ) ), voxeldisplay, SLOT( setField( int ) ) );
			connect( voxeldisplay, SIGNAL( fieldChanged( int ) ), field, SLOT( setValue( int ) ) );

			connect( downSlider, SIGNAL( valueChanged( int ) ), voxeldisplay, SLOT( setSegmentDown( int ) ) );
			connect( voxeldisplay, SIGNAL( segmentDownChanged( int ) ), downSlider, SLOT( setValue( int ) ) );

			connect( upSlider, SIGNAL( valueChanged( int ) ), voxeldisplay, SLOT( setSegmentUp( int ) ) );
			connect( voxeldisplay, SIGNAL( segmentUpChanged( int ) ), upSlider, SLOT( setValue( int ) ) );

		}
		else if( multi( fileName ) )
		{
			buffer = Buffer( fileName ) ;

			connect( layer, SIGNAL( valueChanged( int ) ), this, SLOT( getFile( int ) ) ) ;
			connect( this, SIGNAL( prepareFile( QString ) ), this, SLOT( open( QString ) ) ) ;

//		    open(files.at(0)) ;
		}
	}
}

void MainWindow::getFile( int i )
{
	int index = i ;

	bool updated = buffer.update() ;
	if( index >= buffer.size() || index < 0 )
	{
		index = 0 ;
	}
	if(updated)
		time->setRange( 0, buffer.size() - 1 );
	
	QString name = buffer.at( index ).file ;
	emit prepareFile( name );
}

void MainWindow::getLayer( int i )
{
	int f = buffer.fileInBuffer( triangledisplay->fileName ) ;

	if( f < 0 )
	{

	}
	else
	{
		int index = i ;

		if( index < 0 || i >= buffer.at( f ).size() )
		{
			index = 0 ;
			layer->setValue( index ) ;
			return ;
		}

		open( buffer.at( f ).file ) ;
	}

}

bool MainWindow::triangles( const QString s ) const
{

	QFile file( s ) ;

	if( !file.open( QIODevice::ReadOnly | QIODevice::Text ) )
	{
		return false;
	}

	QTextStream stream;
	stream.setDevice( &file );
	QString type ;
	stream >> type ;
	file.close();

	return ( type == "TRIANGLES" || type == "BIN_TRIANGLES" ) ;
}

bool MainWindow::voxels( const QString s ) const
{

	QFile file( s ) ;

	if( !file.open( QIODevice::ReadOnly | QIODevice::Text ) )
	{
		return false;
	}

	QTextStream stream;
	stream.setDevice( &file );
	QString type ;
	stream >> type ;
	file.close();

	return ( type == "VOXELS" ) ;
}

bool MainWindow::multi( const QString s ) const
{

	QFile file( s ) ;

	if( !file.open( QIODevice::ReadOnly | QIODevice::Text ) )
	{
		return false;
	}

	QTextStream stream;
	stream.setDevice( &file );
	QString type ;
	stream >> type ;
	file.close();

	return ( type == "MULTI" ) ;
}

void MainWindow::disconnect( TriangleGLDrawer *display )
{
	QObject::disconnect( zoom, SIGNAL( valueChanged( int ) ), display, SLOT( setZoom( int ) ) );
	QObject::disconnect( display, SIGNAL( zoomChanged( int ) ), zoom, SLOT( setValue( int ) ) );

	QObject::disconnect( alpha, SIGNAL( valueChanged( int ) ), display, SLOT( setSet( int ) ) );
	QObject::disconnect( display, SIGNAL( setChanged( int ) ), alpha, SLOT( setValue( int ) ) );

	QObject::disconnect( field, SIGNAL( valueChanged( int ) ), display, SLOT( setScale( int ) ) );
	QObject::disconnect( display, SIGNAL( scaleChanged( int ) ), field, SLOT( setValue( int ) ) );

	QObject::disconnect( layer, SIGNAL( valueChanged( int ) ), display, SLOT( setTimePlane( int ) ) );
	QObject::disconnect( display, SIGNAL( timePlaneChanged( int ) ), layer, SLOT( setValue( int ) ) );

	QObject::disconnect( downSlider, SIGNAL( valueChanged( int ) ), display, SLOT( setSegmentDown( int ) ) );
	QObject::disconnect( display, SIGNAL( segmentDownChanged( int ) ), downSlider, SLOT( setValue( int ) ) );

	QObject::disconnect( upSlider, SIGNAL( valueChanged( int ) ), display, SLOT( setSegmentUp( int ) ) );
	QObject::disconnect( display, SIGNAL( segmentUpChanged( int ) ), upSlider, SLOT( setValue( int ) ) );

	QObject::disconnect( printButton, SIGNAL( released() ), display, SLOT( grab() ) );

}

