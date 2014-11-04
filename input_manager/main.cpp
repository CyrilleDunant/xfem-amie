
#include <QApplication>
#include <QFileDialog>
//#include <QtWidgets/QApplication>
//#include <QtWidgets/QFileDialog>
#include <unistd.h>
#include <fstream>
#include <iostream>
#include <string>

int main(int argc, char *argv[]) 
{
	QApplication app(argc, argv);

	QString input = QFileDialog::getOpenFileName(NULL, "Open AMIE input file", "../examples/data/composite","*.ini") ;

	std::ifstream pathFile ;
	pathFile.open("path_to_executable.ini", std::ios::in) ;
	std::string path ;
	getline( pathFile, path ) ;
	std::string command = path+" "+input.toStdString() ;
	#ifdef _WIN32
	std::string winCommand = command ;
	size_t backslash = winCommand.find("/") ;
	while(backslash < std::string::npos)
	{
		std::string left = winCommand.substr(0, backslash) ;
		std::string right = winCommand.substr(backslash+1) ;
		winCommand = left+"\\"+right ;
		backslash = winCommand.find("/") ;
	}
	command = winCommand ;
	#endif
    std::cout << command << std::endl ;
	system( command.c_str() ) ;
	

	return 0 ;
}
