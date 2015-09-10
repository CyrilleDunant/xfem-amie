
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
        if(path.size() == 0) { path = "./2d_composite" ; }
        std::string stdinput = input.toStdString() ;
	std::string command = path+" "+stdinput ;
        std::string dir = std::string() ;
        if(stdinput.find("/") != std::string::npos)
            dir = stdinput.substr( 0, stdinput.find_last_of("/")) ;
	if(input.size() == 0)
		return 0 ;
        bool hasDir = false ;
        for(int i = 1 ; i < argc ; i++)
        {
            if(hasDir)
                command += " " + dir+"/"+std::string(argv[i]) ;
            else
                command += " " + std::string(argv[i]) ;

            hasDir = (std::string(argv[i]) == "--directory") ;   
        }
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
	system( command.c_str() ) ;
	

	return 0 ;
}
