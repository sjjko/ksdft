#ifndef LATEXCOMMENT_H
#define LATEXCOMMENT_H

#include "main.h"
#include <unistd.h>
#include <sys/types.h>
#include <errno.h>
#include <stdio.h>
#include <sys/wait.h>
#include <stdlib.h>
#include "customAssert.h"
#include <texcaller.h>
//include to generate the image folder structure
#include <sys/types.h>
#include <sys/stat.h>

/** \brief insert a latex string comment and write latex document of computations in the code
 */

class latexComment
{
    private:
        std::string _latex;
        std::string _epsTemplateString;

    protected:
        std::string _myFunctionString;  //!< latex string describing the current class

    public:
        latexComment();
        virtual ~latexComment();
        void writeHeader(); //!< write latex document header
        int newLine(string lineString); //!< enter new line into latex document - append newline
        int closeString(); //!< close the latex string by finishing the document
        int emptyLine(); //!< insert a newline
        int subSection(string title); //!< enter a subsection into the latex document
        int subsubSection(string title); //!< enter a subsubsection into the latex document
        int section(string title); //!< enter a section into the latex document
        int writeTheLatexDocument(string DocumentName); //!< write the string into a pdf document
        void startItemize(); //!< start a list to add items successively
        void endItemize(); //!< end the list
        void insertImage(const std::string imageName,const std::string captionString);//!< add image to latex document
        void addItem(std::string stringToAdd);  //!< add an item to the list
        inline void commentMyFunctionAsItem() {addItem(_myFunctionString);}; //!< used by classes to return their latex description
        inline void commentMyFunction() {this->newLine(_myFunctionString);}; //!< used by classes to return their latex description
        inline void setMyFunctionString(std::string inputString) {_myFunctionString=inputString;} //!< used by classes to set the temporary latex string which is concatenated to the global latex string by member functions of latexComment class
        #ifdef INCLUDE_PICTURE_IN_LATEX
        int convertPPMToPS(string fileName);
        #endif // LATEXCOMMENT_H

        string getString();

};

#endif // LATEXCOMMENT_H
