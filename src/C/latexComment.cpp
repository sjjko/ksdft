#include "latexComment.h"

latexComment::latexComment()
{
this->_latex ="";
this->_epsTemplateString="%!PS-Adobe-2.0 \n"
"%!PS-Adobe-2.0 \n"
"%%Title: anyTitle.eps \n"
"%%Creator: gnuplot 4.6 patchlevel 4 \n"
"%%CreationDate: Wed Feb 24 10:19:16 2016 \n"
"%%DocumentFonts: (atend) \n"
"%%BoundingBox: 50 50 402 1058 \n"
"%%Orientation: Portrait \n"
"%%Pages: (atend) \n"
"%%EndComments \n"
"%%BeginProlog \n"
"%%Trailer \n"
"%%DocumentFonts: Times-Roman \n"
"%%Pages: 1 \n";
}

latexComment::~latexComment()
{
    //dtor
}

void latexComment::startItemize()
{
 this->_latex += " \\begin{enumerate} ";
}

void latexComment::insertImage(const std::string imageName,const std::string captionString)
{
    #ifdef WITH_TEX
        this->_latex += " \\begin{figure}[ht] ";
        this->_latex += " \\begin{center} ";
        this->_latex += " \\includegraphics[scale=1.0]{./POSTPROCESSING/EPS/"+imageName+"} ";
        this->_latex += " \\caption{"+captionString+"} ";
        this->_latex += " \\end{center} ";
        this->_latex += " \\end{figure} ";
   #endif
}



void latexComment::endItemize()
{
 this->_latex += " \\end{enumerate} ";
}

void latexComment::addItem(std::string stringToAdd)
{
 this->_latex += " \\item  ";
 this->_latex += stringToAdd;
 this->_latex += " ";
}

void latexComment::writeHeader()
{

#ifdef WITH_TEX

 this->_latex =
    "\\documentclass{article} "
    "\\usepackage{graphicx} "
    "\\begin{document} ";
#endif
}

string latexComment::getString()
{
return this->_latex;
}

int latexComment::emptyLine()
{
    //is a matrix of size: N x Nions G(N,3) * X(3,Na) ---> (N,Na)
    #ifdef WITH_TEX
        this->_latex+=" \\newline ";
    #endif
    return 0;
}

int latexComment::section(string title)
{
    //is a matrix of size: N x Nions G(N,3) * X(3,Na) ---> (N,Na)
    #ifdef WITH_TEX
        this->_latex+=" \\noindent\\makebox[\\linewidth]{\\rule{\\paperwidth}{0.4pt}} ";
        this->_latex+=" \\section{ ";
        this->_latex+=title;
        this->_latex+="} ";
    #endif
    return 0;
}

int latexComment::subSection(string title)
{
    //is a matrix of size: N x Nions G(N,3) * X(3,Na) ---> (N,Na)
    #ifdef WITH_TEX
        this->_latex+=" \\begin{center} \\line(1,0){450} \\end{center} ";
        this->_latex+=" \\subsection{ ";
        this->_latex+=title;
        this->_latex+="} ";
    #endif
    return 0;
}

int latexComment::subsubSection(string title)
{
    //is a matrix of size: N x Nions G(N,3) * X(3,Na) ---> (N,Na)
    #ifdef WITH_TEX
        this->_latex+=" \\begin{center} \\line(1,0){250} \\end{center} ";
        this->_latex+=" \\subsubsection{ ";
        this->_latex+=title;
        this->_latex+="} ";
    #endif
    return 0;
}

int latexComment::newLine(string lineString)
{
    //is a matrix of size: N x Nions G(N,3) * X(3,Na) ---> (N,Na)
    #ifdef WITH_TEX
        this->_latex+=lineString;
        this->_latex+=" \\newline ";
    #endif
    return 0;
}

int latexComment::closeString()
{
#ifdef WITH_TEX
    //this->_latex+=" \\newline ";
    this->_latex+=" \\end{document}";
#endif
return 0;
}

 #ifdef USE_EXTERNAL_LIB_TEXCALLER
int latexComment::writeTheLatexDocument(string DocumentName)
{
        #ifdef WITH_TEX

try {
    std::string pdf;
    std::string info;

    texcaller::convert(pdf, info, this->_latex, "LaTeX", "PDF", 5);
    ofstream myfile(DocumentName);
    if (myfile.is_open())
    {
    myfile << pdf;
    myfile.close();
    }

    std::cout << "Generated PDF of " << pdf.size() << " bytes.";
    std::cout << " Details:" << std::endl << std::endl << info;
} catch (std::domain_error &e) {
    std::cout << "Error: " << e.what() << std::endl;
    return -1;
}

#endif
return 0;
}
#endif // USE_EXTERNAL_LIB_TEXCALLER

#ifdef INCLUDE_PICTURE_IN_LATEX

int latexComment::convertPPMToPS(string fileName)
{
//! \brief convert ppm to post script to include in latex document
//!
//! \param filePath: path to ppm file
//! \param fileName: name of ppm file to convert
//! \return 0 on success
//! fork into two processes, one of which converts the image

    pid_t child_pid;
    child_pid = fork();
    int status;
    int local=0;
   if (child_pid >= 0) /* fork succeeded */
    {
        if (child_pid == 0) /* fork() returns 0 for the child process */
        {
            printf("child process!\n");

            // Increment the local and global variables
            local++;
            //global++;

            printf("child PID =  %d, parent pid = %d\n", getpid(), getppid());
            //printf("\n child's local = %d, child's global = %d\n",local,global);

            const char *fn=fileName.c_str();
            char *current_path;
            getcwd(current_path, 255);
            char *inputFile[]={current_path,(char*)fn,".ppm",(char*)0};
            char *outputFile[]={current_path,(char*)fn,".ps",(char*)0};
            char *argv[] = {"/usr/bin/convert", *inputFile, " -bordercolor White ", *outputFile,0};
            return execv("/usr/bin/convert",argv); // call whoami command

         }
         else /* parent process */
         {
             printf("parent process!\n");
             printf("parent PID =  %d, child pid = %d\n", getpid(), child_pid);
             wait(&status); /* wait for child to exit, and store child's exit status */
             printf("Child exit code: %d\n", WEXITSTATUS(status));

             //The change in local and global variable in child process should not reflect here in parent process.
             //printf("\n Parent'z local = %d, parent's  global = %d\n",local,global);

             printf("Parent says bye!\n");
             //exit(0);  /* parent exits */
         }
    }
    else /* failure */
    {
        perror("fork");
        exit(0);
    }

    return 0;
}

#endif



