#ifndef CASSERT_H
#define CASSERT_H

#include "main.h"
#include "structs.h"

//void cassert(string s,string file,float linie);

//int GLOBAL_verbosityLevel=0;
namespace myFunctions
{
inline void cassert(bool test,bool isCritical,string s,string file,float linie)
{
    /** \brief
     *
     * \param test: condition to test
     * \param s: message if test failed
     * \param file: current source file
     * \param float: current source line
     *
     */

    if(!test)
    {
        cout << s << " \n TEST FAILED " << endl;
        cout << "We are in file " << file << " on line " << linie  << endl;
        if(isCritical)
        {
        cout << "error is critical - exit" << endl;
        exit(EXIT_FAILURE);
        }
        else
        {
        cout << "error is not critical - go on" << endl;
        }
    }
}

inline void verbosity(paramStruct Pa, string s,int verbosityL,string file,float linie)
{

    /** \brief prints message according to global level of verbosity flag set in input file
     *
     * \param s: message to display
     * \param verbosityL: level of verbosity - compared to global one
     * \param file: current source file
     * \param linie: current source line
     */

    //!< simple print function - prints if global verbosity level is greater than given
    if(verbosityL < Pa.globalVL)
    {
        cout << file << " " << linie << " " <<  s  << endl;
    }

}
}



#endif // ASSERT_H
