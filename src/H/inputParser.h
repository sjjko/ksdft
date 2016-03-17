#ifndef INPUTPARSER_H
#define INPUTPARSER_H

#include "main.h"
#include "structs.h"
#include "customAssert.h"
#include "armadillo"

/** \brief class for input data processing
 * uses BOOST property_tree and reads input data ascii file
 */

class inputParser
{

    public:
        inputParser(const string pathToCase, const string inputFilename,struct paramStruct *parameterStruct);
        virtual ~inputParser();
        int parseInput();
        arma::mat readAtomicCoordinates();
        int writeSampleAtomicCoordinateFile();
        int readInput();
        int writeTemplateInputFile();

    private:
        string _inputFilename;
        struct paramStruct* _inputStruct;
        string _coordFilename;
        string _pathToCase;

};

#endif // INPUTPARSER_H
