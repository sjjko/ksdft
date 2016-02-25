#include "inputParser.h"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

/** \brief Class for input data processing
 * uses BOOST property_tree and reads input data ascii file
 */

using namespace myFunctions;

inputParser::inputParser(const string inputFilename,struct paramStruct *parameterStruct)
{
//    /** \brief set up input parser class
//     *
//     * \param inputFilename ... filename of input file
//     * \param parameter struct containing containers for input data
//     *
//     */
    this->_inputFilename=inputFilename;
    this->_inputStruct=parameterStruct;
    if (!std::ifstream(_inputFilename))
    {
        std::cout << "input file does not exist! - write a template one!" << std::endl;
        if(!this->writeTemplateInputFile())
        {
        verbosity("could not write template file - EXIT!",0,__FILE__,__LINE__);
        }
    }
}

inputParser::~inputParser()
{
    //dtor
}

int inputParser::writeTemplateInputFile()
{
//! \brief write a template input file
ofstream myfile("inputFile.param");
  if (myfile.is_open())
  {
    myfile << "[physical]\n";
    myfile << "Z=1\n";
    myfile << "f=2\n";
    myfile << "atomic coordinates file=\"./atoms.param\"";
    myfile << "[geometry]\n";
    myfile << "Nx=20\n";
    myfile << "Ny=25\n";
    myfile << "Nz=30\n";
    myfile << "Rx=6\n";
    myfile << "Ry=6\n";
    myfile << "Rz=6\n";
    myfile << "[numerical]\n";
    myfile << "alpha=0.5e-3\n";
    myfile << "no. iterations steepest descent=250\n";
    myfile << "no. iterations conjugate gradient=350\n";
    myfile << "number of electrons (wavefunctions)=2\n";
    myfile << "[general]\n";
    myfile << "VerbosityLevel=1\n";
    myfile.close();
  }
  else
  {
      cout << "Unable to open file";
      return -1;
  }
  return 0;
}

int inputParser::parseInput()
{
    return this->readInput();
}

arma::mat inputParser::readAtomicCoordinates()
{
  //! \brief read atomic coordinates from input file - name given in input .param file
    //!
    //! \param all internal
    //! \return matrix with coordinates of atoms
    arma::mat returnMatrix;

    cout << "inputParser: read atomic coordinate file: " << this->_coordFilename << endl;
    returnMatrix.load(this->_coordFilename,arma::csv_ascii);
    return returnMatrix;

}


int inputParser::writeSampleAtomicCoordinateFile()
{
  //! \brief if no atomic coordinate file found write a sample one
    //!
    //! \param all internal
    //! \return 1 at success


    _coordFilename="atoms.param";
    ofstream outfile(_coordFilename);

    verbosity("we write an atom file for two single atoms system",1,__FILE__,__LINE__);
    outfile << "0,0,0" << endl;
    outfile << "1,0,0" << endl;
    outfile.close();

    verbosity("finished writing - exit",1,__FILE__,__LINE__);
    return 1;

}

int inputParser::readInput()
{
    using boost::property_tree::ptree;
    ptree pt;

    //namespace pt = boost::property_tree;

    read_ini(this->_inputFilename, pt);

    for (auto& section : pt)
    {
        std::cout << '[' << section.first << "]\n";
        //for (auto& key : section.second)
        //    std::cout << key.first << "=" << key.second.get_value<std::string>() << "\n";
        if(section.first=="physical")
        {
        for (auto& key : section.second)
            {
            if(key.first=="Z") this->_inputStruct->Z=key.second.get_value<double>();
            else if(key.first=="f") this->_inputStruct->f=key.second.get_value<double>();
            else if(key.first=="atomic coordinates filename") this->_coordFilename=key.second.get_value<std::string>();
            }
        }
        else if(section.first=="numerical")
        {
        for (auto& key : section.second)
            {
            if(key.first=="alpha") this->_inputStruct->alpha=key.second.get_value<double>();
            else if(key.first=="no. iterations steepest descent") this->_inputStruct->sdNit=key.second.get_value<int>();
            else if(key.first=="no. iterations conjugate gradient") this->_inputStruct->pccgNit=key.second.get_value<int>();
            else if(key.first=="number of electrons (wavefunctions)") this->_inputStruct->number_of_wavefunctions=key.second.get_value<int>();
            }
        }
        else if(section.first=="geometry")
        {
        this->_inputStruct->R=arma::eye(3,3);
        this->_inputStruct->S=arma::vec(3,fill::zeros);
        for (auto& key : section.second)
            {
            if(key.first=="Nx") this->_inputStruct->S(0)=key.second.get_value<double>();
            if(key.first=="Ny") this->_inputStruct->S(1)=key.second.get_value<int>();
            if(key.first=="Nz") this->_inputStruct->S(2)=key.second.get_value<int>();
            if(key.first=="Rx") this->_inputStruct->R(0,0)=key.second.get_value<double>();
            if(key.first=="Ry") this->_inputStruct->R(1,1)=key.second.get_value<int>();
            if(key.first=="Rz") this->_inputStruct->R(2,2)=key.second.get_value<int>();
            }
        }
        else if(section.first=="general")
        {
        for (auto& key : section.second)
            {
            if(key.first=="VerbosityLevel") this->_inputStruct->globalVL=key.second.get_value<int>();
            }
        }
    this->_inputStruct->prodS=arma::prod(this->_inputStruct->S);

    }

    myFunctions::cassert(((this->_inputStruct->Z>0) && (this->_inputStruct->Z<1000)),ISCRITICAL,"inputParser: Z out of bounds!",__FILE__,__LINE__);
    myFunctions::cassert(((this->_inputStruct->number_of_wavefunctions>0) && (this->_inputStruct->number_of_wavefunctions<500)),
    ISCRITICAL,"inputParser: number of electrons out of bounds!",__FILE__,__LINE__);

    cout << "Z is set to " << this->_inputStruct->Z << endl;
    cout << "The system has " << this->_inputStruct->number_of_wavefunctions << " electrons!" << endl;
    cout << "We compute " << this->_inputStruct->sdNit << " iterations with steepest descent optimizer!" << endl;
    cout << "We compute " << this->_inputStruct->pccgNit << " iterations with conjugate gradient optimizer!" << endl;

    return 0;
}




