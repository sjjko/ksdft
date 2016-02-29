#ifndef INPUTARGUMENTS_H_INCLUDED
#define INPUTARGUMENTS_H_INCLUDED

verbosity(Pa,"start with "+std::to_string(argc)+"input arguments.",2,__FILE__,__LINE__);
myFunctions::cassert(argc==2,NOTCRITICAL,"give case name as first argument to ksdft++! (folder has to exist in ./case directory!)",__FILE__,__LINE__);
Pa.caseName=argv[1];
string pathToCase="./case/"+Pa.caseName;
verbosity(Pa,"check if we have a true case directory at "+pathToCase,2,__FILE__,__LINE__);
string paramFile=pathToCase+"/inputFile.param";
ifstream f(paramFile.c_str());
if(!f.good())
{
f.close();
verbosity(Pa,"maybe we are in bin/Release folder - check this "+pathToCase,2,__FILE__,__LINE__);
string pathToCase="../../case/"+Pa.caseName;
verbosity(Pa,"check if we have a true case directory at "+pathToCase,2,__FILE__,__LINE__);
string paramFile=pathToCase+"/inputFile.param";
ifstream f2(paramFile.c_str());
myFunctions::cassert(f2.good(),ISCRITICAL,"could not find a inputFile.param file in input directory and in directory "+ pathToCase +" - Exit!",__FILE__,__LINE__);
f2.close();
}
verbosity(Pa,"our input file for this run is "+paramFile,2,__FILE__,__LINE__);

#endif // INPUTARGUMENTS_H_INCLUDED
