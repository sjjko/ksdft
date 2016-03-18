#include "timerClass.h"

using boost::timer::cpu_timer;
using boost::timer::cpu_times;
using boost::timer::nanosecond_type;
using namespace myFunctions;

timerClass::timerClass(paramStruct Pa)
{

verbosity(Pa,"Initialize the timer class ",0,__FILE__,__LINE__);

//cpu_times const elapsed_times(timer.elapsed());
//cout << "REPORT TIME: " << elapsed_times.user << " "  << elapsed_times.system << endl;
//exit(-1);



}

void timerClass::getActualTime()
{
//std::cout << format(_timerInstance, _placesForOutput) << std::endl;

}

timerClass::~timerClass()
{
    //dtor
}
