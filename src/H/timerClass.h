#ifndef TIMERCLASS_H
#define TIMERCLASS_H

#include <stdio.h>
#include <boost/timer/timer.hpp>
#include <cmath>
#include "customAssert.h"
#include "structs.h"


class timerClass
{
    public:
        /** Default constructor */
        timerClass(paramStruct Pa);
        /** Default destructor */
        virtual ~timerClass();
        void getActualTime();

    protected:
    private:
    boost::timer::cpu_timer _timerInstance;
    const short _placesForOutput=3;
};

#endif // TIMERCLASS_H
