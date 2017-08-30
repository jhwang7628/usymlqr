#ifndef SIMPLE_TIMER_H
#define SIMPLE_TIMER_H
#include <iostream>
#include <chrono>

//##############################################################################
// Class SimpleTimer
//##############################################################################
template<typename T>
class SimpleTimer
{
    private:
        std::chrono::high_resolution_clock::time_point  _start;
        std::chrono::high_resolution_clock::duration    _elapsed;
        int _count = 0; 

    public:
        SimpleTimer() : _elapsed(0) {}

        inline void Start() { _start = std::chrono::high_resolution_clock::now();}
        inline void Pause() { _elapsed += std::chrono::high_resolution_clock::now() - _start; _count++; }
        inline T Duration() { return (std::chrono::duration_cast<std::chrono::duration<T>>(_elapsed)).count(); }
        inline T DurationAverage() { return _count>0 ? Duration()/(T)_count : -1.0; }

        //// debug methods ////
        inline void DurationDebug()
        {
            std::cout << "elapsed(native) = " << _elapsed.count() << std::endl;
            std::cout << "elapsed(T)      = " << (std::chrono::duration_cast<std::chrono::duration<T>>(_elapsed)).count() << std::endl;
        }
};
//##############################################################################
#endif
