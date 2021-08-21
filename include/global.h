#pragma once
#include <vector>
#include <cmath>
#include <string>
#include <limits>
#include <time.h>

/*
 * to-do list
 * gamma
 * anti-aliasing
 * tri class
 * memory management
 */

const double epsilon = std::numeric_limits<double>::epsilon();

class MyTimer
{
private:
    clock_t tick;

public:
    void start()
    {
        tick = clock();
    }
    void end()
    {
        printf("%fms\n\n", (double)(clock() - tick) / CLOCKS_PER_SEC);
    }
};
