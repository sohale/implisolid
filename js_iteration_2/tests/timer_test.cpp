#include <iostream>
using namespace std;

#include "../js_iteration_1/timer.hpp"

void try_things2()
{
    ;
}

void test_timer(){
    timer t;
    //t.timer_start();
    //timer_start();
    try_things2();
    //t.timer_start();

    auto r = t.stop();
    cout << r;
}

int main(){
    test_timer();
    return 0;
}


