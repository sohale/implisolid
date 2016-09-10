#include "gtest/gtest.h"

#include <iostream>
using namespace std;

#include "../../js_iteration_1/timer.hpp"

void try_things2()
{
    ;
}

TEST(AAA, timer) {
//void test_timer(){
    timer t;
    //t.timer_start();
    //timer_start();
    try_things2();
    //t.timer_start();

    auto r = t.stop();
    cout << r;
}

/*
int main(){
    test_timer();
    cout << std::endl;
    cout << "done" << std::endl;
    return 0;
}
*/

