#include <iostream>
using namespace std;

#include <chrono>

//now sure about the difference between the two classes: steady_clock and system_clock

//auto timer_start_time = chrono::steady_clock::now();
//std::chrono::time_point<qpc_clock, duration> timer_start_time = chrono::steady_clock::now();

//system_clock::time_point
//failed: system_clock::time_point timer_start_time = chrono::steady_clock::now();


//std::chrono::system_clock clock1;
//std::chrono::steady_clock clock2;
//std::chrono::system_clock::time_point timer_start_time = std::chrono::steady_clock::now();

//std::chrono::steady_clock
//: error: no viable conversion from
//'time_point<std::__1::chrono::steady_clock, duration<[...], ratio<[...], 1000000000>>>'
//to
//'time_point<std::__1::chrono::system_clock, duration<[...], ratio<[...], 1000000>>>'

typedef std::chrono::steady_clock::time_point  my_time_point_t;
typedef std::chrono::steady_clock::duration  my_duration_t;
//common_type = ?
// has_facet =?

my_time_point_t timer_start_time = std::chrono::steady_clock::now();


class timer{
    my_time_point_t timer_start_time; // = std::chrono::steady_clock::now();

    void timer_start(){
        timer_start_time = std::chrono::steady_clock::now();

        //time_point

    }

public:
    timer(){
        this->timer_start();
    }

    auto stop(){
        my_time_point_t timer_end_time = chrono::steady_clock::now();
        my_duration_t diff = timer_end_time - timer_start_time;
        cout << endl << "Duration: ";
        double duration = chrono::duration <double, milli> (diff).count();
        cout << duration << " msec" ;
        //cout << " (" << chrono::duration <double, nano> (diff).count() << " ns" << ")";
        cout << endl;
        return duration;
    }
};





/*
class timer
{
    //types:
    //std::chrono::time_point<int>   <--- <int> is wrong
    //std::chrono::system_clock::time_point  <--- did not work
    auto timer_start_time = chrono::steady_clock::now();

    timer(){
        this->timer_start_time = chrono::steady_clock::now();
    }
    void stop()
    {
        auto timer_end_time = chrono::steady_clock::now();
        auto diff = timer_end_time - this->timer_start_time;
        cout << endl << "Duration: ";
        auto duration = chrono::duration <double, milli> (diff).count();
        cout << duration << " msec" ;
        //cout << " (" << chrono::duration <double, nano> (diff).count() << " ns" << ")";
        cout << endl;
        return duration;
    }
};
*/