#pragma once

#include <typeinfo>

#include "implicit_function.hpp"

std::vector<mp5_implicit::implicit_function*> garbage_collector;


void register_new_object(mp5_implicit::implicit_function* p) {
    garbage_collector.push_back(p);
}

void gc_objects() {
   int counter = 0;
   while (!garbage_collector.empty()) {

        mp5_implicit::implicit_function* e = garbage_collector.back(); // garbage_collector.front()
       if (false) {
           loger << "delete(" << e << typeid(*e).name() << ")   "; // << std::endl << std::flush;
       }

       garbage_collector.pop_back();
       delete e;
       counter++;
   }

    if (false) {
        loger << std::endl << std::flush;
        loger << counter << " garbages collected. " << garbage_collector.size()<< std::endl << std::flush;
    }
}
