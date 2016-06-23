
/** **********************************************
 *     my_assert()                               *
 *********************************************** */
/*
#include "boost/assert.hpp"
//#define my_assert(cond, msg)   BOOST_ASSERT_MSG((cond), (cond))

//Usage:  my_assert(execSize == (_oldSize - remaining), "execSize : " << execSize << ", _oldSize : " << _oldSize << ", remaining : " << remaining);

//#if defined ASSERT_ENABLED
    #define my_assert(cond, msg) {\
        if(!(cond))\
        {\
            std::stringstream str;\
            str << msg;\
            BOOST_ASSERT_MSG(cond, str.str().c_str());\
        }\
    }
//#else
//    #define my_assert(...)
//#endif
//#define my_assert(cond)   my_assert(cond, "")

////////////////////////////////////
*/

#include <cassert>

#define my_assert(cond, msg)  assert(cond)
