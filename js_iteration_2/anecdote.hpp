#pragma once
#include <iostream>
#include "boost/multi_array.hpp"

#define USE_ANECDOTE false


// This is not logging for run-time. Just for debugging, and for documenting the code.
#if true  // (ASSERT_USED && USE_ANECDOTE)
    template<typename T>
    inline void anecdote(const T& x, int level=0) {
        // reports, narrates, anecdotes, commentates the results.
        std::cout << x;
    }
    /*
    inline void anecdote_s(int level=0) {
        std::cout << std::endl;
    }
    */

    template<typename T>
    inline void anecdote_c(const boost::multi_array<T,2>& c, int level=0) {
        int size = c.size();
        int ctr = 0;
        for( auto i = c.begin(); i != c.end(); ++i ) {
            for( auto j = i->begin(); j != i->end(); ++j ) {
                std::cout << *j << ", ";
            }
            std::cout << ".     ";
            ctr++;
            if (ctr > 4) {
                std::cout << "... [and " << size-ctr <<" more].."; //<< std::endl;
            }
        }
    }
    template<typename IT>
    inline void anecdote_c(IT begin, IT end, int level=0) {
        std::cout << "[";
        int size = std::distance(begin, end);
        int ctr = 0;
        for( auto i = begin; i != end; ++i ) {
            std::cout << *i;
            auto q = i;
            ++q;
            if (q != end) {
                std::cout << ", ";
            }
            ctr++;
            if (ctr > 4) {
                std::cout << "... [and " << size-ctr <<" more].."; //<< std::endl;
            }
        }
        std::cout << "]";
    }

    /*
    template<typename T>
    inline void anecdote(...) {
        va_list args;
        va_start(args, fmt);

        T v = va_arg(args, T);
    }
    */

#else
    #pragma message("Anecdotes ignored")
    template<typename T>
    inline void anecdote(const T& x, int level=0) {
        // no absolutely nothing.
        T y;
    }
#endif

/*
#include <type_traits>

// Helper to determine whether T is a container.
template<typename T>
struct is_container
{
private:
    // whether it has an iterator
    template<typename C> static char test(typename C::iterator*);  // If C has a C::iterator
    template<typename C> static long test(...);   // otherwise
public:
    enum { value = sizeof(test<T>(0)) == sizeof(char) };  // compiletime?
};
*/

/*
template <typename Container>
typename std::enable_if<is_container<Container>::value,  // the if's condition
                        void>::type
bar(const Container &c, typename Container::value_type const & t)
{
  // Note: no extra check needed for value_type, the check comes for
  //       free in the function signature already.
}
*/
/*
#define ONLY_IF_CONTAINER typename std::enable_if<is_container<Container>::value,  // the if's condition \
                        void>::type
*/