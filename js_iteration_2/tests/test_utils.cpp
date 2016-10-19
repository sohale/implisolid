#include "gtest/gtest.h"

#include <vector>
#include <iostream>

#include "../basic_data_structures.hpp"
#include "../basic_functions.hpp"
using std::vector;
// todo: remove "using" of ublas ('s vector)

template<typename T>
void print(T x)
{
    for (typename T::iterator i = x.begin(); i != x.end(); ++i) {
        std::cout << *i << " ";
    }
    std::cout << std::endl;
}

TEST(subsetrandom, random_subset) {

    std::vector<int> x {1,2,3,4};
    int m = 2;
    chisle_random_subset(x, m);
    print(x);
    EXPECT_EQ(x.size(), m)  << "selecting " << m << " from size " << 4;

    chisle_random_subset(x, m);
    EXPECT_EQ(x.size(), m)  << "selecting " << m << " from size " << m;


    chisle_random_subset(x, 0);
    EXPECT_EQ(x.size(), 0) << "selecting 0 from size 0";

    chisle_random_subset(x, 4);
    EXPECT_EQ(x.size(), 0) << "selecting 4 from size 0";

    int n=5;
    std::vector<int> x100(n);
    for (int i=0;i<n;++i) x100[i]=i;
    for (int m2 = n; m2 >= 0; --m2) {
        chisle_random_subset(x100, m2);
        // print(x100);
        EXPECT_EQ(x100.size(), m2);
    }

}
