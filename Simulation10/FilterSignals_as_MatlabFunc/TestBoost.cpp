#include <boost/version.hpp>
#include "boost/multi_array.hpp"
#include <iostream>
#include <iomanip>

using boost::multi_array;

typedef boost::multi_array_types::extent_range range;
typedef boost::multi_array<double, 3> double3d; // for signaling marginals

int main()
{
    std::cout << "Boost version: "
    << BOOST_VERSION / 100000
    << "."
    << BOOST_VERSION / 100 % 1000
    << "."
    << BOOST_VERSION % 100
    << std::endl;
    
    int W = 2; int M = 2;
    
    double3d eta1(boost::extents[W][M][M]);
    
    for (int m2=0;m2<M;m2++)
        for (int m1=0;m1<M;m1++)
            for (int w1=0;w1<W;w1++){
                int k1=(w1*M+m1)*M+m2;
                eta1[w1][m1][m2]=k1+1;
            }
    
    std::cout << "eta1 = "
    << eta1[1][1][1]
    << std::endl;
    
    return 0;
}
