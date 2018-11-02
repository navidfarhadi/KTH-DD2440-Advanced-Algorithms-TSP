#include <vector>
#include <array>
#include <iostream>
#include "read.hpp"

/*
 * cpp file for the implementation
 * to read in the vertices
 */

void read_in_stdin(std::vector<std::array<double,2>> &v)
{
    std::cout << "Read data in ...\n";
    for(auto &line : v)
    {
        for(auto &d : line)
        {
            std::cin >> d;
            std::cout << d << " ";
        }
        std::cout << "\n";
    }
    //std::cout << (double)(v[1][1]) << "\n";
}