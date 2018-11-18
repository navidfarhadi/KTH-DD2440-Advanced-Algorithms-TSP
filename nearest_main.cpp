#include<iostream>
#include <vector>
#include "read.hpp"
#include "2opt.hpp"
#include "nearest_neighbor.hpp"

int main_a()
{   
    clock_t init_clock = clock();

    int numVertices;
    std::cin >> numVertices;
    std::vector<std::array<double,2>> vertices(numVertices);
    read_in_stdin(vertices);

    std::vector<int> path;
    compute_nearest_neighbor(path,vertices);

    while(double(clock() - init_clock) / CLOCKS_PER_SEC < 1.99)
    {
            twoOpt(path, vertices, init_clock);
    }

    for(int i = 0; i < path.size(); i++)
    {
        std::cout << path[i] << std::endl;
    }

    return 0;
}