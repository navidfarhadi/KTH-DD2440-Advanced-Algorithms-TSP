#include "odd_degree_subgraph.hpp"
#include <iostream>

/*
 * hpp file for the implementation
 * of finding the vertices with an
 * odd degree in a MST and creating
 * the subgraph out of it.
 */

void find_odd_v(std::vector<std::array<int,2>> &e, std::vector<int> &v, int size)
{
    // std::cerr << "find_odd_v\n";
    int deg[size] = {0};
    // increase for every ingoing
    // or outgoing edge
    for(auto &line : e)
    {
        deg[line[0]]++;
        deg[line[1]]++;
    }
    for(int i = 0; i < size; i++)
    {
        // std::cerr << "deg[" << i << "] = " << deg[i] << "\n";
        if(deg[i] %2 == 1)
        {
            // it is an vertex with an odd degree
            // then we include it
            // std::cerr << "included vertex " << i << "\n";
            v.push_back(i);
        }
    }
}

void create_subgraph(std::vector<int> &v, std::vector<std::array<int,2>> &e)
{
    // implementation does not make to much sense
    // just wanted fo fullfil what I had to do
    // std::cerr << "create_subgraph\n";
    for(int i = 0; i < (int)(v.size()); i++)
    {
        for(int j = i + 1; j < (int)(v.size()); j++)
        {
            e.push_back({v[i],v[j]});
            // std::cerr << "edge: " << v[i] << " " << v[j] << "\n";
        }
    }
}