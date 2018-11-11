#ifndef EULER_HAMILTON_CPP
#define EULER_HAMILTON_CPP

#include "euler_hamilton.hpp"
#include <iostream>

void changeToMap(std::vector<std::array<int,2>> &original_graph, std::map<int,std::vector<int>> &new_graph)
{
    for(std::array<int,2> &edge : original_graph)
    {
        int edge_0 = edge[0];
        int edge_1 = edge[1];
        if(new_graph.count(edge_0))
        {
            new_graph[edge_0].push_back(edge_1);
        }
        else
        {
            std::vector<int> temp;
            temp.push_back(edge_1);
            new_graph[edge_0] = temp;
        }

        if(new_graph.count(edge_1))
        {
            new_graph[edge_1].push_back(edge_0);
        }
        else
        {
            std::vector<int> temp;
            temp.push_back(edge_0);
            new_graph[edge_1] = temp;
        }
    }
}

bool hasEulerianCircuit(std::map<int,std::vector<int>> &graph)
{
    for (const auto &p : graph)
    {
        if(p.second.size() % 2 != 0)
        {
            return 0;
        }
    }

    return 1;
}

#endif