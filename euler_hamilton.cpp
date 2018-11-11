#ifndef EULER_HAMILTON_CPP
#define EULER_HAMILTON_CPP

#include "euler_hamilton.hpp"
#include <iostream>
#include <algorithm>

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

void findEulerianCircuit(std::map<int,std::vector<int>> &input_graph, std::vector<int> &eulerian_path)
{
    std::map<int,std::vector<int>> input_graph_copy(input_graph);
    std::vector<int> nodeList;
    int currentNode = 0;
    int currentNodeNeighbor;

    while(true)
    {
        if(input_graph_copy[currentNode].size() == 0 && nodeList.size() == 0)
        {
            break;
        }

        if(input_graph_copy[currentNode].size() == 0)
        {
            eulerian_path.push_back(currentNode);
            currentNode = nodeList.back();
            nodeList.pop_back();
        }
        else
        {
            nodeList.push_back(currentNode);
            currentNodeNeighbor = input_graph_copy[currentNode].back();
            // input_graph_copy[currentNode].pop_back();
            input_graph_copy[currentNode].erase(std::remove(input_graph_copy[currentNode].begin(),input_graph_copy[currentNode].end(),currentNodeNeighbor),input_graph_copy[currentNode].end());
            input_graph_copy[currentNodeNeighbor].erase(std::remove(input_graph_copy[currentNodeNeighbor].begin(),input_graph_copy[currentNodeNeighbor].end(),currentNode),input_graph_copy[currentNodeNeighbor].end());
            currentNode = currentNodeNeighbor;
        }
    }

    eulerian_path.push_back(currentNode);
}

#endif