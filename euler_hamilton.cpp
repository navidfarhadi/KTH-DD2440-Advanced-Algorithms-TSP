#ifndef EULER_HAMILTON_CPP
#define EULER_HAMILTON_CPP

#include "euler_hamilton.hpp"
#include "compute_distance.hpp"
#include <iostream>
#include <unordered_map>
#include <algorithm>

void changeToMap(std::vector<std::array<int,2>> &original_graph, std::unordered_map<int,std::vector<int>> &new_graph)
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

bool hasEulerianCircuit(std::unordered_map<int,std::vector<int>> &graph)
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

void findEulerianCircuit(std::unordered_map<int,std::vector<int>> &input_graph, std::vector<int> &eulerian_path, int node)
{
    std::unordered_map<int,std::vector<int>> input_graph_copy(input_graph);
    std::vector<int> nodeList;
    int currentNode = node;
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
            input_graph_copy[currentNode].pop_back();
            input_graph_copy[currentNodeNeighbor].erase(std::find(input_graph_copy[currentNodeNeighbor].begin(),input_graph_copy[currentNodeNeighbor].end(),currentNode));
            currentNode = currentNodeNeighbor;
        }
    }

    eulerian_path.push_back(currentNode);
    std::reverse(eulerian_path.begin(),eulerian_path.end());
}

void findHamiltonianCircuit(std::vector<int> &hamiltonian_circuit,int n)
{
    bool isNodeVisited[n];
    std::fill(isNodeVisited,isNodeVisited+sizeof(isNodeVisited),0);
    int cost = 0;

    for(int i = 0; i < hamiltonian_circuit.size(); i++)
    {
        if(isNodeVisited[hamiltonian_circuit[i]])
        {
            hamiltonian_circuit[i] = -1;
        }
        else
        {
            isNodeVisited[hamiltonian_circuit[i]] = true;
        }
    }
    hamiltonian_circuit.erase(std::remove(hamiltonian_circuit.begin(),hamiltonian_circuit.end(),-1),hamiltonian_circuit.end());
}

double findTotalCost(std::vector<int> &hamiltonian_circuit, std::vector<std::array<double,2>> &vertices)
{
    double cost = 0;
    int i;

    for(i = 0; i < hamiltonian_circuit.size() - 1; i++)
    {
        cost += compute_rounded_distance(hamiltonian_circuit[i],hamiltonian_circuit[i+1]);
    }

    cost += compute_rounded_distance(hamiltonian_circuit[i],hamiltonian_circuit[0]);

    return cost;
}

#endif
