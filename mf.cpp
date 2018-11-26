/*
 * This is the cpp file for the implementation
 * of the Multiple Fragment Heuristic approach.
 * The goal is to produce an initial tour in
 * a better time than Christofides, but still
 * have low total cost.
 */

#include "mf.hpp"
#include <queue>
#include <array>
#include <algorithm>
#include <iostream>

struct compa
{
    bool operator () (const std::array<int,2> &lhs, const std::array<int,2> &rhs)
    {
        return (compute_rounded_distance(lhs[0],lhs[1]) > compute_rounded_distance(rhs[0],rhs[1]));
    }
}comp;

bool find_tour(std::vector<std::vector<int>> &arr, std::array<int,2> &edge)
{
    int prev = edge[1];
    if(arr[prev].size() == 0)
    {
        return false;
    }
    int curr = arr[prev][0];
    while(true)
    {
        if(arr[curr].size() == 1)
        {
            // only the back reference
            // so we are good
            return false;
        }
        else
        {
            // there are two elements in
            // find back reference and forward reference
            int forward = (arr[curr][0] == prev) ? arr[curr][1] : arr[curr][0];
            if(forward == edge[0])
            {
                //std::cout << "Indirect loop\n";
                return true;
            }
            else
            {
                prev = curr;
                curr = forward;
            }
        }
    }
    return false;
}

void create_inital_tour(int numVertices, std::vector<int> &init_tour)
{
    std::vector<std::vector<int>> arr(numVertices);
    //std::cout << "1\n";
    std::vector<std::array<int,2>> list;
    for(int i = 0; i < numVertices; i++)
    {
        for(int j = 0; j <= i; j++)
        {
            if(i != j)
            {
                list.push_back({i,j});
            }
        }
    }
    //std::cout << "2\n";
    std::sort (list.begin(), list.end(),comp);

    //std::cout << "3\n";


    int counter = 1;
    int starting_node = list.back()[0];
    arr[list.back()[0]].push_back(list.back()[1]);
    arr[list.back()[1]].push_back(list.back()[0]);
    list.pop_back();
    while(counter < numVertices - 1)
    {
        //std::cout << "counter = " << counter << "\n";
        std::array<int,2> edge = list.back();
        //std::cout << "edge " << edge[0] << " - " << edge[1] << "\n";
        //std::cout << "step 4\n";
        if(arr[edge[0]].size() == 1)
        {
            if(arr[edge[0]][0] == edge[1])
            {
                list.pop_back();
                //std::cout << "Direct loop\n";
                continue;
            }
        }
        if(arr[edge[0]].size() < 2 && arr[edge[1]].size() < 2 && !find_tour(arr, edge))
        {
            arr[edge[0]].push_back(edge[1]);
            arr[edge[1]].push_back(edge[0]);
            counter++;
            //std::cout << "IN\n";
        }
        list.pop_back();
    }

    //std::cout << "Start traversing\n";
    for(int i = 0; i < numVertices; i++)
    {
        if(arr[i].size() == 1)
        {
            starting_node = i;
            //std::cout << "starting node " << starting_node << "\n";
            break;
        }
    }

    // traverse into init_tour
    counter = 1;
    int prev = starting_node;
    // has to exist
    int curr = arr[prev][0];
    init_tour.push_back(prev);
    init_tour.push_back(curr);
    while(counter < numVertices - 1)
    {
        int forward = (arr[curr][0] == prev) ? arr[curr][1] : arr[curr][0];
        init_tour.push_back(forward);
        prev = curr;
        curr = forward;
        counter++;
    }


    //std::cout << "END\n";
}