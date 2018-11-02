#include <iostream>
#include <vector>
#include <array>
#include <algorithm>
#include <cstdlib>
#include <tuple>
#include <limits>
#include <cmath>
#include "read.hpp"

/*
 * This is the file that implements the minimum spanning tree
 * algorithm. We use the algorithm of Prim
 */

class MST
{
    public:
        MST();
        void getMST(std::vector<std::array<int,2>> &e);
        std::vector<int> parent;
        std::vector<std::array<double,2>> vertices;
        double computeEdgeDist(int ind1, int ind2);
        double computeDistance(std::array<double,2> v1, std::array<double,2> v2);
};

MST::MST()
{
    std::cout << "MST\n";
}

struct compare
{
    bool operator () (const std::pair<int,double> &lhs, const std::pair<int,double> &rhs)
    {
        return (lhs.second > rhs.second);
    }
};

void MST::getMST(std::vector<std::array<int,2>> &e)
{
    int size = vertices.size();
    double dist[size];
    int pred[size];
    bool visited[size] = {false};
    compare comp;
    // we always start with vertic index 0
    // because it does not make a difference
    
    // distance zero
    dist[0] = 0;
    // all other have max distance
    for(int i = 1; i < size; i++)
    {
        dist[i] = __DBL_MAX__;
    }
    // -2: root node
    // -1: no predecessor up to now
    pred[0] = -2;
    for(int i = 1; i < size; i++)
    {
        pred[i] = -1;
    }
    // boolean value stands for if it is invalid
    // true: invalid
    // false: valid
    std::vector<std::pair<int,double>> pq;

    pq.push_back({0,0});

    std::make_heap(pq.begin(),pq.end(),comp);

    while (!pq.empty())
    {
        std::pop_heap(pq.begin(), pq.end(), comp);
        std::pair<int,double> v = pq.back();
        pq.pop_back();
        if(visited[v.first])
        {
            // we have already visited that one
            // so we can abort the current iteration
            continue;
        }
        else
        {
            // mark it as visited
            visited[v.first] = true;
        }

        for(int w = 0; w < size; w++)
        {
            if(!visited[w])
            {
                double newDist = computeEdgeDist(v.first,w);
                if(newDist < dist[w])
                {
                    pred[w] = v.first;
                    pq.push_back({w,newDist});
                    // we also need to do the changes to the heap
                    std::push_heap(pq.begin(), pq.end(), comp);
                    dist[w] = newDist;
                }
            }
        }
    }
}

double MST::computeEdgeDist(int ind1, int ind2)
{
    return computeDistance(vertices[ind1],vertices[ind2]);
}

// computes the euclidian distance between two points
// no rounding is made as an adjustment for the christofides algorithm
double MST::computeDistance(std::array<double,2> v1, std::array<double,2> v2)
{
    double sum = (v1.at(0) - v2.at(0))*(v1.at(0) - v2.at(0));
    sum += (v1.at(1) - v2.at(1))*(v1.at(1) - v2.at(1));
    sum = std::sqrt(sum);
    return sum;
}


int main()
{
    MST *m = new MST();
    int iters;
    std::cin >> iters;
    std::vector<std::array<double,2>> v(iters);
    read_in_stdin(v);
    return 0;
}