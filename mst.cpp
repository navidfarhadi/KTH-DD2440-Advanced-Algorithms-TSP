#include <iostream>
#include <vector>
#include <array>
#include <algorithm>
#include <cstdlib>
#include <tuple>
#include <limits>
#include <cmath>
#include "read.hpp"
#include "odd_degree_subgraph.hpp"

/*
 * This is the file that implements the minimum spanning tree
 * algorithm. We use the algorithm of Prim
 */

class MST
{
    public:
        MST(std::vector<std::array<double,2>> &vertices);
        void getMST(std::vector<std::array<int,2>> &e);
        std::vector<std::array<double,2>> vertices;
        double computeEdgeDist(int ind1, int ind2);
        double computeDistance(std::array<double,2> v1, std::array<double,2> v2);
        void print_pq(std::vector<std::pair<int,double>> &pq);
};

MST::MST(std::vector<std::array<double,2>> &vertices)
{
    std::cerr << "MST\n";
    this->vertices = vertices;
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
    std::cerr << "start getMST\n";
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
        print_pq(pq);
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

        std::cerr << "EVALUATING: " << v.first << "\n";

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
    for(int i = 0; i < size; i++)
    {
        std::cerr << "pred[" << i << "] == " << pred[i] << "\n"; 
        if(pred[i] >= 0)
        {
            // then it is a valid predecessor
            e.push_back({pred[i],i});
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

void MST::print_pq(std::vector<std::pair<int,double>> &pq)
{
    std::cerr << "print pq\n";
    for(auto &line : pq)
    {
        std::cerr << line.first << " dist: " << line.second << "\n";
    }
}


int main()
{
    int iters;
    std::cin >> iters;
    std::vector<std::array<double,2>> v(iters);
    read_in_stdin(v);
    std::vector<std::array<int,2>> mst_edges;
    MST *m = new MST(v);
    m->getMST(mst_edges);
    for(auto &line : mst_edges)
    {
        std::cout << line[0] << " - " << line[1] << "\n";    
    }
    std::vector<int> vers;
    std::vector<std::array<int,2>> sub_edges;
    find_odd_v(mst_edges,vers,iters);
    create_subgraph(vers, sub_edges);
    return 0;
}