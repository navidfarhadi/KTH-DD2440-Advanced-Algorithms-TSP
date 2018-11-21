#include<iostream>
#include <vector>
#include <map>
#include <climits>
#include "read.hpp"
#include "mst.hpp"
#include "odd_degree_subgraph.hpp"
#include "perfect_matching.hpp"
#include "euler_hamilton.hpp"
#include "k-opt.hpp"


int main()
{   
    clock_t init_clock = clock();

    int numVertices;
    std::cin >> numVertices;
    std::vector<std::array<double,2>> vertices(numVertices);
    read_in_stdin(vertices);

    std::vector<std::array<int,2>> mst_edges;
    MST *m = new MST(vertices);
    m->getMST(mst_edges);

    std::vector<int> oddVertexIndices;
    std::vector<std::array<int,2>> sub_edges;
    find_odd_v(mst_edges,oddVertexIndices,numVertices);
    create_subgraph(oddVertexIndices, sub_edges);

    std::vector<std::array<int,2>> matching_edges;
    find_perfect_matching(sub_edges, vertices, oddVertexIndices, matching_edges);

    std::vector<std::array<int,2>> multigraph;
    multigraph.reserve(mst_edges.size() + matching_edges.size());
    multigraph.insert(multigraph.end(), mst_edges.begin(), mst_edges.end());
    multigraph.insert(multigraph.end(), matching_edges.begin(), matching_edges.end());

    std::map<int,std::vector<int>> new_graph;
    changeToMap(multigraph, new_graph);

    if(hasEulerianCircuit(new_graph))
    {
        std::vector<int> path;
        findEulerianCircuit(new_graph,path, 0);
        findHamiltonianCircuit(path,numVertices);

        double remaining = 2 - (double(clock() - init_clock) / CLOCKS_PER_SEC);
        double twoOptTime = remaining / 4;

        while(double(clock() - init_clock) / CLOCKS_PER_SEC < twoOptTime)
        {
            twoOpt(path, vertices, init_clock);
        }

        while(double(clock() - init_clock) / CLOCKS_PER_SEC < 1.99)
        {
            threeOpt(path, vertices, init_clock);
        }

        for(int i = 0; i < path.size(); i++)
        {
            std::cout << path[i] << std::endl;
        }
    }

    return 0;
}
