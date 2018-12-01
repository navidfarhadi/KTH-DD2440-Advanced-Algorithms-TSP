#include<iostream>
#include <vector>
#include <map>
#include <unordered_map>
#include <climits>
#include "read.hpp"
#include "mst.hpp"
#include "odd_degree_subgraph.hpp"
#include "perfect_matching.hpp"
#include "euler_hamilton.hpp"
#include "k-opt.hpp"
#include "compute_distance.hpp"


int main()
{   
    clock_t init_clock = clock();

    int numVertices;
    std::cin >> numVertices;
    std::vector<std::array<double,2>> vertices(numVertices);
    read_in_stdin(vertices);

    //initialization of the array of distances
    init_dist_array(vertices);

    //Implementation of Christofides algorithm

    //Minimum spanning tree
    std::vector<std::array<int,2>> mst_edges;
    MST *m = new MST(vertices);
    m->getMST(mst_edges);

    //Finding vertices with Odd Degree
    std::vector<int> oddVertexIndices;
    std::vector<std::array<int,2>> sub_edges;
    find_odd_v(mst_edges,oddVertexIndices,numVertices);
    create_subgraph(oddVertexIndices, sub_edges);

    //Perfect matching
    std::vector<std::array<int,2>> matching_edges;
    find_perfect_matching(sub_edges, vertices, oddVertexIndices, matching_edges);

    std::vector<std::array<int,2>> multigraph;
    multigraph.reserve(mst_edges.size() + matching_edges.size());
    multigraph.insert(multigraph.end(), mst_edges.begin(), mst_edges.end());
    multigraph.insert(multigraph.end(), matching_edges.begin(), matching_edges.end());

    std::unordered_map<int,std::vector<int>> new_graph;
    changeToMap(multigraph, new_graph);

    if(hasEulerianCircuit(new_graph))
    {
	//Finding eularian and hamiltonian circuit
        std::vector<int> path;
        findEulerianCircuit(new_graph,path, 0);
        findHamiltonianCircuit(path,numVertices);

	//Optimization k-opt

        while(double(clock() - init_clock) / CLOCKS_PER_SEC < 1.05)
        {
	    twoOpt(path, vertices, init_clock);

            twoHOpt(path, vertices, init_clock);
        }

        while(double(clock() - init_clock) / CLOCKS_PER_SEC < 1.99)
        {
            twoHOpt(path, vertices, init_clock);
            
            threeOpt(path, vertices, init_clock);
        }

        for(int i = 0; i < path.size(); i++)
        {
            std::cout << path[i] << std::endl;
        }
    }

    return 0;
}
