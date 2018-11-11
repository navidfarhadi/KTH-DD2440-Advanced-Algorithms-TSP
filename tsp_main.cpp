#include <iostream>
#include <vector>
#include <map>

#include "read.hpp"
#include "mst.hpp"
#include "odd_degree_subgraph.hpp"
#include "perfect_matching.hpp"
#include "euler_hamilton.hpp"


int main()
{
    int numVertices;
    std::cin >> numVertices;
    std::vector<std::array<double,2>> vertices(numVertices);
    read_in_stdin(vertices);

    //  Minimum Spanning Tree
    std::vector<std::array<int,2>> mst_edges;
    MST *m = new MST(vertices);
    m->getMST(mst_edges);

    //  Odd degree vertex subgraph
    std::vector<int> oddVertexIndices;
    std::vector<std::array<int,2>> sub_edges;
    find_odd_v(mst_edges,oddVertexIndices,numVertices);
    create_subgraph(oddVertexIndices, sub_edges);

    // Min weight perfect matching
    std::vector<std::array<int,2>> matching_edges;
    unsigned int subgraphNumVertices = oddVertexIndices.size();
    find_perfect_matching(sub_edges, vertices, subgraphNumVertices, matching_edges);

    std::vector<std::array<int,2>> multigraph;
    multigraph.reserve(mst_edges.size() + matching_edges.size());
    multigraph.insert(multigraph.end(), mst_edges.begin(), mst_edges.end());
    multigraph.insert(multigraph.end(), matching_edges.begin(), matching_edges.end());

    std::cout << "Perfect matching and MST union edges:\n";
    for (auto &edge : multigraph) {
        std::cout << edge[0] << " - " << edge[1] << "\n";
    }

    // Convert graph to a map of vectors
    std::map<int,std::vector<int>> new_graph;
    changeToMap(multigraph, new_graph);

    std::cout << "size of new graph: " << new_graph.size() << std::endl;

    for (const auto &p : new_graph)
    {
        std::cout << "Node:" << p.first << std::endl << "Edges to: ";
        for(const auto &q : p.second)
        {
            std::cout << q << " ";
        }
        std::cout << std::endl;
    }

    if(hasEulerianCircuit(new_graph))
    {
        std::cout << "Graph has Eulerian circuit" << std::endl;

        std::vector<int> eulerian_path;
        findEulerianCircuit(new_graph,eulerian_path);

        std::cout << "Eulerian circuit: ";
        for(int i = 0; i < eulerian_path.size(); i++)
        {
            std::cout << eulerian_path[i] << " ";
        }
        std::cout << std::endl;
    }
    else
    {
        std::cout << "Graph does NOT have Eulerian circuit" << std::endl;
    }

    return 0;
}