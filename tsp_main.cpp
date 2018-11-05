#include <iostream>

#include "read.hpp"
#include "mst.hpp"
#include "odd_degree_subgraph.hpp"
#include "perfect_matching.hpp"


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

    std::cout << "Perfect matching edges:\n";
    for (auto &edge : matching_edges) {
        std::cout << edge[0] << " - " << edge[1] << "\n";
    }

    std::vector<std::array<int,2>> multigraph;
    multigraph.reserve(mst_edges.size() + matching_edges.size());
    multigraph.insert(multigraph.end(), mst_edges.begin(), mst_edges.end());
    multigraph.insert(multigraph.end(), matching_edges.begin(), matching_edges.end());


    return 0;
}