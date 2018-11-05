#include <iostream>
#include <cmath>
#include <stdio.h>
#include "perfect_matching.hpp"
#include "compute_distance.hpp"

/**
 * Using the implementation of Blossom V by Vladimir Kolmogorov.
 *
 * Vladimir Kolmogorov. "Blossom V: A new implementation of a minimum cost perfect matching algorithm."
 * In Mathematical Programming Computation (MPC), July 2009, 1(1):43-67.
 */
#include "Blossom5/PerfectMatching.h"

/**
 *  find_perfect_matching
 *  @param graph Vector of [v1,v2] edge pairs (vertex indices)
 *  @param vertices Vector of [x,y] coordinates for all vertices (full graph)
 *  @param vertexCount The number of vertices for graph
 *  @param matching The matching to be calculated by reference
 */ 
void find_perfect_matching(std::vector<std::array<int,2>> &graph, std::vector<std::array<double,2>> &vertices, unsigned int vertexCount, std::vector<std::array<int,2>> &matching) {

    unsigned int numEdges = graph.size();

    PerfectMatching *pm = new PerfectMatching(vertexCount, numEdges);
    for (auto &edge : graph) {
        double distance = compute_distance(vertices[edge[0]], vertices[edge[1]]);
        pm->AddEdge(edge[0], edge[1], distance);
    }

    pm->Solve();
    for (unsigned int e = 0; e < numEdges; e++) {
        if (pm->GetSolution(e) == 1) {
            matching.push_back(graph[e]);
        }
    }

    delete pm;
}