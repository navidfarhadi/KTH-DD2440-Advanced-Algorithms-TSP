#include <iostream>
#include <fstream>
#include <cmath>
#include <map>
#include <stdio.h>
#include <unistd.h>
#include <sys/stat.h> 
#include <fcntl.h>
#include "perfect_matching.hpp"
#include "compute_distance.hpp"

/**
 * Using the implementation of Blossom V by Vladimir Kolmogorov.
 *
 * Vladimir Kolmogorov. "Blossom V: A new implementation of a minimum cost perfect matching algorithm."
 * In Mathematical Programming Computation (MPC), July 2009, 1(1):43-67.
 */
#include "Blossom5/PerfectMatching.h"

int suppress_stdout() 
{
    fflush(stdout);
    int ret = dup(1);
    int nullfd = open("/dev/null", O_WRONLY);
    dup2(nullfd, 1);
    close(nullfd);
    return ret;
}

void resume_stdout(int fd) 
{
    fflush(stdout);
    dup2(fd, 1);
    close(fd);
}

/**
 *  find_perfect_matching
 *  @param graph Vector of [v1,v2] edge pairs (vertex indices)
 *  @param vertices Vector of [x,y] coordinates for all vertices (full graph)
 *  @param usedVertices Vector of vertex indices used in graph
 *  @param matching The matching to be calculated by reference
 */ 
void find_perfect_matching(std::vector<std::array<int,2>> &graph, std::vector<std::array<double,2>> &vertices, std::vector<int> &usedVertices, std::vector<std::array<int,2>> &matching) {

    // std::streambuf *original_cout_buffer = std::cout.rdbuf();
    // std::ofstream null_stream("/dev/null");
    // std::cout.rdbuf(null_stream.rdbuf());

    int suppress = suppress_stdout();
    
    unsigned int numVertices = usedVertices.size();
    unsigned int numEdges = graph.size();

    PerfectMatching *pm = new PerfectMatching(numVertices, numEdges);
    // std::cout << "Finding perfect matching, edges:\n";

    // Can't use arbitrary vertex indices, as the implementation
    // needs all sequential 0..N (N = numVertices) indices, therefore
    // we make a mapping.
    std::map<int, int> sequentialIndices;
    for (unsigned int vertex = 0; vertex < usedVertices.size(); vertex++) {
        sequentialIndices.insert({usedVertices.at(vertex), vertex});
    }

    for (auto &edge : graph) {
        double distance = compute_distance(vertices[edge[0]], vertices[edge[1]]);

        int vertA = sequentialIndices.find(edge[0])->second;
        int vertB = sequentialIndices.find(edge[1])->second;
        
        // std::cout << vertA << " -> " << vertB << " = " << distance << std::endl;
        pm->AddEdge(vertA, vertB, distance);
    }

    pm->Solve();
    for (unsigned int e = 0; e < numEdges; e++) {
        if (pm->GetSolution(e) == 1) {
            matching.push_back(graph[e]);
        }
    }

    delete pm;

    resume_stdout(suppress);
}

