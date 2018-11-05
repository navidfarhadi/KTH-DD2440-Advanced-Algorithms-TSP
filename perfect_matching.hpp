#ifndef PERFECT_MATCHING
#define PERFECT_MATCHING

#include <vector>
#include <array>

void find_perfect_matching(std::vector <std::array<int, 2>> &graph, std::vector <std::array<double, 2>> &vertices,
                           unsigned int vertexCount, std::vector <std::array<int, 2>> &matching);

#endif