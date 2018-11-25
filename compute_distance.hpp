#ifndef TSP_COMPUTE_DISTANCE_HPP
#define TSP_COMPUTE_DISTANCE_HPP

#include <array>
#include <vector>

// computes the euclidian distance between two points
// no rounding is made as an adjustment for the christofides algorithm


void init_dist_array(std::vector<std::array<double,2>> v);

double compute_distance(std::array<double,2> v1, std::array<double,2> v2);

int compute_rounded_distance(int i, int j);

#endif //TSP_COMPUTE_DISTANCE_HPP
