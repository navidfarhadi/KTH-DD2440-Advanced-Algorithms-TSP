#ifndef TSP_COMPUTE_DISTANCE_HPP
#define TSP_COMPUTE_DISTANCE_HPP

#include <array>

// computes the euclidian distance between two points
// no rounding is made as an adjustment for the christofides algorithm
double compute_distance(std::array<double,2> v1, std::array<double,2> v2);

#endif //TSP_COMPUTE_DISTANCE_HPP
