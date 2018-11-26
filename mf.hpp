/*
 * This is the hpp file for the implementation
 * of the Multiple Fragment Heuristic approach.
 * The goal is to produce an initial tour in
 * a better time than Christofides, but still
 * have low total cost.
 */

#include "compute_distance.hpp"

void create_inital_tour(int numVertices, std::vector<int> &init_tour);