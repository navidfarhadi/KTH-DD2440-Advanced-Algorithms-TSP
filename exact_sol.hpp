/*
 * This is the hpp file for the implementation of the
 * exact solution of the TSP.
 * This implementation obviously is of a runtime NP
 * and was only for TESTING reasons implemented.
 */

#include <vector>
#include <array>

double computeTotalDist(int array_perms[], int size, std::vector<std::array<double,2>> &v);

double compute_best_distance(int size, std::vector<std::array<double,2>> &v);