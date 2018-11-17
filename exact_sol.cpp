/*
 * This is the hpp file for the implementation of the
 * exact solution of the TSP.
 * This implementation obviously is of a runtime NP
 * and was only for TESTING reasons implemented.
 */

#include "exact_sol.hpp"
#include <iostream>
#include <algorithm>
#include "compute_distance.hpp"

double computeTotalDist(int array_perms[], int size, std::vector<std::array<double,2>> &v)
{
    double totDist = 0.0;
    for(int i = 1; i < size; i++)
    {
        totDist += compute_distance(v.at(array_perms[i - 1]), v.at(array_perms[i]));
    }
    totDist += compute_distance(v.at(array_perms[size - 1]), v.at(array_perms[0]));
    return totDist;
}

double compute_best_distance(int size, std::vector<std::array<double,2>> &v)
{
    int num_perms = 1;
    int array_perms[size];
    //std::cout << "start compute_best_distance\n";
    for(int i = size; i > 0 ; i--)
    {
        num_perms *= i;
        array_perms[i-1] = i-1;
    }
    //std::cout << "num_perms = " << num_perms << "\n";

    double dist = __DBL_MAX__;

    do
    {
        double cur_dist = computeTotalDist(array_perms, size, v);
        if(cur_dist < dist)
        {
            dist = cur_dist;
        }
    } while(std::next_permutation(array_perms, array_perms + size));

    return dist;

}