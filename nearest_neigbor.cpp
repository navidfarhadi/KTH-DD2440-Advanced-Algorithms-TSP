/*
 * This is the cpp file for the implementation
 * of the nearest neighbor algorithm
 */

#include "nearest_neighbor.hpp"
#include "compute_distance.hpp"

void compute_nearest_neighbor(std::vector<int> &path, std::vector<std::array<double,2>> &v)
{
    // we always start with node #0
    // can be changed if needed
    int size = v.size();
    bool visited[size] = {false};
    // the current node
    int cur = 0;
    int counter = 1;
    visited[cur] = true;
    path.push_back(cur);
    while(true)
    {
        double dist = __DBL_MAX__;
        int index = -1;
        for(int i = 0; i < size; i++)
        {
            if(!visited[i])
            {
                double cur_dist = compute_distance(v[cur],v[i]);
                if(cur_dist < dist)
                {
                    dist = cur_dist;
                    index = i;
                }
            }
        }
        cur = index;
        counter++;
        visited[cur] = true;
        path.push_back(cur);
        if(counter == size)
        {
            break;
        }
    }
}