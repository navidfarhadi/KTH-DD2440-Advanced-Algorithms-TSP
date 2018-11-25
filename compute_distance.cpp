#include <cmath>
#include "compute_distance.hpp"
#include <iostream>

static int numVertices;
static int *dist_array;

void init_dist_array(std::vector<std::array<double,2>> v)
{
    numVertices = v.size();
    dist_array = new int[numVertices * numVertices];
    for(int i = 0; i < numVertices; i++)
    {
        for(int j = 0; j <= i; j++)
        {
            double sum = (v[i][0] - v[j][0]) * (v[i][0] - v[j][0]);
            sum += (v[i][1] - v[j][1]) * (v[i][1] - v[j][1]);
            dist_array[(i * numVertices) + j] = (int)(round(std::sqrt(sum)));
        }
    }
}

double compute_distance(std::array<double,2> v1, std::array<double,2> v2) {
    double sum = (v1.at(0) - v2.at(0)) * (v1.at(0) - v2.at(0));
    sum += (v1.at(1) - v2.at(1)) * (v1.at(1) - v2.at(1));
    sum = std::sqrt(sum);
    return sum;
}

int compute_rounded_distance(int i, int j) {
    /*double sum = (v1.at(0) - v2.at(0)) * (v1.at(0) - v2.at(0));
    sum += (v1.at(1) - v2.at(1)) * (v1.at(1) - v2.at(1));
    sum = std::sqrt(sum);
    sum = round(sum);
    return sum;*/
    if (i < j)
    {
        //std::cout << "access: j = " << j << " i = " << i << "\n";
        return  dist_array[(j * numVertices) + i];
    }
    else
    {
        //std::cout << "access: i = " << i << " j = " << j << "\n";
        return  dist_array[(i * numVertices) + j];
    }
}