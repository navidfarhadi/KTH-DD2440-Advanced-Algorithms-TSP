#include "euler_hamilton.hpp"
#include "compute_distance.hpp"
#include "2opt.hpp"
#include <iostream>
#include <algorithm>

void twoOpt(std::vector<int> &circuit, std::vector<std::array<double,2>> &vertices, clock_t &init_clock)
{
    int n = circuit.size();
    int k;
    int l;
    double delta;
    clock_t begin = clock();
    bool flag = 0;

    for (int i=0; i<circuit.size()-1; i++)
    {
		for(int j=i+1; j<circuit.size(); j++)
		{
			k = i-1;
			l = j+1;
			
			if (i==0)
			{
				k=n-1;
			}
			if (j==n-1)
			{
				l=0;
			}

			delta = compute_distance(vertices[circuit[k]], vertices[circuit[j]]) + compute_distance(vertices[circuit[i]], vertices[circuit[l]]) - compute_distance(vertices[circuit[k]], vertices[circuit[i]]) - compute_distance(vertices[circuit[j]], vertices[circuit[l]]);
			//delta = compute_rounded_distance(vertices[circuit[k]], vertices[circuit[j]]) + compute_rounded_distance(vertices[circuit[i]], vertices[circuit[l]]) - compute_rounded_distance(vertices[circuit[k]], vertices[circuit[i]]) - compute_rounded_distance(vertices[circuit[j]], vertices[circuit[l]]);

			if (delta<0)
			{
				std::reverse(std::begin(circuit)+i, std::begin(circuit)+j+1);
			}

			if(double(clock() - init_clock) / CLOCKS_PER_SEC > 1.99)
			{
				flag = 1;
				break;
			}
		}
		if(flag)
		{
			break;
		}
    }
}
