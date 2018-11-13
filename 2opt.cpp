#include "euler_hamilton.hpp"
#include "2opt.hpp"
#include <iostream>
#include <algorithm>

void twoOpt(std::vector<int> &circuit, std::vector<std::array<double,2>> &vertices, clock_t &init_clock)
{
    double best_distance = findTotalCost(circuit, vertices);
    double new_distance;
    std::vector<int> new_circuit(circuit);

	clock_t begin = clock();
	bool flag = 0;

    for (int i=0; i<circuit.size()-1; i++)
    {
		for(int j=i+1; j<circuit.size(); j++)
		{
			new_circuit = circuit;
			std::reverse(std::begin(new_circuit)+i, std::begin(new_circuit)+j);
			new_distance = findTotalCost(new_circuit, vertices);
			if (new_distance < best_distance)
			{
				circuit = new_circuit;
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
