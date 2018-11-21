//#include "euler_hamilton.hpp"
#include "compute_distance.hpp"
#include "k-opt.hpp"
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

void threeOpt(std::vector<int> &circuit, std::vector<std::array<double,2>> &vertices, clock_t &init_clock)
{

	double current;
	bool flag = 0;
	for (int i = 0; i < circuit.size(); i++)
	{
		for (int j = i + 2; j < circuit.size(); j++)
		{
			for (int k = j + 2; k < circuit.size(); k++)
			{
                int casee = 0;
                std::vector<int> h;
				current = compute_rounded_distance(vertices[circuit[i]],vertices[circuit[i+1]])+ compute_rounded_distance(vertices[circuit[j]],vertices[circuit[j + 1]])+ compute_rounded_distance(vertices[circuit[k]],vertices[circuit[k + 1]]);
				if (current > compute_rounded_distance(vertices[circuit[i]],vertices[circuit[i + 1]])+ compute_rounded_distance(vertices[circuit[j]],vertices[circuit[k]])+ compute_rounded_distance(vertices[circuit[j + 1]],vertices[circuit[k + 1]]))
				{
					current = compute_rounded_distance(vertices[circuit[i]],vertices[circuit[i + 1]])+ compute_rounded_distance(vertices[circuit[j]],vertices[circuit[k]])+ compute_rounded_distance(vertices[circuit[j + 1]],vertices[circuit[k + 1]]);
					casee = 1;
				}

				if (current > compute_rounded_distance(vertices[circuit[i]],vertices[circuit[j]])+ compute_rounded_distance(vertices[circuit[i + 1]],vertices[circuit[j + 1]])+ compute_rounded_distance(vertices[circuit[k]],vertices[circuit[k + 1]]))
				{
					current = compute_rounded_distance(vertices[circuit[i]],vertices[circuit[j]])+ compute_rounded_distance(vertices[circuit[i + 1]],vertices[circuit[j + 1]])+ compute_rounded_distance(vertices[circuit[k]],vertices[circuit[k + 1]]);
					casee = 2;
				}

				if (current > compute_rounded_distance(vertices[circuit[i]],vertices[circuit[j]])+ compute_rounded_distance(vertices[circuit[i + 1]],vertices[circuit[k]])+ compute_rounded_distance(vertices[circuit[j + 1]],vertices[circuit[k + 1]]))
				{
					current = compute_rounded_distance(vertices[circuit[i]],vertices[circuit[j]])+ compute_rounded_distance(vertices[circuit[i + 1]],vertices[circuit[k]])+ compute_rounded_distance(vertices[circuit[j + 1]],vertices[circuit[k + 1]]);
					casee = 3;
				}

				if (current > compute_rounded_distance(vertices[circuit[i]],vertices[circuit[j + 1]])+ compute_rounded_distance(vertices[circuit[k]],vertices[circuit[i + 1]])+ compute_rounded_distance(vertices[circuit[j]],vertices[circuit[k + 1]]))
				{
					current = compute_rounded_distance(vertices[circuit[i]],vertices[circuit[j + 1]])+ compute_rounded_distance(vertices[circuit[k]],vertices[circuit[i + 1]])+ compute_rounded_distance(vertices[circuit[j]],vertices[circuit[k + 1]]);
					casee = 4;
				}

				if (current > compute_rounded_distance(vertices[circuit[i]],vertices[circuit[j + 1]])+ compute_rounded_distance(vertices[circuit[k]],vertices[circuit[j]])+ compute_rounded_distance(vertices[circuit[i + 1]],vertices[circuit[k + 1]]))
				{
					current = compute_rounded_distance(vertices[circuit[i]],vertices[circuit[j + 1]])+ compute_rounded_distance(vertices[circuit[k]],vertices[circuit[j]])+ compute_rounded_distance(vertices[circuit[i + 1]],vertices[circuit[k + 1]]);
					casee = 5;
				}

				if (current > compute_rounded_distance(vertices[circuit[i]],vertices[circuit[k]])+ compute_rounded_distance(vertices[circuit[j + 1]],vertices[circuit[i + 1]])+ compute_rounded_distance(vertices[circuit[j]],vertices[circuit[k + 1]]))
				{
					current = compute_rounded_distance(vertices[circuit[i]],vertices[circuit[k]])+ compute_rounded_distance(vertices[circuit[j + 1]],vertices[circuit[i + 1]])+ compute_rounded_distance(vertices[circuit[j]],vertices[circuit[k + 1]]);
					casee = 6;
				}

				if (current > compute_rounded_distance(vertices[circuit[i]],vertices[circuit[k]])+ compute_rounded_distance(vertices[circuit[j + 1]],vertices[circuit[j]])+ compute_rounded_distance(vertices[circuit[i + 1]],vertices[circuit[k + 1]]))
				{
					current = compute_rounded_distance(vertices[circuit[i]],vertices[circuit[k]])+ compute_rounded_distance(vertices[circuit[j + 1]],vertices[circuit[j]])+ compute_rounded_distance(vertices[circuit[i + 1]],vertices[circuit[k + 1]]);
					casee = 7;
				}

				if(casee == 1)
				{
					for(int l=0;l<j+1;l++) h.push_back(circuit[l]);
					for(int l=k;l>=j+1;l--) h.push_back(circuit[l]);
					for(int l=k+1;l<circuit.size();l++) h.push_back(circuit[l]);
				}
				else if(casee == 2)
				{
					for(int l=0;l<i+1;l++) h.push_back(circuit[l]);
					for(int l=j;l>=i+1;l--) h.push_back(circuit[l]);
					for(int l=j+1;l<circuit.size();l++) h.push_back(circuit[l]);
				}
				else if(casee == 3)
				{
					for(int l=0;l<i+1;l++) h.push_back(circuit[l]);
					for(int l=j+1;l>=i+1;l--) h.push_back(circuit[l]);
					for(int l=k;l>j+1;l--) h.push_back(circuit[l]);
					for(int l=k+1;l<circuit.size();l++) h.push_back(circuit[l]);
				}
				else if(casee == 4)
				{
					for(int l=0;l<i+1;l++) h.push_back(circuit[l]);
					for(int l=j+1;l<k+1;l++) h.push_back(circuit[l]);
					for(int l=i+1;l<j+1;l++) h.push_back(circuit[l]);
					for(int l=k+1;l<circuit.size();l++) h.push_back(circuit[l]);
				}
				else if(casee == 5)
				{
					for(int l=0;l<i+1;l++) h.push_back(circuit[l]);
					for(int l=j+1;l<k+1;l++) h.push_back(circuit[l]);
					for(int l=j;l>=i+1;l--) h.push_back(circuit[l]);
					for(int l=k+1;l<circuit.size();l++) h.push_back(circuit[l]);
				}
				else if(casee == 6)
				{
					for(int l=0;l<i+1;l++) h.push_back(circuit[l]);
					for(int l=k;l>=j+1;l--) h.push_back(circuit[l]);
					for(int l=i+1;l<j+1;l++) h.push_back(circuit[l]);
					for(int l=k+1;l<circuit.size();l++) h.push_back(circuit[l]);
				}
				else if(casee == 7)
				{
					for(int l=0;l<i+1;l++) h.push_back(circuit[l]);
					for(int l=k;l>j+1;l--) h.push_back(circuit[l]);
					for(int l=j+1;l>=i+1;l--) h.push_back(circuit[l]);
					for(int l=k+1;l<circuit.size();l++) h.push_back(circuit[l]);
				}
				for(int l=0;l<h.size();l++)
				circuit[l]=h[l];


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
		if(flag)
		{
			break;
		}
	}
}