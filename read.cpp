#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <array>

void read_in_stdin(std::vector<std::array<double,2>> &v)
{
    for(auto &line : v)
    {
        for(auto &d : line)
        {
            std::cin >> d;
            std::cout << d << " ";
        }
        std::cout << "\n";
    }
    std::cout << (double)(v[1][1]) << "\n";
}

// computes the euclidian distance between two points
// no rounding is made as an adjustment for the christofides algorithm
double compute_distance(std::array<double,2> v1, std::array<double,2> v2)
{
    double sum = (v1.at(0) - v2.at(0))*(v1.at(0) - v2.at(0));
    sum += (v1.at(1) - v2.at(1))*(v1.at(1) - v2.at(1));
    sum = std::sqrt(sum);
    return sum;
}

/*int main()
{
    int iters;
    std::cin >> iters;
    std::vector<std::array<double,2>> v(iters);
    read_in_stdin(v);

    std::cout << (double)(v[1][1]) << "\n";
    return 0;
}*/