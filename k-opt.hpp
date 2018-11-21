#ifndef TWOOPT_HPP
#define TWOOPT_HPP

#include <vector>
#include <iostream>
#include <array>

void twoOpt(std::vector<int> &circuit, std::vector<std::array<double,2>> &vertices, clock_t &init_clock);
void twoHOpt(std::vector<int> &circuit, std::vector<std::array<double,2>> &vertices, clock_t &init_clock);
void threeOpt(std::vector<int> &circuit, std::vector<std::array<double,2>> &vertices, clock_t &init_clock);

#endif //TWOOPT__HPP

