#include <cmath>
#include "compute_distance.hpp"

double compute_distance(std::array<double,2> v1, std::array<double,2> v2) {
    double sum = (v1.at(0) - v2.at(0)) * (v1.at(0) - v2.at(0));
    sum += (v1.at(1) - v2.at(1)) * (v1.at(1) - v2.at(1));
    sum = std::sqrt(sum);
    return sum;
}
