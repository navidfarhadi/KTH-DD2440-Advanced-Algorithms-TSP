#include <iostream>
#include <vector>
#include <array>
#include <algorithm>
#include <cstdlib>
#include <tuple>
#include <limits>
#include <cmath>

/*
 * This is the file that implements the minimum spanning tree
 * algorithm. We use the algorithm of Prim
 */

class MST
{
    public:
        std::vector<std::array<double,2>> vertices;

        MST(std::vector<std::array<double,2>> &vertices);
        void getMST(std::vector<std::array<int,2>> &e);
        double computeEdgeDist(int ind1, int ind2);
        void print_pq(std::vector<std::pair<int,double>> &pq);
};
