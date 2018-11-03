#include <vector>
#include <array>


/*
 * hpp file for the implementation
 * of finding the vertices with an
 * odd degree in a MST and creating
 * the subgraph out of it.
 */

void find_odd_v(std::vector<std::array<int,2>> &e, std::vector<int> &v, int size);

void create_subgraph(std::vector<int> &v, std::vector<std::array<int,2>> &e);