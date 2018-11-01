#include <iostream>
#include <vector>
#include <array>
#include <algorithm>

/*
 * This is the file that implements the minimum spanning tree
 * algorithm. We use the algorithm of Kruskal
 */

class MST
{
    public:
        MST();
        std::vector<std::array<int,2>> getMST();
        std::vector<int> parent;
        std::vector<std::array<double,2>> vertices;
        void sortEdges(std::vector<std::array<int,2>> &e);
        int comparator(std::vector<std::array<int,2>> &e1, std::vector<std::array<int,2>> &e2);

};

MST::MST()
{
    std::cout << "MST\n";
}

void MST::sortEdges(std::vector<std::array<int,2>> &e)
{
    //std::sort(e.begin,e.end,comparator);
}

int MST::comparator(std::vector<std::array<int,2>> &e1, std::vector<std::array<int,2>> &e2)
{
    return 0;
}

int main()
{
    MST *m = new MST();
    m->parent[0];
}