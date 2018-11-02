#include <iostream>
#include <vector>
#include <array>
#include <algorithm>
#include <cstdlib>
#include <tuple>

/*
 * This is the file that implements the minimum spanning tree
 * algorithm. We use the algorithm of Prim
 */

class MST
{
    public:
        MST();
        void getMST(std::vector<std::array<int,2>> &e);
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

void MST::getMST(std::vector<std::array<int,2>> &e)
{
    int size = vertices.size();
    double dist[size];
    double pred[size];
    // we always start with vertic index 0
    // because it does not make a difference
    
    // distance zero
    dist[0] = 0;
    // no predecessor
    pred[0] = -1;

    // boolean value stands for if it is invalid
    // true: invalid
    // false: valid
    std::vector<std::tuple<int,double,bool>> pq;

    pq.push_back(std::make_tuple(0,0,false));

    std::make_heap(pq.begin(),pq.end());
}
int main()
{
    MST *m = new MST();
    m->parent[0];
}