#include "mst.hpp"
#include "compute_distance.hpp"


MST::MST(std::vector<std::array<double,2>> &vertices)
{
    std::cerr << "MST\n";
    this->vertices = vertices;
}

struct compare
{
    bool operator () (const std::pair<int,double> &lhs, const std::pair<int,double> &rhs)
    {
        return (lhs.second > rhs.second);
    }
};

void MST::getMST(std::vector<std::array<int,2>> &e)
{
    std::cerr << "start getMST\n";
    int size = vertices.size();
    double dist[size];
    int pred[size];
    bool visited[size] = {false};
    compare comp;
    // we always start with vertic index 0
    // because it does not make a difference

    // distance zero
    dist[0] = 0;
    // all other have max distance
    for(int i = 1; i < size; i++)
    {
        dist[i] = __DBL_MAX__;
    }
    // -2: root node
    // -1: no predecessor up to now
    pred[0] = -2;
    for(int i = 1; i < size; i++)
    {
        pred[i] = -1;
    }
    // boolean value stands for if it is invalid
    // true: invalid
    // false: valid
    std::vector<std::pair<int,double>> pq;

    pq.push_back({0,0});

    std::make_heap(pq.begin(),pq.end(),comp);

    while (!pq.empty())
    {
        print_pq(pq);
        std::pop_heap(pq.begin(), pq.end(), comp);
        std::pair<int,double> v = pq.back();
        pq.pop_back();
        if(visited[v.first])
        {
            // we have already visited that one
            // so we can abort the current iteration
            continue;
        }
        else
        {
            // mark it as visited
            visited[v.first] = true;
        }

        std::cerr << "EVALUATING: " << v.first << "\n";

        for(int w = 0; w < size; w++)
        {
            if(!visited[w])
            {
                double newDist = computeEdgeDist(v.first,w);
                if(newDist < dist[w])
                {
                    pred[w] = v.first;
                    pq.push_back({w,newDist});
                    // we also need to do the changes to the heap
                    std::push_heap(pq.begin(), pq.end(), comp);
                    dist[w] = newDist;
                }
            }
        }
    }
    for(int i = 0; i < size; i++)
    {
        std::cerr << "pred[" << i << "] == " << pred[i] << "\n";
        if(pred[i] >= 0)
        {
            // then it is a valid predecessor
            e.push_back({pred[i],i});
        }
    }
}

double MST::computeEdgeDist(int ind1, int ind2)
{
    return compute_distance(vertices[ind1],vertices[ind2]);
}

void MST::print_pq(std::vector<std::pair<int,double>> &pq)
{
    std::cerr << "print pq\n";
    for(auto &line : pq)
    {
        std::cerr << line.first << " dist: " << line.second << "\n";
    }
}
