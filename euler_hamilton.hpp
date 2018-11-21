#include <iostream>
#include <unordered_map>
#include <vector>

// Seperate chaining hash table to hold graph. May change to a different data structure later.
void changeToMap(std::vector<std::array<int,2>> &original_graph, std::unordered_map<int,std::vector<int>> &new_graph);

// Returns 1 if graph has an Eulerian circuit, 0 if it does not.
bool hasEulerianCircuit(std::unordered_map<int,std::vector<int>> &graph);

void findEulerianCircuit(std::unordered_map<int,std::vector<int>> &input_graph, std::vector<int> &eulerian_path, int node);

void findHamiltonianCircuit(std::vector<int> &hamiltonian_circuit, int n);

double findTotalCost(std::vector<int> &hamiltonian_circuit, std::vector<std::array<double,2>> &vertices);