#pragma once

#include <vector>
#include <list>
#include <Eigen/Core>

using namespace std;
using namespace Eigen;

class Adjacency {
 public:
    Adjacency(const MatrixXi &F, int _nVertices);
    ~Adjacency();

    list<int> getNeighborhood(int s, int maxDepth = 3);

 private:
    int nVertices;
    vector<list<int> > adjList;
    
};
