#pragma once

#include <vector>
#include <list>
#include <Eigen/Core>

using namespace std;
using namespace Eigen;

class Adjacency {
 public:
    Adjacency(const MatrixXi &F, int nVertices);
    ~Adjacency();

    list<int>& getNeighbours(int v){ return adjList[v]; } 

 private:
    vector<list<int> > adjList;
    
};
