#include "adjacency.h"

#include <iterator>
#include <Eigen/SparseCore>

using namespace std;
using namespace Eigen;

Adjacency::~Adjacency(){}

Adjacency::Adjacency(const MatrixXi &F, int nVertices) {
    adjList = vector<list<int> >(nVertices);

    vector<Triplet<bool> > triplets;
    for (int i = 0; i < F.rows(); i++) {
	for (int j = 0; j < 3; j++) {
	    int u = F(i, j);
	    int v = F(i, (j+1)%3);
	    triplets.push_back(Triplet<bool>(u, v, true));
	    triplets.push_back(Triplet<bool>(v, u, true));
	}
    }
    SparseMatrix<bool> adjMatrix(nVertices, nVertices);
    adjMatrix.setFromTriplets(triplets.begin(), triplets.end());

    for (int k = 0; k < adjMatrix.outerSize(); k++) {
	for (SparseMatrix<bool>::InnerIterator it(adjMatrix, k); it; ++it) {
	    adjList[it.row()].push_back(it.col());
	}
    }
}
