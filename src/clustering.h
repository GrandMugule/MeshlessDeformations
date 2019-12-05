#pragma once

#include <vector>
#include <list>
#include <Eigen/Core>

using namespace std;
using namespace Eigen;

class SpectralClustering {
 public:
    SpectralClustering(const MatrixXd &X, const MatrixXi &F, int _k, int _nFeatures = 5);
    ~SpectralClustering(){}

    vector<list<int> >& getClusters(){return clusters;}

 private:
    // Laplace-Beltrami operator
    MatrixXd A;
    MatrixXd W;
    void buildLB(const MatrixXd &X, const MatrixXi &F);

    // solve generalized eigen problem
    // keep nFeatures first eigenvectors
    // and use their coordinates as descriptors
    int nFeatures;
    MatrixXd D;
    void computeDescriptors();

    // perform k-means clustering on descriptors
    int k;
    vector<list<int> > clusters;
    void computeClusters(int maxiter = 50);
    
};
