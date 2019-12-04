#include "clustering.h"

#include <Eigen/Eigenvalues>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <limits>

using namespace std;
using namespace Eigen;

SpectralClustering::SpectralClustering(const MatrixXd &X, const MatrixXi &F, int _k, int _nFeatures) {
    buildLB(X, F);

    nFeatures = _nFeatures;
    computeDescriptors();

    k = _k;
    computeClusters();
}

/*
  Use Laplace-Beltrami operator to compute a descriptor for each vertex.
*/

void SpectralClustering::buildLB(const MatrixXd &X, const MatrixXi &F) {
    int n = X.rows();
    int f = F.rows();

    A = MatrixXd::Zero(n, n);
    W = MatrixXd::Zero(n, n);
    for (int i = 0; i < f; i++) {
	int u, v, w;
	u = F(i, 0); v = F(i, 1); w = F(i, 2);
	
	double d0, d1, d2;
	d0 = (X.row(v) - X.row(u)).norm();
	d1 = (X.row(w) - X.row(v)).norm();
	d2 = (X.row(u) - X.row(w)).norm();

	double s = (d0 + d1 + d2) / 2;
	double area = sqrt(s * (s - d0) * (s - d1) * (s - d2));
	
	A(u, u) += area / 3; A(v, v) += area / 3; A(w, w) += area / 3;

	double w0, w1, w2;
	w0 = (d0*d0 - d1*d1 - d2*d2) / (8 * area);
	w1 = (d1*d1 - d2*d2 - d0*d0) / (8 * area);
	w2 = (d2*d2 - d0*d0 - d1*d1) / (8 * area);

	W(u, v) += w0; W(v, u) += w0; W(u, u) -= w0; W(v, v) -= w0;
	W(v, w) += w1; W(w, v) += w1; W(v, v) -= w1; W(w, w) -= w1;
	W(w, u) += w2; W(u, w) += w2; W(w, w) -= w2; W(u, u) -= w2;
    }
}

void SpectralClustering::computeDescriptors() {
    GeneralizedSelfAdjointEigenSolver<MatrixXd> es(W, A);
    D = es.eigenvectors().block(0, 1, es.eigenvectors().rows(), nFeatures);

    cout << "Eigenvalues of LB operator used in features :" << endl;
    cout << es.eigenvalues().segment(1, nFeatures) << endl;
}

/*
  Perform k-means clustering on descriptor dataset.
*/

void initCentroids(MatrixXd &centroids, const MatrixXd &D) {
    int k = centroids.rows();
    int n = D.rows();

    srand(time(nullptr));
    int r = rand() % k;
    centroids.row(0) = D.row(r);

    for (int j = 1; j < k; j++) {
	// squared distance to closest centroid
	vector<double> sqr_dist(n, numeric_limits<double>::max());
	double sum_dist = 0.;
	for (int i = 0; i < n; i++) {
	    for (int c = 0; c < j; c++) {
		double cur_dist = (centroids.row(c) - D.row(i)).squaredNorm();
		sqr_dist[i] = min(sqr_dist[i], cur_dist);
	    }
	    sum_dist += sqr_dist[i];
	}

	// new centroid chosen randomly using squared distance distribution
	double u = ((double) rand() / RAND_MAX) * sum_dist;
	int r = -1; double cur_sum = 0.;
	while (cur_sum <= u) {
	    r++;
	    cur_sum += sqr_dist[r];
	}
	centroids.row(j) = D.row(r);
    }
}

void updatePartition(vector<int> &partition, const MatrixXd &centroids, const MatrixXd &D) {
    int k = centroids.rows();
    int n = D.rows();

    for (int i = 0; i < n; i++) {
	double min_dist = numeric_limits<double>::max();
	int idx = -1;
	for (int j = 0; j < k; j++) {
	    double cur_dist = (centroids.row(j) - D.row(i)).norm();
	    if (cur_dist < min_dist) {
		min_dist = cur_dist;
		idx = j;
	    }
	}
	partition[i] = idx;
    }
}

void moveCentroids(MatrixXd &centroids, const vector<int> &partition, const MatrixXd &D) {
    int k = centroids.rows();
    int n = D.rows();
    int nFeatures = D.cols();

    vector<int> cluster_size(k, 0);
    MatrixXd vsum = MatrixXd::Zero(k, nFeatures);
    for (int i = 0; i < n; i++) {
	int j = partition[i];
	cluster_size[j]++;
	vsum.row(j) += D.row(i);
    }

    for (int j = 0; j < k; j++) {
	if (cluster_size[j] > 0) {
	    centroids.row(j) = vsum.row(j) / cluster_size[j];
	}
    }
}

void showEnergy(const vector<int> &partition, const MatrixXd &centroids, const MatrixXd &D) {
    int n = D.rows();
    
    double energy = 0.;
    for (int i = 0; i < n; i++) {
	energy += (centroids.row(partition[i]) - D.row(i)).squaredNorm() / n;
    }

    cout << "energy : " << energy << endl;
}

void SpectralClustering::computeClusters(int maxiter) {
    clusters = vector<list<int> >(k);

    MatrixXd centroids(k, nFeatures);
    initCentroids(centroids, D);

    int n = D.rows();
    vector<int> partition(n, 0);
    updatePartition(partition, centroids, D);
    showEnergy(partition, centroids, D);

    for (int i = 0; i < maxiter; i++) {
	moveCentroids(centroids, partition, D);
	updatePartition(partition, centroids, D);
	if (i % 10 == 0) showEnergy(partition, centroids, D);
    }

    for (int i = 0; i < n; i++) {
	clusters[partition[i]].push_back(i);
    }
}
