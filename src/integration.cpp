#include "integration.h"
#include "shapematching.h"

#include <iterator>
#include <cassert>
#include <cmath>

using namespace std;
using namespace Eigen;


/*
  Constructor and setters.
*/

Integration::Integration(MatrixXd& _Xi, MatrixXd& _Xf, float _h, float _alpha){
    Xf = _Xf;
    h = _h;
    alpha = _alpha;
    features = {{Feature::GRAVITY, false},
		{Feature::CLUSTERS, false},
		{Feature::PLASTICITY, false}};

    int n = _Xi.rows();
    X = _Xi;
    V = MatrixXd::Zero(n, 3);
}

void Integration::addFeature(Feature f){
    features[f] = true;
}

void Integration::setGravity(float g){
    assert(features[Feature::GRAVITY]);
    
    gravity = g;
}

void Integration::setClusters(vector<list<int> >& _clusters){
    assert(features[Feature::CLUSTERS] || features[Feature::PLASTICITY]);
    assert(!_clusters.empty());
    
    clusters = _clusters;
    if (features[Feature::PLASTICITY]){
	S = vector<MatrixXd>(clusters.size(), MatrixXd::Identity(3, 3));
    }
}

void Integration::computeDestination(float beta){
    if (features[Feature::PLASTICITY]){
	MatrixXd Xp = X;
	for (int i = 0; i < clusters.size(); i++) {
	    MatrixXd Xc(clusters[i].size(), 3);
	    MatrixXd Xfc(clusters[i].size(), 3);
	    int j = 0;

	    for (list<int>::iterator v = clusters[i].begin(); v != clusters[i].end(); ++v) {
		Xc.row(j) = X.row(*v);
		Xfc.row(j) = Xf.row(*v);
		j++;
	    }

	    ShapeMatching sm(Xc, Xfc, 1., Deformation::LINEAR);
	    MatrixXd D1 = sm.getPureDeformation() - MatrixXd::Identity(3, 3);
	    if (D1.norm() > c_yield){
		S[i] = (MatrixXd::Identity(3, 3) + h * c_creep * D1) * S[i];
		S[i] /= pow(S[i].determinant(), 1/3);
	    }

	    MatrixXd D2 = S[i] - MatrixXd::Identity(3, 3);
	    if (D2.norm() > c_max){
		MatrixXd T = MatrixXd::Identity(3, 3) + c_max * D2.normalized();
		for (list<int>::iterator v = clusters[i].begin(); v != clusters[i].end(); ++v) {
		    Xp.row(*v) = Xp.row(*v) * T.transpose();
		}
	    }
	}
	ShapeMatching sm(Xp, Xf, beta, Deformation::QUADRATIC);
	G = sm.getMatch();
    }
    else {
	ShapeMatching sm(X, Xf, beta, Deformation::QUADRATIC);
	G = sm.getMatch();
    }
}


void Integration::change_destination(MatrixXd& new_Xf) {
	Xf = new_Xf;
	computeDestination();
}


/*
  Integration scheme.
*/

void Integration::performStep(float lambda){
    V *= lambda;
    V += alpha * (G - X) / h;

    if (features[Feature::GRAVITY]) {
	perform_step_gravity();
    }

    if (features[Feature::CLUSTERS]) {
	perform_step_clusters();
    }
    
    X += h * V;
}

void Integration::perform_step_gravity(){
    for (int i = 0; i < V.rows(); i++) {
	V.row(i) -= h*gravity * RowVector3d(0., 1., 0.);
    }
}

void Integration::perform_step_clusters(){
    for (vector<list<int> >::iterator c = clusters.begin(); c != clusters.end(); ++c){
	if (c->empty()) continue;

	MatrixXd Xc(c->size(), 3);
	MatrixXd Xfc(c->size(), 3);
	int i = 0;

	for (list<int>::iterator v = c->begin(); v != c->end(); v++) {
	    Xc.row(i) = X.row(*v);
	    Xfc.row(i) = Xf.row(*v);
	    i++;
	}

	ShapeMatching sm(Xc, Xfc, 0.5, Deformation::LINEAR);
	i = 0;
	for (list<int>::iterator v = c->begin(); v != c->end(); ++v) {
	    V.row(*v) += alpha * (sm.getMatch().row(i) - Xc.row(i)) / h;
	}
    }
}


/*
  Other methods.
*/

bool Integration::check_ground(int axe, double sol,double amortissement) {
	bool contact = false;
	for (int i = 0; i < X.rows(); i++) {
		if (X(i,axe) < sol) {
			std::cout << "Touche le sol" << std::endl;
			X(i,axe) = sol;
			V(i,axe) = amortissement * abs(V(i,axe));
			contact = true;
		}
	}
	return contact;
}
bool Integration::check_height(int axe, double hauteur) {
	bool contact = false;
	for (int i = 0; i < X.rows(); i++) {
		if (X(i, axe) > hauteur) {
			std::cout << "Touche le sol" << std::endl;
			contact = true;
		}
	}
	return contact;
}

bool Integration::check_box(std::map<string,double> box, double amortissement,double epsilon) {
	bool contact = false;
	
	for (int i = 0; i < X.rows(); i++) {
		if (X(i, 0) <= box["x_min"] + epsilon || X(i, 1) <= box["y_min"] + epsilon || X(i, 2) <= box["z_min"] + epsilon || X(i, 0) >= box["x_max"] - epsilon || X(i, 1) >= box["y_max"] - epsilon || X(i, 2) >= box["z_max"] - epsilon) {
			contact = true;
		}
		if (X(i, 0) <= box["x_min"]) {
			X(i, 0) = box["x_min"];
			V(i,0) = - amortissement *V(i,0);
			contact = true;
		}
		if (X(i, 1) <= box["y_min"]) {
			X(i, 1) = box["y_min"];
			V(i, 1) = - amortissement *V(i, 1); 
			contact = true;
		}
		if (X(i, 2) <= box["z_min"]) {
			X(i, 2) = box["z_min"];
			V(i, 2) = - amortissement * V(i, 2);
			contact = true;
		}
		if (X(i, 0) >= box["x_max"]) {
			X(i, 0) = box["x_max"];
			V(i, 0) = - amortissement * V(i, 0);
			contact = true;
		}
		if (X(i, 1) >= box["y_max"]) {
			X(i, 1) = box["y_max"];
			V(i, 1) = - amortissement * V(i, 1);
			contact = true;
		}
		if (X(i, 2) >= box["z_max"]) {
			X(i, 2) = box["z_max"];
			V(i, 2) = - amortissement * V(i, 2);
			contact = true;
		}

	}
	return contact;
}

