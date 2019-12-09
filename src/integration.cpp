#include "integration.h"
#include "shapematching.h"

#include <iterator>
#include <cassert>

using namespace std;
using namespace Eigen;

/*
  Constructor and setters.
*/

Integration::Integration(MatrixXd& _Xi, MatrixXd& _Xf, float _h, float _alpha, Feature _method){
    Xf = _Xf;
    h = _h;
    alpha = _alpha;
    method = _method;

    int n = _Xi.rows();
    X = _Xi;
    V = MatrixXd::Zero(n, 3);
    ShapeMatching sm(X, Xf, 0.5, Deformation::QUADRATIC);
    G = sm.getMatch();
}

void Integration::setGravity(float g){
    assert(method == Feature::GRAVITY);
    gravity = g;
}

void Integration::setClusters(vector<list<int> >& _clusters){
    assert(method == Feature::CLUSTERS);
    assert(!_clusters.empty());
    clusters = _clusters;
}

/*
  Integration scheme.
*/

void Integration::performStep(float lambda){
    switch(method){
    case (Feature::NONE):
	perform_step(lambda);
	break;
    case (Feature::GRAVITY):
	perform_step_gravity(lambda);
	break;
    case (Feature::CLUSTERS):
	perform_step_clusters(lambda);
	break;
    default:
	cout << "Integration error" << endl;
    }
}

void Integration::perform_step(float lambda){
    V *= lambda;
    V += alpha * (G - X) / h;
    
    X += h * V;
}

void Integration::perform_step_gravity(float lambda){
    V *= lambda;
    V += alpha * (G - X) / h;
    
    for (int i = 0; i < V.rows(); i++) {
	V.row(i) -= gravity * RowVector3d(0., 1., 0.);
    }
    
    X += h * V;
}

void Integration::perform_step_clusters(float lambda){
    V *= lambda;
    V += alpha * (G - X) / h;

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

    X += h * V;
}

/*
  Other methods.
*/

void Integration::check_ground(int axe, double sol) {
	for (int i = 0; i < X.rows(); i++) {
		if (X(i,axe) < sol) {
			std::cout << "Touche le sol" << std::endl;
			X(i,axe) = sol;
			V(i,axe) = 0.05 * abs(V(i,axe));
		}
	}
}

bool Integration::check_box(std::map<string,double> box, double amortissement) {
	bool contact = false;
	for (int i = 0; i < X.rows(); i++) {
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

