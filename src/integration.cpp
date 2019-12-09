#include "integration.h"
#include "shapematching.h"
#include <map>

#include <iterator>

using namespace Eigen;

/*
  Constructor and destructor.
*/

Integration::Integration(MatrixXd _Xi, MatrixXd &_Xf, float _h, float _alpha, vector<list<int> >&_clusters)
    : Xf(_Xf), clusters(_clusters)
{
    h = _h;
    alpha = _alpha;

    int n = _Xi.rows();
    X = _Xi;
    V = MatrixXd::Zero(n, 3);
}

/*
  Integration scheme.
*/

void Integration::performStep(float lambda) {
    // loss of energy at each step so that the movement eventually stabilizes
    V *= lambda;

    // match current shape with objective shape
    ShapeMatching sm(X, Xf, 0.5, Deformation::QUADRATIC);
    V += alpha * (sm.getMatch() - X) / h;

	
	if (!clusters.empty()) {
		// match clusters with objective shape
		for (vector<list<int> >::iterator c = clusters.begin(); c != clusters.end(); ++c) {
			if (c->empty()) continue;

			MatrixXd Xc(c->size(), 3);
			MatrixXd Xfc(c->size(), 3);
			int i = 0;

			for (list<int>::iterator v = c->begin(); v != c->end(); ++v) {
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
	
    // update positions
    X += h * V;
}

void Integration::performStep_gravity(float lambda) {
	// loss of energy at each step so that the movement eventually stabilizes
	V *= lambda;

	// match current shape with objective shape
	ShapeMatching sm(X, Xf, 0.5, Deformation::QUADRATIC);
	V += alpha * (sm.getMatch() - X) / h;

	//add gravity
	for (int i = 0; i < V.rows(); i++) {
		V.row(i) += RowVector3d(0, -0.01, 0);
	}
	// update positions
	X += h * V;
}

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

