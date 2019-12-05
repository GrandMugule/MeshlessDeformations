#include "integration.h"
#include "shapematching.h"

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

    // update positions
    X += h * V;
}
