#include "integration.h"

using namespace Eigen;

/*
  Constructor and destructor.
*/

Integration::~Integration() {
    delete X;
    delete V;
}

Integration::Integration(MatrixXd &_Xi, MatrixXd &_Xf, float _h, float _alpha)
    : Xf(_Xf)
{
    h = _h;
    alpha = _alpha;

    int n = _Xi.rows();
    X = new MatrixXd();
    *X = _Xi;
    V = new MatrixXd();
    *V = MatrixXd::Zero(n, 3);
}

/*
  Integration scheme.
*/

void Integration::performStep() {
    int n = X->rows();
    for (int i = 0; i < n; i++) {
	V->row(i) += alpha * (Xf.row(i) - X->row(i)) / h;
	X->row(i) += h * V->row(i);
    }
}
