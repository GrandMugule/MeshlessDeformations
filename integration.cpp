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
    *V += alpha * (Xf - *X) / h;
    *X += h * *V;
}
