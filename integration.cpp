#include "integration.h"

using namespace Eigen;

/*
  Constructor and destructor.
*/

Integration::~Integration() {
	delete X;
	delete V;
}

Integration::Integration(MatrixXd& _Xi, MatrixXd& _Xf, float _h, float _alpha)
	: Xf(_Xf)
{
	h = _h;
	alpha = _alpha;

	int n = _Xi.rows();
	X = new MatrixXd();
	*X = _Xi;
	V = new MatrixXd();
	*V = MatrixXd::Zero(n, 3);

	gravity = false;

}

Integration::Integration(MatrixXd& _Xi, MatrixXd& _Xf, float _h, float _alpha, bool _gravity, int _axis_gravity)
	: Integration(_Xi,  _Xf,  _h,  _alpha)
{
	gravity = _gravity;
	axis_gravity = _axis_gravity;
	sol = _Xf(0, _axis_gravity);
	for (int i = 0; i < _Xf.rows(); i++) {
		float hauteur = _Xf(i, _axis_gravity);
		if (hauteur < sol) {
			sol = hauteur;
		}
	}
}

/*
  Integration scheme.
*/

void Integration::performStep(float lambda) {
	*V *= lambda;
	*V += alpha * (Xf - *X) / h;
    *X += h * *V;
	if (gravity) {
		for (int i = 0; i < X->rows(); i++) {
			if (X->row(i)(axis_gravity) < sol) {
				X->row(i)(axis_gravity) = sol;
				V->row(i)(axis_gravity) = 0.05 * abs(V->row(i)(axis_gravity));
			}
		}
	}
}
