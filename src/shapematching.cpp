#include "shapematching.h"

using namespace Eigen;

/*
  Constructors and destructor.
*/

ShapeMatching::~ShapeMatching() {
    delete G;
}

ShapeMatching::ShapeMatching(MatrixXd &_X0, MatrixXd &_X, VectorXd &_W)
    : X0(_X0), X(_X), W(_W)
{
    beta = 0.;

    int n = X0.rows();
    G = new MatrixXd(n, 3);

    linearDeformation();
}

ShapeMatching::ShapeMatching(MatrixXd &_X0, MatrixXd &_X, VectorXd &_W, float _beta, Deformation method)
    : X0(_X0), X(_X), W(_W)
{
    beta = _beta;

    int n = X0.rows();
    G = new MatrixXd(n, 3);
    
    switch (method) {
    case (Deformation::LINEAR) :
	linearDeformation();
	break;
    case (Deformation::QUADRATIC) :
	quadraticDeformation();
	break;
    default :
	std::cout << "This deformation doesn't exist" << std::endl;
    }
}

/*
  Linear and quadratic deformations.
  Rigid deformations are juste a special case of linear deformations, with beta = 0.
*/

RowVector3d weightedMean(const MatrixXd &V, const VectorXd &weights) {
    RowVector3d m(0, 0, 0);

    double s = weights.sum();
    for (int i = 0; i < V.rows(); i++) {
	m += (weights(i) / s) * V.row(i);
    }

    return m;
}

void ShapeMatching::linearDeformation() {
    int n = X0.rows();

    RowVector3d x0cm = weightedMean(X0, W);
    RowVector3d xcm = weightedMean(X, W);

    MatrixXd Apq = MatrixXd::Zero(3, 3);
    MatrixXd Aqq = MatrixXd::Zero(3, 3);
    for (int i = 0; i < n; i++) {
	MatrixXd p(1, 3), q(1, 3);
	p.row(0) = X.row(i) - xcm;
	q.row(0) = X0.row(i) - x0cm;

	Apq += W(i) * p.transpose() * q;
	Aqq += W(i) * q.transpose() * q;
    }

    MatrixXd A = Apq * Aqq.inverse();
    MatrixXd S = (Apq.transpose() * Apq).sqrt();
    MatrixXd R = Apq * S.inverse();
    MatrixXd L = (beta * A) + ((1 - beta) * R);
    for (int i = 0; i < n; i++) {
	G->row(i) = xcm + ((X0.row(i) - x0cm) * L.transpose());
    }
}

void ShapeMatching::quadraticDeformation() {
    int n = X0.rows();

    RowVector3d x0cm = weightedMean(X0, W);
    RowVector3d xcm = weightedMean(X, W);

    MatrixXd p(n, 3);
    MatrixXd r(n, 3);
    MatrixXd q(n, 9);
    for (int i = 0; i < n; i++) {
	p.row(i) = X.row(i) - xcm;
	r.row(i) = X0.row(i) - x0cm;
	q(i, 0) = r(i, 0); q(i, 1) = r(i, 1); q(i, 2) = r(i, 2);
	q(i, 3) = r(i, 0)*r(i, 0); q(i, 4) = r(i, 1)*r(i, 1); q(i, 5) = r(i, 2)*r(i, 2);
	q(i, 6) = r(i, 0)*r(i, 1); q(i, 7) = r(i, 1)*r(i, 2); q(i, 8) = r(i, 2)*r(i, 0);
    }

    // The article doesn't mention it, but maybe we should center q

    MatrixXd Apr = MatrixXd::Zero(3, 3);
    MatrixXd Apq = MatrixXd::Zero(3, 9);
    MatrixXd Aqq = MatrixXd::Zero(9, 9);
    for (int i = 0; i < n; i++) {
	Apr += W(i) * p.row(i).transpose() * r.row(i);
	Apq += W(i) * p.row(i).transpose() * q.row(i);
	Aqq += W(i) * q.row(i).transpose() * q.row(i);
    }
    
    Aqq += 0.001 * MatrixXd::Identity(9, 9);

    MatrixXd A_quad = Apq * Aqq.inverse();

    MatrixXd S_lin = (Apr.transpose() * Apr).sqrt();
    MatrixXd R_lin = Apr * S_lin.inverse();
    MatrixXd R_quad = MatrixXd::Zero(3, 9);
    R_quad.block(0, 0, 3, 3) = R_lin;

    MatrixXd Q = (beta * A_quad) + ((1 - beta) * R_quad);
    for (int i = 0; i < n; i++) {
	G->row(i) = xcm + (q.row(i) * Q.transpose());
    }
}

