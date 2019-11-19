#include "shapematching.h"

using namespace Eigen;

ShapeMatching::ShapeMatching(const MatrixXd &X0, const MatrixXd &X) {
    // for now only implementing rigid deformation from X0 to X
    n = X0.rows();
    G = new MatrixXd(n, 3);

    RowVector3d x0cm = X0.rowwise().mean();
    RowVector3d xcm = X.rowwise().mean();

    MatrixXd Apq(3, 3) = MatrixXd::Zero(3, 3);
    for (int i = 0; i < n; i++) {
        MatrixXd p(1, 3), q(1, 3);
        p.row(0) = X0.row(i) - x0cm;
        q.row(0) = X.row(i) - xcm;

        Apq += p.transpose() * q;
    }

    MatrixXd S = (Apq.transpose() * Apq).sqrt();
    MatrixXd R = Apq * S.inverse();

    for (int i = 0; i < n; i++) {
        G.row(i) = xcm + R*(X0.row(i) - x0cm);
    }
}

ShapeMatching::~ShapeMatching() {
    delete G;
}

MatrixXd ShapeMatching::getMatch() {
    return *G;
}
