#pragma once

using namespace Eigen;

class ShapeMatching {
public:
    ShapeMatching(const MatrixXd &X0, const MatrixXd &X);
    ~ShapeMatching();

    MatrixXd getMatch();

private:
    MatrixXd* G;
};
