#pragma once

#include <Eigen/Core>
#include <Eigen/SVD>
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/MatrixFunctions>
#include <iostream>

using namespace Eigen;

class ShapeMatching {
public:
    ShapeMatching(const MatrixXd &X0, const MatrixXd &X);
    ~ShapeMatching();

    MatrixXd getMatch();

private:
    MatrixXd* G;
};
