#pragma once

#include <Eigen/Core>
#include <Eigen/SVD>
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/MatrixFunctions>
#include <iostream>

using namespace Eigen;

enum class Deformation {
    RIGID,
    LINEAR,
    QUADRATIC
};

class ShapeMatching {
public:
    // Constructor without beta performs rigid matching
    ShapeMatching(MatrixXd &_X0, MatrixXd &_X, float _beta, Deformation _method);
    ShapeMatching(MatrixXd &_X0, MatrixXd &_X);
    ~ShapeMatching();

    MatrixXd getMatch();
    MatrixXd getPureDeformation();

private:
    // Input
    MatrixXd& X0;
    MatrixXd& X;
    float beta;
    Deformation method;

    // Output
    MatrixXd* G;

    // Pure deformation (not defined in the quadratic case)
    MatrixXd* P;

    // Deformation methods
    void linearDeformation();
    void quadraticDeformation();
};
