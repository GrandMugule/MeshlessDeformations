#pragma once

#include <Eigen/Core>
#include <Eigen/SVD>
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/MatrixFunctions>
#include <iostream>

using namespace Eigen;

enum class Deformation {
    LINEAR,
    QUADRATIC
};

class ShapeMatching {
public:
    // Constructor without beta performs rigid matching
    ShapeMatching(MatrixXd &_X0, MatrixXd &_X, VectorXd &_W, float _beta, Deformation method);
    ShapeMatching(MatrixXd &_X0, MatrixXd &_X, VectorXd &_W);
    ~ShapeMatching();

    MatrixXd getMatch(){ return *G; }

private:
    // Input
    MatrixXd& X0;
    MatrixXd& X;
    VectorXd& W;
    float beta;

    // Output
    MatrixXd* G;

    // Deformations
    void linearDeformation();
    void quadraticDeformation();
};