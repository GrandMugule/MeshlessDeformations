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
    ShapeMatching(MatrixXd &_X0, MatrixXd &_X, float _beta, Deformation method);
    ShapeMatching(MatrixXd &_X0, MatrixXd &_X);
    ~ShapeMatching();

    MatrixXd getMatch(){ return *G; }

private:
    // Input
    MatrixXd& X0;
    MatrixXd& X;
    float beta;

    // Output
    MatrixXd* G;

    // Deformations
    void linearDeformation();
    void quadraticDeformation();
};
