#pragma once

#include <Eigen/Core>
#include <Eigen/SVD>
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/MatrixFunctions>

using namespace Eigen;

class Integration {
 public:
    Integration(MatrixXd &_Xi, MatrixXd &_Xf, float _h, float _alpha);
    ~Integration();

    void performStep();
    MatrixXd currentPosition(){ return *X; }

 private:
    MatrixXd& Xf;
    float h;
    float alpha;

    MatrixXd* X;
    MatrixXd* V;
};
