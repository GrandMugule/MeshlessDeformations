#pragma once

#include <Eigen/Core>
#include <Eigen/SVD>
#include <Eigen/Eigenvalues>
#include <iostream>
#include <unsupported/Eigen/MatrixFunctions>

using namespace Eigen;

class Integration {
 public:
	Integration(MatrixXd& _Xi, MatrixXd& _Xf, float _h, float _alpha);
    Integration(MatrixXd &_Xi, MatrixXd &_Xf, float _h, float _alpha, bool _gravity, int _axis_gravity);
    ~Integration();

    void performStep(float lambda = 0.9);
    MatrixXd currentPosition(){ return *X; }

 private:
    MatrixXd& Xf;
    float h;
    float alpha;


    MatrixXd* X;
    MatrixXd* V;

	//in case of gravity
	bool gravity;
	int axis_gravity;
	float sol;
};
