#pragma once

#include <Eigen/Core>

#include <vector>
#include <list>
#include <map>

using namespace std;
using namespace Eigen;

class Integration {
 public:
    Integration(MatrixXd _Xi, MatrixXd &_Xf, float _h, float _alpha, vector<list<int> > &_clusters);
    ~Integration(){}

    void performStep(float lambda = 0.9);
	void check_ground(int axe, double sol);
	bool check_box(std::map<string, double> box, double amortissement);
	void change_destination(MatrixXd& new_Xf) { Xf = new_Xf; }
    MatrixXd& currentPosition(){ return X; }

 private:
    MatrixXd& Xf;
    float h;
    float alpha;
    vector<list<int> >& clusters;

    MatrixXd X;
    MatrixXd V;
};
