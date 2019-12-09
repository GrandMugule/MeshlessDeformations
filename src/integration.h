#pragma once

#include <Eigen/Core>

#include <vector>
#include <list>
#include <map>

using namespace std;
using namespace Eigen;

enum class Feature {
    NONE,
    GRAVITY,
    CLUSTERS
};

class Integration {
 public:
    Integration(MatrixXd& _Xi, MatrixXd& _Xf, float _h, float _alpha, Feature _method = Feature::NONE);
    ~Integration(){}

    void setGravity(float g);
    void setClusters(vector<list<int> >& _clusters);

    void performStep(float lambda = 0.9);
    MatrixXd& currentPosition(){ return X; }
    
    void check_ground(int axe, double sol);
    bool check_box(map<string, double> box, double amortissement);
    void change_destination(MatrixXd& new_Xf) { Xf = new_Xf; }

 private:
    MatrixXd Xf;
    float h;
    float alpha;
    Feature method;

    MatrixXd X;
    MatrixXd V;
    MatrixXd G;

    float gravity;
    vector<list<int> > clusters;

    void perform_step(float lambda);
    void perform_step_gravity(float lambda);
    void perform_step_clusters(float lambda);
};
