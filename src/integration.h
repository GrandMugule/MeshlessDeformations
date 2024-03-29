#pragma once

#include "shapematching.h"

#include <Eigen/Core>

#include <vector>
#include <list>
#include <map>

using namespace std;
using namespace Eigen;

enum class Feature {
    GRAVITY,
    CLUSTERS,
    PLASTICITY
};

class Integration {
 public:
    Integration(MatrixXd& _Xi, MatrixXd& _Xf, float _h, float _alpha);
    ~Integration(){}

    void addFeature(Feature f);
    void setGravity(int _axe, float g);
    void setClusters(vector<list<int> >& _clusters);
    void setPlasticityCoeffs(float _c_yield, float _c_creep, float _c_max);
    void computeDestination(float beta = 0.5, Deformation deformation = Deformation::QUADRATIC);
    void change_destination(MatrixXd& new_Xf);
    void change_matching(MatrixXd _G);

    void performStep(float lambda = 0.9);
    MatrixXd& currentPosition(){ return X; }
    
    bool check_ground(int axe, double sol, double amortissement);
    bool check_height(int axe, double hauteur);
    bool check_box(map<string, double> box, double amortissement, double epsilon);

 private:
    MatrixXd Xf;
    float h;
    float alpha;
    map<Feature, bool> features;
    
    vector<MatrixXd> S;
    float c_yield = 0.5;
    float c_creep = 0.5;
    float c_max = 0.5;

    MatrixXd X;
    MatrixXd V;
    MatrixXd G;

    float gravity;
	int axe;
    vector<list<int> > clusters;

    void perform_step_gravity();
    void perform_step_clusters();
};
