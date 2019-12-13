#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOFF.h>
#include <igl/unproject_onto_mesh.h>
#include <iostream>
#include <list>
#include <iterator>

#include "shapematching.h"
#include "integration.h"
#include "adjacency.h"
#include "clustering.h"

#include <string>
#include <cassert>

using namespace std;
using namespace Eigen;


/*
  Global variables.
*/

// Initial mesh
MatrixXd X0;
MatrixXi F;
Adjacency* A = nullptr;
MatrixXd C;

// the user can specify these parameters
SpectralClustering* SC = nullptr;
float alpha = 0.1;
float beta = 0.5;
float step = 0.1;
Deformation deformation = Deformation::QUADRATIC;

// Elastic stretching
MatrixXd X;
int axe = 0;
int currentVertex;
list<int> currentNeighborhood;
MatrixXd G;

// Integration scheme
Integration* I = nullptr;


/*
  Parse input command line to initialize global variables.
*/

void init_data(int argc, char *argv[]){
    assert(argc > 1);
    
    // input mesh and its adjacency graph
    igl::readOFF(argv[1], X0, F);
    cout << "Vertices : " << X0.rows() << endl;
    cout << "Faces : " << F.rows() << endl << endl;
    A = new Adjacency(F, X0.rows());

    // rescale intput mesh
    double scale = (X0.colwise().maxCoeff() - X0.colwise().minCoeff()).norm();
    X0 *= 10 / scale;


    // parse input
    for (int i = 2; i < argc; i += 2) {
	string s(argv[i]); string t(argv[i+1]);
	if (s.compare("--clusters") == 0){
	    SC = new SpectralClustering(X0, F, stoi(t));
	    continue;
	}
	if (s.compare("--alpha") == 0){
	    alpha = stof(t);
	    continue;
	}
	if (s.compare("--beta") == 0){
	    beta = stof(t);
	    continue;
	}
	if (s.compare("--step") == 0){
	    step = stof(t);
	    continue;
	}
	if (s.compare("--deformation") == 0){
	    if (t.compare("r") == 0) deformation = Deformation::RIGID;
	    else if (t.compare("l") == 0) deformation = Deformation::LINEAR;
	    else if (t.compare("q") == 0) deformation = Deformation::QUADRATIC;
	    continue;
	}
    }

    // set color
    C = MatrixXd(X0.rows(), 3);
    if (SC == nullptr) {
	for (int i = 0; i < C.rows(); i++) {
	    C.row(i) << 1.0, 1.0, 0.0;
	}
    }
    else {
	for (vector<list<int> >::iterator c = SC->getClusters().begin(); c != SC->getClusters().end(); ++c) {
	    if (c->empty()) continue;
	    RowVector3d color = (RowVector3d::Random() + RowVector3d::Constant(1.)) / 2;
	    for (list<int>::iterator v = c->begin(); v != c->end(); ++v) {
		C.row(*v) = color;
	    }
	}
    }

    X = X0;
    G = X0;
}


/*
  Callbacks : 
  - mouse_down is used to select a vertex
  - kew_down is used to move the selected vertex and start/stop the animation
  - pre_draw is called when animation is launched to perform steps in the integration scheme.
*/

bool mouse_down(igl::opengl::glfw::Viewer& viewer, int button, int modifier) {
    int fid, vid;
    Vector3f bc;
    double x = viewer.current_mouse_x;
    double y = viewer.core().viewport(3) - viewer.current_mouse_y;
    Vector2f mouse_position(x, y);
    if (igl::unproject_onto_mesh(mouse_position, viewer.core().view, viewer.core().proj, viewer.core().viewport, G, F, fid, bc)) {
	// find nearest vertex
        int bcid;
        bc.maxCoeff(&bcid);
        vid = F(fid, bcid);
        MatrixXd P(1, 3);
        P.row(0) = G.row(vid);
	
	// update current vertex and neighborhood
	currentVertex = vid;
	currentNeighborhood = A->getNeighborhood(currentVertex);

	// add a red dot on the viewer
        viewer.data().add_points(P, RowVector3d(1, 0, 0));
        return true;
    }
    return false;
}

bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifier) {
    std::cout << "pressed Key: " << key << " " << (unsigned int)key << std::endl;
    if (key == '1') {
	axe = 0;
	cout << "modification axe : " << key << " " << (unsigned int)key << endl;
	return true;
    }
    if (key == '2') {
	axe = 1;
	cout << "modification axe : " << key << " " << (unsigned int)key << endl;
        return true;
    }
    if (key == '3') {
	axe = 2;
	cout << "modification axe : " << key << " " << (unsigned int)key << endl;
	return true;
    }
    if ((unsigned int)key == 6 || (unsigned int)key == 7) { //6 : touche fleche de droite, 7 : touche fleche de gauche
	RowVector3d oldPoint = G.row(currentVertex);
	RowVector3d newPoint = X.row(currentVertex);
	if ((unsigned int)key == 6) { //fleche de droite
	    newPoint(axe) += 0.1;
	}
	else {
	    newPoint(axe) -= 0.1;
	}
	X.row(currentVertex) = newPoint;
	viewer.data().clear();
	viewer.data().add_points(newPoint, RowVector3d(0, 1, 0));
	viewer.data().add_points(oldPoint, RowVector3d(1, 0, 0));
	viewer.data().add_edges(oldPoint, newPoint, Eigen::RowVector3d(0, 0, 1));
	viewer.data().set_mesh(G, F);
	viewer.data().set_colors(C);
	return true;
    }
    if ((unsigned int)key == 32) { //touche espace
        // update all vertices in neighborhood
        RowVector3d delta = X.row(currentVertex) - G.row(currentVertex);
        for (list<int>::iterator it = currentNeighborhood.begin(); it != currentNeighborhood.end(); ++it) {
	    X.row(*it) += delta;
	}
	//update G
        ShapeMatching sm = (deformation == Deformation::RIGID) ? ShapeMatching(X0, X) : ShapeMatching(X0, X, beta, deformation);
	G = sm.getMatch();
	viewer.data().clear();
	viewer.data().set_mesh(G, F);
	viewer.data().set_colors(C);
	return true;
    }
    if ((unsigned int)key == 'D') {
    if (SC == nullptr) {
	I = new Integration(G, X0, step, alpha);
    }
    else {
	I = new Integration(G, X0, step, alpha);
	I->addFeature(Feature::CLUSTERS);
	I->setClusters(SC->getClusters());
    }
    I->change_matching(X0);
    cout << "Animation is running..." << endl;
    viewer.core().is_animating = true;
    return true;
    }
    if ((unsigned int)key == 'S') {
	std::cout << "Animation stopped" << std::endl;
	G = I->currentPosition();
	viewer.core().is_animating = false;
	return true;
    }
    return false;
}

bool pre_draw(igl::opengl::glfw::Viewer& viewer) {
    if (viewer.core().is_animating) {
	I->performStep();
	viewer.data().clear();
	viewer.data().set_mesh(I->currentPosition(), F);
	viewer.data().set_colors(C);
    }
    return false;
}


/*
  Main function.
*/

int main(int argc, char *argv[]) {
    init_data(argc, argv);

    // initialize viewer
    igl::opengl::glfw::Viewer viewer; // create the 3d viewer
    viewer.callback_key_down = &key_down; // for dealing with keyboard events
    viewer.callback_mouse_down = &mouse_down;
    viewer.core().is_animating = false;
    viewer.callback_pre_draw = &pre_draw; // to perform animation steps
    viewer.data().set_mesh(G, F); // load a face-based representation of the input 3d shape
    viewer.data().set_colors(C);
    viewer.launch(); // run the editor
}
