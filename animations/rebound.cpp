#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOFF.h>
#include <igl/unproject_onto_mesh.h>
#include <iostream>
#include <ostream>
#include <list>
#include <iterator>
#include <map>

#include "shapematching.h"
#include "integration.h"
#include "adjacency.h"
#include "clustering.h"

//NOT FINISHED YET

using namespace Eigen; // to use the classes provided by Eigen library

// Initial mesh
MatrixXd X0;
MatrixXi F;


// Elastic stretching
MatrixXd X;
int axe = 0;
MatrixXd G;
MatrixXd Xf;

// Integration scheme
Integration* I;

//Rebound
bool contact = false;
double amortissement = 0.05;
map<string, double> box;

bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier) {
	std::cout << "pressed Key: " << key << " " << (unsigned int)key << std::endl;


	if ((unsigned int)key == 'D') {
		std::cout << "Animation is running..." << std::endl;
		vector<list<int> > vectorNull;
		I = new Integration(G, X0, 0.1, 0.01, vectorNull);
		viewer.core().is_animating = true;
		return true;
	}

	if ((unsigned int)key == 'C') {
		std::cout << "Change destination" << std::endl;
		for (int i = 0; i < Xf.rows(); i++) {
			RowVector3d translation(0, 0, 1);
			Xf.row(i) += translation;
		}

		I->change_destination(Xf);
		return true;
	}

	if ((unsigned int)key == 'S') {
		std::cout << "Fin integration" << std::endl;
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
		//contact = I->check_box(box, amortissement);
	}
	return false;
}



int main(int argc, char* argv[]) {
	// initialize input mesh

	if (argc < 2) {
		igl::readOFF("../../data/bunny.off", X0, F);
	}
	else {
		igl::readOFF(argv[1], X0, F);
	}

	//  print the number of mesh elements
	std::cout << "Vertices: " << X0.rows() << std::endl;
	std::cout << "Faces:    " << F.rows() << std::endl;
	std::cout << std::endl;

	// elastic stretching matrices
	X = X0;
	G = X0;
	Xf = X0;



	// initialize viewer
	igl::opengl::glfw::Viewer viewer; // create the 3d viewer
	viewer.callback_key_down = &key_down; // for dealing with keyboard events
	viewer.core().is_animating = false;
	viewer.callback_pre_draw = &pre_draw; // to perform animation steps
	viewer.data().set_mesh(X, F); // load a face-based representation of the input 3d shape
	viewer.launch(); // run the editor
}
