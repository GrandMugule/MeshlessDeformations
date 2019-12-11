#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOFF.h>
#include <igl/unproject_onto_mesh.h>
#include <iostream>
#include <ostream>
#include <list>
#include <iterator>

#include "shapematching.h"
#include "integration.h"
#include "adjacency.h"
#include "clustering.h"

using namespace Eigen; // to use the classes provided by Eigen library

// Initial mesh
MatrixXd X0;
MatrixXi F;
Adjacency* A;
SpectralClustering* SC;
MatrixXd C;

// Elastic stretching
MatrixXd X;
int axe = 0;
float sol = 0.;
MatrixXd G;

// Integration scheme
Integration* I;

//rebond
double hauteur = 0.4;
bool contact_g = false;
bool previous_contact_g = false;
bool contact_h = false;
bool previous_contact_h = false;
int etape = 0;

bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier) {
	/*
	Resumé des touches utilisees :
	- X,Y,Z : selectionne l'axe pour les déplacements
	- flèche droite, flèche gauche : déplace le point selectionné selon l'axe
	- espace : actualise la forme en effectuant le shapematching
	- S (Stop) : arrête l'intégration. Permet de recommencer en déplacement d'autres points
	- M (Montée) : translate l'ensemble de la forme dans la direction de l'axe (NB attention sur clavier azerty c'est la virgule)
	- C (Chute) : simule une chute libre de l'objet qui revient à sa position initiale, avec un sol et des rebonds
	*/
	std::cout << "pressed Key: " << key << " " << (unsigned int)key << std::endl;
	if (key == 'X') {
		axe = 0;
		std::cout << "modification axe : " << key << " " << (unsigned int)key << std::endl;
		return true;
	}
	if (key == 'Y') {
		axe = 1;
		std::cout << "modification axe : " << key << " " << (unsigned int)key << std::endl;
		return true;
	}
	if (key == 'Z') {
		axe = 2;
		std::cout << "modification axe : " << key << " " << (unsigned int)key << std::endl;
		return true;
	}


	if ((unsigned int)key == 'S') {
		std::cout << "Fin integration" << std::endl;
		G = I->currentPosition();
		viewer.core().is_animating = false;
		return true;
	}
	if ((unsigned int)key == 'M') { //M comme Montée
		std::cout << "Montée de l'objet selon l'axe choisi" << std::endl;
		RowVector3d offset(0, 0, 0);
		offset(axe) = hauteur;
		for (int i = 0; i < G.rows(); i++) {
			G.row(i) += offset;
		}
		viewer.data().clear();
		viewer.data().set_mesh(G, F);
		return true;
	}
	if ((unsigned int)key == 'C') { //C comme Chute
		std::cout << "Initialisation chute libre en direction de l'axe choisi" << std::endl;
		sol = X0.colwise().minCoeff()(axe);
		I = new Integration(G, X0, 0.1, 0.03);
		viewer.core().is_animating = true;
		return true;
	}
	return false;
}

void rebond() {
	contact_g = I->check_ground(axe, sol);
	if (contact_g && !previous_contact_g) {
		MatrixXd Xf = X0;
		hauteur = 0.5 * hauteur;
		RowVector3d offset(0, 0, 0);
		offset(axe) = hauteur;
		for (int i = 0; i < G.rows(); i++) {
			Xf.row(i) += offset;
		}
		I->change_destination(Xf);
	}
	if (contact_g) {
		previous_contact_g = true;
	}
	else {
		previous_contact_g = false;
	}
	contact_h = I->check_height(axe, hauteur);
	if (contact_h && !previous_contact_h) {
		I->change_destination(X0);
	}
	if (contact_h) {
		previous_contact_h = true;
	}
	else {
		previous_contact_h = false;
	}
}

bool pre_draw(igl::opengl::glfw::Viewer& viewer) {
	if (viewer.core().is_animating) {
		I->performStep();
		rebond();
		viewer.data().clear();
		viewer.data().set_mesh(I->currentPosition(), F);
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

	// initialize adjacency graph
	A = new Adjacency(F, X0.rows());

	// perform clustering on initial shape and give a different colour to each cluster
	std::cout << "Performing spectral clustering on input mesh..." << std::endl;
	SC = new SpectralClustering(X0, F, 10);
	std::cout << "Done" << std::endl;
	std::cout << std::endl;

	C = MatrixXd(X0.rows(), 3);
	for (vector<list<int> >::iterator c = SC->getClusters().begin(); c != SC->getClusters().end(); ++c) {
		RowVector3d color = (RowVector3d::Random() + RowVector3d::Constant(1.)) / 2;
		for (list<int>::iterator v = c->begin(); v != c->end(); ++v) {
			C.row(*v) = color;
		}
	}

	// elastic stretching matrices
	X = X0;
	G = X0;

	// initialize viewer
	igl::opengl::glfw::Viewer viewer; // create the 3d viewer
	viewer.callback_key_down = &key_down; // for dealing with keyboard events
	viewer.core().is_animating = false;
	viewer.callback_pre_draw = &pre_draw; // to perform animation steps
	viewer.data().set_mesh(G, F); // load a face-based representation of the input 3d shape
	viewer.launch(); // run the editor
}
