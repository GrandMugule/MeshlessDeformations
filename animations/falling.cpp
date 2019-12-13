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
MatrixXd C;

// Elastic falling
MatrixXd X;
int axe = 0;
float sol = 0.;
MatrixXd G;

// the user can specify these parameters
SpectralClustering* SC = nullptr;
float alpha = 0.005;
float beta = 0.5;
float step = 0.1;
double amortissement = 0.05;

// Integration scheme
Integration* I;

//gravity 
bool Gravity = false;
double gravity = 10.;


//plasticity
bool plasticity = false;
int periode = 24; //periode à laquelle on réactualise le shapematching
int compteur = 0;

void init_data(int argc, char* argv[]) {
	assert(argc > 1);

	// input mesh and its adjacency graph
	igl::readOFF(argv[1], X0, F);
	cout << "Vertices : " << X0.rows() << endl;
	cout << "Faces : " << F.rows() << endl << endl;

	// rescale intput mesh
	double scale = (X0.colwise().maxCoeff() - X0.colwise().minCoeff()).norm();
	X0 *= 10 / scale;


	// parse input
	for (int i = 2; i < argc; i += 2) {
		string s(argv[i]); string t(argv[i + 1]);
		if (s.compare("--clusters") == 0) {
			SC = new SpectralClustering(X0, F, stoi(t));
			continue;
		}
		if (s.compare("--alpha") == 0) {
			alpha = stof(t);
			continue;
		}
		if (s.compare("--beta") == 0) {
			beta = stof(t);
			continue;
		}
		if (s.compare("--step") == 0) {
			step = stof(t);
			continue;
		}
		if (s.compare("--gravity") == 0) {
			gravity = stof(t);
			Gravity = true;
			continue;
		}
		if (s.compare("--plasticity") == 0) {
			periode = stof(t);
			plasticity = true;
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
		X = I->currentPosition();
		viewer.core().is_animating = false;
		return true;
	}
	if ((unsigned int)key == 'M') { //M comme Montée
		std::cout << "Montée de l'objet selon l'axe choisi" << std::endl;
		RowVector3d offset(0, 0, 0);
		offset(axe) = 2.;
		for (int i = 0; i < G.rows(); i++) {
			X.row(i) += offset;
		}
		viewer.data().clear();
		viewer.data().set_mesh(X, F);
		viewer.data().set_colors(C);
		return true;
	}
	if ((unsigned int)key == 'C') { //C comme Chute
		//Initialisation de I en fonction des paramètres
		std::cout << "Initialisation chute libre en direction de l'axe choisi" << std::endl;
		sol = X0.colwise().minCoeff()(axe);
		if (SC == nullptr) {
			I = new Integration(X, X0, step, alpha);
		}
		else {
			I = new Integration(X, X0, step, alpha);
			I->addFeature(Feature::CLUSTERS);
			I->setClusters(SC->getClusters());
		}
		if (Gravity) {
			I->addFeature(Feature::GRAVITY);
			I->setGravity(gravity);
		}
		if (plasticity) {
			I->addFeature(Feature::PLASTICITY);
		}
		I->computeDestination(beta);
		viewer.core().is_animating = true;
		return true;
	}
	return false;
}



bool pre_draw(igl::opengl::glfw::Viewer& viewer) {
	if (viewer.core().is_animating) {
		I->performStep();
		//Gestion du contact avec le sol
		bool contact = I->check_ground(axe,sol,amortissement);		
		//Actualisation du matching si on a la plasticite
		if (plasticity && compteur % periode == 0) {
			I->computeDestination();
		}
		compteur += 1;
		viewer.data().clear();
		viewer.data().set_mesh(I->currentPosition(), F);
		viewer.data().set_colors(C);
	}
	return false;
}



int main(int argc, char* argv[]) {
	init_data(argc, argv);

	// initialize viewer
	igl::opengl::glfw::Viewer viewer; // create the 3d viewer
	viewer.callback_key_down = &key_down; // for dealing with keyboard events
	viewer.core().is_animating = false;
	viewer.callback_pre_draw = &pre_draw; // to perform animation steps
	viewer.data().set_mesh(X, F); // load a face-based representation of the input 3d shape
	viewer.data().set_colors(C);
	viewer.launch(); // run the editor
}
