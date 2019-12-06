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
int currentVertex;
list<int> currentNeighborhood;
MatrixXd G;

// Integration scheme
Integration* I;

bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifier) {
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
	if ((unsigned int)key == 6 || (unsigned int)key == 7) { //6 : touche fleche de droite, 7 : touche fleche de gauche
		RowVector3d oldPoint = G.row(currentVertex);
		RowVector3d newPoint = X.row(currentVertex);
		if ((unsigned int)key == 6) { //fleche de droite
		        newPoint(axe) += 0.1;
			std::cout << "offset : " << 0.1 << " " << std::endl;
		}
		else {
		        newPoint(axe) -= 0.1;
			std::cout << "offset : " << -0.1 << " " << std::endl;
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
	        G = ShapeMatching(X0, X, 0.5, Deformation::QUADRATIC).getMatch();
		viewer.data().clear();
		viewer.data().set_mesh(G, F);
		viewer.data().set_colors(C);
		return true;
	}
	if ((unsigned int)key == 'D') {
		std::cout << "Animation is running..." << std::endl;
		I = new Integration(G, X0, 0.1, 0.1, SC->getClusters());
		viewer.core().is_animating = true;
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



bool pre_draw(igl::opengl::glfw::Viewer& viewer) {
    if (viewer.core().is_animating) {
	I->performStep();
	viewer.data().clear();
	viewer.data().set_mesh(I->currentPosition(), F);
	viewer.data().set_colors(C);
    }
    return false;
}



int main(int argc, char *argv[]) {
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

    // initialize adjacency graph
    std::cout << "Computing adjacency graph of input mesh..." << std::endl;
    A = new Adjacency(F, X0.rows());
    std::cout << "Done" << std::endl;
    std::cout << std::endl;

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
    viewer.callback_mouse_down = &mouse_down;
    viewer.core().is_animating = false;
    viewer.callback_pre_draw = &pre_draw; // to perform animation steps
    viewer.data().set_mesh(G, F); // load a face-based representation of the input 3d shape
    viewer.data().set_colors(C);
    viewer.launch(); // run the editor
}
