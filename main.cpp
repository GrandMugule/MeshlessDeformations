#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOFF.h>
#include <igl/unproject_onto_mesh.h>
#include <iostream>
#include <ostream>

#include "shapematching.h"

using namespace Eigen; // to use the classes provided by Eigen library

MatrixXd X0;
MatrixXd X;
MatrixXd G;
MatrixXi F;
int axe = 0;
int currentVertex;

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
		MatrixXd Oldpoint = G.row(currentVertex);
		MatrixXd Newpoint = X.row(currentVertex);
		if ((unsigned int)key == 6) { //fleche de droite
			X(currentVertex, axe) += 0.1;
			std::cout << "offset : " << 0.1 << " " << std::endl;
		}
		else {
			X(currentVertex, axe) -= 0.1;
			std::cout << "offset : " << -0.1 << " " << std::endl;
		}
		viewer.data().clear();
		viewer.data().add_points(Newpoint, RowVector3d(0, 1, 0));
		viewer.data().add_points(Oldpoint, RowVector3d(1, 0, 0));
		viewer.data().add_edges(Oldpoint, Newpoint, Eigen::RowVector3d(0, 0, 1));
		viewer.data().set_mesh(G, F);
		return true;
	}
	if ((unsigned int)key == 32) { //touche espace
		//update la forme de G
		//update(X0,X,G)
		G = ShapeMatching(X0, X).getMatch();
		viewer.data().clear();
		viewer.data().set_mesh(G, F);
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
		//update la forme de G si on touche un nouveau point
		// update(X0,X,G);
        int bcid;
        bc.maxCoeff(&bcid);
        vid = F(fid, bcid);
        MatrixXd P(1, 3);
        P.row(0) = G.row(vid);
		currentVertex = vid;
        viewer.data().add_points(P, RowVector3d(1, 0, 0));
        return true;
    }
    return false;
}



int main(int argc, char *argv[]) {
    if (argc < 2) {
        igl::readOFF("../data/octagon.off", X0, F); // default input mesh
    }
    else {
        igl::readOFF(argv[1], X0, F); // input mesh given in command line
    }
	X = X0;
	G = X0;
    //  print the number of mesh elements
    std::cout << "Vertices: " << X.rows() << std::endl;
    std::cout << "Faces:    " << F.rows() << std::endl;

    igl::opengl::glfw::Viewer viewer; // create the 3d viewer
    viewer.callback_key_down = &key_down; // for dealing with keyboard events
    viewer.callback_mouse_down = &mouse_down;
    viewer.data().set_mesh(G, F); // load a face-based representation of the input 3d shape
    viewer.launch(); // run the editor
}
