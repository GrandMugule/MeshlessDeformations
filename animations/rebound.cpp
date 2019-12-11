#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOFF.h>
#include <igl/unproject_onto_mesh.h>
#include <iostream>
#include <ostream>
#include <list>
#include <iterator>
#include <map>
#include <stdlib.h>

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
bool previous_contact = false;
double amortissement = 0.9;
double epsilon = 0.05;
map<string, double> box;
RowVector3d PointXf;
int face_destination = 0;
int etape = 0;

void random_destination() {
	
	int random = rand() % 2 +1;
	face_destination = (face_destination + random) % 3;
	PointXf(face_destination)= 1.5* (rand() % 2) * (box["x_max"] -box["x_min"]) + box["x_min"];
	PointXf((face_destination+1)%3) = (box["x_max"] + box["x_min"]) /2. ;
	PointXf((face_destination + 2) % 3) = (box["x_max"] + box["x_min"]) /2.;
	/*
	PointXf((face_destination + 1) % 3) = (box["x_max"] + box["x_min"]) / 2.;
	PointXf((face_destination + 2) % 3) = (box["x_max"] + box["x_min"]) / 2.;
	if (etape % 2 == 0) {
		PointXf(face_destination) = 1.5 * box["x_max"];
	}
	else {
		PointXf(face_destination) = 1.5 * box["x_min"];
	}
	etape += 1;
	*/
	RowVector3d xfcm = Xf.colwise().mean();
	MatrixXd Rotation(3, 3);
	Rotation = Matrix3d::Zero();
	Rotation(face_destination, face_destination) = 1.;
	Rotation((face_destination+2)%3, (face_destination +1) % 3) = 1.;
	Rotation((face_destination+1) % 3, (face_destination +2) % 3) = - 1.;
	for (int i = 0; i < Xf.rows(); i++) {
		Xf.row(i) -= xfcm;
		Xf.row(i) = Xf.row(i) * Rotation.transpose();
		Xf.row(i) += PointXf;
	}

}

void view_box(igl::opengl::glfw::Viewer & viewer) {
	RowVector3d Point1(box["x_max"], box["y_max"], box["z_max"]);
	RowVector3d Point2(box["x_max"], box["y_max"], box["z_min"]);
	RowVector3d Point3(box["x_max"], box["y_min"], box["z_max"]);
	RowVector3d Point4(box["x_max"], box["y_min"], box["z_min"]);
	RowVector3d Point5(box["x_min"], box["y_max"], box["z_max"]);
	RowVector3d Point6(box["x_min"], box["y_max"], box["z_min"]);
	RowVector3d Point7(box["x_min"], box["y_min"], box["z_max"]);
	RowVector3d Point8(box["x_min"], box["y_min"], box["z_min"]);

	viewer.data().add_points(Point1, RowVector3d(0, 1, 0));
	viewer.data().add_points(Point2, RowVector3d(0, 1, 0));
	viewer.data().add_points(Point3, RowVector3d(0, 1, 0));
	viewer.data().add_points(Point4, RowVector3d(0, 1, 0));
	viewer.data().add_points(Point5, RowVector3d(0, 1, 0));
	viewer.data().add_points(Point6, RowVector3d(0, 1, 0));
	viewer.data().add_points(Point7, RowVector3d(0, 1, 0));
	viewer.data().add_points(Point8, RowVector3d(0, 1, 0));

	viewer.data().add_edges(Point1, Point2, Eigen::RowVector3d(0, 0, 1));
	viewer.data().add_edges(Point1, Point3, Eigen::RowVector3d(0, 0, 1));
	viewer.data().add_edges(Point1, Point5, Eigen::RowVector3d(0, 0, 1));
	viewer.data().add_edges(Point2, Point6, Eigen::RowVector3d(0, 0, 1));
	viewer.data().add_edges(Point2, Point4, Eigen::RowVector3d(0, 0, 1));
	viewer.data().add_edges(Point3, Point4, Eigen::RowVector3d(0, 0, 1));
	viewer.data().add_edges(Point3, Point7, Eigen::RowVector3d(0, 0, 1));
	viewer.data().add_edges(Point4, Point8, Eigen::RowVector3d(0, 0, 1));
	viewer.data().add_edges(Point5, Point6, Eigen::RowVector3d(0, 0, 1));
	viewer.data().add_edges(Point5, Point7, Eigen::RowVector3d(0, 0, 1));
	viewer.data().add_edges(Point6, Point8, Eigen::RowVector3d(0, 0, 1));
	viewer.data().add_edges(Point7, Point8, Eigen::RowVector3d(0, 0, 1));
}

bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier) {
	std::cout << "pressed Key: " << key << " " << (unsigned int)key << std::endl;


	if ((unsigned int)key == 'D') {
		std::cout << "Animation is running..." << std::endl;
		I = new Integration(G, Xf, 0.001, 0.001);
		viewer.core().is_animating = true;
		return true;
	}

	if ((unsigned int)key == 'C') {
		std::cout << "Change destination" << std::endl;
		random_destination();
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
		I->performStep(0.95);
		viewer.data().clear();
		view_box(viewer);
		viewer.data().set_mesh(I->currentPosition(), F);
		contact = I->check_box(box, amortissement,epsilon);
		if (contact && !previous_contact) {
			random_destination();
			I->change_destination(Xf);
		}
		previous_contact = contact;
		/*
		if (contact) {
			random_destination();
			I->change_destination(Xf);
		}
		*/
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

	//def of the box
	box["x_min"] = -1.;
	box["x_max"] = 1.;
	box["y_min"] = -1.;
	box["y_max"] = 1.;
	box["z_min"] = -1.;
	box["z_max"] = 1.;



	// initialize viewer
	igl::opengl::glfw::Viewer viewer; // create the 3d viewer
	view_box(viewer);
	viewer.callback_key_down = &key_down; // for dealing with keyboard events
	viewer.core().is_animating = false;
	viewer.callback_pre_draw = &pre_draw; // to perform animation steps
	viewer.data().set_mesh(X, F); // load a face-based representation of the input 3d shape
	viewer.launch(); // run the editor
}