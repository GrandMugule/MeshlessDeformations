#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOFF.h>
#include <iostream>
#include <ostream>

using namespace Eigen; // to use the classes provided by Eigen library

MatrixXd V; // matrix storing the vertex coordinates of the input mesh (n rows, 3 columns)
MatrixXi F; // incidence relations between faces and edges (f rows, 3 columns)

bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifier) {
    std::cout << "pressed Key: " << key << " " << (unsigned int)key << std::endl;

    return false;
}

int main(int argc, char *argv[]) {
    igl::readOFF("../data/octagon.off", V, F); // load input mesh in off format

    //  print the number of mesh elements
    std::cout << "Vertices: " << V.rows() << std::endl;
    std::cout << "Faces:    " << F.rows() << std::endl;

    igl::opengl::glfw::Viewer viewer; // create the 3d viewer
    viewer.callback_key_down = &key_down; // for dealing with keyboard events
    viewer.data().set_mesh(V, F); // load a face-based representation of the input 3d shape
    viewer.launch(); // run the editor
}
