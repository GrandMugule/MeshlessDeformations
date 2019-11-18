#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOFF.h>
#include <igl/unproject_onto_mesh.h>
#include <iostream>
#include <ostream>

using namespace Eigen; // to use the classes provided by Eigen library

MatrixXd V;
MatrixXi F;
MatrixXd C;

bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifier) {
    std::cout << "pressed Key: " << key << " " << (unsigned int)key << std::endl;

    return false;
}

bool mouse_down(igl::opengl::glfw::Viewer& viewer, int button, int modifier) {
    int fid, vid;
    Vector3f bc;
    double x = viewer.current_mouse_x;
    double y = viewer.core().viewport(3) - viewer.current_mouse_y;
    Vector2f mouse_position(x, y);
    if (igl::unproject_onto_mesh(mouse_position, viewer.core().view, viewer.core().proj, viewer.core().viewport, V, F, fid, bc)) {
        int bcid;
        bc.maxCoeff(&bcid);
        vid = F(fid, bcid);
        MatrixXd P(1, 3);
        P.row(0) = V.row(vid);
        viewer.data().add_points(P, RowVector3d(1, 0, 0));
        return true;
    }
    return false;
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        igl::readOFF("../data/octagon.off", V, F); // default input mesh
    }
    else {
        igl::readOFF(argv[1], V, F); // input mesh given in command line
    }

    //  print the number of mesh elements
    std::cout << "Vertices: " << V.rows() << std::endl;
    std::cout << "Faces:    " << F.rows() << std::endl;

    C = MatrixXd::Constant(F.rows(), 3, 1);

    igl::opengl::glfw::Viewer viewer; // create the 3d viewer
    viewer.callback_key_down = &key_down; // for dealing with keyboard events
    viewer.callback_mouse_down = &mouse_down;
    viewer.data().set_mesh(V, F); // load a face-based representation of the input 3d shape
    viewer.data().set_colors(C);
    viewer.launch(); // run the editor
}
