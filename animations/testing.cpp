#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOFF.h>

#include <iterator>

#include "clustering.h"

using namespace std;
using namespace Eigen;

MatrixXd V;
MatrixXi F;
MatrixXd C;

int main(int argc, char *argv[]) {
    igl::readOFF("../data/bunny.off", V, F);

    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(V, F);

    int nClusters = 30;
    SpectralClustering sc(V, F, nClusters);
    vector<list<int> > clusters = sc.getClusters();
    C = MatrixXd(V.rows(), 3);
    for (int i = 0; i < nClusters; i++) {
	RowVector3d color = (RowVector3d::Random() + RowVector3d::Constant(1.)) / 2;
	for (list<int>::iterator it = clusters[i].begin(); it != clusters[i].end(); ++it) {
	    C.row(*it) = color;
	}
    }
    viewer.data().set_colors(C);
    
    viewer.launch();
}
