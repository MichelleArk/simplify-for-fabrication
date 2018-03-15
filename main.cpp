#include "edges.h"
#include "euler_characteristic.h"
#include <igl/read_triangle_mesh.h>
#include <igl/viewer/Viewer.h>
#include "compute_normal_sets.h"
#include <normalSet.h>
#include <vector>

int main(int argc, char *argv[])
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  // Load in a mesh
  igl::read_triangle_mesh(argc>1 ? argv[1] : "../shared/data/bunny.off", V, F);


  std::vector<NormalSet> normal_sets = compute_normal_sets(F, V);
  // Create a libigl Viewer object
  igl::viewer::Viewer viewer;
  // Set the vertices and faces for the viewer
  viewer.data.set_mesh(V, F);

  Eigen::MatrixXd C = Eigen::MatrixXd::Constant(F.rows(),3,1);
  // each normal set has a different color
  for(std::vector<NormalSet>::iterator set = normal_sets.begin(); set != normal_sets.end(); set++){
    std::set<int> normal_set = (*set).face_set;
    double r = ((double) rand() / (RAND_MAX)); double g = ((double) rand() / (RAND_MAX)); double b = ((double) rand() / (RAND_MAX));
    std::set<int>::iterator face;
    for (face = normal_set.begin(); face != normal_set.end(); ++face){
        int face_idx = *face;

        C.row(face_idx) << Eigen::RowVector3d(r,g,b);
    }
  }

  viewer.data.set_colors(C);
  // Launch a viewer instance
  viewer.launch();
  return 0;
}
