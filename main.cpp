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
  //igl::read_triangle_mesh(argc>1 ? argv[1] : "../shared/data/cube.obj", V, F);
  igl::read_triangle_mesh(argc>1 ? argv[1] : "../shared/data/max-face-low-res.obj", V, F);

  std::vector<NormalSet> normal_sets = compute_normal_sets(F, V);

  // Straighten edges
  Eigen::MatrixXd newV;
  Eigen::MatrixXi newF;
  Eigen::MatrixXd P1, P2;
  straightenEdges(V, F, normal_sets, newV, newF, P1, P2);

  // Create a libigl Viewer object
  igl::viewer::Viewer viewer;
  // Set the vertices and faces for the viewer
  viewer.data.set_mesh(newV, newF);
  // viewer.data.add_points(newV, Eigen::RowVector3d(1,0,0));
  viewer.data.add_edges(P1, P2, Eigen::RowVector3d(0.54, 0.47, 0.39)); // brown for edges

  // white faces for joints
  Eigen::MatrixXd C = Eigen::MatrixXd::Constant(newF.rows(),3,1);
  // // each normal set has a different color
  // for(std::vector<NormalSet>::iterator set = normal_sets.begin(); set != normal_sets.end(); set++){
  //   std::set<int> normal_set = (*set).face_set;
  //   double r = ((double) rand() / (RAND_MAX)); double g = ((double) rand() / (RAND_MAX)); double b = ((double) rand() / (RAND_MAX));
  //   std::set<int>::iterator face;
  //   for (face = normal_set.begin(); face != normal_set.end(); ++face){
  //       int face_idx = *face;
  //       C.row(face_idx) << Eigen::RowVector3d(r,g,b);
  //   }
  // }
  viewer.data.set_colors(C);
  viewer.core.show_lines = false; // don't show wireframe on joints
  // Launch a viewer instance
  viewer.launch();
  return 0;
}
