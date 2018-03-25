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
  Eigen::VectorXd Cost;
  straightenEdges(V, F, normal_sets, newV, newF, P1, P2, Cost);

  // Create a libigl Viewer object
  igl::viewer::Viewer viewer;
  // Set the vertices and faces for the viewer
  //viewer.data.set_mesh(V, F);
  viewer.data.set_mesh(newV, newF);
  // viewer.data.add_points(newV, Eigen::RowVector3d(1,0,0));
  viewer.data.add_edges(P1, P2, Eigen::RowVector3d(0.54, 0.47, 0.39)); // brown for edges

  // white faces for joints
  Eigen::MatrixXd C = Eigen::MatrixXd::Constant(newF.rows(),3,1);
  //Eigen::MatrixXd C = Eigen::MatrixXd(newF.rows(),3);
  for(int sphere_idx = 0; sphere_idx < newF.rows() / 20; sphere_idx++){
    // lookup the cost of that sphere
    double cost = Cost(sphere_idx);
    for(int f_idx = 0; f_idx < 20; f_idx++){
      C.row((sphere_idx*20)+f_idx) *= cost;
    }
  }

  viewer.data.set_colors(C);
  viewer.core.show_lines = false; // don't show wireframe on joints
  // Launch a viewer instance
  viewer.launch();
  return 0;
}
