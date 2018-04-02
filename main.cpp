#include "Edge.h"
#include "euler_characteristic.h"
#include <igl/read_triangle_mesh.h>
#include <igl/viewer/Viewer.h>
#include "compute_normal_sets.h"
#include <normalSet.h>
#include <vector>
#include <igl/unproject_onto_mesh.h>
#include <igl/upsample.h>
#include <igl/qslim.h>
#include <igl/doublearea.h>
#include <igl/readOFF.h>
#include <fstream>

void view_straightened_mesh(Eigen::MatrixXd &V, Eigen::MatrixXi &F, std::vector<NormalSet> &normal_sets, igl::viewer::Viewer &viewer, Eigen::MatrixXd& newVCenters, Eigen::MatrixXi& Edges);
void color_normal_sets(Eigen::MatrixXd &V, Eigen::MatrixXi &F, std::vector<NormalSet> &normal_sets, igl::viewer::Viewer &viewer);
void preprocess_mesh(Eigen::MatrixXd &V, Eigen::MatrixXi &F);
Eigen::Vector3d randcolor();

int main(int argc, char *argv[])
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  Eigen::MatrixXd originalV;
  Eigen::MatrixXi originalF;
  // Keep V, F used for merging
  Eigen::MatrixXd mergeV;
  Eigen::MatrixXi mergeF;

  // Load in a mesh
  //igl::read_triangle_mesh(argc>1 ? argv[1] : "../shared/data/cube.obj", V, F);
  //igl::read_triangle_mesh(argc>1 ? argv[1] : "../shared/data/max-face-low-res.obj", V, F);
  igl::read_triangle_mesh(argc>1 ? argv[1] : "../shared/data/bunny.off", V, F);
  //igl::read_triangle_mesh(argc > 1 ? argv[1] : "../shared/data/cylinder.obj", V, F);
  //igl::read_triangle_mesh(argc > 1 ? argv[1] : "../shared/data/spot_triangulated.obj", V, F);
  //igl::read_triangle_mesh(argc > 1 ? argv[1] : "../shared/data/spot_control_mesh.obj", V, F);
  //igl::read_triangle_mesh("../shared/data/teapot-low.obj", V, F);
  //igl::read_triangle_mesh("../shared/data/teapot.obj", V, F);
  int num_regions = 80;

  // Keep the original mesh
  originalV = V;
  originalF = F;

  // Compute normal and double area of each face
  Eigen::MatrixXd N;
  igl::per_face_normals(V, F, Eigen::Vector3d(1,1,1).normalized(), N);
  Eigen::VectorXd doubA;
  igl::doublearea(V, F, doubA);

  // Initialization
  std::vector<NormalSet> normal_sets;
  std::set<int> painted_faces;
  NormalSet painting_set;
  std::vector<int> fid;

  // Launch a viewer instance
  igl::viewer::Viewer viewer;

  // Set the vertices and faces for the viewer
  viewer.data.set_mesh(V, F);
  // Color matrix and vector
  Eigen::MatrixXd C = Eigen::MatrixXd::Constant(F.rows(), 3, 1);
  Eigen::Vector3d painting_color;

  // final outputs for writing
  Eigen::MatrixXd newVCenters;
  Eigen::MatrixXi E;

  bool painting = false;
  viewer.callback_key_pressed =
    [&](igl::viewer::Viewer & viewer, unsigned char key, int mod)->bool
   {
    switch(key)
    {
      default:
        return false;
      case 'p':
      {
    		Eigen::MatrixXd PC = Eigen::MatrixXd::Constant(F.rows(), 3, 1);
    		viewer.data.set_colors(PC);
    		painting = !painting; // toggle
    		painting_color = randcolor();
        painting_set = NormalSet(true);
        break;
      }
  	  case 'b': {
  		  if (painting == true && fid.size()) {
  			  int rm_fid = fid.back();
  			  painting_set.erase(rm_fid, N.row(rm_fid), doubA(rm_fid) / 2.0);
  			  painted_faces.erase(rm_fid);
  			  C.row(rm_fid) = Eigen::Vector3d::Constant(1);
  			  fid.pop_back();
  			  viewer.data.set_colors(C);
  		  }
  		  break;
  	  }
      case 'd':
      {
  		if (painting == true) {
  		    normal_sets.push_back(painting_set);
  		    std::cout << "Done Painting Set " << painting_set.id << std::endl;
    			for (std::set<int>::iterator iter = painting_set.face_set.begin(); iter != painting_set.face_set.end(); iter++) {
    			    std::cout << *iter << std::endl;
    			}
    			std::cout << "normal_sets size: " << normal_sets.size() << std::endl;

    			// initialize new painting_set
    			painting_color = randcolor();
    			painting_set = NormalSet(true);
    		}
          break;
      }
      case 'w':
      {
        std::ofstream myfile;
        myfile.open("output_bunny.txt");
        myfile << "V: " << newVCenters.rows() << "\n";
        for(int v_idx = 0; v_idx < newVCenters.rows(); v_idx++){
          myfile << V(v_idx, 0) << ", " << V(v_idx, 1) << ", " << V(v_idx,2) << "\n";
        }
        myfile << "E: " << E.rows() << "\n";
        for(int e_idx = 0; e_idx < E.rows(); e_idx++){
          myfile << E(e_idx, 0) << ", " << E(e_idx, 1) << "\n";
        }
        myfile.close();
        break;
      }
      case 'n':
      {
        compute_normal_sets(F, V, normal_sets, painted_faces);
		    color_normal_sets(V, F, normal_sets, viewer);

        mergeV = V;
        mergeF = F;
        painting = false;
        break;
      }
      case 's':
      {
        std::cout << "num_regions: " << num_regions << std::endl;
        mergeNormalSets(mergeV, mergeF, normal_sets, num_regions);
        color_normal_sets(mergeV, mergeF, normal_sets, viewer);

        view_straightened_mesh(V, F, normal_sets, viewer, newVCenters, E);
        num_regions = fmax(num_regions - 5, 1); // don't go negative
        painting = false;

        break;
      }
      case 'q':
      {
        preprocess_mesh(V, F);
        viewer.data.clear();
        viewer.data.set_mesh(V, F);

    		// re-initialize everything
    		C = Eigen::MatrixXd::Constant(F.rows(), 3, 1);
    		igl::per_face_normals(V, F, Eigen::Vector3d(1, 1, 1).normalized(), N);
    		igl::doublearea(V, F, doubA);
    		painting = false;
    		break;
      }
  	  case 'r':
  	  {
  		  V = originalV;
  		  F = originalF;

  		  viewer.data.clear();
  		  viewer.core.show_lines = true;
  		  viewer.data.set_mesh(V, F);
  		  // Clean everything
  		  normal_sets.clear();
  		  painted_faces.clear();
  		  painting_set.clearSet();
        num_regions = 80;
  		  C = Eigen::MatrixXd::Constant(F.rows(), 3, 1);
  		  painting = false;
  		  break;
  	  }
    }
    return true;
  };

  painting_color = randcolor();
  viewer.callback_mouse_down =
	  [&V, &F, &painting, &C, &normal_sets, &painting_set, &painting_color, &N, &painted_faces, &doubA, &fid]
  (igl::viewer::Viewer& viewer, int, int)->bool
  {
	  if (painting) {
		  int cur_fid;
		  Eigen::Vector3f bc;
		  // Cast a ray in the view direction starting from the mouse position
		  double x = viewer.current_mouse_x;
		  double y = viewer.core.viewport(3) - viewer.current_mouse_y;
		  if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core.view * viewer.core.model, viewer.core.proj, viewer.core.viewport, V, F, cur_fid, bc)) {
			  painting_set.addToSet(cur_fid, N.row(cur_fid), doubA(cur_fid) / 2.0);
			  painted_faces.insert(cur_fid);
			  C.row(cur_fid) = painting_color;
			  fid.push_back(cur_fid);
			  viewer.data.set_colors(C);
		  }
	  }
	  return false; // regular 3D manipulations
  };

  //Instructions for viewer
  std::cout<<R"(
(step 0: press 'q' preprocess mesh by decimating it)
step 1: press 'p' activate painting mode
		- click on face to paint face
		- press 'd' to start a new group
		- press 'b' to unpaint the last painted face
step 2: press 'n' compute normal set grouping (based on similar normal)
step 3: press 's' merge normal set
step 4: press 's' show simplified result
press 'r' reset
  )";

  viewer.launch();
  return 0;
}

Eigen::Vector3d randcolor() {
	double r = ((double)rand() / (RAND_MAX));
	double g = ((double)rand() / (RAND_MAX));
	double b = ((double)rand() / (RAND_MAX));
	return Eigen::Vector3d(r, g, b);
}

void preprocess_mesh(Eigen::MatrixXd &V, Eigen::MatrixXi &F){
  // upsample then downsample with qslim
  // Eigen::MatrixXd NV;
  // Eigen::MatrixXi NF;
  // igl::upsample(V, F, NV, NF);
  Eigen::VectorXi J, I;
  igl::qslim(V, F, F.rows()/2, V, F, J, I);
}

void color_normal_sets(Eigen::MatrixXd &V, Eigen::MatrixXi &F, std::vector<NormalSet> &normal_sets, igl::viewer::Viewer &viewer){
  // Set the vertices and faces for the viewer
  viewer.data.clear();
  viewer.data.set_mesh(V, F);

  Eigen::MatrixXd C = Eigen::MatrixXd::Constant(F.rows(),3,1);
  // each normal set has a different color
  for(std::vector<NormalSet>::iterator set = normal_sets.begin(); set != normal_sets.end(); set++){
    std::set<int> normal_set = (*set).face_set;
	  Eigen::Vector3d set_color = randcolor();
    std::set<int>::iterator face;
    for (face = normal_set.begin(); face != normal_set.end(); ++face){
        int face_idx = *face;
        C.row(face_idx) = set_color;
    }
  }

  viewer.data.set_colors(C);
}

void view_straightened_mesh(Eigen::MatrixXd &V, Eigen::MatrixXi &F, std::vector<NormalSet> &normal_sets, igl::viewer::Viewer &viewer, Eigen::MatrixXd& newVCenters, Eigen::MatrixXi& Edges){
  // Straighten edges
  Eigen::MatrixXd newV;
  Eigen::MatrixXi newF;
  Eigen::MatrixXd P1, P2;
  Eigen::VectorXd Cost;
  std::cout << "gets here" << std::endl;
  straightenEdges(V, F, normal_sets, newV, newF, P1, P2, Cost, newVCenters, Edges);

  // Create a libigl Viewer object
  viewer.data.clear();
  viewer.data.set_mesh(newV, newF);
  // viewer.data.add_points(newV, Eigen::RowVector3d(1,0,0));
  //viewer.data.add_edges(P1, P2, Eigen::RowVector3d(0.54, 0.47, 0.39)); // brown for edges
  viewer.data.add_edges(P1, P2, Eigen::RowVector3d(0,0,0));

  // white faces for joints
  Eigen::MatrixXd C = Eigen::MatrixXd::Constant(newF.rows(),3,1);
  viewer.data.set_colors(C);
  viewer.core.show_lines = false; // don't show wireframe on joints
  std::cout << "Num Joints: " << newV.rows() / 12 << std::endl;
  std::cout << "Num Rods: " << P1.rows() << std::endl;
}
