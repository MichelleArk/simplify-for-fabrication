#include "compute_normal_sets.h"
#include <iostream>
#include <igl/doublearea.h>

typedef Eigen::Triplet<int> T;

void compute_normal_sets(
  const Eigen::MatrixXi &F,
  const Eigen::MatrixXd &V,
  std::vector<NormalSet>& all_normal_sets, std::set<int>& painted_faces)
{
  std::map<std::string, int> edge_to_face = preprocess_edge_to_face(F);

  // Compute per-face normals
  Eigen::MatrixXd N;
  igl::per_face_normals(V, F, Eigen::Vector3d(1,1,1).normalized(), N);

  // Compute doublearea
  Eigen::VectorXd doubA;
  igl::doublearea(V, F, doubA);
  double total_area = doubA.sum()/2.0;


  // Initialize visited set to empty
  std::set<int> visited;
  // For obtaining next seed
  int unseen_face_idx = 0;

  // BFS on entire F graph
  std::cout << "# Total Faces: " << F.rows() << std::endl;
  while(visited.size() != F.rows()){
    std::cout << "# Faces visited: " << visited.size() << std::endl;
    // FIND NEW SEED: Get a 'random' face_idx from to_visit
    // while unseen_face_idx is in visited or in painted, increment it
    while((visited.find(unseen_face_idx) != visited.end()) || (painted_faces.find(unseen_face_idx) != painted_faces.end())){
      unseen_face_idx += 1;
    }
    int face_idx = unseen_face_idx;
    if(face_idx >= F.rows()){
      break;
    }

    // Initialize set with unseen face in it
    NormalSet normalSet(face_idx, N.row(face_idx), doubA(face_idx)/2.0);

    // Initialize Q with random unseen face on it
    std::queue<int> Q;
    Q.push(face_idx);

    // BFS on a single connected componenet
    std::set<int> already_seen;
    // push onto visited for cycle detection / avoidance
    visited.insert(face_idx);

    while(!Q.empty()){
      // Get current face
      int curr_face = Q.front();
      Q.pop();

      // Iterate over current face's neighbours
      std::vector<int> neighbours = get_neighbours(F, curr_face, edge_to_face);
      for(std::vector<int>::iterator itr = neighbours.begin(); itr != neighbours.end(); itr++){
        int neighbour = *itr;
        Eigen::Vector3d neighbour_normal = N.row(neighbour);
        // Check that neighbour has similar normal and has not been on the queue already
        if(similar_normals(neighbour_normal, normalSet.avg_normal) && !(visited.find(neighbour) != visited.end()) ){
          Q.push(neighbour);
          if(painted_faces.find(neighbour) == painted_faces.end()){ // neighbor is not painted
            normalSet.addToSet(neighbour, neighbour_normal, doubA(neighbour)/2.0);
          } else{ // neighbor was painted - consider the normal but do not add to set, since it is already claimed by painted region
            normalSet.updateAvgNormal(neighbour_normal);
          }
          //push onto visited for cycle detection / avoidance
          visited.insert(neighbour);
        }
      }
    }
    all_normal_sets.push_back(normalSet);
  }
  for (std::vector<NormalSet>::iterator set = all_normal_sets.begin(); set != all_normal_sets.end(); set++) {
    (*set).computeBoundary(F, V);
  }
}

bool similar_normals(
  Eigen::Vector3d n1,
  Eigen::Vector3d n2)
{
  double threshold = 0.95; // 30 degrees
  n1.normalize();
  n2.normalize();
  return n1.dot(n2) > threshold;
}

std::vector<int> get_neighbours(
  const Eigen::MatrixXi &F,
  int f_idx,
  std::map<std::string,int> edge_to_f)
{
  std::vector<int> neighbours;

  int v0 = F(f_idx, 0); int v1 = F(f_idx, 1); int v2 = F(f_idx, 2);

  std::string key10 = std::to_string(v1) + "," + std::to_string(v0);
  std::string key21 = std::to_string(v2) + "," + std::to_string(v1);
  std::string key02 = std::to_string(v0) + "," + std::to_string(v2);

  // add to neighbours if key is found
  if(edge_to_f.find(key10) != edge_to_f.end()){
    neighbours.push_back(edge_to_f[key10]);
  }
  if(edge_to_f.find(key21) != edge_to_f.end()){
    neighbours.push_back(edge_to_f[key21]);
  }
  if(edge_to_f.find(key02) != edge_to_f.end()){
    neighbours.push_back(edge_to_f[key02]);
  }
  return neighbours;
}

std::map<std::string, int> preprocess_edge_to_face(const Eigen::MatrixXi &F )
{
  std::map<std::string, int> edge_to_f;
  int num_faces = F.rows();
  for(int f_idx= 0; f_idx < num_faces; f_idx++){
    int v0 = F(f_idx, 0); int v1 = F(f_idx, 1); int v2 = F(f_idx, 2);
    // Create string keys
    std::string key01 = std::to_string(v0) + "," + std::to_string(v1);
    std::string key12 = std::to_string(v1) + "," + std::to_string(v2);
    std::string key20 = std::to_string(v2) + "," + std::to_string(v0);
    // Add mappings to edge_to_f
    edge_to_f.insert(make_pair(key01,f_idx));
    edge_to_f.insert(make_pair(key12,f_idx));
    edge_to_f.insert(make_pair(key20,f_idx));
  }
  return edge_to_f;
}

bool sharedBoundary(
  Eigen::VectorXi bnd1,
  Eigen::VectorXi bnd2,
  bool set1_painted,
  bool set2_painted,
  std::vector<int> &endpoints,
  std::set<int> &foundSharedVertices,
  Eigen::MatrixXd &V,
  double &shared_bnd_length)
{
	shared_bnd_length = 0;
	std::set<int> localFoundSharedVertices;
	for (int i = 0; i < bnd1.size(); i++) {
		// Look for bnd1(i) inside of bnd2 if bnd1(i) is not already in a shared boundary
		if (foundSharedVertices.find(bnd1(i)) == foundSharedVertices.end() && localFoundSharedVertices.find(bnd1(i)) == localFoundSharedVertices.end()) {
			for (int j = 0; j < bnd2.size(); j++) {
				if (bnd1(i) == bnd2(j)) {
					int endpoint1 = bnd1(i);
					int endpoint2 = bnd1(i);
					int f = 1;
					int prev_endpoint = endpoint1;
					// fixed code, please double check if it works (I checked with the cube and it works)
					while (f < bnd1.size() && bnd1((i + f) % bnd1.size()) == bnd2((j - f + bnd2.size()) % bnd2.size())) {
						endpoint2 = bnd2((j - f + bnd2.size()) % bnd2.size());
						shared_bnd_length += ((V.row(prev_endpoint) - V.row(endpoint2)).norm());
						prev_endpoint = endpoint2;
						foundSharedVertices.insert(endpoint2);
						localFoundSharedVertices.insert(endpoint2);
						if(set2_painted || set1_painted) {
							endpoints.push_back(endpoint2);
						}
						f++;
					}
					// go backwards
					int b = 1;
					prev_endpoint = endpoint1;
					while (b < bnd1.size() && bnd1((i - b + bnd1.size()) % bnd1.size()) == bnd2((j + b) % bnd2.size())) {
						endpoint1 = bnd2((j + b) % bnd2.size());
						shared_bnd_length += ((V.row(prev_endpoint) - V.row(endpoint1)).norm());
						prev_endpoint = endpoint1;
						foundSharedVertices.insert(endpoint1);
						localFoundSharedVertices.insert(endpoint1);
						if(set2_painted || set1_painted){
							endpoints.push_back(endpoint1);
						}
						b++;
					}
					if (endpoint1 != endpoint2){
						endpoints.push_back(endpoint1);
  					if(!set2_painted){ // to avoid pushing back endpoint2 twice
  						endpoints.push_back(endpoint2);
  					}
						foundSharedVertices.erase(endpoint1);
						foundSharedVertices.erase(endpoint2);
						//i = fmin(i+f+1, bnd1.size()-1); // double check
					}
				}
			}
		}
	}
  return (endpoints.size() > 0);
}

// indicies into V of new vertices, get min cost at some v to remove
void computeAngleDeficitPerVertex(
  std::vector<int> newVertices,
  Eigen::MatrixXd V,
  std::vector<NormalSet> normal_sets,
  Eigen::VectorXd &C,
  int& min_cost_vid,
  double& min_cost)
{
  // initialize all costs to be 2pi
  C = Eigen::VectorXd::Constant(newVertices.size(), 2 * M_PI);

  int cost_idx = 0;
  for (std::vector<int>::iterator v_itr= newVertices.begin(); v_itr != newVertices.end(); ++v_itr){
    int v_idx = *v_itr;
    // find every boundary loop that v_idx is part of, and compute internal angle in radians
    for (std::vector<NormalSet>::iterator n_iter = normal_sets.begin(); n_iter != normal_sets.end(); n_iter++) {
			NormalSet set = *n_iter;
      // loop over simplified_bnd of set
      int bnd_size = set.simplified_bnd.size();
      for(int i = 0; i < bnd_size; i++){
        if(set.simplified_bnd(i) == v_idx){
          // make triangle from i-1, i, i+1;
          int b_idx;
          if(i-1 == -1){
            b_idx = set.simplified_bnd(bnd_size-1);
          }else{
            b_idx = set.simplified_bnd(i-1);
          }
          int c_idx = set.simplified_bnd((i+1) % bnd_size);
          // get length of edges
          double bc = (V.row(b_idx) - V.row(c_idx)).norm();
          double cv = (V.row(c_idx) - V.row(v_idx)).norm();
          double vb = (V.row(v_idx) - V.row(b_idx)).norm();

          double bc_squared = bc * bc;
          double cv_squared = cv * cv;
          double vb_squared = vb * vb;
          // get internal angle in radians
          C(cost_idx) -= acos(
            (cv_squared + vb_squared - bc_squared)
            / (2 * cv * vb));
          // break out of loop that searches for v_idx in this simplified bnd
          break;
        }
      }
    }
    double v_cost = abs(C(cost_idx));
    C(cost_idx) = v_cost;
    if(min_cost_vid == -1 || v_cost < min_cost){
      min_cost_vid = v_idx;
      min_cost = v_cost;
    }
    cost_idx++;
  }
}
// Remove unimportant vertex
void removeVertex(
  std::vector<int> &newVertices,
  std::vector<NormalSet> &normal_sets,
  int v_idx)
{
	std::vector<int>::iterator position = std::find(newVertices.begin(), newVertices.end(), v_idx);
	if (position != newVertices.end()) {
		newVertices.erase(position);
	}
	for (std::vector<NormalSet>::iterator n_iter = normal_sets.begin(); n_iter != normal_sets.end(); ++n_iter) {
		NormalSet set = *n_iter;
		for (int i = 0; i < set.simplified_bnd.size(); i++) {
			Eigen::VectorXi new_bnd(set.simplified_bnd.size() - 1);
			bool update = false;
			int new_idx = 0;
			for (int bnd_idx = 0; bnd_idx < set.simplified_bnd.size(); bnd_idx++) {
				if (set.simplified_bnd(bnd_idx) == v_idx) {
					update = true;
				}
				else {
					if(new_idx < new_bnd.size())
						new_bnd(new_idx) = set.simplified_bnd(bnd_idx);
					new_idx++;
				}
			}
			if (update == true) {
				(*n_iter).simplified_bnd.resize(new_bnd.size());
				(*n_iter).simplified_bnd = new_bnd;
			}
		}
	}
}
