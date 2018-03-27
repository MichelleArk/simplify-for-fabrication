#include "compute_normal_sets.h"
#include <iostream>

typedef Eigen::Triplet<int> T;

void compute_normal_sets(const Eigen::MatrixXi &F, const Eigen::MatrixXd &V, std::vector<NormalSet>& all_normal_sets)
{
  // initialize normal sets
  //std::vector<NormalSet> all_normal_sets;

  std::map<std::string, int> edge_to_face = preprocess_edge_to_face(F);

  Eigen::MatrixXd N;
  igl::per_face_normals(V, F, Eigen::Vector3d(1,1,1).normalized(), N);

  std::set<int> visited;
  int unseen_face_idx = 0;

  std::cout << "# Total Faces: " << F.rows() << std::endl;
  // BFS on entire F graph
  //int id = 0;
  while(visited.size() != F.rows()){
    std::cout << "# Faces visited: " << visited.size() << std::endl;
    // Get a 'random' face_idx from to_visit
    // while unseen_face_idx is in visited, increment it
    while(visited.find(unseen_face_idx) != visited.end()){
      unseen_face_idx += 1;
    }
    int face_idx = unseen_face_idx;
    // Initialize set with unseen face in it
    NormalSet normalSet(face_idx, N.row(face_idx));
	  //id++;

    // Initialize Q with random unseen face on it
    std::queue<int> Q;
    Q.push(face_idx);

    // BFS on a single 'connected' componenet
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
        if( similar_normals(neighbour_normal, normalSet.avg_normal) && !(visited.find(neighbour) != visited.end()) ){
          Q.push(neighbour);
          normalSet.addToSet(neighbour, neighbour_normal);

          //push onto visited for cycle detection / avoidance
          visited.insert(neighbour);
        }
      }
    }
    all_normal_sets.push_back(normalSet);
  }
  //return all_normal_sets;
}

bool similar_normals(Eigen::Vector3d n1, Eigen::Vector3d n2)
{
  double threshold = 0.866; // 30 degrees
  //double threshold = 0.8;
  n1.normalize();
  n2.normalize();
  return n1.dot(n2) > threshold;
}

std::vector<int> get_neighbours(const Eigen::MatrixXi &F, int f_idx, std::map<std::string,int> edge_to_f)
{
  std::vector<int> neighbours;

  int v0 = F(f_idx, 0); int v1 = F(f_idx, 1); int v2 = F(f_idx, 2);

  std::string key10 = std::to_string(v1) + "," + std::to_string(v0);
  std::string key21 = std::to_string(v2) + "," + std::to_string(v1);
  std::string key02 = std::to_string(v0) + "," + std::to_string(v2);

  int n1 = edge_to_f[key10];
  int n2 = edge_to_f[key21];
  int n3 = edge_to_f[key02];

  // TODO: double-check that this is the right way to check if n_i was found
  if(n1 != 0){
    neighbours.push_back(n1);
  }
  if(n2 != 0){
    neighbours.push_back(n2);
  }
  if(n3 != 0){
    neighbours.push_back(n3);
  }
  return neighbours;
}

std::map<std::string, int> preprocess_edge_to_face( const Eigen::MatrixXi &F )
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

bool sharedBoundary(Eigen::VectorXi bnd1, Eigen::VectorXi bnd2, std::vector<int> &endpoints, std::set<int> &foundSharedVertices)
{
	for (int i = 0; i < bnd1.size(); i++) {
		// Look for bnd1(i) inside of bnd2 if bnd1(i) is not already in a shared boundary
		if (foundSharedVertices.find(bnd1(i)) == foundSharedVertices.end()) {
			for (int j = 0; j < bnd2.size(); j++) {
				if (bnd1(i) == bnd2(j)) {
					int endpoint1 = bnd1(i);
					int endpoint2 = bnd1(i);
					int f = 1;
					// fixed code, please double check if it works (I checked with the cube and it works)
					while (f < bnd1.size() && bnd1((i + f) % bnd1.size()) == bnd2((j - f + bnd2.size()) % bnd2.size())) {
						endpoint2 = bnd2((j - f + bnd2.size()) % bnd2.size());
						foundSharedVertices.insert(endpoint2);
						f++;
					}
					if (endpoint1 != endpoint2){
						endpoints.push_back(endpoint1);
						endpoints.push_back(endpoint2);
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

void straightenEdges(Eigen::MatrixXd &V, Eigen::MatrixXi &F, std::vector<NormalSet> &normal_sets, Eigen::MatrixXd &newV, Eigen::MatrixXi &newF, Eigen::MatrixXd &P1, Eigen::MatrixXd &P2, Eigen::VectorXd& Cost)
{
	// Initialize boundaries
	for (std::vector<NormalSet>::iterator set = normal_sets.begin(); set != normal_sets.end(); set++) {
		(*set).computeBoundary(F);
	}

	// Find longest shared boundaries
	std::set<int> newVertices;
	std::set<int> foundSharedVertices;
	int F_size = 0;
	for (std::vector<NormalSet>::iterator iter1 = normal_sets.begin(); iter1 != normal_sets.end(); iter1++) {
		std::set<int> boundingVertices;
		NormalSet set1 = *iter1;
		for (std::vector<NormalSet>::iterator iter2 = normal_sets.begin(); iter2 != normal_sets.end(); iter2++) {
			NormalSet set2 = *iter2;
			if (set1.id != set2.id) {
				// Find shared boundary
				std::vector<int> endpoints;
				if (sharedBoundary(set1.bnd, set2.bnd, endpoints, foundSharedVertices)) {
					for (int i = 0; i < endpoints.size(); i++) {
						boundingVertices.insert(endpoints[i]);
						newVertices.insert(endpoints[i]);
					}
				}
			}
		}
		if (boundingVertices.size() > 2) {// Fixed by adding this condition, still need to figure out why there are cases = 2?
			F_size += (boundingVertices.size() - 2);
      //can't directly update the original bnd becasue we still need to use it for comparison (iter2)
			set1.simplifyBoundary(boundingVertices);
      *iter1 = set1;
		}
	}

  // make newVertices a vector, not a set
  std::vector<int> newVerticesVector;
  newVerticesVector.assign(newVertices.begin(), newVertices.end());
  // newV and newF store V, F to represent icosahedrons centered at newVertices
  createApproxSpheres(newVerticesVector, V, newV, newF);
  // P1, P2 store edges between centers in newVertices based on normal_sets
  connectApproxSpheres(normal_sets, V, P1, P2);
  // compute cost per vertex
  computeRemovalCostPerVertex(newVerticesVector, V, normal_sets, Cost);

  // update normal sets and vertices whose cost are not qualified
  std::vector<int> newVerticesCopy(newVerticesVector);
  int Cost_idx = 0;
  for (std::vector<int>::iterator iter = newVerticesCopy.begin(); iter != newVerticesCopy.end(); iter++) {
    int cur_V = *iter;
    if (Cost(Cost_idx) < 0.6) {// doublecheck the threshold
      removeVertex(newVerticesVector, normal_sets, cur_V);
      createApproxSpheres(newVerticesVector, V, newV, newF);
      connectApproxSpheres(normal_sets, V, P1, P2);
      std::cout << newF.size() << std::endl;
    }
    Cost_idx++;
  }
  computeRemovalCostPerVertex(newVerticesVector, V, normal_sets, Cost);
}

void createApproxSpheres(std::vector<int> icoCenters, Eigen::MatrixXd &V, Eigen::MatrixXd &newV, Eigen::MatrixXi &newF)
{
  Eigen::MatrixXd icoV; Eigen::MatrixXi icoF;
  igl::read_triangle_mesh("../shared/data/icosahedron.obj", icoV, icoF);
  int v_step = icoV.rows(); int f_step = icoF.rows();
  icoV *= 2; // scale
  newV.resize(v_step * icoCenters.size(),3);
  newF.resize(f_step * icoCenters.size(),3);

  std::vector<int>::iterator iter;
  int v_i = 0; int f_i = 0;
  // for each center, add the appropriate vertices and faces
  int num_centers_seen = 0;
  for (iter = icoCenters.begin(); iter != icoCenters.end(); ++iter) {
  	int c_idx = *iter;
  	Eigen::Vector3d center = V.row(c_idx);
    //Joint new_Joint(c_idx, center, v_i, f_i);
    for(int v = 0; v < v_step; v++){
      //icoV vertices translated by center
      newV(v_i,0) = center(0) + icoV(v,0);
      newV(v_i,1) = center(1) + icoV(v,1);
      newV(v_i,2) = center(2) + icoV(v,2);
      v_i++;
    }
    for(int f = 0; f < f_step; f++){
      //update indices to be global
      newF(f_i,0) = icoF(f,0) + (num_centers_seen * v_step);
      newF(f_i,1) = icoF(f,1) + (num_centers_seen * v_step);
      newF(f_i,2) = icoF(f,2) + (num_centers_seen * v_step);
      f_i++;
    }
    num_centers_seen += 1;
  }
}

void connectApproxSpheres(std::vector<NormalSet> &normal_sets, Eigen::MatrixXd &V, Eigen::MatrixXd &P1, Eigen::MatrixXd &P2)
{
  // This is very hacky and needs to be written better. twice as slow as it should be
  int E_size = 0;
  Eigen::MatrixXi seenEdges = Eigen::MatrixXi::Zero(V.rows(), V.rows());

  typedef Eigen::Triplet< double > Triplet;
  std::vector< double > P1triplets; // to store connecting triangles
  std::vector< double > P2triplets;

  for (std::vector<NormalSet>::iterator iter = normal_sets.begin(); iter != normal_sets.end(); iter++) {
    NormalSet cur_set = *iter;
    for (int i = 0; i < cur_set.simplified_bnd.size(); i++) {
      int v1 = cur_set.simplified_bnd(i);
      int v2 = cur_set.simplified_bnd((i+1) % cur_set.simplified_bnd.size());

      if(!seenEdges(v1, v2)){
        P1triplets.push_back(V(v1,0));
        P1triplets.push_back(V(v1,1));
        P1triplets.push_back(V(v1,2));

        P2triplets.push_back(V(v2,0));
        P2triplets.push_back(V(v2,1));
        P2triplets.push_back(V(v2,2));

        // Mark the edge as seen in either direction
        seenEdges(v1, v2) = 1;
        seenEdges(v2, v1) = 1;
        E_size += 1;
      }
    }
  }

  P1.resize(E_size, 3);
  P2.resize(E_size, 3);
  for(int e_idx = 0; e_idx < (E_size*3); e_idx+=3){
    P1.row(e_idx/3) = Eigen::Vector3d(P1triplets[e_idx], P1triplets[e_idx+1], P1triplets[e_idx+2]);
    P2.row(e_idx/3) = Eigen::Vector3d(P2triplets[e_idx], P2triplets[e_idx+1], P2triplets[e_idx+2]);
  }
}

// indicies into V of new vertices
void computeRemovalCostPerVertex(std::vector<int> newVertices, Eigen::MatrixXd V, std::vector<NormalSet> normal_sets, Eigen::VectorXd &C){
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
    C(cost_idx) = abs(C(cost_idx));
    cost_idx++; // could be a point of error..
  }
}

// Remove unimportant vertex
void removeVertex(std::vector<int> &newVertices, std::vector<NormalSet> &normal_sets, int v_idx) {
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
