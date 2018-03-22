#include "compute_normal_sets.h"
#include <iostream>

typedef Eigen::Triplet<int> T;

std::vector<NormalSet> compute_normal_sets(const Eigen::MatrixXi &F, const Eigen::MatrixXd &V)
{
  // initialize normal sets
  std::vector<NormalSet> all_normal_sets;

  std::map<std::string, int> edge_to_face = preprocess_edge_to_face(F);

  Eigen::MatrixXd N;
  igl::per_face_normals(V, F, Eigen::Vector3d(1,1,1).normalized(), N);

  std::set<int> visited;
  int unseen_face_idx = 0;

  std::cout << "# Total Faces: " << F.rows() << std::endl;
  // BFS on entire F graph
  int id = 0;
  while(visited.size() != F.rows()){
    std::cout << "# Faces visited: " << visited.size() << std::endl;
    // Get a 'random' face_idx from to_visit
    // while unseen_face_idx is in visited, increment it
    while(visited.find(unseen_face_idx) != visited.end()){
      unseen_face_idx += 1;
    }
    int face_idx = unseen_face_idx;
    // Initialize set with unseen face in it
    NormalSet normalSet(face_idx, N.row(face_idx), id);
	id++;

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
  return all_normal_sets;
}

bool similar_normals(Eigen::Vector3d n1, Eigen::Vector3d n2){
  double threshold = 0.866; // 30 degrees
  //double threshold = 0.8;

  // TODO: sanity check this works ok
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

std::map<std::string, int> preprocess_edge_to_face( const Eigen::MatrixXi &F ){
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

bool sharedBoundary(Eigen::VectorXi bnd1, Eigen::VectorXi bnd2, std::vector<int> &endpoints, std::set<int> &foundSharedVertices) {
	for (int i = 0; i < bnd1.size(); i++) {
		// Look for bnd1(i) inside of bnd2 if bnd1(i) is not already in a shared boundary
		if (foundSharedVertices.find(bnd1(i)) == foundSharedVertices.end()) {
			for (int j = 0; j < bnd2.size(); j++) {
				if (bnd1(i) == bnd2(j)) {
					//std::cout << "checking" << std::endl;
					int endpoint1 = bnd1(i);
					int endpoint2 = bnd1(i);
					int f = 1;
					// fixed code, please double check if it works (I checked with the cube and it works)
					while (f < bnd1.size() && bnd1((i + f) % bnd1.size()) == bnd2((j - f + bnd2.size()) % bnd2.size())) {
						endpoint2 = bnd2((j - f + bnd2.size()) % bnd2.size());
						foundSharedVertices.insert(endpoint2);
						f++;
					}
					// Original code
					/*
					while (i + f < bnd1.size() && j - f >= 0 && bnd1(i + f) == bnd2(j - f)) {
						std::cout << "checking 1" << std::endl;
						endpoint2 = bnd2(j - f);
						foundSharedVertices.insert(endpoint2);
						f++;
					}
					int b = 1;
					while (i - b >= 0 && j + b < bnd2.size() && bnd1(i - b) == bnd2(j + b)) {
						std::cout << "checking 2" << std::endl;
						endpoint1 = bnd2(j + b);
						foundSharedVertices.insert(endpoint1);
						b++;
					}
					*/
					if (endpoint1 != endpoint2){
						//std::cout << "checking 3" << std::endl;
						endpoints.push_back(endpoint1);
						endpoints.push_back(endpoint2);
						foundSharedVertices.erase(endpoint1);
						foundSharedVertices.erase(endpoint2);
						//hasSharedBoundary = true;
						//i = f; // double check
					}
					//else {
						//std::cout << "checking erase" << std::endl;
						//foundSharedVertices.erase(endpoint1);
					//}
				}
			}
		}
	}
	return (endpoints.size() > 0);
		//hasSharedBoundary;
}


void straightenEdges(Eigen::MatrixXd &V, Eigen::MatrixXi &F, std::vector<NormalSet> &normal_sets, Eigen::MatrixXd &newV, Eigen::MatrixXi &newF) {
	// Initialize boundaries
	for (std::vector<NormalSet>::iterator set = normal_sets.begin(); set != normal_sets.end(); set++) {
		std::set<int> normal_set = (*set).face_set;
		Eigen::MatrixXi F_set(normal_set.size(), 3);
		std::set<int>::iterator face;
		int i = 0;
		for (face = normal_set.begin(); face != normal_set.end(); ++face) {
			int face_idx = *face;
			F_set.row(i) = F.row(face_idx);
			i++;
		}
		Eigen::VectorXi bnd;
		igl::boundary_loop(F_set, bnd);
		(*set).addBoundary(bnd);
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
			//std::cout << "set1 id" << std::endl;
			//std::cout << set1.id << std::endl;
			//std::cout << "set2 id" << std::endl;
			//std::cout << set2.id << std::endl;
			if (set1.id != set2.id) {
				// Find shared boundary
				std::vector<int> endpoints;
				//std::cout << "here 1" << std::endl;
				if (sharedBoundary(set1.bnd, set2.bnd, endpoints, foundSharedVertices)) {
					//std::cout << "here 2" << std::endl;
					for (int i = 0; i < endpoints.size(); i++) {
						//std::cout<< endpoints[i] <<std::endl;
						boundingVertices.insert(endpoints[i]);
						newVertices.insert(endpoints[i]);
					}
					//F_size += (boundingVertices.size() - 2);
				}
			}
		}
		if (boundingVertices.size() > 2) {// Fixed by adding this condition, still need to figure out why there are cases = 2?
			F_size += (boundingVertices.size() - 2);
			(*iter1).simplifyBoundary(boundingVertices);
			//std::cout << "check here" << std::endl;
			//std::cout << (*iter1).simplified_bnd << std::endl; // Add another attribute to the class to store simplified_bnd;
															   //can't directly update the original bnd becasue we still need to use it for comparison (iter2)
		}
	}

	std::cout << "FindAllNewVertices" << std::endl;
	// Update V
	newV.resize(newVertices.size(), 3);
	std::set<int>::iterator veriter;
	int i = 0;
	for (veriter = newVertices.begin(); veriter != newVertices.end(); ++veriter) {
		int v_idx = *veriter;
		newV.row(i) = V.row(v_idx);
		i++;
	}

	std::cout << "updatedV" << std::endl;
	// Triangulate to make new F and update face_set in normalSet
	newF.resize(F_size, 3);
	// This is the final version but it still does not quite work; also the coloring of faces is another problem (I ran it with release mode to bypass that problem)
  int F_idx = 0;
	for (std::vector<NormalSet>::iterator iter = normal_sets.begin(); iter != normal_sets.end(); iter++) {
		NormalSet cur_set = *iter;
		cur_set.clearSet();
		for (int i = 2; i < cur_set.simplified_bnd.size(); i++) { // Is it safe??
			if (cur_set.simplified_bnd.size() > 0) {
				newF(F_idx, 0) = cur_set.simplified_bnd(0);
				newF(F_idx, 1) = cur_set.simplified_bnd(i - 1);
				newF(F_idx, 2) = cur_set.simplified_bnd(i);
        cur_set.addToSet(F_idx, Eigen::Vector3d::Zero()); // 0 normal vector for now
				F_idx++;
			}
		}
    *iter = cur_set;
	}

	// One with holes
	/*
	for (std::set<std::vector<int>>::iterator iter = bnd_sets.begin(); iter != bnd_sets.end(); iter++) {
		std::vector<int> cur_bnd = *iter;
		std::cout << "here" << std::endl;
		std::cout << cur_bnd.size() << std::endl;
		for (int i = 2; i < cur_bnd.size(); i++) { // Is it safe??
			newF(F_idx, 0) = cur_bnd[0];
			newF(F_idx, 1) = cur_bnd[i - 1];
			newF(F_idx, 2) = cur_bnd[i];
			//cur_set.face_set.insert(F_idx);
			F_idx++;
			//std::cout << F_idx << std::endl;
		}
	}*/

	// Does not work
	/*
	for (std::vector<NormalSet>::iterator iter = normal_sets.begin(); iter != normal_sets.end(); iter++) {
		NormalSet cur_set = *iter;
		cur_set.face_set.clear();
		std::cout << "here" << std::endl;
		std::cout << cur_set.bnd.size() << std::endl;
		for (int i = 2; i < cur_set.bnd.size(); i++) { // Is it safe??
			newF(F_idx, 0) = cur_set.bnd(0);
			newF(F_idx, 1) = cur_set.bnd(i - 1);
			newF(F_idx, 2) = cur_set.bnd(i);
			cur_set.face_set.insert(F_idx);
			F_idx++;
			//std::cout << F_idx << std::endl;
		}
	}
	*/
	std::cout << "updatednewF" << std::endl;
	//F = newF; // Check if it works?
}
