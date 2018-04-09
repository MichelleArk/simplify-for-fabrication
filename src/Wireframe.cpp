#include "Wireframe.h"

Wireframe::Wireframe(Eigen::MatrixXd V, Eigen::MatrixXi F){
  _V = V;
  _F = F;
}

void Wireframe::Update(std::vector<NormalSet> new_regions){
  regions = new_regions;
  straightenEdges();
}

void Wireframe::straightenEdges()
{
	// Initialize boundaries
	for (std::vector<NormalSet>::iterator set = regions.begin(); set != regions.end(); set++) {
		(*set).computeBoundary(_F, _V);
	}

	// Find longest shared boundaries
	std::set<int> newVertices;
	std::set<int> foundSharedVertices;
	for (std::vector<NormalSet>::iterator iter1 = regions.begin(); iter1 != regions.end(); iter1++) {
		std::set<int> boundingVertices;
		NormalSet set1 = *iter1;
		for (std::vector<NormalSet>::iterator iter2 = regions.begin(); iter2 != regions.end(); iter2++) {
			NormalSet set2 = *iter2;
			if (set1.id != set2.id) {
				// Find shared boundary
				std::vector<int> endpoints;
				double shared_bnd_length;
				if (sharedBoundary(set1.bnd, set2.bnd, set1.painted, set2.painted, endpoints, foundSharedVertices, _V, shared_bnd_length)) {
					for (int i = 0; i < endpoints.size(); i++) {
						boundingVertices.insert(endpoints[i]);
					}
				}
			}
		}
		if (boundingVertices.size() > 2) {// Fixed by adding this condition, still need to figure out why there are cases = 2?
      //can't directly update the original bnd becasue we still need to use it for comparison (iter2)
			set1.simplifyBoundary(boundingVertices);
      *iter1 = set1;
      // prevents adding 'floating' bounding vertices of size 2
      for (std::set<int>::iterator b_iter = boundingVertices.begin(); b_iter != boundingVertices.end(); b_iter++) {
        newVertices.insert(*b_iter);
      }
		} else{
      // no bounding vertices => this normal set is fully contained by another normal set
      (*iter1).simplified_bnd.resize(0);
      (*iter1).bnd.resize(0);
    }
	}

  // make newVertices a vector, not a set
  newVerticesVector.clear();
  newVerticesVector.assign(newVertices.begin(), newVertices.end());
  // newV and newF store V, F to represent icosahedrons centered at newVertices
  createApproxSpheres(); // some points in normal_sets, not in newVertices vector
  // P1, P2 store edges between centers in newVertices based on normal_sets
  connectApproxSpheres();
}

void Wireframe::createApproxSpheres()
{
  Eigen::MatrixXd icoV; Eigen::MatrixXi icoF;
  igl::read_triangle_mesh("../shared/data/icosahedron.obj", icoV, icoF);
  int v_step = icoV.rows(); int f_step = icoF.rows();
  //icoV *= 0.001; // scale for bunny
  icoV *= 0.002;
  //icoV *= 2; // scale for face
  //icoV *= 0.01; // scale for spot
  //icoV *= 0.2; // scale for teapot
  wire_V.resize(v_step * newVerticesVector.size(),3);
  wire_F.resize(f_step * newVerticesVector.size(),3);

  std::vector<int>::iterator iter;
  int v_i = 0; int f_i = 0;
  // for each center, add the appropriate vertices and faces
  int num_centers_seen = 0;
  for (iter = newVerticesVector.begin(); iter != newVerticesVector.end(); ++iter) {
  	int c_idx = *iter;
  	Eigen::Vector3d center = _V.row(c_idx);
    //Joint new_Joint(c_idx, center, v_i, f_i);
    for(int v = 0; v < v_step; v++){
      //icoV vertices translated by center
      wire_V(v_i,0) = center(0) + icoV(v,0);
      wire_V(v_i,1) = center(1) + icoV(v,1);
      wire_V(v_i,2) = center(2) + icoV(v,2);
      v_i++;
    }
    for(int f = 0; f < f_step; f++){
      //update indices to be global
      wire_F(f_i,0) = icoF(f,0) + (num_centers_seen * v_step);
      wire_F(f_i,1) = icoF(f,1) + (num_centers_seen * v_step);
      wire_F(f_i,2) = icoF(f,2) + (num_centers_seen * v_step);
      f_i++;
    }
    num_centers_seen += 1;
  }
}

void Wireframe::connectApproxSpheres()
{
  // This is very hacky and needs to be written better. twice as slow as it should be
  int E_size = 0;
  Eigen::MatrixXi seenEdges = Eigen::MatrixXi::Zero(_V.rows(), _V.rows());
  Eigen::MatrixXi seenEdgesRelnewV = Eigen::MatrixXi::Zero(newVerticesVector.size(), newVerticesVector.size());

  typedef Eigen::Triplet< double > Triplet;
  std::vector< double > P1triplets; // to store connecting triangles
  std::vector< double > P2triplets;

  for (std::vector<NormalSet>::iterator iter = regions.begin(); iter != regions.end(); iter++) {
    NormalSet cur_set = *iter;
    for (int i = 0; i < cur_set.simplified_bnd.size(); i++) {
      int v1 = cur_set.simplified_bnd(i);
      int v2 = cur_set.simplified_bnd((i+1) % cur_set.simplified_bnd.size());

      bool v1found = false;
      int v1foundidx = -1;
      bool v2found = false;
      int v2foundidx = -1;

      for(int v_idx = 0; v_idx < newVerticesVector.size(); v_idx++){
        if(newVerticesVector[v_idx] == v1){
          v1found = true;
          v1foundidx = v_idx;
        }
        if(newVerticesVector[v_idx] == v2){
          v2found = true;
          v2foundidx = v_idx;
        }
      }
      if(!v1found || !v2found){
        std::cout << "PROBLEM" << std::endl;
      }

      if(!seenEdges(v1, v2) && !seenEdges(v2, v1)){
        P1triplets.push_back(_V(v1,0));
        P1triplets.push_back(_V(v1,1));
        P1triplets.push_back(_V(v1,2));

        P2triplets.push_back(_V(v2,0));
        P2triplets.push_back(_V(v2,1));
        P2triplets.push_back(_V(v2,2));

        // Mark the edge as seen in either direction
        seenEdges(v1, v2) = 1;
        seenEdges(v2, v1) = 1;
        E_size += 1;

        // Mark as seen in rel indexing of icoCenters, ie newV
        seenEdgesRelnewV(v1foundidx, v2foundidx) = 1;
        seenEdgesRelnewV(v2foundidx, v1foundidx) = 1;
      }
    }
  }

  wire_P1.resize(E_size, 3);
  wire_P2.resize(E_size, 3);
  for(int e_idx = 0; e_idx < (E_size*3); e_idx+=3){
    wire_P1.row(e_idx/3) = Eigen::Vector3d(P1triplets[e_idx], P1triplets[e_idx+1], P1triplets[e_idx+2]);
    wire_P2.row(e_idx/3) = Eigen::Vector3d(P2triplets[e_idx], P2triplets[e_idx+1], P2triplets[e_idx+2]);
  }

  // write out to final format
  newVCenters.resize(newVerticesVector.size(), 3);
  for(int i = 0; i < newVerticesVector.size(); i++){
    newVCenters.row(i) = _V.row(newVerticesVector[i]);
  }

  // final format
  newEdges.resize(E_size, 2);
  int e_idx = 0;
  for(int i = 0; i < seenEdgesRelnewV.rows(); i++){
    for(int j = i; j < seenEdgesRelnewV.rows(); j++){
      if(seenEdgesRelnewV(i,j) == 1){
        newEdges.row(e_idx) << i, j;
        e_idx++;
      }
    }
  }
}
