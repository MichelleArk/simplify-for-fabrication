#include "RAG.h"

RAG::RAG(std::vector<NormalSet> normal_sets, Eigen::MatrixXd &V, Eigen::MatrixXi &F) {
  regions = normal_sets;
  _V = V; _F = F;
  // Compute double area and total area
	igl::doublearea(V, F, _doubA);
  _total_area = _doubA.sum() / 2.0;
  int num_regions = 80;

  // edge_costs initialization
  edge_costs.resize(normal_sets.size(), normal_sets.size());
  edge_costs.setZero();
  typedef Eigen::Triplet< double > Triplet;
  std::vector< Triplet > triplets;
	// Compute edge cost if find shared boundaries
	std::set<int> foundSharedVertices;
	for (std::vector<NormalSet>::iterator iter1 = normal_sets.begin(); iter1 != normal_sets.end(); iter1++) {
		NormalSet set1 = *iter1;
    if(!set1.painted){ // don't merge painted sets
  		for (std::vector<NormalSet>::iterator iter2 = normal_sets.begin(); iter2 != normal_sets.end(); iter2++) {
  			NormalSet set2 = *iter2;
        if(!set2.painted){
    			if (set1.id != set2.id) {
    				// Find shared boundary
    				std::vector<int> endpoints;
    				double shared_bnd_length;
    				if (set1.id < set2.id) { // To avoid adding ij and ji twice
    					if (sharedBoundary(set1.bnd, set2.bnd, set1.painted, set2.painted, endpoints, foundSharedVertices, V, shared_bnd_length)) {
                double area_weight = fmin(set1.area, set2.area) / _total_area;
                if(area_weight < ((F.rows() / (double) num_regions) / 4.0) / F.rows()){ // one region is less than 4 faces
                  area_weight = 0.0000001;
                }else if (fmax(set1.area, set2.area) / _total_area > ((F.rows() / (double) num_regions) * 2.0) / F.rows()){
                  area_weight = 20;
                }else{
                  area_weight = 10;
                }
    						double perimeter_weight = fmin(set1.perimeter, set2.perimeter) / shared_bnd_length;
    						set1.avg_normal.normalize();
    						set2.avg_normal.normalize();
    						double normal_weight = (1.0 - set1.avg_normal.dot(set2.avg_normal));
                double edge_cost = area_weight * (normal_weight*normal_weight) * perimeter_weight;

                // Add mappings to edge_to_f
                triplets.push_back(Triplet(set1.id, set2.id, edge_cost));
    					}
    				}
    			}
        }
  		}
    }
	}
  edge_costs.setFromTriplets(triplets.begin(), triplets.end());
}

void RAG::MergeMinCostRegions(int num_to_remove, int num_regions){
  std::vector<NormalSet>::iterator set_i, set_j;
  for(int regions_removed = 0; regions_removed < num_to_remove; regions_removed++){
    // Find regions to merge
    GetMinCostRegions(set_i, set_j);

    // Get set sizes before any updates
    int seti_size = (*set_i).face_set.size();
    int setj_size = (*set_j).face_set.size();

    int seti_id = (*set_i).id;
    int setj_id = (*set_j).id;
    // Add set_j's faces to set_i's faces
    for (std::set<int>::iterator set_iter = (*set_j).face_set.begin(); set_iter != (*set_j).face_set.end(); set_iter++) {
      (*set_i).face_set.insert(*set_iter);
    }

    // Update set_i's normal
    (*set_i).avg_normal = ((*set_i).avg_normal*seti_size + (*set_j).avg_normal*setj_size) / (seti_size + setj_size);

    // Update set_i's area
    (*set_i).area += (*set_j).area;

    // Remove set_j
    regions.erase(set_j);
    _removed_regions.insert(setj_id);

    // Update all boundaries
    for (std::vector<NormalSet>::iterator set = regions.begin(); set != regions.end(); set++) {
      (*set).computeBoundary(_F, _V);
    }

    // Update cost at all neighbors of set_i and set_j with set_i
    //UpdateEdgeCosts(seti_id, setj_id);
    UpdateEdgeCosts(num_regions);
  }
}

void RAG::GetMinCostRegions(std::vector<NormalSet>::iterator & set_i, std::vector<NormalSet>::iterator & set_j){
  int set1_id, set2_id;
  double min_cost = -1;
  for (int k=0; k < edge_costs.outerSize(); ++k){
    for (Eigen::SparseMatrix<double>::InnerIterator it(edge_costs,k); it; ++it){
      // row and col have not been removed
      if(_removed_regions.find(it.row()) == _removed_regions.end() && _removed_regions.find(it.col()) == _removed_regions.end()){
        double edge_cost = edge_costs.coeffRef(it.row(), it.col());
        if(min_cost == -1 || (edge_cost < min_cost && edge_cost != -1)){
          min_cost = edge_cost;
          set1_id = it.row();
          set2_id = it.col();
        }
      }
    }
  }

  // get pointers to normal sets with indices set1_id and set2_id
  for (std::vector<NormalSet>::iterator iter = regions.begin(); iter != regions.end(); iter++) {
		NormalSet set = *iter;
    if(set.id == set1_id){
      set_i = iter;
    }
    if(set.id == set2_id){
      set_j = iter;
    }
  }
  std::cout << "merging " << set1_id << " with " << set2_id << " with cost " << min_cost << std::endl;
}

void RAG::UpdateEdgeCosts(int num_regions){
  edge_costs.setZero();

  typedef Eigen::Triplet< double > Triplet;
  std::vector< Triplet > triplets;
	// Compute edge cost if find shared boundaries
	std::set<int> foundSharedVertices;
	for (std::vector<NormalSet>::iterator iter1 = regions.begin(); iter1 != regions.end(); iter1++) {
		NormalSet set1 = *iter1;
    if(!set1.painted){ // don't merge painted sets
  		for (std::vector<NormalSet>::iterator iter2 = regions.begin(); iter2 != regions.end(); iter2++) {
  			NormalSet set2 = *iter2;
        if(!set2.painted){
    			if (set1.id != set2.id) {
    				// Find shared boundary
    				std::vector<int> endpoints;
    				double shared_bnd_length;
    				if (set1.id < set2.id) { // To avoid adding ij and ji twice
    					if (sharedBoundary(set1.bnd, set2.bnd, set1.painted, set2.painted, endpoints, foundSharedVertices, _V, shared_bnd_length)) {
                double area_weight = fmin(set1.area, set2.area) / _total_area;
                if(area_weight < ((_F.rows() / (double) num_regions) / 4.0) / _F.rows()){ // one region is less than 4 faces
                  area_weight = 0.0000001;
                }else if (fmax(set1.area, set2.area) / _total_area > ((_F.rows() / (double) num_regions) * 2.0) / _F.rows()){
                  area_weight = 20;
                }else{
                  area_weight = 10;
                }
    						double perimeter_weight = fmin(set1.perimeter, set2.perimeter) / shared_bnd_length;
    						set1.avg_normal.normalize();
    						set2.avg_normal.normalize();
    						double normal_weight = (1.0 - set1.avg_normal.dot(set2.avg_normal));
                double edge_cost = area_weight * (normal_weight*normal_weight) * perimeter_weight;

                // Add mappings to edge_cost
                triplets.push_back(Triplet(set1.id, set2.id, edge_cost));
    					}
    				}
    			}
        }
  		}
    }
	}
  edge_costs.setFromTriplets(triplets.begin(), triplets.end());
}
void RAG::UpdateEdgeCosts(int set_remaining, int set_removed){
  // Handle by case
  for (int k=0; k < edge_costs.outerSize(); ++k){
    for (Eigen::SparseMatrix<double>::InnerIterator it(edge_costs,k); it; ++it){
      int set1_id = it.row();
      int set2_id = it.col();

      if(set1_id == set_remaining && set2_id != set_removed){ // set2 is neighbor of set that was kept
        if(_removed_regions.find(set2_id) == _removed_regions.end()){
          double new_edge_cost = GetEdgeCostBetweenRegions(set2_id, set_remaining);
          edge_costs.coeffRef(set1_id, set2_id) = new_edge_cost;
        }
      }
      if(set2_id == set_remaining && set1_id != set_removed){ // set1 is neighbor of set that was kept
        if(_removed_regions.find(set1_id) == _removed_regions.end()){
          double new_edge_cost = GetEdgeCostBetweenRegions(set1_id, set_remaining);
          edge_costs.coeffRef(set1_id, set2_id) = new_edge_cost;
        }
      }
      if(set1_id == set_removed){ // set2 is neighbor of set that was removed
        if(set2_id != set_remaining && _removed_regions.find(set2_id) == _removed_regions.end()){
          double new_edge_cost = GetEdgeCostBetweenRegions(set2_id, set_remaining);
          edge_costs.coeffRef(set2_id, set_remaining) = new_edge_cost;
        }
      }
      if(set2_id == set_removed){ // set1 is neighbor of set that was removed
        if(set1_id != set_remaining && _removed_regions.find(set1_id) == _removed_regions.end()){
          double new_edge_cost = GetEdgeCostBetweenRegions(set1_id, set_remaining);
          edge_costs.coeffRef(set1_id, set_remaining) = new_edge_cost;
        }
      }
    }
  }
}

double RAG::GetEdgeCostBetweenRegions(int region1_id, int region2_id){
  // get pointers to normal sets with indices set1_id and set2_id
  NormalSet set1, set2;
  for (std::vector<NormalSet>::iterator iter = regions.begin(); iter != regions.end(); iter++) {
    NormalSet set = *iter;
    if(set.id == region1_id){
      set1 = set;
    }
    if(set.id == region2_id){
      set2 = set;
    }
  }

  double edge_cost = -1; // stupid
  // Find shared boundary
  std::vector<int> endpoints;
  double shared_bnd_length;
  std::set<int> foundSharedVertices;
  if (sharedBoundary(set1.bnd, set2.bnd, set1.painted, set2.painted, endpoints, foundSharedVertices, _V, shared_bnd_length)) {
    double area_weight = fmin(set1.area, set2.area) / _total_area;
    if(area_weight < ((_F.rows() / (double) regions.size()) / 4.0) / _F.rows()){ // one region is less than 4 faces
      area_weight = 0.0000001;
    }else if (fmax(set1.area, set2.area) / _total_area > ((_F.rows() / (double) regions.size()) * 2.0) / _F.rows()){
      area_weight = 20;
    }else{
      area_weight = 10;
    }
    double perimeter_weight = fmin(set1.perimeter, set2.perimeter) / shared_bnd_length;
    set1.avg_normal.normalize();
    set2.avg_normal.normalize();
    double normal_weight = (1.0 - set1.avg_normal.dot(set2.avg_normal));
    edge_cost = area_weight * (normal_weight*normal_weight) * perimeter_weight;
  }
  return edge_cost;
}
