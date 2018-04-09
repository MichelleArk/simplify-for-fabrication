#include "RAG.h"

RAG::RAG(std::vector<NormalSet> normal_sets, Eigen::MatrixXd &V, Eigen::MatrixXi &F) {
  regions = normal_sets;
  _V = V; _F = F;
  // Compute double area and total area
	igl::doublearea(V, F, _doubA);
  _total_area = _doubA.sum() / 2.0;
  // initialize edge_costs
  edge_costs.resize(normal_sets.size(), normal_sets.size());
  edge_costs.setZero();
  UpdateEdgeCosts(80);
}

void RAG::MergeMinCostRegions(int num_to_remove, int num_regions){
  std::vector<NormalSet>::iterator set_i, set_j;
  for(int regions_removed = 0; regions_removed < num_to_remove; regions_removed++){
    // Find regions to merge
    set_i = set_to_merge1;
    set_j = set_to_merge2;
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

    // Update all boundaries
    for (std::vector<NormalSet>::iterator set = regions.begin(); set != regions.end(); set++) {
      (*set).computeBoundary(_F, _V);
    }

    // Remove set_j
    regions.erase(set_j);
    _removed_regions.insert(setj_id);

    UpdateEdgeCosts(num_regions);
  }
}

void RAG::UpdateEdgeCosts(int num_regions){
  edge_costs.setZero();

  double min_edge_cost = -1;

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
                double area_threshold =  ((_F.rows() / (double) num_regions) / 4.0) / _F.rows();
                if(area_weight < area_threshold){ // one region is less than 4 faces
                  area_weight = 0.0000001;
                }else{
                  area_weight = 1;
                }
    						double perimeter_weight = fmin(set1.perimeter, set2.perimeter) / shared_bnd_length;
    						set1.avg_normal.normalize();
    						set2.avg_normal.normalize();
    						double normal_weight = (1.0 - set1.avg_normal.dot(set2.avg_normal));
                double edge_cost = area_weight * (normal_weight*normal_weight) * perimeter_weight;
                if(min_edge_cost == -1 || edge_cost < min_edge_cost){
                  set_to_merge1 = iter1;
                  set_to_merge2 = iter2;
                  min_edge_cost = edge_cost;
                }
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
