#pragma once

#ifndef RAG_H
#define RAG_H

#include <Eigen/Core>
#include "normalSet.h"
#include "compute_normal_sets.h"
#include <igl/doublearea.h>
#include <igl/gaussian_curvature.h>

class RAG {
public:
  // constructor
  RAG(){};
	RAG(std::vector<NormalSet> normal_sets, Eigen::MatrixXd &V, Eigen::MatrixXi &F);

  void MergeMinCostRegions(int num_to_remove, int num_regions);
  void UpdateEdgeCosts(int num_regions);
  void UpdateEdgeCostsSmartly(int set_remaining, int set_removed, std::set<int> neighbors_to_update);
  void GetMinCostRegions(std::vector<NormalSet>::iterator & set_i, std::vector<NormalSet>::iterator & set_j);
  double GetEdgeCostBetweenRegions(int region1_id, int region2_id);


  std::vector<NormalSet> regions;
  Eigen::SparseMatrix<double> edge_costs;

private:
  Eigen::VectorXd _doubA;
  double _total_area;
  Eigen::MatrixXd _V;
  Eigen::MatrixXi _F;

  std::set<int> _removed_regions;
};

#endif
