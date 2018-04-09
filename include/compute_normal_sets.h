#ifndef COMPUTE_NORMAL_SETS_H
#define COMPUTE_NORMAL_SETS_H

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
#include <set>
#include <map>
#include <queue>
#include <igl/per_face_normals.h>
#include <igl/read_triangle_mesh.h>
#include "normalSet.h"

void compute_normal_sets(
  const Eigen::MatrixXi &F,
  const Eigen::MatrixXd &V,
  std::vector<NormalSet>& all_normal_sets,
  std::set<int>& visited);

std::vector<int> get_neighbours(
  const Eigen::MatrixXi &F,
  int f_idx,
  std::map<std::string,int> edge_to_f);

bool similar_normals(
  Eigen::Vector3d n1,
  Eigen::Vector3d n2);

std::map<std::string, int> preprocess_edge_to_face(const Eigen::MatrixXi &F);

bool sharedBoundary(
  Eigen::VectorXi bnd1,
  Eigen::VectorXi bnd2,
  bool set1_painted,
  bool set2_painted,
  std::vector<int> &endpoints,
  std::set<int> &foundSharedVertices,
  Eigen::MatrixXd &V,
  double &shared_bnd_length);

// These two functions aren't currently used. But they will be helpful in implementing the future work
void computeAngleDeficitPerVertex(
  std::vector<int> newVertices,
  Eigen::MatrixXd V,
  std::vector<NormalSet> normal_sets,
  Eigen::VectorXd &C,
  int& min_cost_vid,
  double& min_cost);

void removeVertex(
  std::vector<int> &newVertices,
  std::vector<NormalSet> &normal_sets,
  int v_idx);

#endif
