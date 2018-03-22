#ifndef COMPUTE_NORMAL_SETS_H
#define COMPUTE_NORMAL_SETS_H

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
#include <set>
#include <map>
#include <queue>
#include <igl/per_face_normals.h>
#include "normalSet.h"
#include <igl/boundary_loop.h>

// Compute the Euler Characteristic of a given triangle mesh.
//
// Inputs:
//   F  #F by 3 list of triangle indices into some vertex list V
// Returns Euler Characteristic as an integer
std::vector<NormalSet> compute_normal_sets( const Eigen::MatrixXi &F, const Eigen::MatrixXd &V);

std::vector<int> get_neighbours(const Eigen::MatrixXi &F, int f_idx, std::map<std::string,int> edge_to_f);

bool similar_normals(Eigen::Vector3d n1, Eigen::Vector3d n2);

std::map<std::string, int> preprocess_edge_to_face(const Eigen::MatrixXi &F);

bool sharedBoundary(Eigen::VectorXi bnd1, Eigen::VectorXi bnd2, std::vector<int> &endpoints, std::set<int> &foundSharedVertices);

void straightenEdges(Eigen::MatrixXd &V, Eigen::MatrixXi &F, std::vector<NormalSet> &normal_sets, Eigen::MatrixXd &newV, Eigen::MatrixXi &newF);
#endif
