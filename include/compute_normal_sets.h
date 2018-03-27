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
#include <igl/read_triangle_mesh.h>

void compute_normal_sets( const Eigen::MatrixXi &F, const Eigen::MatrixXd &V, std::vector<NormalSet>& all_normal_sets);

std::vector<int> get_neighbours(const Eigen::MatrixXi &F, int f_idx, std::map<std::string,int> edge_to_f);

bool similar_normals(Eigen::Vector3d n1, Eigen::Vector3d n2);

std::map<std::string, int> preprocess_edge_to_face(const Eigen::MatrixXi &F);

bool sharedBoundary(Eigen::VectorXi bnd1, Eigen::VectorXi bnd2, std::vector<int> &endpoints, std::set<int> &foundSharedVertices);

// TODO: move cost outta here
void straightenEdges(Eigen::MatrixXd &V, Eigen::MatrixXi &F, std::vector<NormalSet> &normal_sets, Eigen::MatrixXd &newV, Eigen::MatrixXi &newF, Eigen::MatrixXd &P1, Eigen::MatrixXd &P2, Eigen::VectorXd &Cost);

void createApproxSpheres(std::vector<int> icoCenters, Eigen::MatrixXd &V, Eigen::MatrixXd &newV, Eigen::MatrixXi &newF);

void connectApproxSpheres(std::vector<NormalSet> &normal_sets, Eigen::MatrixXd &V, Eigen::MatrixXd &P1, Eigen::MatrixXd &P2);

void computeRemovalCostPerVertex(std::vector<int> newVertices, Eigen::MatrixXd V, std::vector<NormalSet> normal_sets, Eigen::VectorXd &C);

void removeVertex(std::vector<int> &newVertices, std::vector<NormalSet> &normal_sets, int v_idx);
#endif
