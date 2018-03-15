#pragma once

#include <Eigen/Core>
#include <set>

class NormalSet {
  public:
    NormalSet();
    NormalSet(int face_idx, Eigen::Vector3d normal, int id);
    void addToSet(int face_idx, Eigen::Vector3d normal);
	void NormalSet::addBoundary(Eigen::VectorXi boundary);

	int id;
	std::set<int> face_set;
    Eigen::Vector3d avg_normal;
	Eigen::VectorXi bnd;
	std::set<int> bndset;
};
