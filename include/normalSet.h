#pragma once

#include <Eigen/Core>
#include <set>

class NormalSet {
  public:
    NormalSet();
    NormalSet(int face_idx, Eigen::Vector3d normal, int id);
    void addToSet(int face_idx, Eigen::Vector3d normal);
	void addBoundary(Eigen::VectorXi boundary);
	void simplifyBoundary(std::set<int> new_bnd);

	int id;
	std::set<int> face_set;
    Eigen::Vector3d avg_normal;
	Eigen::VectorXi bnd;
	Eigen::VectorXi simplified_bnd;
};
