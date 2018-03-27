#pragma once

#include <Eigen/Core>
#include <set>
#include <igl/boundary_loop.h>

class NormalSet {
  protected:
    static int current_id;
  public:
    NormalSet();
    NormalSet(int face_idx, Eigen::Vector3d normal);
    void addToSet(int face_idx, Eigen::Vector3d normal);
    void clearSet();
  	void computeBoundary(Eigen::MatrixXi &F);
  	void simplifyBoundary(std::set<int> new_bnd);

	int id;
	std::set<int> face_set;
    Eigen::Vector3d avg_normal;
	Eigen::VectorXi bnd;
	Eigen::VectorXi simplified_bnd;
};
