#pragma once

#include <Eigen/Core>
#include <set>
#include <igl/boundary_loop.h>

class NormalSet {
  protected:
    static int current_id;
  public:
    NormalSet(bool is_painted = false);
    NormalSet(int face_idx, Eigen::Vector3d normal, double cur_area);
    void addToSet(int face_idx, Eigen::Vector3d normal, double cur_area);
    void updateAvgNormal(Eigen::Vector3d normal);
    void clearSet();
	  void erase(int face_idx, Eigen::Vector3d normal, double cur_area);
	  void computeBoundary(const Eigen::MatrixXi &F, const Eigen::MatrixXd &V);
  	void simplifyBoundary(std::set<int> new_bnd);

  	int id;
  	bool painted;
  	std::set<int> face_set;
  	Eigen::Vector3d avg_normal;
  	Eigen::VectorXi bnd;
  	Eigen::VectorXi simplified_bnd;
  	double area;
  	double perimeter;

    // for RAG
    std::set<int> neighbors;
};
