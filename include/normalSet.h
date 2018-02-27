#pragma once

#include <Eigen/Core>
#include <set>

class NormalSet {
  public:
    NormalSet();
    NormalSet(int face_idx, Eigen::Vector3d normal);
    void addToSet(int face_idx, Eigen::Vector3d normal);

    std::set<int> face_set;
    Eigen::Vector3d avg_normal;
};
