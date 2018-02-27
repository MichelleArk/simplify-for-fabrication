#include "normalSet.h"

NormalSet::NormalSet(){
  avg_normal = Eigen::Vector3d::Zero();
}

NormalSet::NormalSet(int face_idx, Eigen::Vector3d normal) {
  face_set.insert(face_idx);
  avg_normal = normal;
}

void NormalSet::addToSet(int face_idx, Eigen::Vector3d normal) {
  // Update avg_normal
  int num_faces = face_set.size();
  avg_normal = ((avg_normal * num_faces) + normal) / (num_faces + 1);
  // Add new face to set
  face_set.insert(face_idx);
}
