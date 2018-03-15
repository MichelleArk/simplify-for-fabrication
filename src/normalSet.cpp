#include "normalSet.h"

NormalSet::NormalSet(){
  avg_normal = Eigen::Vector3d::Zero();
}

NormalSet::NormalSet(int face_idx, Eigen::Vector3d normal, int set_id) {
  face_set.insert(face_idx);
  avg_normal = normal;
  id = set_id;
}

void NormalSet::addToSet(int face_idx, Eigen::Vector3d normal) {
  // Update avg_normal
  int num_faces = face_set.size();
  avg_normal = ((avg_normal * num_faces) + normal) / (num_faces + 1);
  // Add new face to set
  face_set.insert(face_idx);
}

void NormalSet::addBoundary(Eigen::VectorXi boundary) {
	bnd = boundary;
	for (int i = 0; i < bnd.size(); i++) {
		bndset.insert(bnd(i));
	}
}