#include "normalSet.h"
#include <iostream>

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

void NormalSet::clearSet(){
  face_set.clear();
  avg_normal = Eigen::Vector3d::Zero();
}

void NormalSet::addBoundary(Eigen::VectorXi boundary) {
	bnd = boundary;
}

void NormalSet::simplifyBoundary(std::set<int> new_bnd) {
	simplified_bnd.resize(new_bnd.size());
	int idx = 0;
	for (int i = 0; i < bnd.size(); i++) {
		std::set<int>::iterator new_iter;
		for (new_iter = new_bnd.begin(); new_iter != new_bnd.end(); ++new_iter) {
			if (bnd(i) == *new_iter) {
				//std::cout << "inside simplify" << std::endl;
				simplified_bnd(idx) = bnd(i);
				//std::cout << simplified_bnd(idx) << std::endl;
				idx++;
			}
		}
	}
	//std::cout << "what you want to check" << std::endl;
	//std::cout << simplified_bnd << std::endl;
}
