#include "normalSet.h"
#include <iostream>

int NormalSet::current_id = -1;

NormalSet::NormalSet(bool is_painted){
  avg_normal = Eigen::Vector3d::Zero();
  painted = is_painted;
  id = current_id;
  current_id++;
}

NormalSet::NormalSet(int face_idx, Eigen::Vector3d normal) {
  face_set.insert(face_idx);
  avg_normal = normal;
  id = current_id;
  current_id++;
  painted = false;
}

void NormalSet::addToSet(int face_idx, Eigen::Vector3d normal) {
  // Update avg_normal
  //int num_faces = face_set.size();
  //avg_normal = ((avg_normal * num_faces) + normal) / (num_faces + 1);
  updateAvgNormal(normal);
  // Add new face to set
  face_set.insert(face_idx);
}

void NormalSet::updateAvgNormal(Eigen::Vector3d new_normal){
  int num_faces = face_set.size();
  avg_normal = ((avg_normal * num_faces) + new_normal) / (num_faces + 1);
}
void NormalSet::clearSet(){
  face_set.clear();
  avg_normal = Eigen::Vector3d::Zero();
}

void NormalSet::computeBoundary(Eigen::MatrixXi &F) {
  std::set<int> normal_set = face_set;
  Eigen::MatrixXi F_set(normal_set.size(), 3);
  std::set<int>::iterator face;
  int i = 0;
  for (face = normal_set.begin(); face != normal_set.end(); ++face) {
  	int face_idx = *face;
  	F_set.row(i) = F.row(face_idx);
  	i++;
  }
  Eigen::VectorXi boundary;
  igl::boundary_loop(F_set, boundary);
  bnd = boundary;
}

void NormalSet::simplifyBoundary(std::set<int> new_bnd) {
	simplified_bnd.resize(new_bnd.size());
	int idx = 0;
	for (int i = 0; i < bnd.size(); i++) {
		std::set<int>::iterator new_iter;
		for (new_iter = new_bnd.begin(); new_iter != new_bnd.end(); ++new_iter) {
			if (bnd(i) == *new_iter) {
				simplified_bnd(idx) = bnd(i);
				idx++;
			}
		}
	}
}
