#include "normalSet.h"
#include <iostream>

int NormalSet::current_id = -1;

NormalSet::NormalSet(bool is_painted){
  avg_normal = Eigen::Vector3d::Zero();
  painted = is_painted;
  id = current_id;
  current_id++;
  area = 0;
  perimeter = 0;
}

NormalSet::NormalSet(int face_idx, Eigen::Vector3d normal, double cur_area) {
  face_set.insert(face_idx);
  avg_normal = normal;
  id = current_id;
  current_id++;
  painted = false;
  area = cur_area;
  perimeter = 0;
}

void NormalSet::addToSet(int face_idx, Eigen::Vector3d normal, double cur_area) {
  // Update avg_normal
  updateAvgNormal(normal);
  // Add new face to set
  face_set.insert(face_idx);
  area += cur_area;
}

void NormalSet::erase(int face_idx, Eigen::Vector3d normal, double cur_area) {
	// Get current face_set size
	int num_faces = face_set.size();
	// Update avg_normal
	avg_normal = (avg_normal*num_faces - normal) / (num_faces - 1);
	// Remove last inserted face
	face_set.erase(face_idx);
	area -= cur_area;
}

void NormalSet::updateAvgNormal(Eigen::Vector3d new_normal){
  int num_faces = face_set.size();
  avg_normal = ((avg_normal * num_faces) + new_normal) / (num_faces + 1);
}

void NormalSet::clearSet(){
  face_set.clear();
  avg_normal = Eigen::Vector3d::Zero();
  current_id = 0;
}

void NormalSet::computeBoundary(Eigen::MatrixXi &F, Eigen::MatrixXd &V) {
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
  for (int j = 0; j < bnd.size(); j++) {
	  perimeter += ((V.row(bnd(j)) - V.row(bnd((j + 1) % bnd.size()))).norm());
  }
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
