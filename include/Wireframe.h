#pragma once

#ifndef WIREFRAME_H
#define WIREFRAME_H

#include <Eigen/Core>
#include "normalSet.h"
#include "compute_normal_sets.h"

class Wireframe {
public:
  // constructor
  Wireframe(){};
	Wireframe(Eigen::MatrixXd V, Eigen::MatrixXi F);

  void Update(std::vector<NormalSet> new_regions);
  void straightenEdges();
  void createApproxSpheres();
  void connectApproxSpheres();

  std::vector<NormalSet> regions;
  Eigen::MatrixXd wire_V; // for viewer.set_mesh()
  Eigen::MatrixXi wire_F;
  Eigen::MatrixXd wire_P1; // for viewer.set_edges()
  Eigen::MatrixXd wire_P2;

  //indices into _V of remaining newVertices
  std::vector<int> newVerticesVector;
  // Final output for writing
  Eigen::MatrixXd newVCenters;
  Eigen::MatrixXi newEdges;

private:
  Eigen::MatrixXd _V; //original mesh
  Eigen::MatrixXi _F; // original mesh
};

#endif
