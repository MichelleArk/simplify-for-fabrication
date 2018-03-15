#include "compute_normal_sets.h"
#include <iostream>

typedef Eigen::Triplet<int> T;

std::vector<NormalSet> compute_normal_sets(const Eigen::MatrixXi &F, const Eigen::MatrixXd &V)
{
  // initialize normal sets
  std::vector<NormalSet> all_normal_sets;

  std::map<std::string, int> edge_to_face = preprocess_edge_to_face(F);

  Eigen::MatrixXd N;
  igl::per_face_normals(V, F, Eigen::Vector3d(1,1,1).normalized(), N);

  std::set<int> visited;
  int unseen_face_idx = 0;

  std::cout << "# Total Faces: " << F.rows() << std::endl;
  // BFS on entire F graph
  while(visited.size() != F.rows()){
    std::cout << "# Faces visited: " << visited.size() << std::endl;
    // Get a 'random' face_idx from to_visit
    // while unseen_face_idx is in visited, increment it
    while(visited.find(unseen_face_idx) != visited.end()){
      unseen_face_idx += 1;
    }
    int face_idx = unseen_face_idx;
    // Initialize set with unseen face in it
    NormalSet normalSet(face_idx, N.row(face_idx));

    // Initialize Q with random unseen face on it
    std::queue<int> Q;
    Q.push(face_idx);

    // BFS on a single 'connected' componenet
    std::set<int> already_seen;
    // push onto visited for cycle detection / avoidance
    visited.insert(face_idx);

    while(!Q.empty()){
      // Get current face
      int curr_face = Q.front();
      Q.pop();

      // Iterate over current face's neighbours
      std::vector<int> neighbours = get_neighbours(F, curr_face, edge_to_face);
      for(std::vector<int>::iterator itr = neighbours.begin(); itr != neighbours.end(); itr++){
        int neighbour = *itr;
        Eigen::Vector3d neighbour_normal = N.row(neighbour);
        // Check that neighbour has similar normal and has not been on the queue already
        if( similar_normals(neighbour_normal, normalSet.avg_normal) && !(visited.find(neighbour) != visited.end()) ){
          Q.push(neighbour);
          normalSet.addToSet(neighbour, neighbour_normal);

          //push onto visited for cycle detection / avoidance
          visited.insert(neighbour);
        }
      }
    }
    all_normal_sets.push_back(normalSet);
  }
  return all_normal_sets;
}

bool similar_normals(Eigen::Vector3d n1, Eigen::Vector3d n2){
  double threshold = 0.866; // 30 degrees
  //double threshold = 0.8;

  // TODO: sanity check this works ok
  n1.normalize();
  n2.normalize();
  return n1.dot(n2) > threshold;
}

std::vector<int> get_neighbours(
  const Eigen::MatrixXi &F,
  int f_idx,
  std::map<std::string,int> edge_to_f)
{
  std::vector<int> neighbours;

  int v0 = F(f_idx, 0); int v1 = F(f_idx, 1); int v2 = F(f_idx, 2);

  std::string key10 = std::to_string(v1) + "," + std::to_string(v0);
  std::string key21 = std::to_string(v2) + "," + std::to_string(v1);
  std::string key02 = std::to_string(v0) + "," + std::to_string(v2);

  int n1 = edge_to_f[key10];
  int n2 = edge_to_f[key21];
  int n3 = edge_to_f[key02];

  // TODO: double-check that this is the right way to check if n_i was found
  if(n1 != 0){
    neighbours.push_back(n1);
  }
  if(n2 != 0){
    neighbours.push_back(n2);
  }
  if(n3 != 0){
    neighbours.push_back(n3);
  }
  return neighbours;
}

std::map<std::string, int> preprocess_edge_to_face( const Eigen::MatrixXi &F ){
  std::map<std::string, int> edge_to_f;
  int num_faces = F.rows();
  for(int f_idx= 0; f_idx < num_faces; f_idx++){
    int v0 = F(f_idx, 0); int v1 = F(f_idx, 1); int v2 = F(f_idx, 2);
    // Create string keys
    std::string key01 = std::to_string(v0) + "," + std::to_string(v1);
    std::string key12 = std::to_string(v1) + "," + std::to_string(v2);
    std::string key20 = std::to_string(v2) + "," + std::to_string(v0);
    // Add mappings to edge_to_f
    edge_to_f.insert(make_pair(key01,f_idx));
    edge_to_f.insert(make_pair(key12,f_idx));
    edge_to_f.insert(make_pair(key20,f_idx));
  }
  return edge_to_f;
}
