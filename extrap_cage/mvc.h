#include <Eigen/Core>

// Implementation of the algorithm described in the paper "Mean Value Coordinates
// for Closed Triangular Meshes" by Tao Ju, Scott Schaefer and Joe Warren.
// Given a mesh (V, F) and a list of points X, compute the mean value coordinates
// of each point X with respect to the mesh.
//
// Inputs:
//   V: #V by 3 matrix of vertices
//   F: #F by 3 matrix of faces
//   X: n by 3 matrix of query points
//   W: n by #V matrix of weights

void mvc(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::MatrixXd & X,
  Eigen::MatrixXd & W);