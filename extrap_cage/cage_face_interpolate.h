#include <Eigen>

// Given a mesh (V, F), and a coarse mesh (V_C, F_C) exterior to (V, F), return 
// a matrix W such W(i, j) measures the extent to which vertex i on the coarse mesh
// is "influenced" by face j on the *fine* mesh.
//
// Inputs:
//   V    #V by 3 matrix of vertices for fine mesh
//   F    #F by 3 matrix of faces for fine mesh
//   V_C  #V_C by 3 matrix of vertices for coarse mesh
//   F_C  #F_C by 3 matrix of faces for coarse mesh
// Outpts:
//   W    #V_C by F matrix of weights with each row summing to one.

void cage_face_interpolate(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXd & F,
  const Eigen::MatrixXd & V_C,
  const Eigen::MatrixXd & F_C,
  Eigen::MatrixXd & W);