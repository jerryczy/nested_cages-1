#include <Eigen/Core>

// Given two meshes having a one-to-one vertex and face correspondence, return a
// list of matrices and transformations such that the ith such pair is an affine
// transformation that takes the ith vertex of the first mesh to the ith vertex of
// the second.
//
// Inputs:
//   V1  #V by 3 matrix of vertex positions for the first mesh
//   V2  #V by 3 matrix of vertex positions for the second mesh
//   F   #F by 3 matrix of face positions common to both meshes
// Outputs:
//   VA   #V by 9 matrix where the ith row of A stores the 12 coefficients
//        corresponding to the affine transformation for the ith vertex

void recover_affine_transformations(
  const Eigen::MatrixXd & V1,
  const Eigen::MatrixXd & V2,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & VA);