#include <Eigen>

// Given two meshes having a one-to-one vertex and face correspondence, return a
// list of matrices and transformations such that the ith such pair is an affine
// transformation that takes the ith face of the first mesh to the ith face of the
// second.
//
// Inputs:
//   V1  #V by 3 matrix of vertex positions for the first mesh
//   V2  #V by 3 matrix of vertex positions for the second mesh
//   F   #F by 3 matrix of face positions common to both meshes
// Outputs:
//   FA   #F by 9 matrix where the ith row of A is the (col-major) flattened 3x3 
//       matrix corresponding to the linear transformation for the ith face
//   FT   #F by 3 matrix where the ith row of T is the translation vector for
//       the ith face

void recover_face_transformations(
  const Eigen::MatrixXd & V1,
  const Eigen::MatrixXd & V2,
  const Eigen::MatrixXd & F,
  Eigen::MatrixXd & FA,
  Eigen::MatrixXd & FT);