#include "recover_affine_transformations.h"

void recover_face_transformations(
  const Eigen::MatrixXd & V1,
  const Eigen::MatrixXd & V2,
  const Eigen::MatrixXd & F,
  Eigen::MatrixXd & FA,
  Eigen::MatrixXd & FT)
{
  FA.resize(F.rows(), 9);
  FT.resize(F.rows(), 3);
  
  for (int i = 0; i < F.rows(); i++) {
    Eigen::Matrix3d T1, T2;
    for (int j = 0; j < 3; j++) {
      T1.col(j) = V1.row(F(i, j)).transpose();
      T2.col(j) = V2.row(F(i, j)).transpose();
    }
    
    // find the "best" A and T s.t. T2 = A*T1 + T
    // ...
  }
  return tfms;
}