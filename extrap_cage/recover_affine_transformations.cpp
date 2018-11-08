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
    
    // for now do something simple: align the centroid of both triangles with
    // origin and find the linear transformation taking one centered triangle
    // to the other; TODO: is this a good idea?
    
    Eigen::Vector3d c1 = (1/3.0) * (T1.row(0) + T1.row(1) + T1.row(2));
    Eigen::Vector3d c2 = (1/3.0) * (T2.row(0) + T2.row(1) + T2.row(2));
    
    // find the matrix A such that A*(T1-c1) = (T2-c2)
    FA.row(i) = ((T2 - c2) * (T1 - c1).inverse()).resize(1, 9);
    FT.row(i) = c2 - c1;
  }
}