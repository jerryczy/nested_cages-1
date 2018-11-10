#include "recover_affine_transformations.h"
#include <Eigen/LU>
#include <Eigen/QR>

void recover_affine_transformations(
  const Eigen::MatrixXd & V1,
  const Eigen::MatrixXd & V2,
  const Eigen::MatrixXi & F,
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
    
    Eigen::Vector3d c1 = (1/3.0) * (T1.col(0) + T1.col(1) + T1.col(2));
    Eigen::Vector3d c2 = (1/3.0) * (T2.col(0) + T2.col(1) + T2.col(2));
    
    // find the matrix A such that A*(T1-c1) = (T2-c2)
    //Eigen::MatrixXd A = (T2.colwise() - c2) * (T1.colwise() - c1).inverse();
    Eigen::MatrixXd A = (T1.colwise() - c1).colPivHouseholderQr().solve(T2.colwise() - c2);
    A.resize(1, 9);
    FA.row(i) = A;
    FT.row(i) = c2 - c1;


    // Affine transformation, need stack four points to solve linear transformation
    // use three vertices and center point
    /*
    |x'|   |a b c|   |x|   |j|     |x y z 0 0 0 0 0 0 1 0 0|   |a|   |x'|
    |y'| = |d e f| * |y| + |k| ==> |0 0 0 x y z 0 0 0 0 1 0| * |-| = |y'|
    |x'|   |g h i|   |z|   |l|     |0 0 0 0 0 0 x y z 0 0 1|   |l|   |z'|
    */
    /*
    // code, put in comment for now
    #include <Eigen/SVD>
    Eigen::MatrixXd T1 = Eigen::MatrixXd::Zeros(12, 12);
    Eigen::MatrixXd T2(3, 4);
    Eigen::Vector3d c1 = Eigen::Vector3d::Zeros(3);
    Eigen::Vector3d c2 = Eigen::Vector3d::Zeros(3);
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        T1.block(3*j + k, 3*k, 3, 1) = V1.row(F(i, j));
        T1(3*j + k, 9 + k) = 1;
      }
      T2.col(j) = V2.row(F(i, j)).transpose();
      c1 += (1 / 3.0) * V1.row(F(i, j));
      c2 += (1 / 3.0) * V2.row(F(i, j));
    }
    for (int k = 0; k < 3; k++) {
        T1.block(9 + k, 3*k, 3, 1) = c1;
        T1(9 + k, 9 + k) = 1;
      }
    T2.col(4) = c2.transpose();
    T2.resize(12, 1)
    
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(T1, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::VectorXd A = svd.solve(T2);
    FA.row(i) = A.head(9);
    FT.row(i) = A.tail(3);
    */
  }
}