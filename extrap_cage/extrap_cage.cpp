#include "extrap_cage.h"
#include "cage_vertex_interpolate.h"
#include "recover_affine_transformations.h"

void extrap_cage(
  const Eigen::MatrixXd & V1,
  const Eigen::MatrixXd & V2,
  const Eigen::MatrixXi & F,
  const Eigen::MatrixXd & V1_C,
  const Eigen::MatrixXi & F_C,
  Eigen::MatrixXd & V2_C)
{
  V2_C = Eigen::MatrixXd::Zero(V1_C.rows(), 3);
  
  Eigen::MatrixXd W;
  cage_vertex_interpolate(V1, F, V1_C, F_C, W);
  
  Eigen::MatrixXd VA;
  recover_affine_transformations(V1, V2, F, VA);
  
  for (int vi = 0; vi < V1_C.rows(); vi++) {
    Eigen::MatrixXd A(3, 3);
    Eigen::Vector3d t;
    A << VA.block(vi, 0, 1, 9);
    t << VA.block(vi, 9, 1, 3);
    for (int vci = 0; vci < V1_C.rows(); vci++) {
      V2_C.row(vci) += W(vci, vi) * A * (V1_C.row(vci) - V1.row(vi)) + t + V1.row(vi);
    }
  }
}