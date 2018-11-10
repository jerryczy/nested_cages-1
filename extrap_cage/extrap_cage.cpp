#include "extrap_cage.h"
#include "cage_face_interpolate.h"
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
  cage_face_interpolate(V1, F, V1_C, F_C, W);
  
  Eigen::MatrixXd FA, FT;
  recover_affine_transformations(V1, V2, F, FA, FT);
  
  for (int f = 0; f < F.rows(); f++) {
    Eigen::MatrixXd A(3, 3);
    A << FA.row(f);
    Eigen::RowVector3d t = FT.row(f);
    
    for (int vc = 0; vc < V1_C.rows(); vc++) {
      V2_C.row(vc) += W(vc, f) * (V1_C.row(vc) * A.transpose() + t);
    }
  }
}