#include "extrap_cage.h"

void extrap_cage(
  const Eigen::MatrixXd & V1,
  const Eigen::MatrixXd & V2,
  const Eigen::MatrixXd & F,
  const Eigen::MatrixXd & V1_C,
  const Eigen::MatrixXd & F_C,
  Eigen::MatrixXd & V2_C)
{
  V2_C.resize(V1_C.rows(), 3);
  
  Eigen::VectorXd & W;
  cage_face_interpolate(V1, F, V1_C, F_C, W);
  
  Eigen::MatrixXd FA, FT;
  recover_face_transformations(V1, F, V2, F2, FA, FT);
  
  for (int f = 0; f < F.rows(); f++) {  
    Eigen::Matrix3d A = FA.row(f).resize(3, 3);
    Eigen::Vector3d t = FT.row(f);
    
    for (int vc = 0; vc < V1_C.rows(); vc++) {
      V2_C.row(vc) += W(vc, f) * (V1_C.row(vc) * A.transpose() + t);
    }
  }
}