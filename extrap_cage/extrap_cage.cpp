#include "extrap_cage.h"
#include "mvc.h"

void extrap_cage(
  const Eigen::MatrixXd & V1,
  const Eigen::MatrixXd & V2,
  const Eigen::MatrixXi & F,
  const Eigen::MatrixXd & V1_C,
  const Eigen::MatrixXi & F_C,
  Eigen::MatrixXd & V2_C)
{
  Eigen::MatrixXd W(V1_C.rows(), 3);
  mvc(V1, F, V1_C, W);
  V2_C = W * V2;
}