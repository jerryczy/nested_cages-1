#include "extrap_cage.h"
#include "mvc.h"
#include <Eigen/QR>

void extrap_cage(
  const Eigen::MatrixXd & V1,
  const Eigen::MatrixXd & V2,
  const Eigen::MatrixXi & F,
  const Eigen::MatrixXd & V1_C,
  const Eigen::MatrixXi & F_C,
  Eigen::MatrixXd & V2_C)
{
  Eigen::MatrixXd W;
  mvc(V1_C, F_C, V1, W);
  assert((V1 - W * V1_C).maxCoeff() < 0.01);
  V2_C = (W.transpose() * W).ldlt().solve(W.transpose() * V2);
}