#include "recover_affine_transformations.h"
#include <Eigen/LU>
#include <Eigen/QR>
#include <igl/adjacency_list.h>

void recover_affine_transformations(
  const Eigen::MatrixXd & V1,
  const Eigen::MatrixXd & V2,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & VA)
{
  VA.resize(F.rows(), 9);
  
  std::vector<std::vector<int>> adj;
  igl::adjacency_list(F, adj);
  
  for (int i = 0; i < V1.rows(); i++) {
    Eigen::MatrixXd X(3*adj[i].size(), 12);
    Eigen::Vector3d Y(3*adj[i].size());
    for (int j = 0; j < X.rows(); j++) {
      X.block(j, 4*(j % 3), 1, 4) << V1.row(adj[i][j / 3]) - V1.row(i), 1;
      Y[j] = V2(adj[i][j / 3], j % 3) - V2(i, j % 3);
    }
    
    VA.row(i) << (X.transpose() * X).ldlt().solve(X.transpose() * Y);
  }
}