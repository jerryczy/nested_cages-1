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
  VA.resize(V1.rows(), 12);
  
  std::vector<std::vector<int>> adj;
  igl::adjacency_list(F, adj);
  
  for (int i = 0; i < V1.rows(); i++) {
    Eigen::MatrixXd X(3*adj[i].size(), 12);
    Eigen::VectorXd Y(3*adj[i].size());
    for (int j = 0; j < adj[i].size(); j++) {
      Eigen::RowVector3d x = V1.row(adj[i][j]);
      Eigen::RowVector3d y = V2.row(adj[i][j]);
      X.block(3*j, 0, 3, 12) << x, 0, 0, 0, 0, 0, 0, 1, 0, 0,
                                0, 0, 0, x, 0, 0, 0, 0, 1, 0,
                                0, 0, 0, 0, 0, 0, x, 0, 0, 1;
      Y.block(3*j, 0, 3, 1) << Eigen::Vector3d(y);
    }
    
    VA.row(i) = (X.transpose() * X).ldlt().solve(X.transpose() * Y);
  }
}