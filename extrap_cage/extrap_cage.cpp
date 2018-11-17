#include "extrap_cage.h"
#include <igl/point_mesh_squared_distance.h>
#include <igl/barycentric_coordinates.h>

void extrap_cage(
  const Eigen::MatrixXd & V1,
  const Eigen::MatrixXd & V2,
  const Eigen::MatrixXi & F,
  const Eigen::MatrixXd & V1_C,
  const Eigen::MatrixXi & F_C,
  Eigen::MatrixXd & V2_C)
{
  V2_C = V1_C.replicate(1, 1);
  Eigen::MatrixXd V2_C_inc = Eigen::MatrixXd::Zero(V2_C.rows(), 3);
  
  Eigen::VectorXd sqrD;
  Eigen::VectorXi I;
  Eigen::MatrixXd C;
  igl::point_mesh_squared_distance(V1, V1_C, F_C, sqrD, I, C);
  
  Eigen::VectorXd total_vertex_weights = Eigen::VectorXd::Zero(V2_C.rows());
  
  for (int vi = 0; vi < V1.rows(); vi++) {
    Eigen::RowVector3i coarse_face = F_C.row(I[vi]);
    Eigen::RowVector3d coarse_point = C.row(vi);
    
    Eigen::RowVector3d a = V1_C.row(coarse_face[0]);
    Eigen::RowVector3d b = V1_C.row(coarse_face[1]);
    Eigen::RowVector3d c = V1_C.row(coarse_face[2]);
    Eigen::RowVector3d l;
    igl::barycentric_coordinates(coarse_point, a, b, c, l);
    
    V2_C_inc.row(coarse_face[0]) += (1 - l[0]) * (V2.row(vi) - V1.row(vi));
    V2_C_inc.row(coarse_face[1]) += (1 - l[1]) * (V2.row(vi) - V1.row(vi));
    V2_C_inc.row(coarse_face[2]) += (1 - l[2]) * (V2.row(vi) - V1.row(vi));
    
    total_vertex_weights[coarse_face[0]] += 1 - l[0];
    total_vertex_weights[coarse_face[1]] += 1 - l[1];
    total_vertex_weights[coarse_face[2]] += 1 - l[2];
  }
  
  for (int vci = 0; vci < V2_C.rows(); vci++) {
    V2_C_inc.row(vci) /= total_vertex_weights[vci];
  }
  
  V2_C += V2_C_inc;
}