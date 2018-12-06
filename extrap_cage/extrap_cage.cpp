#include "extrap_cage.h"
#include <igl/per_vertex_normals.h>
#include <igl/ray_mesh_intersect.h>
#include <igl/barycentric_to_global.h>
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
  
  Eigen::MatrixXd N1, N2;
  igl::per_vertex_normals(V1, F, N1);
  igl::per_vertex_normals(V2, F, N2);
  
  Eigen::VectorXd total_vertex_weights = Eigen::VectorXd::Zero(V2_C.rows());
  
  for (int vi = 0; vi < V1.rows(); vi++) {
    igl::Hit hit;
    igl::ray_mesh_intersect(V1.row(vi), N1.row(vi), V1_C, F_C, hit);
    Eigen::MatrixXd bc(1, 3);
    bc << hit.id, hit.u, hit.v;
    Eigen::RowVector3i coarse_face = F_C.row(hit.id);
    Eigen::RowVector3d coarse_point = igl::barycentric_to_global(V1_C, F_C, bc);
    
    Eigen::RowVector3d a = V1_C.row(coarse_face[0]);
    Eigen::RowVector3d b = V1_C.row(coarse_face[1]);
    Eigen::RowVector3d c = V1_C.row(coarse_face[2]);
    Eigen::RowVector3d l;
    igl::barycentric_coordinates(coarse_point, a, b, c, l);
    
    Eigen::RowVector3d t = V2.row(vi) - V1.row(vi) + N2.row(vi) * (coarse_point - V1.row(vi)).norm();
    V2_C_inc.row(coarse_face[0]) += l[0] * t;
    V2_C_inc.row(coarse_face[1]) += l[1] * t;
    V2_C_inc.row(coarse_face[2]) += l[2] * t;
    
    total_vertex_weights[coarse_face[0]] += l[0];
    total_vertex_weights[coarse_face[1]] += l[1];
    total_vertex_weights[coarse_face[2]] += l[2];
  }
  
  for (int vci = 0; vci < V2_C.rows(); vci++) {
    V2_C_inc.row(vci) /= total_vertex_weights[vci];
  }
  
  V2_C += V2_C_inc;
}