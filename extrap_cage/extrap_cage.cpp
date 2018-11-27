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
  
  Eigen::MatrixXd N;
  igl::per_vertex_normals(V1, F, N);
  
  Eigen::VectorXd total_vertex_weights = Eigen::VectorXd::Zero(V2_C.rows());
  
  for (int vi = 0; vi < V1.rows(); vi++) {
    igl::Hit hit;
    igl::ray_mesh_intersect(V1.row(vi), N.row(vi), V1_C, F_C, hit);
    Eigen::MatrixXd bc(1, 3);
    bc << hit.id, hit.u, hit.v;
    Eigen::RowVector3i coarse_face = F_C.row(hit.id);
    Eigen::RowVector3d coarse_point = igl::barycentric_to_global(V1_C, F_C, bc);
    
    Eigen::RowVector3d a = V1_C.row(coarse_face[0]);
    Eigen::RowVector3d b = V1_C.row(coarse_face[1]);
    Eigen::RowVector3d c = V1_C.row(coarse_face[2]);
    Eigen::RowVector3d l;
    igl::barycentric_coordinates(coarse_point, a, b, c, l);
    
    V2_C_inc.row(coarse_face[0]) += l[0] * (V2.row(vi) - V1.row(vi));
    V2_C_inc.row(coarse_face[1]) += l[1] * (V2.row(vi) - V1.row(vi));
    V2_C_inc.row(coarse_face[2]) += l[2] * (V2.row(vi) - V1.row(vi));
    
    total_vertex_weights[coarse_face[0]] += l[0];
    total_vertex_weights[coarse_face[1]] += l[1];
    total_vertex_weights[coarse_face[2]] += l[2];
  }
  
  for (int vci = 0; vci < V2_C.rows(); vci++) {
    V2_C_inc.row(vci) /= total_vertex_weights[vci];
  }
  
  V2_C += V2_C_inc;
  
  // patchwork: find portruding translations and translate corresponding triangles
  for (int vi = 0; vi < V1.rows(); vi++) {
    if ((V2.row(vi) - V1.row(vi)).dot(N.row(vi)) > 0) {
      igl::Hit hit;
      igl::ray_mesh_intersect(V1.row(vi), V2.row(vi) - V1.row(vi), V2_C, F_C, hit);
      Eigen::MatrixXd bc(1, 3);
      bc << hit.id, hit.u, hit.v;
      Eigen::RowVector3i coarse_face = F_C.row(hit.id);
      Eigen::RowVector3d coarse_point = igl::barycentric_to_global(V2_C, F_C, bc);
      
      if ((V2.row(vi) - V1.row(vi)).norm() > (coarse_point - V1.row(vi)).norm()) {
        V2_C.row(coarse_face[0]) += V2.row(vi) - coarse_point;
        V2_C.row(coarse_face[1]) += V2.row(vi) - coarse_point;
        V2_C.row(coarse_face[2]) += V2.row(vi) - coarse_point;
      }
    }
  }
}