#include "extrap_cage.h"
#include "recover_affine_transformations.h"
#include <igl/per_vertex_normals.h>
#include <igl/ray_mesh_intersect.h>
#include <igl/barycentric_to_global.h>
#include <igl/barycentric_coordinates.h>
#include <igl/svd3x3.h>

void extrap_cage(
  const Eigen::MatrixXd & V1,
  const Eigen::MatrixXd & V2,
  const Eigen::MatrixXi & F,
  const Eigen::MatrixXd & V1_C,
  const Eigen::MatrixXi & F_C,
  Eigen::MatrixXd & V2_C)
{
  V2_C = Eigen::MatrixXd::Zero(V1_C.rows(), 3);
  
  Eigen::MatrixXd N;
  igl::per_vertex_normals(V1, F, N);
  
  Eigen::MatrixXd VA;
  recover_affine_transformations(V1, V2, F, VA);
  
  Eigen::MatrixXd VA_C(V2_C.rows(), 9);
  
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

    VA_C.row(coarse_face[0]) += l[0] * VA.row(vi);
    VA_C.row(coarse_face[1]) += l[1] * VA.row(vi);
    VA_C.row(coarse_face[2]) += l[2] * VA.row(vi);
    
    total_vertex_weights[coarse_face[0]] += l[0];
    total_vertex_weights[coarse_face[1]] += l[1];
    total_vertex_weights[coarse_face[2]] += l[2];
  }
  
  for (int vci = 0; vci < V2_C.rows(); vci++) {
    Eigen::Matrix3d A;
    Eigen::RowVector3d t;
    A << VA_C.block(vci, 0, 1, 3), VA_C.block(vci, 3, 1, 3), VA_C.block(vci, 6, 1, 3);
    t << VA_C.block(vci, 9, 1, 3);
    
    Eigen::Matrix3d U, V;
    Eigen::Vector3d S;
    igl::svd3x3(A, U, S, V);
    Eigen::Matrix3d Sigma;
    Sigma << 1, 0, 0,
             0, 1, 0,
             0, 0, (U * V.transpose()).determinant();
    Eigen::Matrix3d R = U * Sigma * V.transpose();
    
    V2_C.row(vci) = (1 / total_vertex_weights[vci]) * V1_C.row(vci) * R.transpose() + t;
  }
}