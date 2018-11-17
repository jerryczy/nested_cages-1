#define SAMPLE_SIZE 100

#include "cage_vertex_interpolate.h"
#include "float.h"
#include <vector>
#include <igl/point_mesh_squared_distance.h>

void cage_vertex_interpolate(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::MatrixXd & V_C,
  const Eigen::MatrixXi & F_C,
  Eigen::MatrixXd & W)
{
  W = Eigen::MatrixXd::Zero(V_C.rows(), V.rows());
  
  Eigen::VectorXd sqrD, I;
  Eigen::MatrixXd C;
  igl::point_mesh_squared_distance(V_C, V, F, sqrD, I, C);

  for (int vci = 0; vci < V_C.rows(); vci++) {
    Eigen::RowVector3d v = V_C.row(vci);
    Eigen::RowVector3d closest_point = C.row(vci);
    Eigen::RowVector3i closest_face = F.row(I[vci]);
    Eigen::RowVector3d a = V.row(closest_face[0]);
    Eigen::RowVector3d b = V.row(closest_face[1]);
    Eigen::RowVector3d c = V.row(closest_face[2]);
    Eigen::RowVector3d l;
    igl::barycentric_coordinates(v, a, b, c, l);
    W(vci, closest_face[0]) = 1 - l[0];
    W(vci, closest_face[1]) = 1 - l[1];
    W(vci, closest_face[2]) = 1 - l[2];
  }
}