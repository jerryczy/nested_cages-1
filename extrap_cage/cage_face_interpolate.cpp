#define SAMPLE_SIZE 100

#include "cage_face_interpolate.h"
#include "float.h"
#include <vector>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/doublearea.h>
#include <igl/point_mesh_squared_distance.h>

void cage_face_interpolate(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::MatrixXd & V_C,
  const Eigen::MatrixXi & F_C,
  Eigen::MatrixXd & W)
{
  W = Eigen::MatrixXd::Zero(V_C.rows(), F.rows());
  
  std::vector<std::vector<int>> V_CF, V_CI;
  igl::vertex_triangle_adjacency(V_C.rows(), F_C, V_CF, V_CI);
  
  Eigen::VectorXd dblA;
  igl::doublearea(V_C, F_C, dblA);
  
  Eigen::MatrixXd P(SAMPLE_SIZE, 3);

  for (int vci = 0; vci < V_C.rows(); vci++) {
  
    std::vector<int> vci_faces = V_CF[vci];
    std::vector<int> vci_face_inds = V_CI[vci];
    
    // calculate cumulative areas of incident faces (TODO: use voronoi areas)
    double total_area = 0;
    std::vector<double> cum_area;
    int n = vci_faces.size();
    for (int i = 0; i < n; i++) {
      double area = dblA(vci_faces[i]);
      total_area += area;
      cum_area.push_back(total_area);
    }

    for (int s = 0; s < SAMPLE_SIZE; s++) {
    
      // pick a face incident on vc;
      double idx = total_area * (std::rand() / double(RAND_MAX));
      int fci;
      for (int i = 0; i < n; i++) {
        if (cum_area[i] > idx) {
          fci = i;
          break;
        }
      }
      int fc = vci_faces[fci];
      int ic = vci_face_inds[fci];
      
      // obtain coordinates of face
      Eigen::Vector3d v1 = V_C.row(vci);
      Eigen::Vector3d v2 = V_C.row(F_C(fc, (ic + 1) % 3));
      Eigen::Vector3d v3 = V_C.row(F_C(fc, (ic + 2) % 3));
    
      // pick a random point in the voronoi region corresponding to vc
      double alpha = std::rand() / (2 * double(RAND_MAX));
      double beta = std::rand() / (2 * double(RAND_MAX));
      if (alpha + beta > 1) {
        alpha = 1 - alpha;
        beta = 1 - beta;
      }
      P.row(s) = v1 + alpha * (v2 - v1) + beta * (v3 - v1);   
    }
    
    Eigen::VectorXd I(SAMPLE_SIZE), sqrD;
    Eigen::MatrixXd C(SAMPLE_SIZE, 3);
    igl::point_mesh_squared_distance(P, V, F, sqrD, I, C);
    
    // update weights according to the closest face to each sampled point
    for (int s = 0; s < I.rows(); s++) {
      W(vci, I[s]) += 1.0 / SAMPLE_SIZE;
    }
  }
}