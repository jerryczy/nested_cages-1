#define SAMPLE_SIZE 100

#include "cage_face_interpolate.h"
#include <igl/point_simplex_squared_distance.h>

void cage_face_interpolate(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXd & F,
  const Eigen::MatrixXd & V_C,
  const Eigen::MatrixXd & F_C,
  Eigen::MatrixXd & W)
{
  W = Eigen::VectorXd::Zeros(V_C.rows(), F.rows());
  
  std::vector<std::vector<int>> V_CF, V_CI;
  igl::vertex_triangle_adjacency(V_C.rows(), F_C, V_CF, V_CI);
  Eigen::VectorXd dblA;
  igl::doublearea(V_C, F_C, dblA);

  for (int vci = 0; vci < V_C.rows(); vci++) {
  
    std::vector vci_faces = V_CF[vci];
    std::vector vci_face_inds = V_CI[vci];

    for (int s = 0; s < SAMPLE_SIZE; s++) {
    
      // pick a face incident on vc;
      double total_area;
      std::vector<double> cum_area;
      int n = vci_faces.size();
      for (int i = 0; i < n; i++) {
        double area = dblA(vci_faces[i]);
        total_area += area;
        cum_area.push_back(total_area)
      }
      int idx = std::rand() % std::floor(total_area);
      int fci;
      for (int i = n-1; i >= 0; i--) {
        if (idx <= cum_area(i)) {
          fci = i;
          break;
        }
      }
      int fc = vci_faces[fci];
      int ic = vci_face_inds[fci];
      Eigen::Vector3d v1 = V_C.row(vci);
      Eigen::Vector3d v2 = V_C.row(F_C(fc, (ic + 1) % 3));
      Eigen::Vector3d v3 = V_C.row(F_C(fc, (ic + 2) % 3));
    
      // pick a random point in this voronoi region
      double alpha = std::rand() / (2 * double(RAND_MAX));
      double beta = std::rand() / (2 * double(RAND_MAX));
      if (alpha + beta > 1) {
        alpha = 1 - alpha;
        beta = 1 - beta;
      }
      Eigen::Vector3d x = v1 + alpha * (v2 - v1) + beta * (v3 - v1);
      
      // find mesh face closest to x; TODO: there's probably a more efficient way
      double smallest_dist = DBL_MAX;
      int closest_face = 0;
      for (int f = 0; f < F.rows(); f++) {
        double dist, c;
        igl::point_simplex_squared_distance(x, V, F, f, &dist, &c);
        if (dist < smallest_dist) {
          smallest_dist = dist;
          closest_face = f;
        }
      }
      
      W(vci, closest_face) += 1.0 / SAMPLE_SIZE;
      
    }
  }
}