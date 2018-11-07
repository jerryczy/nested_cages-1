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
  
  for (int vci = 0; vci < V_C.rows(); vci++) {
  
    std::vector vci_faces = V_CF[vci];
    std::vector vci_face_inds = V_CI[vci];

    for (int s = 0; s < SAMPLE_SIZE; s++) {
    
      // pick a face incident on vc; TODO: probability must be based on face areas
      int fci = rand() % vci_faces.size();
      int fc = vci_faces[fci];
      int ic = vci_face_inds[ic];
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
        igl::point_simplex_squared_distance(x, V, F, i, &dist, &c);
        if (dist < smallest_dist) {
          smallest_dist = dist;
          closest_face = f;
        }
      }
      
      W(vci, closest_face) += 1.0 / SAMPLE_SIZE;
      
    }
  }
}