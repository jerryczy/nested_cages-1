// Meshfix include
#define MESHFIX_WITH_EIGEN
#include <meshfix.h>

// Our header files
#include "io.h"
#include "cgal.h"
#include "flow.h"
#include "reinflate.h"
#include "extrap_cage/extrap_cage.h"

// libigl includes
#include <igl/doublearea.h>
#include <igl/writeOBJ.h>
#include <igl/read_triangle_mesh.h>
#include <igl/copyleft/cgal/remesh_self_intersections.h>
#include <igl/copyleft/cgal/polyhedron_to_mesh.h>
#include <igl/copyleft/cgal/mesh_to_polyhedron.h>
#include <igl/copyleft/cgal/intersect_other.h>
#include <igl/opengl/glfw/Viewer.h>
 
// useful namespaces
using namespace Eigen;
using namespace igl;
using namespace std;

// For debuggin'
int at(
  Eigen::MatrixXi & M,
  const int i,
  const int j)
{
  return M(i,j);
}

void view_cages(igl::opengl::glfw::Viewer &viewer) {
  // based on tutorial code (107)
  int k = viewer.data_list.size();
  std::map<int, Eigen::RowVector3d> colors;
  for (int i = 0; i < k; i++)
    colors.emplace(i, 0.5*Eigen::RowVector3d::Random().array() + 0.5);
  int last_selected = -1;
  viewer.callback_key_down =
    [&](igl::opengl::glfw::Viewer &, unsigned int key, int mod)
  {
    if(key == GLFW_KEY_UP)
    {
      viewer.selected_data_index = (viewer.selected_data_index + 1) % k;
      return true;
    }
    else if(key == GLFW_KEY_DOWN)
    {
      viewer.selected_data_index = (k + viewer.selected_data_index - 1) % k;
      return true;
    }
    return false;
  };
  viewer.callback_pre_draw =
    [&](igl::opengl::glfw::Viewer &)
  {
    if (last_selected != viewer.selected_data_index)
    {
      for (auto &data : viewer.data_list)
      {
        data.set_colors(colors[data.id]);
      }
      viewer.data_list[viewer.selected_data_index].set_colors(Eigen::RowVector3d(1., 1., 1.));
      last_selected = viewer.selected_data_index;
    }
    return false;
  };
  viewer.launch();
}

int main(int argc, char * argv[])
{
  using namespace igl::copyleft::cgal;

  if (argc==1)
  {
    cout << R"(Usage: 

    ./nested_cages input.off q L(1) L(2) ... L(k) EnergyExpansion EnergyFinal output

input: the program accepts files in the following formats: .off, .obj, .ply, .stl, .wrl .mesh

output: cages will be saved as output_1.obj, output_2.obj, ..., output_k.obj

q is the quadrature order for the shrinking flow

L(1) > L(2) > ... > L(k) is the number of faces for each cage.
If L(k) is followed by 'r' the initial decimation for this cage will be regular
(adaptive if no 'r').
Each L(k) can be replace by a filed with an input decimation.

EnergyExpansion is the energy to be minimized for the re-inflation
Energies implemented: None, DispStep, DispInitial, Volume, SurfARAP, VolARAP 

EnergyFinal is the energy to be minimized after the re-inflation (additional processing)
Energies implemented: None, DispStep, DispInitial, Volume, SurfARAP, VolARAP
)";
    return EXIT_FAILURE;
  }
  
  igl::opengl::glfw::Viewer viewer;
  
  if (strncmp(argv[1], "view", 4) == 0) {
    char *filename;
    char *suffix;
    asprintf(&filename, "%s%s", argv[2], "_0.obj");
    std::ifstream testfile(filename);
    int i = 0;
    while (testfile) {
      viewer.load_mesh_from_file(filename);
      asprintf(&suffix,"_%d.obj",++i);
      asprintf(&filename, "%s%s", argv[2], suffix);
      testfile = std::ifstream(filename);
    }
    view_cages(viewer);
    return 0;
  }
  
  // TODO: move this functionality to a separate file
  if (strncmp(argv[1], "extrap", 4) == 0) {
    MatrixXd V1, V2, V1_C, V2_C_original, V2_C;
    MatrixXi F, F1_C, F2_C_original;
    if (!read_triangle_mesh(argv[2],V1,F) ||
        !read_triangle_mesh(argv[3],V2,F) ||
        !read_triangle_mesh(argv[4],V1_C,F1_C) ||
        !read_triangle_mesh(argv[5],V2_C_original,F2_C_original)){
      cout << "unable to read input file"  << endl;
      return 0;
    }
    extrap_cage(V1, V2, F, V1_C, F1_C, V2_C);

    viewer.data().set_mesh(V1, F);
    viewer.data().set_face_based(true);
    viewer.append_mesh();
    viewer.data().set_mesh(V2, F);
    viewer.data().set_face_based(true);
    viewer.append_mesh();
    viewer.data().set_mesh(V1_C, F1_C);
    viewer.data().set_face_based(true);
    //viewer.append_mesh();
    //viewer.data().set_mesh(V2_C_original, F2_C_original);
    //viewer.data().set_face_based(true);
    viewer.append_mesh();
    viewer.data().set_mesh(V2_C, F1_C);
    viewer.data().set_face_based(true);
    
    // TODO: this is not the best interface but will do for now
    view_cages(viewer);
    return 0;
  }
    
  // number of layers
  int k = argc-6;
  #ifdef VERBOSE_DEBUG
    cout << "number of layers = " << k << endl;
  #endif

  // quadrature order
  int quad_order = atoi(argv[2]);

  // read input mesh
  Surface_mesh M; 
  MatrixXd V0;
  MatrixXi F0;
  if (!read_triangle_mesh(argv[1],V0,F0)){
    cout << "unable to read input file"  << endl;
    return 0;
  }
  // convert to CGAL format
  mesh_to_polyhedron(V0,F0,M);

  // output input mesh as level output_0.off
  char* filename; 
  char *suffix;
  if (asprintf(&filename, "%s%s", argv[argc-1], "_0.obj")!=-1){
    writeOBJ(filename,V0,F0);
  } else {
    cout << "unable to allocate space for output file name"  << endl;
    return 0;
  }

  // First fine mesh is the input mesh
  MatrixXd V = V0;
  MatrixXi F = F0;
  
  // add initial mesh to viewer
  viewer.data().set_mesh(V, F);
  viewer.data().set_face_based(true);

  // vector where each entry is the number of faces for eachj level
  int L[k];
  // standard CGAL decimation is adaptive
  bool adaptive = true;

  // declare input decimation
  Surface_mesh M_hat;
  MatrixXd V_coarse;
  MatrixXi F_coarse;

  // Loop over levels
  for(int i = 0;i<k;i++){
    // check if argv[i+3] is a valid file. If so, it will be used as input decimation 
    std::ifstream is_file(argv[i+3]);
    if (is_file)
    {
      if (!read_triangle_mesh(argv[i+3],V_coarse,F_coarse)){
        cout << "error: input decimation is an existing file, but couldn't be read"  << endl;
        return 0;
      }
      // convert to CGAL format
      mesh_to_polyhedron(V_coarse,F_coarse,M_hat);
    } 
    // otherwise the input decimation is computed by decimating the previous layer
    // with CGAL decimation (regular or adaptive)
    else
    {
      // first check if last charcater of argv[i+2] is r. If it is, drop the 'r' and adaptive = false
      adaptive = remove_all_chars_and_count(argv[i+3], 'r')==0;
      // check if argv[i+2] is a valid integer (throw an error if it is not)
      if (!legal_int(argv[i+3])){
        cout << "you have to pass integer values or valid input deimatations"  << endl;
        cout << "the invalid argument you have passed is " << argv[i+3] << endl;
        return 0;
      }
      // specified number of faces for this cage
      L[i] = atoi(argv[i+3]);
      // value to pass to the decimator
      float ratio = (1.*L[i])/(1.*M.size_of_facets()); 
      // the previously computed will be decimated
      M_hat = M;
      // decimate
      decimate_CGAL(&M_hat,ratio,adaptive);
    }

    // Convert decimation to LibIGL/Eigen format
    polyhedron_to_mesh(M_hat,V_coarse,F_coarse); 
    // Parameters to call function to check for decimation's self-intersections
    RemeshSelfIntersectionsParam params;
    params.detect_only = true;
    params.first_only = true;
    MatrixXd tempV;
    MatrixXi tempF;
    MatrixXi IF;
    VectorXi J;
    VectorXi IM;
    remesh_self_intersections(V_coarse,F_coarse,params,tempV,tempF,IF,J,IM);
    // If input coarse mesh self-intersect, remove self-intersections with Meshfix (to-do)
    if (IF.rows()>0){
      #ifdef VERBOSE_DEBUG
    	  cout << i+1 << "-th input decimation self-intersects. Fixing with Meshfix " << endl;
      #endif
      cout << "Polishing M" << i+1 << "..." << endl;
      meshfix(MatrixXd(V_coarse),MatrixXi(F_coarse),V_coarse,F_coarse);
      cout << "Success!" << endl;
      remesh_self_intersections(V_coarse,F_coarse,params,tempV,tempF,IF,J,IM);
      if (IF.rows()==0)
      {
        #ifdef VERBOSE_DEBUG
          cout << "Meshfix succesfully removed self-intersections" << endl;  
        #endif
      }
      else{
        #ifdef VERBOSE_DEBUG
          cout << "Meshfix wasn't able to remove all self-intersections. Quitting..." << endl; 
        #endif
    	  return 0;
      }
    }

	  // calculate triangle areas for initial mesh (will be used
	  // to define the integral at _every_ step - the metric is fixed)
	  VectorXd area_0;
	  doublearea(V,F,area_0);
	  area_0 = 0.5*area_0;
    // Precompute matrix that convert gradients at quadrature points to gradients at mesh vertices
  	SparseMatrix<double> A_qv;
    gradQ_to_gradV(V, F, area_0, quad_order, A_qv);
  	// Flow M inside M_hat and save the result to a stack M of flow meshes
    stack<MatrixXd> H;
    cout << "Flowing M" << i << " inside M" << i+1 << "..." << endl;
    if (!flow_fine_inside_coarse(V,F,V_coarse,F_coarse,A_qv,H))
    {
      cout << "Flow failed to take fine mesh inside coarse mesh after 1000 iterations. Quitting" << endl;
      return 0;
    }
    cout << "Success!" << endl;

    // Reinflate and output to cage to C
    MatrixXd C;
    cout << "Reinflating M" << i << ", pushing M" << i+1 << "..." << endl;
    reinflate(H,F,V_coarse,F_coarse,argv[argc-3],argv[argc-2],C);
    cout << "Success!" << endl;

    // sanity check: cage should never self-intersect at this stage
    remesh_self_intersections(C,F_coarse,params,tempV,tempF,IF,J,IM);
    if (IF.rows()>0){
      cout << i+1 << "-th output cage self-intersects. ERROR! Quitting...  " << endl;  
      return 0;
    }

    // sanity check: cage should never intersect input decimation
    intersect_other(C,F_coarse,V,F,true,IF);
    if (IF.rows()>0){
      cout << i+1 << "-th output cage intersect previous cage. ERROR! Quitting...  " << endl;  
      return 0;
    }

    // output cage is the input for the next level
    M.clear();
    mesh_to_polyhedron(C,F_coarse,M); 
    V = C;
    F = F_coarse;
    
    // add mesh to viewer
    viewer.append_mesh();
    viewer.data().set_mesh(V, F);
    viewer.data().set_face_based(true);

    // Output cage to file output_i.obj
    if ((asprintf(&suffix,"_%d.obj",i+1)!=-1) && (asprintf(&filename, "%s%s", argv[argc-1], suffix)!=-1)){
      writeOBJ(filename,C,F_coarse);
    } else {
      cout << "unable to allocate space for output file name"  << endl;
      return 0;
    }
    // back to adaptive decimation (standard)
    adaptive = true;
  }
  
  view_cages(viewer);
  
  free(filename);
  return 1;
}
