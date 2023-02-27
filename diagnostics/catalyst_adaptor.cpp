#include "catalyst_adaptor.hpp"

#if USE_CATALYST
#include "catalyst.hpp"
#include "catalyst_conduit_blueprint.hpp"
#include <iostream>
#include <mpi.h>
#include <sstream>

extern "C" {
void catalyst_adaptor() {}

void catalyst_adaptor_initialise(char *a_catalyst_scripts) {
  conduit_cpp::Node node;

  // Pass script to Catalyst
  int iscript = 0;
  std::stringstream catalyst_scripts_stream(a_catalyst_scripts);
  std::string script;
  while (std::getline(catalyst_scripts_stream, script, ':')) {
    node["catalyst/scripts/script" + std::to_string(iscript++)].set_string(
        script);
  }

  // Initialize Catalyst
  catalyst_status err1 = catalyst_initialize(conduit_cpp::c_node(&node));
  if (err1 != catalyst_status_ok) {
    std::cerr << "Failed to initialize Catalyst: " << err1 << std::endl;
  }

  conduit_cpp::Node catalyst_info;
  catalyst_status err2 = catalyst_about(conduit_cpp::c_node(&catalyst_info));
  if (err2 != catalyst_status_ok) {
    std::cerr << "Failed to get Catalyst implementation info: " << err2
              << std::endl;
  }
  int my_id;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
  if (my_id == 0) {
    std::cout << "**********************************\n";
    std::cout << "* Initializing Catalyst pipeline *\n";
    std::cout << "**********************************\n";
    std::cout << "Catalyst info:";
    if (catalyst_info.has_path("catalyst/tpl/conduit/license")) {
      // remove the long license string from the output
      catalyst_info.remove("catalyst/tpl/conduit/license");
    }
    std::string catalyst_info_yaml = catalyst_info.to_yaml();
    std::cout << catalyst_info_yaml;
  }
}

void catalyst_adaptor_execute(int *a_step_index, double *a_time) {
  int my_id;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_id);

  if (my_id == 0) {
    std::cout << "*******************************\n";
    std::cout << "* Executing Catalyst pipeline *\n";
    std::cout << "*******************************\n";
    std::cout << "Step " << *a_step_index << " (t_now = " << *a_time << ")";
    std::cout << std::endl;
  }

  // use variable names consistent with JOREK
  int nsub = 0;
  int n_elements = 0;
  int n_scalars = 0;

  // We want only rank 0 to pass the grid to Catalyst
  if (my_id == 0)
    catalyst_get_params(&nsub, &n_elements, &n_scalars);
  int nnos = nsub * nsub * n_elements; // Number of nodes in the Catalyst grid
  int nel = (nsub - 1) * (nsub - 1) *
            n_elements; // Number of cells in the Catalyst grid
  int nnoel = 4;        // Number of points needed to define a cell

  std::vector<std::string> scalar_names(n_scalars);
  if (my_id == 0) {
    char name[36];
    for (int iscalar = 0; iscalar < n_scalars; ++iscalar) {
      int f_iscalar = iscalar + 1; // Fortran is 1-indexed
      catalyst_get_scalar_name(name, &f_iscalar);
      scalar_names[iscalar] = name;
    }
  }

  conduit_cpp::Node exec_params;

  // add time/cycle information
  auto state = exec_params["catalyst/state"];
  state["timestep"].set(*a_step_index);
  state["time"].set(*a_time);
  state["multiblock"].set(0);

  // Add channels.
  // We only have 1 channel here. Let's name it 'grid'.
  auto channel = exec_params["catalyst/channels/grid"];

  // Use Conduit Mesh Blueprint
  channel["type"].set("mesh");

  // now create the mesh.
  auto mesh = channel["data"];

  mesh["coordsets/coords/type"].set("explicit");

  std::vector<float> coords_R(nnos);
  std::vector<float> coords_Z(nnos);
  std::vector<int> cell_points(nnoel * nel);

  if (my_id == 0) {
    catalyst_interp_grid(&nnos, &nel, coords_R.data(), coords_Z.data(),
                         cell_points.data());
  }

  // set coordinates of nodes
  mesh["coordsets/coords/type"].set("explicit");
  mesh["coordsets/coords/values/r"].set_external(coords_R);
  mesh["coordsets/coords/values/z"].set_external(coords_Z);

  // set topology
  mesh["topologies/mesh/type"].set("unstructured");
  mesh["topologies/mesh/coordset"].set("coords");
  mesh["topologies/mesh/elements/shape"].set("quad");
  mesh["topologies/mesh/elements/connectivity"].set_external(cell_points);

  // store all scalars in this 2D array
  std::vector<std::vector<float>> scalars(n_scalars);

  if (my_id == 0) {
    // add scalars
    auto fields = mesh["fields"];
    for (int iscalar = 0; iscalar < n_scalars; ++iscalar) {
      scalars[iscalar].resize(nnos);
      int f_iscalar = iscalar + 1; // Fortran is 1-indexed
      catalyst_interp_scalar(scalars[iscalar].data(), &f_iscalar);
      // auto scalar = fields[scalar_names[iscalar]];
      std::string scalar_path = scalar_names[iscalar] + "/";
      fields[scalar_path + std::string("association")].set("vertex");
      fields[scalar_path + std::string("topology")].set("mesh");
      fields[scalar_path + std::string("volume_dependent")].set("false");
      fields[scalar_path + std::string("values")].set_external(
          scalars[iscalar]);
    }
  }

  conduit_cpp::Node mesh_info;
  if (!conduit_cpp::BlueprintMesh::verify(mesh, mesh_info)) {
    std::cerr << "Grid does not satisfy Conduit mesh blueprint" << std::endl;
    mesh_info.print();
  }

  catalyst_status err = catalyst_execute(conduit_cpp::c_node(&exec_params));
  if (err != catalyst_status_ok) {
    std::cerr << "Failed to execute Catalyst: " << err << std::endl;
  }

  if (my_id == 0)
    std::cout << std::endl;
}

void catalyst_adaptor_finalise() {
  conduit_cpp::Node node;

  // Finalize Catalyst
  catalyst_status err = catalyst_finalize(conduit_cpp::c_node(&node));
  if (err != catalyst_status_ok) {
    std::cerr << "Failed to finalize Catalyst : " << err << std::endl;
  }
}
}

#endif