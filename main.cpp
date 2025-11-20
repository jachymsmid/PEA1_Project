#include <iostream>
#include <iomanip>
#include <sstream>
#include "functions.h"

//

using NumericalScheme = Lax_Wendroff;

int main()
{
  // domain definition
  const RealNumber x_min = - 1.5, x_max = 1.5, y_min = - 1.5, y_max = 1.5;
  // number of grid points in the x, y direction respectively
  const int N_x = 64, N_y = 64;
  const RealNumber dx = (x_max-x_min) / N_x, dy = (y_max-y_min) / N_y;
  const RealNumber T = 2.f;
  RealNumber t = 0.f;
  const std::string sim_name;
  RealNumber dt = CFL< NumericalScheme>(dx, dy);

  SimulationInfo sim_info(dt, dx, dy, x_min, y_min, x_max, y_max, N_x, N_y);
  

  // initialize the mesh, U and file_name
  Mesh mesh(sim_info); 
  std::string file_name;

  mesh.construct_regular_grid(sim_info);

  InitialConditions<My_Initial_Conditions >(mesh);
  std::cout << "Imposing initial conditions..." << std::endl;

  mesh.write_data("sim/output_t_0.00000.vti");

  // std::cout << "Enter simulation name (folder with the same name will be created): ";
  // std::cin << sim_name;

  // main loop
  while (t <= T)
  {
    t += dt;
    NumericalSolver< NumericalScheme >(mesh, sim_info);
    BoundaryConditions< Zeros >(mesh);

    // this is kinda awkward
    std::ostringstream fn;
    fn << "sim/output_t_" << std::fixed << std::setprecision(5) << t << ".vti";
    file_name = fn.str();

    mesh.write_data(file_name);
  }
  return 0;
}
