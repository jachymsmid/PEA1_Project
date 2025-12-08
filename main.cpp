#include <iostream>
#include <iomanip>
#include <sstream>
#include "utils/Mesh.hpp"
#include "utils/SimulationInfo.hpp"
#include "utils/NumericalSchemes.hpp"
#include "utils/InitialConditions.hpp"
#include "utils/BoundaryConditions.hpp"

// TODO:  - [x] better mesh representation
//            - [x] change the array of structer paradigm to structure of array
//                - [x] how to acces the data?
//            - [x] copy constructor
//            - [x] implement flattened arrays 
//               - [o] maybe Z-curves?
//                  - [o] conditions for Z-curves usage
//        - [ ] irregular grid method - too elaborate?
//        - [x] upwind method
//        - [x] lax-friedrichs method
//        - [x] lax-wendroff method
//            - [ ] calculate the spatial step in the function - preparation for irregular grid
//            - [x] clearer and simpler code
//        - [x] SimulationInfo struct
//            - [x] constructor - is it working??
//        - [x] how to implement the numerical solver
//            - [x] using the function template
//            - [o] using a pointer to a function
//        - [x] how to implement boundary/initial conditions?
//            - [x] again function templates
//        - [x] VTK library for data storage and visualization
//        - [o] error handling
//        - [o] generalize the schemes?

using RealNumber = float;
using NumericalScheme_ = Lax_Wendroff< RealNumber >;
using Mesh_ = Mesh< RealNumber >;
using InitialConditions_ = MyInitialConditions< RealNumber, Mesh_ >;
using BoundaryConditions_ = Zeros< RealNumber >;

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
  RealNumber dt = CFL< RealNumber, NumericalScheme_ >(dx, dy);

  SimulationInfo< RealNumber > sim_info(dt, dx, dy, x_min, y_min, x_max, y_max, N_x, N_y);
  

  // initialize the mesh, U and file_name
  Mesh_ mesh( sim_info ); 
  std::string file_name;

  mesh.construct_regular_grid( sim_info );

  InitialConditions< InitialConditions_, Mesh_ >( mesh );

  mesh.write_data( "sim/output_t_0.00000.vti" );

  // main loop
  while (t <= T)
  {
    t += dt;
    NumericalSolver< RealNumber, NumericalScheme_ >( mesh, sim_info );
    BoundaryConditions< RealNumber, BoundaryConditions_ >( mesh );

    // this is kinda awkward
    std::ostringstream fn;
    fn << "sim/output_t_" << std::fixed << std::setprecision(5) << t << ".vti";
    file_name = fn.str();

    mesh.write_data(file_name);
  }
  return 0;
}
