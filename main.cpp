#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <iomanip>
#include <sstream>

// TODO:  - [ ] better mesh representation
//            - [x] copy constructor
//            - [x] implement flattened arrays 
//               - [ ] maybe Z curves?
//            - [ ] irregular grid method - too elaborate?
//        - [ ] upwind method
//        - [ ] lax-friedrichs method
//        - [ ] lax-wendroff method
//            - [ ] calculate the spatial step in the function - preparation for irregular grid
//            - [ ] clearer and simpler code
//        - [x] SimulationInfo struct
//            - [x] constructor - is it working??
//        - [x] how to implement the numerical solver
//            - [x] using the function template
//            - [o] using a pointer to a function
//        - [ ] how to implement boundary/initial conditions?

using RealNumber = float;

// struct to hold infrmation about the simulation
struct SimulationInfo
{
  RealNumber step_t, step_x, step_y, x_min, y_min, x_max, y_max;
  int number_x, number_y;

  // constructor
  SimulationInfo( const RealNumber step_t,
                  const RealNumber step_x,
                  const RealNumber step_y,
                  const RealNumber x_min,
                  const RealNumber y_min,
                  const RealNumber x_max,
                  const RealNumber y_max,
                  const int number_x,
                  const int number_y)
                : step_t(step_t),
                  step_x(step_x),
                  step_y(step_y),
                  x_min(x_min),
                  y_min(y_min),
                  x_max(x_max),
                  y_max(y_max),
                  number_x(number_x),
                  number_y(number_y) {}
};

// struct to hold cartesian coordinates of a point and a value at that point  
struct DataPoint
{
  RealNumber x, y, value;
};

// mesh class
//  - glorified array of DataPoints
//  - flattened array  
//  - with copy constructor
//  - methods
//    - write_data
//    - getRows
//    - getCols
//    - print_data
class Mesh {
    std::vector<DataPoint> data;
    size_t rows, cols;

public:
    // constructor
    Mesh(SimulationInfo sim_info)
        : data(sim_info.number_x * sim_info.number_y), rows(sim_info.number_x), cols(sim_info.number_y) {}

    // copy constructor
    Mesh(Mesh &mesh)
        : data(mesh.data), rows(mesh.rows), cols(mesh.cols) {}

    DataPoint& operator()(size_t i, size_t j)
    {
        return data[i * cols + j];
    }
    const DataPoint& operator()(size_t i, size_t j) const
    {
        return data[i * cols + j];
    }

    // getters for number of columns and rows
    size_t getCols() { return cols; }
    size_t getRows() { return rows; }

    // method for regular grid construction
    void construct_grid(SimulationInfo sim_info)
    {
      for (size_t i = 0; i < rows; ++i)
      {
        for (size_t j = 0; j < cols; ++j)
        {
            data[i * cols + j].x = sim_info.x_min + i * sim_info.step_x;
            data[i * cols + j].y = sim_info.y_min + j * sim_info.step_y;
        }
      }
    }

    // write the data to a file
    void write_data(std::string file_name)
    {
      std::ofstream file;
      file.open(file_name);

      if (file)
      {
        for (size_t i = 0; i < rows; i++)
        {
          for (size_t j = 0; j < cols;j++)
          {
            // we want to end every line with just a \n character
            if (j == cols - 1)
            {
              file << data[i*cols + j].value << std::endl;
            }
            else
            {
              file << data[i*cols + j].value << ","; 
            }
          }
        }
      }
      file.close();
    }

    // print the mesh to the command line, for troubleshooting
    void print_mesh(Mesh &mesh)
    {
      for (size_t i = 0; i < mesh.rows; i++)
      {
        for (size_t j = 0; j < mesh.cols; j++)
        {
          std::cout << mesh(i,j).value << " ";
        }
        std::cout << std::endl;
      }
    }
};

// numerical solvers structs
// lax-friedrichs scheme
struct Lax_Friedrichs
{
  static void func()
  {
  }
};

// lax-wendroff scheme
struct Lax_Wendroff
{
  static void func(Mesh &mesh, SimulationInfo sim_info)
  {
    Mesh U_n(mesh);
    RealNumber lx, ly, dx, dy, dt;
    dx = sim_info.step_x;
    dy = sim_info.step_y;
    dt = sim_info.step_t;
    lx = dt / dx;
    for(size_t i = 1; i < mesh.getRows() - 1; i++)
    {
      for(size_t j = 1; j < mesh.getCols() - 1; j++)
      {
        ly = 2*mesh(i,j).x*dt/dy;
        mesh(i,j).value = U_n(i,j).value
                - lx / 2 * ( U_n(i + 1,j).value - U_n(i - 1,j).value )
                - ly / 2 * ( U_n(i,j + 1).value - U_n(i,j - 1).value )
                + pow(lx,2.0) / 2 * ( U_n(i + 1,j).value - 2 * U_n(i,j).value + U_n(i - 1,j).value )
                + pow(ly,2.0) / 2 * (U_n(i,j + 1).value - 2 * U_n(i,j).value + U_n(i,j - 1).value )
                + ly * lx / 4 * (U_n(i + 1,j + 1).value - U_n(i - 1,j + 1).value - U_n(i + 1,j - 1).value + U_n(i - 1,j - 1).value );
      }
    }
  }
};

// upwind scheme
struct Upwind
{
  static void func()
  {
  }
};

// general function to call specific numerical scheme
template <typename T>
void NumericalSolver(Mesh &mesh, SimulationInfo sim_info)
{
  T::func(mesh, sim_info);
}



// impose initial conditions
void initial_conditions(Mesh &mesh)
{
  for (size_t i = 0; i < mesh.getRows(); i++)
  {
    for (size_t j = 0; j < mesh.getCols(); j++)
    {
      mesh(i,j).value = 100 * exp( - ( pow(mesh(i,j).x + 1, 2.0)+pow(mesh(i,j).y, 2.0) ) / 0.01);
    }
  }
}

// impose boundary conditions
void boundary_conditions(Mesh &mesh)
{
  for (size_t i = 0; i < mesh.getRows(); ++i)
  {
    for (size_t j = 0; j < mesh.getCols(); ++j)
    {
      if (i == 0 || i == mesh.getRows() - 1 || j == 0 || j == mesh.getCols() - 1)
      {
        mesh(i,j).value = 0; 
      }
    }
  }
}



// numerical schemes
void lax_friedrichs()
{
}

void upwind()
{
}

// lax-wendroff scheme


RealNumber cfl_lax_wendroff(RealNumber dx, RealNumber dy)
{
  return 1.0 / (1.0 / dx + 1.0 / 2.0 * 1.5 * dy);
}

int main()
{
  // domain definition
  const RealNumber x_min = - 1.5, x_max = 1.5, y_min = - 1.5, y_max = 1.5;
  // number of grid points in the x, y direction respectively
  const int N_x = 10, N_y = 10;
  const RealNumber dx = (x_max-x_min) / N_x, dy = (y_max-y_min) / N_y;
  const RealNumber T = 2.f;
  RealNumber t = 0.f;
  const std::string sim_name;
  RealNumber dt = cfl_lax_wendroff(dx, dy);

  SimulationInfo sim_info(dt, dx, dy, x_min, y_min, x_max, y_max, N_x, N_y);
  

  // initialize the mesh, U and file_name
  std::cout << "Simulation info initialization..." << std::endl;
  Mesh mesh(sim_info); 
  std::cout << "Mesh initialization..." << std::endl;
  std::string file_name;

  mesh.construct_grid(sim_info);
  std::cout << "Grid construction..." << std::endl;

  initial_conditions(mesh);
  std::cout << "Imposing initial conditions..." << std::endl;

  mesh.write_data("sim/output_t_0.00000.csv");

  // std::cout << "Enter simulation name (folder with the same name will be created): ";
  // std::cin << sim_name;

  // main loop
  while (t <= T)
  {
    t += dt;
    NumericalSolver< Lax_Wendroff >(mesh, sim_info);
    boundary_conditions(mesh);

    // this is kinda awkward
    std::ostringstream fn;
    fn << "sim/output_t_" << std::fixed << std::setprecision(5) << t << ".csv";
    file_name = fn.str();

    mesh.write_data(file_name);
  }
  return 0;
}
