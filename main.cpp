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
//            - [ ] irregular grid method
//        - [ ] upwind method
//        - [ ] lax-friedrichs method
//        - [ ] lax-wendroff method
//            - [ ] calculate the spatial step in the function
//        - [ ] SimulationInfo struct
//            - [x] constructor
//            - [ ] constructor for default values - should there be one?
//        - [ ] how to implement the numerical solver
//        - [ ] how to implement boundary/initial conditions?

using RealNumber = float;

// struct to hold infrmation about the simulation
struct SimulationInfo
{
  RealNumber step_t, step_x, step_y, x_min, y_min, x_max, y_max;
  int number_x, number_y;

  // constructor
  SimulationInfo( RealNumber step_t,
                  RealNumber step_x,
                  RealNumber step_y,
                  RealNumber x_min,
                  RealNumber y_min,
                  RealNumber x_max,
                  RealNumber y_max,
                  int number_x,
                  int number_y)
                : step_t(step_t),
                  step_x(step_x),
                  step_y(step_y),
                  x_min(x_min),
                  y_min(y_min),
                  x_max(x_max),
                  y_max(y_max),
                  number_x(number_x),
                  number_y() {}
};

// struct to hold cartesian coordinates and a value at that point  
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
class mesh {
    std::vector<DataPoint> data;
    size_t rows, cols;

public:
    // constructor
    mesh(size_t rows, size_t cols)
        : data(rows * cols), rows(rows), cols(cols) {}

    // copy constructor
    mesh(mesh &Mesh)
        : data(Mesh.data), rows(Mesh.rows), cols(Mesh.cols) {}

    DataPoint& operator()(size_t i, size_t j)
    {
        return data[i * cols + j];
    }
    const DataPoint& operator()(size_t i, size_t j) const
    {
        return data[i * cols + j];
    }

    size_t rowCount() const { return rows; }
    size_t colCount() const { return cols; }

    // method for regular grid construction
    void construct_grid(RealNumber dx, RealNumber dy, RealNumber x_min, RealNumber y_min)
    {
      for (size_t i = 0; i < rows; ++i)
      {
        for (size_t j = 0; j < cols; ++j)
        {
            data[i * cols + j].x = x_min + i * dx;
            data[i * cols + j].y = y_min + j * dy;
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
            file << data[i*cols + j].value << " "; 
          }
          file << std::endl;
        }
      }
      file.close();
    }
};

// impose initial conditions
void initial_conditions(mesh &Mesh)
{
  for (size_t i = 0; i < Mesh.rowCount(); i++)
  {
    for (size_t j = 0; j < Mesh.colCount(); j++)
    {
      Mesh(i,j).value = 100*exp(-(pow(Mesh(i,j).x+1,2.0)+pow(Mesh(i,j).y,2.0))/0.01);
    }
  }
}

// impose boundary conditions
void boundary_conditions(mesh &Mesh)
{
  for (size_t i = 0; i < Mesh.rowCount(); ++i)
  {
    for (size_t j = 0; j < Mesh.colCount(); ++j)
    {
      if (i == 0 || i == Mesh.rowCount() - 1 || j == 0 || j == Mesh.colCount() - 1)
      {
        Mesh(i,j).value = 0; 
      }
    }
  }
}

void print_mesh(mesh &Mesh)
{
  for (size_t i = 0; i < Mesh.rowCount(); i++)
  {
    for (size_t j = 0; j < Mesh.colCount(); j++)
    {
      std::cout << Mesh(i,j).value << " ";
    }
    std::cout << std::endl;
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
void lax_wendroff(mesh &Mesh, RealNumber dx, RealNumber dy, RealNumber dt)
{
  mesh U_n(Mesh);
  RealNumber lx, ly;
  lx = dt/dx;
  for(size_t i = 1; i < Mesh.rowCount()-1; i++)
  {
    for(size_t j = 1; j < Mesh.colCount()-1; j++)
    {
      ly = 2*Mesh(i,j).x*dt/dy;
      Mesh(i,j).value = U_n(i,j).value
              - lx/2*(U_n(i+1,j).value - U_n(i-1,j).value)
              - ly/2*(U_n(i,j+1).value - U_n(i,j-1).value)
              + pow(lx,2.0)/2*(U_n(i+1,j).value - 2*U_n(i,j).value + U_n(i-1,j).value)
              + pow(ly,2.0)/2*(U_n(i,j+1).value - 2*U_n(i,j).value + U_n(i,j-1).value)
              + ly*lx/4*(U_n(i+1,j+1).value - U_n(i-1,j+1).value - U_n(i+1,j-1).value + U_n(i-1,j-1).value);
    }
  }
}

RealNumber cfl_lax_wendroff(RealNumber dx, RealNumber dy)
{
  return 1.0/(1.0/dx+1.0/2.0*1.5*dy);
}

int main()
{
  // domain definition
  const RealNumber x_min = -1.5, x_max = 1.5, y_min = -1.5, y_max = 1.5;
  // number of grid points in the x, y direction respectively
  const int N_x = 0, N_y = 0;
  const RealNumber dx = (x_max-x_min)/N_x, dy = (y_max-y_min)/N_y;
  const RealNumber T = 2.f;
  RealNumber t = 0.f;
  const std::string sim_name;
  SimulationInfo sim_info();
  

  // initialize the mesh, U and file_name
  mesh Mesh(N_x, N_y); 
  std::string file_name;

  Mesh.construct_grid(dx, dy, x_min, y_min);
  initial_conditions(Mesh);
  Mesh.write_data("sim/output_t_0.00000.csv");
  RealNumber dt = cfl_lax_wendroff(dx, dy);

  // std::cout << "Enter simulation name (folder with the same name will be created): ";
  // std::cin << sim_name;

  // main loop
  while (t <= T)
  {
    t += dt;
    lax_wendroff(Mesh, sim_info);
    boundary_conditions(Mesh);
    
    // this is kinda awkward
    std::ostringstream fn;
    fn << "sim/output_t_" << std::fixed << std::setprecision(5) << t << ".csv";
    file_name = fn.str();

    Mesh.write_data(file_name);
  }
  return 0;
}
