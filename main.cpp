#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <iomanip>
#include <sstream>

// TODO:  - [ ] implement better mesh representation
//        - [x] implement flattened arrays 
//        - [ ] upwind method
//        - [ ] lax-friedrichs method
//        - [x] lax-wendroff method

using RealNumber = float;

// point in 2D Cartesian coordinates
struct Point
{
  RealNumber x, y;
};

// matrix class for faster 2d vector implementation
class Matrix {
    std::vector<RealNumber> data;
    size_t rows, cols;

public:
    Matrix(size_t rows, size_t cols)
        : data(rows * cols), rows(rows), cols(cols) {}

    RealNumber& operator()(size_t i, size_t j)
    {
        return data[i * cols + j];
    }
    const RealNumber& operator()(size_t i, size_t j) const
    {
        return data[i * cols + j];
    }

    size_t rowCount() const { return rows; }
    size_t colCount() const { return cols; }
};

// impose initial conditions
void initial_conditions(Matrix &U, std::vector<std::vector<Point>> &Mesh)
{
  for (size_t i = 0; i < U.rowCount(); i++)
  {
    for (size_t j = 0; j < U.colCount(); j++)
    {
      U(i,j) = 100*exp(-(pow(Mesh[i][j].x+1,2.0)+pow(Mesh[i][j].y,2.0))/0.01);
    }
  }
}

// impose boundary conditions
void boundary_conditions(Matrix &U)
{
  for (size_t i = 0; i < U.rowCount(); ++i)
  {
    for (size_t j = 0; j < U.colCount(); ++j)
    {
      if (i == 0 || i == U.rowCount() - 1 || j == 0 || j == U.colCount() - 1)
      {
        U(i,j) = 0; 
      }
    }
  }
}

// construct the mesh
void mesh(RealNumber dx, RealNumber dy, RealNumber x_min, RealNumber y_min, std::vector<std::vector<Point>> &Mesh)
{
  for (size_t i = 0; i < Mesh.size(); ++i)
  {
    for (size_t j = 0; j < Mesh[i].size(); ++j)
    {
        Mesh[i][j].x = x_min + i * dx;
        Mesh[i][j].y = y_min + j * dy;
    }
  }
}

// print mesh for troubleshooting reasons
void print_mesh(std::vector<std::vector<Point>> &Mesh)
{
  for (size_t i = 0; i < Mesh.size(); ++i)
  {
    for (size_t j = 0; j < Mesh[i].size(); ++j)
    {
      std::cout << "(" << Mesh[i][j].x << ", " << Mesh[i][j].y << ") ";
    }
    std::cout << std::endl;
  }
}

void print_field(Matrix &U)
{
  for (size_t i = 0; i < U.rowCount(); i++)
  {
    for (size_t j = 0; j < U.colCount(); j++)
    {
      std::cout << U(i,j) << " ";
    }
    std::cout << std::endl;
  }
}

// write the data to a file
void write_data(Matrix &U, std::string file_name)
{
  std::ofstream file;
  file.open(file_name);

  if (file)
  {
    for (size_t i = 0; i < U.rowCount(); i++)
    {
      for (size_t j = 0; j < U.colCount(); j++)
      {
        file << U(i,j) << " "; 
      }
      file << std::endl;
    }
  }
  file.close();
}

// numerical schemes
void lax_friedrichs()
{
}

void upwind()
{
}

void lax_wendroff(Matrix &U, std::vector<std::vector<Point>> &Mesh, RealNumber dx, RealNumber dy, RealNumber dt)
{
  Matrix U_n = U;
  RealNumber lx, ly;
  lx = dt/dx;
  for(size_t i = 1; i < U.rowCount()-1; i++)
  {
    for(size_t j = 1; j < U.colCount()-1; j++)
    {
      ly = 2*Mesh[i][j].x*dt/dy;
      U(i,j) = U_n(i,j)
              - lx/2*(U_n(i+1,j) - U_n(i-1,j))
              - ly/2*(U_n(i,j+1) - U_n(i,j-1))
              + pow(lx,2.0)/2*(U_n(i+1,j) - 2*U_n(i,j) + U_n(i-1,j))
              + pow(ly,2.0)/2*(U_n(i,j+1) - 2*U_n(i,j) + U_n(i,j-1))
              + ly*lx/4*(U_n(i+1,j+1) - U_n(i-1,j+1) - U_n(i+1,j-1) + U_n(i-1,j-1));
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

  // initialize the mesh, U and file_name
  Matrix U(N_x, N_y); 
  std::vector<std::vector<Point>> Mesh (N_x, std::vector<Point>(N_y));
  std::string file_name;

  mesh(dx, dy, x_min, y_min, Mesh);
  initial_conditions(U, Mesh);
  write_data(U, "sim/output_t_0.00000.dat");
  RealNumber dt = cfl_lax_wendroff(dx, dy);

  // std::cout << "Enter simulation name (folder with the same name will be created): ";
  // std::cin << sim_name;

  // main loop
  while (t <= T)
  {
    t += dt;
    lax_wendroff(U, Mesh, dx, dy, dt);
    boundary_conditions(U);
    std::ostringstream fn;
    fn << "sim/output_t_" << std::fixed << std::setprecision(5) << t << ".dat";
    file_name = fn.str();
    write_data(U, file_name);
  }
  return 0;
}
