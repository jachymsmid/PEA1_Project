#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <iomanip>
#include <sstream>


using RealNumber = float;

// point in 2D Cartesian coordinates
struct Point
{
  RealNumber x, y;
};


// impose initial conditions
void initial_conditions(std::vector<std::vector< RealNumber >> &U, std::vector<std::vector<Point>> &Mesh)
{
  for (size_t i = 0; i < U.size(); i++)
  {
    for (size_t j = 0; j < U[i].size(); j++)
    {
      U[i][j] = 100*exp(-(pow(Mesh[i][j].x+1,2.0)+pow(Mesh[i][j].y,2.0))/0.01);
    }
  }
}

// impose boundary conditions
void boundary_conditions(std::vector<std::vector<RealNumber>> &U)
{
  for (size_t i = 0; i < U.size(); ++i)
  {
    for (size_t j = 0; j < U[i].size(); ++j)
    {
      if (i == 0 || i == U.size() - 1 || j == 0 || j == U[i].size() - 1)
      {
        U[i][j] = 0; 
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

void print_field(std::vector<std::vector<RealNumber>> &U)
{
  for (size_t i = 0; i < U.size(); i++)
  {
    for (size_t j = 0; j < U[i].size(); j++)
    {
      std::cout << U[i][j] << " ";
    }
    std::cout << std::endl;
  }
}

// write the data to a file
void write_data(std::vector<std::vector<RealNumber>> &U, std::string file_name)
{
  std::ofstream file;
  file.open(file_name);

  if (file)
  {
    for (size_t i = 0; i < U.size(); i++)
    {
      for (size_t j = 0; j < U[i].size(); j++)
      {
        file << U[i][j] << " "; 
      }
      file << std::endl;
    }
  }
  file.close();
}

// numerical schemes
void lax_friedrichs(){}

void upwind(){}

void lax_wendroff(std::vector<std::vector<RealNumber>> &U, std::vector<std::vector<Point>> &Mesh)
{
  for(size_t i = 0; i < U.size()-1; i++)
  {
    for(size_t j = 0; j < U[i].size()-1; j++)
    {
      U[i][j] = 2*i*j;
    }
  }
}

RealNumber cfl_lax_wendroff(RealNumber dx, RealNumber dy)
{
  return 1/(dx*dy);
}

int main()
{
  // domain definition
  const RealNumber x_min = -1.5, x_max = 1.5, y_min = -1.5, y_max = 1.5;
  // number of grid points in the x, y direction respectively
  const int N_x = 10, N_y = 10;
  const RealNumber dx = (x_max-x_min)/N_x, dy = (y_max-y_min)/N_y;
  const int T = 2;
  RealNumber t = 0.f;

  // initialize the mesh, U and file_name
  std::vector<std::vector<RealNumber>> U(N_x, std::vector<RealNumber>(N_y)); 
  std::vector<std::vector<Point>> Mesh (N_x, std::vector<Point>(N_y));
  std::string file_name;
  std::ostringstream fn;

  mesh(dx, dy, x_min, y_min, Mesh);
  initial_conditions(U, Mesh);
  write_data(U, "sim/t_0.dat");
  RealNumber dt = cfl_lax_wendroff(dx, dy);

  // main loop
  while (t < T-1)
  {
    t += dt; 
    lax_wendroff(U, Mesh);
    boundary_conditions(U);
    fn << "sim/output_t_" << std::fixed << std::setprecision(5) << t << ".dat";
    file_name = fn.str(); 
    write_data(U, file_name);
  }

  return 0;
}
