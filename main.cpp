#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>

// point in 2D Cartesian coordinates
struct Point
{
  double x, y;
};



// impose initial conditions
void initial_conditions(std::vector<std::vector<double>> &U, std::vector<std::vector<Point>> &Mesh, int N_x, int N_y)
{
  for (int i = 0; i < N_x; i++)
  {
    for (int j = 0; j < N_y; j++)
    {
      U[i][j] = 100*exp((-pow(Mesh[i][j].x+1,2.0)+pow(Mesh[i][j].y,2.0))/0.01);
    }
  }
}

// impose boundary conditions
void boundary_conditions(std::vector<std::vector<double>> &U, int N_x, int N_y)
{
  for (int i = 0; i < N_x; i++)
  {
    U[i][0] = 0.0;
    U[i][-1] = 0.0;
  }
  for (int j = 0; j <N_y; j++)
  {
    U[0][j] = 0.0;
    U[-1][j] = 0.0;
  }
}

// construct the mesh
void mesh(int N_x, int N_y, double x_min, double x_max, double y_min, double y_max, std::vector<std::vector<Point>> &Mesh)
{
  double dx = (x_max - x_min)/(N_x-1);
  double dy = (y_max - y_min)/(N_y-1);
  for (int i = 0; i < N_x; ++i)
  {
    for (int j = 0; j < N_y; ++j)
    {
        Mesh[i][j].x = x_min + i * dx;
        Mesh[i][j].y = y_min + j * dy;
    }
  }
}

// print mesh for troubleshooting reasons
void print_mesh(int N_x, int N_y, std::vector<std::vector<Point>> &Mesh)
{
  for (int i = 0; i < N_x; ++i)
  {
    for (int j = 0; j < N_y; ++j)
    {
      std::cout << "(" << Mesh[i][j].x << ", " << Mesh[i][j].y << ") ";
    }
    std::cout << std::endl;
  }
}

void print_field(std::vector<std::vector<double>> &U, int N_x, int N_y)
{
  for (int i = 0; i < N_x; i++)
  {
    for (int j = 0; j < N_y; j++)
    {
      std::cout << U[i][j] << " ";
    }
    std::cout << std::endl;
  }
}

// write the data to a file
void write_data(){}

// numerical schemes
void lax_friedrichs(){}
void lax_wendorf(){}
void upwind(){}

int main()
{
  // domain definition
  const double x_min = 0, x_max = 2, y_min = 0, y_max = 2;
  // number of grid points in the x, y direction respectively
  const int N_x = 5, N_y = 5;

  
  std::vector<std::vector<Point>> Mesh (N_x, std::vector<Point>(N_y));
  mesh(N_x, N_y, x_min, x_max, y_min, y_max, Mesh);
  //print_mesh(N_x, N_y, Mesh);

  std::vector<std::vector<double>> U(N_x, std::vector<double>(N_y));
  initial_conditions(U, Mesh, N_x, N_y);
  print_field(U, N_x, N_y);

  return 0;
}
