#include "functions.h"
#include <iostream>
#include <fstream>
#include <cmath>

// SimulationInfo

SimulationInfo::SimulationInfo(
                  const RealNumber step_t,
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
                  number_y(number_y)
  {
    std::cout << "Simulation info initialization..." << std::endl;
  }

// Mesh

Mesh::Mesh( SimulationInfo sim_info )
        : data(sim_info.number_x * sim_info.number_y),
          x_cord(sim_info.number_x * sim_info.number_y),
          y_cord(sim_info.number_x * sim_info.number_y),
          rows(sim_info.number_x), cols(sim_info.number_y)
{
  std::cout << "Mesh initialization..." << std::endl;
}

Mesh::Mesh( const Mesh &mesh )
        : data(mesh.data),
          x_cord(mesh.x_cord),
          y_cord(mesh.y_cord),
          rows(mesh.rows),
          cols(mesh.cols) {}

RealNumber& Mesh::value_ref( size_t i, size_t j )
{
  if ( i > rows || j > cols )
    {
      throw std::out_of_range("Index out of bounds");
    }

    return data[i * cols + j];
}

const RealNumber& Mesh::value( size_t i,size_t j ) const
{
  if ( i > rows || j > cols )
  {
    throw std::out_of_range("Index out of bounds");
  }

  return data[i * cols + j];
}

const RealNumber& Mesh::x( size_t i,size_t j ) const
{
  if ( i > rows || j > cols )
  {
    throw std::out_of_range("Index out of bounds");
  }

  return x_cord[i * cols + j];
}

const RealNumber& Mesh::y( size_t i,size_t j ) const
{
  if ( i > rows || j > cols )
  {
    throw std::out_of_range("Index out of bounds");
  }

  return y_cord[i * cols + j];
}

size_t Mesh::getCols() const { return cols; }
size_t Mesh::getRows() const { return rows; }

void Mesh::construct_regular_grid( SimulationInfo sim_info )
{
  for (size_t i = 0; i < rows; ++i)
  {
    for (size_t j = 0; j < cols; ++j)
    {
        x_cord[i * cols + j] = sim_info.x_min + i * sim_info.step_x;
        y_cord[i * cols + j] = sim_info.y_min + j * sim_info.step_y;
    }
  }
  std::cout << "Regular grid constructed..." << std::endl;
}

void Mesh::write_data( std::string file_name )
{
  std::ofstream file;
  file.open(file_name);

  if (file)
  {
    for (size_t i = 0; i < rows; i++)
    {
      for (size_t j = 0; j < cols; j++)
      {
        file << x_cord[i * cols + j] << "," << y_cord[i * cols + j] << "," << data[i * cols + j] << std::endl;
      }
    }
  }
  file.close();
}

// ------------------numerical schemes-----------------------------

// Lax-Friedrichs

void Lax_Friedrichs::solve( Mesh &mesh, SimulationInfo sim_info )
{
  Mesh U_n(mesh); 
  RealNumber dx, dy, dt, p = 1.0f;

  dx = sim_info.step_x;
  dy = sim_info.step_y;
  dt = sim_info.step_t;

  for(size_t i = 1; i < mesh.getRows() - 1; i++)
  {
      for(size_t j = 1; j < mesh.getCols() - 1; j++)
      {

          RealNumber x_j = mesh.x(i, j); 

          mesh.value_ref(i,j) = 
                1.0f/4.0f * ( U_n.value(i + 1, j) + U_n.value(i - 1, j) + U_n.value(i, j + 1) + U_n.value(i, j - 1))
              - dt * (( U_n.value(i + 1, j) - U_n.value(i - 1, j)) / ( 2.0f * dx ) // u_x
                  - 2.0f * p * x_j * ( U_n.value(i, j + 1) - U_n.value(i, j - 1)) / ( 2.0f * dy ) 
                );
      }
  }
}

RealNumber Lax_Friedrichs::cfl( RealNumber dx, RealNumber dy )
{
    RealNumber dt_limit_x = dx / 2.0f;
    RealNumber dt_limit_y = dy / 6.0f; 
    
    return std::min(dt_limit_x, dt_limit_y);
}

// Lax-Wendroff

void Lax_Wendroff::solve(Mesh &mesh, SimulationInfo sim_info)
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
      ly = 2*mesh.x(i,j)*dt/dy;
      mesh.value_ref(i,j) = U_n.value(i,j)
              - lx / 2 * ( U_n.value(i + 1,j) - U_n.value(i - 1,j) )
              - ly / 2 * ( U_n.value(i,j + 1) - U_n.value(i,j - 1) )
              + pow(lx,2.0) / 2 * ( U_n.value(i + 1,j) - 2 * U_n.value(i,j) + U_n.value(i - 1,j) )
              + pow(ly,2.0) / 2 * (U_n.value(i,j + 1) - 2 * U_n.value(i,j) + U_n.value(i,j - 1) )
              + ly * lx / 4 * (U_n.value(i + 1,j + 1) - U_n.value(i - 1,j + 1) - U_n.value(i + 1,j - 1) + U_n.value(i - 1,j - 1) );
    }
  }
}

RealNumber Lax_Wendroff::cfl(RealNumber dx, RealNumber dy)
{
  return 1.0 / (1.0 / dx + 1.0 / 2.0 * 1.5 * dy);
}

// Upwind

void Upwind::solve(Mesh &mesh, SimulationInfo sim_info)
{
    Mesh U_n(mesh);  // uložení stavu z předchozího kroku
    RealNumber dx = sim_info.step_x;
    RealNumber dy = sim_info.step_y;
    RealNumber dt = sim_info.step_t;

    // rychlosti šíření
    RealNumber a = 1.0f;   // směr proudění ve směru x
    RealNumber b = 1.0f;   // směr proudění ve směru y

    for (size_t i = 1; i < mesh.getRows() - 1; i++)
    {
        for (size_t j = 1; j < mesh.getCols() - 1; j++)
        {
            mesh.value_ref(i, j) =
                U_n.value(i, j)
                - a * (dt / dx) * (U_n.value(i, j) - U_n.value(i - 1, j))
                - b * (dt / dy) * (U_n.value(i, j) - U_n.value(i, j - 1));
        }
    }
}

RealNumber Upwind::cfl(RealNumber dx, RealNumber dy)
{
    RealNumber a = 1.0f, b = 1.0f;
    //podmínka: dt <= 1 / ( |a|/dx + |b|/dy )
    return 1.0f / (fabs(a)/dx + fabs(b)/dy);
}

// initial conditions

void My_Initial_Conditions::impose( Mesh &mesh )
{
  for (size_t i = 0; i < mesh.getRows(); i++)
  {
    for (size_t j = 0; j < mesh.getCols(); j++)
    {
      mesh.value_ref(i,j) = 100 * exp( - ( pow(mesh.x(i,j) + 1, 2.0) + pow(mesh.y(i,j), 2.0) ) / 0.01);
    }
  }
}

// boundary conditions

void Zeros::impose( Mesh &mesh )
{
  for (size_t i = 0; i < mesh.getRows(); ++i)
  {
    for (size_t j = 0; j < mesh.getCols(); ++j)
    {
      if (i == 0 || i == mesh.getRows() - 1 || j == 0 || j == mesh.getCols() - 1)
      {
        mesh.value_ref(i,j) = 0; 
      }
    }
  }
}
