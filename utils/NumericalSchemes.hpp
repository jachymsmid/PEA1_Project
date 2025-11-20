#pragma once

#include <cstddef>

template< class RealNumber, class Mesh, class SimulationInfo >
struct Lax_Friedrichs
{
  void solve( Mesh &mesh, SimulationInfo sim_info)
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
  RealNumber cfl( const SimulationInfo &sim_info )
  {
    // will the compiler create new variables here? or will it substitute?
    RealNumber dx = sim_info.step_x;
    RealNumber dy = sim_info.step_y;
    RealNumber a = sim_info.x_speed;
    RealNumber b = sim_info.y_speed;

    RealNumber dt_limit_x = dx / a;
    RealNumber dt_limit_y = dy / b; 

    RealNumber min = dt_limit_x < dt_limit_y ? dt_limit_x : dt_limit_y;
    
    return min;
  }
};

template< class RealNumber, class Mesh, class SimulationInfo >
struct Lax_Wendroff
{
  void solve( Mesh &mesh, SimulationInfo sim_info)
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
  // wrong
  RealNumber cfl( RealNumber dx, RealNumber dy )
  {
    return 1.0 / (1.0 / dx + 1.0 / 2.0 * 1.5 * dy);
  }
};

template< class RealNumber, class Mesh, class SimulationInfo >
struct Upwind
{
  static void solve( Mesh &mesh, SimulationInfo sim_info);
  static RealNumber cfl( RealNumber dx, RealNumber dy );
};

