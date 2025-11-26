#pragma once

#include "Mesh.hpp"
#include "GeneralVector.hpp"
#include "SimulationInfo.hpp"
#include <cmath>

template < class RealNumber >
struct Lax_Friedrichs
{
  static void solve( Mesh< RealNumber > &mesh, SimulationInfo< RealNumber > sim_info)
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
  static RealNumber cfl( RealNumber dx, RealNumber dy )
  {
    RealNumber dt_limit_x = dx / 2.0f;
    RealNumber dt_limit_y = dy / 6.0f; 
    
    return std::min(dt_limit_x, dt_limit_y);
  }
};

template < class RealNumber >
struct Lax_Wendroff
{
  static void solve( Mesh< RealNumber > &mesh, SimulationInfo< RealNumber > sim_info)
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

  static RealNumber cfl( RealNumber dx, RealNumber dy )
  {
    return 1.0 / (1.0 / dx + 1.0 / 2.0 * 1.5 * dy);
  }

};

template < class RealNumber >
struct Upwind
{
  static void solve( Mesh< RealNumber > &mesh, SimulationInfo< RealNumber > sim_info)
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
  static RealNumber cfl( RealNumber dx, RealNumber dy )
  {
    RealNumber a = 1.0f, b = 1.0f;
    //podmínka: dt <= 1 / ( |a|/dx + |b|/dy )
    return 1.0f / (fabs(a)/dx + fabs(b)/dy);
  }
};

template < class RealNumber, class Scheme >
void NumericalSolver( Mesh<RealNumber > &mesh, SimulationInfo< RealNumber > sim_info )
{
  Scheme::solve( mesh, sim_info );
}

template < class RealNumber, class Scheme >
RealNumber CFL( RealNumber dx, RealNumber dy )
{
  Scheme::cfl( dx, dy );
}
