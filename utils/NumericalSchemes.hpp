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
    RealNumber dt = sim_info.step_t;
    RealNumber a = sim_info.x_speed;   // advekce ve směru x
    RealNumber b = sim_info.y_speed;   // advekce ve směru y

    for(size_t i = 1; i < mesh.getRows() - 1; i++)
    {
        for(size_t j = 1; j < mesh.getCols() - 1; j++)
        {
            RealNumber dx_local = mesh.x(i+1,j) - mesh.x(i-1,j);
            RealNumber dy_local = mesh.y(i,j+1) - mesh.y(i,j-1);
            RealNumber x_j = mesh.x(i, j); 

            RealNumber avg = (U_n.value(i + 1, j) + U_n.value(i - 1, j) + U_n.value(i, j + 1) + U_n.value(i, j - 1)) / RealNumber(4.0);

            mesh.value_ref(i, j) = avg - dt * (a * (U_n.value(i+1,j) - U_n.value(i-1,j)) / dx_local 
                                             + b * (U_n.value(i,j+1) - U_n.value(i,j-1)) / dy_local);

        }
    }
  }
// CFL podmínka pro maximální krok
static RealNumber CFL(Mesh<RealNumber> &mesh, SimulationInfo<RealNumber> sim_info)
{
    RealNumber max_dx = 0;
    RealNumber max_dy = 0;

    for (size_t i = 1; i < mesh.getRows() - 1; i++)
    {
        for (size_t j = 1; j < mesh.getCols() - 1; j++)
        {
            RealNumber dx_local = mesh.x(i+1, j) - mesh.x(i-1, j);
            RealNumber dy_local = mesh.y(i, j+1) - mesh.y(i, j-1);

            if (std::abs(dx_local) > max_dx)
                max_dx = std::abs(dx_local);

            if (std::abs(dy_local) > max_dy)
                max_dy = std::abs(dy_local);
        }
    }

    RealNumber a = sim_info.x_speed;
    RealNumber b = sim_info.y_speed;

    RealNumber dt = RealNumber(1.0) /
                   (std::abs(a) / max_dx + std::abs(b) / max_dy);

    return dt;
}
};

template < class RealNumber >
struct Lax_Wendroff
{
  static void solve( Mesh< RealNumber > &mesh, SimulationInfo< RealNumber > sim_info)
  {
    Mesh U_n(mesh);
    RealNumber dt = sim_info.step_t;

    RealNumber a = sim_info.x_speed; // rychlost ve směru x
    RealNumber b = sim_info.y_speed; // rychlost ve směru y

    for(size_t i = 1; i < mesh.getRows() - 1; i++)
    {
      for(size_t j = 1; j < mesh.getCols() - 1; j++)
      {
        RealNumber dx_local = mesh.x(i + 1, j) - mesh.x(i - 1, j);
        RealNumber dy_local = mesh.y(i, j + 1) - mesh.y(i, j - 1);

        RealNumber lx = a * dt / dx_local;
        RealNumber ly = b * dt / dy_local;

        mesh.value_ref(i,j) = U_n.value(i,j)
                - lx / 2 * ( U_n.value(i + 1,j) - U_n.value(i - 1,j) )
                - ly / 2 * ( U_n.value(i,j + 1) - U_n.value(i,j - 1) )
                + pow(lx,2.0) / 2 * (U_n.value(i + 1,j) - 2 * U_n.value(i,j) + U_n.value(i - 1,j) )
                + pow(ly,2.0) / 2 * (U_n.value(i,j + 1) - 2 * U_n.value(i,j) + U_n.value(i,j - 1) )
                + pow(ly,2.0) / 4 * (U_n.value(i + 1,j + 1) - U_n.value(i - 1,j + 1) - U_n.value(i + 1,j - 1) + U_n.value(i - 1,j - 1) );
      }
    }
  }

static RealNumber cfl(Mesh<RealNumber> &mesh, SimulationInfo<RealNumber> sim_info)
{
    RealNumber max_dx = 0;
    RealNumber max_dy = 0;

    for (size_t i = 1; i < mesh.getRows() - 1; i++)
    {
        for (size_t j = 1; j < mesh.getCols() - 1; j++)
        {
            RealNumber dx_local = mesh.x(i+1, j) - mesh.x(i-1, j);
            RealNumber dy_local = mesh.y(i, j+1) - mesh.y(i, j-1);

            if (std::abs(dx_local) > max_dx)
                max_dx = std::abs(dx_local);

            if (std::abs(dy_local) > max_dy)
                max_dy = std::abs(dy_local);
        }
    }

    RealNumber a = sim_info.x_speed;
    RealNumber b = sim_info.y_speed;

    RealNumber dt = RealNumber(1.0) /
                   (std::abs(a) / max_dx + std::abs(b) / max_dy);

    return dt;
}
};

template < class RealNumber >
struct Upwind
{
  static void solve( Mesh< RealNumber > &mesh, SimulationInfo< RealNumber > sim_info)
  {
    Mesh U_n(mesh);  // uložení stavu z předchozího kroku
    RealNumber dt = sim_info.step_t;

    // rychlosti šíření
    RealNumber a = sim_info.x_speed;   // směr proudění ve směru x
    RealNumber b = sim_info.y_speed;   // směr proudění ve směru y

    for (size_t i = 1; i < mesh.getRows() - 1; i++)
    {
        for (size_t j = 1; j < mesh.getCols() - 1; j++)
        {
            RealNumber dx_local = mesh.x(i,j)-mesh.x(i-1,j);
            RealNumber dy_local = mesh.y(i,j)-mesh.y(i,j-1);
            mesh.value_ref(i, j) =
                U_n.value(i, j)
                - a * (dt / dx_local) * (U_n.value(i, j) - U_n.value(i - 1, j))
                - b * (dt / dy_local) * (U_n.value(i, j) - U_n.value(i, j - 1));
        }
    }
  }
  

static RealNumber cfl(Mesh<RealNumber> &mesh, SimulationInfo<RealNumber> sim_info)
{
    RealNumber max_dx = 0;
    RealNumber max_dy = 0;

    for (size_t i = 1; i < mesh.getRows() - 1; i++)
    {
        for (size_t j = 1; j < mesh.getCols() - 1; j++)
        {
            RealNumber dx_local = mesh.x(i+1, j) - mesh.x(i-1, j);
            RealNumber dy_local = mesh.y(i, j+1) - mesh.y(i, j-1);

            if (std::abs(dx_local) > max_dx)
                max_dx = std::abs(dx_local);

            if (std::abs(dy_local) > max_dy)
                max_dy = std::abs(dy_local);
        }
    }

    RealNumber a = sim_info.x_speed;
    RealNumber b = sim_info.y_speed;

    RealNumber dt = RealNumber(1.0) /
                   (std::abs(a) / max_dx + std::abs(b) / max_dy);

    return dt;
}
};

template < class RealNumber, class Scheme >
void NumericalSolver( Mesh<RealNumber > &mesh, SimulationInfo< RealNumber > sim_info )
{
  Scheme::solve( mesh, sim_info );
}

template < class RealNumber, class Scheme >
RealNumber CFL(Mesh<RealNumber> &mesh, SimulationInfo<RealNumber> sim_info)
{
  return Scheme::cfl(mesh, sim_info);
}
