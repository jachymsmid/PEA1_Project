#include "functions.h"
#include <iostream>
#include <fstream>
#include <cmath>

#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkXMLImageDataWriter.h>

// SimulationInfo




// ------------------numerical schemes-----------------------------

// Lax-Friedrichs

void Lax_Friedrichs::solve( Mesh &mesh, SimulationInfo sim_info )


RealNumber Lax_Friedrichs::cfl( RealNumber dx, RealNumber dy )


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
