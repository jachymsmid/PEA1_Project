#pragma once

#include <string>
#include <vector>
using RealNumber = float;


// numerical schemes

struct Lax_Friedrichs
{
  static void solve( Mesh &mesh, SimulationInfo sim_info);
  static RealNumber cfl( RealNumber dx, RealNumber dy );
};

struct Lax_Wendroff
{
  static void solve( Mesh &mesh, SimulationInfo sim_info);
  static RealNumber cfl( RealNumber dx, RealNumber dy );
};

struct Upwind
{
  static void solve( Mesh &mesh, SimulationInfo sim_info);
  static RealNumber cfl( RealNumber dx, RealNumber dy );
};

template < typename T >
void NumericalSolver( Mesh &mesh, SimulationInfo sim_info );
#include "templates/NumericalSolver.tpp"

template < typename T >
RealNumber CFL( RealNumber dx, RealNumber dy );
#include "templates/CFL.tpp"

// initial conditions

struct My_Initial_Conditions
{
  static void impose( Mesh &mesh );
};

template < typename T >
void InitialConditions( Mesh &mesh );
#include "templates/InitialConditions.tpp"

// boundary conditions

struct Zeros
{
  static void impose( Mesh &mesh );
};

template < typename T >
void BoundaryConditions( Mesh &mesh );
#include "templates/BoundaryConditions.tpp"

