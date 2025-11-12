#pragma once

#include <string>
#include <vector>
using RealNumber = float;

struct SimulationInfo
{
  RealNumber step_t, step_x, step_y, x_min, y_min, x_max, y_max;
  int number_x, number_y;

  // constructor
  SimulationInfo( const RealNumber step_t,
                  const RealNumber step_x,
                  const RealNumber step_y,
                  const RealNumber x_min,
                  const RealNumber y_min,
                  const RealNumber x_max,
                  const RealNumber y_max,
                  const int number_x,
                  const int number_y);
};

class Mesh
{
private:

  std::vector< RealNumber > data;
  std::vector< RealNumber > x_cord;
  std::vector< RealNumber > y_cord;
  size_t rows, cols;

public:

  Mesh( SimulationInfo sim_info ); 

  Mesh( const Mesh &mesh );

  RealNumber &value_ref( size_t i, size_t j );

  const RealNumber &value( size_t i, size_t j ) const;

  const RealNumber &x( size_t i, size_t j ) const;

  const RealNumber &y( size_t i, size_t j ) const;

  size_t getCols() const;
  size_t getRows() const;

  void construct_regular_grid( SimulationInfo sim_info );

  void write_data( std::string file_name );

};

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

