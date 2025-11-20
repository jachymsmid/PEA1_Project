#pragma once

#include <string>
#include <vector>
using RealNumber = float;


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

