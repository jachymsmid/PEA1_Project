#pragma once

#include <cmath>

template < class RealNumber, class Mesh >
struct MyInitialConditions
{
  static void impose( Mesh &mesh )
  {
    for (size_t i = 0; i < mesh.getRows(); i++)
    {
      for (size_t j = 0; j < mesh.getCols(); j++)
      {
        mesh.value_ref(i,j) = 100 * exp( - ( pow(mesh.x(i,j) + 1, 2.0) + pow(mesh.y(i,j), 2.0) ) / 0.01);
      }
    }
  }
};


template < class Mesh, class initialConditions >
void InitialConditions( Mesh &mesh )
{
  initialConditions::impose( mesh );
}
