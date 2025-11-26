#pragma once

#include "Mesh.hpp"

template < class RealNumber >
struct Zeros
{
  void impose( Mesh< RealNumber > &mesh )
  {
    for (size_t i = 0; i < mesh.getRows(); ++i)
    {
      for (size_t j = 0; j < mesh.getCols(); ++j)
      {
        // this condition is wrong because of different indexing
        if (i == 0 || i == mesh.getRows() - 1 || j == 0 || j == mesh.getCols() - 1)
        {
          mesh.value_ref(i,j) = 0; 
        }
      }
    }
  }
};

template < class RealNumber, class initialConditions >
void InitialConditions( Mesh< RealNumber > &mesh )
{
  initialConditions::impose( mesh );
}
