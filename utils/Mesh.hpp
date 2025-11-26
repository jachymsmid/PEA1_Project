#pragma once

#include "GeneralVector.hpp"
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include "SimulationInfo.hpp"


template < class RealNumber >
class Mesh
{
private:

  std::vector< RealNumber > data;
  std::vector< RealNumber > x_cord;
  std::vector< RealNumber > y_cord;
  std::size_t rows, cols;

public:

  Mesh( SimulationInfo< RealNumber > sim_info )
    : data(sim_info.number_x * sim_info.number_y),
      x_cord(sim_info.number_x * sim_info.number_y),
      y_cord(sim_info.number_x * sim_info.number_y),
      rows(sim_info.number_x),
      cols(sim_info.number_y) {}

  Mesh( const Mesh &mesh )
    : data(mesh.data),
      x_cord(mesh.x_cord),
      y_cord(mesh.y_cord),
      rows(mesh.rows),
      cols(mesh.cols) {}

  RealNumber &value_ref( std::size_t i, std::size_t j ) { return data[i * cols + j]; }

  const RealNumber &value( size_t i, size_t j ) const { return data[i*cols + j]; }

  const RealNumber &x( size_t i, size_t j ) const { return x_cord[i*cols + j]; }

  const RealNumber &y( size_t i, size_t j ) const { return y_cord[i*cols + j]; }

  size_t getCols() const { return cols; }
  size_t getRows() const { return rows; }


  void construct_regular_grid( SimulationInfo< RealNumber > sim_info )
  {
    for (size_t i = 0; i < rows; ++i)
    {
      for (size_t j = 0; j < cols; ++j)
      {
          x_cord[i * cols + j] = sim_info.x_min + i * sim_info.step_x;
          y_cord[i * cols + j] = sim_info.y_min + j * sim_info.step_y;
      }
    }
  }

  void write_data( std::string file_name )
  {
    const int pocPoli = 1; //<------------------ TADY

    int pocBunek = (cols)*(rows);

    std::ofstream file;
    file.open( file_name );

    file << "# vtk DataFile Version 3.0" << std::endl;
    file << "vtk output" << std::endl;
    file << "ASCII" << std::endl;
    file << "DATASET STRUCTURED_GRID" << std::endl;
    // why are dimensions cols + 1 ??
    file << "DIMENSIONS "<< cols + 1 << " " << rows + 1 << " 1" << std::endl;
    file << "POINTS " << ( cols + 1 )*( rows + 1) << " float" << std::endl;

    for( size_t i=0; i < cols + 1; i++)
      for(float j=0; j < rows + 1; j++)
        file << x_cord[i] << " " << y_cord[j] << " 0" << std::endl;

    file << "CELL_DATA " << pocBunek << std::endl;
    file << "FIELD FieldData " << pocPoli << std::endl;

    for ( size_t k = 0; k < pocPoli; k++ )
    {
      file << "value " << k << pocBunek << " float" << std::endl;
      for( size_t i = 0; i < cols; i++)
        for( size_t j = 0; j < rows; j++)
            file << data[i*cols+j] << std::endl;
    }

    file.close();
  }

};
